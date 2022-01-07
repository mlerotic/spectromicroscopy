# 
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
# 
#   Copyright (C) 2011 Mirna Lerotic, 2nd Look
#   http://2ndlookconsulting.com
#   License: GNU GPL v3
#
#   Mantis is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   any later version.
#
#   Mantis is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details <http://www.gnu.org/licenses/>.


from __future__ import division
from __future__ import print_function

import numpy as np
import scipy as sp
import scipy.interpolate
import scipy.ndimage
import h5py
import datetime
import os

from .file_plugins import file_stk
from .file_plugins import file_sdf
from .file_plugins import file_xrm
from .file_plugins import file_ncb
from .file_plugins import file_dataexch_hdf5

from . import data_struct


# ----------------------------------------------------------------------
class data:
    def __init__(self, data_struct):
        self.data_struct = data_struct
        self.i0_dwell = None
        self.i0data = np.zeros(1)
        self.n_ev = 0
        self.n_theta = 0
        self.shifts = []
        self.stack4D = None

    # ----------------------------------------------------------------------
    def new_data(self):
        self.n_cols = 0
        self.n_rows = 0
        self.n_ev = 0

        self.x_dist = 0
        self.y_dist = 0

        self.x_start = 0
        self.x_stop = 0
        self.y_start = 0
        self.y_stop = 0
        self.x_pxsize = 0
        self.y_pxsize = 0
        self.squarepx = True

        self.i0_dwell = None

        self.ev = 0
        self.absdata = 0

        self.i0data = np.zeros(1)
        self.evi0 = 0

        self.od = 0
        self.od3d = 0

        self.xshifts = 0
        self.yshifts = 0
        self.shifts = []

        self.stack4D = None
        self.n_theta = 0
        self.theta = 0
        self.od4D = 0

        self.data_struct.spectromicroscopy.normalization.white_spectrum = None
        self.data_struct.spectromicroscopy.normalization.white_spectrum_energy = None
        self.data_struct.spectromicroscopy.normalization.white_spectrum_energy_units = None

        self.data_struct.spectromicroscopy.optical_density = None

    # ----------------------------------------------------------------------
    def read_stk_i0(self, filename, extension):
        if extension == '.xas':
            file_stk.read_stk_i0_xas(self, filename)
        elif extension == '.csv':
            file_stk.read_stk_i0_csv(self, filename)

        self.calculate_optical_density()
        self.fill_h5_struct_normalization()

    # ----------------------------------------------------------------------
    def read_sdf_i0(self, filename):
        file_sdf.read_sdf_i0(self, filename)
        self.calculate_optical_density()

        self.fill_h5_struct_normalization()

    # ----------------------------------------------------------------------
    def read_xrm_ReferenceImages(self, filenames):

        self.calculate_optical_density_from_refimgs(filenames)

        self.fill_h5_struct_normalization()

    # ----------------------------------------------------------------------
    def read_h54D(self, filename):

        file_dataexch_hdf5.read(filename, self)

        if self.data_struct.spectromicroscopy.normalization.white_spectrum is not None:
            self.calculate_optical_density()
            self.fill_h5_struct_normalization()

        self.scale_bar()

    # ----------------------------------------------------------------------
    def read_ncb4D(self, filenames):
        self.new_data()
        file_ncb.read_ncb4D(self, filenames)
        now = datetime.datetime.now()

        self.data_struct.implements = 'information:exchange:spectromicroscopy'
        self.data_struct.version = '1.0'

        self.data_struct.information.file_creation_datetime = now.strftime("%Y-%m-%dT%H:%M")
        self.data_struct.information.comment = 'Converted in Mantis'

        self.data_struct.exchange.data = self.stack4D
        self.data_struct.exchange.data_signal = 1
        self.data_struct.exchange.data_axes = 'x:y:energy:theta'

        self.data_struct.exchange.theta = np.array(self.theta)
        self.data_struct.exchange.theta_units = 'degrees'

        self.data_struct.exchange.x = self.x_dist
        self.data_struct.exchange.y = self.y_dist

        self.scale_bar()

    # ----------------------------------------------------------------------
    def read_ncb4Denergy(self, filename):

        f = open(str(filename), 'rU')

        elist = []

        for line in f:
            if line.startswith('*'):
                if 'Common name' in line:
                    spectrum_common_name = line.split(':')[-1].strip()

            else:
                e, = [float(x) for x in line.split()]
                elist.append(e)

        self.ev = np.array(elist)

        f.close()

        self.n_ev = self.ev.size
        self.data_struct.exchange.energy = self.ev
        self.data_struct.exchange.energy_units = 'ev'

    # ----------------------------------------------------------------------
    def read_dpt(self, filename):
        self.new_data()

        n_rows = 11
        n_cols = 8

        imgstack = np.zeros((n_rows, n_cols))

        f = open(str(filename), 'r')

        elist = []

        for line in f:
            if line.startswith("*"):
                pass
            else:
                x = line.split(',')
                e = float(x[0])
                x = x[1:]
                data = []
                for i in range(len(x)):
                    data.append(float(x[i]))
                elist.append(e)
                data = np.array(data)
                data = np.reshape(data, (n_rows, n_cols), order='F')
                imgstack = np.dstack((imgstack, data))

        imgstack = imgstack[:, :, 1:]

        f.close()

        self.n_cols = imgstack.shape[0]
        self.n_rows = imgstack.shape[1]
        self.n_ev = imgstack.shape[2]

        pixelsize = 1
        # Since we do not have a scanning microscope we fill the x_dist and y_dist from pixel_size
        self.x_dist = np.arange(np.float(self.n_cols)) * pixelsize
        self.y_dist = np.arange(np.float(self.n_rows)) * pixelsize

        self.ev = np.array(elist)

        msec = np.ones((self.n_ev))

        self.data_dwell = msec

        self.absdata = imgstack

        # Check if the energies are consecutive, if they are not sort the data
        sort = 0
        for i in range(self.n_ev - 1):
            if self.ev[i] > self.ev[i + 1]:
                sort = 1
                break
        if sort == 1:
            sortind = np.argsort(self.ev)
            self.ev = self.ev[sortind]
            self.absdata = self.absdata[:, :, sortind]

        #         self.original_n_cols = imgstack.shape[0]
        #         self.original_n_rows = imgstack.shape[1]
        #         self.original_n_ev = imgstack.shape[2]
        #         self.original_ev = self.ev.copy()
        #         self.original_absdata = self.absdata.copy()

        self.fill_h5_struct_from_stk()

        self.scale_bar()

        # Fix the normalization
        self.evi0 = self.ev.copy()
        self.i0data = np.ones(self.n_ev)

        self.i0_dwell = self.data_dwell

        self.fill_h5_struct_normalization()

        # Optical density does not have to be calculated - use raw data

        self.od3d = self.absdata.copy()

        self.od = np.reshape(self.od3d, (n_rows * n_cols, self.n_ev), order='F')

    # ----------------------------------------------------------------------
    def fill_h5_struct_from_stk(self):

        now = datetime.datetime.now()

        self.data_struct.implements = 'information:exchange:spectromicroscopy'
        self.data_struct.version = '1.0'

        self.data_struct.information.file_creation_datetime = now.strftime("%Y-%m-%dT%H:%M")
        self.data_struct.information.comment = 'Converted in Mantis'

        self.data_struct.exchange.data = self.absdata
        self.data_struct.exchange.data_signal = 1
        self.data_struct.exchange.data_axes = 'x:y:energy'

        self.data_struct.exchange.energy = self.ev
        self.data_struct.exchange.energy_units = 'ev'

        self.data_struct.exchange.x = self.x_dist
        self.data_struct.exchange.y = self.y_dist

    # ----------------------------------------------------------------------
    def fill_h5_struct_normalization(self):

        self.data_struct.spectromicroscopy.normalization.white_spectrum = self.i0data
        self.data_struct.spectromicroscopy.normalization.white_spectrum_energy = self.evi0
        self.data_struct.spectromicroscopy.normalization.white_spectrum_energy_units = 'eV'

        if self.stack4D is None:
            self.data_struct.spectromicroscopy.optical_density = self.od
        else:
            self.data_struct.spectromicroscopy.optical_density = self.od4D

    # ----------------------------------------------------------------------
    def calc_histogram(self):
        # calculate average flux for each pixel
        self.averageflux = np.mean(self.absdata, axis=2)
        self.histogram = self.averageflux

        return

    # ----------------------------------------------------------------------
    def i0_from_histogram(self, i0_indices):

        self.evi0hist = self.ev.copy()

        # i0_indices = np.where((fluxmin<=self.averageflux)&(self.averageflux<=fluxmax))

        self.evi0 = self.ev.copy()

        self.i0_dwell = self.data_dwell

        if self.stack4D is None:
            self.i0datahist = np.zeros((self.n_ev))
            self.i0data = self.i0datahist
            if np.any(i0_indices):
                invnumel = 1. / self.averageflux[i0_indices].shape[0]
                for ie in range(self.n_ev):
                    thiseng_abs = self.absdata[:, :, ie]
                    self.i0datahist[ie] = np.sum(thiseng_abs[i0_indices]) * invnumel

            self.calculate_optical_density()

        else:
            self.i0datahist = np.zeros((self.n_ev, self.n_theta))
            self.i0data = self.i0datahist
            self.od4D = np.zeros((self.n_cols, self.n_rows, self.n_ev, self.n_theta))

            if np.any(i0_indices):
                invnumel = 1. / self.averageflux[i0_indices].shape[0]
            else:
                return

            for i in range(self.n_theta):
                for ie in range(self.n_ev):
                    thiseng_abs = self.stack4D[:, :, ie, i]
                    self.i0datahist[ie, i] = np.sum(thiseng_abs[i0_indices]) * invnumel

            self.calculate_optical_density_4D()

        self.fill_h5_struct_normalization()

        return

    # ----------------------------------------------------------------------
    def UsePreNormalizedData(self):

        self.evi0 = self.ev.copy()
        self.i0data = np.ones(self.n_ev)

        self.i0_dwell = self.data_dwell

        self.od = np.empty((self.n_cols, self.n_rows, self.n_ev))
        for i in range(self.n_ev):
            self.od[:, :, i] = self.absdata[:, :, i]

        self.od3d = self.od.copy()

        n_pixels = self.n_cols * self.n_rows
        # Optical density matrix is rearranged into n_pixelsxn_ev
        self.od = np.reshape(self.od, (n_pixels, self.n_ev), order='F')

        if self.stack4D is not None:
            self.od4D = self.stack4D.copy()

        self.fill_h5_struct_normalization()

        return

    # ----------------------------------------------------------------------
    def set_i0(self, i0data, evdata):

        self.evi0 = evdata
        self.i0data = i0data

        self.i0_dwell = self.data_dwell

        self.calculate_optical_density()

        self.fill_h5_struct_normalization()

        return

    # ----------------------------------------------------------------------
    def reset_i0(self):

        self.i0_dwell = None

        self.i0data = 0
        self.evi0 = 0

        self.od = 0
        self.od3d = 0

        self.data_struct.spectromicroscopy.normalization.white_spectrum = None
        self.data_struct.spectromicroscopy.normalization.white_spectrum_energy = None
        self.data_struct.spectromicroscopy.normalization.white_spectrum_energy_units = None

        self.data_struct.spectromicroscopy.optical_density = None

    # ----------------------------------------------------------------------
    # Normalize the data: calculate optical density matrix D
    def calculate_optical_density(self):

        if self.stack4D is not None:
            self.calculate_optical_density_4D()
            return

        n_pixels = self.n_cols * self.n_rows
        self.od = np.empty((self.n_cols, self.n_rows, self.n_ev))

        # little hack to deal with rounding errors
        self.evi0[self.evi0.size - 1] += 0.001
        self.evi0[0] -= 0.001
        if len(self.evi0) > 3: # >3 is needed to avoid boundary error!
            fi0int = scipy.interpolate.interp1d(self.evi0.astype(np.double), self.i0data.astype(np.double),
                                                kind='cubic', bounds_error=False, fill_value=0.0)
        elif len(self.evi0) > 1: # use linear interpolation when there are fewer points
            fi0int = scipy.interpolate.interp1d(self.evi0.astype(np.double), self.i0data.astype(np.double),
                                                bounds_error=False, fill_value=0.0)
        else: # use constant value when only a single value is available
            fi0int = lambda x: self.i0data.astype(np.double)
        
        i0 = fi0int(self.ev)

        if (self.data_dwell is not None) and (self.i0_dwell is not None):
            i0 = i0 * (self.data_dwell / self.i0_dwell)

        # zero out all negative values in the image stack
        negative_indices = np.where(self.absdata <= 0)
        if negative_indices:
            self.absdata[negative_indices] = 0.01

        for i in range(self.n_ev):
            self.od[:, :, i] = - np.log(self.absdata[:, :, i] / i0[i])

        # clean up the result
        nan_indices = np.where(np.isfinite(self.od) == False)
        if nan_indices:
            self.od[nan_indices] = 0

        self.od3d = self.od.copy()

        # Optical density matrix is rearranged into n_pixelsxn_ev
        self.od = np.reshape(self.od, (n_pixels, self.n_ev), order='F')

        return

    # ----------------------------------------------------------------------
    # Normalize the data: calculate optical density matrix D
    def calculate_optical_density_4D(self):

        n_pixels = self.n_cols * self.n_rows
        self.od4D = np.zeros((self.n_cols, self.n_rows, self.n_ev, self.n_theta))

        # little hack to deal with rounding errors
        self.evi0[self.evi0.size - 1] += 0.001

        self.i0data = np.array(self.i0data)
        i0dims = self.i0data.shape

        for ith in range(self.n_theta):
            self.od = np.empty((self.n_cols, self.n_rows, self.n_ev))

            if len(i0dims) == 2:
                self.i0data = self.i0datahist[:, ith]

            if len(self.evi0) > 2:
                fi0int = scipy.interpolate.interp1d(self.evi0, self.i0data, kind='cubic', bounds_error=False,
                                                    fill_value=0.0)
            else:
                fi0int = scipy.interpolate.interp1d(self.evi0, self.i0data, bounds_error=False, fill_value=0.0)
            i0 = fi0int(self.ev)

            if (self.data_dwell is not None) and (self.i0_dwell is not None):
                i0 = i0 * (self.data_dwell / self.i0_dwell)

            # zero out all negative values in the image stack
            negative_indices = np.where(self.stack4D <= 0)
            if negative_indices:
                self.stack4D[negative_indices] = 0.01

            for i in range(self.n_ev):
                self.od[:, :, i] = - np.log(self.stack4D[:, :, i, ith] / i0[i])

            # clean up the result
            nan_indices = np.where(np.isfinite(self.od) == False)
            if nan_indices:
                self.od[nan_indices] = 0

            self.od4D[:, :, :, ith] = self.od[:, :, :]

        self.od3d = self.od.copy()

        # Optical density matrix is rearranged into n_pixelsxn_ev
        self.od = np.reshape(self.od, (n_pixels, self.n_ev), order='F')

        return

    # ----------------------------------------------------------------------
    # Normalize the data: calculate optical density matrix D
    def calculate_optical_density_from_refimgs(self, files):

        n_pixels = self.n_cols * self.n_rows
        self.od = np.empty((self.n_cols, self.n_rows, self.n_ev))

        # zero out all negative values in the image stack
        negative_indices = np.where(self.absdata <= 0)
        if negative_indices:
            self.absdata[negative_indices] = 0.01

        # Load reference images
        refimgs = np.empty((self.n_cols, self.n_rows, self.n_ev))
        refimgs_ev = []
        for i in range(len(files)):
            ncols, nrows, iev, imgdata = file_xrm.read_xrm_fileinfo(files[i], readimgdata=True)
            refimgs[:, :, i] = np.reshape(imgdata, (ncols, nrows), order='F')
            refimgs_ev.append(iev)

        # Check if the energies are consecutive, if they are not sort the data
        consec = 0
        for i in range(len(refimgs_ev) - 1):
            if refimgs_ev[i] > refimgs_ev[i + 1]:
                consec = 1
                break
        if consec == 1:
            sortind = np.argsort(refimgs_ev)
            refimgs_ev = refimgs_ev[sortind]
            refimgs = refimgs[:, :, refimgs_ev]

        for i in range(self.n_ev):
            if self.ev[i] != refimgs_ev[i]:
                print('Error, wrong reference image energy')
                return

            self.od[:, :, i] = - np.log(self.absdata[:, :, i] / refimgs[:, :, i])

        # clean up the result
        nan_indices = np.where(np.isfinite(self.od) == False)
        if nan_indices:
            self.od[nan_indices] = 0

        self.od3d = self.od.copy()

        # Optical density matrix is rearranged into n_pixelsxn_ev
        self.od = np.reshape(self.od, (n_pixels, self.n_ev), order='F')

        self.evi0 = refimgs_ev
        self.i0data = np.ones((self.n_ev))
        self.i0_dwell = self.data_dwell

        return

    # ----------------------------------------------------------------------
    def scale_bar(self):
        self.x_start = np.min(self.x_dist)
        self.x_stop = np.max(self.x_dist)
        self.x_pxsize = np.round(np.abs(self.x_stop - self.x_start) / (self.n_cols - 1),
                                 5)  # um per px in y direction, "-1" because stop-start is 1 px shorter than n_rows

        self.y_start = np.min(self.y_dist)
        self.y_stop = np.max(self.y_dist)
        self.y_pxsize = np.round(np.abs(self.y_stop - self.y_start) / (self.n_rows - 1),
                                 5)  # um per px in y direction, "-1" because stop-start is 1 px shorter than n_rows

        if self.x_pxsize == self.y_pxsize:
            self.squarepx = True
        else:
            self.squarepx = False
        bar_microns = 0.2 * np.abs(self.x_stop - self.x_start)

        if bar_microns >= 10.:
            bar_microns = 10. * int(0.5 + 0.1 * int(0.5 + bar_microns))
            bar_string = str(int(0.01 + bar_microns)).strip()
        elif bar_microns >= 1.:
            bar_microns = float(int(0.5 + bar_microns))
            if bar_microns == 1.:
                bar_string = '1'
            else:
                bar_string = str(int(0.01 + bar_microns)).strip()
        else:
            bar_microns = np.maximum(0.1 * int(0.5 + 10 * bar_microns), 0.1)
            bar_string = str(bar_microns).strip()

        self.scale_bar_string = bar_string

        self.scale_bar_pixels_x = int(0.5 + float(self.n_cols) *
                                      float(bar_microns) / float(abs(self.x_stop - self.x_start)))

        self.scale_bar_pixels_y = int(0.01 * self.n_rows)

        if self.scale_bar_pixels_y < 2:
            self.scale_bar_pixels_y = 2

    # ----------------------------------------------------------------------
    def write_xas(self, filename, evdata, data):
        f = open(filename, 'w')
        print('*********************  X-ray Absorption Data  ********************', file=f)
        print('*', file=f)
        print('* Formula: ', file=f)
        print('* Common name: ', file=f)
        print('* Edge: ', file=f)
        print('* Acquisition mode: ', file=f)
        print('* Source and purity: ', file=f)
        print('* Comments: Stack list ROI ""', file=f)
        print('* Delta eV: ', file=f)
        print('* Min eV: ', file=f)
        print('* Max eV: ', file=f)
        print('* Y axis: ', file=f)
        print('* Contact person: ', file=f)
        print('* Write date: ', file=f)
        print('* Journal: ', file=f)
        print('* Authors: ', file=f)
        print('* Title: ', file=f)
        print('* Volume: ', file=f)
        print('* Issue number: ', file=f)
        print('* Year: ', file=f)
        print('* Pages: ', file=f)
        print('* Booktitle: ', file=f)
        print('* Editors: ', file=f)
        print('* Publisher: ', file=f)
        print('* Address: ', file=f)
        print('*--------------------------------------------------------------', file=f)
        for ie in range(self.n_ev):
            print('\t {0:06.2f}, {1:06f}'.format(evdata[ie], data[ie]), file=f)

        f.close()

        return

    # ----------------------------------------------------------------------
    def write_csv(self, filename, evdata, data, cname=''):
        f = open(filename, 'w')
        print('*********************  X-ray Absorption Data  ********************', file=f)
        print('*', file=f)
        print('* Formula: ', file=f)
        print('* Common name: {0}'.format(cname), file=f)
        print('* Edge: ', file=f)
        print('* Acquisition mode: ', file=f)
        print('* Source and purity: ', file=f)
        print('* Comments: Stack list ROI ""', file=f)
        print('* Delta eV: ', file=f)
        print('* Min eV: ', file=f)
        print('* Max eV: ', file=f)
        print('* Y axis: ', file=f)
        print('* Contact person: ', file=f)
        print('* Write date: ', file=f)
        print('* Journal: ', file=f)
        print('* Authors: ', file=f)
        print('* Title: ', file=f)
        print('* Volume: ', file=f)
        print('* Issue number: ', file=f)
        print('* Year: ', file=f)
        print('* Pages: ', file=f)
        print('* Booktitle: ', file=f)
        print('* Editors: ', file=f)
        print('* Publisher: ', file=f)
        print('* Address: ', file=f)
        print('*--------------------------------------------------------------', file=f)
        for ie in range(self.n_ev):
            print('{0:06.2f}, {1:06f}'.format(evdata[ie], data[ie]), file=f)

        f.close()

        return

    # ----------------------------------------------------------------------
    # Read x-ray absorption spectrum
    def read_xas(self, filename):

        spectrum_common_name = ' '

        f = open(str(filename), 'rU')

        elist = []
        ilist = []

        for line in f:
            if line.startswith('*'):
                if 'Common name' in line:
                    spectrum_common_name = line.split(':')[-1].strip()

            else:
                e, i = [float(x) for x in line.split()]
                elist.append(e)
                ilist.append(i)

        spectrum_evdata = np.array(elist)
        spectrum_data = np.array(ilist)

        f.close()

        if spectrum_evdata[-1] < spectrum_evdata[0]:
            spectrum_evdata = spectrum_evdata[::-1]
            spectrum_data = spectrum_data[::-1]

        if spectrum_common_name == ' ':
            spectrum_common_name = os.path.splitext(os.path.basename(str(filename)))[0]

        return spectrum_evdata, spectrum_data, spectrum_common_name

    # ----------------------------------------------------------------------
    # Read x-ray absorption spectrum
    def read_txt(self, filename):

        spectrum_common_name = os.path.splitext(os.path.basename(str(filename)))[0]

        f = open(str(filename), 'rU')

        elist = []
        ilist = []

        for line in f:
            if line.startswith('%'):
                pass
            else:
                e, i = [float(x) for x in line.split()]
                elist.append(e)
                ilist.append(i)

        spectrum_evdata = np.array(elist)
        spectrum_data = np.array(ilist)

        f.close()

        if spectrum_evdata[-1] < spectrum_evdata[0]:
            spectrum_evdata = spectrum_evdata[::-1]
            spectrum_data = spectrum_data[::-1]

        return spectrum_evdata, spectrum_data, spectrum_common_name

    # ----------------------------------------------------------------------
    # Read x-ray absorption spectrum
    def read_csv(self, filename):

        spectrum_common_name = ' '

        f = open(str(filename), 'rU')

        elist = []
        ilist = []

        # Check the first character of the line and skip if not a number
        allowedchars = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '-', '.']

        for line in f:
            if line.startswith('*'):
                if 'Common name' in line:
                    spectrum_common_name = line.split(':')[-1].strip()
            elif line[0] not in allowedchars:
                continue
            else:
                e, i = [float(x) for x in line.split(',')]
                elist.append(e)
                ilist.append(i)

        spectrum_evdata = np.array(elist)
        spectrum_data = np.array(ilist)

        f.close()

        if spectrum_evdata[-1] < spectrum_evdata[0]:
            spectrum_evdata = spectrum_evdata[::-1]
            spectrum_data = spectrum_data[::-1]

        if spectrum_common_name == ' ':
            spectrum_common_name = os.path.splitext(os.path.basename(str(filename)))[0]

        return spectrum_evdata, spectrum_data, spectrum_common_name

    # ----------------------------------------------------------------------
    # Register images using Fourier Shift Theorem
    # EdgeEnhancement: 0 = no edge enhacement; 1 = sobel; 2 = prewitt
    def register_images(self, ref_image, image2, have_ref_img_fft=False, edge_enhancement=0):

        if have_ref_img_fft == False:
            if edge_enhancement == 1:
                self.ref_fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(scipy.ndimage.filters.sobel(ref_image))))
            elif edge_enhancement == 2:
                self.ref_fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(scipy.ndimage.filters.prewitt(ref_image))))
            else:
                self.ref_fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(ref_image)))

        if edge_enhancement == 1:
            img2_fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(scipy.ndimage.filters.sobel(image2))))
        if edge_enhancement == 2:
            img2_fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(scipy.ndimage.filters.prewitt(image2))))
        else:
            img2_fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(image2)))

        fr = (self.ref_fft * img2_fft.conjugate()) / (np.abs(self.ref_fft) * np.abs(img2_fft))
        fr = np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(fr)))
        fr = np.abs(fr)

        shape = ref_image.shape

        xc, yc = np.unravel_index(np.argmax(fr), shape)

        # Limit the search to 1 pixel border
        if xc == 0:
            xc = 1
        if xc == shape[0] - 1:
            xc = shape[0] - 2

        if yc == 0:
            yc = 1
        if yc == shape[1] - 1:
            yc = shape[1] - 2

        # Use peak fit to find the shifts
        xpts = [xc - 1, xc, xc + 1]
        ypts = fr[xpts, yc]
        xf, fit = self.peak_fit(xpts, ypts)

        xpts = [yc - 1, yc, yc + 1]
        ypts = fr[xc, xpts]
        yf, fit = self.peak_fit(xpts, ypts)

        xshift = xf - np.float(shape[0]) / 2.0
        yshift = yf - np.float(shape[1]) / 2.0

        return xshift, yshift, fr

    # ----------------------------------------------------------------------
    # Apply image registration
    def apply_image_registration(self, image, xshift, yshift):

        shape = image.shape
        nx = shape[0]
        ny = shape[1]

        outofboundariesval = np.sum(image) / float(nx * ny)
        shifted_img = scipy.ndimage.interpolation.shift(image, [xshift, yshift],
                                                        mode='constant',
                                                        cval=outofboundariesval)

        return shifted_img

    # ----------------------------------------------------------------------
    # Apply image registration
    def crop_registed_images(self, images, min_xshift, max_xshift, min_yshift, max_yshift):

        # if the image is moved to the right (positive) we need to crop the left side
        xleft = int(np.ceil(max_xshift))
        if xleft < 0:
            xleft = 0
        # if the image is moved to the left (negative) we need to crop the right side
        xright = int(np.floor(self.n_cols + min_xshift))
        if xright > (self.n_cols):
            xright = int(self.n_cols)

        ybottom = int(np.ceil(max_yshift))
        if ybottom < 0:
            ybottom = 0
        ytop = int(np.floor(self.n_rows + min_yshift))
        if ytop > (self.n_rows):
            ytop = int(self.n_rows)

        if self.stack4D is not None:
            cropped_stack = images[xleft:xright, ybottom:ytop, :, :]
        else:
            cropped_stack = images[xleft:xright, ybottom:ytop, :]

        return cropped_stack, xleft, xright, ybottom, ytop

    # ----------------------------------------------------------------------
    # Quadratic peak fit: Fits the 3 data pairs to y=a+bx+cx^2, returning fit=[a,b,c]'
    #  and xpeak at position of inflection'
    def peak_fit(self, x, y):

        y1y0 = y[1] - y[0]
        y2y0 = y[2] - y[0]
        x1x0 = np.float(x[1] - x[0])
        x2x0 = np.float(x[2] - x[0])
        x1x0sq = np.float(x[1] * x[1] - x[0] * x[0])
        x2x0sq = np.float(x[2] * x[2] - x[0] * x[0])

        c_num = y2y0 * x1x0 - y1y0 * x2x0
        c_denom = x2x0sq * x1x0 - x1x0sq * x2x0

        if c_denom == 0:
            print('Divide by zero error')
            return

        c = c_num / np.float(c_denom)
        if x1x0 == 0:
            print('Divide by zero error')
            return

        b = (y1y0 - c * x1x0sq) / np.float(x1x0)
        a = y[0] - b * x[0] - c * x[0] * x[0]

        fit = [a, b, c]
        if c == 0:
            xpeak = 0.
            print('Cannot find xpeak')
            return
        else:
            # Constrain the fit to be within these three points.
            xpeak = -b / (2.0 * c)
            if xpeak < x[0]:
                xpeak = np.float(x[0])
            if xpeak > x[2]:
                xpeak = np.float(x[2])

        return xpeak, fit

    # -----------------------------------------------------------------------------
    # Despike image using Enhanced Lee Filter
    def despike(self, image, leefilt_percent=50.0):

        fimg = self.lee_filter(image)

        leefilt_max = np.amax(fimg)
        threshold = (1. + 0.01 * leefilt_percent) * leefilt_max

        datadim = np.int32(image.shape)

        ncols = datadim[0].copy()
        nrows = datadim[1].copy()

        spikes = np.where(image > threshold)
        n_spikes = fimg[spikes].shape[0]

        result_img = image.copy()

        if n_spikes > 0:

            xsp = spikes[0]
            ysp = spikes[1]
            for i in range(n_spikes):
                ix = xsp[i]
                iy = ysp[i]
                print(ix, iy)
                if ix == 0:
                    ix1 = 1
                    ix2 = 2
                elif ix == (ncols - 1):
                    ix1 = ncols - 2
                    ix2 = ncols - 3
                else:
                    ix1 = ix - 1
                    ix2 = ix + 1

                if iy == 0:
                    iy1 = 1
                    iy2 = 2
                elif iy == (nrows - 1):
                    iy1 = nrows - 2
                    iy2 = nrows - 3
                else:
                    iy1 = iy - 1
                    iy2 = iy + 1

                print(result_img[ix, iy])
                result_img[ix, iy] = 0.25 * (image[ix1, iy] + image[ix2, iy] +
                                             image[ix, iy1] + image[ix, iy2])
                print(result_img[ix, iy])

        return result_img

    # -----------------------------------------------------------------------------
    # Lee filter
    def lee_filter(self, image):

        nbox = 5  # The size of the filter box is 2N+1.  The default value is 5.
        sig = 5.0  # Estimate of the standard deviation.  The default is 5.

        delta = int((nbox - 1) / 2)  # width of window

        datadim = np.int32(image.shape)

        n_cols = datadim[0].copy()
        n_rows = datadim[1].copy()

        Imean = np.zeros((n_cols, n_rows))
        scipy.ndimage.filters.uniform_filter(image, size=nbox, output=Imean)

        Imean2 = Imean ** 2

        # variance
        z = np.empty((n_cols, n_rows))

        for l in range(delta, n_cols - delta):
            for s in range(delta, n_rows - delta):
                z[l, s] = np.sum((image[l - delta:l + delta, s - delta:s + delta] - Imean[l, s]) ** 2)

        z = z / float(nbox ** 2 - 1.0)

        z = (z + Imean2) / float(1.0 + sig ** 2) - Imean2

        ind = np.where(z < 0)
        n_ind = z[ind].shape[0]
        if n_ind > 0:
            z[ind] = 0

        lf_image = Imean + (image - Imean) * (z / (Imean2 * sig ** 2 + z))

        return lf_image
