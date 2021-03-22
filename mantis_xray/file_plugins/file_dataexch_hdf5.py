# -*- coding: utf-8 -*-
#
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
#
#   Copyright (C) 2011 Mirna Lerotic - 2nd Look, Nicholas Schwarz, Benjamin Watts
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

# File format details can be found on
# https://confluence.aps.anl.gov/display/NX/Data+Exchange+Basics


from __future__ import division
from __future__ import print_function

import os
import numpy as np
import h5py

from mantis_xray import data_struct
verbose = 0

title = 'Exchange'
extension = ['*.hdf','*.hdf5']
read_types = ['spectrum','image','stack']
write_types = ['spectrum','image','stack','results']


def identify(filename):
    try:
        # Open HDF5 file
         aps_format = False
         # Open HDF5 file
         f = h5py.File(filename, 'r')
         #Check is exchange is in the file
         if 'exchange' in f:
             aps_format = True
         f.close()
         return aps_format
    except:
        return False


def read_old( filepath, data_stk):
    data_stk.read_h5(filepath)
    return

def GetFileStructure(FileName):
    return None


#----------------------------------------------------------------------
def read(filename, stack_object, selection=None, *args, **kwargs):
    """ToDo: .attrs returns a list in py2 and a set in py3! .attrs and .keys() need to be wrapped in list()"""
    have4d = 0

    #new_stack = data_stack.data(data_struct)
    # Open HDF5 file
    f = h5py.File(filename, 'r')


    # Read basic definitions
    ds = f['implements']
    data_struct.implements = ds[...]
    ds = f['version']
    data_struct.version = ds[...]

    stack_object.have_dimscale = 0

    #Information HDF5 group
    if 'information' in f:
        informationGrp = f['information']
        if 'title' in informationGrp:
            title = informationGrp['title']
            data_struct.information.title = title[...]
        if 'comment' in informationGrp:
            com = informationGrp['comment']
            data_struct.information.comment = com[...]
        if 'file_creation_datetime' in informationGrp:
            fcdt = informationGrp['file_creation_datetime']
            data_struct.information.file_creation_datetime = fcdt[...]

        # /ids
        if 'ids' in informationGrp:
            idsGrp = informationGrp['ids']
            if 'proposal' in idsGrp:
                prop = idsGrp['proposal']
                data_struct.information.ids.proposal = prop[...]
            if 'acitivity' in idsGrp:
                act = idsGrp['activity']
                data_struct.information.ids.activity = act[...]
            if 'esaf' in idsGrp:
                esaf = idsGrp['esaf']
            data_struct.information.ids.esaf = esaf[...]

        # /experimenter
        if 'experimenter' in informationGrp:
            experimenterGrp = informationGrp['experimenter']
            if 'name' in experimenterGrp:
                name = experimenterGrp['name']
                data_struct.information.experimenter.name = name[...]
            if 'role' in experimenterGrp:
                role = experimenterGrp['role']
                data_struct.information.experimenter.role = role[...]
            if 'affiliation' in experimenterGrp:
                affiliation = experimenterGrp['affiliation']
                data_struct.information.experimenter.affiliation = affiliation[...]
            if 'address' in experimenterGrp:
                address = experimenterGrp['address']
                data_struct.information.experimenter.address = address[...]
            if 'phone' in experimenterGrp:
                phone = experimenterGrp['phone']
                data_struct.information.experimenter.phone = phone[...]
            if 'email' in experimenterGrp:
                email = experimenterGrp['email']
                data_struct.information.experimenter.email = email[...]
            if 'facility_user_id' in experimenterGrp:
                facility_user_id = experimenterGrp['facility_user_id']
                data_struct.information.experimenter.facility_user_id = facility_user_id[...]

        # /sample
        if 'sample' in informationGrp:
            sampleGrp = informationGrp['sample']
            if 'name' in sampleGrp:
                name = sampleGrp['name']
                data_struct.information.sample.name = name[...]
            if 'description' in sampleGrp:
                description = sampleGrp['description']
                data_struct.information.sample.description = description[...]
            if 'preparation_datetime' in sampleGrp:
                preparation_datetime = sampleGrp['preparation_datetime']
                data_struct.information.sample.preparation_datetime = preparation_datetime[...]
            if 'chemical_formula' in sampleGrp:
                cf = sampleGrp['chemical_formula']
                data_struct.information.sample.chemical_formula = cf[...]
            if 'environment' in sampleGrp:
                env = sampleGrp['environment']
                data_struct.information.sample.environment = env[...]
            if 'temperature' in sampleGrp:
                temp = sampleGrp['temperature']
                data_struct.information.sample.temperature = temp[...]
                data_struct.information.sample.temperature_units = temp.attrs['units']
            if 'pressure' in sampleGrp:
                press = sampleGrp['pressure']
                data_struct.information.sample.pressure = press[...]
                data_struct.information.sample.pressure_units = press.attrs['units']

        # /objective
        if 'objective' in informationGrp:
            objectiveGrp = informationGrp['objective']
            if 'manufacturer' in objectiveGrp:
                man = objectiveGrp['manufacturer']
                data_struct.information.objective.manufacturer = man[...]
            if 'model' in objectiveGrp:
                model = objectiveGrp['model']
                data_struct.information.objective.model = model[...]
            if 'comment' in objectiveGrp:
                com = objectiveGrp['comment']
                data_struct.information.objective.comment = com[...]
            if 'magnification' in objectiveGrp:
                mag = objectiveGrp['magnification']
                data_struct.information.objective.magnification = mag[...]

        # /scintillator
        if 'scintillator' in informationGrp:
            scintGrp = informationGrp['scintillator']
            if 'name' in scintGrp:
                name = scintGrp['name']
                data_struct.information.scintillator.name = name[...]
            if 'type' in scintGrp:
                type = scintGrp['type']
                data_struct.information.scintillator.type = type[...]
            if 'comment' in scintGrp:
                com = scintGrp['comment']
                data_struct.information.scintillator.comment = com[...]
            if 'scintillating_thickness' in scintGrp:
                sct = scintGrp['scintillating_thickness']
                data_struct.information.scintillator.scintillating_thickness = sct[...]
                data_struct.information.scintillator.scintillating_thickness_units = sct.attrs['units']
            if 'substrate_thickness' in scintGrp:
                subt = scintGrp['substrate_thickness']
                data_struct.information.scintillator.substrate_thickness = subt[...]
                data_struct.information.scintillator.substrate_thickness_units = subt.attrs['units']


        # /facility
        if 'facility' in informationGrp:
            facilityGrp = informationGrp['facility']
            if 'name' in facilityGrp:
                name = facilityGrp['name']
                data_struct.information.facility.name = name[...]
            if 'beamline' in facilityGrp:
                bl = facilityGrp['beamline']
                data_struct.information.facility.beamline = bl[...]


        # /accelerator
        if 'accelerator' in informationGrp:
            acceleratorGrp = informationGrp['accelerator']
            if 'ring_current' in acceleratorGrp:
                rc = acceleratorGrp['ring_current']
                data_struct.information.accelerator.ring_current = rc[...]
                data_struct.information.accelerator.ring_current_units = rc.attrs['units']
            if 'primary_beam_energy' in acceleratorGrp:
                pbe = acceleratorGrp['primary_beam_energy']
                data_struct.information.accelerator.primary_beam_energy = pbe[...]
                data_struct.information.accelerator.primary_beam_energy_units = pbe.attrs['units']
            if 'monostripe' in acceleratorGrp:
                ms = acceleratorGrp['monostripe']
                data_struct.information.accelerator.monostripe = ms[...]


        # /detector
        if 'detector' in informationGrp:
            detectorGrp = informationGrp['detector']
            if 'manufacturer' in detectorGrp:
                manf = detectorGrp['manufacturer']
                data_struct.information.detector.manufacturer = manf[...]
            if 'model' in detectorGrp:
                model = detectorGrp['model']
                data_struct.information.detector.model = model[...]
            if 'serial_number' in detectorGrp:
                sn = detectorGrp['serial_number']
                data_struct.information.detector.serial_number = sn[...]
            if 'bit_depth' in detectorGrp:
                bd = detectorGrp['bit_depth']
                data_struct.information.detector.bit_depth = bd[...]
            if 'operating_temperature' in detectorGrp:
                ot = detectorGrp['operating_temperature']
                data_struct.information.detector.operating_temperature = ot[...]
                data_struct.information.detector.operating_temperature_units = ot.attrs['units']
            if 'exposure_time' in detectorGrp:
                et = detectorGrp['exposure_time']
                data_struct.information.detector.exposure_time = et[...]
                data_struct.information.detector.exposure_time_units = et.attrs['units']
            if 'frame_rate' in detectorGrp:
                fr = detectorGrp['frame_rate']
                data_struct.information.detector.frame_rate = fr[...]

            if 'pixel_size' in detectorGrp:
                psGrp = detectorGrp['pixel_size']
                hor = psGrp['horizontal']
                data_struct.information.detector.pixel_size.horizontal = hor[...]
                data_struct.information.detector.pixel_size.horizontal_units = hor.attrs['units']
                ver = psGrp['vertical']
                data_struct.information.detector.pixel_size.vertical = ver[...]
                data_struct.information.detector.pixel_size.vertical_units = ver.attrs['units']

            if 'dimensions' in detectorGrp:
                dimGrp = detectorGrp['dimenstions']
                h = dimGrp['horizontal']
                data_struct.information.detector.dimensions.horizontal = h[...]
                v = dimGrp['vertical']
                data_struct.information.detector.dimensions.vertical = v[...]

            if 'binning' in detectorGrp:
                binningGrp = detectorGrp['binning']
                h = binningGrp['horizontal']
                data_struct.information.detector.binning.horizontal = h[...]
                v = binningGrp['vertical']
                data_struct.information.detector.binning.vertical = v[...]

            if 'axis_directions' in detectorGrp:
                axisdirGrp = detectorGrp['axis_directions']
                h = axisdirGrp['horizontal']
                data_struct.information.detector.axis_directions.horizontal = h[...]
                v = axisdirGrp['vertical']
                data_struct.information.detector.axis_directions.vertical = v[...]

            if 'roi' in detectorGrp:
                roiGrp = detectorGrp['roi']
                x1 = roiGrp['x1']
                data_struct.information.detector.roi.x1 = x1[...]
                y1 = roiGrp['y1']
                data_struct.information.detector.roi.y1 = y1[...]
                x2 = roiGrp['x2']
                data_struct.information.detector.roi.x2 = x2[...]
                y2 = roiGrp['x2']
                data_struct.information.detector.roi.y2 = y2[...]


    #Exchange HDF5 Group
    if 'exchange' in f:
        exchangeGrp = f['exchange']
        if 'title' in exchangeGrp:
            title = exchangeGrp['title']
            data_struct.exchange.title = title[...]
        if 'comment' in exchangeGrp:
            com = exchangeGrp['comment']
            data_struct.exchange.comment = com[...]
        if 'data_collection_datetime' in exchangeGrp:
            dcdt = exchangeGrp['data_collection_datetime']
            data_struct.exchange.data_collection_datetime = dcdt[...]

        if data_struct.version != '1.0':
            #load old version data
            # /exchange/detector
#                 detectorGrp = exchangeGrp['detector']
#                 dsdata = detectorGrp['data']
#                 data_struct.exchange.data = dsdata[...]

            #detectorGrp = exchangeGrp['detector']
            dsdata = exchangeGrp['data']
            data_struct.exchange.data = dsdata[...]

            if 'axes' in dsdata.attrs:
                data_struct.exchange.data_axes = dsdata.attrs['axes']

                if 'x' in exchangeGrp:
                    dsx = exchangeGrp['x']
                    data_struct.exchange.x = dsx[...]
                    data_struct.have_dimscale = 1
                if 'y' in exchangeGrp:
                    dsy = exchangeGrp['y']
                    data_struct.exchange.y = dsy[...]
#
#                 if 'x' in detectorGrp:
#                     dsx = detectorGrp['x']
#                     data_struct.exchange.x = dsx[...]
#                     stack_object.have_dimscale = 1
#
#                 if 'y' in detectorGrp:
#                     dsy = detectorGrp['y']
#                     data_struct.exchange.y = dsy[...]

            dseng = exchangeGrp['energy']
            data_struct.exchange.energy = dseng[...]
            data_struct.exchange.energy_units = dseng.attrs['units']

        else:
            # Version 1.0
            data = exchangeGrp['data']
            data_struct.exchange.data = data[...]
            if 'signal' in data.attrs:
                data_struct.exchange.data_signal = data.attrs['signal']
            if 'description' in data.attrs:
                data_struct.exchange.data_description = data.attrs['description']
            if 'units' in data.attrs:
                data_struct.exchange.data_units = data.attrs['units']
            if 'detector' in data.attrs:
                data_struct.exchange.data_detector = data.attrs['detector']

            if 'axes' in data.attrs:
                data_struct.exchange.data_axes = data.attrs['axes']
                axes_list = data_struct.exchange.data_axes.split(':')


                #Axes list can be arbitrary but for spectromicroscopy it is always x:y:z
                for i in axes_list:
                    ax = exchangeGrp[i]
                    if i == 'x':
                        data_struct.exchange.x = ax[...]
                        if 'units' in ax.attrs:
                            data_struct.exchange.x_units = ax.attrs['units']
                        stack_object.have_dimscale = 1
                    if i == 'y':
                        data_struct.exchange.y = ax[...]
                        if 'units' in ax.attrs:
                            data_struct.exchange.y_units = ax.attrs['units']
                    if i == 'z':
                        data_struct.exchange.z = ax[...]
                        if 'units' in ax.attrs:
                            data_struct.exchange.z_units = ax.attrs['units']

            energy = exchangeGrp['energy']
            data_struct.exchange.energy = energy[...]
            data_struct.exchange.energy_units= energy.attrs['units']


        if 'theta' in exchangeGrp:
            th = exchangeGrp['theta']
            data_struct.exchange.theta = th[...]
            data_struct.exchange.theta_units = th.attrs['units']
            have4d = 1

        if 'white_data' in exchangeGrp:
            wd = exchangeGrp['white_data']
            data_struct.exchange.white_data = wd[...]
            data_struct.exchange.white_data_units = wd.attrs['units']
        if 'dark_data' in exchangeGrp:
            dd = exchangeGrp['dark_data']
            data_struct.exchange.dark_data = dd[...]
            data_struct.exchange.dark_data_units = dd.attrs['units']
        if 'rotation' in exchangeGrp:
            rot = exchangeGrp['rotation']
            data_struct.exchange.rotation = rot[...]



    # Spectromicroscopy HDF5 group
    if 'spectromicroscopy' in f:
        spectromicroscopyGrp = f['spectromicroscopy']
        if 'positions' in spectromicroscopyGrp:
            pos = spectromicroscopyGrp['positions']
            data_struct.spectromicroscopy.positions = pos[...]
            data_struct.spectromicroscopy.positions_units = pos.attrs['units']
            data_struct.spectromicroscopy.positions_names = pos.attrs['names']

#         if 'optical_density' in spectromicroscopyGrp:
#             od = spectromicroscopyGrp['optical_density']
#             stack_object.data_struct.spectromicroscopy.optical_density = od[...]



        if 'normalization' in spectromicroscopyGrp:
            normGrp = spectromicroscopyGrp['normalization']
            ws = normGrp['white_spectrum']
            #Check if white_spectrum is an array
            wspectrum = ws[...]
            try:
                data_struct.spectromicroscopy.normalization.white_spectrum = ws[...]
                if 'units' in ws.attrs:
                    data_struct.spectromicroscopy.normalization.white_spectrum_units = ws.attrs['units']
                wse = normGrp['white_spectrum_energy']
                data_struct.spectromicroscopy.normalization.white_spectrum_energy = wse[...]
                if 'units' in ws.attrs:
                    data_struct.spectromicroscopy.normalization.white_spectrum_energy_units = wse.attrs['units']
            except:
                pass

    # Close
    f.close()

    if have4d == 0:
        stack_object.absdata = data_struct.exchange.data


    else:
        stack_object.stack4D = data_struct.exchange.data
        stack_object.theta = data_struct.exchange.theta
        stack_object.n_theta = len(stack_object.theta)
        stack_object.absdata = stack_object.stack4D[:,:,:,0]


    datadim = np.int32(stack_object.absdata.shape)


    stack_object.n_cols = datadim[0].copy()
    stack_object.n_rows =  datadim[1].copy()
    stack_object.ev = data_struct.exchange.energy

    stack_object.n_ev = np.int32(stack_object.ev.shape[0]).copy()

    npixels = stack_object.n_cols*stack_object.n_rows*stack_object.n_ev


    #Check if the energies are consecutive, if they are not sort the data
    consec = 0
    for i in range(stack_object.n_ev - 1):
        if stack_object.ev[i] > stack_object.ev[i+1]:
            consec = 1
            break
    if (consec == 1) and (have4d == 0):
        if verbose == 1: print("sort the energy data")
        sortind = np.argsort(stack_object.ev)
        stack_object.ev = stack_object.ev[sortind]
        stack_object.absdata = stack_object.absdata[:,:,sortind]
        #Save sorted energies:
        stack_object.data_struct.exchange.data = stack_object.absdata
        stack_object.data_struct.exchange.energy = stack_object.ev


    if stack_object.have_dimscale == 1:
        stack_object.x_dist = data_struct.exchange.x
        stack_object.y_dist = data_struct.exchange.y
    else:
        stack_object.x_dist = range(stack_object.n_cols)
        stack_object.y_dist = range(stack_object.n_rows)


    stack_object.i0data = data_struct.spectromicroscopy.normalization.white_spectrum
    stack_object.evi0 = data_struct.spectromicroscopy.normalization.white_spectrum_energy

    stack_object.data_dwell = np.ones((stack_object.n_ev))
    stack_object.i0_dwell = np.ones((stack_object.n_ev))


    if stack_object.data_struct.spectromicroscopy.normalization.white_spectrum is not None:
        stack_object.calculate_optical_density()
        stack_object.fill_h5_struct_normalization()


    if verbose == 1:
        print('filename ', filename)
        #print 'File creation date ', data_struct.file_creation_datetime
        print('Data array shape: ', stack_object.absdata.shape)
        print('n columns ', stack_object.n_cols)
        print('n_rows ', stack_object.n_rows)
        print('n_ev ', stack_object.n_ev)
        print('ev array ', stack_object.ev)
        print('x dist ', stack_object.x_dist)
        print('y_dist ', stack_object.y_dist)
        #print 'type ', type
        #print 'i0 data ', stack_object.i0data
        #print 'evi0 ', stack_object.evi0


    return

#----------------------------------------------------------------------
def write(filename, data_object, data_type):
    """Switchyard for writing different types of data."""
    if data_type in ['stack']:
        write_h5(filename, data_object.data_struct)

#----------------------------------------------------------------------
def write_h5(filename, data_struct):
    # Open HDF5 file
    f = h5py.File(filename, 'w')

    # This is an implementation of version 1.0
    data_struct.implements = 'information:exchange:spectromicroscopy'
    data_struct.version = '1.0'

    # HDF5 Root Group
    ds = f.create_dataset('implements', data = data_struct.implements)
    ds = f.create_dataset('version', data = data_struct.version)


    #Information HDF5 group
    if 'information' in dir(data_struct):
        informationGrp = f.create_group('information')
        if data_struct.information.title is not None:
            ds = informationGrp.create_dataset('title', data = data_struct.information.title)
        if data_struct.information.comment is not None:
            ds = informationGrp.create_dataset('comment', data = str(data_struct.information.comment))
        if data_struct.information.file_creation_datetime is not None:
            ds = informationGrp.create_dataset('file_creation_datetime', data = str(data_struct.information.file_creation_datetime))

        # /ids
        idsGrp = informationGrp.create_group('ids')
        have_ids = 0
        if data_struct.information.ids.proposal is not None:
            ds = idsGrp.create_dataset('proposal', data = data_struct.information.ids.proposal)
            have_ids = 1
        if data_struct.information.ids.activity is not None:
            ds = idsGrp.create_dataset('activity', data = data_struct.information.ids.activity)
            have_ids = 1
        if data_struct.information.ids.esaf is not None:
            ds = idsGrp.create_dataset('esaf', data = data_struct.information.ids.esaf)
            have_ids = 1
        if have_ids == 0:
            del informationGrp['ids']



        # /experimenter
        experimenterGrp = informationGrp.create_group('experimenter')
        have_exp = 0
        if data_struct.information.experimenter.name is not None:
            ds = experimenterGrp.create_dataset('name', data = str(data_struct.information.experimenter.name))
            have_exp = 1
        if data_struct.information.experimenter.role is not None:
            ds = experimenterGrp.create_dataset('role', data = data_struct.information.experimenter.role)
            have_exp = 1
        if data_struct.information.experimenter.affiliation is not None:
            ds = experimenterGrp.create_dataset('affiliation', data = data_struct.information.experimenter.affiliation)
            have_exp = 1
        if data_struct.information.experimenter.address is not None:
            ds = experimenterGrp.create_dataset('address', data = data_struct.information.experimenter.address)
            have_exp = 1
        if data_struct.information.experimenter.phone is not None:
            ds = experimenterGrp.create_dataset('phone', data = data_struct.information.experimenter.phone)
            have_exp = 1
        if data_struct.information.experimenter.email is not None:
            ds = experimenterGrp.create_dataset('email', data = data_struct.information.experimenter.email)
            have_exp = 1
        if data_struct.information.experimenter.facility_user_id is not None:
            ds = experimenterGrp.create_dataset('facility_user_id', data = data_struct.information.experimenter.facility_user_id)
            have_exp = 1
        if have_exp == 0:
            del informationGrp['experimenter']


        # /sample
        sampleGrp = informationGrp.create_group('sample')
        have_samp = 0
        if  data_struct.information.sample.name is not None:
            ds = sampleGrp.create_dataset('name', data = str(data_struct.information.sample.name))
            have_samp = 1
        if  data_struct.information.sample.description is not None:
            ds = sampleGrp.create_dataset('description', data = data_struct.information.sample.description)
            have_samp = 1
        if  data_struct.information.sample.preparation_datetime is not None:
            ds = sampleGrp.create_dataset('preparation_datetime', data = data_struct.information.sample.preparation_datetime)
            have_samp = 1
        if  data_struct.information.sample.chemical_formula is not None:
            ds = sampleGrp.create_dataset('chemical_formula', data = data_struct.information.chemical_formula)
            have_samp = 1
        if  data_struct.information.sample.environment is not None:
            ds = sampleGrp.create_dataset('environment', data = data_struct.information.sample.environment)
            have_samp = 1
        if  data_struct.information.sample.temperature is not None:
            ds = sampleGrp.create_dataset('temperature', data = data_struct.information.sample.temperature)
            ds.attrs['units'] = data_struct.information.sample.temperature_units
            have_samp = 1
        if  data_struct.information.sample.pressure is not None:
            ds = sampleGrp.create_dataset('pressure', data = data_struct.information.sample.pressure)
            ds.attrs['units'] = data_struct.information.sample.pressure_units
            have_samp = 1
        if have_samp == 0:
            del informationGrp['sample']


        # /objective
        objectiveGrp = informationGrp.create_group('objective')
        have_obj = 0
        if data_struct.information.objective.manufacturer is not None:
            ds = objectiveGrp.create_dataset('manufacturer', data = data_struct.information.objective.manufacturer)
            have_obj = 1
        if data_struct.information.objective.model is not None:
            ds = objectiveGrp.create_dataset('model', data = data_struct.information.objective.model)
            have_obj = 1
        if data_struct.information.objective.comment is not None:
            ds = objectiveGrp.create_dataset('comment', data = data_struct.information.objective.comment)
            have_obj = 1
        if data_struct.information.objective.magnification is not None:
            ds = objectiveGrp.create_dataset('magnification', data = data_struct.information.objective.magnification)
            have_obj = 1
        if have_obj == 0:
            del informationGrp['objective']



        # /scintillator

        scintGrp = informationGrp.create_group('scintillator')
        have_scint = 0
        if data_struct.information.scintillator.name is not None:
            ds = scintGrp.create_dataset('name', data = data_struct.information.scintillator.name)
            have_scint = 1
        if data_struct.information.scintillator.type is not None:
            ds = scintGrp.create_dataset('type', data = data_struct.information.scintillator.type)
            have_scint = 1
        if data_struct.information.scintillator.comment is not None:
            ds = scintGrp.create_dataset('comment', data = data_struct.information.scintillator.comment)
            have_scint = 1
        if data_struct.information.scintillator.scintillating_thickness is not None:
            ds = scintGrp.create_dataset('scintillating_thickness', data = data_struct.information.scintillator.scintillating_thickness)
            ds.attrs['units'] = data_struct.information.scintillator.scintillating_thickness_units
            have_scint = 1
        if data_struct.information.scintillator.substrate_thickness is not None:
            ds = scintGrp.create_dataset('substrate_thickness', data = data_struct.information.scintillator.substrate_thickness)
            ds.attrs['units'] = data_struct.information.scintillator.substrate_thickness_units
            have_scint = 1
        if have_scint == 0:
            del informationGrp['scintillator']


        # /facility
        facilityGrp = informationGrp.create_group('facility')
        have_fac = 0
        if data_struct.information.facility.name is not None:
            ds = facilityGrp.create_dataset('name', data = data_struct.information.facility.name)
            have_fac = 1
        if data_struct.information.facility.beamline is not None:
            ds = facilityGrp.create_dataset('beamline', data = data_struct.information.facility.beamline)
            have_fac = 1
        if have_fac == 0:
            del informationGrp['facility']


        # /accelerator
        acceleratorGrp = informationGrp.create_group('accelerator')
        have_ac = 0
        if data_struct.information.accelerator.ring_current is not None:
            ds = acceleratorGrp.create_dataset('ring_current', data = data_struct.information.accelerator.ring_current)
            ds.attrs['units'] = data_struct.information.accelerator.ring_current_units
            have_ac = 1
        if data_struct.information.accelerator.primary_beam_energy is not None:
            ds = acceleratorGrp.create_dataset('primary_beam_energy', data = data_struct.information.accelerator.primary_beam_energy)
            ds.attrs['units'] = data_struct.information.accelerator.primary_beam_energy_units
            have_ac = 1
        if data_struct.information.accelerator.monostripe is not None:
            ds = acceleratorGrp.create_dataset('monostripe', data = data_struct.information.accelerator.monostripe)
            have_ac = 1
        if have_ac == 0:
            del informationGrp['accelerator']




        # /detector
        detectorGrp = informationGrp.create_group('detector')
        have_detc = 0
        if data_struct.information.detector.manufacturer is not None:
            ds = detectorGrp.create_dataset('manufacturer', data = data_struct.information.detector.manufacturer)
            have_detc = 1
        if data_struct.information.detector.model is not None:
            ds = detectorGrp.create_dataset('model', data = data_struct.information.detector.model)
            have_detc = 1
        if data_struct.information.detector.serial_number is not None:
            ds = detectorGrp.create_dataset('serial_number', data = data_struct.information.detector.serial_number)
            have_detc = 1
        if data_struct.information.detector.bit_depth is not None:
            ds = detectorGrp.create_dataset('bit_depth', data = data_struct.information.detector.bit_depth)
            have_detc = 1
        if data_struct.information.detector.operating_temperature is not None:
            ds = detectorGrp.create_dataset('operating_temperature', data = data_struct.information.detector.operating_temperature)
            ds.attrs['units'] = data_struct.information.detector.operating_temperature_units
            have_detc = 1
        if data_struct.information.detector.exposure_time is not None:
            ds = detectorGrp.create_dataset('exposure_time', data = data_struct.information.detector.exposure_time)
            ds.attrs['units'] = data_struct.information.detector.exposure_time_units
            have_detc = 1
        if data_struct.information.detector.frame_rate is not None:
            ds = detectorGrp.create_dataset('frame_rate', data = data_struct.information.detector.frame_rate)
            have_detc = 1

        if data_struct.information.detector.pixel_size.horizontal is not None:
            psGrp = detectorGrp.create_group('pixel_size')
            ds = psGrp.create_dataset('horizontal', data = data_struct.information.detector.pixel_size.horizontal)
            ds.attrs['units'] = data_struct.information.detector.pixel_size.horizontal_units
            ds = psGrp.create_dataset('vertical', data = data_struct.information.detector.pixel_size.vertical)
            ds.attrs['units'] = data_struct.information.detector.pixel_size.vertical_units
            have_detc = 1

        if data_struct.information.detector.dimensions.horizontal is not None:
            dimGrp = detectorGrp.create_group('dimensions')
            ds = dimGrp.create_dataset('horizontal', data = data_struct.information.detector.dimensions.horizontal)
            ds = dimGrp.create_dataset('vertical', data = data_struct.information.detector.dimensions.vertical)
            have_detc = 1

        if data_struct.information.detector.binning.horizontal is not None:
            binningGrp = detectorGrp.create_group('binning')
            ds = binningGrp.create_dataset('horizontal', data = data_struct.information.detector.binning.horizontal)
            ds = binningGrp.create_dataset('vertical', data = data_struct.information.detector.binning.vertical)
            have_detc = 1

        if data_struct.information.detector.axis_directions.horizontal is not None:
            axisdirGrp = detectorGrp.create_group('axis_directions')
            ds = axisdirGrp.create_dataset('horizontal', data = data_struct.information.detector.axis_directions.horizontal)
            ds = axisdirGrp.create_dataset('vertical', data = data_struct.information.detector.axis_directions.vertical)
            have_detc = 1

        if data_struct.information.detector.roi.x1 is not None:
            roiGrp = detectorGrp.create_group('roi')
            ds = roiGrp.create_dataset('x1', data = data_struct.information.detector.roi.x1)
            ds = roiGrp.create_dataset('y1', data = data_struct.information.detector.roi.y1)
            ds = roiGrp.create_dataset('x2', data = data_struct.information.detector.roi.x2)
            ds = roiGrp.create_dataset('y2', data = data_struct.information.detector.roi.y2)
            have_detc = 1

        if have_detc == 0:
            del informationGrp['detector']


    # exchange definition

    # exchange HDF5 group
    # /exchange
    if 'exchange' in dir(data_struct):
        exchangeGrp = f.create_group("exchange")
        if data_struct.exchange.title is not None:
            ds = exchangeGrp.create_dataset('title', data = data_struct.exchange.title)
        if data_struct.exchange.comment is not None:
            ds = exchangeGrp.create_dataset('comment', data = data_struct.exchange.comment)
        if data_struct.exchange.data_collection_datetime is not None:
            ds = exchangeGrp.create_dataset('data_collection_datetime', data = data_struct.exchange.data_collection_datetime)

        # /exchange/data
        ds_data = exchangeGrp.create_dataset('data', data = data_struct.exchange.data)
        if data_struct.exchange.data_signal is not None:
            ds_data.attrs['signal'] = data_struct.exchange.data_signal
        if data_struct.exchange.data_description is not None:
            ds_data.attrs['description'] = data_struct.exchange.data_description
        if data_struct.exchange.data_units is not None:
            ds_data.attrs['units'] = data_struct.exchange.data_units
        if data_struct.exchange.data_detector is not None:
            ds_data.attrs['detector'] = data_struct.exchange.data_detector
        if data_struct.exchange.data_axes is not None:
            ds_data.attrs['axes'] = data_struct.exchange.data_axes

        if data_struct.exchange.x is not None:
            ds = exchangeGrp.create_dataset('x', data = data_struct.exchange.x)
        if data_struct.exchange.x_units is not None:
            ds.attrs['units'] = data_struct.exchange.x_units
        if data_struct.exchange.y is not None:
            ds = exchangeGrp.create_dataset('y', data = data_struct.exchange.y)
        if data_struct.exchange.y_units is not None:
            ds.attrs['units'] = data_struct.exchange.y_units
        if data_struct.exchange.z is not None:
            ds = exchangeGrp.create_dataset('z', data = data_struct.exchange.z)
        if data_struct.exchange.z_units is not None:
            ds.attrs['units'] = data_struct.exchange.z_units

        if data_struct.exchange.energy is not None:
            ds = exchangeGrp.create_dataset('energy', data = data_struct.exchange.energy)
        if data_struct.exchange.energy_units is not None:
            ds.attrs['units'] = data_struct.exchange.energy_units

        if data_struct.exchange.theta is not None:
            ds = exchangeGrp.create_dataset('theta', data = data_struct.exchange.theta)
            ds.attrs['units'] = data_struct.exchange.theta_units

        # /exchange/white_data
        if data_struct.exchange.white_data is not None:
            ds = exchangeGrp.create_dataset('white_data', data = data_struct.exchange.white_data)
            ds.attrs['units'] = data_struct.exchange.white_data_units

        # /exchange/dark_data
        if data_struct.exchange.dark_data is not None:
            ds = exchangeGrp.create_dataset('dark_data', data = data_struct.exchange.dark_data)
            ds.attrs['units'] = data_struct.exchange.dark_data_units

        if data_struct.exchange.rotation is not None:
            ds = exchangeGrp.create_dataset('rotation', data = data_struct.exchange.rotation)



    # Spectromicroscopy HDF5 group
    if 'spectromicroscopy' in dir(data_struct):
        spectromicroscopyGrp = f.create_group('spectromicroscopy')
        if data_struct.spectromicroscopy.positions is not None:
            ds = spectromicroscopyGrp.create_dataset('positions', data = data_struct.spectromicroscopy.positions)
        if data_struct.spectromicroscopy.positions_units is not None:
            ds.attrs['units'] = data_struct.spectromicroscopy.positions_units
        if data_struct.spectromicroscopy.positions_names is not None:
            ds.attrs['names'] = data_struct.spectromicroscopy.positions_names
        if data_struct.spectromicroscopy.xshifts is not None:
            ds = spectromicroscopyGrp.create_dataset('xshifts', data = data_struct.spectromicroscopy.xshifts)
        if data_struct.spectromicroscopy.yshifts is not None:
            ds = spectromicroscopyGrp.create_dataset('yshifts', data = data_struct.spectromicroscopy.yshifts)
        if data_struct.spectromicroscopy.optical_density is not None:
            ds = spectromicroscopyGrp.create_dataset('optical_density', data = data_struct.spectromicroscopy.optical_density)

        # /spectromicroscopy/normalization
        normalizationGrp = spectromicroscopyGrp.create_group('normalization')
        if data_struct.spectromicroscopy.normalization.white_spectrum is not None:
            ds = normalizationGrp.create_dataset('white_spectrum',
                                                data = data_struct.spectromicroscopy.normalization.white_spectrum)
            if data_struct.spectromicroscopy.normalization.white_spectrum_units is not None:
                ds.attrs['units'] = data_struct.spectromicroscopy.normalization.white_spectrum_units
            ds = normalizationGrp.create_dataset('white_spectrum_energy',
                                                    data = data_struct.spectromicroscopy.normalization.white_spectrum_energy)
            if data_struct.spectromicroscopy.normalization.white_spectrum_energy_units is not None:
                ds.attrs['units'] = data_struct.spectromicroscopy.normalization.white_spectrum_energy_units



        else:
            del spectromicroscopyGrp['normalization']


    # Close
    f.close()


#----------------------------------------------------------------------
def write_results_h5(filename, data_struct, anlz):

    test_file = 0
    #Check if file exists
    try:
        # Open HDF5 file
        f = h5py.File(filename, 'r')
        test_file = 1
        f.close()
    except:
        pass
    #Try to save new hdf5 file
    if test_file == 0:
        #try:
        write_h5(filename, data_struct)
        #except:
        #    print 'Error: Could not open nor create HDF5 file ', filename
        #    return -1


    # Open HDF5 file
    f = h5py.File(filename, 'r+')
    # Read basic definitions
    ds = f['implements']
    implements = ds[...]

    if 'Mantis' in f:
        mantisGrp = f['Mantis']
    else:
        del f['implements']
        implements = implements+':Mantis'
        ds = f.create_dataset('implements', data = implements)
        mantisGrp = f.create_group('Mantis')


    if anlz.pca_calculated == 1:
        if 'pca' in mantisGrp:
            del mantisGrp['pca']
            pcaGrp = mantisGrp.create_group('pca')
        else:
            pcaGrp = mantisGrp.create_group('pca')
        ds_data = pcaGrp.create_dataset('pca_images', data = anlz.pcaimages)
        ds_data = pcaGrp.create_dataset('pca_eigenvectors', data = anlz.eigenvecs)
        ds_data = pcaGrp.create_dataset('pca_eigenvalues', data = anlz.eigenvals)


    if anlz.clusters_calculated == 1:
        if 'cluster_analysis' in mantisGrp:
            del mantisGrp['cluster_analysis']
            caGrp = mantisGrp.create_group('cluster_analysis')
        else:
            caGrp = mantisGrp.create_group('cluster_analysis')
        ds_data = caGrp.create_dataset('cluster_indices', data = anlz.cluster_indices)
        ds_data = caGrp.create_dataset('cluster_spectra', data = anlz.clusterspectra)
        ds_data = caGrp.create_dataset('cluster_distances', data = anlz.cluster_distances)

    if anlz.tspectrum_loaded == 1:
        if 'spectral_maps' in mantisGrp:
            del mantisGrp['spectral_maps']
            caGrp = mantisGrp.create_group('spectral_maps')
        else:
            caGrp = mantisGrp.create_group('spectral_maps')
        ds_data = caGrp.create_dataset('raw_maps', data = anlz.target_svd_maps)
        ds_data = caGrp.create_dataset('fitted_maps', data = anlz.target_pcafit_maps)
        ds_data = caGrp.create_dataset('raw_spectra', data = anlz.target_spectra)
        ds_data = caGrp.create_dataset('fitted_spectra', data = anlz.target_pcafit_spectra)

    f.close()


    return 0





#----------------------------------------------------------------------
class h5:
    def __init__(self):
        pass


#----------------------------------------------------------------------
    def check_h5_format(self, filename):

        aps_format = False
        # Open HDF5 file
        f = h5py.File(filename, 'r')

        #Check is exchange is in the file
        #Exchange HDF5 Group
        if 'exchange' in f:
            aps_format = True

        f.close()

        return aps_format

#----------------------------------------------------------------------
    def read_h5(self, filename):#, data_struct):

        have4d = 0

        # Open HDF5 file
        f = h5py.File(filename, 'r')


        # Read basic definitions
        ds = f['implements']
        data_struct.implements = ds[...]
        ds = f['version']
        data_struct.version = ds[...]


        self.have_dimscale = 0


        #Information HDF5 group
        if 'information' in f:
            informationGrp = f['information']
            if 'title' in informationGrp:
                title = informationGrp['title']
                data_struct.information.title = title[...]
            if 'comment' in informationGrp:
                com = informationGrp['comment']
                data_struct.information.comment = com[...]
            if 'file_creation_datetime' in informationGrp:
                fcdt = informationGrp['file_creation_datetime']
                data_struct.information.file_creation_datetime = fcdt[...]

            # /ids
            if 'ids' in informationGrp:
                idsGrp = informationGrp['ids']
                if 'proposal' in idsGrp:
                    prop = idsGrp['proposal']
                    data_struct.information.ids.proposal = prop[...]
                if 'acitivity' in idsGrp:
                    act = idsGrp['activity']
                    data_struct.information.ids.activity = act[...]
                if 'esaf' in idsGrp:
                    esaf = idsGrp['esaf']
                data_struct.information.ids.esaf = esaf[...]

            # /experimenter
            if 'experimenter' in informationGrp:
                experimenterGrp = informationGrp['experimenter']
                if 'name' in experimenterGrp:
                    name = experimenterGrp['name']
                    data_struct.information.experimenter.name = name[...]
                if 'role' in experimenterGrp:
                    role = experimenterGrp['role']
                    data_struct.information.experimenter.role = role[...]
                if 'affiliation' in experimenterGrp:
                    affiliation = experimenterGrp['affiliation']
                    data_struct.information.experimenter.affiliation = affiliation[...]
                if 'address' in experimenterGrp:
                    address = experimenterGrp['address']
                    data_struct.information.experimenter.address = address[...]
                if 'phone' in experimenterGrp:
                    phone = experimenterGrp['phone']
                    data_struct.information.experimenter.phone = phone[...]
                if 'email' in experimenterGrp:
                    email = experimenterGrp['email']
                    data_struct.information.experimenter.email = email[...]
                if 'facility_user_id' in experimenterGrp:
                    facility_user_id = experimenterGrp['facility_user_id']
                    data_struct.information.experimenter.facility_user_id = facility_user_id[...]

            # /sample
            if 'sample' in informationGrp:
                sampleGrp = informationGrp['sample']
                if 'name' in sampleGrp:
                    name = sampleGrp['name']
                    data_struct.information.sample.name = name[...]
                if 'description' in sampleGrp:
                    description = sampleGrp['description']
                    data_struct.information.sample.description = description[...]
                if 'preparation_datetime' in sampleGrp:
                    preparation_datetime = sampleGrp['preparation_datetime']
                    data_struct.information.sample.preparation_datetime = preparation_datetime[...]
                if 'chemical_formula' in sampleGrp:
                    cf = sampleGrp['chemical_formula']
                    data_struct.information.sample.chemical_formula = cf[...]
                if 'environment' in sampleGrp:
                    env = sampleGrp['environment']
                    data_struct.information.sample.environment = env[...]
                if 'temperature' in sampleGrp:
                    temp = sampleGrp['temperature']
                    data_struct.information.sample.temperature = temp[...]
                    data_struct.information.sample.temperature_units = temp.attrs['units']
                if 'pressure' in sampleGrp:
                    press = sampleGrp['pressure']
                    data_struct.information.sample.pressure = press[...]
                    data_struct.information.sample.pressure_units = press.attrs['units']

            # /objective
            if 'objective' in informationGrp:
                objectiveGrp = informationGrp['objective']
                if 'manufacturer' in objectiveGrp:
                    man = objectiveGrp['manufacturer']
                    data_struct.information.objective.manufacturer = man[...]
                if 'model' in objectiveGrp:
                    model = objectiveGrp['model']
                    data_struct.information.objective.model = model[...]
                if 'comment' in objectiveGrp:
                    com = objectiveGrp['comment']
                    data_struct.information.objective.comment = com[...]
                if 'magnification' in objectiveGrp:
                    mag = objectiveGrp['magnification']
                    data_struct.information.objective.magnification = mag[...]

            # /scintillator
            if 'scintillator' in informationGrp:
                scintGrp = informationGrp['scintillator']
                if 'name' in scintGrp:
                    name = scintGrp['name']
                    data_struct.information.scintillator.name = name[...]
                if 'type' in scintGrp:
                    type = scintGrp['type']
                    data_struct.information.scintillator.type = type[...]
                if 'comment' in scintGrp:
                    com = scintGrp['comment']
                    data_struct.information.scintillator.comment = com[...]
                if 'scintillating_thickness' in scintGrp:
                    sct = scintGrp['scintillating_thickness']
                    data_struct.information.scintillator.scintillating_thickness = sct[...]
                    data_struct.information.scintillator.scintillating_thickness_units = sct.attrs['units']
                if 'substrate_thickness' in scintGrp:
                    subt = scintGrp['substrate_thickness']
                    data_struct.information.scintillator.substrate_thickness = subt[...]
                    data_struct.information.scintillator.substrate_thickness_units = subt.attrs['units']


            # /facility
            if 'facility' in informationGrp:
                facilityGrp = informationGrp['facility']
                if 'name' in facilityGrp:
                    name = facilityGrp['name']
                    data_struct.information.facility.name = name[...]
                if 'beamline' in facilityGrp:
                    bl = facilityGrp['beamline']
                    data_struct.information.facility.beamline = bl[...]


            # /accelerator
            if 'accelerator' in informationGrp:
                acceleratorGrp = informationGrp['accelerator']
                if 'ring_current' in acceleratorGrp:
                    rc = acceleratorGrp['ring_current']
                    data_struct.information.accelerator.ring_current = rc[...]
                    data_struct.information.accelerator.ring_current_units = rc.attrs['units']
                if 'primary_beam_energy' in acceleratorGrp:
                    pbe = acceleratorGrp['primary_beam_energy']
                    data_struct.information.accelerator.primary_beam_energy = pbe[...]
                    data_struct.information.accelerator.primary_beam_energy_units = pbe.attrs['units']
                if 'monostripe' in acceleratorGrp:
                    ms = acceleratorGrp['monostripe']
                    data_struct.information.accelerator.monostripe = ms[...]


            # /detector
            if 'detector' in informationGrp:
                detectorGrp = informationGrp['detector']
                if 'manufacturer' in detectorGrp:
                    manf = detectorGrp['manufacturer']
                    data_struct.information.detector.manufacturer = manf[...]
                if 'model' in detectorGrp:
                    model = detectorGrp['model']
                    data_struct.information.detector.model = model[...]
                if 'serial_number' in detectorGrp:
                    sn = detectorGrp['serial_number']
                    data_struct.information.detector.serial_number = sn[...]
                if 'bit_depth' in detectorGrp:
                    bd = detectorGrp['bit_depth']
                    data_struct.information.detector.bit_depth = bd[...]
                if 'operating_temperature' in detectorGrp:
                    ot = detectorGrp['operating_temperature']
                    data_struct.information.detector.operating_temperature = ot[...]
                    data_struct.information.detector.operating_temperature_units = ot.attrs['units']
                if 'exposure_time' in detectorGrp:
                    et = detectorGrp['exposure_time']
                    data_struct.information.detector.exposure_time = et[...]
                    data_struct.information.detector.exposure_time_units = et.attrs['units']
                if 'frame_rate' in detectorGrp:
                    fr = detectorGrp['frame_rate']
                    data_struct.information.detector.frame_rate = fr[...]

                if 'pixel_size' in detectorGrp:
                    psGrp = detectorGrp['pixel_size']
                    hor = psGrp['horizontal']
                    data_struct.information.detector.pixel_size.horizontal = hor[...]
                    data_struct.information.detector.pixel_size.horizontal_units = hor.attrs['units']
                    ver = psGrp['vertical']
                    data_struct.information.detector.pixel_size.vertical = ver[...]
                    data_struct.information.detector.pixel_size.vertical_units = ver.attrs['units']

                if 'dimensions' in detectorGrp:
                    dimGrp = detectorGrp['dimenstions']
                    h = dimGrp['horizontal']
                    data_struct.information.detector.dimensions.horizontal = h[...]
                    v = dimGrp['vertical']
                    data_struct.information.detector.dimensions.vertical = v[...]

                if 'binning' in detectorGrp:
                    binningGrp = detectorGrp['binning']
                    h = binningGrp['horizontal']
                    data_struct.information.detector.binning.horizontal = h[...]
                    v = binningGrp['vertical']
                    data_struct.information.detector.binning.vertical = v[...]

                if 'axis_directions' in detectorGrp:
                    axisdirGrp = detectorGrp['axis_directions']
                    h = axisdirGrp['horizontal']
                    data_struct.information.detector.axis_directions.horizontal = h[...]
                    v = axisdirGrp['vertical']
                    data_struct.information.detector.axis_directions.vertical = v[...]

                if 'roi' in detectorGrp:
                    roiGrp = detectorGrp['roi']
                    x1 = roiGrp['x1']
                    data_struct.information.detector.roi.x1 = x1[...]
                    y1 = roiGrp['y1']
                    data_struct.information.detector.roi.y1 = y1[...]
                    x2 = roiGrp['x2']
                    data_struct.information.detector.roi.x2 = x2[...]
                    y2 = roiGrp['x2']
                    data_struct.information.detector.roi.y2 = y2[...]


        #Exchange HDF5 Group
        if 'exchange' in f:
            exchangeGrp = f['exchange']
            if 'title' in exchangeGrp:
                title = exchangeGrp['title']
                data_struct.exchange.title = title[...]
            if 'comment' in exchangeGrp:
                com = exchangeGrp['comment']
                data_struct.exchange.comment = com[...]
            if 'data_collection_datetime' in exchangeGrp:
                dcdt = exchangeGrp['data_collection_datetime']
                data_struct.exchange.data_collection_datetime = dcdt[...]

            if data_struct.version != '1.0':
                #load old version data
                # /exchange/detector
#                 detectorGrp = exchangeGrp['detector']
#                 dsdata = detectorGrp['data']
#                 data_struct.exchange.data = dsdata[...]

                #detectorGrp = exchangeGrp['detector']
                dsdata = exchangeGrp['data']
                data_struct.exchange.data = dsdata[...]

                if 'axes' in dsdata.attrs:
                    data_struct.exchange.data_axes = dsdata.attrs['axes']


#
#                 if 'x' in detectorGrp:
#                     dsx = detectorGrp['x']
#                     data_struct.exchange.x = dsx[...]
#                     self.have_dimscale = 1
#
#                 if 'y' in detectorGrp:
#                     dsy = detectorGrp['y']
#                     data_struct.exchange.y = dsy[...]

                dseng = exchangeGrp['energy']
                data_struct.exchange.energy = dseng[...]
                #data_struct.exchange.energy_units = dseng.attrs['units']

            else:
                # Version 1.0
                data = exchangeGrp['data']
                data_struct.exchange.data = data[...]
                if 'signal' in data.attrs:
                    data_struct.exchange.data_signal = data.attrs['signal']
                if 'description' in data.attrs:
                    data_struct.exchange.data_description = data.attrs['description']
                if 'units' in data.attrs:
                    data_struct.exchange.data_units = data.attrs['units']
                if 'detector' in data.attrs:
                    data_struct.exchange.data_detector = data.attrs['detector']

                if 'axes' in data.attrs:
                    data_struct.exchange.data_axes = data.attrs['axes']
                    axes_list = data_struct.exchange.data_axes.split(':')

                    #Axes list can be arbitrary but for spectromicroscopy it is always x:y:z
                    for i in axes_list:
                        ax = exchangeGrp[i]
                        if i == 'x':
                            data_struct.exchange.x = ax[...]
                            if 'units' in ax.attrs:
                                data_struct.exchange.x_units = ax.attrs['units']
                            self.have_dimscale = 1
                        if i == 'y':
                            data_struct.exchange.y = ax[...]
                            if 'units' in ax.attrs:
                                data_struct.exchange.y_units = ax.attrs['units']
                        if i == 'z':
                            data_struct.exchange.z = ax[...]
                            if 'units' in ax.attrs:
                                data_struct.exchange.z_units = ax.attrs['units']

                energy = exchangeGrp['energy']
                data_struct.exchange.energy = energy[...]
                data_struct.exchange.energy_units= energy.attrs['units']

            if 'theta' in exchangeGrp:
                th = exchangeGrp['theta']
                data_struct.exchange.theta = th[...]
                data_struct.exchange.theta_units = th.attrs['units']
                have4d = 1

            if 'white_data' in exchangeGrp:
                wd = exchangeGrp['white_data']
                data_struct.exchange.white_data = wd[...]
                data_struct.exchange.white_data_units = wd.attrs['units']
            if 'dark_data' in exchangeGrp:
                dd = exchangeGrp['dark_data']
                data_struct.exchange.dark_data = dd[...]
                data_struct.exchange.dark_data_units = dd.attrs['units']
            if 'rotation' in exchangeGrp:
                rot = exchangeGrp['rotation']
                data_struct.exchange.rotation = rot[...]



        # Spectromicroscopy HDF5 group
        if 'spectromicroscopy' in f:
            spectromicroscopyGrp = f['spectromicroscopy']
            if 'positions' in spectromicroscopyGrp:
                pos = spectromicroscopyGrp['positions']
                data_struct.spectromicroscopy.positions = pos[...]
                data_struct.spectromicroscopy.positions_units = pos.attrs['units']
                data_struct.spectromicroscopy.positions_names = pos.attrs['names']

#             if 'optical_density' in spectromicroscopyGrp:
#                 od = spectromicroscopyGrp['optical_density']
#                 self.data_struct.spectromicroscopy.optical_density = od[...]



            if 'normalization' in spectromicroscopyGrp:
                normGrp = spectromicroscopyGrp['normalization']
                ws = normGrp['white_spectrum']
                #Check if white_spectrum is an array
                wspectrum = ws[...]
                try:
                    data_struct.spectromicroscopy.normalization.white_spectrum = ws[...]
                    if 'units' in ws.attrs:
                        data_struct.spectromicroscopy.normalization.white_spectrum_units = ws.attrs['units']
                    wse = normGrp['white_spectrum_energy']
                    data_struct.spectromicroscopy.normalization.white_spectrum_energy = wse[...]
                    if 'units' in ws.attrs:
                        data_struct.spectromicroscopy.normalization.white_spectrum_energy_units = wse.attrs['units']
                except:
                    pass

        # Close
        f.close()

        if have4d == 0:
            self.absdata = data_struct.exchange.data
        else:
            self.stack4D = data_struct.exchange.data
            self.theta = data_struct.exchange.theta
            self.n_theta = len(self.theta)
            self.absdata = self.stack4D[:,:,:,0]


        datadim = np.int32(self.absdata.shape)


        self.n_cols = datadim[0].copy()
        self.n_rows =  datadim[1].copy()
        self.ev = data_struct.exchange.energy

        self.n_ev = np.int32(self.ev.shape[0]).copy()

        npixels = self.n_cols*self.n_rows*self.n_ev


        #Check if the energies are consecutive, if they are not sort the data
        consec = 0
        for i in range(self.n_ev - 1):
            if self.ev[i] > self.ev[i+1]:
                consec = 1
                break
        if consec == 1:
            if verbose == 1: print("sort the energy data")
            sortind = np.argsort(self.ev)
            self.ev = self.ev[sortind]
            self.absdata = self.absdata[:,:,sortind]
            #Save sorted energies:
            self.data_struct.exchange.data = self.absdata
            self.data_struct.exchange.energy = self.ev


        if self.have_dimscale == 1:
            self.x_dist = data_struct.exchange.x
            self.y_dist = data_struct.exchange.y
        else:
            self.x_dist = range(self.n_cols)
            self.y_dist = range(self.n_rows)


        self.i0data = data_struct.spectromicroscopy.normalization.white_spectrum
        self.evi0 = data_struct.spectromicroscopy.normalization.white_spectrum_energy

        self.data_dwell = np.ones((self.n_ev))
        self.i0_dwell = np.ones((self.n_ev))



        if verbose == 1:
            print('filename ', filename)
            #print 'File creation date ', data_struct.file_creation_datetime
            print('Data array shape: ', self.absdata.shape)
            print('n columns ', self.n_cols)
            print('n_rows ', self.n_rows)
            print('n_ev ', self.n_ev)
            print('ev array ', self.ev)
            print('x dist ', self.x_dist)
            print('y_dist ', self.y_dist)
            #print 'type ', type
            #print 'i0 data ', self.i0data
            #print 'evi0 ', self.evi0


        return

#----------------------------------------------------------------------
    def write_h5(self, filename, data_struct):

        # Open HDF5 file
        f = h5py.File(filename, 'w')

        # This is an implementation of version 1.0
        data_struct.implements = 'information:exchange:spectromicroscopy'
        data_struct.version = '1.0'

        # HDF5 Root Group
        ds = f.create_dataset('implements', data = data_struct.implements)
        ds = f.create_dataset('version', data = data_struct.version)


        #Information HDF5 group
        informationGrp = f.create_group('information')
        if data_struct.information.title is not None:
            ds = informationGrp.create_dataset('title', data = data_struct.information.title)
        if data_struct.information.comment is not None:
            ds = informationGrp.create_dataset('comment', data = str(data_struct.information.comment))
        if data_struct.information.file_creation_datetime is not None:
            ds = informationGrp.create_dataset('file_creation_datetime', data = str(data_struct.information.file_creation_datetime))

        # /ids
        idsGrp = informationGrp.create_group('ids')
        have_ids = 0
        if data_struct.information.ids.proposal is not None:
            ds = idsGrp.create_dataset('proposal', data = data_struct.information.ids.proposal)
            have_ids = 1
        if data_struct.information.ids.activity is not None:
            ds = idsGrp.create_dataset('activity', data = data_struct.information.ids.activity)
            have_ids = 1
        if data_struct.information.ids.esaf is not None:
            ds = idsGrp.create_dataset('esaf', data = data_struct.information.ids.esaf)
            have_ids = 1
        if have_ids == 0:
            del informationGrp['ids']



        # /experimenter
        experimenterGrp = informationGrp.create_group('experimenter')
        have_exp = 0
        if data_struct.information.experimenter.name is not None:
            ds = experimenterGrp.create_dataset('name', data = str(data_struct.information.experimenter.name))
            have_exp = 1
        if data_struct.information.experimenter.role is not None:
            ds = experimenterGrp.create_dataset('role', data = data_struct.information.experimenter.role)
            have_exp = 1
        if data_struct.information.experimenter.affiliation is not None:
            ds = experimenterGrp.create_dataset('affiliation', data = data_struct.information.experimenter.affiliation)
            have_exp = 1
        if data_struct.information.experimenter.address is not None:
            ds = experimenterGrp.create_dataset('address', data = data_struct.information.experimenter.address)
            have_exp = 1
        if data_struct.information.experimenter.phone is not None:
            ds = experimenterGrp.create_dataset('phone', data = data_struct.information.experimenter.phone)
            have_exp = 1
        if data_struct.information.experimenter.email is not None:
            ds = experimenterGrp.create_dataset('email', data = data_struct.information.experimenter.email)
            have_exp = 1
        if data_struct.information.experimenter.facility_user_id is not None:
            ds = experimenterGrp.create_dataset('facility_user_id', data = data_struct.information.experimenter.facility_user_id)
            have_exp = 1
        if have_exp == 0:
            del informationGrp['experimenter']


        # /sample
        sampleGrp = informationGrp.create_group('sample')
        have_samp = 0
        if  data_struct.information.sample.name is not None:
            ds = sampleGrp.create_dataset('name', data = str(data_struct.information.sample.name))
            have_samp = 1
        if  data_struct.information.sample.description is not None:
            ds = sampleGrp.create_dataset('description', data = data_struct.information.sample.description)
            have_samp = 1
        if  data_struct.information.sample.preparation_datetime is not None:
            ds = sampleGrp.create_dataset('preparation_datetime', data = data_struct.information.sample.preparation_datetime)
            have_samp = 1
        if  data_struct.information.sample.chemical_formula is not None:
            ds = sampleGrp.create_dataset('chemical_formula', data = data_struct.information.chemical_formula)
            have_samp = 1
        if  data_struct.information.sample.environment is not None:
            ds = sampleGrp.create_dataset('environment', data = data_struct.information.sample.environment)
            have_samp = 1
        if  data_struct.information.sample.temperature is not None:
            ds = sampleGrp.create_dataset('temperature', data = data_struct.information.sample.temperature)
            ds.attrs['units'] = data_struct.information.sample.temperature_units
            have_samp = 1
        if  data_struct.information.sample.pressure is not None:
            ds = sampleGrp.create_dataset('pressure', data = data_struct.information.sample.pressure)
            ds.attrs['units'] = data_struct.information.sample.pressure_units
            have_samp = 1
        if have_samp == 0:
            del informationGrp['sample']


        # /objective
        objectiveGrp = informationGrp.create_group('objective')
        have_obj = 0
        if data_struct.information.objective.manufacturer is not None:
            ds = objectiveGrp.create_dataset('manufacturer', data = data_struct.information.objective.manufacturer)
            have_obj = 1
        if data_struct.information.objective.model is not None:
            ds = objectiveGrp.create_dataset('model', data = data_struct.information.objective.model)
            have_obj = 1
        if data_struct.information.objective.comment is not None:
            ds = objectiveGrp.create_dataset('comment', data = data_struct.information.objective.comment)
            have_obj = 1
        if data_struct.information.objective.magnification is not None:
            ds = objectiveGrp.create_dataset('magnification', data = data_struct.information.objective.magnification)
            have_obj = 1
        if have_obj == 0:
            del informationGrp['objective']



        # /scintillator

        scintGrp = informationGrp.create_group('scintillator')
        have_scint = 0
        if data_struct.information.scintillator.name is not None:
            ds = scintGrp.create_dataset('name', data = data_struct.information.scintillator.name)
            have_scint = 1
        if data_struct.information.scintillator.type is not None:
            ds = scintGrp.create_dataset('type', data = data_struct.information.scintillator.type)
            have_scint = 1
        if data_struct.information.scintillator.comment is not None:
            ds = scintGrp.create_dataset('comment', data = data_struct.information.scintillator.comment)
            have_scint = 1
        if data_struct.information.scintillator.scintillating_thickness is not None:
            ds = scintGrp.create_dataset('scintillating_thickness', data = data_struct.information.scintillator.scintillating_thickness)
            ds.attrs['units'] = data_struct.information.scintillator.scintillating_thickness_units
            have_scint = 1
        if data_struct.information.scintillator.substrate_thickness is not None:
            ds = scintGrp.create_dataset('substrate_thickness', data = data_struct.information.scintillator.substrate_thickness)
            ds.attrs['units'] = data_struct.information.scintillator.substrate_thickness_units
            have_scint = 1
        if have_scint == 0:
            del informationGrp['scintillator']


        # /facility
        facilityGrp = informationGrp.create_group('facility')
        have_fac = 0
        if data_struct.information.facility.name is not None:
            ds = facilityGrp.create_dataset('name', data = data_struct.information.facility.name)
            have_fac = 1
        if data_struct.information.facility.beamline is not None:
            ds = facilityGrp.create_dataset('beamline', data = data_struct.information.facility.beamline)
            have_fac = 1
        if have_fac == 0:
            del informationGrp['facility']


        # /accelerator
        acceleratorGrp = informationGrp.create_group('accelerator')
        have_ac = 0
        if data_struct.information.accelerator.ring_current is not None:
            ds = acceleratorGrp.create_dataset('ring_current', data = data_struct.information.accelerator.ring_current)
            ds.attrs['units'] = data_struct.information.accelerator.ring_current_units
            have_ac = 1
        if data_struct.information.accelerator.primary_beam_energy is not None:
            ds = acceleratorGrp.create_dataset('primary_beam_energy', data = data_struct.information.accelerator.primary_beam_energy)
            ds.attrs['units'] = data_struct.information.accelerator.primary_beam_energy_units
            have_ac = 1
        if data_struct.information.accelerator.monostripe is not None:
            ds = acceleratorGrp.create_dataset('monostripe', data = data_struct.information.accelerator.monostripe)
            have_ac = 1
        if have_ac == 0:
            del informationGrp['accelerator']




        # /detector
        detectorGrp = informationGrp.create_group('detector')
        have_detc = 0
        if data_struct.information.detector.manufacturer is not None:
            ds = detectorGrp.create_dataset('manufacturer', data = data_struct.information.detector.manufacturer)
            have_detc = 1
        if data_struct.information.detector.model is not None:
            ds = detectorGrp.create_dataset('model', data = data_struct.information.detector.model)
            have_detc = 1
        if data_struct.information.detector.serial_number is not None:
            ds = detectorGrp.create_dataset('serial_number', data = data_struct.information.detector.serial_number)
            have_detc = 1
        if data_struct.information.detector.bit_depth is not None:
            ds = detectorGrp.create_dataset('bit_depth', data = data_struct.information.detector.bit_depth)
            have_detc = 1
        if data_struct.information.detector.operating_temperature is not None:
            ds = detectorGrp.create_dataset('operating_temperature', data = data_struct.information.detector.operating_temperature)
            ds.attrs['units'] = data_struct.information.detector.operating_temperature_units
            have_detc = 1
        if data_struct.information.detector.exposure_time is not None:
            ds = detectorGrp.create_dataset('exposure_time', data = data_struct.information.detector.exposure_time)
            ds.attrs['units'] = data_struct.information.detector.exposure_time_units
            have_detc = 1
        if data_struct.information.detector.frame_rate is not None:
            ds = detectorGrp.create_dataset('frame_rate', data = data_struct.information.detector.frame_rate)
            have_detc = 1

        if data_struct.information.detector.pixel_size.horizontal is not None:
            psGrp = detectorGrp.create_group('pixel_size')
            ds = psGrp.create_dataset('horizontal', data = data_struct.information.detector.pixel_size.horizontal)
            ds.attrs['units'] = data_struct.information.detector.pixel_size.horizontal_units
            ds = psGrp.create_dataset('vertical', data = data_struct.information.detector.pixel_size.vertical)
            ds.attrs['units'] = data_struct.information.detector.pixel_size.vertical_units
            have_detc = 1

        if data_struct.information.detector.dimensions.horizontal is not None:
            dimGrp = detectorGrp.create_group('dimensions')
            ds = dimGrp.create_dataset('horizontal', data = data_struct.information.detector.dimensions.horizontal)
            ds = dimGrp.create_dataset('vertical', data = data_struct.information.detector.dimensions.vertical)
            have_detc = 1

        if data_struct.information.detector.binning.horizontal is not None:
            binningGrp = detectorGrp.create_group('binning')
            ds = binningGrp.create_dataset('horizontal', data = data_struct.information.detector.binning.horizontal)
            ds = binningGrp.create_dataset('vertical', data = data_struct.information.detector.binning.vertical)
            have_detc = 1

        if data_struct.information.detector.axis_directions.horizontal is not None:
            axisdirGrp = detectorGrp.create_group('axis_directions')
            ds = axisdirGrp.create_dataset('horizontal', data = data_struct.information.detector.axis_directions.horizontal)
            ds = axisdirGrp.create_dataset('vertical', data = data_struct.information.detector.axis_directions.vertical)
            have_detc = 1

        if data_struct.information.detector.roi.x1 is not None:
            roiGrp = detectorGrp.create_group('roi')
            ds = roiGrp.create_dataset('x1', data = data_struct.information.detector.roi.x1)
            ds = roiGrp.create_dataset('y1', data = data_struct.information.detector.roi.y1)
            ds = roiGrp.create_dataset('x2', data = data_struct.information.detector.roi.x2)
            ds = roiGrp.create_dataset('y2', data = data_struct.information.detector.roi.y2)
            have_detc = 1

        if have_detc == 0:
            del informationGrp['detector']


        # exchange definition

        # exchange HDF5 group
        # /exchange
        exchangeGrp = f.create_group("exchange")
        if data_struct.exchange.title is not None:
            ds = exchangeGrp.create_dataset('title', data = data_struct.exchange.title)
        if data_struct.exchange.comment is not None:
            ds = exchangeGrp.create_dataset('comment', data = data_struct.exchange.comment)
        if data_struct.exchange.data_collection_datetime is not None:
            ds = exchangeGrp.create_dataset('data_collection_datetime', data = data_struct.exchange.data_collection_datetime)

        # /exchange/data
        ds_data = exchangeGrp.create_dataset('data', data = data_struct.exchange.data)
        if data_struct.exchange.data_signal is not None:
            ds_data.attrs['signal'] = data_struct.exchange.data_signal
        if data_struct.exchange.data_description is not None:
            ds_data.attrs['description'] = data_struct.exchange.data_description
        if data_struct.exchange.data_units is not None:
            ds_data.attrs['units'] = data_struct.exchange.data_units
        if data_struct.exchange.data_detector is not None:
            ds_data.attrs['detector'] = data_struct.exchange.data_detector
        if data_struct.exchange.data_axes is not None:
            ds_data.attrs['axes'] = data_struct.exchange.data_axes

        if data_struct.exchange.x is not None:
            ds = exchangeGrp.create_dataset('x', data = data_struct.exchange.x)
        if data_struct.exchange.x_units is not None:
            ds.attrs['units'] = data_struct.exchange.x_units
        if data_struct.exchange.y is not None:
            ds = exchangeGrp.create_dataset('y', data = data_struct.exchange.y)
        if data_struct.exchange.y_units is not None:
            ds.attrs['units'] = data_struct.exchange.y_units
        if data_struct.exchange.z is not None:
            ds = exchangeGrp.create_dataset('z', data = data_struct.exchange.z)
        if data_struct.exchange.z_units is not None:
            ds.attrs['units'] = data_struct.exchange.z_units

        if data_struct.exchange.energy is not None:
            ds = exchangeGrp.create_dataset('energy', data = data_struct.exchange.energy)
            ds.attrs['units'] = data_struct.exchange.energy_units

        if data_struct.exchange.theta is not None:
            ds = exchangeGrp.create_dataset('theta', data = data_struct.exchange.theta)
            ds.attrs['units'] = data_struct.exchange.theta_units

        # /exchange/white_data
        if data_struct.exchange.white_data is not None:
            ds = exchangeGrp.create_dataset('white_data', data = data_struct.exchange.white_data)
            ds.attrs['units'] = data_struct.exchange.white_data_units

        # /exchange/dark_data
        if data_struct.exchange.dark_data is not None:
            ds = exchangeGrp.create_dataset('dark_data', data = data_struct.exchange.dark_data)
            ds.attrs['units'] = data_struct.exchange.dark_data_units

        if data_struct.exchange.rotation is not None:
            ds = exchangeGrp.create_dataset('rotation', data = data_struct.exchange.rotation)



        # Spectromicroscopy HDF5 group
        spectromicroscopyGrp = f.create_group('spectromicroscopy')
        if data_struct.spectromicroscopy.positions is not None:
            ds = spectromicroscopyGrp.create_dataset('positions', data = data_struct.spectromicroscopy.positions)
        if data_struct.spectromicroscopy.positions_units is not None:
            ds.attrs['units'] = data_struct.spectromicroscopy.positions_units
        if data_struct.spectromicroscopy.positions_names is not None:
            ds.attrs['names'] = data_struct.spectromicroscopy.positions_names
        if data_struct.spectromicroscopy.xshifts is not None:
            ds = spectromicroscopyGrp.create_dataset('xshifts', data = data_struct.spectromicroscopy.xshifts)
        if data_struct.spectromicroscopy.yshifts is not None:
            ds = spectromicroscopyGrp.create_dataset('yshifts', data = data_struct.spectromicroscopy.yshifts)
        if self.data_struct.spectromicroscopy.optical_density is not None:
            ds = spectromicroscopyGrp.create_dataset('optical_density', data = self.data_struct.spectromicroscopy.optical_density)

        # /spectromicroscopy/normalization
        normalizationGrp = spectromicroscopyGrp.create_group('normalization')
        if data_struct.spectromicroscopy.normalization.white_spectrum is not None:
            ds = normalizationGrp.create_dataset('white_spectrum',
                                             data = data_struct.spectromicroscopy.normalization.white_spectrum)
            if data_struct.spectromicroscopy.normalization.white_spectrum_units is not None:
                ds.attrs['units'] = data_struct.spectromicroscopy.normalization.white_spectrum_units
            ds = normalizationGrp.create_dataset('white_spectrum_energy',
                                                 data = data_struct.spectromicroscopy.normalization.white_spectrum_energy)
            if data_struct.spectromicroscopy.normalization.white_spectrum_energy_units is not None:
                ds.attrs['units'] = data_struct.spectromicroscopy.normalization.white_spectrum_energy_units



        else:
            del spectromicroscopyGrp['normalization']


        # Close
        f.close()


#----------------------------------------------------------------------
    def write_results_h5(self, filename, data_struct, anlz):

        test_file = 0
        #Check if file exists
        try:
            # Open HDF5 file
            f = h5py.File(filename, 'r')
            test_file = 1
            f.close()
        except:
            pass
        #Try to save new hdf5 file
        if test_file == 0:
            #try:
            self.write_h5(filename, data_struct)
            #except:
            #    print 'Error: Could not open nor create HDF5 file ', filename
            #    return -1


        # Open HDF5 file
        f = h5py.File(filename, 'r+')
        # Read basic definitions
        ds = f['implements']
        implements = ds[...]

        if 'Mantis' in f:
            mantisGrp = f['Mantis']
        else:
            del f['implements']
            implements = implements+':Mantis'
            ds = f.create_dataset('implements', data = implements)
            mantisGrp = f.create_group('Mantis')


        if anlz.pca_calculated == 1:
            if 'pca' in mantisGrp:
                del mantisGrp['pca']
                pcaGrp = mantisGrp.create_group('pca')
            else:
                pcaGrp = mantisGrp.create_group('pca')
            ds_data = pcaGrp.create_dataset('pca_images', data = anlz.pcaimages)
            ds_data = pcaGrp.create_dataset('pca_eigenvectors', data = anlz.eigenvecs)
            ds_data = pcaGrp.create_dataset('pca_eigenvalues', data = anlz.eigenvals)


        if anlz.clusters_calculated == 1:
            if 'cluster_analysis' in mantisGrp:
                del mantisGrp['cluster_analysis']
                caGrp = mantisGrp.create_group('cluster_analysis')
            else:
                caGrp = mantisGrp.create_group('cluster_analysis')
            ds_data = caGrp.create_dataset('cluster_indices', data = anlz.cluster_indices)
            ds_data = caGrp.create_dataset('cluster_spectra', data = anlz.clusterspectra)
            ds_data = caGrp.create_dataset('cluster_distances', data = anlz.cluster_distances)

        if anlz.tspectrum_loaded == 1:
            if 'spectral_maps' in mantisGrp:
                del mantisGrp['spectral_maps']
                caGrp = mantisGrp.create_group('spectral_maps')
            else:
                caGrp = mantisGrp.create_group('spectral_maps')
            ds_data = caGrp.create_dataset('raw_maps', data = anlz.target_svd_maps)
            ds_data = caGrp.create_dataset('fitted_maps', data = anlz.target_pcafit_maps)
            ds_data = caGrp.create_dataset('raw_spectra', data = anlz.target_spectra)
            ds_data = caGrp.create_dataset('fitted_spectra', data = anlz.target_pcafit_spectra)

        f.close()


        return 0








