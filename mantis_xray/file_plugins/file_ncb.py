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
import os
import re

verbose = 0

title = 'Ncb'
extension = ['*.ncb']
read_types = ['image', 'stack']
write_types = ['image', 'stack', 'results']


def identify(filename):
    try:
        basename, extension = os.path.splitext(filename)
        dat_fn = basename + '.dat'
        f = open(str(dat_fn), 'r')
        f.close()
        return True
    except:
        return False


def GetFileStructure(FileName):
    return None


# ----------------------------------------------------------------------
def read(filename, self, selection=None, *args, **kwargs):
    basename, extension = os.path.splitext(filename)
    dat_fn = basename + '.dat'

    f = open(str(dat_fn), 'r')
    lines = f.readlines()
    temp = lines[0].split()

    self.n_cols = int(temp[0])
    self.n_rows = int(temp[1])
    scale = float(temp[2])

    # print self.n_cols, self.n_rows, scale

    temp = lines[1].split()
    x_start = float(temp[0])
    x_stop = float(temp[1])

    temp = lines[2].split()
    y_start = float(temp[0])
    y_stop = float(temp[1])

    self.n_ev = int(lines[3])

    # print self.n_ev

    self.ev = np.zeros((self.n_ev))
    self.data_dwell = np.zeros((self.n_ev))
    filename_list = []
    for i in range(self.n_ev):
        self.ev[i] = float(lines[4 + i])

    # print self.ev

    for i in range(self.n_ev):
        self.ev[i] = float(lines[4 + i])

    try:
        for i in range(self.n_ev):
            temp = lines[4 + self.n_ev + i].split()

            filename_list.append(temp[0])
            self.data_dwell[i] = float(temp[2])
    except:
        self.data_dwell = np.ones((self.n_ev))

    f.close()

    if scale < 0:
        dataformat = np.float32
    else:
        dataformat = np.int16

    # print 'data format', dataformat

    f = open(str(filename), 'rb')
    big_array = np.fromfile(f, dataformat, self.n_cols * self.n_rows * self.n_ev)
    f.close()

    if scale > 0 and scale != 1:
        image_stack = big_array.astype(np.float) / scale
        print("data rescaled by ", 1. / scale)
    else:
        image_stack = big_array.astype(np.float)

    if (x_start > x_stop):
        image_stack = image_stack[::-1, :, :]
        t = x_start
        x_start = x_stop
        x_stop = t
        print('x data reversed')

    if (y_start > y_stop):
        image_stack = image_stack[:, ::-1, :]
        t = y_start
        y_start = y_stop
        y_stop = t
        print('y data reversed')

    xstep = (x_stop - x_start) / (self.n_cols - 1)
    self.x_dist = np.arange(x_start, x_stop + xstep, xstep)

    ystep = (y_stop - y_start) / (self.n_rows - 1)
    self.y_dist = np.arange(y_start, y_stop + ystep, ystep)

    self.absdata = np.empty((self.n_cols, self.n_rows, self.n_ev))

    self.absdata = np.reshape(image_stack, (self.n_cols, self.n_rows, self.n_ev), order='F')

    self.fill_h5_struct_from_stk()

    return


# ----------------------------------------------------------------------
def read_ncb4D(self, filenames):
    if verbose: print('First energy stack:', filenames[0])

    data_files = natural_sort(filenames)

    neng = len(data_files)
    if verbose: print('Number of energies ', neng)

    self.stack4D = []
    theta = []

    for i in range(neng):
        # Read the data
        thisstack, thistheta = read_ncb_data(self, data_files[i])

        # Check if we have negative angles, if yes convert to 0-360deg
        for ith in range(len(thistheta)):
            if thistheta[ith] < 0:
                thistheta[ith] = thistheta[ith] + 360

        self.stack4D.append(thisstack)
        theta.append(thistheta)

        dims = thisstack.shape

    self.stack4D = np.array(self.stack4D)
    self.stack4D = np.transpose(self.stack4D, axes=(1, 2, 0, 3))

    self.theta = theta[0]
    self.n_theta = len(theta)

    self.absdata = self.stack4D[:, :, :, 0]

    self.data_dwell = np.ones((neng))

    self.fill_h5_struct_from_stk()

    return


# ----------------------------------------------------------------------
def read_ncb_data(self, filename):
    basename, extension = os.path.splitext(filename)
    dat_fn = basename + '.dat'

    f = open(str(dat_fn), 'r')
    lines = f.readlines()
    temp = lines[0].split()

    self.n_cols = int(temp[0])
    self.n_rows = int(temp[1])
    scale = float(temp[2])

    if verbose:
        print('self.n_cols, self.n_rows, scale', self.n_cols, self.n_rows, scale)

    temp = lines[1].split()
    x_start = float(temp[0])
    x_stop = float(temp[1])

    temp = lines[2].split()
    y_start = float(temp[0])
    y_stop = float(temp[1])

    self.n_theta = int(lines[3])

    angles = np.zeros((self.n_theta))
    self.data_dwell = np.zeros((self.n_theta))
    filename_list = []
    for i in range(self.n_theta):
        angles[i] = float(lines[4 + i])

    for i in range(self.n_theta):
        angles[i] = float(lines[4 + i])

    try:
        for i in range(self.n_theta):
            temp = lines[4 + self.n_ev + i].split()

            filename_list.append(temp[0])
            self.data_dwell[i] = float(temp[2])
    except:
        self.data_dwell = np.ones((self.n_theta))

    f.close()

    if scale < 0:
        dataformat = np.float32
    else:
        dataformat = np.int16

    f = open(str(filename), 'rb')
    big_array = np.fromfile(f, dataformat, self.n_cols * self.n_rows * self.n_theta)
    f.close()

    if scale > 0 and scale != 1:
        image_stack = big_array.astype(np.float) / scale
        if verbose: print("data rescaled by ", 1. / scale)
    else:
        image_stack = big_array.astype(np.float)

    if (x_start > x_stop):
        image_stack = image_stack[::-1, :, :]
        t = x_start
        x_start = x_stop
        x_stop = t
        if verbose: print('x data reversed')

    if (y_start > y_stop):
        image_stack = image_stack[:, ::-1, :]
        t = y_start
        y_start = y_stop
        y_stop = t
        if verbose: print('y data reversed')

    xstep = (x_stop - x_start) / (self.n_cols - 1)
    self.x_dist = np.arange(x_start, x_stop + xstep, xstep)

    ystep = (y_stop - y_start) / (self.n_rows - 1)
    self.y_dist = np.arange(y_start, y_stop + ystep, ystep)

    image_stack = np.reshape(image_stack, (self.n_cols, self.n_rows, self.n_theta), order='F')

    return image_stack, angles


# ----------------------------------------------------------------------
#  This procedure writes  a whole stack (3d (E,x,y) array) to a binary file
#  with associated *.dat file to track paramaters
def write(filename, stack, data_type):  # ,norm):
    #
    # def operate_on_Narray(A, B, function):
    #     try:
    #         return [operate_on_Narray(a, b, function) for a, b in zip(A, B)]
    #     except TypeError as e:
    #         # Not iterable
    #         return function(A, B)
    # import matplotlib.pyplot as plt
    # def export_figure_matplotlib(name,arr, dpi=100, resize_fact=1, plt_show=False):
    #     fig = plt.figure(frameon=False)
    #     fig.set_size_inches(arr.shape[1] / dpi, arr.shape[0] / dpi)
    #     ax = plt.Axes(fig, [0., 0., 1., 1.])
    #     ax.set_axis_off()
    #     fig.add_axes(ax)
    #     ax.imshow(arr, cmap='gray', interpolation='none')
    #     if plt_show:
    #         plt.show()
    #     else:
    #         plt.savefig(name+'.tif', dpi=(dpi * resize_fact))
    #         plt.close()
    print('Writing .ncb stack:', filename)

    basename, extension = os.path.splitext(filename)
    dat_fn = basename + '.dat'
    filename = basename + '.ncb'

    image_stack = np.transpose(stack.absdata, axes=(1, 0, 2))
    '''Experimental OD filtering and edge substraction'''
    # if norm:
    #
    #     image_stack = np.transpose(stack.od3d, axes=(1,0,2))
    #     #image_stack = stack.od3d.copy()
    #     OD_filter = np.amax(image_stack, 2)
    #     #returning the ODmask as tif
    #     mask = np.where((OD_filter > 2.5), 0, np.ones(np.shape(image_stack[:,:,0])))
    #     export_figure_matplotlib(basename+"od_mask",mask)
    #
    #     for i in range(np.size(image_stack, 2)):
    #         image_stack[:,:,i] = np.where((OD_filter > 2.5),0, image_stack[:,:,i])
    #     #image_stack = np.where(image_stack > 2,0, image_stack) #OD > 2 filter
    #     e0 = 0
    #     num = 5
    #     ##pre-edge
    #     ave = image_stack[:,:,e0]
    #     for i in range(num-1):
    #         ave = ave + image_stack[:,:,e0+i+1]
    #     ave = ave / num
    #
    #     #export pre-edge_average for contour plot.
    #     pre = np.transpose(stack.absdata, axes=(1,0,2))[:,:,e0]
    #     for i in range(num-1):
    #         pre = pre + np.transpose(stack.absdata, axes=(1,0,2))[:,:,e0+i+1]
    #     pre = pre / num
    #     export_figure_matplotlib(basename,pre)
    #
    #     for i in range(np.size(image_stack, 2)):
    #         image_stack[:,:,i] = image_stack[:,:,i]-ave
    #
    #     #post-edge
    #     ave = image_stack[:,:,-1]
    #     for i in range(num-1):
    #         ave = ave + image_stack[:,:,e0-i-2]
    #     ave = (ave / num) + 0.0001
    #     #print(ave)
    #     for i in range(np.size(image_stack, 2)):
    #         image_stack[:,:,i] = image_stack[:,:,i] - operate_on_Narray(ave, image_stack[:,:,i], lambda a, b: a/((1+ np.exp(-b+0.5*a)**(25*1/a))))

    image_stack = np.reshape(image_stack, (stack.n_cols * stack.n_rows * stack.n_ev), order='F')

    tmax = np.amax(image_stack)
    tmin = np.amin(image_stack)
    test = np.amax([abs(tmin), abs(tmax)])
    scale = 1.
    if test > 3e4 or test < 1e3:
        scale = 10. ** (3 - np.fix(np.math.log10(test)))

    # Save image stack to a binary .dat file
    f = open(str(filename), 'wb')
    if scale != 1.0:
        print('Scaling the data by ', scale)
        saveddata = image_stack * scale
        saveddata.astype(np.int16).tofile(f)
    else:
        scale = -1.0
        image_stack.astype(np.float32).tofile(f)
    f.close()

    # print 'imagedims', image_stack.shape

    f = open(str(dat_fn), 'w')

    print('\t%d\t%d\t%.6f' % (stack.n_rows, stack.n_cols, scale), file=f)

    x_start = stack.y_dist[0]
    x_stop = stack.y_dist[-1]
    if x_start != 0.:
        x_stop = x_stop - x_start
        x_start = 0.

    print('\t%.6f\t%.6f' % (x_start, x_stop), file=f)

    y_start = stack.x_dist[0]
    y_stop = stack.x_dist[-1]
    if y_start != 0.:
        y_stop = y_stop - y_start
        y_start = 0.
    print('\t%.6f\t%.6f' % (y_start, y_stop), file=f)

    print('\t%d' % (stack.n_ev), file=f)

    for i in range(stack.n_ev):
        print('\t%.6f' % (stack.ev[i]), file=f)

    for i in range(stack.n_ev):
        thisstr = 'image' + str(i + 1)
        print('%s\t%.6f\t%.6f' % (thisstr, stack.ev[i], stack.data_dwell[i]), file=f)

    f.close()

    return


# ----------------------------------------------------------------------
def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


# ----------------------------------------------------------------------
class Cncb:
    def __init__(self):
        pass

    # ----------------------------------------------------------------------
    def read_ncb(self, filename):

        print('Reading stack:', filename)

        basename, extension = os.path.splitext(filename)
        dat_fn = basename + '.dat'

        f = open(str(dat_fn), 'r')
        lines = f.readlines()
        temp = lines[0].split()

        self.n_cols = int(temp[0])
        self.n_rows = int(temp[1])
        scale = float(temp[2])

        # print self.n_cols, self.n_rows, scale

        temp = lines[1].split()
        x_start = float(temp[0])
        x_stop = float(temp[1])

        temp = lines[2].split()
        y_start = float(temp[0])
        y_stop = float(temp[1])

        self.n_ev = int(lines[3])

        # print self.n_ev

        self.ev = np.zeros((self.n_ev))
        self.data_dwell = np.zeros((self.n_ev))
        filename_list = []
        for i in range(self.n_ev):
            self.ev[i] = float(lines[4 + i])

        # print self.ev

        for i in range(self.n_ev):
            self.ev[i] = float(lines[4 + i])

        try:
            for i in range(self.n_ev):
                temp = lines[4 + self.n_ev + i].split()

                filename_list.append(temp[0])
                self.data_dwell[i] = float(temp[2])
        except:
            self.data_dwell = np.ones((self.n_ev))

        f.close()

        if scale < 0:
            dataformat = np.float32
        else:
            dataformat = np.int16

        # print 'data format', dataformat

        f = open(str(filename), 'rb')
        big_array = np.fromfile(f, dataformat, self.n_cols * self.n_rows * self.n_ev)
        f.close()

        if scale > 0 and scale != 1:
            image_stack = big_array.astype(np.float) / scale
            print("data rescaled by ", 1. / scale)
        else:
            image_stack = big_array.astype(np.float)

        if (x_start > x_stop):
            image_stack = image_stack[::-1, :, :]
            t = x_start
            x_start = x_stop
            x_stop = t
            print('x data reversed')

        if (y_start > y_stop):
            image_stack = image_stack[:, ::-1, :]
            t = y_start
            y_start = y_stop
            y_stop = t
            print('y data reversed')

        xstep = (x_stop - x_start) / (self.n_cols - 1)
        self.x_dist = np.arange(x_start, x_stop + xstep, xstep)

        ystep = (y_stop - y_start) / (self.n_rows - 1)
        self.y_dist = np.arange(y_start, y_stop + ystep, ystep)

        self.absdata = np.empty((self.n_cols, self.n_rows, self.n_ev))

        self.absdata = np.reshape(image_stack, (self.n_cols, self.n_rows, self.n_ev), order='F')

        return

    # ----------------------------------------------------------------------
    #  This procedure writes  a whole stack (3d (E,x,y) array) to a binary file
    #  with associated *.dat file to track paramaters
    def write_ncb(self, filename, data_struct):

        print('Writing .ncb stack:', filename)

        basename, extension = os.path.splitext(filename)
        dat_fn = basename + '.dat'

        image_stack = np.transpose(self.absdata, axes=(1, 0, 2))

        # image_stack = self.absdata.copy()

        image_stack = np.reshape(image_stack, (self.n_cols * self.n_rows * self.n_ev), order='F')

        tmax = np.amax(image_stack)
        tmin = np.amin(image_stack)
        test = np.amax([abs(tmin), abs(tmax)])
        scale = 1.
        if test > 3e4 or test < 1e3:
            scale = 10. ** (3 - np.fix(np.math.log10(test)))

        # Save image stack to a binary .dat file
        f = open(str(filename), 'wb')
        if scale != 1.0:
            print('Scaling the data by ', scale)
            saveddata = image_stack * scale
            saveddata.astype(np.int16).tofile(f)
        else:
            scale = -1.0
            image_stack.astype(np.float32).tofile(f)
        f.close()

        # print 'imagedims', image_stack.shape

        f = open(str(dat_fn), 'w')

        print('\t%d\t%d\t%.6f' % (self.n_rows, self.n_cols, scale), file=f)
        # print('\t%d\t%d\t%.6f' %(self.n_cols, self.n_rows, scale)

        x_start = self.y_dist[0]
        x_stop = self.y_dist[-1]
        if x_start != 0.:
            x_stop = x_stop - x_start
            x_start = 0.

        print('\t%.6f\t%.6f' % (x_start, x_stop), file=f)

        y_start = self.x_dist[0]
        y_stop = self.x_dist[-1]
        if y_start != 0.:
            y_stop = y_stop - y_start
            y_start = 0.
        print('\t%.6f\t%.6f' % (y_start, y_stop), file=f)

        print('\t%d' % (self.n_ev), file=f)

        for i in range(self.n_ev):
            print('\t%.6f' % (self.ev[i]), file=f)

        for i in range(self.n_ev):
            thisstr = 'image' + str(i + 1)
            print('%s\t%.6f\t%.6f' % (thisstr, self.ev[i], self.data_dwell[i]), file=f)

        f.close()

        return

    # ----------------------------------------------------------------------
    def read_ncb_data(self, filename):

        # print 'Reading stack:', filename

        basename, extension = os.path.splitext(filename)
        dat_fn = basename + '.dat'

        f = open(str(dat_fn), 'r')
        lines = f.readlines()
        temp = lines[0].split()

        self.n_cols = int(temp[0])
        self.n_rows = int(temp[1])
        scale = float(temp[2])

        if verbose:
            print('self.n_cols, self.n_rows, scale', self.n_cols, self.n_rows, scale)

        temp = lines[1].split()
        x_start = float(temp[0])
        x_stop = float(temp[1])

        temp = lines[2].split()
        y_start = float(temp[0])
        y_stop = float(temp[1])

        self.n_theta = int(lines[3])

        # print self.n_ev

        angles = np.zeros((self.n_theta))
        self.data_dwell = np.zeros((self.n_theta))
        filename_list = []
        for i in range(self.n_theta):
            angles[i] = float(lines[4 + i])

        # print self.ev

        for i in range(self.n_theta):
            angles[i] = float(lines[4 + i])

        try:
            for i in range(self.n_theta):
                temp = lines[4 + self.n_ev + i].split()

                filename_list.append(temp[0])
                self.data_dwell[i] = float(temp[2])
        except:
            self.data_dwell = np.ones((self.n_theta))

        f.close()

        if scale < 0:
            dataformat = np.float32
        else:
            dataformat = np.int16

        # print 'data format', dataformat

        f = open(str(filename), 'rb')
        big_array = np.fromfile(f, dataformat, self.n_cols * self.n_rows * self.n_theta)
        f.close()

        if scale > 0 and scale != 1:
            image_stack = big_array.astype(np.float) / scale
            if verbose: print("data rescaled by ", 1. / scale)
        else:
            image_stack = big_array.astype(np.float)

        if (x_start > x_stop):
            image_stack = image_stack[::-1, :, :]
            t = x_start
            x_start = x_stop
            x_stop = t
            if verbose: print('x data reversed')

        if (y_start > y_stop):
            image_stack = image_stack[:, ::-1, :]
            t = y_start
            y_start = y_stop
            y_stop = t
            if verbose: print('y data reversed')

        xstep = (x_stop - x_start) / (self.n_cols - 1)
        self.x_dist = np.arange(x_start, x_stop + xstep, xstep)

        ystep = (y_stop - y_start) / (self.n_rows - 1)
        self.y_dist = np.arange(y_start, y_stop + ystep, ystep)

        image_stack = np.reshape(image_stack, (self.n_cols, self.n_rows, self.n_theta), order='F')

        return image_stack, angles

    # ----------------------------------------------------------------------
    def read_ncb4D(self, filenames):

        if verbose: print('First energy stack:', filenames[0])

        data_files = natural_sort(filenames)

        neng = len(data_files)
        if verbose: print('Number of energies ', neng)

        self.stack4D = []
        theta = []

        for i in range(neng):
            # Read the data
            thisstack, thistheta = self.read_ncb_data(data_files[i])

            # Check if we have negative angles, if yes convert to 0-360deg
            for ith in range(len(thistheta)):
                if thistheta[ith] < 0:
                    thistheta[ith] = thistheta[ith] + 360

            self.stack4D.append(thisstack)
            theta.append(thistheta)

            dims = thisstack.shape

        self.stack4D = np.array(self.stack4D)
        self.stack4D = np.transpose(self.stack4D, axes=(1, 2, 0, 3))

        self.theta = theta[0]
        self.n_theta = len(theta)

        self.absdata = self.stack4D[:, :, :, 0]

        self.data_dwell = np.ones((neng))

        return
