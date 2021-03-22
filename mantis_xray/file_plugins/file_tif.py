from __future__ import print_function
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

import numpy as np
from PIL import Image
import os


title = 'Tiff'
extension = ['*.tif','*.tiff']
read_types = ['image','stack']
write_types = ['image','stack','results']


def identify(filename):
    try:
        # Open Tiff file
        img = Image.open(filename)
        return True
    except:
        return False

def GetFileStructure(FileName):
    return None


#----------------------------------------------------------------------
def read(filename, self, selection=None, *args, **kwargs):

    img = Image.open(filename)

    imgstack = []

    i = 0
    while True:
        try:
            img.seek(i)
            i += 1
            imgstack.append(np.array((img)))
        except EOFError:
            break

    imgstack = np.array((imgstack))

    imgstack = np.transpose(imgstack, axes=(1,2,0))

    self.n_cols = imgstack.shape[0]
    self.n_rows = imgstack.shape[1]
    self.n_ev = imgstack.shape[2]

    print('File {0} dims: [{1},{2},{3}]'.format(filename, self.n_cols, self.n_rows, self.n_ev))

    pixelsize = 1
    #Since we do not have a scanning microscope we fill the x_dist and y_dist from pixel_size
    self.x_dist = np.arange(np.float(self.n_cols))*pixelsize
    self.y_dist = np.arange(np.float(self.n_rows))*pixelsize

    #Read energies from file
    basename, extension = os.path.splitext(filename)
    engfilename = basename+'.txt'
    f = open(str(engfilename),'r')

    elist = []

    for line in f:
        if line.startswith("*"):
            pass
        else:
            e = float(line)
            elist.append(e)

    self.ev = np.array(elist)

    f.close()


    msec = np.ones((self.n_ev))

    self.data_dwell = msec

    self.absdata = imgstack

    #Check if the energies are consecutive, if they are not sort the data
    sort = 0
    for i in range(self.n_ev - 1):
        if self.ev[i] > self.ev[i+1]:
            sort = 1
            break
    if sort == 1:
        sortind = np.argsort(self.ev)
        self.ev = self.ev[sortind]
        self.absdata = self.absdata[:,:,sortind]


    self.fill_h5_struct_from_stk()

#----------------------------------------------------------------------
def write(filename, data_object, data_type):
    """Switchyard for writing different types of data."""
    if data_type in ['stack']:
        write_tif(filename, data_object.absdata, data_object.ev)

#----------------------------------------------------------------------
def write_tif(filename, data, energies = []):

    dims = data.shape
    for i in range(dims[2]):
        basename, extension = os.path.splitext(filename)
        thisfn = basename + '_' + str(i+1) + extension

        img1 = Image.fromarray(data[:,:,i])
        img1.save(thisfn)


    if len(energies) > 0:
        thisfn = basename + '.txt'
        f = open(thisfn, 'w')
        for i in range(len(energies)):
            print('%.6f' %(energies[i]), file=f)
        f.close()



