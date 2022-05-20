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


#-----------------------------------------------------------------------
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

def read_tif_info(filename):
    img = Image.open(filename)
    arr = np.array(img)
    n_cols = arr.shape[1]
    n_rows = arr.shape[0]
    return n_cols, n_rows

#-----------------------------------------------------------------------
def read_tif_list(self, filelist, filepath, ds):
    from PyQt5.QtCore import Qt
    n_ev = len(filelist)
    ev = []
    cols = []
    rows =[]
    for i in range(n_ev):
        row = self.filelist.findItems(filelist[i], Qt.MatchContains)[0].row()
        cols.append(int(self.filelist.item(row,1).text()))
        rows.append(int(self.filelist.item(row,2).text()))
        ev.append(float(self.filelist.item(row,3).text()))
    assert rows.count(rows[0]) == len(rows), 'x dimension not equally sized'
    assert cols.count(cols[0]) == len(cols), 'y dimension not equally sized'
    n_rows = rows[0]
    n_cols = cols[0]

    absdata = np.zeros((n_cols, n_rows, n_ev))

    for i in range(n_ev):
        fn = filelist[i]
        filename = os.path.join(filepath, fn)
        img = Image.open(filename)
        imagestack = np.rot90(np.array(img),3)
        absdata[:, :, i] = np.reshape(imagestack, (n_cols, n_rows), order='C')
    # Since we do not know the pixel size, x_dist and y_dist is derived from pixel number
    pixelsize = 1
    x_dist = np.arange(np.float(n_cols))*pixelsize
    y_dist = np.arange(np.float(n_rows))*pixelsize
    # Since unknown, set dwell time to 1
    msec = 1
    data_dwell = np.ones((n_ev))*msec
    #
    # Fill the data structure with data:
    ds.implements = 'information:exchange:spectromicroscopy'
    ds.version = '1.0'
    ds.information.comment = 'Converted from .tif file list in MANTiS',

    import datetime
    now = datetime.datetime.now()
    ds.information.file_creation_datetime = now.strftime("%Y-%m-%dT%H:%M")

    ds.information.experimenter.name = ''
    ds.information.sample.name = ''

    ds.exchange.data = absdata
    ds.exchange.data_signal = 1
    ds.exchange.data_axes = 'x:y'

    ds.exchange.energy = np.array(ev)
    ds.exchange.energy_units = 'eV'

    ds.exchange.x = x_dist
    ds.exchange.x_units = 'um'
    ds.exchange.y = y_dist
    ds.exchange.y_units = 'um'

    ds.spectromicroscopy.data_dwell = data_dwell
    return

def write(filename, data_object, data_type):
    """Switchyard for writing different types of data."""
    if data_type in ['stack']:
        write_tif(filename, data_object.absdata, data_object.ev)

#-----------------------------------------------------------------------
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



