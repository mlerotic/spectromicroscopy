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

import numpy
import os


title = 'text table'
extension = ['*.csv','*.txt','*.xas']
read_types = ['spectrum','stack']
write_types = ['spectrum','stack']


def identify(filename):
    try:
        with open(filename, 'tr') as check_file:  # try open file in text mode and read a single line
            check_file.readline()
            return True
    except:
        return False

def GetFileStructure(FileName):
    return None


#----------------------------------------------------------------------
def read(filename, self, selection=None, *args, **kwargs):

    with open(str(filename),'r') as f:
        Line = f.readline().split()
    if Line == ['#X', '#Y', '#Wave', '#Intensity']:
        read_stack(filename, self, **kwargs)
    else:
        read_spectrum(filename, self, **kwargs)
    
    return

def read_stack(filename, self, selection=None, *args, **kwargs):
    """
    This reads the text version of an SPC hyperspectral map from a Raman measurement.
    """

    f = open(str(filename),'r')
    
    xlist = []
    ylist = []    
    wlist = []
    ilist = []    

    for line in f:
        if line.startswith(("*","%","#")):
            pass
        else:
            x, y, w, i = [float (x) for x in line.split()] 
            xlist.append(x)
            ylist.append(y)
            wlist.append(w)
            ilist.append(i)
    
    self.x_dist = numpy.unique(xlist)
    self.y_dist = numpy.unique(ylist)
    self.ev = numpy.unique(wlist)

    self.n_cols = len(self.x_dist)
    self.n_rows = len(self.y_dist)
    self.n_ev = len(self.ev)

    self.data_dwell = numpy.ones((self.n_ev))*1.0
    self.absdata = numpy.empty((self.n_cols,self.n_rows, self.n_ev))

    self.absdata = numpy.flip(numpy.transpose(numpy.reshape(ilist, (self.n_rows, self.n_cols, self.n_ev), order='C'), axes=[1,0,2]), axis=2)

    self.fill_h5_struct_from_stk()


    return

def read_spectrum(filename, self, selection=None, *args, **kwargs):

    f = open(str(filename),'r')
    
    elist = []
    ilist = []    

    for line in f:
        if line.startswith(("*","%","#")):
            pass
        else:
            e, i = [float (x) for x in line.split(',')] 
            elist.append(e)
            ilist.append(i)
            
    self.evi0 = numpy.array(elist)
    self.i0data = numpy.array(ilist) 
            
    f.close()
    
    self.i0_dwell = None
       
    return
    
#----------------------------------------------------------------------
def write(filename, data_object, data_type):
    """Switchyard for writing different types of data."""
    if data_type in ['spectrum']:
        write_spectrum(filename, data_object.absdata, data_object.ev)
    elif data_type in ['stack']:
        write_spectrum(filename, numpy.average(data_object.absdata,axis=(0,1)), energies=data_object.ev, title='I0 Spectrum')

#----------------------------------------------------------------------
def write_spectrum(filename, data, energies, title=None ):
    if title is None:
        title = 'Spectrum'
    
    with open(filename, 'w') as f:
        print('********************* '+title+'  ********************', file=f)
        print('*', file=f)
        print('* ev, intensity', file=f)
        for ie in range(len(energies)):
            print('{0:06.2f}, {1:06f}'.format(energies[ie], data[ie]), file=f)



