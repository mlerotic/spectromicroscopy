# -*- coding: utf-8 -*-
# 
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
# 
#   Copyright (C) 2015 Benjamin Watts, Paul Scherrer Institute
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


import sys, os, numpy, h5py
from collections import OrderedDict
from PyQt4 import QtGui

title = 'NXstxm'
extension = ['*.hdf','*.hdf5','*.nxs']
read_types = ['spectrum','image','stack']
write_types = ['spectrum','image','stack']

def identify(filename):
    try:
        # Open HDF5 file
        f = h5py.File(filename, 'r')
        # Count valid entries
        n_regions = 0
        for entry in f:
            if 'NX_class' in f[entry].attrs and f[entry].attrs['NX_class'] == 'NXentry':
                if 'definition' in f[entry] and f[entry]['definition'][0] == 'NXstxm':
                    n_regions += 1
        f.close()
        return n_regions>0 # return true if file contains at least one NXstxm entry
    except:
        return False

def read(FileName,stack_object,selection=(0,0)):
    D = GetFileStructure(FileName)
    entry = D.keys()[selection[0]]
    detector = D[entry].keys()[selection[1]]
    F = h5py.File(FileName, 'r')
    if 'energy' in list(F[entry][detector]):
        stack_object.ev = numpy.array(F[entry][detector]['energy'])
    elif 'photon_energy' in list(F[entry][detector]):
        stack_object.ev = numpy.array(F[entry][detector]['photon_energy'])
    else:
        print "Can't find photon energy!"
    stack_object.x_dist = numpy.array(F[entry][detector]['sample_x'])
    stack_object.y_dist = numpy.array(F[entry][detector]['sample_y'])
    stack_object.data_dwell = numpy.array(F[entry][detector]['count_time'])
    stack_object.n_cols = len(stack_object.x_dist)
    stack_object.n_rows = len(stack_object.y_dist)
    stack_object.n_ev = len(stack_object.ev)
    if 'axes' in F[entry][detector].attrs: # Specification correct
        axes_list = list(F[entry][detector].attrs['axes'])
        axes_order = [axes_list.index('sample_x'),axes_list.index('sample_y'),axes_list.index('energy')]
    else: # Old version from before the specification was finalised
        if 'energy' in list(F[entry][detector]):
            energy_axis = F[entry][detector]['energy'].attrs['axis']
        elif 'photon_energy' in list(F[entry][detector]):
            energy_axis = F[entry][detector]['photon_energy'].attrs['axis']
        else:
            print "Can't find photon energy!"
        axes_order = [F[entry][detector]['sample_x'].attrs['axis']-1,F[entry][detector]['sample_y'].attrs['axis']-1,energy_axis-1]
    stack_object.absdata = numpy.transpose(numpy.array(F[entry][detector]['data']),axes=axes_order)
    
    
    F.close()
    
    stack_object.fill_h5_struct_from_stk()


def GetFileStructure(FileName):
    F = h5py.File(FileName, 'r')
    D = OrderedDict()
    for entry in F:
        if 'NX_class' in F[entry].attrs and F[entry].attrs['NX_class'] == 'NXentry':
            D[entry] = OrderedDict()
            D[entry].norm_data = OrderedDict()
            D[entry].definition = None
            D[entry].scan_type = None
            D[entry].data_shape = None
            D[entry].data_axes = None
            for data in F[entry]:
                if 'NX_class' in F[entry][data].attrs and F[entry][data].attrs['NX_class'] == 'NXdata':
                    D[entry][data] = OrderedDict()
                    #print "should collect more info in each NXdata group"
                elif 'NX_class' in F[entry][data].attrs and F[entry][data].attrs['NX_class'] == 'NXmonitor':
                    D[entry].norm_data[data] = OrderedDict()
            if len(D[entry].norm_data) == 0:
                D[entry].norm_data = None
            if 'definition' in F[entry]:
                D[entry].definition = F[entry]['definition'][0]
            if len(D[entry].keys()) > 0:
                channel_zero = D[entry].keys()[0]
                if 'stxm_scan_type' in F[entry][channel_zero]:
                    D[entry].scan_type = F[entry][channel_zero]['stxm_scan_type'][0]
                signal_name = 'data'
                if 'signal' in F[entry][channel_zero].attrs:
                    signal_name = F[entry][channel_zero].attrs['signal']
                if signal_name in F[entry][channel_zero]:
                    D[entry].data_shape = F[entry][channel_zero][signal_name].shape
                if 'axes' in F[entry][channel_zero].attrs:
                    D[entry].data_axes = F[entry][channel_zero].attrs['axes']
    F.close()
    if len(D) == 0:
        return None
    else:
        return D

