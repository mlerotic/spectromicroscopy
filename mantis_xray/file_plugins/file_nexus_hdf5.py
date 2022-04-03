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


from __future__ import print_function
import sys, os, numpy, h5py
from collections import OrderedDict
#from PyQt5 import QtGui

title = 'NXstxm'
extension = ['*.hdf','*.hdf5','*.nxs']
read_types = ['spectrum','image','stack']
write_types = []#'spectrum','image','stack']

def perhaps_decode(value):
	try:
		return value.decode("utf-8")	# h5py <3.0
	except AttributeError:
		return value					# h5py >=3.0

def identify(filename):
    try:
        # Open HDF5 file
        f = h5py.File(filename, 'r')
        # Count valid entries
        n_regions = 0
        for entry in f:
            if 'NX_class' in list(f[entry].attrs) and f[entry].attrs['NX_class'] in [b'NXentry', u'NXentry']:
                if 'definition' in list(f[entry]) and f[entry]['definition'][0] in [b'NXstxm', u'NXstxm']:
                    n_regions += 1
        f.close()
        return n_regions>0 # return true if file contains at least one NXstxm entry
    except:
        return False

def read(FileName,stack_object,selection=(0,0), *args, **kwargs):
    """Todo: Add support for single images!"""
    D = GetFileStructure(FileName)
    entry = list(D.keys())[selection[0]]
    detector = list(D[entry].keys())[selection[1]] #[counter0, ...]
    F = h5py.File(FileName, 'r')
    if 'energy' in list(F[entry][detector]):
        stack_object.ev = numpy.array(F[entry][detector]['energy'])
    elif 'photon_energy' in list(F[entry][detector]):
        stack_object.ev = numpy.array(F[entry][detector]['photon_energy'])
    else:
        print("Can't find photon energy!")
    stack_object.x_dist = numpy.array(F[entry][detector]['sample_x'])
    stack_object.y_dist = numpy.array(F[entry][detector]['sample_y'])
    stack_object.data_dwell = numpy.array(F[entry][detector]['count_time'])
    stack_object.n_cols = len(stack_object.x_dist)
    stack_object.n_rows = len(stack_object.y_dist)
    stack_object.n_ev = len(stack_object.ev)
    if 'axes' in list(F[entry][detector].attrs): # Specification correct
        axes_list = [item.decode('UTF-8') for item in F[entry][detector].attrs['axes']]
        if 'line_position' in axes_list:
            axes_order = [axes_list.index('line_position'),axes_list.index('line_position'),axes_list.index('energy')] #linescan
        elif 'energy' in axes_list:
            axes_order = [axes_list.index('sample_x'),axes_list.index('sample_y'),axes_list.index('energy')] #stack
        else:
            axes_order = [axes_list.index('sample_x'),axes_list.index('sample_y')] #image
    else: # Old version from before the specification was finalised
        if 'energy' in list(F[entry][detector]):
            try:
                energy_axis = F[entry][detector]['energy'].attrs['axis']
            except:
                print("Only stacks are supported!")
        elif 'photon_energy' in list(F[entry][detector]):
            energy_axis = F[entry][detector]['photon_energy'].attrs['axis']
        else:
            print("Can't find photon energy!")
        axes_order = [F[entry][detector]['sample_x'].attrs['axis']-1,F[entry][detector]['sample_y'].attrs['axis']-1,energy_axis-1]
    signal_name = 'data'
    if 'signal' in list(F[entry][detector].attrs):
        signal_name = perhaps_decode(F[entry][detector].attrs['signal'])
    if axes_order[0] == axes_order[1]: #i.e. if linescan
        temp = numpy.transpose(numpy.array(F[entry][detector][signal_name]),axes=axes_order[1:])
        stack_object.absdata = numpy.tile(temp,(temp.shape[0],1,1))
    elif len(axes_order)<3: #i.e. if image
        stack_object.absdata = numpy.expand_dims(numpy.transpose(numpy.array(F[entry][detector][signal_name]),axes=axes_order),2)
    else:
        stack_object.absdata = numpy.transpose(numpy.array(F[entry][detector][signal_name]),axes=axes_order)

    F.close()

    stack_object.fill_h5_struct_from_stk()


def GetFileStructure(FileName):
    """ToDo: Currently, the file will be opened two times. Maybe a solution like in the sdf-plugin would be better."""
    F = h5py.File(FileName, 'r')
    D = OrderedDict()
    for entry in F:
        if 'NX_class' in list(F[entry].attrs) and F[entry].attrs['NX_class'] in [b'NXentry', u'NXentry']:
            D[entry] = OrderedDict()
            D[entry].norm_data = OrderedDict()
            D[entry].definition = None
            D[entry].scan_type = None
            D[entry].data_shape = None
            D[entry].data_axes = None
            for data in F[entry]:
                if 'NX_class' in list(F[entry][data].attrs) and F[entry][data].attrs['NX_class'] in [b'NXdata', u'NXdata']:
                    D[entry][data] = OrderedDict()
                    #print "should collect more info in each NXdata group"
                elif 'NX_class' in list(F[entry][data].attrs) and F[entry][data].attrs['NX_class'] in [b'NXmonitor', u'NXmonitor']:
                    D[entry].norm_data[data] = OrderedDict()
            if len(D[entry].norm_data) == 0:
                D[entry].norm_data = None
            if 'definition' in list(F[entry]):
                D[entry].definition = F[entry]['definition'][0].decode("utf-8")
            if len(D[entry].keys()) > 0:
                channel_zero = list(D[entry].keys())[0]
                if 'stxm_scan_type' in list(F[entry][channel_zero]):
                    D[entry].scan_type = F[entry][channel_zero]['stxm_scan_type'][0].decode("utf-8")
                signal_name = 'data'
                if 'signal' in list(F[entry][channel_zero].attrs):
                    signal_name = perhaps_decode(F[entry][channel_zero].attrs['signal'])
                if signal_name in list(F[entry][channel_zero]):
                    D[entry].data_shape = F[entry][channel_zero][signal_name].shape
                if 'axes' in list(F[entry][channel_zero].attrs):
                    D[entry].data_axes = [perhaps_decode(item) for item in F[entry][channel_zero].attrs['axes']]
    F.close()
    if len(D) == 0:
        return None
    else:
        return D

