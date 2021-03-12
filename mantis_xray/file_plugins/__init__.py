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

"""
The file_plugins system is exposed to general code through the functions defined here in __init__.py:

identify(filename)                  : Returns an instance of the plugin that claims to deal with the file at the URL 'filename'.
GetFileStructure(filename)          : Returns a structure describing the internal organisation of the file, indicating sets of data available to choose from.
load(filename,stack_object,..)      : Loads data from the URL 'filename' into the object (data_stack type) 'stack_object'. The plugin used can be stated or determined automatically, using 'identify'.

Further functions for writing files via the plugins need to be written yet. To access the system, you should import the module ('import file_plugins') and then access the above functions as attributes of the module (e.g. 'file_plugins.load('data.hdf5',data_stk)' ).

Each file plugin should be included here in the 'file_plugins' directory. Each plugin should define the following:

title                           : A short string naming the plugin.
extension                       : A list of strings indicating the file extensions that the plugin handles (e.g. ['*.txt']).
read_types                      : A list of strings indicating the data types that the plugin will read (e.g. ['spectrum','image','stack']).
write_types                     : A list of strings indicating the data types that the plugin will write (e.g. ['spectrum','image','stack']).
identify(filename)              : Returns boolean indicating if the plugin can read the file at URL 'filename'.
GetFileStructure(filename)      : Returns a structure describing the internal organisation of the file, indicating sets of data available to choose from.
read(filename,stack_object,..)  : Loads data from the URL 'filename' into the object (data_stack type) 'stack_object'.

"""
from __future__ import print_function

import pkgutil, imp, os, sys
import numpy
from .. import data_stack

verbose = True

# These variables declare the options that each plugin can claim the ability to handle
actions = ['read','write']
data_types = ['spectrum','image','stack','results']

# Go through the directory and try to load each plugin
plugins = []

for m in pkgutil.iter_modules(path=__path__):
    if verbose: print("Loading file plugin:", m[1], ".", end=' ')
    try:
        details = imp.find_module(m[1],__path__)
        # check if there is a read() function in plugin
        if 'read' in dir(imp.load_module(m[1],*details)):
            plugins.append(imp.load_module(m[1],*details))
            if verbose: print("("+plugins[-1].title+") Success!")
        else:
            if verbose: print('Not a valid plugin - skipping.')

    except ImportError as e:
        if verbose: print("prerequisites not satisfied:", e)


# if getattr(sys, 'frozen', False):
#     module_names = ['file_dataexch_hdf5', 'file_ncb', 'file_nexus_hdf5', 'file_sdf', 'file_stk', 'file_tif', 'file_xrm']
#     for m in module_names:
#         if verbose: print "Loading file plugin:", m, "...",
#         try:
#
#
#             details = imp.find_module(m)
#             # check if there is a read() function in plugin
#             if 'read' in dir(imp.load_module(m,*details)):
#                 plugins.append(imp.load_module(m,*details))
#                 if verbose: print "("+plugins[-1].title+") Success!"
#             else:
#                 if verbose: print 'Not a valid plugin - skipping.'
#
#         except ImportError as e:
#             if verbose: print "prerequisites not satisfied:", e




# Go through set of plugins and assemble lists of supported file types for each action and data type
supported_filters = dict([a,dict([t,[]] for t in data_types)] for a in actions)
supported_plugins = dict([a,dict([t,[]] for t in data_types)] for a in actions)
filter_list = dict([a,dict([t,[]] for t in data_types)] for a in actions)
for P in plugins:
    for action in actions:
        for data_type in data_types:
            if data_type in getattr(P,action+'_types'):
                filter_list[action][data_type].append(P.title+' ('+' '.join(P.extension)+')')
                supported_plugins[action][data_type].append(P)
                for ext in P.extension:
                    if ext not in supported_filters[action][data_type]:
                        supported_filters[action][data_type].append(ext)
for data_type in data_types:
    filter_list['read'][data_type] = ['Supported Formats ('+' '.join(supported_filters['read'][data_type])+')']+filter_list['read'][data_type]
    filter_list['read'][data_type].append('All files (*.*)')


def load(filename, stack_object=None, plugin=None, selection=None, json=None):
    """
    Pass the load command over to the appropriate plugin so that it can import data from the named file.
    """
    if plugin is None:
        plugin = identify(filename)
    if plugin is None:
        return None
    else:
        print("load", filename, "with the", plugin.title, "plugin.")
        if selection is None:
            plugin.read(filename, stack_object, selection, json)
        elif len(selection) == 1:
            plugin.read(filename, stack_object, selection[0], json)
        else:
            plugin.read(filename,stack_object,selection[0],json)
            temp_stack = data_stack.data(stack_object.data_struct)
            full_stack = stack_object.absdata.copy()
            for s in selection[1:]:
                plugin.read(filename,temp_stack,s)
                if full_stack.shape[1] > temp_stack.absdata.shape[1]:
                    temp_stack.absdata = numpy.pad(temp_stack.absdata,((0,0),(0,full_stack.shape[1]-temp_stack.absdata.shape[1]),(0,0)), mode='constant',constant_values=0)
                elif full_stack.shape[1] < temp_stack.absdata.shape[1]:
                    full_stack = numpy.pad(full_stack,((0,0),(0,temp_stack.absdata.shape[1]-full_stack.shape[1]),(0,0)), mode='constant',constant_values=0)
                full_stack = numpy.vstack((full_stack,temp_stack.absdata))
            stack_object.absdata = full_stack
            stack_object.x_dist = numpy.arange(full_stack.shape[0])
            stack_object.y_dist = numpy.arange(full_stack.shape[1])
            stack_object.n_cols = len(stack_object.x_dist)
            stack_object.n_rows = len(stack_object.y_dist)
            return

def save(filename, data_object, data_type, plugin=None):
    """
    Pass the save command over to the appropriate plugin so that it can write data to the named file.
    """
    print("save", filename, "with the", plugin.title, "plugin.")
    plugin.write(filename, data_object, data_type)

def GetFileStructure(filename, plugin=None):
    """
    Use the plugin to skim-read the file and return the structure of the data.
    Returns None if there is only a single data array (i.e. no choices to be made).
    """
    if plugin is None:
        plugin = identify(filename)
    if plugin is None:
        return None
    else:
        print("get info from", filename, "with the", plugin.title, "plugin.")
        FileInfo = plugin.GetFileStructure(filename)
        #if FileInfo is not None:
            #print len(FileInfo), len(FileInfo[next(iter(FileInfo))])
            #print FileInfo
        return FileInfo

def identify(filename):
    """
    Cycle through plugins until finding one that claims to understand the file format.
    First it tries those claiming corresponding file extensions, followed by all other plugins until an appropriate plugin is found.
    """
    print("Identifying file:", filename, "...", end=' ')
    ext = os.path.splitext(filename)[1]
    flag = [True]*len(plugins)
    for i,P in enumerate(plugins):
        if '*'+ext in P.extension:
            if P.identify(filename):
                return P
            elif flag[i] == True: #if plugin returns False, e.g. dataexch_hdf5 does not match, try the next plugin and find the same extension
                flag[i] = False
                continue
            else:
                break
    print("Error! unknown file type.")
    return None

