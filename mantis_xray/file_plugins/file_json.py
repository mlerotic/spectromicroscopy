# -*- coding: utf-8 -*-
#
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
#
#   Copyright (C) 2015 Benjamin Watts, Paul Scherrer Institut
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
import json
import re, numpy, sys
from os.path import splitext, isfile
from collections import OrderedDict


title = 'JSON'
extension = ['*.json', '*.hdr']
read_types = ['spectrum','image','stack']
write_types = []
def identify(filename):
    try:
        if isfile(splitext(filename)[0] + '.json'):
            JS = JS_FileLoader(splitext(filename)[0] + '.json')
            return JS.num_regions>0 # return true if file contains at least one region
        else:
            return False
    except:
        return False

def GetFileStructure(FileName):
    JS = JS_FileLoader(FileName)
    if JS.num_regions<2 and JS.num_channels<2 and isfile(splitext(FileName)[0] + '.json'): #if json file exists, skip datachoicedialog
        print("JSON-file found.")
        return None # exit if only one choice
    D = OrderedDict()
    for i,R in enumerate(['Region_'+str(r) for r in range(JS.num_regions)]):
        D[R] = OrderedDict()
        D[R].definition = 'JSON'
        D[R].scan_type = JS.js['ScanDefinition']['Type']
        D[R].data_shape = JS.data_size[i]
        D[R].data_axes = [JS.js['ScanDefinition']['Regions'][1]['PAxis']['Name'],JS.js['ScanDefinition']['Regions'][1]['QAxis']['Name'],JS.js['ScanDefinition']['StackAxis']['Name']]
        for ch in range(1,JS.num_channels+1):
            D[R][JS.js['Channels'][ch]['Name']] = OrderedDict()
    return D



#----------------------------------------------------------------------
class JS_FileLoader:
  """Parse .hdr file for metadata."""
  js = []
  f = []
  def __init__(self, fileName,identify=False):
    fileName = (splitext(fileName)[0] + '.json')
    if not JS_FileLoader.js or JS_FileLoader.f != fileName: # prevent class from fetching and parsing the *.hdr file multiple times. Use the class attribute "hdr" instead if available and check if a new file is loaded.
        # Parse the JSON file
        with open(fileName) as f:
            JS_FileLoader.js = json.load(f)
        JS_FileLoader.f = fileName
    else:
        None
    self.js = JS_FileLoader.js
    if 'ScanDefinition' in self.js:
      self.num_regions = int(self.js['ScanDefinition']['Regions'][0])
      self.num_channels = int(self.js['Channels'][0])
      self.file_path, self.file_ext = splitext(fileName)
      self.data_size = self.parse_DataSize()
      self.data_names = self.parseDataNames()
    else:
      self.num_regions = 0
      self.num_channels = 0

#----------------------------------------------------------------------
  def parseDataNames(self): #ToDo: Assigning Scan types by DataFlags is prone to errors. Better detect automatically via (HDR.hdr['ScanDefinition']['Regions'][0]) etc.
    """Figure out names for the .xsp or .xim files that contain the actual data, then check that the files actually exist, printing warnings if they don't."""
    DataNames = []#Regions
    DataNames2 = []#Channels
    DataNames3 = []#Energies
    Alphabet = 'abcdefghijklmnopqrstuvwxyz'
    DataFlag = self.js['ScanDefinition']['Flags']
    if DataFlag in ['Spectra','Multi-Region Spectra']:
      for num_R in range(self.num_regions):
        DataNames2 = []
        for num_Ch in range(self.num_channels):
          DataNames2.append([self.file_path+'_'+str(num_R)+'.xsp'])
        DataNames.append(DataNames2)
    elif DataFlag == 'Image':
      DataNames2 = []
      for num_Ch in range(self.num_channels):
        DataNames2.append([self.file_path+'_'+Alphabet[num_Ch]+'.xim'])
      DataNames.append(DataNames2)
    elif DataFlag in ['Multi-Region Image']:
      for num_R in range(self.num_regions):
          DataNames2 = []
          for num_Ch in range(self.num_channels):
              DataNames2.append([self.file_path + '_' + Alphabet[num_Ch] + str(num_R) + '.xim'])
          DataNames.append(DataNames2)
    elif DataFlag in ['Multi-Region Image Stack']:
        for num_Ch in range(self.num_channels):
            DataNames2 = []
            for num_R in range(self.num_regions):
                DataNames3 = []
                for num_E in range(self.data_size[0][2]):
                    DataNames3.append(
                        self.file_path + '_' + Alphabet[num_Ch] + str(num_E).zfill(3) + str(num_R) + '.xim')
                DataNames2.append(DataNames3)
            DataNames = [DataNames2]
            # print(DataNames[0])
    elif DataFlag == 'Image Stack':
      DataNames2 = []
      for num_Ch in range(self.num_channels):
        DataNames3 = []
        for num_E in range(self.data_size[0][2]):
          DataNames3.append(self.file_path+'_'+Alphabet[num_Ch]+str(num_E).zfill(3)+'.xim')
        DataNames2.append(DataNames3)
      DataNames = [DataNames2]
    else:
      print("WARNING! Unknown flag:", DataFlag)
    #for num_R in range(len(DataNames)):#Check that names correspond to existing files
      #for num_Ch in range(len(DataNames[num_R])):
        #for num_E in range(len(DataNames[num_R][num_Ch])):
          #if exists(DataNames[num_R][num_Ch][num_E]) == False:
            #print "WARNING! Data file doesn't exist:", DataNames[num_R][num_Ch][num_E]
    return DataNames

#----------------------------------------------------------------------
  def parse_DataSize(self):
    """Calculate data array size. This is useful for making sure all of the lists of data are the correct length."""
    DataSize = []
    for R_num in range(self.num_regions):
      DataSize.append([1,1,1])# [PAxis,QAxis,StackAxis] (switch to [X1,X2,E] later)]
      DataSize[R_num][0] = int(self.js['ScanDefinition']['Regions'][R_num+1]['PAxis']['Points'][0])
      if 'QAxis' in self.js['ScanDefinition']['Regions'][R_num+1] and 'Points' in self.js['ScanDefinition']['Regions'][R_num+1]['QAxis']:
        DataSize[R_num][1] = int(self.js['ScanDefinition']['Regions'][R_num+1]['QAxis']['Points'][0])
      if 'StackAxis' in self.js['ScanDefinition'] and 'Points' in self.js['ScanDefinition']['StackAxis']:
        DataSize[R_num][2] = int(self.js['ScanDefinition']['StackAxis']['Points'][0])
      if self.js['ScanDefinition']['Type'] in ['NEXAFS Point Scan','NEXAFS Line Scan']:
        DataSize[R_num] = [DataSize[R_num][1],1,DataSize[R_num][0]]#switch to [X1,X2,E] format
#        DataSize[R_num] = [1,DataSize[R_num][1],DataSize[R_num][0]]#also works, but might be problematic for finding number of spatial points
    return DataSize

#----------------------------------------------------------------------
def read(filename, self, selection=None, *args, **kwargs):
    JS = JS_FileLoader(filename)
    allowed = ['Image Stack', 'Image', 'Multi-Region Image Stack']
    if JS.js['ScanDefinition']['Flags'] in allowed:
        self.x_dist = numpy.array([float(i) for i in JS.js['ScanDefinition']['Regions'][selection[0]+1]['PAxis']['Points'][1:] ])
        self.y_dist = numpy.array([float(i) for i in JS.js['ScanDefinition']['Regions'][selection[0]+1]['QAxis']['Points'][1:] ])
        self.ev = numpy.array([float(i) for i in JS.js['ScanDefinition']['StackAxis']['Points'][1:] ])
        self.n_cols = len(self.x_dist)
        self.n_rows = len(self.y_dist)
        self.n_ev = len(self.ev)

        msec = float(JS.js['ScanDefinition']['Dwell'])
        self.data_dwell = numpy.ones((self.n_ev))*msec

        imagestack = numpy.empty((self.n_cols,self.n_rows,self.n_ev), numpy.int32)
        for i in range(len(JS.data_names[selection[1]][selection[0]])):
            try:
                imagestack[:,:,i] = numpy.loadtxt(JS.data_names[selection[1]][selection[0]][i], numpy.int32).T
            except ValueError:
                print("Aborted stack or XIMs with inconsistent dimensions.")
                imagestack[:,:,i] = numpy.nan
            except IOError:
                print("Image file not found.")
                imagestack[:,:,i] = numpy.nan
        self.absdata = numpy.empty((self.n_cols,self.n_rows, self.n_ev))

        self.absdata = numpy.reshape(imagestack, (self.n_cols,self.n_rows, self.n_ev), order='F')

        self.fill_h5_struct_from_stk()

    else:
        print("Unknown Format")

    return

#----------------------------------------------------------------------
def read_js_i0(self, filename, *args, **kwargs):
    JS = JS_FileLoader(filename)

    if 'ScanType' in JS.js['ScanDefinition'] and JS.js['ScanDefinition']['ScanType'] == 'Spectra':
        Energies = JS.js['ScanDefinition']['Regions'][1]['PAxis']['Points'][1:]
        tempimage = numpy.loadtxt(JS.data_names[0][0][0], numpy.float32)
        Data = tempimage[:,1]
    elif JS.js['ScanDefinition']['Type'] == 'NEXAFS Line Scan':
        Energies = JS.js['ScanDefinition']['Regions'][1]['PAxis']['Points'][1:]
        tempimage = numpy.loadtxt(JS.data_names[0][0][0], numpy.int32)
        Data = numpy.mean(tempimage,axis=0)
    else:# Image Stack
        Energies = JS.js['ScanDefinition']['StackAxis']['Points'][1:]
        tempimage = numpy.empty((JS.data_size[0][0],JS.data_size[0][1]), numpy.int32)
        Data = numpy.empty((JS.data_size[0][2]), numpy.int32)
        for i in range(len(JS.data_names[0][0])):
            tempimage = numpy.loadtxt(JS.data_names[0][0][i], numpy.int32)
            Data[i] = numpy.mean(tempimage)


    msec = float(JS.js['ScanDefinition']['Dwell'])#shouldn't this be needed?
    self.i0_dwell = msec
    self.evi0 = numpy.array([float(i) for i in Energies])
    self.i0data = Data
    return