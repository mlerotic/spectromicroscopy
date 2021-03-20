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
from os.path import splitext, join, dirname, isfile
from collections import OrderedDict


title = 'SDF'
extension = ['*.hdr']
read_types = ['spectrum','image','stack']
write_types = []
def identify(filename):
    try:
        if isfile(splitext(filename)[0] + '.json'):
            print("JSON-file found. No need to fetch again.")
            return False
        else:
            HDR = HDR_FileParser(filename)
            return HDR.num_regions > 0  # return true if file contains at least one region
    except:
        return False

def GetFileStructure(FileName):
    HDR = HDR_FileParser(FileName)
    if HDR.num_regions<2 and HDR.num_channels<2 and isfile(splitext(FileName)[0] + '.json'): #if json file exists, skip datachoicedialog
        return None # exit if only one choice
    D = OrderedDict()
    for i,R in enumerate(['Region_'+str(r) for r in range(HDR.num_regions)]):
        D[R] = OrderedDict()
        D[R].definition = 'SDF'
        D[R].scan_type = HDR.hdr['ScanDefinition']['Type']
        D[R].data_shape = HDR.data_size[i]
        D[R].data_axes = [HDR.hdr['ScanDefinition']['Regions'][1]['PAxis']['Name'],HDR.hdr['ScanDefinition']['Regions'][1]['QAxis']['Name'],HDR.hdr['ScanDefinition']['StackAxis']['Name']]
        for ch in range(1,HDR.num_channels+1):
            D[R][HDR.hdr['Channels'][ch]['Name']] = OrderedDict()
    return D



#----------------------------------------------------------------------
class HDR_FileParser:
  """Parse .hdr file for metadata."""
  hdr = []
  f = []
  def __init__(self, fileName,identify=False):
    if not HDR_FileParser.hdr or HDR_FileParser.f != fileName: # prevent class from fetching and parsing the *.hdr file multiple times. Use the class attribute "hdr" instead if available and check if a new file (not "f") is loaded.
        # compile some regular expressions
        self.MatchReStruct = re.compile('[\s\{\}\(\)=";]')
        self.MatchReArray = re.compile('[,\s\{\(\);]')
        self.__file = open(fileName)
        # Parse the HDR file
        HDR_FileParser.hdr = self.parseStructure()
        HDR_FileParser.f = fileName
        self.__file.close()
    else:
        None
    self.hdr = HDR_FileParser.hdr
    if 'ScanDefinition' in self.hdr:
      self.num_regions = int(self.hdr['ScanDefinition']['Regions'][0])
      self.num_channels = int(self.hdr['Channels'][0])
      self.file_path, self.file_ext = splitext(fileName)
      self.data_size = self.parse_DataSize()
      self.data_names = self.parseDataNames()
    else:
      self.num_regions = 0
      self.num_channels = 0



#----------------------------------------------------------------------
  def parseStructure(self):
    """.hdr files consist of structures and arrays. This routine sorts through the structure parts."""
    Structure = {}
    BuildWord=''
    BeforeEq=True
    QuotedWord=False
    raw = self.__file.read(1)
    while len(raw) > 0:#until we reach the end of the file
      matched = self.MatchReStruct.match(raw)
      if matched == None:
        BuildWord+=raw
      elif matched.group() == '"':
        QuotedWord= not QuotedWord
      elif QuotedWord==True:
        BuildWord+=raw
      elif matched.group() == '=':
        FieldName=BuildWord
        BuildWord=''
        BeforeEq=False
      elif matched.group() == ';':
        try:  # convert numbers into ints or floats or end up with a string
            Structure[FieldName] = int(BuildWord)
        except ValueError:
            try:
                Structure[FieldName] = float(BuildWord)
            except ValueError:
                Structure[FieldName] = BuildWord
        except TypeError:
            Structure[FieldName] = BuildWord
        BuildWord=''
        BeforeEq=True
      elif matched.group() == '{':
        #Must be after =
        BuildWord = self.parseStructure()
      elif matched.group() == '}':
        #break loop and return dictionary
        break
      elif matched.group() == '(':
        #Must be after =
        BuildWord= self.parseArray()
      elif matched.group() == ')':
        #This should not happen
        print(') in structure')
      raw = self.__file.read(1)
    return Structure

#----------------------------------------------------------------------
  def parseArray(self):
    """.hdr files consist of structures and arrays. This rountine sorts through the array parts."""
    Array = []
    BuildWord=''
    raw = self.__file.read(1)
    while len(raw) > 0:#until we reach the end of the file
      matched = self.MatchReArray.match(raw)
      if matched == None:
        BuildWord+=raw
      elif matched.group() == ',':
        if len(BuildWord) > 0:
            try:    #convert numbers in array into ints or floats
                Array.append(int(BuildWord))
            except ValueError:
                Array.append(float(BuildWord))
            BuildWord=''
      elif matched.group() == ';':
        print('; in array')
      elif matched.group() == '{':
        Array.append(self.parseStructure())
      elif matched.group() == '(':
        Array.append(self.parseArray())
      elif matched.group() == ')':
        if len(BuildWord) > 0:
            try:  #convert numbers in array into ints or floats
                Array.append(int(BuildWord))
            except ValueError:
                try:
                    Array.append(float(BuildWord))
                except ValueError: # There is an error in the HDF5toSDF conversion script at the SLS/PolLux which appends arrays with a superfluous "}". Obviously, this char cannot be converted to int or float. We therefore just skip it here.
                    pass
        break
      raw = self.__file.read(1)
    return Array

#----------------------------------------------------------------------
  def parseDataNames(self): #ToDo: Assigning Scan types by DataFlags is prone to errors. Better detect automatically via (HDR.hdr['ScanDefinition']['Regions'][0]) etc.
    """Figure out names for the .xsp or .xim files that contain the actual data, then check that the files actually exist, printing warnings if they don't."""
    DataNames = []#Regions
    DataNames2 = []#Channels
    DataNames3 = []#Energies
    Alphabet = 'abcdefghijklmnopqrstuvwxyz'
    DataFlag = self.hdr['ScanDefinition']['Flags']
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
          DataNames2.append([self.file_path+'_'+Alphabet[num_Ch]+str(num_R)+'.xim'])
        DataNames.append(DataNames2)
    elif DataFlag in ['Multi-Region Image Stack']:
      for num_Ch in range(self.num_channels):
        DataNames2 = []
        for num_R in range(self.num_regions):
            DataNames3 = []
            for num_E in range(self.data_size[0][2]):
                DataNames3.append(self.file_path + '_' + Alphabet[num_Ch] + str(num_E).zfill(3)+str(num_R)+ '.xim')
            DataNames2.append(DataNames3)
        DataNames = [DataNames2]
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
      DataSize[R_num][0] = int(self.hdr['ScanDefinition']['Regions'][R_num+1]['PAxis']['Points'][0])
      if 'QAxis' in self.hdr['ScanDefinition']['Regions'][R_num+1] and 'Points' in self.hdr['ScanDefinition']['Regions'][R_num+1]['QAxis']:
        DataSize[R_num][1] = int(self.hdr['ScanDefinition']['Regions'][R_num+1]['QAxis']['Points'][0])
      if 'StackAxis' in self.hdr['ScanDefinition'] and 'Points' in self.hdr['ScanDefinition']['StackAxis']:
        DataSize[R_num][2] = int(self.hdr['ScanDefinition']['StackAxis']['Points'][0])
      if self.hdr['ScanDefinition']['Type'] in ['NEXAFS Point Scan','NEXAFS Line Scan']:
        DataSize[R_num] = [DataSize[R_num][1],1,DataSize[R_num][0]]#switch to [X1,X2,E] format
#        DataSize[R_num] = [1,DataSize[R_num][1],DataSize[R_num][0]]#also works, but might be problematic for finding number of spatial points
    return DataSize




#----------------------------------------------------------------------
#----------------------------------------------------------------------
def read(filename, self, selection=None, JSONstatus=None):
    HDR = HDR_FileParser(filename)
    if JSONstatus:
        with open(splitext(filename)[0] + '.json', 'w') as outfile:
            json.dump(HDR.hdr, outfile, indent=4, sort_keys=True, ensure_ascii=True)
            print("JSON-file written at "+ splitext(filename)[0] + '.json')
    allowed =['Image Stack','Image','Multi-Region Image Stack']
    if HDR.hdr['ScanDefinition']['Flags'] in allowed:

        self.x_dist = numpy.array([float(i) for i in HDR.hdr['ScanDefinition']['Regions'][selection[0]+1]['PAxis']['Points'][1:] ])
        self.y_dist = numpy.array([float(i) for i in HDR.hdr['ScanDefinition']['Regions'][selection[0]+1]['QAxis']['Points'][1:] ])
        self.ev = numpy.array([float(i) for i in HDR.hdr['ScanDefinition']['StackAxis']['Points'][1:] ])

        self.n_cols = len(self.x_dist)
        self.n_rows = len(self.y_dist)
        self.n_ev = len(self.ev)

        msec = float(HDR.hdr['ScanDefinition']['Dwell'])
        self.data_dwell = numpy.ones((self.n_ev))*msec

        imagestack = numpy.empty((self.n_cols,self.n_rows,self.n_ev), numpy.int32)
        for i in range(len(HDR.data_names[selection[1]][selection[0]])):
            try:
                imagestack[:,:,i] = numpy.loadtxt(HDR.data_names[selection[1]][selection[0]][i], numpy.int32).T
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
def read_sdf_i0(self, filename):
    HDR = HDR_FileParser(filename)

    if 'ScanType' in HDR.hdr['ScanDefinition'] and HDR.hdr['ScanDefinition']['ScanType'] == 'Spectra':
        Energies = HDR.hdr['ScanDefinition']['Regions'][1]['PAxis']['Points'][1:]
        tempimage = numpy.loadtxt(HDR.data_names[0][0][0], numpy.float32)
        Data = tempimage[:,1]
    elif HDR.hdr['ScanDefinition']['Type'] == 'NEXAFS Line Scan':
        Energies = HDR.hdr['ScanDefinition']['Regions'][1]['PAxis']['Points'][1:]
        tempimage = numpy.loadtxt(HDR.data_names[0][0][0], numpy.int32)
        Data = numpy.mean(tempimage,axis=0)
    else:# Image Stack
        Energies = HDR.hdr['ScanDefinition']['StackAxis']['Points'][1:]
        tempimage = numpy.empty((HDR.data_size[0][0],HDR.data_size[0][1]), numpy.int32)
        Data = numpy.empty((HDR.data_size[0][2]), numpy.int32)
        for i in range(len(HDR.data_names[0][0])):
            tempimage = numpy.loadtxt(HDR.data_names[0][0][i], numpy.int32)
            Data[i] = numpy.mean(tempimage)

    msec = float(HDR.hdr['ScanDefinition']['Dwell'])#shouldn't this be needed?
    self.i0_dwell = msec
    self.evi0 = numpy.array([float(i) for i in Energies])
    self.i0data = Data
    return