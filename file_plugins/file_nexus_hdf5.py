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

def read(stack_object,FileName,entry=0,channel=0):
    D = GetFileInfo(FileName)
    first_entry = next(iter(D))
    first_detector = next(iter(D[first_entry]))
    #print D
    F = h5py.File(FileName, 'r')
    stack_object.ev = numpy.array(F[first_entry][first_detector]['energy'])
    stack_object.x_dist = numpy.array(F[first_entry][first_detector]['sample_x'])
    stack_object.y_dist = numpy.array(F[first_entry][first_detector]['sample_y'])
    stack_object.data_dwell = numpy.array(F[first_entry][first_detector]['count_time'])
    stack_object.n_cols = len(stack_object.x_dist)
    stack_object.n_rows = len(stack_object.y_dist)
    stack_object.n_ev = len(stack_object.ev)
    if 'axes' in F[first_entry][first_detector].attrs: # Specification correct
        axes_list = list(F[first_entry][first_detector].attrs['axes'])
        axes_order = [axes_list.index('sample_x'),axes_list.index('sample_y'),axes_list.index('energy')]
    else: # Old version from before the specification was finalised
        axes_order = [F[first_entry][first_detector]['sample_x'].attrs['axis']-1,F[first_entry][first_detector]['sample_y'].attrs['axis']-1,F[first_entry][first_detector]['energy'].attrs['axis']-1]
    stack_object.absdata = numpy.transpose(numpy.array(F[first_entry][first_detector]['data']),axes=axes_order)
    
    
    F.close()


def GetFileInfo(FileName):
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

class NeXusLoadDialog(QtGui.QWidget):

    def __init__(self):
        super(NeXusLoadDialog, self).__init__()
        self.path = None
        self.filename = None
        self.h5data = h5data()
        self.initUI()

    def initUI(self):
        #self.setGeometry(300, 300, 290, 150)
        self.setWindowTitle('Load NeXus HDF5')
        
        # A vertical box layout containing rows
        self.MainSizer = QtGui.QVBoxLayout()
        hbox1 = QtGui.QHBoxLayout()
        hbox2 = QtGui.QHBoxLayout()
        hbox3 = QtGui.QHBoxLayout()
        hbox4 = QtGui.QHBoxLayout()
        
        # Add the widgets to the first row
        hbox1.addWidget(QtGui.QLabel("Path:"))
        if self.path is None:
            self.Path_text = QtGui.QLabel("")
        else:
            self.Path_text = QtGui.QLabel(self.path)
        hbox1.addWidget(self.Path_text,stretch=1)
        self.MainSizer.addLayout(hbox1)
        
        # Add the widgets to the second row
        hbox2.addWidget(QtGui.QLabel('File:'))
        self.File_text = QtGui.QLineEdit(self)
        hbox2.addWidget(self.File_text,stretch=1)
        browse_button = QtGui.QPushButton('Browse...')
        browse_button.clicked.connect(self.OnBrowse)
        hbox2.addWidget(browse_button)
        self.MainSizer.addLayout(hbox2)
        
        # Add the widgets to the third row - dynamic set of widgets to display file info
        self.Entry_info = EntryInfoBox()
        self.MainSizer.addWidget(self.Entry_info)
        
        # Add widgets for the fourth row - just OK, cancel buttons
        hbox4.addStretch(1)
        self.button_ok = QtGui.QPushButton('Accept')
        self.button_ok.clicked.connect(self.OnAccept)
        self.button_ok.setEnabled(False)
        hbox4.addWidget(self.button_ok)
        button_cancel = QtGui.QPushButton('Cancel')
        button_cancel.clicked.connect(self.close)
        hbox4.addWidget(button_cancel)
        self.MainSizer.addLayout(hbox4)
        
        # Set self.MainSizer as the layout for the window
        self.setLayout(self.MainSizer)
        self.show()

    def OnAccept(self):
        #print self.h5data.check_h5_format(os.path.join(self.path,self.filename))
        return

    def OnBrowse(self):
        if self.path is None:
            start_path = ''
        else:
            start_path = self.path
        FileChoice = QtGui.QFileDialog.getOpenFileName(self, "Choose a file", start_path, filter='HDF (*.hdf5);;*.*')
        if FileChoice == '':
            return
        self.path, self.filename = os.path.split(str(FileChoice))
        
        self.Path_text.setText(self.path)
        self.File_text.setText(self.filename)
        self.contents = self.h5data.GetFileInfo(str(FileChoice))
        self.Entry_info.UpdateInfo(self.contents)
        

class InfoItem(QtGui.QHBoxLayout):
    '''Row of widgets describing the contents of a file entry that appear within the EntryInfoBox'''
    def __init__(self,name,contents,allowed_application_definitions):
        super(InfoItem, self).__init__()
        self.valid_flag = self.checkValidity(contents,allowed_application_definitions)
        self.populate(name,contents)
        
    def populate(self,name,contents):
        self.radioButton = QtGui.QRadioButton(name+' ('+str(contents.definition)+':'+str(contents.scan_type)+')')
        self.addWidget(self.radioButton)
        self.addStretch(1)
        self.addWidget(QtGui.QLabel('Channels:'))
        channel_combobox = QtGui.QComboBox()
        for c in contents:
            channel_combobox.addItem(c)
        self.addWidget(channel_combobox)
        self.addStretch(1)
        Data_Size_Label = QtGui.QLabel('Points: '+str(contents.data_shape))
        if contents.data_axes is not None:
            Data_Size_Label.setToolTip(' | '.join(contents.data_axes))
        self.addWidget(Data_Size_Label)
        if not self.valid_flag:
            self.radioButton.setEnabled(False)
            channel_combobox.setEnabled(False)
    
    def setChecked(self,value=True):
        self.radioButton.setChecked(value)

    def isEnabled(self):
        return self.radioButton.isEnabled()

    def isValid(self):
        return self.valid_flag

    def checkValidity(self,contents,allowed_application_definitions):
        '''Check if data file entry appears to have the correct structure.'''
        if contents.definition in allowed_application_definitions and contents.scan_type is not None and contents.data_shape is not None and contents.data_axes is not None:
            return True
        else:
            return False

class EntryInfoBox(QtGui.QGroupBox):
    '''Widgets giving a summary of the contents of a file via an InfoItem object per file entry.'''
    def __init__(self):
        super(EntryInfoBox, self).__init__('File Summary')
        self.vbox = QtGui.QVBoxLayout()
        self.setLayout(self.vbox)
        self.valid_entry_flag = False
        
    def ClearAll(self):
        # Cycle through children and mark for deletion
        while self.vbox.count():
            row = self.vbox.takeAt(0)
            while row.count():
                item = row.takeAt(0)
                if item.widget() is not None:
                    item.widget().deleteLater()
        self.valid_entry_flag = False
        
    def UpdateInfo(self, FileContents):
        self.ClearAll()
        # Now popoulate widgets representing file contents
        for entry in FileContents:
            entryCheckBox = InfoItem(entry,FileContents[entry],self.parent().h5data.allowed_application_definitions)
            self.vbox.addLayout(entryCheckBox)
            if self.valid_entry_flag is False and entryCheckBox.isValid() is True:
                entryCheckBox.setChecked(True)
                self.valid_entry_flag = True
                self.parent().button_ok.setEnabled(True)
        

