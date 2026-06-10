
from PyQt5 import QtWidgets, QtCore, QtGui
import os
from ... import file_plugins

class File_GUI_Class():
    """
    Ask user to choose file and then use an appropriate plugin to read and return a data structure.
    """
    def __init__(self):
        self.last_path = dict([a,dict([t,os.getcwd()] for t in file_plugins.data_types)] for a in file_plugins.actions)
        self.last_filter = dict([a,dict([t,0] for t in file_plugins.data_types)] for a in file_plugins.actions)
        self.supported_filters = file_plugins.supported_filters
        self.filter_list = file_plugins.filter_list
        self.option_write_json = False
        self.option_norm_ringcurrent = True
        #print(self.filter_list)

    def SelectFile(self,action,data_type):
        #print(action,data_type)
        dlg=QtWidgets.QFileDialog(None)
        dlg.setWindowTitle('Choose File')
        dlg.setViewMode(QtWidgets.QFileDialog.Detail)
        #dlg.setOption(QtWidgets.QFileDialog.DontUseNativeDialog)
        if action == "write":
            dlg.setAcceptMode(QtWidgets.QFileDialog.AcceptSave)
        dlg.setDirectory(self.last_path[action][data_type])
        dlg.setNameFilters(self.filter_list[action][data_type])
        dlg.selectNameFilter(self.filter_list[action][data_type][self.last_filter[action][data_type]])
        if dlg.exec_(): #if not cancelled
            self.last_path[action][data_type] = os.path.split(str(dlg.selectedFiles()[0]))[0]
            chosen_plugin = None
            if action == 'read':
                checklist = self.filter_list[action][data_type][1:-1] #take into account the extra "Supported" and "All" filter entries
            else:
                checklist = self.filter_list[action][data_type]
            for i,filt in enumerate(checklist):
                if filt==dlg.selectedNameFilter():
                    chosen_plugin = file_plugins.supported_plugins[action][data_type][i]
                    break
            if chosen_plugin is not None:
                self.last_filter[action][data_type] = i
            return (str(dlg.selectedFiles()[0]),chosen_plugin)
        else:
            return (None,None)


#----------------------------------------------------------------------
    class DataChoiceDialog(QtWidgets.QDialog):

        def __init__(self,filepath=None,filestruct=None,plugin=None):
            # Explicitly calling QDialog's init
            super(File_GUI_Class.DataChoiceDialog, self).__init__() 
            self.filepath = filepath
            self.selection = None
            self.setWindowTitle('Choose Dataset')

            # A vertical box layout containing rows
            self.MainSizer = QtWidgets.QVBoxLayout()
            hbox1 = QtWidgets.QHBoxLayout()
            hbox2 = QtWidgets.QHBoxLayout()
            hbox3 = QtWidgets.QHBoxLayout()
            hbox4 = QtWidgets.QHBoxLayout()

            # Add the widgets to the first row
            hbox1.addWidget(QtWidgets.QLabel("Path:"))
            self.Path_text = QtWidgets.QLabel("")
            hbox1.addWidget(self.Path_text,stretch=1)
            self.MainSizer.addLayout(hbox1)

            # Add the widgets to the second row
            hbox2.addWidget(QtWidgets.QLabel('File:'))
            self.File_text = QtWidgets.QLineEdit(self)
            hbox2.addWidget(self.File_text,stretch=1)
            browse_button = QtWidgets.QPushButton('Browse...')
            browse_button.clicked.connect(self.OnBrowse)
            hbox2.addWidget(browse_button)
            self.MainSizer.addLayout(hbox2)

            # Add the widgets to the third row - dynamic set of widgets to display file info
            self.Entry_info = File_GUI_Class.EntryInfoBox()
            self.MainSizer.addWidget(self.Entry_info)

            # Add widgets for the fourth row - just OK, cancel buttons
            self.jsoncheck = QtWidgets.QCheckBox("Write Data Structure to JSON-file (for .HDR only!)")
            self.jsoncheck.stateChanged.connect(lambda: self.setChecked("jsoncheck",self.jsoncheck.isChecked()))
            hbox4.addWidget(self.jsoncheck)
            self.jsoncheck.setChecked(File_GUI.option_write_json)
            if plugin and plugin.title not in ['SDF']:
                self.jsoncheck.hide()

            hbox4.addStretch(1)
            self.button_ok = QtWidgets.QPushButton('Accept')
            self.button_ok.clicked.connect(self.OnAccept)
            self.button_ok.setEnabled(False)
            hbox4.addWidget(self.button_ok)
            button_cancel = QtWidgets.QPushButton('Cancel')
            button_cancel.clicked.connect(self.close)
            hbox4.addWidget(button_cancel)
            self.MainSizer.addLayout(hbox4)

            # Set self.MainSizer as the layout for the window
            self.setLayout(self.MainSizer)
            self.setModal(True)
            #self.show()

            if filepath is not None:
                path, filename = os.path.split(str(self.filepath))
                self.Path_text.setText(path)
                self.File_text.setText(filename)
                if filestruct is None:
                    self.contents = file_plugins.GetFileStructure(str(self.filepath))
                    #print(self.contents)
                else:
                    self.contents = filestruct
                self.Entry_info.UpdateInfo(self.contents)

        def setChecked(self,name,value=True):
            if name == "jsoncheck":
                File_GUI.option_write_json = value
                #self.jsoncheck.setChecked(value)
        def OnAccept(self):
            self.selection = self.Entry_info.GetSelection()
            if len(self.selection) == 0:
                self.close() #accepting with no regions selected is the same as clicking 'cancel'
            else:
                self.accept()

        def OnBrowse(self):
            filepath, plugin = File_GUI.SelectFile('read','stack')
            if filepath is None:
                return
            path, filename = os.path.split(str(filepath))
            self.filepath = str(filepath)
            self.Path_text.setText(path)
            self.File_text.setText(filename)
            self.contents = file_plugins.GetFileStructure(str(filepath),plugin=plugin)
            if self.contents is None:
                self.selection = [(0,0)]
                self.accept()
            else:
                self.Entry_info.UpdateInfo(self.contents)


    #-------------------------------------
    class InfoItem(QtWidgets.QHBoxLayout):
        '''Row of widgets describing the contents of a file entry that appear within the EntryInfoBox'''
        def __init__(self,name,contents,allowed_application_definitions,index):
            super(File_GUI_Class.InfoItem, self).__init__()
            self.valid_flag = self.checkValidity(contents,allowed_application_definitions)
            self.index = index
            self.populate(name,contents)

        def populate(self,name,contents):
            self.checkbox = QtWidgets.QCheckBox(name+' ('+str(contents.definition)+':'+str(contents.scan_type)+')')
            self.checkbox.clicked.connect(self.setChecked)
            self.addWidget(self.checkbox)
            self.addStretch(1)
            self.addWidget(QtWidgets.QLabel('Channels:'))
            self.channel_combobox = QtWidgets.QComboBox()
            for c in contents:
                self.channel_combobox.addItem(c)
            self.addWidget(self.channel_combobox)
            self.addStretch(1)
            Data_Size_Label = QtWidgets.QLabel('Points: '+str(contents.data_shape))
            if contents.data_axes is not None:
                Data_Size_Label.setToolTip(' | '.join(contents.data_axes))
            self.addWidget(Data_Size_Label)
            if not self.valid_flag:
                self.checkbox.setEnabled(False)
                self.channel_combobox.setEnabled(False)

            Data_Normalize_Label = QtWidgets.QLabel('Normalize data to: ')
            self.addWidget(Data_Normalize_Label)
            self.norm_combobox = QtWidgets.QComboBox()
            self.norm_combobox.addItem('none')
            if contents.norm_data is not None:
                for normkey in contents.norm_data.keys():
                    self.norm_combobox.addItem(normkey)
            # Only add fallback option for NXstxm files
            if contents.definition == 'NXstxm':
                self.norm_combobox.addItem('ring current (fallback)')
            preferences = ["control", "ringcurrent","ring current (fallback)"]
            for pref in preferences:
                index = self.norm_combobox.findText(pref)
                if index != -1:
                    self.norm_combobox.setCurrentIndex(index)
                    break

            self.addWidget(self.norm_combobox)
            self.addStretch(1)

        def setChecked(self,value=True):
            self.checkbox.setChecked(value)

        def isEnabled(self):
            return self.checkbox.isEnabled()

        def isValid(self):
            return self.valid_flag

        def checkValidity(self,contents,allowed_application_definitions):
            '''Check if data file entry appears to have the correct structure.'''
            return True
            #if contents.definition in allowed_application_definitions and contents.scan_type is not None and contents.data_shape is not None and contents.data_axes is not None:
                #return True
            #else:
                #return False

        def GetStatus(self):
            return (self.checkbox.isChecked(),self.channel_combobox.currentIndex(),self.norm_combobox.currentIndex())

    #---------------------------------------
    class EntryInfoBox(QtWidgets.QGroupBox):
        '''Widgets giving a summary of the contents of a file via an InfoItem object per file entry.'''
        def __init__(self):
            super(File_GUI_Class.EntryInfoBox, self).__init__('File Summary')
            self.vbox = QtWidgets.QVBoxLayout()
            self.setLayout(self.vbox)
            self.valid_entry_flag = False
            self.selection = (0,0)

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
            try:
                # Assume FileContents works like a dict or list of entries
                # Adjust iteration based on real FileContents structure
                for i,entry in enumerate(FileContents):
                    entryCheckBox = File_GUI_Class.InfoItem(entry,FileContents[entry],['NXstxm'],i)
                    self.vbox.addLayout(entryCheckBox)
                    if self.valid_entry_flag is False and entryCheckBox.isValid() is True:
                        entryCheckBox.setChecked(True)
                        self.valid_entry_flag = True
                        self.parent().button_ok.setEnabled(True)
            except:
                pass

        def GetSelection(self):
            selection = []
            for i in range(self.vbox.count()):
                status = self.vbox.itemAt(i).GetStatus()
                if status[0]:
                    selection.append((i,status[1],status[2])) # (region,detector,normalize)
            return selection

File_GUI = File_GUI_Class()
