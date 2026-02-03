
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtCore import Qt
import os
import re
import numpy as np
from ...file_plugins import file_xrm, file_bim, file_tif # file_sm_netcdf imported conditionally

class StackListFrame(QtWidgets.QDialog):

    def __init__(self, parent, filepath, com, stack, data_struct):
        QtWidgets.QWidget.__init__(self, parent)

        self.parent = parent

        self.data_struct = data_struct
        self.stk = stack
        self.common = com

        self.resize(600, 500)
        self.setWindowTitle('Stack File List')

        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)

        self.filepath = str(filepath)

        self.have1st = 0
        self.havelast = 0
        self.file1st = ' '
        self.filelast = ' '

        self.filetype = ''

        vbox = QtWidgets.QVBoxLayout()

        self.textt = QtWidgets.QLabel(self)
        self.textt.setText('Select first stack file')

        vbox.addStretch(1)
        vbox.addWidget(self.textt)


        self.filelist = QtWidgets.QTableWidget()
        self.filelist.setMinimumHeight(450)
        self.filelist.setColumnCount(4)
        self.filelist.setHorizontalHeaderLabels(('File list', 'X', 'Y', 'eV'))
        self.filelist.setShowGrid(False)
        self.filelist.verticalHeader().setVisible(False)

        self.filelist.setColumnWidth(0,400)
        self.filelist.setColumnWidth(1,50)
        self.filelist.setColumnWidth(2,50)
        self.filelist.setColumnWidth(3,50)

        self.filelist.setRowCount(0)

        self.filelist.cellClicked.connect(self.OnFileList)
        self.filelist.horizontalHeader().setSortIndicatorShown(True)
        self.filelist.horizontalHeader().sectionClicked.connect(self.OnSort)

        vbox.addWidget(self.filelist)
        vbox.addStretch(1)

        self.tc_first = QtWidgets.QLabel(self)
        self.tc_first.setText('First stack file: ')
        self.tc_last = QtWidgets.QLabel(self)
        self.tc_last.setText('Last stack file: ')


        vbox.addWidget(self.tc_first)
        vbox.addWidget(self.tc_last)
        vbox.addStretch(1)


        hbox = QtWidgets.QHBoxLayout()


        self.button_accept = QtWidgets.QPushButton('Accept')
        self.button_accept.setEnabled(False)
        self.button_accept.clicked.connect( self.OnAccept)
        hbox.addWidget(self.button_accept)

        button_cancel = QtWidgets.QPushButton('Cancel')
        button_cancel.clicked.connect( self.close)
        hbox.addWidget(button_cancel)

        vbox.addLayout(hbox)

        vbox.addStretch(1)


        self.setLayout(vbox)

        self.ShowFileList()

# ----------------------------------------------------------------------
    def OnSort(self):
        self.sm_files = [self.filelist.item(i,0).text() for i in range(self.filelist.rowCount())]
        self.have1st = 0
        self.havelast = 0
        self.file1st = ' '
        self.filelast = ' '
        self.button_accept.setEnabled(False)
        self.tc_first.setText('First stack file: ')
        self.tc_last.setText('Last stack file: ')
#----------------------------------------------------------------------
    def OnFileList(self, row, column):


        if (self.have1st == 1) and (self.havelast==1):
            self.have1st = 0
            self.havelast = 0
            self.button_accept.setEnabled(False)
            self.tc_first.setText('First stack file: ')
            self.tc_last.setText('Last stack file: ')

        item = self.filelist.item(row, 0)
        fn =  item.text()
        if self.have1st == 0:
            self.tc_first.setText('First stack file: ' + fn)
            self.textt.setText('Select last stack file')
            self.file1st = fn
            self.have1st = 1
        elif self.havelast == 0:
            self.tc_last.setText('Last stack file: ' + fn)
            self.textt.setText('Select first stack file')
            self.filelast = fn
            self.button_accept.setEnabled(True)
            self.havelast = 1
#----------------------------------------------------------------------
    def ShowFileList(self):

        filepath = str(self.filepath)
        self.sm_files = [x for x in os.listdir(filepath) if x.endswith('.sm')]

        if self.sm_files:

            self.filetype = 'sm'

            try:
                from netCDF4 import Dataset
                from ...file_plugins import file_sm_netcdf

            except:
                QtWidgets.QMessageBox.warning(self, 'Error', "Could not import netCDF4 library.")
                return

            count = 0

            for i in range(len(self.sm_files)):

                filename = self.sm_files[i]
                thisfile = os.path.join(filepath, filename)

                filever, ncols, nrows, iev = file_sm_netcdf.read_sm_header(thisfile)

                if filever > 0:
                    self.filelist.insertRow(count)
                    self.filelist.setRowHeight(count,20)

                    self.filelist.setItem(count, 0, QtWidgets.QTableWidgetItem(filename))
                    self.filelist.setItem(count, 1, QtWidgets.QTableWidgetItem(str(ncols)))
                    self.filelist.setItem(count, 2, QtWidgets.QTableWidgetItem(str(nrows)))
                    self.filelist.setItem(count, 3, QtWidgets.QTableWidgetItem('{0:5.2f}'.format(iev)))

                    count += 1
                else:
                    continue




        self.xrm_files = [x for x in os.listdir(filepath) if x.endswith('.xrm')]


        if self.xrm_files:

            self.filetype = 'xrm'

            count = 0

            for i in range(len(self.xrm_files)):

                filename = self.xrm_files[i]
                thisfile = os.path.join(filepath, filename)

                ncols, nrows, iev = file_xrm.read_xrm_fileinfo(thisfile)

                if ncols > 0:
                    self.filelist.insertRow(count)
                    self.filelist.setRowHeight(count,20)

                    self.filelist.setItem(count, 0, QtWidgets.QTableWidgetItem(filename))
                    self.filelist.setItem(count, 1, QtWidgets.QTableWidgetItem(str(ncols)))
                    self.filelist.setItem(count, 2, QtWidgets.QTableWidgetItem(str(nrows)))
                    self.filelist.setItem(count, 3, QtWidgets.QTableWidgetItem('{0:5.2f}'.format(iev)))

                    count += 1


            self.sm_files = self.xrm_files


        self.bim_files = [x for x in os.listdir(filepath) if x.endswith('.bim')]


        if self.bim_files:

            self.filetype = 'bim'



            count = 0

            for i in range(len(self.bim_files)):
                #print sm_files
                filename = self.bim_files[i]
                thisfile = os.path.join(filepath, filename)

                ncols, nrows, iev = file_bim.read_bim_info(thisfile)

                if ncols >0 :
                    self.filelist.insertRow(count)
                    self.filelist.setRowHeight(count,20)

                    self.filelist.setItem(count, 0, QtWidgets.QTableWidgetItem(filename))
                    self.filelist.setItem(count, 1, QtWidgets.QTableWidgetItem(str(ncols)))
                    self.filelist.setItem(count, 2, QtWidgets.QTableWidgetItem(str(nrows)))
                    self.filelist.setItem(count, 3, QtWidgets.QTableWidgetItem('{0:5.2f}'.format(iev)))

                    count += 1


            self.sm_files = self.bim_files

        self.tif_files = [x for x in os.listdir(filepath) if x.endswith('.tif')]
        if self.tif_files:

            self.filetype = 'tif'

            for i in range(len(self.tif_files)):

                filename = self.tif_files[i]
                thisfile = os.path.join(filepath, filename)

                ncols, nrows = file_tif.read_tif_info(thisfile)
                #auto-read energies when the following syntax applies "<str>_XXX.XeV_XX.tif"
                fnlist = filename.split('_')
                try:
                    ind =[ m for m, j in enumerate(fnlist) if re.search('\deV', j)][0]
                    iev = float(fnlist[ind][:-2])
                except IndexError:
                    iev = i
                self.filelist.insertRow(i)
                self.filelist.setRowHeight(i, 20)

                self.filelist.setItem(i, 0, QtWidgets.QTableWidgetItem(filename))
                self.filelist.setItem(i, 1, QtWidgets.QTableWidgetItem(str(ncols)))
                self.filelist.setItem(i, 2, QtWidgets.QTableWidgetItem(str(nrows)))
                self.filelist.setItem(i, 3, QtWidgets.QTableWidgetItem('{0:5.2f}'.format(iev)))

            self.sm_files = self.tif_files

        self.filelist.setSortingEnabled(True)
        return


#----------------------------------------------------------------------
    def OnAccept(self, evt):

        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))

        if self.common.stack_loaded == 1:
            self.parent.new_stack_refresh()
            self.stk.new_data()

        ind1st = self.sm_files.index(self.file1st)
        indlast = self.sm_files.index(self.filelast)
        if indlast < ind1st:
            filelist = self.sm_files[indlast:ind1st+1]
            filelist.reverse()
        else:
            filelist = self.sm_files[ind1st:indlast+1]


        if self.filetype == 'sm':
            from ...file_plugins import file_sm_netcdf
            file_sm_netcdf.read_sm_list(self, filelist, self.filepath, self.data_struct)
        elif self.filetype == 'xrm':
            file_xrm.read_xrm_list(self, filelist, self.filepath, self.data_struct)
        elif self.filetype == 'bim':
            file_bim.read_bim_list(self, filelist, self.filepath, self.data_struct)
        elif self.filetype == 'tif':
            file_tif.read_tif_list(self, filelist, self.filepath, self.data_struct)
        else:
            print('Wrong file type')
            return


        #fill the gui structure data
        self.stk.absdata = self.data_struct.exchange.data

        datadim = np.int32(self.stk.absdata.shape)
        self.stk.n_cols = datadim[0].copy()
        self.stk.n_rows =  datadim[1].copy()
        self.stk.ev = self.data_struct.exchange.energy
        self.stk.n_ev = np.int32(self.stk.ev.shape[0]).copy()

        npixels = self.stk.n_cols*self.stk.n_rows*self.stk.n_ev


        self.stk.x_dist = self.data_struct.exchange.x
        self.stk.y_dist = self.data_struct.exchange.y

        self.stk.data_dwell = self.data_struct.spectromicroscopy.data_dwell


        self.stk.fill_h5_struct_from_stk()

        self.stk.setScale()


        #self.parent.page1.iev = int(self.stk.n_ev/3) #Is this correct?

        self.parent.ix = int(self.stk.n_cols/2)
        self.parent.iy = int(self.stk.n_rows/2)

        self.common.stack_loaded = 1

        self.parent.page0.absimgfig.loadNewImage()
        directory = os.path.dirname(str(self.filepath))
        self.parent.ShowInfo(os.path.basename(str(self.filepath)), directory)
        # self.page1.ResetDisplaySettings()
        self.parent.page1.absimgfig.loadNewImageWithROI()
        self.parent.page1.button_multicrop.setText('Crop stack 3D...')
        # print(x,y), (self.ix,self.iy), self.stk.absdata.shape
        self.parent.page1.specfig.ClearandReload()

        # self.parent.refresh_widgets()
        # self.parent.page1.ResetDisplaySettings()
        # self.parent.page1.filename = filelist[0]
        # # self.parent.page1.textctrl.setText(filelist[0])
        #
        # self.parent.page0.slider_eng.setRange(0,self.stk.n_ev-1)
        # #self.parent.page0.iev = int(self.stk.n_ev/2)
        # self.parent.page0.slider_eng.setValue(self.parent.page1.iev)
        #
        # self.parent.page1.slider_eng.setRange(0,self.stk.n_ev-1)
        # #self.parent.page1.iev = self.stk.n_ev/2
        # self.parent.page1.slider_eng.setValue(self.parent.page1.iev)
        #
        # self.parent.page1.specfig.loadNewSpectrum()
        # self.parent.page1.absimgfig.loadNewImageWithROI()
        #
        #
        # self.parent.ShowInfo(filelist[0], self.filepath)
        self.parent.page5.updatewidgets()
        self.parent.refresh_widgets()
        QtWidgets.QApplication.restoreOverrideCursor()
        self.close()
