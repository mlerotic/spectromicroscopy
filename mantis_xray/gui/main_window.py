import os
import numpy as np
import sys
import getopt
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtCore import Qt
from ..analysis import data_struct
from ..analysis import data_stack
from ..analysis import analyze
from ..analysis import nnma
from ..core.constants import Winsizex, Winsizey
from .. import __version__ as version
from ..core.state import common
from ..helpers import resource_path
from ..file_plugins import file_dataexch_hdf5
from .. import file_plugins
from .pages import (
    PageLoadData, PageStack, PagePCACluster,
    PageSpectral, PagePeakID, PageXrayPeakFitting,
    PageNNMA, PageTomo
)
from .utils.file_gui import File_GUI
from .dialogs import StackListFrame, AboutFrame

class MainFrame(QtWidgets.QMainWindow):

    def __init__(self):
        super(MainFrame, self).__init__()
        self.active_widget = None
        self.initUI()



#----------------------------------------------------------------------
    def initUI(self):

        self.data_struct = data_struct.h5()
        self.stk = data_stack.data(self.data_struct)
        self.stk.ask_fix_callback = self.ask_fix_duplicates_callback
        self.anlz = analyze.analyze(self.stk)
        self.nnma = nnma.nnma(self.stk)
        self.common = common()


        self.resize(Winsizex, Winsizey)
        self.setWindowTitle('Mantis v.{0}'.format(version))

        self.initToolbar()

        self.controlc = QtWidgets.QShortcut(QtGui.QKeySequence("Ctrl+C"), self)
        self.controlc.activated.connect(self.OnCopyWidget)

        ico = QtGui.QIcon(resource_path(os.path.join('images','logo-2l-32.ico')))
        self.setWindowIcon(ico)

        self.tabs = QtWidgets.QTabWidget()
        self.tab_tomo = PageTomo(self.common, self.data_struct, self.stk, self.anlz)

        # create the tab windows with descriptive attribute names
        self.tab_load = PageLoadData(self.common, self.data_struct, self.stk)
        self.tab_prep = PageStack(self.common, self.data_struct, self.stk)
        self.tab_clus = PagePCACluster(self.common, self.data_struct, self.stk, self.anlz)
        self.tab_spec = PageSpectral(self.common, self.data_struct, self.stk, self.anlz)
        self.tab_peak = PagePeakID(self.common, self.data_struct, self.stk, self.anlz)
        self.tab_pfit = PageXrayPeakFitting(self.common, self.data_struct, self.stk, self.anlz)
        self.tab_nnma = PageNNMA(self.common, self.data_struct, self.stk, self.anlz, self.nnma)

        self.tabs.addTab(self.tab_load, "Load Data")
        self.tabs.addTab(self.tab_prep, "Preprocess Data")
        self.tabs.addTab(self.tab_clus, "PCA && Cluster Analysis")
        self.tabs.addTab(self.tab_spec, "Spectral Maps")
        self.tabs.addTab(self.tab_nnma, "NNMA Analysis")
        self.tabs.addTab(self.tab_peak, "Peak ID")
        self.tabs.addTab(self.tab_pfit, "XrayPeakFitting")

        self.tabs.addTab(self.tab_tomo, "Tomography")
        self.tabs.setTabEnabled(self.tabs.indexOf(self.tab_tomo), True)
        self._set_tomo_tab_style(True)

#        if showmaptab:
#            self.tab_map = PageMap(self.common, self.data_struct, self.stk)
#            self.tabs.addTab(self.tab_map, "Image Maps")

        if sys.platform == 'win32':
            self.tabs.setMinimumHeight(750)
        else:
            self.tabs.setMinimumHeight(400)


        # print Qt colours
        BackgroundColour = self.palette().color(self.palette().Background)
        DarkBackgroundFlag = (BackgroundColour.red()+BackgroundColour.green()+BackgroundColour.blue()) < 375
        if DarkBackgroundFlag:
            self.tabs.tabBar().setTabTextColor(0, QtGui.QColor('lightgreen'))
            self.tabs.tabBar().setTabTextColor(1, QtGui.QColor('lightgreen'))
            self.tabs.tabBar().setTabTextColor(2, QtGui.QColor('tomato'))
            self.tabs.tabBar().setTabTextColor(3, QtGui.QColor('tomato'))
            self.tabs.tabBar().setTabTextColor(4, QtGui.QColor('orchid'))
            self.tabs.tabBar().setTabTextColor(5, QtGui.QColor('orchid'))
        else:
            self.tabs.tabBar().setTabTextColor(0, QtGui.QColor('green'))
            self.tabs.tabBar().setTabTextColor(1, QtGui.QColor('green'))
            self.tabs.tabBar().setTabTextColor(2, QtGui.QColor('darkRed'))
            self.tabs.tabBar().setTabTextColor(3, QtGui.QColor('darkRed'))
            self.tabs.tabBar().setTabTextColor(4, QtGui.QColor('purple'))
            self.tabs.tabBar().setTabTextColor(5, QtGui.QColor('purple'))
        tomo_index = self.tabs.indexOf(self.tab_tomo)
        self._set_tomo_tab_style(tomo_index != -1 and self.tabs.isTabEnabled(tomo_index))
        #if showmaptab:
        #    self.tabs.tabBar().setTabTextColor(9, QtGui.QColor('darkblue'))



        # Only add "expert" pages if option "--key" is given in command line
        try:
            options, extraParams = getopt.getopt(sys.argv[1:], '', ['wx', 'batch', 'nnma', 'ica', 'keyeng'])
        except:
            print('Error - wrong command line option used. Available options are --wx, --batch and --nnma')
            return

#         for opt, arg in options:
#             if opt in '--nnma':
#                 if verbose: print "Running with NNMA."
#                 self.tab_nnma = PageNNMA(self.common, self.data_struct, self.stk, self.anlz, self.nnma)
#                 self.tabs.addTab(self.tab_nnma, "NNMA Analysis")

        layout = QtWidgets.QVBoxLayout()

        layout.addWidget(self.tabs)
        #self.setCentralWidget(self.tabs)
        self.ShowInfo('File name', os.getcwd())
        self.tab_load.FnOverlayCheckBox.stateChanged.connect(self.syncFnOverlayCheckBox)
        self.tab_prep.FnOverlayCheckBox.stateChanged.connect(self.syncFnOverlayCheckBox)
        self.tab_clus.FnOverlayCheckBox.stateChanged.connect(self.syncFnOverlayCheckBox)

        self.scrollArea = QtWidgets.QScrollArea()
        #self.scrollArea.setBackgroundRole(QtGui.QPalette.Dark)
        self.scrollArea.setWidget(self.tabs)
        #self.scrollArea.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        #self.scrollArea.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.scrollArea.setWidgetResizable(True)
        self.setCentralWidget(self.scrollArea)

        screen = QtWidgets.QDesktopWidget().screenGeometry()
        center = screen.center()

        winrect = self.frameGeometry()
        winrect.moveCenter(center)
        self.move(winrect.topLeft())

        if screen.height() <= 1080 - 50 | screen.width() <= 1920:
            self.showMaximized()

        self.show()
        if sys.platform == "darwin":
            self.raise_()


    #----------------------------------------------------------------------

    def OnCopyWidget(self):
        if self.active_widget:
            self.active_widget.OnCopy()

    #----------------------------------------------------------------------

    def initToolbar(self):

        self.actionOpen = QtWidgets.QAction(self)
        self.actionOpen.setObjectName('actionOpen')
        self.actionOpen.setIcon(QtGui.QIcon(resource_path(os.path.join('images','document-open.png'))))
        self.toolbar = self.addToolBar('actionOpen')
        self.toolbar.addAction(self.actionOpen)
        self.actionOpen.triggered.connect(self.LoadStack)

        self.actionOpenSL = QtWidgets.QAction(self)
        self.actionOpenSL.setObjectName('actionOpenSL')
        self.actionOpenSL.setIcon(QtGui.QIcon(resource_path(os.path.join('images','open-sl.png'))))
        self.toolbar.addAction(self.actionOpenSL)
        self.actionOpenSL.triggered.connect(self.BuildStack)

        self.actionSave = QtWidgets.QAction(self)
        self.actionSave.setObjectName('actionSave')
        self.actionSave.setIcon(QtGui.QIcon(resource_path(os.path.join('images','media-floppy.png'))))
        self.toolbar.addAction(self.actionSave)
        self.actionSave.triggered.connect(self.onSaveResultsToH5)
        self.actionSave.setEnabled(False)

        self.actionInfo = QtWidgets.QAction(self)
        self.actionInfo.setObjectName('actionInfo')
        self.actionInfo.setIcon(QtGui.QIcon(resource_path(os.path.join('images','help-browser.png'))))
        self.toolbar.addAction(self.actionInfo)
        self.actionInfo.triggered.connect(self.onAbout)

#----------------------------------------------------------------------
    def LoadStack(self):
        """
        Browse for a stack file:
        """

        filepath, plugin = File_GUI.SelectFile('read','stack')
        if filepath is not None:
            if plugin is None: # auto-assign appropriate plugin
                plugin = file_plugins.identify(filepath)
            FileStruct = file_plugins.GetFileStructure(filepath, plugin=plugin)
            FileInternalSelection = [(0,0)]
            if FileStruct is not None:
                dlg = File_GUI.DataChoiceDialog(filepath=filepath, filestruct=FileStruct, plugin=plugin)
                if not dlg.exec_():
                    return # do nothing if GUI is cancelled
                FileInternalSelection = dlg.selection
                if dlg.filepath != filepath:
                    filepath = dlg.filepath
                    plugin = file_plugins.identify(dlg.filepath)
            if plugin is None:
                QtWidgets.QMessageBox.warning(self, 'Error!', "Unknown file type")
                return

            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))

            if self.common.stack_loaded == 1:
                self.new_stack_refresh()
                self.stk.new_data()
                self.anlz.delete_data()
            try:    #if checkboxes exist, return checked/unchecked
                JSONconvert = dlg.jsoncheck.isChecked()
                #ringnorm = dlg.ringnormcheck.isChecked()
            except UnboundLocalError: #if checkboxes missing, i.e. for prenormalized data, *.ncb, etc.
                print("DataChoiceDialog skipped")
                JSONconvert = None
                #ringnorm = None
            file_plugins.load(filepath, stack_object=self.stk, plugin=plugin, selection=FileInternalSelection,json=JSONconvert)
            
            # Check for duplicates immediately after load
            self.check_energy_duplicates()
            
            directory = os.path.dirname(str(filepath))


            self.tab_prep.filename = os.path.basename(str(filepath))


            #Update widgets
            x=self.stk.n_cols
            y=self.stk.n_rows
            self.tab_prep.imgrgb = np.zeros(x*y*3,dtype = "uint8")
            self.tab_prep.maxval = np.amax(self.stk.absdata)


            self.ix = int(x/2)
            self.iy = int(y/2)

            self.tab_prep.ix = self.ix
            self.tab_prep.iy = self.iy

            #self.iev = 0
            #self.tab_load.slider_eng.setRange(0,self.stk.n_ev-1)
            #self.tab_load.iev = self.iev
            #self.tab_load.slider_eng.setValue(self.iev)

            #self.tab_prep.slider_eng.setRange(0,self.stk.n_ev-1)
            #self.tab_prep.iev = self.iev
            #self.tab_prep.slider_eng.setValue(self.iev)
            #if showmaptab:
            #    self.tab_map.Clear()
            #    self.tab_map.slider_eng.setRange(0,self.stk.n_ev-1)
            self.stk.setScale()
            self.common.stack_loaded = 1
            self.common.path = directory
            self.common.fntocaption = self.tab_load.FnOverlayCheckBox.isChecked()

            if self.stk.data_struct.spectromicroscopy.normalization.white_spectrum is not None:
                if self.stk.calculate_optical_density() is False:
                    print("Normalization aborted due to duplicates.")
                    self.common.i0_loaded = 0
                else:
                    self.stk.fill_h5_struct_normalization()
                    self.common.i0_loaded = 1

            #self.tab_load.Clear()
            self.ShowInfo(self.tab_prep.filename, directory)
            self.tab_load.absimgfig.loadNewImage()
            #self.tab_prep.ResetDisplaySettings()
            self.tab_prep.absimgfig.loadNewImageWithROI()
            self.tab_prep.button_multicrop.setText('Crop stack 3D...')
            #print(x,y), (self.ix,self.iy), self.stk.absdata.shape
            self.tab_prep.specfig.ClearandReload()
            #self.tab_prep.textctrl.setText(self.tab_prep.filename)

            self.tab_peak.updatewidgets()

            QtWidgets.QApplication.restoreOverrideCursor()

            #if showmaptab:
            #    self.tab_map.Clear()
            #    self.tab_map.loadData()
        self.refresh_widgets()

#-----------------------------------------------------------------------
    def BuildStack(self):
        """
        Browse for .sm files
        """

        #try:
        if True:
            directory = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose a directory", '', QtWidgets.QFileDialog.ShowDirsOnly|QtWidgets.QFileDialog.ReadOnly )



            if directory == '':
                return

            self.common.path = directory

            stackframe = StackListFrame(self, directory, self.common, self.stk, self.data_struct)
            stackframe.show()

#         except:
#             print 'Error could not build stack list.'
#             self.common.stack_loaded = 0
#             self.common.i0_loaded = 0
#             self.new_stack_refresh()
#             self.refresh_widgets()
#
#             QtWidgets.QMessageBox.warning(self,'Error',"Error could not build stack list")
#             import sys; print sys.exc_info()


#-----------------------------------------------------------------------
    def LoadStack4D(self):

        try:
        #if True:

            wildcard = "Supported 4D formats (*.hdf5 *.ncb);;HDF5 files (*.hdf5);;NCB files (*.ncb);;"

            filepath, _filter = QtWidgets.QFileDialog.getOpenFileNames(self, 'Open files', '', wildcard)

            if filepath == '':
                return

            filenames = []
            for name in filepath:
                filenames.append(str(name))


            directory =  os.path.dirname(str(filenames[0]))
            self.tab_prep.filename = os.path.basename(filenames[0])

            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
            basename, extension = os.path.splitext(self.tab_prep.filename)
            extension = extension.strip()

            self.common.path = directory
            self.common.filename = self.tab_prep.filename

            if extension == '.hdf5':
                if self.common.stack_loaded == 1:
                    self.new_stack_refresh()
                    self.stk.new_data()
                    self.anlz.delete_data()
                self.stk.read_h54D(filenames[0])
                
                self.check_energy_duplicates()

                print('Finished reading 4D stack', filenames[0])

            elif extension == '.ncb':
                if self.common.stack_loaded == 1:
                    self.new_stack_refresh()
                    self.stk.new_data()
                    #self.stk.data_struct.delete_data()
                    self.anlz.delete_data()
                self.stk.read_ncb4D(filenames)
                #Get energy list
                wildcard =  "Text file (*.txt);;"
                engfilepath, _filter = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', directory, wildcard)
                self.stk.read_ncb4Denergy(str(engfilepath))


            #Update widgets
            if not hasattr(self.stk, 'theta') or type(self.stk.theta) == int:
                raise TypeError("Not a 4D stack")
            self.common.stack_4d = 1
            self._set_tomo_tab_enabled(True)

            x=self.stk.n_cols
            y=self.stk.n_rows
            self.tab_prep.imgrgb = np.zeros(x*y*3,dtype = "uint8")
            self.tab_prep.maxval = np.amax(self.stk.absdata)


            self.ix = int(x/2)
            self.iy = int(y/2)

            self.tab_prep.ix = self.ix
            self.tab_prep.iy = self.iy

            self.iev = int(self.stk.n_ev/2)
            #self.tab_load.slider_eng.setRange(0,self.stk.n_ev-1)
            self.tab_load.iev = self.iev
            #self.tab_load.slider_eng.setValue(self.iev)

            #self.tab_prep.slider_eng.setRange(0,self.stk.n_ev-1)
            self.tab_prep.iev = self.iev
            #self.tab_prep.slider_eng.setValue(self.iev)

            self.tab_load.slider_theta.setVisible(True)
            #self.tab_load.tc_imagetheta.setVisible(True)
            self.itheta = 0
            self.tab_load.slider_theta.setRange(0,self.stk.n_theta-1)
            self.tab_load.itheta = self.itheta
            self.tab_load.slider_theta.setValue(self.itheta)
            #self.tab_load.tc_imagetheta.setText("4D Data Angle: "+str(self.stk.theta[self.itheta]))


            self.tab_prep.slider_theta.setVisible(True)
            #self.tab_prep.tc_imagetheta.setVisible(True)
            self.tab_prep.slider_theta.setRange(0,self.stk.n_theta-1)
            self.tab_prep.itheta = self.itheta
            self.tab_prep.slider_theta.setValue(self.itheta)
            #self.tab_prep.tc_imagetheta.setText("4D Data Angle: "+str(self.stk.theta[self.itheta]))

            self.tab_clus.slider_theta.setVisible(True)
            self.tab_clus.slider_theta.setEnabled(True)
            self.tab_clus.slider_theta.setRange(0, self.stk.n_theta - 1)
            self.tab_clus.itheta = self.itheta
            self.tab_clus.slider_theta.setValue(self.itheta)

            self.tab_clus.button_calcpca4D.setVisible(True) if hasattr(self.tab_clus, 'button_calcpca4D') else None
            self.tab_spec.button_calc4d.setVisible(True)


            self.common.stack_loaded = 1

            if self.stk.data_struct.spectromicroscopy.normalization.white_spectrum is not None:
                self.common.i0_loaded = 1

            self.tab_load.absimgfig.loadNewImage()
            self.ShowInfo(self.tab_prep.filename, directory)
            #self.tab_prep.ResetDisplaySettings()
            self.tab_prep.absimgfig.loadNewImageWithROI()
            self.tab_prep.button_multicrop.setText('Crop stack 4D...')
            self.tab_prep.specfig.ClearandReload()
            # self.tab_prep.textctrl.setText(self.tab_prep.filename)

            self.tab_peak.updatewidgets()

            QtWidgets.QApplication.restoreOverrideCursor()

            #if showmaptab:
            #    self.tab_map.Clear()
            #    self.tab_map.loadData()
        except:

            self.common.stack_loaded = 0
            self.common.i0_loaded = 0
            self.new_stack_refresh()
            self.tab_prep.button_multicrop.setText('Crop stack 3D/4D...')
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, 'Error', 'Image stack not loaded.')

            import sys
            print(sys.exc_info())


        self.refresh_widgets()

#----------------------------------------------------------------------
    def OnSaveProcessedStack(self, event):
        self.SaveProcessedStack()


#----------------------------------------------------------------------
    def SaveProcessedStack(self):

        """
        Export processed stack to file
        """
        filepath, plugin = File_GUI.SelectFile('write','stack')
        if filepath is not None and plugin is not None:
            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
            try:
                file_plugins.save(filepath, self.stk, 'stack', plugin=plugin)
                QtWidgets.QApplication.restoreOverrideCursor()
            except:
                QtWidgets.QApplication.restoreOverrideCursor()
                QtWidgets.QMessageBox.warning(self, 'Error', 'Could not save processed stack file.')




        return


#----------------------------------------------------------------------
    def onSaveResultsToH5(self, event):
        self.SaveResultsToH5()

#----------------------------------------------------------------------
    def SaveResultsToH5(self):

        """
        Browse for .hdf5 file
        """

        try:

            #print self.data_struct.exchange.data.shape
            #print self.data_struct.exchange.data_axes

            wildcard = "HDF5 files (*.hdf5)"

            filepath, _filter = QtWidgets.QFileDialog.getSaveFileName(self, 'Save as .hdf5', '', wildcard)

            filepath = str(filepath)
            if filepath == '':
                return


            directory =  os.path.dirname(str(filepath))
            self.tab_prep.filename =  os.path.basename(str(filepath))

            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))

            self.common.path = directory
            self.common.filename = self.tab_prep.filename


            file_dataexch_hdf5.write_results_h5(filepath, self.data_struct, self.anlz)

            QtWidgets.QApplication.restoreOverrideCursor()

        except:

            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, 'Error', 'Could not save HDF5 file.')


        self.refresh_widgets()

        return

#----------------------------------------------------------------------
    def onAbout(self):
        dlg = AboutFrame(self)
        dlg.exec_()

# ----------------------------------------------------------------------
    def ShowInfo(self, filename, filepath):
        for i in range(self.tabs.count()):
            page = self.tabs.widget(i)
            tc_file = page.findChild(QtWidgets.QLabel, "tc_file")
            tc_path = page.findChild(QtWidgets.QLabel, "tc_path")
            if tc_file:
                tc_file.setText(filename)
            if tc_path:
                tc_path.setText(filepath)

# ----------------------------------------------------------------------

    def syncFnOverlayCheckBox(self, state):
        self.common.fntocaption = state
        for i in range(self.tabs.count()):
            tab = self.tabs.widget(i)
            if hasattr(tab, 'FnOverlayCheckBox'):
                fig = getattr(tab, 'absimgfig', None) or getattr(tab, 'climgfig', None)
                if fig:
                    pca = getattr(tab, 'calcpca', None)
                    fig.setCaption(setpcalabel=pca)
                spec = getattr(tab, 'specfig', None)
                if spec:
                    spec.setCaption()
                tab.FnOverlayCheckBox.stateChanged.disconnect()
                tab.FnOverlayCheckBox.setChecked(state)
                tab.FnOverlayCheckBox.stateChanged.connect(self.syncFnOverlayCheckBox)


#----------------------------------------------------------------------
    def refresh_widgets(self):

        if self.common.stack_loaded == 0:

            self.tab_prep.FnOverlayCheckBox.setDisabled(True)
            self.tab_load.FnOverlayCheckBox.setDisabled(True)
            self.tab_clus.FnOverlayCheckBox.setDisabled(True)
            self.tab_prep.button_i0ffile.setEnabled(False)
            self.tab_prep.button_i0.setEnabled(False)
            self.tab_prep.button_artefacts.setEnabled(False)
            self.tab_prep.button_prenorm.setEnabled(False)
            self.tab_prep.button_refimgs.setEnabled(False)
            self.tab_prep.button_multicrop.setEnabled(False)
            #self.tab_prep.button_limitev.setEnabled(False)
            #self.tab_prep.button_subregion.setEnabled(False)
            self.tab_prep.button_darksig.setEnabled(False)
            self.tab_prep.button_save.setEnabled(False)
            self.tab_prep.button_savestack.setEnabled(False)
            self.tab_prep.button_align.setEnabled(False)
            self.tab_prep.button_meanflux.setEnabled(False)
            self.tab_prep.button_slideshow.setEnabled(False)
            self.tab_prep.button_spectralROI.setEnabled(False)
            self.tab_prep.button_ROIdosecalc.setEnabled(False)
            self.tab_prep.button_lockspectrum.setEnabled(False)
            self.tab_prep.button_clearlastroi.setEnabled(False)
            self.tab_prep.button_mergeroi.setEnabled(False)
            self.tab_prep.button_subtractroi.setEnabled(False)
            self.tab_prep.button_clearspecfig.setEnabled(False)
            self.tab_prep.ROIShapeBox.setEnabled(False)
            self.tab_prep.ROIvisibleCheckBox.setEnabled(False)
            #self.tab_prep.button_resetdisplay.setEnabled(False)
            #self.tab_prep.button_despike.setEnabled(False)
            #self.tab_prep.button_displaycolor.setEnabled(False)
            self.actionSave.setEnabled(False)

        else:
            self.tab_prep.FnOverlayCheckBox.setDisabled(False)
            self.tab_load.FnOverlayCheckBox.setDisabled(False)
            self.tab_clus.FnOverlayCheckBox.setDisabled(False)
            self.tab_prep.button_i0ffile.setEnabled(True)
            self.tab_prep.button_i0.setEnabled(True)
            self.tab_prep.button_artefacts.setEnabled(True)
            self.tab_prep.button_prenorm.setEnabled(True)
            self.tab_prep.button_multicrop.setEnabled(True)
            #self.tab_prep.button_limitev.setEnabled(True)
            if self.common.stack_4d == 0:
                self.tab_prep.button_refimgs.setEnabled(True)
                #self.tab_prep.button_subregion.setEnabled(True)
                self.tab_prep.button_darksig.setEnabled(True)
            else:
                self.tab_prep.button_refimgs.setEnabled(False)
                #self.tab_prep.button_subregion.setEnabled(False)
                self.tab_prep.button_darksig.setEnabled(False)
            self.tab_load.button_save.setEnabled(True)
            self.tab_prep.button_save.setEnabled(True)
            self.tab_load.pb_copy_img.setEnabled(True)
            self.tab_prep.pb_copy_img.setEnabled(True)
            self.tab_prep.pb_copy_specimg.setEnabled(True)
            self.tab_prep.button_savestack.setEnabled(True)
            self.tab_prep.button_align.setEnabled(True)
            self.tab_prep.button_meanflux.setEnabled(True)
            self.tab_prep.button_slideshow.setEnabled(True)
            self.tab_prep.button_spectralROI.setEnabled(True)
            self.tab_prep.button_ROIdosecalc.setEnabled(True)
            self.tab_prep.button_lockspectrum.setEnabled(True)
            self.tab_prep.button_clearlastroi.setEnabled(True)
            self.tab_prep.button_mergeroi.setEnabled(True)
            self.tab_prep.button_subtractroi.setEnabled(True)
            self.tab_prep.button_clearspecfig.setEnabled(True)
            self.tab_prep.ROIShapeBox.setEnabled(True)
            self.tab_prep.ROIvisibleCheckBox.setEnabled(True)
            #self.tab_prep.button_resetdisplay.setEnabled(True)
            #self.tab_prep.button_despike.setEnabled(True)
            #self.tab_prep.button_displaycolor.setEnabled(True)
            self.actionSave.setEnabled(True)


        if self.tab_tomo is not None:
            self._set_tomo_tab_enabled(bool(self.common.stack_loaded and self.common.stack_4d))

        if self.common.i0_loaded == 0:
            self.tab_prep.specfig.I0Cancel()
            self.tab_prep.button_showi0.setEnabled(False)
            self.tab_prep.button_showi0.setChecked(False)
            self.tab_clus.button_calcpca.setEnabled(False)
            self.tab_spec.button_loadtspec.setEnabled(False)
            self.tab_spec.button_addflat.setEnabled(False)
            self.tab_peak.button_save.setEnabled(False)
            if self.tab_nnma:
                self.tab_nnma.button_calcnnma.setEnabled(False)
                self.tab_nnma.button_muroi.setEnabled(False)
                self.tab_nnma.button_mufile.setEnabled(False)
                self.tab_nnma.button_murand.setEnabled(False)
            if self.tab_tomo is not None:
                self.tab_tomo.button_engdata.setEnabled(False)
        else:
            self.tab_prep.button_showi0.setEnabled(True)
            self.tab_clus.button_calcpca.setEnabled(True)
            self.tab_spec.button_loadtspec.setEnabled(True)
            self.tab_spec.button_addflat.setEnabled(True)
            self.tab_peak.button_save.setEnabled(True)
            if self.tab_nnma:
                self.tab_nnma.button_calcnnma.setEnabled(True)
                self.tab_nnma.button_muroi.setEnabled(True)
                self.tab_nnma.button_mufile.setEnabled(True)
                self.tab_nnma.button_murand.setEnabled(True)
            if self.tab_tomo is not None:
                self.tab_tomo.button_engdata.setEnabled(True)

        if self.common.pca_calculated == 0:
            self.tab_clus.button_savepca.setEnabled(False)
            self.tab_clus.button_toggle_pca_plot.setEnabled(False)
            self.tab_clus.npcaspin.setEnabled(False)
            self.tab_clus.button_movepcup.setEnabled(False)
            self.tab_clus.button_calcca.setEnabled(False)
            self.tab_clus.nclusterspin.setEnabled(False)
            self.tab_clus.remove1stpcacb.setEnabled(False)
            self.tab_clus.cb_splitclusters.setEnabled(False)
            self.tab_clus.ntc_pcscaling_2.setEnabled(False)
            self.tab_spec.rb_fit.setEnabled(False)
        else:
            self.tab_clus.button_calcca.setEnabled(True)
            self.tab_clus.nclusterspin.setEnabled(True)
            self.tab_clus.remove1stpcacb.setEnabled(True)
            self.tab_clus.cb_splitclusters.setEnabled(True)
            self.tab_clus.ntc_pcscaling_2.setEnabled(True)
            self.tab_spec.rb_fit.setEnabled(True)

            # In PCA_new cluster mode, lock the PCA block (except Calculate PCA).
            if self.tab_clus.cluster_display_mode:
                self.tab_clus.button_savepca.setEnabled(False)
                self.tab_clus.button_toggle_pca_plot.setEnabled(False)
                self.tab_clus.button_movepcup.setEnabled(False)
                self.tab_clus.npcaspin.setEnabled(False)
            else:
                self.tab_clus.button_savepca.setEnabled(True)
                self.tab_clus.button_toggle_pca_plot.setEnabled(True)
                self.tab_clus.npcaspin.setEnabled(True)
                # Keep Move PC Up in sync with current plot view state.
                self.tab_clus.UpdatePCAPlotView()


        if self.common.cluster_calculated == 0:
            self.tab_spec.button_addclspec.setEnabled(False)
            self.tab_peak.button_addclspec.setEnabled(False)
            self.tab_clus.button_scatterplots.setEnabled(False)
            self.tab_clus.button_savecluster.setEnabled(False)
            self.tab_clus.showallspectracb.setEnabled(False)
            self.tab_clus.button_errormap.setEnabled(False)
            if self.tab_nnma:
                self.tab_nnma.button_mucluster.setEnabled(False)
        else:
            self.tab_spec.button_addclspec.setEnabled(True)
            self.tab_peak.button_addclspec.setEnabled(True)
            self.tab_clus.button_scatterplots.setEnabled(True)
            self.tab_clus.button_savecluster.setEnabled(True)
            self.tab_clus.showallspectracb.setEnabled(True)
            self.tab_clus.button_errormap.setEnabled(True)
            if self.tab_nnma:
                self.tab_nnma.button_mucluster.setEnabled(True)


        if self.common.spec_anl_calculated == 0:
            self.tab_spec.button_removespec.setEnabled(False)
            self.tab_spec.button_movespdown.setEnabled(False)
            self.tab_spec.button_movespup.setEnabled(False)
            self.tab_spec.button_save.setEnabled(False)
            self.tab_spec.button_showrgb.setEnabled(False)
            #self.tab_spec.button_histogram.setEnabled(False)
            self.tab_spec.button_calc4d.setEnabled(False)

        else:
            self.tab_spec.button_removespec.setEnabled(True)
            self.tab_spec.button_movespdown.setEnabled(True)
            self.tab_spec.button_movespup.setEnabled(True)
            self.tab_spec.button_save.setEnabled(True)
            self.tab_spec.button_showrgb.setEnabled(True)
            #self.tab_spec.button_histogram.setEnabled(True)
            self.tab_spec.button_calc4d.setEnabled(True)


        if self.tab_tomo is not None:
            if self.common.spec_anl4D_calculated == 0:
                self.tab_tomo.button_spcomp.setEnabled(False)
            else:
                self.tab_tomo.button_spcomp.setEnabled(True)


        if self.tab_pfit != None:
            if self.common.cluster_calculated == 0:
                self.tab_pfit.button_addclspec.setEnabled(False)
            else:
                self.tab_pfit.button_addclspec.setEnabled(True)

            if self.common.xpf_loaded == 1:
                self.tab_pfit.slider_spec.setEnabled(True)
            else:
                self.tab_pfit.slider_spec.setEnabled(False)


        # ToDo: RestoreResetDisplaysetting functionality # Really needed?
        #self.tab_prep.ResetDisplaySettings()



#-----------------------------------------------------------------------
    def new_stack_refresh(self):


        self.common.i0_loaded = 0
        self.common.pca_calculated = 0
        self.common.cluster_calculated = 0
        self.common.spec_anl_calculated = 0
        self.common.xpf_loaded = 0
        self.common.stack_4d = 0
        self.common.pca4D_calculated = 0
        self.common.spec_anl4D_calculated = 0

        self.refresh_widgets()

        #page 0
        self.tab_load.slider_theta.setVisible(False)
        #self.tab_load.tc_imagetheta.setVisible(False)

        #page 1
        self.tab_prep.button_i0.disconnect()
        self.tab_prep.button_i0.setText("Select I0")
        self.tab_prep.button_i0.clicked.connect(self.tab_prep.specfig.OnI0Histogram)
        #self.tab_prep.rb_flux.setChecked(True)
        #self.tab_prep.rb_od.setChecked(False)
        #self.tab_prep.showflux = True

        #fig = self.tab_prep.specfig
        #fig.clf()
        #self.tab_prep.SpectrumPanel.draw()
        #self.tab_prep.tc_spec.setText("Spectrum at point: ")

        #fig = self.tab_prep.absimgfig
        #fig.clf()
        #self.tab_prep.AbsImagePanel.draw()
        #self.tab_prep.tc_imageeng.setText("Image at energy: ")

        # self.tab_prep.textctrl.setText(' ')

        self.tab_prep.ResetDisplaySettings()
        #page 1
        self.tab_prep.slider_theta.setVisible(False)

        # page 2_1 (PCA & Cluster Analysis)
        self.tab_clus.reset_after_new_stack()
        self.tab_clus.slider_theta.setVisible(False)
        self.tab_clus.slider_theta.setEnabled(False)

        #page 4
        self.tab_spec.ClearWidgets()

        #page 5
        fig = self.tab_peak.kespecfig
        fig.clf()
        self.tab_peak.KESpecPan.draw()
        fig = self.tab_peak.absimgfig
        fig.clf()
        self.tab_peak.AbsImagePanel.draw()

        #page 6
        if self.tab_pfit != None:
            fig = self.tab_pfit.Specfig
            fig.clf()
            self.tab_pfit.SpectrumPanel.draw()
            self.tab_pfit.slider_spec.setEnabled(False)

        # tab_tomo
        if self.tab_tomo is not None:
            self._set_tomo_tab_enabled(False)

    def _set_tomo_tab_enabled(self, enabled):
        if self.tab_tomo is None:
            return
        tomo_index = self.tabs.indexOf(self.tab_tomo)
        if tomo_index == -1:
            return
        self.tabs.setTabEnabled(tomo_index, enabled)
        self._set_tomo_tab_style(enabled)
        if not enabled:
            if self.tabs.currentWidget() == self.tab_tomo:
                self.tabs.setCurrentWidget(self.tab_prep)
            self.tab_tomo.NewStackClear()

    def _set_tomo_tab_style(self, enabled):
        if self.tab_tomo is None:
            return
        tomo_index = self.tabs.indexOf(self.tab_tomo)
        if tomo_index == -1:
            return
        if enabled:
            if self.palette().color(self.palette().Background).red() + self.palette().color(self.palette().Background).green() + self.palette().color(self.palette().Background).blue() < 375:
                color = QtGui.QColor('dodgerblue')
            else:
                color = QtGui.QColor('darkblue')
        else:
            color = QtGui.QColor('gray')
        self.tabs.tabBar().setTabTextColor(tomo_index, color)

    def check_energy_duplicates(self):
        if hasattr(self.stk, 'ev'):
            ev = self.stk.ev
            # Basic check if it's an array/list and not empty scalar
            if hasattr(ev, '__len__') and len(ev) > 1:
                if len(ev) != len(np.unique(ev)):
                    # Temporarily restore cursor for interaction
                    QtWidgets.QApplication.restoreOverrideCursor()
                    
                    msg = QtWidgets.QMessageBox()
                    msg.setIcon(QtWidgets.QMessageBox.Question)
                    msg.setWindowTitle("Duplicate Energies Detected")
                    msg.setText("The loaded stack contains duplicate energy values.")
                    msg.setInformativeText("This causes errors during analysis.\n"
                                           "Do you want to fix this by shifting duplicates by 1E-5 eV?")
                    
                    msg.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
                    msg.setDefaultButton(QtWidgets.QMessageBox.Yes)
                    
                    ret = msg.exec_()
                    
                    if ret == QtWidgets.QMessageBox.Yes:
                        self.stk.ev = self.stk.handle_duplicates(self.stk.ev)
            
            # Check for unsorted energies
            if hasattr(ev, '__len__') and len(ev) > 1:
                # Re-fetch ev in case it was modified
                ev = self.stk.ev
                if not np.all(np.diff(ev) >= 0):
                    # Temporarily restore cursor
                    QtWidgets.QApplication.restoreOverrideCursor()
                    
                    msg = QtWidgets.QMessageBox()
                    msg.setIcon(QtWidgets.QMessageBox.Question)
                    msg.setWindowTitle("Unsorted Energies Detected")
                    msg.setText("The loaded stack has unsorted energy values.")
                    msg.setInformativeText("Energies must be in ascending order for correct analysis.\n"
                                           "Sort automatically?")
                    
                    msg.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
                    msg.setDefaultButton(QtWidgets.QMessageBox.Yes)
                    
                    ret = msg.exec_()
                    
                    if ret == QtWidgets.QMessageBox.Yes:
                        sort_idx = np.argsort(ev)
                        self.stk.sort_by_indices(sort_idx)
                    
                    # Set cursor back to wait
                    QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))

    def ask_fix_duplicates_callback(self, message):
        """
        Callback used by data_stack to ask user for confirmation via GUI.
        """        
        cursor = QtWidgets.QApplication.overrideCursor()
        if cursor:
            QtWidgets.QApplication.restoreOverrideCursor()
            
        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Question)
        msg.setWindowTitle("Warning!")
        msg.setText(message)
        msg.setInformativeText("Fix automatically to proceed?")
        
        msg.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        msg.setDefaultButton(QtWidgets.QMessageBox.Yes)
        
        ret = msg.exec_()
        
        result = False
        if ret == QtWidgets.QMessageBox.Yes:
            result = True
            
        if cursor:
             QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
             
        return result
