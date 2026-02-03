
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
from ..core.constants import Winsizex, Winsizey, showtomotab
from .. import __version__ as version
from ..core.state import common
from ..helpers import resource_path
from ..file_plugins import file_dataexch_hdf5
from .. import file_plugins
from .pages import (
    PageLoadData, PageStack, PagePCA, PagePCACluster,
    PageCluster, PageSpectral, PagePeakID, PageXrayPeakFitting,
    PageNNMA, PageTomo
)
from .utils.file_gui import File_GUI
from .dialogs import AboutFrame, StackListFrame

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

        # create the page windows as tabs
        self.page0 = PageLoadData(self.common, self.data_struct, self.stk)
        self.page1 = PageStack(self.common, self.data_struct, self.stk)
        self.page2_1 = PagePCACluster(self.common, self.data_struct, self.stk, self.anlz)
        self.page2 = PagePCA(self.common, self.data_struct, self.stk, self.anlz)
        self.page3 = PageCluster(self.common, self.data_struct, self.stk, self.anlz)
        self.page4 = PageSpectral(self.common, self.data_struct, self.stk, self.anlz)
        self.page5 = PagePeakID(self.common, self.data_struct, self.stk, self.anlz)
        self.page6 = PageXrayPeakFitting(self.common, self.data_struct, self.stk, self.anlz)
        self.page7 = PageNNMA(self.common, self.data_struct, self.stk, self.anlz, self.nnma)


        self.tabs.addTab(self.page0,"Load Data")
        self.tabs.addTab(self.page1, "Preprocess Data")
        self.tabs.addTab(self.page2, "PCA")
        self.tabs.addTab(self.page2_1,"PCA_new")
        self.tabs.addTab(self.page3,"Cluster Analysis")
        self.tabs.addTab(self.page4,"Spectral Maps")
        self.tabs.addTab(self.page7, "NNMA Analysis")
        self.tabs.addTab(self.page5,"Peak ID")
        self.tabs.addTab(self.page6, "XrayPeakFitting")

        if showtomotab:
            self.page8 = PageTomo(self.common, self.data_struct, self.stk, self.anlz)
            self.tabs.addTab(self.page8, "Tomography")

#        if showmaptab:
#            self.page9 = PageMap(self.common, self.data_struct, self.stk)
#            self.tabs.addTab(self.page9, "Image Maps")

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
            self.tabs.tabBar().setTabTextColor(4, QtGui.QColor('tomato'))
            self.tabs.tabBar().setTabTextColor(5, QtGui.QColor('tomato'))
            self.tabs.tabBar().setTabTextColor(6, QtGui.QColor('orchid'))
            self.tabs.tabBar().setTabTextColor(7, QtGui.QColor('orchid'))
            if showtomotab:
                self.tabs.tabBar().setTabTextColor(8, QtGui.QColor('dodgerblue'))
        else:
            self.tabs.tabBar().setTabTextColor(0, QtGui.QColor('green'))
            self.tabs.tabBar().setTabTextColor(1, QtGui.QColor('green'))
            self.tabs.tabBar().setTabTextColor(2, QtGui.QColor('darkRed'))
            self.tabs.tabBar().setTabTextColor(3, QtGui.QColor('darkRed'))
            self.tabs.tabBar().setTabTextColor(4, QtGui.QColor('darkRed'))
            self.tabs.tabBar().setTabTextColor(5, QtGui.QColor('darkRed'))
            self.tabs.tabBar().setTabTextColor(6, QtGui.QColor('purple'))
            self.tabs.tabBar().setTabTextColor(7, QtGui.QColor('purple'))
            if showtomotab:
                self.tabs.tabBar().setTabTextColor(8, QtGui.QColor('darkblue'))
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
#                 self.page7 = PageNNMA(self.common, self.data_struct, self.stk, self.anlz, self.nnma)
#                 self.tabs.addTab(self.page7, "NNMA Analysis")

        layout = QtWidgets.QVBoxLayout()

        layout.addWidget(self.tabs)
        #self.setCentralWidget(self.tabs)
        self.ShowInfo('File name', os.getcwd())
        self.page0.FnOverlayCheckBox.stateChanged.connect(self.syncFnOverlayCheckBox)
        self.page1.FnOverlayCheckBox.stateChanged.connect(self.syncFnOverlayCheckBox)
        self.page2_1.FnOverlayCheckBox.stateChanged.connect(self.syncFnOverlayCheckBox)

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


            self.page1.filename = os.path.basename(str(filepath))


            #Update widgets
            x=self.stk.n_cols
            y=self.stk.n_rows
            self.page1.imgrgb = np.zeros(x*y*3,dtype = "uint8")
            self.page1.maxval = np.amax(self.stk.absdata)


            self.ix = int(x/2)
            self.iy = int(y/2)

            self.page1.ix = self.ix
            self.page1.iy = self.iy

            #self.iev = 0
            #self.page0.slider_eng.setRange(0,self.stk.n_ev-1)
            #self.page0.iev = self.iev
            #self.page0.slider_eng.setValue(self.iev)

            #self.page1.slider_eng.setRange(0,self.stk.n_ev-1)
            #self.page1.iev = self.iev
            #self.page1.slider_eng.setValue(self.iev)
            #if showmaptab:
            #    self.page9.Clear()
            #    self.page9.slider_eng.setRange(0,self.stk.n_ev-1)
            self.stk.setScale()
            self.common.stack_loaded = 1
            self.common.path = directory
            self.common.fntocaption = self.page0.FnOverlayCheckBox.isChecked()

            if self.stk.data_struct.spectromicroscopy.normalization.white_spectrum is not None:
                if self.stk.calculate_optical_density() is False:
                    print("Normalization aborted due to duplicates.")
                    self.common.i0_loaded = 0
                else:
                    self.stk.fill_h5_struct_normalization()
                    self.common.i0_loaded = 1

            #self.page0.Clear()
            self.ShowInfo(self.page1.filename, directory)
            self.page0.absimgfig.loadNewImage()
            #self.page1.ResetDisplaySettings()
            self.page1.absimgfig.loadNewImageWithROI()
            self.page1.button_multicrop.setText('Crop stack 3D...')
            #print(x,y), (self.ix,self.iy), self.stk.absdata.shape
            self.page1.specfig.ClearandReload()
            #self.page1.textctrl.setText(self.page1.filename)

            self.page5.updatewidgets()

            QtWidgets.QApplication.restoreOverrideCursor()

            #if showmaptab:
            #    self.page9.Clear()
            #    self.page9.loadData()
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
            self.page1.filename = os.path.basename(filenames[0])

            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
            basename, extension = os.path.splitext(self.page1.filename)
            extension = extension.strip()

            self.common.path = directory
            self.common.filename = self.page1.filename

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

            x=self.stk.n_cols
            y=self.stk.n_rows
            self.page1.imgrgb = np.zeros(x*y*3,dtype = "uint8")
            self.page1.maxval = np.amax(self.stk.absdata)


            self.ix = int(x/2)
            self.iy = int(y/2)

            self.page1.ix = self.ix
            self.page1.iy = self.iy

            self.iev = int(self.stk.n_ev/2)
            #self.page0.slider_eng.setRange(0,self.stk.n_ev-1)
            self.page0.iev = self.iev
            #self.page0.slider_eng.setValue(self.iev)

            #self.page1.slider_eng.setRange(0,self.stk.n_ev-1)
            self.page1.iev = self.iev
            #self.page1.slider_eng.setValue(self.iev)

            self.page0.slider_theta.setVisible(True)
            #self.page0.tc_imagetheta.setVisible(True)
            self.itheta = 0
            self.page0.slider_theta.setRange(0,self.stk.n_theta-1)
            self.page0.itheta = self.itheta
            self.page0.slider_theta.setValue(self.itheta)
            #self.page0.tc_imagetheta.setText("4D Data Angle: "+str(self.stk.theta[self.itheta]))


            self.page1.slider_theta.setVisible(True)
            #self.page1.tc_imagetheta.setVisible(True)
            self.page1.slider_theta.setRange(0,self.stk.n_theta-1)
            self.page1.itheta = self.itheta
            self.page1.slider_theta.setValue(self.itheta)
            #self.page1.tc_imagetheta.setText("4D Data Angle: "+str(self.stk.theta[self.itheta]))

            self.page2.button_calcpca4D.setVisible(True)
            self.page4.button_calc4d.setVisible(True)


            self.common.stack_loaded = 1

            if self.stk.data_struct.spectromicroscopy.normalization.white_spectrum is not None:
                self.common.i0_loaded = 1

            self.page0.absimgfig.loadNewImage()
            self.ShowInfo(self.page1.filename, directory)
            #self.page1.ResetDisplaySettings()
            self.page1.absimgfig.loadNewImageWithROI()
            self.page1.button_multicrop.setText('Crop stack 4D...')
            self.page1.specfig.ClearandReload()
            # self.page1.textctrl.setText(self.page1.filename)

            self.page5.updatewidgets()

            QtWidgets.QApplication.restoreOverrideCursor()

            #if showmaptab:
            #    self.page9.Clear()
            #    self.page9.loadData()
        except:

            self.common.stack_loaded = 0
            self.common.i0_loaded = 0
            self.new_stack_refresh()
            self.page1.button_multicrop.setText('Crop stack 3D/4D...')
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
            self.page1.filename =  os.path.basename(str(filepath))

            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))

            self.common.path = directory
            self.common.filename = self.page1.filename


            file_dataexch_hdf5.write_results_h5(filepath, self.data_struct, self.anlz)

            QtWidgets.QApplication.restoreOverrideCursor()

        except:

            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, 'Error', 'Could not save HDF5 file.')


        self.refresh_widgets()

        return

#----------------------------------------------------------------------
    def onAbout(self):

        #self.popup = AboutFrame(self)
        #self.popup.show()
        pass

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

            self.page1.FnOverlayCheckBox.setDisabled(True)
            self.page0.FnOverlayCheckBox.setDisabled(True)
            self.page2_1.FnOverlayCheckBox.setDisabled(True)
            self.page1.button_i0ffile.setEnabled(False)
            self.page1.button_i0.setEnabled(False)
            self.page1.button_artefacts.setEnabled(False)
            self.page1.button_prenorm.setEnabled(False)
            self.page1.button_refimgs.setEnabled(False)
            self.page1.button_multicrop.setEnabled(False)
            #self.page1.button_limitev.setEnabled(False)
            #self.page1.button_subregion.setEnabled(False)
            self.page1.button_darksig.setEnabled(False)
            self.page1.button_save.setEnabled(False)
            self.page1.button_savestack.setEnabled(False)
            self.page1.button_align.setEnabled(False)
            self.page1.button_meanflux.setEnabled(False)
            self.page1.button_slideshow.setEnabled(False)
            self.page1.button_spectralROI.setEnabled(False)
            self.page1.button_lockspectrum.setEnabled(False)
            self.page1.button_clearlastroi.setEnabled(False)
            self.page1.button_mergeroi.setEnabled(False)
            self.page1.button_subtractroi.setEnabled(False)
            self.page1.button_clearspecfig.setEnabled(False)
            self.page1.ROIShapeBox.setEnabled(False)
            self.page1.ROIvisibleCheckBox.setEnabled(False)
            #self.page1.button_resetdisplay.setEnabled(False)
            #self.page1.button_despike.setEnabled(False)
            #self.page1.button_displaycolor.setEnabled(False)
            self.actionSave.setEnabled(False)

        else:
            self.page1.FnOverlayCheckBox.setDisabled(False)
            self.page0.FnOverlayCheckBox.setDisabled(False)
            self.page2_1.FnOverlayCheckBox.setDisabled(False)
            self.page1.button_i0ffile.setEnabled(True)
            self.page1.button_i0.setEnabled(True)
            self.page1.button_artefacts.setEnabled(True)
            self.page1.button_prenorm.setEnabled(True)
            self.page1.button_multicrop.setEnabled(True)
            #self.page1.button_limitev.setEnabled(True)
            if self.common.stack_4d == 0:
                self.page1.button_refimgs.setEnabled(True)
                #self.page1.button_subregion.setEnabled(True)
                self.page1.button_darksig.setEnabled(True)
            else:
                self.page1.button_refimgs.setEnabled(False)
                #self.page1.button_subregion.setEnabled(False)
                self.page1.button_darksig.setEnabled(False)
            self.page0.button_save.setEnabled(True)
            self.page1.button_save.setEnabled(True)
            self.page0.pb_copy_img.setEnabled(True)
            self.page1.pb_copy_img.setEnabled(True)
            self.page1.pb_copy_specimg.setEnabled(True)
            self.page1.button_savestack.setEnabled(True)
            self.page1.button_align.setEnabled(True)
            self.page1.button_meanflux.setEnabled(True)
            self.page1.button_slideshow.setEnabled(True)
            self.page1.button_spectralROI.setEnabled(True)
            self.page1.button_lockspectrum.setEnabled(True)
            self.page1.button_clearlastroi.setEnabled(True)
            self.page1.button_mergeroi.setEnabled(True)
            self.page1.button_subtractroi.setEnabled(True)
            self.page1.button_clearspecfig.setEnabled(True)
            self.page1.ROIShapeBox.setEnabled(True)
            self.page1.ROIvisibleCheckBox.setEnabled(True)
            #self.page1.button_resetdisplay.setEnabled(True)
            #self.page1.button_despike.setEnabled(True)
            #self.page1.button_displaycolor.setEnabled(True)
            self.actionSave.setEnabled(True)


        if self.common.i0_loaded == 0:
            self.page1.specfig.I0Cancel()
            self.page1.button_showi0.setEnabled(False)
            self.page1.button_showi0.setChecked(False)
            #self.page1.rb_flux.setEnabled(False)
            #self.page1.rb_od.setEnabled(False)
            #self.page1.button_reseti0.setEnabled(False)
            #self.page1.button_saveod.setEnabled(False)
            self.page2.button_calcpca.setEnabled(False)
            self.page2_1.button_calcpca.setEnabled(False)
            self.page2.button_calcpca4D.setEnabled(False)
            self.page4.button_loadtspec.setEnabled(False)
            self.page4.button_addflat.setEnabled(False)
            self.page5.button_save.setEnabled(False)
            if self.page7:
                self.page7.button_calcnnma.setEnabled(False)
                self.page7.button_muroi.setEnabled(False)
                self.page7.button_mufile.setEnabled(False)
                self.page7.button_murand.setEnabled(False)
            if showtomotab:
                self.page8.button_engdata.setEnabled(False)
        else:
            self.page1.button_showi0.setEnabled(True)
            #self.page1.rb_flux.setEnabled(True)
            #self.page1.rb_od.setEnabled(True)
            #self.page1.button_reseti0.setEnabled(True)
            #self.page1.button_saveod.setEnabled(True)
            self.page2.button_calcpca.setEnabled(True)
            self.page2_1.button_calcpca.setEnabled(True)
            self.page2.button_calcpca4D.setEnabled(True)
            self.page4.button_loadtspec.setEnabled(True)
            self.page4.button_addflat.setEnabled(True)
            self.page5.button_save.setEnabled(True)
            if self.page7:
                self.page7.button_calcnnma.setEnabled(True)
                self.page7.button_muroi.setEnabled(True)
                self.page7.button_mufile.setEnabled(True)
                self.page7.button_murand.setEnabled(True)
            if showtomotab:
                self.page8.button_engdata.setEnabled(True)

        if self.common.pca_calculated == 0:
            self.page2.button_savepca.setEnabled(False)
            self.page2_1.button_savepca.setEnabled(False)
            self.page2.slidershow.setEnabled(False)
            self.page2.button_movepcup.setEnabled(False)
            self.page3.button_calcca.setEnabled(False)
            self.page4.rb_fit.setEnabled(False)
        else:
            self.page2.button_savepca.setEnabled(True)
            self.page2_1.button_savepca.setEnabled(True)
            self.page2.slidershow.setEnabled(True)
            self.page2.button_movepcup.setEnabled(True)
            self.page3.button_calcca.setEnabled(True)
            self.page4.rb_fit.setEnabled(True)

        if self.common.cluster_calculated == 0:
            self.page3.button_scatterplots.setEnabled(False)
            self.page3.button_savecluster.setEnabled(False)
            self.page3.slidershow.setEnabled(False)
            self.page4.button_addclspec.setEnabled(False)
            self.page5.button_addclspec.setEnabled(False)
            if self.page7:
                self.page7.button_mucluster.setEnabled(False)
        else:
            self.page3.button_scatterplots.setEnabled(True)
            self.page3.button_savecluster.setEnabled(True)
            self.page3.slidershow.setEnabled(True)
            self.page4.button_addclspec.setEnabled(True)
            self.page5.button_addclspec.setEnabled(True)
            if self.page7:
                self.page7.button_mucluster.setEnabled(True)


        if self.common.spec_anl_calculated == 0:
            self.page4.button_removespec.setEnabled(False)
            self.page4.button_movespdown.setEnabled(False)
            self.page4.button_movespup.setEnabled(False)
            self.page4.button_save.setEnabled(False)
            self.page4.button_showrgb.setEnabled(False)
            #self.page4.button_histogram.setEnabled(False)
            self.page4.button_calc4d.setEnabled(False)

        else:
            self.page4.button_removespec.setEnabled(True)
            self.page4.button_movespdown.setEnabled(True)
            self.page4.button_movespup.setEnabled(True)
            self.page4.button_save.setEnabled(True)
            self.page4.button_showrgb.setEnabled(True)
            #self.page4.button_histogram.setEnabled(True)
            self.page4.button_calc4d.setEnabled(True)


        if showtomotab:
            if self.common.spec_anl4D_calculated == 0:
                self.page8.button_spcomp.setEnabled(False)
            else:
                self.page8.button_spcomp.setEnabled(True)


        if self.page6 != None:
            if self.common.cluster_calculated == 0:
                self.page6.button_addclspec.setEnabled(False)
            else:
                self.page6.button_addclspec.setEnabled(True)

            if self.common.xpf_loaded == 1:
                self.page6.slider_spec.setEnabled(True)
            else:
                self.page6.slider_spec.setEnabled(False)


        # ToDo: RestoreResetDisplaysetting functionality # Really needed?
        #self.page1.ResetDisplaySettings()



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
        self.page0.slider_theta.setVisible(False)
        #self.page0.tc_imagetheta.setVisible(False)

        #page 1
        self.page1.button_i0.disconnect()
        self.page1.button_i0.setText("Select I0")
        self.page1.button_i0.clicked.connect(self.page1.specfig.OnI0Histogram)
        #self.page1.rb_flux.setChecked(True)
        #self.page1.rb_od.setChecked(False)
        #self.page1.showflux = True

        #fig = self.page1.specfig
        #fig.clf()
        #self.page1.SpectrumPanel.draw()
        #self.page1.tc_spec.setText("Spectrum at point: ")

        #fig = self.page1.absimgfig
        #fig.clf()
        #self.page1.AbsImagePanel.draw()
        #self.page1.tc_imageeng.setText("Image at energy: ")

        # self.page1.textctrl.setText(' ')

        self.page1.ResetDisplaySettings()
        #page 0
        self.page1.slider_theta.setVisible(False)
        #self.page1.tc_imagetheta.setVisible(False)

        #page 2
        fig = self.page2.pcaevalsfig
        fig.clf()
        self.page2.PCAEvalsPan.draw()

        fig = self.page2.pcaimgfig
        fig.clf()
        self.page2.PCAImagePan.draw()

        fig = self.page2.pcaspecfig
        fig.clf()
        self.page2.PCASpecPan.draw()

        self.page2.vartc.setText('0%')
        self.page2.npcaspin.setValue(1)
        self.page2.tc_PCAcomp.setText("PCA component ")
        self.page2.text_pcaspec.setText("PCA spectrum ")

        self.page2.selpca = 1
        self.page2.numsigpca = 2
        self.page2.slidershow.setValue(self.page2.selpca)

        self.page2.slider_theta.setVisible(False)
        self.page2.tc_imagetheta.setVisible(False)
        self.page2.button_calcpca4D.setVisible(False)

        #page 3
        fig = self.page3.clusterimgfig
        fig.clf()
        self.page3.ClusterImagePan.draw()

        fig = self.page3.clusterindvimgfig
        fig.clf()
        self.page3.ClusterIndvImagePan.draw()

        fig = self.page3.clusterspecfig

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
        fig.clf()
        self.page3.ClusterSpecPan.draw()

        fig = self.page3.clusterdistmapfig
        fig.clf()
        self.page3.ClusterDistMapPan.draw()

        self.page3.selcluster = 1
        self.page3.slidershow.setValue(self.page3.selcluster)
        self.page3.numclusters = 5
        self.page3.nclusterspin.setValue(self.page3.numclusters)
        self.page3.tc_cluster.setText("Cluster ")
        self.page3.tc_clustersp.setText("Cluster spectrum")
        self.page3.wo_1st_pca = 0
        self.page3.remove1stpcacb.setChecked(False)


        #page 4
        self.page4.ClearWidgets()

        #page 5
        fig = self.page5.kespecfig
        fig.clf()
        self.page5.KESpecPan.draw()
        fig = self.page5.absimgfig
        fig.clf()
        self.page5.AbsImagePanel.draw()


        #page 6
        if self.page6 != None:
            fig = self.page6.Specfig
            fig.clf()
            self.page6.SpectrumPanel.draw()
            self.page6.slider_spec.setEnabled(False)

        #page8
        if showtomotab:
            self.page8.NewStackClear()
