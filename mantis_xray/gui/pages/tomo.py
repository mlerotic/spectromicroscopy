
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
import os
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.widgets import LassoSelector
import matplotlib
from ...analysis import tomo_reconstruction
from ...file_plugins import file_ncb
from ...core.constants import PlotH
from ..dialogs.roi_histogram import ROIHistogram
from ..dialogs.plot_frame import PlotFrame

# Helper rebin function (originally in mantis_qt.py)
def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

""" ------------------------------------------------------------------------------------------------"""
class PageTomo(QtWidgets.QWidget):
    def __init__(self, common, data_struct, stack, anlz):
        super(PageTomo, self).__init__()

        self.initUI(common, data_struct, stack, anlz)

#----------------------------------------------------------------------
    def initUI(self, common, data_struct, stack, anlz):

        self.data_struct = data_struct
        self.stack = stack
        self.com = common
        self.anlz = anlz

        self.tr = tomo_reconstruction.Ctomo(self.stack)
        self.theta = []
        self.tomodata = None

        self.algonames = ['Compressed Sensing', 'SIRT']
        self.algo = 0

        self.fulltomorecdata = []
        self.ncomponents = 0

        self.tomo_calculated = 0
        self.full_tomo_calculated = 0
        self.energiesloaded = 0
        self.preview_loaded = False
        self.single_data_is_volume = False
        self.single_volume_data = None

        self.datanames = []

        self.haveROI = 0
        self.ROIarray = []
        self.ROIvol = []

        self.select1 = 0

        self.icomp = 0
        self.islice = 0

        self.maxIters = 10
        self.beta = 0.5
        self.engpar = 0.001
        self.useengreg = 0
        self.samplethick = 0
        self.nonnegconst = 1

        self.nprocessors = 6


        #panel 1
        sizer1 = QtWidgets.QGroupBox('Tomo Data')
        vbox1 = QtWidgets.QVBoxLayout()


        self.button_spcomp = QtWidgets.QPushButton('Load Tomo Data for Spectral Components')
        self.button_spcomp.clicked.connect( self.OnLoadTomoComponents)
        self.button_spcomp.setEnabled(False)
        vbox1.addWidget(self.button_spcomp)

        self.button_engdata = QtWidgets.QPushButton('Load Tomo Data for each Energy')
        self.button_engdata.clicked.connect( self.OnLoadTomoEng)
        self.button_engdata.setEnabled(False)
        vbox1.addWidget(self.button_engdata)

        self.button_expdata = QtWidgets.QPushButton('Export Tomo Data as .mrc')
        self.button_expdata.clicked.connect( self.OnExportData)
        self.button_expdata.setEnabled(False)
        vbox1.addWidget(self.button_expdata)

        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.HLine)
        line.setFrameShadow(QtWidgets.QFrame.Sunken)

        vbox1.addStretch(1)
        vbox1.addWidget(line)
        vbox1.addStretch(1)

        self.button_loadmrc = QtWidgets.QPushButton('Load Single Tomo Dataset')
        self.button_loadmrc.clicked.connect( self.OnLoadSingleMrc)
        vbox1.addWidget(self.button_loadmrc)

        sizer1.setLayout(vbox1)

        #panel 2
        sizer2 = QtWidgets.QGroupBox('Tomo Reconstruction')
        vbox2 = QtWidgets.QVBoxLayout()


        self.button_calc1 = QtWidgets.QPushButton( 'Calculate One Dataset')
        self.button_calc1.clicked.connect( self.OnCalcTomo1)
        self.button_calc1.setEnabled(False)
        vbox2.addWidget(self.button_calc1)

        hbox21 = QtWidgets.QHBoxLayout()
        tc1 = QtWidgets.QLabel(self)
        hbox21.addWidget(tc1)
        tc1.setText('Choose Tomo Dataset:')
        self.combonames = QtWidgets.QComboBox(self)
        self.combonames.activated[int].connect(self.OnSelect1Comp)
        hbox21.addWidget(self.combonames)
        vbox2.addLayout(hbox21)

        hbox22 = QtWidgets.QHBoxLayout()
        tc2 = QtWidgets.QLabel(self)
        hbox22.addWidget(tc2)
        tc2.setText('Binning:  ')
        self.combobin = QtWidgets.QComboBox(self)
        self.combobin.addItems(['1','2','4','8'])
        hbox22.addWidget(self.combobin)
        vbox2.addLayout(hbox22)

        self.button_calcall = QtWidgets.QPushButton( 'Calculate All Datasets')
        self.button_calcall.clicked.connect( self.OnCalcTomoFull)
        self.button_calcall.setEnabled(False)
        vbox2.addWidget(self.button_calcall)


        self.button_save = QtWidgets.QPushButton( 'Save as .mrc')
        self.button_save.clicked.connect( self.OnSave)
        self.button_save.setEnabled(False)
        vbox2.addWidget(self.button_save)

        self.button_saveall = QtWidgets.QPushButton( 'Save All as .mrc')
        self.button_saveall.clicked.connect( self.OnSaveAll)
        self.button_saveall.setEnabled(False)
        vbox2.addWidget(self.button_saveall)

        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.HLine)
        line.setFrameShadow(QtWidgets.QFrame.Sunken)


        vbox2.addStretch(1)
        vbox2.addWidget(line)
        vbox2.addStretch(1)

#         hbox21 = QtWidgets.QHBoxLayout()
#         tc1 = QtWidgets.QLabel(self)
#         hbox21.addWidget(tc1)
#         tc1.setText('Choose Tomo Dataset: ')
        self.comboalgos = QtWidgets.QComboBox(self)
        self.comboalgos.activated[int].connect(self.OnSelectAlgo)
        self.comboalgos.addItems(self.algonames)
        vbox2.addWidget(self.comboalgos)

        hbox22 = QtWidgets.QHBoxLayout()
        text1 = QtWidgets.QLabel(self)
        text1.setText('Number of iterations')
        hbox22.addWidget(text1)
        hbox22.addStretch(1)
        self.ntc_niterations = QtWidgets.QLineEdit(self)
        self.ntc_niterations.setFixedWidth(65)
        self.ntc_niterations.setValidator(QtGui.QIntValidator(1, 99999, self))
        self.ntc_niterations.setAlignment(QtCore.Qt.AlignRight)
        self.ntc_niterations.setText(str(self.maxIters))
        hbox22.addWidget(self.ntc_niterations)
        vbox2.addLayout(hbox22)

        hbox23 = QtWidgets.QHBoxLayout()
        self.tc_par = QtWidgets.QLabel(self)
        self.tc_par.setText("CS Parameter Beta")
        hbox23.addWidget(self.tc_par)
        self.ntc_beta = QtWidgets.QLineEdit(self)
        self.ntc_beta.setFixedWidth(65)
        self.ntc_beta.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.ntc_beta.setAlignment(QtCore.Qt.AlignRight)
        hbox23.addStretch(1)
        self.ntc_beta.setText(str(self.beta))
        hbox23.addWidget(self.ntc_beta)

        vbox2.addLayout(hbox23)


        hbox25 = QtWidgets.QHBoxLayout()

        self.cb_ereg = QtWidgets.QCheckBox('CS Energy Reg Parameter', self)
        self.cb_ereg.setChecked(False)
        self.cb_ereg.stateChanged.connect(self.OnCBEngReg)
        hbox25.addWidget(self.cb_ereg)

        self.ntc_ereg = QtWidgets.QLineEdit(self)
        self.ntc_ereg.setFixedWidth(65)
        self.ntc_ereg.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.ntc_ereg.setAlignment(QtCore.Qt.AlignRight)
        hbox25.addStretch(1)
        self.ntc_ereg.setText(str(self.engpar))
        self.ntc_ereg.setEnabled(False)
        hbox25.addWidget(self.ntc_ereg)

        vbox2.addLayout(hbox25)


        hbox24 = QtWidgets.QHBoxLayout()
        text1 = QtWidgets.QLabel(self)
        text1.setText('Sample Thickness')
        hbox24.addWidget(text1)
        hbox24.addStretch(1)
        self.ntc_samplethick = QtWidgets.QLineEdit(self)
        self.ntc_samplethick.setFixedWidth(65)
        self.ntc_samplethick.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.ntc_samplethick.setAlignment(QtCore.Qt.AlignRight)
        self.ntc_samplethick.setText(str(self.samplethick))
        hbox24.addWidget(self.ntc_samplethick)
        vbox2.addLayout(hbox24)

        self.cb_nneg = QtWidgets.QCheckBox('Non-Negativity Constraints', self)
        self.cb_nneg.setChecked(True)
        self.cb_nneg.stateChanged.connect(self.OnNonNegConst)
        vbox2.addWidget(self.cb_nneg)

        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.HLine)
        line.setFrameShadow(QtWidgets.QFrame.Sunken)
        vbox2.addStretch(1)
        vbox2.addWidget(line)
        vbox2.addStretch(1)

        vbox2.addWidget(QtWidgets.QLabel('Multiprocessing:'))

        hbox26 = QtWidgets.QHBoxLayout()
        text1 = QtWidgets.QLabel(self)
        text1.setText('Number of processors')
        hbox26.addWidget(text1)
        hbox26.addStretch(1)
        self.ntc_processors = QtWidgets.QLineEdit(self)
        self.ntc_processors.setFixedWidth(65)
        self.ntc_processors.setValidator(QtGui.QIntValidator(1, 1000, self))
        self.ntc_processors.setAlignment(QtCore.Qt.AlignRight)
        self.ntc_processors.setText(str(self.nprocessors))
        hbox26.addWidget(self.ntc_processors)
        vbox2.addLayout(hbox26)

        sizer2.setLayout(vbox2)

        #panel 3
        sizer3 = QtWidgets.QGroupBox('ROI')
        vbox3 = QtWidgets.QVBoxLayout()


        self.button_roi = QtWidgets.QPushButton( 'Select ROI')
        self.button_roi.clicked.connect(self.OnSelectROI)
        self.button_roi.setEnabled(False)
        vbox3.addWidget(self.button_roi)

        self.button_roihist = QtWidgets.QPushButton( 'Histogram ROI selection...')
        self.button_roihist.clicked.connect(self.OnROIHistogram)
        self.button_roihist.setEnabled(False)
        vbox3.addWidget(self.button_roihist)


        self.button_roispec = QtWidgets.QPushButton( 'Show ROI Spectrum')
        self.button_roispec.clicked.connect(self.OnShowROISpec)
        self.button_roispec.setEnabled(False)
        vbox3.addWidget(self.button_roispec)

        self.button_roidel = QtWidgets.QPushButton( 'Reset ROI')
        self.button_roidel.clicked.connect(self.OnResetROI)
        self.button_roidel.setEnabled(False)
        vbox3.addWidget(self.button_roidel)

        self.button_saveroi = QtWidgets.QPushButton( 'Save ROI as .mrc')
        self.button_saveroi.clicked.connect(self.OnSaveROI)
        self.button_saveroi.setEnabled(False)
        vbox3.addWidget(self.button_saveroi)

        self.button_loadroi = QtWidgets.QPushButton( 'Load ROI from .mrc')
        self.button_loadroi.clicked.connect(self.OnLoadROI)
        self.button_loadroi.setEnabled(False)
        vbox3.addWidget(self.button_loadroi)

        sizer3.setLayout(vbox3)





        #panel 5

        vbox5 = QtWidgets.QVBoxLayout()
        vbox5.addStretch(1)

        self.tc_imagecomp = QtWidgets.QLabel(self)
        self.tc_imagecomp.setText("Dataset: ")
        vbox5.addWidget(self.tc_imagecomp)


        gridsizer5 = QtWidgets.QGridLayout()
        gridsizer5.setSpacing(5)

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.absimgfig = Figure((PlotH*0.9, PlotH*0.9))

        self.AbsImagePanel = FigureCanvas(self.absimgfig)
        self.AbsImagePanel.setParent(self)
        self.cid1 = self.AbsImagePanel.mpl_connect('button_press_event', self.OnPointImage)


        fbox.addWidget(self.AbsImagePanel)
        frame.setLayout(fbox)
        gridsizer5.addWidget(frame, 1, 1, QtCore .Qt. AlignLeft)


        self.slider_slice = QtWidgets.QScrollBar(QtCore.Qt.Vertical)
        self.slider_slice.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_slice.valueChanged[int].connect(self.OnScrollSlice)
        self.slider_slice.setRange(0, 100)

        gridsizer5.addWidget(self.slider_slice, 1, 0, QtCore .Qt. AlignLeft)


        self.slider_comp = QtWidgets.QScrollBar(QtCore.Qt.Horizontal)
        self.slider_comp.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_comp.valueChanged[int].connect(self.OnScrollComp)
        self.slider_comp.setRange(0, 100)
        self.slider_comp.setEnabled(True)
        self.tc_comp = QtWidgets.QLabel(self)
        self.tc_comp.setText("Component: ")
        hbox51 = QtWidgets.QHBoxLayout()
        hbox51.addWidget(self.tc_comp)
        hbox51.addWidget(self.slider_comp)
        gridsizer5.addLayout(hbox51, 0, 1)



        frame2 = QtWidgets.QFrame()
        frame2.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox2 = QtWidgets.QHBoxLayout()

        self.absimgfig2 = Figure((PlotH*0.9, PlotH*0.9))

        self.AbsImagePanel2 = FigureCanvas(self.absimgfig2)
        self.AbsImagePanel2.setParent(self)


        fbox2.addWidget(self.AbsImagePanel2)
        frame2.setLayout(fbox2)


        gridsizer5.addWidget(frame2, 2, 1, QtCore .Qt. AlignLeft)


        frame3 = QtWidgets.QFrame()
        frame3.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox3 = QtWidgets.QHBoxLayout()

        self.absimgfig3 = Figure((PlotH*0.9, PlotH*0.9))

        self.AbsImagePanel3 = FigureCanvas(self.absimgfig3)
        self.AbsImagePanel3.setParent(self)


        fbox3.addWidget(self.AbsImagePanel3)
        frame3.setLayout(fbox3)
        gridsizer5.addWidget(frame3, 2, 2, QtCore .Qt. AlignLeft)

        gridsizer5.addWidget(sizer3, 1, 2, QtCore .Qt. AlignCenter)



        vbox5.addLayout(gridsizer5)
        vbox5.addStretch(1)


        vboxtop = QtWidgets.QVBoxLayout()

        hboxtop = QtWidgets.QHBoxLayout()
        vboxt1 = QtWidgets.QVBoxLayout()
        vboxt1.addStretch (1)
        vboxt1.addWidget(sizer1)
        vboxt1.addStretch (1)
        vboxt1.addWidget(sizer2)
        vboxt1.addStretch (1)
        #vboxt1.addWidget(sizer3)
        #vboxt1.addStretch (1)

        hboxtop.addStretch (5)
        hboxtop.addLayout(vboxt1)
        hboxtop.addStretch (5)
        hboxtop.addLayout(vbox5)
        hboxtop.addStretch (5)




        vboxtop.addStretch (5)
        vboxtop.addLayout(hboxtop)
        vboxtop.addStretch (9)

        vboxtop.setContentsMargins(20,20,20,20)
        self.setLayout(vboxtop)



    def OnLoadTomoEng(self, event):

        self.fulltomorecdata = []
        self.tomo_calculated = 0
        self.full_tomo_calculated = 0
        self.single_data_is_volume = False
        self.single_volume_data = None

        self.NewStackClear()

        self.button_save.setEnabled(False)
        self.tc_comp.setText('Component: ')
        self.slider_comp.setEnabled(False)


        self.tomodata = self.stack.od4d

        self.theta = self.stack.theta
        self.n_cols = self.stack.n_cols
        self.n_rows = self.stack.n_rows

        self.datanames = []
        for i in range(self.stack.n_ev):
            self.datanames.append(str(self.stack.ev[i]))

        self.ncomponents = self.stack.n_ev

        self.combonames.clear()
        self.combonames.addItems(self.datanames)

        self.tc_imagecomp.setText("Dataset: Energies")

        self.button_calcall.setEnabled(True)
        self.button_calc1.setEnabled(True)
        self.button_roi.setEnabled(False)
        self.energiesloaded = 1

        self.button_expdata.setEnabled(True)

        self.preview_loaded = True

        self.islice = 0
        self.slider_slice.setRange(0, self.tomodata.shape[3] - 1)
        self.slider_slice.setValue(self.islice)
        self.ShowImage()


#----------------------------------------------------------------------
    def OnLoadTomoComponents(self, event):

        self.fulltomorecdata = []
        self.tomo_calculated = 0
        self.full_tomo_calculated = 0
        self.single_data_is_volume = False
        self.single_volume_data = None

        self.NewStackClear()


        self.button_save.setEnabled(False)
        self.tc_comp.setText('Component: ')
        self.slider_comp.setEnabled(False)

        self.ncomponents = self.anlz.n_target_spectra

        self.tomodata = np.zeros((self.stack.n_cols, self.stack.n_rows, self.ncomponents, self.stack.n_theta))


        for i in range (self.stack.n_theta):
            if self.window().tab_spec.showraw == True:
                self.tomodata[:,:,:,i] = self.anlz.target_svd_maps4D[i]
            else:
                self.tomodata[:,:,:,i] = self.anlz.target_pcafit_maps[i]

        self.theta = self.stack.theta
        self.n_cols = self.stack.n_cols
        self.n_rows = self.stack.n_rows

        self.datanames = []
        for i in range(self.ncomponents):
            self.datanames.append(str(self.anlz.tspec_names[i]))


        self.combonames.clear()
        self.combonames.addItems(self.datanames)

        self.tc_imagecomp.setText("Dataset: Spectral Components")

        self.button_calcall.setEnabled(True)
        self.button_calc1.setEnabled(True)
        self.button_roi.setEnabled(False)
        self.energiesloaded = 0

        self.button_expdata.setEnabled(True)

        self.preview_loaded = True

        self.islice = 0
        self.slider_slice.setRange(0, self.tomodata.shape[3] - 1)
        self.slider_slice.setValue(self.islice)
        self.ShowImage()

#----------------------------------------------------------------------
    def OnLoadSingleMrc(self, event):


        wildcard = "Supported 4D formats (*.mrc *.ali *.ncb);;Mrc files (*.mrc *.ali);;NCB files (*.ncb);;"


        OpenFileName, _filter = QtWidgets.QFileDialog.getOpenFileName(self, 'Load Tomo Dataset', '', wildcard,
                                                         None, QtWidgets.QFileDialog.DontUseNativeDialog)

        OpenFileName = str(OpenFileName)

        if OpenFileName == '':
            return

        basename, extension = os.path.splitext(OpenFileName)
        extension = extension.lower()

        if extension in ('.mrc', '.ali'):

            data = tomo_reconstruction.load_mrc(OpenFileName)
            dims = data.shape
            self.single_volume_data = None

            # Read projection angles from an optional text file.
            wildcard = "Angle files (*.txt *.csv *.ang *.dat *.*);;"
            OpenFileName2, _filter = QtWidgets.QFileDialog.getOpenFileName(self, 'Load Angle data', '', wildcard,
                                                             None, QtWidgets.QFileDialog.DontUseNativeDialog)

            OpenFileName2 = str(OpenFileName2)
            tlist = []
            angle_error = None
            if OpenFileName2 != '':
                _ang_ext = os.path.splitext(OpenFileName2)[1].lower()
                if _ang_ext in ('.mrc', '.ali', '.ncb'):
                    angle_error = 'Selected angle file looks like tomography data, not text angles.'
                else:
                    try:
                        with open(OpenFileName2, 'r', encoding='utf-8') as f:
                            for line in f:
                                line = line.strip()
                                if line == '' or line.startswith('*') or line.startswith('#'):
                                    continue
                                # Support one value per line or comma-separated files.
                                tokens = line.replace(',', ' ').split()
                                for token in tokens:
                                    tlist.append(float(token))
                    except (UnicodeDecodeError, ValueError, OSError) as exc:
                        angle_error = str(exc)

            if len(tlist) == 0:
                # No angle file selected (or invalid/empty): assume this is a slice volume.
                self.single_data_is_volume = True
                # Saved MANTiS recon volumes are written as data.T, so restore canonical
                # display/reconstruction orientation (x, y, z) when no angles are provided.
                canonical_volume = np.asarray(data, dtype=np.float32).T.copy()
                self.single_volume_data = canonical_volume
                self.theta = np.arange(canonical_volume.shape[2], dtype=float)
                data = canonical_volume
                QtWidgets.QMessageBox.information(
                    self,
                    'Assuming slice volume',
                    'No valid angle file selected. Loading as slice volume.'
                )
                if angle_error is not None:
                    QtWidgets.QMessageBox.warning(
                        self,
                        'Invalid angle file',
                        f'Could not read angle file as text values:\n{angle_error}'
                    )
            else:
                self.theta = np.array(tlist)
                self.single_data_is_volume = False
                self.single_volume_data = None


            ntheta = len(self.theta)
            if ntheta != dims[2]:
                data = np.swapaxes(data, 0, 2)


            print('angle num:', ntheta)
            dims = data.shape
            print('Data shape', dims)
            self.n_cols = dims[0]
            self.n_rows = dims[1]
            self.tomodata = np.zeros((dims[0], dims[1], 1, dims[2]))
            self.tomodata[:,:,0,:] = data

        elif extension == '.ncb' :

            data, thetalist = file_ncb.read_ncb_data(self, OpenFileName)

            dims = data.shape
            self.theta = np.array(thetalist)

            ntheta = len(self.theta)
            if ntheta != dims[2]:
                data = np.swapaxes(data, 0, 2)


            dims = data.shape
            self.n_cols = dims[0]
            self.n_rows = dims[1]
            self.tomodata = np.zeros((dims[0], dims[1], 1, dims[2]))
            self.tomodata[:,:,0,:] = data
            self.single_data_is_volume = False
            self.single_volume_data = None

        else:
            QtWidgets.QMessageBox.warning(
                self,
                'Unsupported format',
                'Please choose a .mrc, .ali, or .ncb tomography dataset.'
            )
            return


        self.fulltomorecdata = []
        self.tomo_calculated = 0
        self.full_tomo_calculated = 0

        # NewStackClear resets mode flags; preserve the single-load interpretation.
        _single_data_is_volume = self.single_data_is_volume
        _single_volume_data = self.single_volume_data


        self.NewStackClear()

        self.single_data_is_volume = _single_data_is_volume
        self.single_volume_data = _single_volume_data

        self.button_save.setEnabled(False)
        self.tc_comp.setText('Component: ')
        self.slider_comp.setEnabled(False)


        basename = os.path.basename(str(OpenFileName))
        self.datanames = []


        self.datanames.append(str(basename))

        self.ncomponents = 1

        self.combonames.clear()
        #self.combonames.addItems(self.datanames)

        self.tc_imagecomp.setText("Dataset:" + OpenFileName)

        self.button_calcall.setEnabled(False)
        self.button_calc1.setEnabled(not self.single_data_is_volume)
        self.button_roi.setEnabled(False)
        self.energiesloaded = 0

        self.button_expdata.setEnabled(True)

        self.preview_loaded = True

        if self.single_data_is_volume and self.single_volume_data is not None:
            # Match normal lower-panel behavior by using the loaded volume as display tomorec.
            self.tr.tomorec = np.asarray(self.single_volume_data, dtype=np.float32).copy()
            self.tomo_calculated = 1
            self.full_tomo_calculated = 0
            self.nslices = self.tr.tomorec.shape[2]

        # Preview the loaded projections immediately; reconstruction can come later.
        self.islice = 0
        self.slider_slice.setRange(0, dims[2] - 1)
        self.slider_slice.setValue(self.islice)
        self.ShowImage()


#----------------------------------------------------------------------
    def OnCalcTomoFull(self, event):

        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))

        value = self.ntc_niterations.text()
        self.maxIters = int(value)
        value = self.ntc_beta.text()
        self.beta = float(value)
        value = self.ntc_samplethick.text()
        self.samplethick = int(value)
        self.nprocessors = int(self.ntc_processors.text())

        self.fulltomorecdata = []

        binningfactor = int (self.combobin.currentText())
        dims = self.tomodata[:,:,0,:].shape


        if (self.useengreg == 1) and (self.algo == 0):

            value = self.ntc_ereg.text()
            beta2 = float(value)

            print('Calculate initial reconstructions')

            initrecs = []
            for i in range(self.ncomponents):

                print('Progress ',i+1,' / ',self.ncomponents)

                if binningfactor > 1:
                    shape = (int(dims[0]/binningfactor), int(dims[1]/binningfactor))
                    projdata = np.zeros((shape[0], shape[1], dims[2]))
                    #print 'Binning factor:', binningfactor
                    #print 'Binned data dims', shape
                    for j in range(dims[2]):
                        projdata[:,:,j] = rebin(self.tomodata[0:shape[0]*binningfactor,0:shape[1]*binningfactor,i,j], shape)
                else:
                    projdata = self.tomodata[:,:,i,:]


                self.tr.calc_tomo(projdata,
                                   self.theta,
                                   self.maxIters,
                                   self.beta,
                                   0,
                                   algorithm = self.algo,
                                   nonnegconst = self.nonnegconst,
                                   nprocessors = self.nprocessors)

                initrecs.append(np.swapaxes(np.array(self.tr.tomorec.copy()), 0, 1))

            for i in range(self.ncomponents):

                print('Progress ',i+1,' / ',self.ncomponents)

                if binningfactor > 1:
                    shape = (int(dims[0]/binningfactor), int(dims[1]/binningfactor))
                    projdata = np.zeros((shape[0], shape[1], dims[2]))
                    #print 'Binning factor:', binningfactor
                    #print 'Binned data dims', shape
                    for j in range(dims[2]):
                        projdata[:,:,j] = rebin(self.tomodata[0:shape[0]*binningfactor,0:shape[1]*binningfactor,i,j], shape)
                else:
                    projdata = self.tomodata[:,:,i,:]

                self.tr.calc_tomo(projdata,
                                   self.theta,
                                   self.maxIters,
                                   self.beta,
                                   self.samplethick,
                                   algorithm = 2,
                                   x0=initrecs,
                                   comp = i, beta2=beta2,
                                   nonnegconst = self.nonnegconst,
                                   nprocessors = self.nprocessors)

                self.fulltomorecdata.append(self.tr.tomorec.copy())


        else:

            for i in range(self.ncomponents):

                print('Progress ',i+1,' / ',self.ncomponents)

                if binningfactor > 1:
                    shape = (int(dims[0]/binningfactor), int(dims[1]/binningfactor))
                    projdata = np.zeros((shape[0], shape[1], dims[2]))
                    #print 'Binning factor:', binningfactor
                    #print 'Binned data dims', shape
                    for j in range(dims[2]):
                        projdata[:,:,j] = rebin(self.tomodata[0:shape[0]*binningfactor,0:shape[1]*binningfactor,i,j], shape)
                else:
                    projdata = self.tomodata[:,:,i,:]

                self.tr.calc_tomo(projdata,
                                   self.theta,
                                   self.maxIters,
                                   self.beta,
                                   self.samplethick,
                                   algorithm = self.algo,
                                   nonnegconst = self.nonnegconst,
                                   nprocessors = self.nprocessors)

                self.fulltomorecdata.append(self.tr.tomorec.copy())



        self.tr.tomorec = self.fulltomorecdata[self.icomp]



        self.full_tomo_calculated = 1
        self.tomo_calculated = 1

        dims = self.tr.tomorec.shape
        self.n_cols = dims[0]
        self.n_rows = dims[1]

        self.nslices = dims[2]

        self.islice = int(dims[2]/2)

        self.slider_slice.setRange(0, dims[2]-1)

        self.ROIvol =  [[]]* dims[2]
        self.ROIarray = np.zeros((dims[0], dims[1], dims[2]))

        self.tc_comp.setText('Component: '+self.datanames[self.icomp])

        self.slider_slice.setValue(self.islice)
        self.slider_comp.setRange(0, self.ncomponents-1)
        self.slider_comp.setEnabled(True)
        self.slider_comp.setValue(self.icomp)
        self.button_save.setEnabled(True)
        self.button_saveall.setEnabled(True)
        self.button_roi.setEnabled(True)
        self.button_roihist.setEnabled(True)
        self.button_loadroi.setEnabled(True)


        self.ShowImage()

        QtWidgets.QApplication.restoreOverrideCursor()

#----------------------------------------------------------------------
    def OnCalcTomo1(self, event):
        if self.single_data_is_volume:
            QtWidgets.QMessageBox.information(
                self,
                'Slice volume mode',
                'This dataset was loaded without an angle file. Slices and derived projections are shown automatically; calculation is not required.'
            )
            return

        if self.tomodata is None:
            QtWidgets.QMessageBox.warning(self, 'No tomo data loaded',
                                          'Load a tomography dataset before calculating.')
            return

        if self.tomodata.ndim != 4 or self.tomodata.shape[2] <= 0:
            QtWidgets.QMessageBox.warning(self, 'Invalid tomo data',
                                          'Loaded tomography data has an unexpected shape.')
            return

        if len(self.theta) != self.tomodata.shape[3]:
            QtWidgets.QMessageBox.warning(self, 'Angle mismatch',
                                          'Number of projection angles does not match dataset.')
            return

        if self.select1 < 0 or self.select1 >= self.tomodata.shape[2]:
            self.select1 = 0

        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
        try:
            value = self.ntc_niterations.text()
            self.maxIters = int(value)
            value = self.ntc_beta.text()
            self.beta = float(value)
            value = self.ntc_samplethick.text()
            self.samplethick = int(value)
            self.nprocessors = int(self.ntc_processors.text())

            dims = self.tomodata[:,:,self.select1,:].shape
            print('Data dims = ', dims)

            projdata = self.tomodata[:,:,self.select1,:]

            binningfactor = int(self.combobin.currentText())

            if binningfactor > 1:
                shape = (int(dims[0]/binningfactor), int(dims[1]/binningfactor))
                projdata = np.zeros((shape[0], shape[1], dims[2]))
                print('Binning factor:', binningfactor)
                print('Binned data dims', shape)
                for i in range(dims[2]):
                    projdata[:,:,i] = rebin(
                        self.tomodata[0:shape[0]*binningfactor, 0:shape[1]*binningfactor, self.select1, i],
                        shape
                    )

            self.tr.calc_tomo(projdata,
                               self.theta,
                               self.maxIters,
                               self.beta,
                               self.samplethick,
                               algorithm = self.algo,
                               nonnegconst=self.nonnegconst,
                               nprocessors = self.nprocessors)

            self.tomo_calculated = 1
            self.full_tomo_calculated = 0

            dims = self.tr.tomorec.shape
            self.nslices = dims[2]
            self.n_cols = dims[0]
            self.n_rows = dims[1]

            self.ROIvol = [[]]* dims[2]
            self.ROIarray = np.zeros((dims[0], dims[1], dims[2]))
            self.button_roispec.setEnabled(False)

            self.islice = int(dims[2]/2)
            self.slider_slice.setValue(self.islice)
            self.slider_slice.setRange(0, dims[2]-1)

            self.tc_comp.setText('Component: '+self.datanames[self.select1])
            self.slider_comp.setRange(0, 0)
            self.button_roi.setEnabled(True)
            self.button_save.setEnabled(True)
            self.button_roispec.setEnabled(False)
            self.button_roi.setEnabled(True)
            self.button_loadroi.setEnabled(True)
            self.button_roihist.setEnabled(True)

            self.ShowImage()

        except ValueError:
            QtWidgets.QMessageBox.warning(self, 'Invalid reconstruction parameters',
                                          'Iterations, beta, thickness, and processors must be valid numbers.')
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()


#----------------------------------------------------------------------
    def OnSelect1Comp(self, value):
        item = value
        self.select1 = item


#----------------------------------------------------------------------
    def OnSelectAlgo(self, value):
        item = value
        self.algo = item

        # 0 - CS
        if self.algo == 0:
            self.tc_par.setEnabled(True)
            self.ntc_beta.setEnabled(True)
            self.cb_ereg.setEnabled(True)
            self.ntc_ereg.setEnabled(True)
        else:
            self.tc_par.setEnabled(False)
            self.ntc_beta.setEnabled(False)
            self.cb_ereg.setEnabled(False)
            self.ntc_ereg.setEnabled(False)


 #----------------------------------------------------------------------
    def OnCBEngReg(self, state):

        if state == QtCore.Qt.Checked:
            self.useengreg = 1
            self.ntc_ereg.setEnabled(True)
        else:
            self.useengreg = 0
            self.ntc_ereg.setEnabled(False)

 #----------------------------------------------------------------------
    def OnNonNegConst(self, state):

        if state == QtCore.Qt.Checked:
            self.nonnegconst = 1
        else:
            self.nonnegconst = 0


#----------------------------------------------------------------------
    def OnScrollSlice(self, value):
        self.islice = value

        self.ShowImage()


#----------------------------------------------------------------------
    def OnScrollComp(self, value):
        self.icomp = value

        if self.full_tomo_calculated == 0:
            return

        self.tr.tomorec = self.fulltomorecdata[self.icomp]

        self.tc_comp.setText('Component: '+self.datanames[self.icomp])

        self.ShowImage()


#----------------------------------------------------------------------
    def OnPointImage(self, evt):


        if self.energiesloaded == 0 or self.full_tomo_calculated == 0:
            return

        x = evt.xdata
        y = evt.ydata


        if (x == None) or (y == None):
            return


        ix = int(np.floor(x))
        iy = self.n_rows-1-int(np.floor(y))

        if ix<0 :
            ix=0
        if ix>self.n_cols-1 :
            ix=self.n_cols-1
        if iy<0 :
            iy=0
        if iy>self.n_rows-1 :
            iy=self.n_rows-1

        spectrum = []

        for i in range(self.ncomponents):
            spectrum.append(self.fulltomorecdata[i][ix,iy,self.islice])

        title = 'Point [{0:d}, {0:d}, {0:d}]'.format(ix,iy,self.islice)

        plot = PlotFrame(self, self.stack.ev, spectrum, title=title)
        plot.show()

#----------------------------------------------------------------------
    def OnExportData(self, event):


        wildcard = "Mrc files (*.mrc);;"

        SaveFileName, _filter = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Tomo Reconstructions', '', wildcard)

        SaveFileName = str(SaveFileName)
        if SaveFileName == '':
            return

        if os.path.splitext(SaveFileName)[1] == '':
            SaveFileName += '.mrc'


        basename, extension = os.path.splitext(SaveFileName)

        for i in range(self.ncomponents):
            data = self.tomodata[:,:,i,:]

            savefn = basename + '_TiltS_'+self.datanames[i]+extension

            self.tr.save_mrc(savefn, data)

        savefn2 = basename + '_Angles.txt'

        f = open(str(savefn2),'wt')

        for i in range(len(self.theta)):
            print(self.theta[i], file=f)

        f.close()


#----------------------------------------------------------------------
    def OnSave(self, event):


        wildcard = "Mrc files (*.mrc);;"

        SaveFileName, _filter = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Tomo Reconstructions', '', wildcard)

        SaveFileName = str(SaveFileName)
        if SaveFileName == '':
            return

        if os.path.splitext(SaveFileName)[1] == '':
            SaveFileName += '.mrc'


        data = self.tr.tomorec

        self.tr.save_mrc(SaveFileName, data.T)


#----------------------------------------------------------------------
    def OnSaveAll(self, event):


        wildcard = "Mrc files (*.mrc);;"

        SaveFileName, _filter = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Tomo Reconstructions', '', wildcard)

        SaveFileName = str(SaveFileName)
        if SaveFileName == '':
            return

        basename, extension = os.path.splitext(SaveFileName)
        if extension == '':
            extension = '.mrc'



        for i in range(self.ncomponents):
            data = self.fulltomorecdata[i]

            savefn = basename + '_'+self.datanames[i]+extension

            self.tr.save_mrc(savefn, data.T)


#----------------------------------------------------------------------
    def OnSaveROI(self, event):


        wildcard = "Mrc files (*.mrc);;"

        SaveFileName, _filter = QtWidgets.QFileDialog.getSaveFileName(self, 'Save ROI Selection', '', wildcard)

        SaveFileName = str(SaveFileName)
        if SaveFileName == '':
            return

        if os.path.splitext(SaveFileName)[1] == '':
            SaveFileName += '.mrc'


        data = self.ROIarray

        self.tr.save_mrc(SaveFileName, data)


#----------------------------------------------------------------------
    def OnLoadROI(self, event):


        wildcard = "Mrc files (*.mrc);;"

        OpenFileName, _filter = QtWidgets.QFileDialog.getOpenFileName(self, 'Load ROI Selection', '', wildcard)

        OpenFileName = str(OpenFileName)
        if OpenFileName == '':
            return


        self.ROIarray = tomo_reconstruction.load_mrc(OpenFileName)


        for i in range(self.nslices):

            ROIpix = np.ma.array(self.ROIarray[:,:,i])

            ROIpix_masked =  np.ma.masked_values(ROIpix, 0)
            self.ROIvol[i] = ROIpix_masked


        self.haveROI = 1

        if self.full_tomo_calculated == 1:
            self.button_roispec.setEnabled(True)

        self.button_roidel.setEnabled(True)
        self.button_saveroi.setEnabled(True)

        self.ShowImage()

#----------------------------------------------------------------------

    def OnSelectLasso(self,verts):
#         self.lasso.disconnect_events()
#         path = matplotlib.path.Path(verts)
#         self.ind = np.nonzero(path.contains_points(self.xys))[0]
#         print 'Selected '+str(len(self.ind))+' points'
#
#         indices = path.contains_points(self.xys)
#
#         mask = np.array([ path.contains_point((j,i)) for i in range(self.stack.n_rows) for j in range(self.stack.n_cols)]).reshape(self.stack.n_cols,self.stack.n_rows)
#
#         self.ROIvol[:,:,self.islice] = mask
#
#         print np.sum(self.ROIvol)

        path = matplotlib.path.Path(verts)
        #find pixels inside the polygon

        ROIpix = np.zeros((self.n_cols,self.n_rows))

        for i in range(self.n_cols):
            for j in range(self.n_rows):
                Pinside = path.contains_point((i,j))
                if Pinside == True:
                    ROIpix[i, self.n_rows-1-j] = 255



        ROIpix = np.ma.array(ROIpix)

        ROIpix_masked =  np.ma.masked_values(ROIpix, 0)


        self.ROIvol[self.islice] = ROIpix_masked
        self.ROIarray[:,:,self.islice] = ROIpix[:]

        self.haveROI = 1

        self.ShowImage()



#----------------------------------------------------------------------
    def OnSelectROI(self, event):

        self.AbsImagePanel.mpl_disconnect(self.cid1)
        if self.full_tomo_calculated == 1 and self.energiesloaded == 1:
            self.button_roispec.setEnabled(True)

        self.button_roidel.setEnabled(True)
        self.button_roihist.setEnabled(True)
        self.button_saveroi.setEnabled(True)

        lineprops = dict(color='red', linestyle='-', linewidth = 1, alpha=1)


        self.lasso = LassoSelector(self.axes, onselect=self.OnSelectLasso, useblit=False, props=lineprops)

#----------------------------------------------------------------------
    def OnROIHistogram(self, event):
        #self.window().Hide()
        image = self.tr.tomorec[:,:,self.islice].copy()
        n_cols, n_rows = image.shape
        histogram = ROIHistogram(self, image, n_cols, n_rows)
        histogram.show()

#----------------------------------------------------------------------
    def OnResetROI(self, event):

        self.button_roispec.setEnabled(False)
        self.button_roidel.setEnabled(False)
        self.button_roihist.setEnabled(True)
        self.button_saveroi.setEnabled(False)


        self.ROIvol =  [[]]*self.nslices

        self.haveROI = 0

        self.cid1 = self.AbsImagePanel.mpl_connect('button_press_event', self.OnPointImage)

        self.ShowImage()

#----------------------------------------------------------------------
    def CalcROISpectrum(self):

        ROIspectrum = np.zeros((self.stack.n_ev))

        for ie in range(self.stack.n_ev):
            indices = np.where(self.ROIarray == 255)
            numroipix = self.ROIarray[indices].shape[0]
            if numroipix > 0:
                roivoxels = self.fulltomorecdata[ie]
                ROIspectrum[ie] = np.sum(roivoxels[indices])/numroipix

        return ROIspectrum


#----------------------------------------------------------------------
    def OnShowROISpec(self):

        spectrum = self.CalcROISpectrum()

        title = 'ROI spectrum'

        plot = PlotFrame(self, self.stack.ev, spectrum, title=title)
        plot.show()

#----------------------------------------------------------------------
    def ShowImage(self):

        if self.tomo_calculated == 0:
            if not self.preview_loaded or self.tomodata is None:
                return
            if self.single_data_is_volume and self.single_volume_data is not None:
                angle = min(max(self.islice, 0), self.single_volume_data.shape[2] - 1)
                image = self.single_volume_data[:, :, angle].copy()
            else:
                comp = min(max(self.select1, 0), self.tomodata.shape[2] - 1)
                angle = min(max(self.islice, 0), self.tomodata.shape[3] - 1)
                image = self.tomodata[:,:,comp,angle].copy()
        else:
            if self.single_data_is_volume and self.single_volume_data is not None:
                angle = min(max(self.islice, 0), self.single_volume_data.shape[2] - 1)
                image = self.single_volume_data[:, :, angle].copy()
            else:
                image = self.tr.tomorec[:,:,self.islice].copy()

        fig = self.absimgfig
        fig.clf()
        fig.add_axes(((0.0,0.0,1.0,1.0)))
        axes = fig.gca()
        fig.patch.set_alpha(1.0)

        im = axes.imshow(np.rot90(image), cmap=matplotlib.colormaps["gray"])

        if self.haveROI == 1:
            if len(self.ROIvol[self.islice]) > 0:
                im_red = axes.imshow(np.rot90(self.ROIvol[self.islice]), cmap=matplotlib.colormaps["autumn"])

#         if self.window().tab_prep.show_scale_bar == 1:
#             #Show Scale Bar
#             if self.com.white_scale_bar == 1:
#                 sbcolor = 'white'
#             else:
#                 sbcolor = 'black'
#             startx = int(self.stk.n_cols*0.05)
#             starty = self.stk.n_rows-int(self.stk.n_rows*0.05)-self.stk.scale_bar_pixels_y
#             um_string = ' $\mathrm{\mu m}$'
#             microns = '$'+self.stk.scale_bar_string+' $'+um_string
#             axes.text(self.stk.scale_bar_pixels_x+startx+1,starty+1, microns, horizontalalignment='left', verticalalignment='center',
#                       color = sbcolor, fontsize=14)
#             #Matplotlib has flipped scales so I'm using rows instead of cols!
#             p = matplotlib.patches.Rectangle((startx,starty), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
#                                    color = sbcolor, fill = True)
#             axes.add_patch(p)


        axes.axis("off")
        self.AbsImagePanel.draw()
        self.axes = axes

        self.xys = np.dstack(np.meshgrid(np.arange(self.n_cols), np.arange(self.n_rows))).reshape(-1,2)


        if self.tomo_calculated == 0:
            # No reconstruction yet: keep the secondary panels clear.
            self.absimgfig2.clf()
            self.AbsImagePanel2.draw()
            self.absimgfig3.clf()
            self.AbsImagePanel3.draw()
            return


        #Show orthogonal slices
        dims = self.tr.tomorec.shape
        image2 = self.tr.tomorec[:,int(dims[1]/2),:]

        fig = self.absimgfig2
        fig.clf()
        fig.add_axes(((0.0,0.0,1.0,1.0)))
        axes2 = fig.gca()
        fig.patch.set_alpha(1.0)

        im = axes2.imshow(np.rot90(image2), cmap=matplotlib.colormaps["gray"])

        axes2.axis("off")
        self.AbsImagePanel2.draw()


        image3 = self.tr.tomorec[int(dims[0]/2),:,:].T

        fig = self.absimgfig3
        fig.clf()
        fig.add_axes(((0.0,0.0,1.0,1.0)))
        axes3 = fig.gca()
        fig.patch.set_alpha(1.0)

        im = axes3.imshow(np.rot90(image3), cmap=matplotlib.colormaps["gray"])

        axes3.axis("off")
        self.AbsImagePanel3.draw()


#----------------------------------------------------------------------
    def MakeHistogramROI(self, histmin, histmax):

        # Keep ROI containers aligned with current reconstruction dimensions.
        rec_dims = self.tr.tomorec.shape
        if np.shape(self.ROIarray) != rec_dims:
            self.ROIarray = np.zeros(rec_dims)
            self.ROIvol = [[]] * rec_dims[2]
            self.nslices = rec_dims[2]

        self.n_cols = rec_dims[0]
        self.n_rows = rec_dims[1]

        for i in range(self.nslices):

            slice_data = self.tr.tomorec[:,:,i]
            hist_indices = np.where((histmin < slice_data) & (slice_data < histmax))

            ROIpix = np.zeros(slice_data.shape)
            ROIpix[hist_indices] = 255

            ROIpix = np.ma.array(ROIpix)

            ROIpix_masked = np.ma.masked_values(ROIpix, 0)


            self.ROIvol[i] = ROIpix_masked
            self.ROIarray[:,:,i] = ROIpix[:]

        self.haveROI = 1

        if self.full_tomo_calculated == 1:
            self.button_roispec.setEnabled(True)

        self.button_roidel.setEnabled(True)
        self.button_roihist.setEnabled(True)
        self.button_saveroi.setEnabled(True)

        self.ShowImage()


#----------------------------------------------------------------------
    def NewStackClear(self):

        self.tr = tomo_reconstruction.Ctomo(self.stack)

        self.fulltomorecdata = []
        self.ncomponents = 0

        self.tomo_calculated = 0
        self.full_tomo_calculated = 0
        self.energiesloaded = 0
        self.single_data_is_volume = False
        self.single_volume_data = None

        self.datanames = []

        self.ROIvol = []
        self.ROIarray = []
        self.haveROI = 0


        fig = self.absimgfig
        fig.clf()
        self.AbsImagePanel.draw()

        fig = self.absimgfig2
        fig.clf()
        self.AbsImagePanel2.draw()


        fig = self.absimgfig3
        fig.clf()
        self.AbsImagePanel3.draw()

        self.button_spcomp.setEnabled(False)
        self.button_engdata.setEnabled(False)
        self.button_expdata.setEnabled(False)
        self.button_calc1.setEnabled(False)
        self.button_calcall.setEnabled(False)
        self.button_save.setEnabled(False)
        self.button_roi.setEnabled(False)
        self.button_roidel.setEnabled(False)
        self.button_roihist.setEnabled(False)
        self.button_saveroi.setEnabled(False)
        self.button_roispec.setEnabled(False)

        self.tc_imagecomp.setText("Dataset: ")
        self.tc_comp.setText('Component: ')

        self.slider_comp.setEnabled(False)

        self.combonames.clear()
    # ... (Truncated here for brevity in this extraction step, assuming other methods are copied by valid implementation logic)
