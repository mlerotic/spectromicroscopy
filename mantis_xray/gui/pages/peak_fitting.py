
import os
import numpy as np
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtCore import Qt
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from ..dialogs.fit_params import FitParams
from ...core.constants import PlotH, PlotW

class PageXrayPeakFitting(QtWidgets.QWidget):
    def __init__(self, common, data_struct, stack, anlz):
        super(PageXrayPeakFitting, self).__init__()

        self.initUI(common, data_struct, stack, anlz)

#----------------------------------------------------------------------
    def initUI(self, common, data_struct, stack, anlz):


        self.com = common
        self.data_struct = data_struct
        self.stk = stack
        self.anlz = anlz

        self.i_spec = 1

        self.base = 0.0
        self.stepfitparams = np.zeros((6))
        self.gauss_fp_a = np.zeros((12))
        self.gauss_fp_m = np.zeros((12))
        self.gauss_fp_s = np.zeros((12))

        self.nsteps = []
        self.npeaks = []
        self.spectrumfitted = []
        self.initdone = []
        self.firstinit = []

        self.fits = []
        self.fits_sep = []

        self.showevlines = 0

        self.plotcolors = [(238/255.,42/255.,36/255.), (15/255.,135/255.,15/255.), (190/255.,0,190/255.), (0,0,190/255.),(190/255.,190/255.,0), (0,88/255.,0),
                            (132/255.,132/255.,255/255.), (158/255.,79/255.,70/255.), (0,132/255.,150/255.), (150/255.,211/255.,79/255.),(247/255.,158/255.,220/255.), (211/255.,17/255.,255/255.)]

        self.peak_engs = [285.0, 286.5, 288.5, 289.5, 290.5]
        self.peak_names = ['aromatic', 'phenolic', 'carboxyl', 'alkyl', 'carbonyl']

        vboxL = QtWidgets.QVBoxLayout()
        vboxR = QtWidgets.QVBoxLayout()
        hboxT = QtWidgets.QHBoxLayout()


        #panel
        sizer = QtWidgets.QGroupBox('Fit Parameters')
        vbox = QtWidgets.QVBoxLayout()

        font = QtGui.QFont()
        font.setPointSize(7)

        fgs1 = QtWidgets.QGridLayout()


        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Position [eV]')
        textctrl.setFont(font)
        fgs1.addWidget(textctrl, 0, 1)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Height [OD]')
        textctrl.setFont(font)
        fgs1.addWidget(textctrl, 0, 2)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Sigma [eV]')
        textctrl.setFont(font)
        fgs1.addWidget(textctrl, 0, 3)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Base')
        textctrl.setFont(font)
        fgs1.addWidget(textctrl, 16, 0)


        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Step 1')
        fgs1.addWidget(textctrl, 1, 0)

        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Step 2')
        fgs1.addWidget(textctrl, 2, 0)

        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Height [OD]')
        textctrl.setFont(font)
        fgs1.addWidget(textctrl, 3, 2)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Position [eV]')
        textctrl.setFont(font)
        fgs1.addWidget(textctrl, 3, 1)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Sigma [eV]')
        textctrl.setFont(font)
        fgs1.addWidget(textctrl, 3, 3)


        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 1')
        fgs1.addWidget(textctrl, 4, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 2')
        fgs1.addWidget(textctrl, 5, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 3')
        fgs1.addWidget(textctrl, 6, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 4')
        fgs1.addWidget(textctrl, 7, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 5')
        fgs1.addWidget(textctrl, 8, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 6')
        fgs1.addWidget(textctrl, 9, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 7')
        fgs1.addWidget(textctrl, 10, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 8')
        fgs1.addWidget(textctrl, 11, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 9')
        fgs1.addWidget(textctrl, 12, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 10')
        fgs1.addWidget(textctrl, 13, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 11')
        fgs1.addWidget(textctrl, 14, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 12')
        fgs1.addWidget(textctrl, 15, 0)


        self.le_sa1 = QtWidgets.QLineEdit(self)
        self.le_sa1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sa1.setAlignment(QtCore.Qt.AlignRight)
        self.le_sa1.setMaximumWidth(65)
        self.le_sa1.setText(str(self.stepfitparams[0]))
        fgs1.addWidget(self.le_sa1, 1, 1)

        self.le_sp1 = QtWidgets.QLineEdit(self)
        self.le_sp1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sp1.setAlignment(QtCore.Qt.AlignRight)
        self.le_sp1.setMaximumWidth(65)
        self.le_sp1.setText(str(self.stepfitparams[1]))
        fgs1.addWidget(self.le_sp1, 1, 2)

        self.le_sf1 = QtWidgets.QLineEdit(self)
        self.le_sf1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sf1.setAlignment(QtCore.Qt.AlignRight)
        self.le_sf1.setMaximumWidth(65)
        self.le_sf1.setText(str(self.stepfitparams[2]))
        fgs1.addWidget(self.le_sf1, 1, 3)

        self.le_base = QtWidgets.QLineEdit(self)
        self.le_base.setValidator(QtGui.QDoubleValidator(-99999, 99999, 2, self))
        self.le_base.setAlignment(QtCore.Qt.AlignRight)
        self.le_base.setMaximumWidth(65)
        self.le_base.setText(str(self.base))
        fgs1.addWidget(self.le_base, 16, 1)


        self.le_sa2 = QtWidgets.QLineEdit(self)
        self.le_sa2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sa2.setAlignment(QtCore.Qt.AlignRight)
        self.le_sa2.setMaximumWidth(65)
        self.le_sa2.setText(str(self.stepfitparams[3]))
        fgs1.addWidget(self.le_sa2, 2, 1)

        self.le_sp2 = QtWidgets.QLineEdit(self)
        self.le_sp2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sp2.setAlignment(QtCore.Qt.AlignRight)
        self.le_sp2.setMaximumWidth(65)
        self.le_sp2.setText(str(self.stepfitparams[4]))
        fgs1.addWidget(self.le_sp2, 2, 2)

        self.le_sf2 = QtWidgets.QLineEdit(self)
        self.le_sf2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sf2.setAlignment(QtCore.Qt.AlignRight)
        self.le_sf2.setMaximumWidth(65)
        self.le_sf2.setText(str(self.stepfitparams[5]))
        fgs1.addWidget(self.le_sf2, 2, 3)



        self.le_amplitudes = []
        self.le_mu = []
        self.le_sigma = []

        self.le_ga1 = QtWidgets.QLineEdit(self)
        self.le_ga1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga1.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga1.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga1, 4, 2)
        self.le_amplitudes.append(self.le_ga1)

        self.le_gm1 = QtWidgets.QLineEdit(self)
        self.le_gm1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm1.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm1.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm1, 4, 1)
        self.le_mu.append(self.le_gm1)

        self.le_gs1 = QtWidgets.QLineEdit(self)
        self.le_gs1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs1.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs1.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs1, 4, 3)
        self.le_sigma.append(self.le_gs1)

        self.le_ga2 = QtWidgets.QLineEdit(self)
        self.le_ga2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga2.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga2.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga2, 5, 2)
        self.le_amplitudes.append(self.le_ga2)

        self.le_gm2 = QtWidgets.QLineEdit(self)
        self.le_gm2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm2.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm2.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm2, 5, 1)
        self.le_mu.append(self.le_gm2)

        self.le_gs2 = QtWidgets.QLineEdit(self)
        self.le_gs2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs2.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs2.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs2, 5, 3)
        self.le_sigma.append(self.le_gs2)

        self.le_ga3 = QtWidgets.QLineEdit(self)
        self.le_ga3.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga3.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga3.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga3, 6, 2)
        self.le_amplitudes.append(self.le_ga3)

        self.le_gm3 = QtWidgets.QLineEdit(self)
        self.le_gm3.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm3.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm3.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm3, 6, 1)
        self.le_mu.append(self.le_gm3)

        self.le_gs3 = QtWidgets.QLineEdit(self)
        self.le_gs3.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs3.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs3.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs3, 6, 3)
        self.le_sigma.append(self.le_gs3)

        self.le_ga4 = QtWidgets.QLineEdit(self)
        self.le_ga4.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga4.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga4.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga4, 7, 2)
        self.le_amplitudes.append(self.le_ga4)

        self.le_gm4 = QtWidgets.QLineEdit(self)
        self.le_gm4.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm4.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm4.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm4, 7, 1)
        self.le_mu.append(self.le_gm4)

        self.le_gs4 = QtWidgets.QLineEdit(self)
        self.le_gs4.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs4.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs4.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs4, 7, 3)
        self.le_sigma.append(self.le_gs4)

        self.le_ga5 = QtWidgets.QLineEdit(self)
        self.le_ga5.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga5.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga5.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga5, 8, 2)
        self.le_amplitudes.append(self.le_ga5)

        self.le_gm5 = QtWidgets.QLineEdit(self)
        self.le_gm5.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm5.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm5.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm5, 8, 1)
        self.le_mu.append(self.le_gm5)

        self.le_gs5 = QtWidgets.QLineEdit(self)
        self.le_gs5.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs5.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs5.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs5, 8, 3)
        self.le_sigma.append(self.le_gs5)

        self.le_ga6 = QtWidgets.QLineEdit(self)
        self.le_ga6.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga6.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga6.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga6, 9, 2)
        self.le_amplitudes.append(self.le_ga6)

        self.le_gm6 = QtWidgets.QLineEdit(self)
        self.le_gm6.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm6.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm6.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm6, 9, 1)
        self.le_mu.append(self.le_gm6)

        self.le_gs6 = QtWidgets.QLineEdit(self)
        self.le_gs6.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs6.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs6.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs6, 9, 3)
        self.le_sigma.append(self.le_gs6)

        self.le_ga7 = QtWidgets.QLineEdit(self)
        self.le_ga7.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga7.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga7.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga7, 10, 2)
        self.le_amplitudes.append(self.le_ga7)

        self.le_gm7 = QtWidgets.QLineEdit(self)
        self.le_gm7.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm7.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm7.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm7, 10, 1)
        self.le_mu.append(self.le_gm7)

        self.le_gs7 = QtWidgets.QLineEdit(self)
        self.le_gs7.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs7.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs7.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs7, 10, 3)
        self.le_sigma.append(self.le_gs7)

        self.le_ga8 = QtWidgets.QLineEdit(self)
        self.le_ga8.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga8.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga8.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga8, 11, 2)
        self.le_amplitudes.append(self.le_ga8)

        self.le_gm8 = QtWidgets.QLineEdit(self)
        self.le_gm8.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm8.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm8.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm8, 11, 1)
        self.le_mu.append(self.le_gm8)

        self.le_gs8 = QtWidgets.QLineEdit(self)
        self.le_gs8.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs8.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs8.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs8, 11, 3)
        self.le_sigma.append(self.le_gs8)

        self.le_ga9 = QtWidgets.QLineEdit(self)
        self.le_ga9.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga9.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga9.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga9, 12, 2)
        self.le_amplitudes.append(self.le_ga9)

        self.le_gm9 = QtWidgets.QLineEdit(self)
        self.le_gm9.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm9.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm9.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm9, 12, 1)
        self.le_mu.append(self.le_gm9)

        self.le_gs9 = QtWidgets.QLineEdit(self)
        self.le_gs9.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs9.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs9.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs9, 12, 3)
        self.le_sigma.append(self.le_gs9)

        self.le_ga10 = QtWidgets.QLineEdit(self)
        self.le_ga10.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga10.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga10.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga10, 13, 2)
        self.le_amplitudes.append(self.le_ga10)

        self.le_gm10 = QtWidgets.QLineEdit(self)
        self.le_gm10.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm10.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm10.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm10, 13, 1)
        self.le_mu.append(self.le_gm10)

        self.le_gs10 = QtWidgets.QLineEdit(self)
        self.le_gs10.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs10.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs10.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs10, 13, 3)
        self.le_sigma.append(self.le_gs10)

        self.le_ga11 = QtWidgets.QLineEdit(self)
        self.le_ga11.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga11.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga11.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga11, 14, 2)
        self.le_amplitudes.append(self.le_ga11)

        self.le_gm11 = QtWidgets.QLineEdit(self)
        self.le_gm11.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm11.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm11.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm11, 14, 1)
        self.le_mu.append(self.le_gm11)

        self.le_gs11 = QtWidgets.QLineEdit(self)
        self.le_gs11.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs11.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs11.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs11, 14, 3)
        self.le_sigma.append(self.le_gs11)

        self.le_ga12 = QtWidgets.QLineEdit(self)
        self.le_ga12.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga12.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga12.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga12, 15, 2)
        self.le_amplitudes.append(self.le_ga12)

        self.le_gm12 = QtWidgets.QLineEdit(self)
        self.le_gm12.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm12.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm12.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm12, 15, 1)
        self.le_mu.append(self.le_gm12)

        self.le_gs12 = QtWidgets.QLineEdit(self)
        self.le_gs12.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs12.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs12.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs12, 15, 3)
        self.le_sigma.append(self.le_gs12)

        for i in range(12):
            self.le_amplitudes[i].setText(str(self.gauss_fp_a[i]))
            self.le_mu[i].setText(str(self.gauss_fp_m[i]))
            self.le_sigma[i].setText(str(self.gauss_fp_s[i]))

        vbox.addLayout(fgs1)
        sizer.setLayout(vbox)



        #panel 2
        vbox2 = QtWidgets.QVBoxLayout()

        self.tc_spec = QtWidgets.QLabel(self)
        self.tc_spec.setText("X-ray Spectrum: ")
        hbox11 = QtWidgets.QHBoxLayout()

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()
        self.Specfig = Figure((PlotW*1.25, PlotH*1.25))
        self.SpectrumPanel = FigureCanvas(self.Specfig)
        fbox.addWidget(self.SpectrumPanel)
        frame.setLayout(fbox)

        self.slider_spec = QtWidgets.QScrollBar(QtCore.Qt.Vertical)
        self.slider_spec.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_spec.setEnabled(False)
        self.slider_spec.valueChanged[int].connect(self.OnSScroll)
        self.slider_spec.setRange(1, 5)

        hbox11.addWidget(frame)
        hbox11.addWidget(self.slider_spec)

        vbox2.addWidget(self.tc_spec)
        vbox2.addLayout(hbox11)



        #panel 3
        sizer3 = QtWidgets.QGroupBox('Load spectra')
        vbox3 = QtWidgets.QVBoxLayout()

        self.button_addclspec = QtWidgets.QPushButton('Add cluster spectra')
        self.button_addclspec.clicked.connect( self.OnAddClusterSpectra)
        self.button_addclspec.setEnabled(False)
        vbox3.addWidget(self.button_addclspec)

        self.button_loadtspec = QtWidgets.QPushButton('Load spectra')
        self.button_loadtspec.clicked.connect(self.OnSpecFromFile)
        vbox3.addWidget(self.button_loadtspec)

        self.button_save = QtWidgets.QPushButton('Save images...')
        self.button_save.clicked.connect( self.OnSave)
        self.button_save.setEnabled(False)
        vbox3.addWidget(self.button_save)

        sizer3.setLayout(vbox3)



        #panel 4
        sizer4 = QtWidgets.QGroupBox('Fit settings')
        vbox4 = QtWidgets.QVBoxLayout()


        hbox41 = QtWidgets.QHBoxLayout()
        text1 = QtWidgets.QLabel(self)
        text1.setText('Number of steps')
        hbox41.addWidget(text1)
        self.nstepsspin = QtWidgets.QSpinBox()
        self.nstepsspin.setMinimumWidth(85)
        self.nstepsspin.setRange(0,2)
        self.nstepsspin.setValue(0)
        self.nstepsspin.valueChanged[int].connect(self.OnNstepsspin)
        hbox41.addWidget(text1)
        hbox41.addWidget(self.nstepsspin)
        vbox4.addLayout(hbox41)

        hbox42 = QtWidgets.QHBoxLayout()
        text2 = QtWidgets.QLabel(self)
        text2.setText('Number of peaks')
        hbox42.addWidget(text2)
        self.ngaussspin = QtWidgets.QSpinBox()
        self.ngaussspin.setMinimumWidth(85)
        self.ngaussspin.setRange(0,12)
        self.ngaussspin.setValue(0)
        self.ngaussspin.valueChanged[int].connect(self.OnNgaussspin)
        hbox42.addWidget(text2)
        hbox42.addWidget(self.ngaussspin)
        vbox4.addLayout(hbox42)


#         self.button_initfitparams = QtWidgets.QPushButton('Initialize Fit Parameters')
#         self.button_initfitparams.clicked.connect( self.OnInitFitParams)
#         self.button_initfitparams.setEnabled(False)
#         vbox4.addWidget(self.button_initfitparams)

        self.button_fitspec = QtWidgets.QPushButton('Fit Spectrum')
        self.button_fitspec.clicked.connect( self.OnFitSpectrum)
        self.button_fitspec.setEnabled(False)
        vbox4.addWidget(self.button_fitspec)

        sizer4.setLayout(vbox4)


        #panel 5

        sizer5 = QtWidgets.QGroupBox('Peak Positions [eV]')
        vbox5 = QtWidgets.QVBoxLayout()

        self.tc_peakid = []
        for i in range(12):
            tc = QtWidgets.QLabel(self)
            self.tc_peakid.append(tc)
            vbox5.addWidget(tc)

        vbox5.addStretch(1)

        self.cb_showengs = QtWidgets.QCheckBox('Show eV lines', self)
        self.cb_showengs.stateChanged.connect(self.OnShowEngs)
        vbox5.addWidget(self.cb_showengs)

        sizer5.setLayout(vbox5)


        vboxR.addLayout(vbox2)
        vboxR.addWidget(sizer5)

        vboxL.addWidget(sizer3)
        vboxL.addWidget(sizer4)
        vboxL.addWidget(sizer)

        hboxT.setContentsMargins(20,20,20,20)

        hboxT.addStretch(1)
        hboxT.addLayout(vboxL)
        hboxT.addStretch(1)
        hboxT.addLayout(vboxR)
        hboxT.addStretch(1)
        self.setLayout(hboxT)


#----------------------------------------------------------------------
    def OnSpecFromFile(self, event):
        #ToDo: Clean-up needed. One dialog for all pages is sufficient.

        try:
        #if True:
            wildcard = "ASCII files (*.csv *.txt *.xas);;CSV spectrum (*.csv);;TXT spectrum (*.txt);;XAS spectrum (*.xas)"
            directory = self.com.path
            filepath, _filter = QtWidgets.QFileDialog.getOpenFileName(self, 'Choose Spectrum file', directory, wildcard)

            filepath = str(filepath)
            if filepath == '':
                return

            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))

            self.anlz.load_xraypeakfit_spectra(filename=filepath)
            if not self.anlz.xrayfitsp_to_add:
                raise IndexError
            for i in range(self.anlz.xrayfitsp_to_add):
                self.com.xpf_loaded = 1
                self.spectrumfitted.append(0)
                self.firstinit.append(1)
                self.initdone.append(0)
                self.nsteps.append(1)
                self.npeaks.append(4)
                self.fits.append(0)
                self.fits_sep.append(0)

            curr_idx = self.anlz.n_xrayfitsp - 1
            self.i_spec = self.anlz.n_xrayfitsp

            self.nstepsspin.setValue(self.nsteps[self.i_spec-1])
            self.ngaussspin.setValue(self.npeaks[self.i_spec-1])
            self.slider_spec.setMaximum(self.i_spec)
            self.loadSpectrum()
            self.updatewidgets()
            self.slider_spec.setValue(self.i_spec)
            QtWidgets.QApplication.restoreOverrideCursor()
        except IndexError:
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, 'Error', 'Spectrum file not loaded. Data invalid or zero.')

        except TypeError:
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, 'Error', 'Spectrum file not loaded. Unknown format.')

        except Exception as error:
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, 'Error', "Unexpected Error: "+ str(error) +" "+ str(type(error)))

        self.window().refresh_widgets()


#----------------------------------------------------------------------
    def OnAddClusterSpectra(self, event):

        try:

            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))


            for i in range(self.anlz.nclusters):
                self.anlz.load_xraypeakfit_clusterspectrum(i)

                self.spectrumfitted.append(0)
                self.firstinit.append(1)
                self.initdone.append(0)
                self.nsteps.append(1)
                self.npeaks.append(4)
                self.fits.append(0)
                self.fits_sep.append(0)


            self.com.xpf_loaded = 1
            self.i_spec = self.anlz.n_xrayfitsp
            self.slider_spec.setMaximum(self.anlz.n_xrayfitsp)
            self.slider_spec.setValue(self.i_spec)
            self.nstepsspin.setValue(self.nsteps[self.i_spec-1])
            self.ngaussspin.setValue(self.npeaks[self.i_spec-1])


            self.loadSpectrum()
            #self.ShowSpectraList()

            QtWidgets.QApplication.restoreOverrideCursor()

        except:
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, 'Error', 'Cluster spectra not loaded.')



        self.window().refresh_widgets()
        self.updatewidgets()

#----------------------------------------------------------------------
    def OnSScroll(self, value):


        sel = value
        self.i_spec = sel

        #self.tc_speclist.setCurrentRow(self.i_spec-1)


        if self.anlz.n_xrayfitsp > 0:
            self.loadSpectrum()
            self.ShowFitParams()


#----------------------------------------------------------------------
    def OnSpectraListClick(self):
        item = self.tc_speclist.currentRow()

        sel = item
        self.i_spec = sel+1

        if self.com.xpf_loaded == 1:
            self.loadSpectrum()
            self.slider_spec.setValue(self.i_spec)


#----------------------------------------------------------------------
    def OnFitSpectrum(self, event):

        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))

        self.SetFitParams()

        fitted_spectrum, separate_peaks = self.anlz.fit_spectrum(self.i_spec-1, self.nsteps[self.i_spec-1], self.npeaks[self.i_spec-1])

        self.fits[self.i_spec-1] = fitted_spectrum
        self.fits_sep[self.i_spec-1] = separate_peaks

        self.spectrumfitted[self.i_spec-1] = 1

        self.loadSpectrum()

        QtWidgets.QApplication.restoreOverrideCursor()
        self.updatewidgets()



#----------------------------------------------------------------------
    def OnSave(self, event):


        #Save images
        wildcard = "Portable Network Graphics (*.png);;Adobe PDF Files (*.pdf);; SVG (*.svg)"

        fileName = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Fit Plot', '', wildcard)

        fileName = str(fileName)
        if fileName == '':
            return

        path, ext = os.path.splitext(fileName)
        ext = ext[1:].lower()


        if ext != 'png' and ext != 'pdf' and ext != 'svg':
            error_message = (
                  'Only the PNG, PDF and SVG image formats are supported.\n'
                 'A file extension of `png\' or `pdf\' must be used.')

            QtWidgets.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % error_message)
            return



        try:

            matplotlib.rcParams['pdf.fonttype'] = 42

            fig = self.Specfig
            fig.savefig(fileName)


        except IOError as e:
            if e.strerror:
                err = e.strerror
            else:
                err = e


            QtWidgets.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err)


        #Save text file with fit info
        textfilepath = path+'_'+self.anlz.xfspec_names[self.i_spec-1]+'_fitinfo.txt'
        f = open(textfilepath, 'w')
        print('*********************  Fit Results  ********************', file=f)

        print('\n', file=f)
        print('Base:\t\t'+'{0:04.3f}'.format(self.anlz.xfitpars[self.i_spec-1].base), file=f)

        print('\n', file=f)
        text = 'Step Inflection Point [eV]:\t'
        for i in range(self.nsteps[self.i_spec-1]):
            text = text + '{0:04.3f}'.format(self.anlz.xfitpars[self.i_spec-1].stepfitparams[i*3]) + ',  '
        print(text, file=f)

        text = 'Step Height:\t'
        for i in range(self.nsteps[self.i_spec-1]):
            text = text + '{0:04.3f}'.format(self.anlz.xfitpars[self.i_spec-1].stepfitparams[i*3+1]) + ',  '
        print(text, file=f)

        text = 'Step FWHM:\t'
        for i in range(self.nsteps[self.i_spec-1]):
            text = text + '{0:04.3f}'.format(self.anlz.xfitpars[self.i_spec-1].stepfitparams[i*3+2]) + ',  '
        print(text, file=f)

        print('\n', file=f)
        text = 'Peak Positions:\t'
        for i in range(self.npeaks[self.i_spec-1]):
            text = text + '{0:04.3f}'.format(self.anlz.xfitpars[self.i_spec-1].gauss_fp_m[i]) + ',  '
        print(text, file=f)

        text = 'Peak Sigma:\t'
        for i in range(self.npeaks[self.i_spec-1]):
            text = text + '{0:04.3f}'.format(self.anlz.xfitpars[self.i_spec-1].gauss_fp_s[i]) + ',  '
        print(text, file=f)

        text = 'Peak Amplitude:\t'
        for i in range(self.npeaks[self.i_spec-1]):
            text = text + '{0:04.3f}'.format(self.anlz.xfitpars[self.i_spec-1].gauss_fp_a[i]) + ',  '
        print(text, file=f)

        f.close()

        return

#----------------------------------------------------------------------
    def OnShowFitParams(self, event):

        fitparams = [self.anlz.xfitpars[self.i_spec-1].base,
                     self.anlz.xfitpars[self.i_spec-1].stepfitparams,
                     self.anlz.xfitpars[self.i_spec-1].gauss_fp_a,
                     self.anlz.xfitpars[self.i_spec-1].gauss_fp_m,
                     self.anlz.xfitpars[self.i_spec-1].gauss_fp_s]

        fitparswin = FitParams(self.window(), 'Fit Parameters', fitparams, True, False)
        fitparswin.show()

#----------------------------------------------------------------------
    def OnNstepsspin(self, value):
        num = value
        self.nsteps[self.i_spec-1] = num

#----------------------------------------------------------------------
    def OnNgaussspin(self, value):
        num = value
        self.npeaks[self.i_spec-1] = num


#----------------------------------------------------------------------
    def OnShowEngs(self, state):

        if state == QtCore.Qt.Checked:
            self.showevlines = 1
        else: self.showevlines = 0

        if self.anlz.n_xrayfitsp > 0:
            self.loadSpectrum()


#----------------------------------------------------------------------
    def updatewidgets(self):


        self.button_fitspec.setEnabled(True)

        for i in range(12):
            self.tc_peakid[i].setText('')

        peaknames = []
        if len(self.npeaks) > 0:
            for i in range(self.npeaks[self.i_spec-1]):
                i_eng=(np.abs(self.peak_engs-self.anlz.xfitpars[self.i_spec-1].gauss_fp_m[i])).argmin()
                diff =   np.abs(self.peak_engs[i_eng]-self.anlz.xfitpars[self.i_spec-1].gauss_fp_m[i])
                if np.abs(diff) < 0.5:
                    peaknames.append(self.peak_names[i_eng] + '\t({0:01.1f})'.format(diff))
                else:
                    peaknames.append('unknown')



        if self.spectrumfitted[self.i_spec-1] == 1:
            self.button_save.setEnabled(True)


            for i in range(self.npeaks[self.i_spec-1]):
                text = '\t{0:04.3f}'.format(self.anlz.xfitpars[self.i_spec-1].gauss_fp_m[i])
                self.tc_peakid[i].setText(text + '\t' + peaknames[i])
                self.tc_peakid[i].setStyleSheet('color: rgb({0:}, {1}, {2})'.format(int(255*self.plotcolors[i][0]),
                                                                                 int(255*self.plotcolors[i][1]),
                                                                                 int(255*self.plotcolors[i][2])))
        else:
            self.button_save.setEnabled(False)


        self.ShowFitParams()


#----------------------------------------------------------------------
    def loadSpectrum(self):

        spectrum = self.anlz.xrayfitspectra[self.i_spec-1, :]


        fig = self.Specfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()


        line1 = axes.plot(self.stk.ev,spectrum, color='black')


        if self.spectrumfitted[self.i_spec-1] == 1:

            offset = self.fits_sep[self.i_spec-1][0]
            line3 = axes.plot(self.stk.ev, offset, color = 'blue')
            for i in range(1, 1+self.nsteps[self.i_spec-1]):
                y = self.fits_sep[self.i_spec-1][i]+offset
                line3 = axes.plot(self.stk.ev, y, color = 'green')
            for i in range(0, self.npeaks[self.i_spec-1]):
                y = self.fits_sep[self.i_spec-1][i+self.nsteps[self.i_spec-1]+1]+offset
                line3 = axes.plot(self.stk.ev, y, color = self.plotcolors[i])

            line2 = axes.plot(self.stk.ev,self.fits[self.i_spec-1], color='red')

        lines = axes.get_lines()
        self.colors = []
        for i, line in enumerate(lines):
            self.colors.append(line.get_color())


        if self.showevlines:
            peak_engs = self.window().tab_peak.peak_engs
            for i in range(len(peak_engs)):
                if (peak_engs[i]>self.stk.ev[0]) and (peak_engs[i]<self.stk.ev[-1]):
                    axes.axvline(x=peak_engs[i], color = 'g', alpha=0.5)


        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Optical Density')

        self.SpectrumPanel.draw()


        self.tc_spec.setText("X-ray Spectrum: " +
                               self.anlz.xfspec_names[self.i_spec-1])


        self.nstepsspin.setValue(self.nsteps[self.i_spec-1])
        self.ngaussspin.setValue(self.npeaks[self.i_spec-1])




        self.window().refresh_widgets()
        self.updatewidgets()


#----------------------------------------------------------------------
    def ShowSpectraList(self):

        self.tc_speclist.clear()

        for i in range(self.anlz.n_xrayfitsp):
            self.tc_speclist.addItem(self.anlz.xfspec_names[i])




#----------------------------------------------------------------------
    def SetFitParams(self):

        self.base = float(self.le_base.text())


        self.stepfitparams[0] = float(self.le_sa1.text())
        self.stepfitparams[1] = float(self.le_sp1.text())
        self.stepfitparams[2] = float(self.le_sf1.text())
        self.stepfitparams[3] = float(self.le_sa2.text())
        self.stepfitparams[4] = float(self.le_sp2.text())
        self.stepfitparams[5] = float(self.le_sf2.text())

        for i in range(12):

            self.gauss_fp_a[i] = float(self.le_amplitudes[i].text())
            self.gauss_fp_m[i] = float(self.le_mu[i].text())
            self.gauss_fp_s[i] = float(self.le_sigma[i].text())


        self.anlz.set_init_fit_params(self.i_spec-1, self.base, self.stepfitparams,
                                       self.gauss_fp_a, self.gauss_fp_m, self.gauss_fp_s)



#----------------------------------------------------------------------
    def ShowFitParams(self):

        fp = self.anlz.xfitpars[self.i_spec-1]
        self.le_base.setText('{0:.2f}'.format(fp.base))

        self.le_sa1.setText('{0:.2f}'.format(fp.stepfitparams[0]))
        self.le_sp1.setText('{0:.2f}'.format(fp.stepfitparams[1]))
        self.le_sf1.setText('{0:.2f}'.format(fp.stepfitparams[2]))
        self.le_sa2.setText('{0:.2f}'.format(fp.stepfitparams[3]))
        self.le_sp2.setText('{0:.2f}'.format(fp.stepfitparams[4]))
        self.le_sf2.setText('{0:.2f}'.format(fp.stepfitparams[5]))

        for i in range(12):
            self.le_amplitudes[i].setText('{0:.2f}'.format(fp.gauss_fp_a[i]))
            self.le_mu[i].setText('{0:.2f}'.format(fp.gauss_fp_m[i]))
            self.le_sigma[i].setText('{0:.2f}'.format(fp.gauss_fp_s[i]))

