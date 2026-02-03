
import os
import numpy as np
from PIL import Image
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtCore import Qt
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable

from ..dialogs.save import SaveWinP5
from ...core.constants import PlotH, PlotW

class PageNNMA(QtWidgets.QWidget):
    def __init__(self, common, data_struct, stack, anlz, nnma):
        super(PageNNMA, self).__init__()

        self.initUI(common, data_struct, stack, anlz, nnma)

#----------------------------------------------------------------------
    def initUI(self, common, data_struct, stack, anlz, nnma):


        self.com = common
        self.data_struct = data_struct
        self.stk = stack
        self.anlz = anlz
        self.nnma = nnma

        self.i_map = 0
        self.show_scale_bar = 0
        self.nnmacalculated = 0

        self.nComponents = 5
        self.maxIters = 1000
        self.deltaErrorThreshold = 1e-3
        self.initMatrices = 'Random'

        self.lambdaSparse = 0.5
        self.lambdaClusterSim = 10.
        self.lambdaSmooth = 0.

        #panel 1
        sizer1 = QtWidgets.QGroupBox('NNMA analysis')
        hbox1 = QtWidgets.QHBoxLayout()
        vbox1 = QtWidgets.QVBoxLayout()

        vbox1.addStretch(1)
        self.button_calcnnma = QtWidgets.QPushButton('Calculate NNMA')
        self.button_calcnnma.clicked.connect( self.OnCalcNNMA)
        self.button_calcnnma.setEnabled(False)
        vbox1.addWidget(self.button_calcnnma)
        self.button_muroi = QtWidgets.QPushButton('Load initial ROI spectra')
        self.button_muroi.clicked.connect( self.OnLoadROISpectra)
        self.button_muroi.setEnabled(False)
        self.button_mucluster = QtWidgets.QPushButton('Load initial cluster spectra')
        self.button_mucluster.clicked.connect( self.OnLoadClusterSpectra)
        self.button_mucluster.setEnabled(False)
        self.button_mufile = QtWidgets.QPushButton('Load initial standard spectra')
        self.button_mufile.clicked.connect( self.OnLoadStandardSpectra)
        self.button_mufile.setEnabled(False)
        self.button_murand = QtWidgets.QPushButton('Load initial random spectra')
        self.button_murand.clicked.connect( self.OnLoadRandomSpectra)
        self.button_murand.setEnabled(False)
        self.tc_initspectra = QtWidgets.QLabel(self)
        self.tc_initspectra.setText('Initial Spectra: ' + self.initMatrices)
        self.button_savennma = QtWidgets.QPushButton('Save NNMA Results...')
        self.button_savennma.clicked.connect( self.OnSave)
        self.button_savennma.setEnabled(False)
        vbox1.addWidget(self.button_savennma)

        vbox1.addStretch(1)


        sizer11 = QtWidgets.QGroupBox('NNMA settings')
        vbox11 = QtWidgets.QVBoxLayout()

        hbox11 = QtWidgets.QHBoxLayout()
        text1 = QtWidgets.QLabel(self)
        text1.setText('Number of components k')
        hbox11.addWidget(text1)
        hbox11.addStretch(1)
        self.ncompspin = QtWidgets.QSpinBox()
        self.ncompspin.setRange(2,20)
        self.ncompspin.setValue(self.nComponents)
        self.ncompspin.valueChanged[int].connect(self.OnNNMAspin)
        hbox11.addWidget(self.ncompspin)

        vbox11.addLayout(hbox11)

        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.HLine)
        line.setFrameShadow(QtWidgets.QFrame.Sunken)
        vbox11.addStretch(1)
        vbox11.addWidget(line)
        vbox11.addStretch(1)

        text1 = QtWidgets.QLabel(self)
        text1.setText('Regularization parameters')
        vbox11.addWidget(text1)

        hbox15a = QtWidgets.QHBoxLayout()
        tc1 = QtWidgets.QLabel(self)
        tc1.setText("\tSpectra similarity")
        hbox15a.addWidget(tc1)
        self.ntc_lamsim = QtWidgets.QLineEdit(self)
        self.ntc_lamsim.setFixedWidth(65)
        self.ntc_lamsim.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.ntc_lamsim.setAlignment(QtCore.Qt.AlignRight)
        hbox15a.addStretch(1)
        self.ntc_lamsim.setText(str(self.lambdaClusterSim))
        hbox15a.addWidget(self.ntc_lamsim)

        vbox11.addLayout(hbox15a)

        hbox15b = QtWidgets.QHBoxLayout()
        tc1 = QtWidgets.QLabel(self)
        tc1.setText("\tSpectra smoothness")
        hbox15b.addWidget(tc1)
        self.ntc_lamsmooth = QtWidgets.QLineEdit(self)
        self.ntc_lamsmooth.setFixedWidth(65)
        self.ntc_lamsmooth.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.ntc_lamsmooth.setAlignment(QtCore.Qt.AlignRight)
        hbox15b.addStretch(1)
        self.ntc_lamsmooth.setText(str(self.lambdaSmooth))
        hbox15b.addWidget(self.ntc_lamsmooth)

        vbox11.addLayout(hbox15b)

        hbox15c = QtWidgets.QHBoxLayout()
        tc1 = QtWidgets.QLabel(self)
        tc1.setText("\tSparseness")
        hbox15c.addWidget(tc1)
        self.ntc_lamspar = QtWidgets.QLineEdit(self)
        self.ntc_lamspar.setFixedWidth(65)
        self.ntc_lamspar.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.ntc_lamspar.setAlignment(QtCore.Qt.AlignRight)
        hbox15c.addStretch(1)
        self.ntc_lamspar.setText(str(self.lambdaSparse))
        hbox15c.addWidget(self.ntc_lamspar)

        vbox11.addLayout(hbox15c)

        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.HLine)
        line.setFrameShadow(QtWidgets.QFrame.Sunken)
        vbox11.addStretch(1)
        vbox11.addWidget(line)
        vbox11.addStretch(1)

        hbox12 = QtWidgets.QHBoxLayout()
        text1 = QtWidgets.QLabel(self)
        text1.setText('Number of iterations')
        hbox12.addWidget(text1)
        hbox12.addStretch(1)
        self.ntc_niterations = QtWidgets.QLineEdit(self)
        self.ntc_niterations.setFixedWidth(65)
        self.ntc_niterations.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.ntc_niterations.setAlignment(QtCore.Qt.AlignRight)
        self.ntc_niterations.setText(str(self.maxIters))
        hbox12.addWidget(self.ntc_niterations)

        vbox11.addLayout(hbox12)


        hbox14a = QtWidgets.QHBoxLayout()
        tc1 = QtWidgets.QLabel(self)
        tc1.setText("Delta Error Threshold")
        hbox14a.addWidget(tc1)
        self.ntc_dthresh = QtWidgets.QLineEdit(self)
        self.ntc_dthresh.setFixedWidth(65)
        self.ntc_dthresh.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.ntc_dthresh.setAlignment(QtCore.Qt.AlignRight)
        hbox14a.addStretch(1)
        self.ntc_dthresh.setText(str(self.deltaErrorThreshold))
        hbox14a.addWidget(self.ntc_dthresh)

        vbox11.addLayout(hbox14a)
        sizer11.setLayout(vbox11)

        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.HLine)
        line.setFrameShadow(QtWidgets.QFrame.Sunken)


        vbox1.addStretch(1)
        vbox1.addWidget(line)
        vbox1.addStretch(1)
        vbox1.addWidget(self.button_mucluster)
        vbox1.addWidget(self.button_muroi)
        vbox1.addWidget(self.button_mufile)
        vbox1.addWidget(self.button_murand)
        vbox1.addWidget(self.tc_initspectra)


        hbox1.addLayout(vbox1)
        hbox1.addWidget(sizer11)

        sizer1.setLayout(hbox1)


        #panel 2
        vbox2 = QtWidgets.QVBoxLayout()

        self.tc_NNMAcomp = QtWidgets.QLabel(self)
        self.tc_NNMAcomp.setText("NNMA thickness map ")
        vbox2.addWidget(self.tc_NNMAcomp)

        hbox21 = QtWidgets.QHBoxLayout()

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()
        self.nnmaimgfig = Figure((PlotH*1.10, PlotH))
        self.NNMAImagePan = FigureCanvas(self.nnmaimgfig)
        self.NNMAImagePan.setParent(self)

        fbox.addWidget(self.NNMAImagePan)
        frame.setLayout(fbox)
        hbox21.addWidget(frame)

        self.slidershow = QtWidgets.QScrollBar(QtCore.Qt.Vertical)
        self.slidershow.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slidershow.setRange(0,20)
        #self.slidershow.setEnabled(False)
        self.slidershow.valueChanged[int].connect(self.OnScroll)


        hbox21.addWidget(self.slidershow)

        vbox2.addLayout(hbox21)


        #panel 3
        vbox3 = QtWidgets.QVBoxLayout()

        self.tc_nnmaspec = QtWidgets.QLabel(self)
        self.tc_nnmaspec.setText("NNMA spectrum ")
        vbox3.addWidget(self.tc_nnmaspec)

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.nnmaspecfig = Figure((PlotW, PlotH))
        self.NNMASpecPan = FigureCanvas(self.nnmaspecfig)
        self.NNMASpecPan.setParent(self)

        fbox.addWidget(self.NNMASpecPan)
        frame.setLayout(fbox)
        vbox3.addWidget(frame)


        #panel 4
        vbox4 = QtWidgets.QVBoxLayout()

        text4 = QtWidgets.QLabel(self)
        text4.setText("Cost function")
        vbox4.addWidget(text4)

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.costffig = Figure((PlotW*0.7, PlotH*0.7))
        self.CostFPan = FigureCanvas(self.costffig)
        self.CostFPan.setParent(self)
        #self.CostFPan.mpl_connect('button_press_event', self.OnPointEvalsImage)

        fbox.addWidget(self.CostFPan)
        frame.setLayout(fbox)
        vbox4.addWidget(frame)


        #panel 5
        sizer5 = QtWidgets.QGroupBox('Display')
        vbox5 = QtWidgets.QVBoxLayout()

        hbox51 = QtWidgets.QHBoxLayout()
        self.cb_inputsp = QtWidgets.QCheckBox('Overlay input spectra', self)
        self.cb_inputsp.stateChanged.connect(self.ShowSpectrum)
        hbox51.addWidget(self.cb_inputsp)
        vbox5.addLayout(hbox51)

        hbox52 = QtWidgets.QHBoxLayout()
        self.cb_showallsp = QtWidgets.QCheckBox('Show all spectra', self)
        self.cb_showallsp.stateChanged.connect(self.ShowSpectrum)
        hbox52.addWidget(self.cb_showallsp)
        vbox5.addLayout(hbox52)

        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.HLine)
        line.setFrameShadow(QtWidgets.QFrame.Sunken)

        vbox5.addWidget(line)

        hbox52 = QtWidgets.QHBoxLayout()
        self.cb_totalcost = QtWidgets.QCheckBox('Total cost', self)
        self.cb_totalcost.setChecked(True)
        self.cb_totalcost.stateChanged.connect(self.ShowCostFunction)
        hbox52.addWidget(self.cb_totalcost)
        vbox5.addLayout(hbox52)

        hbox53 = QtWidgets.QHBoxLayout()
        self.cb_simlcost = QtWidgets.QCheckBox('Similarity cost', self)
        self.cb_simlcost.stateChanged.connect(self.ShowCostFunction)
        hbox53.addWidget(self.cb_simlcost)
        vbox5.addLayout(hbox53)

        hbox54 = QtWidgets.QHBoxLayout()
        self.cb_smoothcost = QtWidgets.QCheckBox('Smoothness cost', self)
        self.cb_smoothcost.stateChanged.connect(self.ShowCostFunction)
        hbox54.addWidget(self.cb_smoothcost)
        vbox5.addLayout(hbox54)

        hbox54 = QtWidgets.QHBoxLayout()
        self.cb_sparcost = QtWidgets.QCheckBox('Sparsness cost', self)
        self.cb_sparcost.stateChanged.connect(self.ShowCostFunction)
        hbox54.addWidget(self.cb_sparcost)
        vbox5.addLayout(hbox54)

        sizer5.setLayout(vbox5)



        vboxtop = QtWidgets.QVBoxLayout()
        hboxtopL = QtWidgets.QHBoxLayout()
        vboxtopL = QtWidgets.QVBoxLayout()
        vboxtopL.addWidget(sizer1)
        hboxtopL.addLayout(vboxtopL)
        hboxtopL.addWidget(sizer5)
        hboxtopL.addStretch(1)

        gridsizertop = QtWidgets.QGridLayout()
        gridsizertop.setContentsMargins(15,0,0,0)


        gridsizertop.addLayout(hboxtopL, 0, 0, QtCore .Qt. AlignLeft)
        gridsizertop.addLayout(vbox4, 0, 1, QtCore .Qt. AlignLeft)

        gridsizertop.addLayout(vbox3, 1, 0, QtCore .Qt. AlignCenter)
        gridsizertop.addLayout(vbox2, 1, 1, QtCore .Qt. AlignLeft)

        vboxtop.addStretch(1)
        vboxtop.addLayout(gridsizertop)
        vboxtop.addStretch(1)
        self.setLayout(vboxtop)

    # ----------------------------------------------------------------------
    def OnLoadROISpectra(self, event):
        #self.ncompspin.setValue(xxxx n rois)
        self.ncompspin.setEnabled(False)
        plot_num = len(self.window().page1.specfig.pi.items)
        kNNMA = 0
        for i in range(plot_num - 1):
            spec = self.window().page1.specfig.pi.items[i + 1].yData
            if not np.any(spec):  # exclude spectra with only zeros
                continue
            kNNMA += 1
            if kNNMA == 1:
                muInit = spec
            elif kNNMA > 1:
                muInit = np.vstack((muInit, spec))
        if kNNMA < 2:
            QtWidgets.QMessageBox.warning(self, 'Error', 'Select at least two ROIs.')
            return
        self.ncompspin.setValue(kNNMA)
        self.ncompspin.setEnabled(False)
        self.nnma.setROISpectra(muInit.transpose())

        self.window().refresh_widgets()

        self.initMatrices = 'ROI'

        self.tc_initspectra.setText('Initial Spectra: ' + self.initMatrices)

    #----------------------------------------------------------------------
    def OnLoadClusterSpectra(self, event):
        self.ncompspin.setValue(self.window().page3.numclusters)
        self.ncompspin.setEnabled(False)
        self.nnma.setClusterSpectra(self.anlz.clusterspectra)

        self.window().refresh_widgets()

        self.initMatrices = 'Cluster'

        self.tc_initspectra.setText('Initial Spectra: ' + self.initMatrices)

#----------------------------------------------------------------------
    def OnLoadRandomSpectra(self, event):
        self.ncompspin.setEnabled(True)
        self.initMatrices = 'Random'

        self.lambdaClusterSim = 0.0
        self.ntc_lamsim.setText(str(self.lambdaClusterSim))

        self.tc_initspectra.setText('Initial Spectra: ' + self.initMatrices)

#----------------------------------------------------------------------
    def OnLoadStandardSpectra(self, event):

        #if True:
        try:
            wildcard = "ASCII files (*.csv *.txt *.xas);;CSV spectrum (*.csv);;TXT spectrum (*.txt);;XAS spectrum (*.xas)"
            directory = self.com.path
            filepath, _filter = QtWidgets.QFileDialog.getOpenFileNames(self, 'Choose Spectrum file(s)', directory, wildcard)

            #filepath, _filter = QtWidgets.QFileDialog.getOpenFileNames(self, 'Choose files containing Spectra', '', wildcard,
            #                                                           None, QtWidgets.QFileDialog.DontUseNativeDialog)
            if filepath == '':
                return

            self.filenames = []
            for name in filepath:
                self.filenames.append(str(name))
            muInit = []
            kNNMA = 0

            for f in self.filenames:
                spectrum_common_names, fitspectra = self.anlz.load_filter_remap_spectra(f)

                try:
                    muInit = np.vstack((muInit, fitspectra))
                except ValueError:
                    muInit  = fitspectra
                kNNMA += len(spectrum_common_names)

            if kNNMA < 2:
                QtWidgets.QMessageBox.warning(self, 'Error', 'Select at least two spectra.')
                return
            self.ncompspin.setValue(kNNMA)
            self.ncompspin.setEnabled(False)
            self.nnma.setStandardsSpectra(muInit.transpose())
            self.initMatrices = 'Standards'
            self.lambdaClusterSim = 0.0
            self.ntc_lamsim.setText(str(self.lambdaClusterSim))

            self.tc_initspectra.setText('Initial Spectra: ' + self.initMatrices)
            QtWidgets.QApplication.restoreOverrideCursor()

        except:
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, 'Error', 'Spectra not loaded.')

        self.window().refresh_widgets()

#----------------------------------------------------------------------
    def OnCalcNNMA(self, event):

        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))

        value = self.ntc_niterations.text()
        self.maxIters = int(value)
        value = self.ntc_dthresh.text()
        self.deltaErrorThreshold = float(value)

        value = self.ntc_lamspar.text()
        self.lambdaSparse = float(value)
        value = self.ntc_lamsim.text()
        self.lambdaClusterSim = float(value)
        value = self.ntc_lamsmooth.text()
        self.lambdaSmooth = float(value)




        self.nnma.setParameters(kNNMA = self.nComponents,
                                maxIters = self.maxIters,
                                deltaErrorThreshold = self.deltaErrorThreshold,
                                initMatrices = self.initMatrices,
                                lambdaSparse = self.lambdaSparse,
                                lambdaClusterSim = self.lambdaClusterSim,
                                lambdaSmooth = self.lambdaSmooth)

        self.nnma.calcNNMA(initmatrices = self.initMatrices)

        self.nnmacalculated = 1
        self.i_map = 0
        self.slidershow.setValue(self.i_map)
        self.slidershow.setMaximum(self.nComponents-1)
        self.ShowMaps()
        self.ShowSpectrum()
        self.ShowCostFunction()
        self.updatewidgets()

        QtWidgets.QApplication.restoreOverrideCursor()



#----------------------------------------------------------------------
    def OnNNMAspin(self, value):
        num = value
        self.nComponents = num

#----------------------------------------------------------------------
    def OnScroll(self, value):
        self.i_map = value
        if self.nnmacalculated == 1:
            self.ShowMaps()
            self.ShowSpectrum()

#----------------------------------------------------------------------
    def OnSave(self, event):

        savewin = SaveWinP5(self.window())
        savewin.show()


#----------------------------------------------------------------------
    def Save(self, filename, path, spec_png = True, spec_pdf = False, spec_svg = False, spec_csv = False,
             map_png = True, map_pdf = False, map_svg = False, map_tif = False,
             costf_png = True, costf_pdf = False, costf_svg = False):


        self.SaveFileName = os.path.join(path,filename)

        try:
        #if True:
            matplotlib.rcParams['pdf.fonttype'] = 42
            if costf_png:
                ext = 'png'
                fileName_evals = self.SaveFileName+"_CostFunction."+ext

                fig = self.costffig
                fig.savefig(fileName_evals)

            if costf_pdf:
                ext = 'pdf'
                fileName_evals = self.SaveFileName+"_CostFunction."+ext

                fig = self.costffig
                fig.savefig(fileName_evals)

            if costf_svg:
                ext = 'svg'
                fileName_evals = self.SaveFileName+"_CostFunction."+ext

                fig = self.costffig
                fig.savefig(fileName_evals)


            from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
            matplotlib.rcParams['pdf.fonttype'] = 42

            ext = 'png'

            colors=['#FF0000','#000000','#FFFFFF']
            spanclrmap=matplotlib.colors.LinearSegmentedColormap.from_list('spancm',colors)

            if map_png:
                for i in range(self.nComponents):

                    mapimage = self.nnma.tRecon[i, :,:]

                    fig = matplotlib.figure.Figure(figsize =(PlotH*1.15, PlotH))
                    _canvas = FigureCanvas(fig)
                    axes = fig.gca()

                    divider = make_axes_locatable(axes)
                    ax_cb = divider.new_horizontal(size="3%", pad=0.03)
                    fig.add_axes(ax_cb)
                    axes.set_position([0.03,0.03,0.8,0.94])

                    min_val = np.min(mapimage)
                    max_val = np.max(mapimage)
                    bound = np.max((np.abs(min_val), np.abs(max_val)))

                    if self.show_scale_bar == 1:
                        um_string = ' $\mathrm{\mu m}$'
                        microns = '$'+self.stk.scale_bar_string+' $'+um_string
                        axes.text(self.stk.scale_bar_pixels_x+10,self.stk.n_cols-9, microns, horizontalalignment='left', verticalalignment='center',
                                  color = 'white', fontsize=14)
                        #Matplotlib has flipped scales so I'm using rows instead of cols!
                        p = matplotlib.patches.Rectangle((5,self.stk.n_cols-10), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                                               color = 'white', fill = True)
                        axes.add_patch(p)

                    im = axes.imshow(np.rot90(mapimage), cmap=spanclrmap, vmin = -bound, vmax = bound)
                    _cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)

                    axes.axis("off")

                    fileName_img = self.SaveFileName+"_NNMA_" +str(i+1)+"."+ext
                    fig.savefig(fileName_img, pad_inches = 0.0)

            if spec_png:
                for i in range(self.nComponents):

                    fig = matplotlib.figure.Figure(figsize =(PlotW, PlotH))
                    _canvas = FigureCanvas(fig)
                    fig.add_axes((0.15,0.15,0.75,0.75))
                    axes = fig.gca()

                    nspectrum = self.nnma.muRecon[:,i]

                    line1 = axes.plot(self.stk.ev,nspectrum)
                    lcolor = line1[0].get_color()

                    if self.cb_inputsp.isChecked():
                        initspec = self.nnma.muinit[:,self.i_map]
                        _line2 = axes.plot(self.stk.ev,initspec, color=lcolor, linestyle = '--', label = 'Fit')

                    axes.set_xlabel('Photon Energy [eV]')
                    axes.set_ylabel('Absorption coefficient [a.u.]')

                    fileName_spec = self.SaveFileName+"_NNMAspectrum_" +str(i+1)+"."+ext
                    fig.savefig(fileName_spec)

            if spec_csv:
                for i in range(self.nComponents):
                    nspectrum = self.nnma.muRecon[:,i]
                    fileName_spec = self.SaveFileName+"_NNMAspectrum_" +str(i+1)+".csv"
                    cname = "NNMAspectrum_" +str(i+1)
                    self.stk.write_csv(fileName_spec, self.stk.ev, nspectrum, cname = cname)


            ext = 'pdf'

            if map_pdf:
                for i in range(self.nComponents):

                    mapimage = self.nnma.tRecon[i, :,:]

                    fig = matplotlib.figure.Figure(figsize =(PlotH*1.15, PlotH))
                    _canvas = FigureCanvas(fig)
                    axes = fig.gca()

                    divider = make_axes_locatable(axes)
                    ax_cb = divider.new_horizontal(size="3%", pad=0.03)
                    fig.add_axes(ax_cb)
                    axes.set_position([0.03,0.03,0.8,0.94])

                    min_val = np.min(mapimage)
                    max_val = np.max(mapimage)
                    bound = np.max((np.abs(min_val), np.abs(max_val)))

                    if self.show_scale_bar == 1:
                        um_string = ' $\mathrm{\mu m}$'
                        microns = '$'+self.stk.scale_bar_string+' $'+um_string
                        axes.text(self.stk.scale_bar_pixels_x+10,self.stk.n_cols-9, microns, horizontalalignment='left', verticalalignment='center',
                                  color = 'white', fontsize=14)
                        #Matplotlib has flipped scales so I'm using rows instead of cols!
                        p = matplotlib.patches.Rectangle((5,self.stk.n_cols-10), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                                               color = 'white', fill = True)
                        axes.add_patch(p)

                    im = axes.imshow(np.rot90(mapimage), cmap=spanclrmap, vmin = -bound, vmax = bound)
                    _cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)

                    axes.axis("off")

                    fileName_img = self.SaveFileName+"_NNMA_" +str(i+1)+"."+ext
                    fig.savefig(fileName_img, pad_inches = 0.0)

            if spec_pdf:
                for i in range(self.nComponents):


                    fig = matplotlib.figure.Figure(figsize =(PlotW, PlotH))
                    _canvas = FigureCanvas(fig)
                    fig.add_axes((0.15,0.15,0.75,0.75))
                    axes = fig.gca()

                    nspectrum = self.nnma.muRecon[:,i]

                    line1 = axes.plot(self.stk.ev,nspectrum)
                    lcolor = line1[0].get_color()

                    if self.cb_inputsp.isChecked():
                        initspec = self.nnma.muinit[:,self.i_map]
                        _line2 = axes.plot(self.stk.ev,initspec, color=lcolor, linestyle = '--', label = 'Fit')

                    axes.set_xlabel('Photon Energy [eV]')
                    axes.set_ylabel('Absorption coefficient [a.u.]')

                    fileName_spec = self.SaveFileName+"_NNMAspectrum_" +str(i+1)+"."+ext
                    fig.savefig(fileName_spec)

            ext = 'svg'

            if map_svg:
                for i in range(self.nComponents):

                    mapimage = self.nnma.tRecon[i, :,:]

                    fig = matplotlib.figure.Figure(figsize =(PlotH*1.15, PlotH))
                    _canvas = FigureCanvas(fig)
                    axes = fig.gca()

                    divider = make_axes_locatable(axes)
                    ax_cb = divider.new_horizontal(size="3%", pad=0.03)
                    fig.add_axes(ax_cb)
                    axes.set_position([0.03,0.03,0.8,0.94])

                    min_val = np.min(mapimage)
                    max_val = np.max(mapimage)
                    bound = np.max((np.abs(min_val), np.abs(max_val)))

                    if self.show_scale_bar == 1:
                        um_string = ' $\mathrm{\mu m}$'
                        microns = '$'+self.stk.scale_bar_string+' $'+um_string
                        axes.text(self.stk.scale_bar_pixels_x+10,self.stk.n_cols-9, microns, horizontalalignment='left', verticalalignment='center',
                                  color = 'white', fontsize=14)
                        #Matplotlib has flipped scales so I'm using rows instead of cols!
                        p = matplotlib.patches.Rectangle((5,self.stk.n_cols-10), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                                               color = 'white', fill = True)
                        axes.add_patch(p)

                    im = axes.imshow(np.rot90(mapimage), cmap=spanclrmap, vmin = -bound, vmax = bound)
                    _cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)

                    axes.axis("off")

                    fileName_img = self.SaveFileName+"_NNMA_" +str(i+1)+"."+ext
                    fig.savefig(fileName_img, pad_inches = 0.0)

            if spec_svg:
                for i in range(self.nComponents):


                    fig = matplotlib.figure.Figure(figsize =(PlotW, PlotH))
                    _canvas = FigureCanvas(fig)
                    fig.add_axes((0.15,0.15,0.75,0.75))
                    axes = fig.gca()

                    nspectrum = self.nnma.muRecon[:,i]

                    line1 = axes.plot(self.stk.ev,nspectrum)
                    lcolor = line1[0].get_color()

                    if self.cb_inputsp.isChecked():
                        initspec = self.nnma.muinit[:,self.i_map]
                        _line2 = axes.plot(self.stk.ev,initspec, color=lcolor, linestyle = '--', label = 'Fit')

                    axes.set_xlabel('Photon Energy [eV]')
                    axes.set_ylabel('Absorption coefficient [a.u.]')

                    fileName_spec = self.SaveFileName+"_NNMAspectrum_" +str(i+1)+"."+ext
                    fig.savefig(fileName_spec)


            if map_tif:
                for i in range(self.nComponents):
                    mapimage = self.nnma.tRecon[i, :,:]
                    fileName_img = self.SaveFileName+"_NNMA_" +str(i+1)+".tif"
                    img1 = Image.fromarray(mapimage)
                    img1.save(fileName_img)


        except IOError as e:
            if e.strerror:
                err = e.strerror
            else:
                err = e

            QtWidgets.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err)
#


#----------------------------------------------------------------------
    def ShowMaps(self):


        mapimage = self.nnma.tRecon[self.i_map, :,:]


        colors=['#FF0000','#000000','#FFFFFF']

        spanclrmap=matplotlib.colors.LinearSegmentedColormap.from_list('spancm',colors)

        fig = self.nnmaimgfig
        fig.clf()

        axes = fig.gca()

        divider = make_axes_locatable(axes)
        ax_cb = divider.new_horizontal(size="3%", pad=0.03)

        fig.add_axes(ax_cb)

        axes.set_position([0.03,0.03,0.8,0.94])


        min_val = np.min(mapimage)
        max_val = np.max(mapimage)
        bound = np.max((np.abs(min_val), np.abs(max_val)))

        if self.show_scale_bar == 1:
            um_string = ' $\mathrm{\mu m}$'
            microns = '$'+self.stk.scale_bar_string+' $'+um_string
            axes.text(self.stk.scale_bar_pixels_x+10,self.stk.n_cols-9, microns, horizontalalignment='left', verticalalignment='center',
                      color = 'white', fontsize=14)
            #Matplotlib has flipped scales so I'm using rows instead of cols!
            p = matplotlib.patches.Rectangle((5,self.stk.n_cols-10), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                                   color = 'white', fill = True)
            axes.add_patch(p)

        im = axes.imshow(np.rot90(mapimage), cmap=spanclrmap, vmin = -bound, vmax = bound)
        cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)

        axes.axis("off")
        self.NNMAImagePan.draw()



#----------------------------------------------------------------------
    def ShowSpectrum(self):

        if self.nnmacalculated == 0:
            return

        fig = self.nnmaspecfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()


        if self.cb_showallsp.isChecked():
            for i in range(self.nComponents):
                nspectrum = self.nnma.muRecon[:,i]

                line1 = axes.plot(self.stk.ev,nspectrum)
                lcolor = line1[0].get_color()

                if self.cb_inputsp.isChecked():
                    initspec = self.nnma.muinit[:,i]
                    _line2 = axes.plot(self.stk.ev,initspec, color=lcolor, linestyle = '--', label = 'Fit')

        else:

            nspectrum = self.nnma.muRecon[:,self.i_map]

            line1 = axes.plot(self.stk.ev,nspectrum)
            lcolor = line1[0].get_color()

            if self.cb_inputsp.isChecked():
                initspec = self.nnma.muinit[:,self.i_map]
                line2 = axes.plot(self.stk.ev,initspec, color=lcolor, linestyle = '--', label = 'Fit')



        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Absorption coefficient [a.u.]')

        self.NNMASpecPan.draw()

        self.tc_nnmaspec.setText("NNMA spectrum " + str(self.i_map+1))


#----------------------------------------------------------------------
    def ShowCostFunction(self):

        if self.nnmacalculated == 0:
            return

        costTotal = self.nnma.costFnArray[:,0]
        costSparse = self.nnma.costFnArray[:,1]
        costClusterSim = self.nnma.costFnArray[:,2]
        costSmooth = self.nnma.costFnArray[:,3]


        fig = self.costffig
        fig.clf()
        fig.add_axes((0.20,0.15,0.75,0.75))
        axes = fig.gca()

        if self.cb_totalcost.isChecked():
            _line1 = axes.plot(costTotal, label='Total')

        if self.cb_sparcost.isChecked():
            _line2 = axes.plot(costSparse, label='Sparseness ')

        if self.cb_simlcost.isChecked():
            _line3 = axes.plot(costClusterSim, label='Similarity')

        if self.cb_smoothcost.isChecked():
            _line4 = axes.plot(costSmooth, label='Smoothness')


        fontP = matplotlib.font_manager.FontProperties()
        fontP.set_size('small')

        axes.legend(loc=1, prop = fontP)


        axes.set_xlabel('Iteration')
        axes.set_ylabel('Cost Function')

        self.CostFPan.draw()


#----------------------------------------------------------------------
    def updatewidgets(self):

        if self.nnmacalculated == 0:
            self.button_savennma.setEnabled(False)
        else:
            self.button_savennma.setEnabled(True)
