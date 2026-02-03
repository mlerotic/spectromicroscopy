
import os
import numpy as np
from PIL import Image
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtCore import Qt
import matplotlib
import matplotlib.cm
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable

from ..dialogs.save import SaveWinP4
from ..dialogs.composite_rgb import ShowCompositeRBGmap 
from ..dialogs.histogram import ShowMapHistogram

from ...core.constants import PlotH, PlotW

class PageSpectral(QtWidgets.QWidget):
    def __init__(self, common, data_struct, stack, anlz):
        super(PageSpectral, self).__init__()

        self.initUI(common, data_struct, stack, anlz)

#----------------------------------------------------------------------
    def initUI(self, common, data_struct, stack, anlz):

        self.data_struct = data_struct
        self.stk = stack
        self.com = common
        self.anlz = anlz

        self.i_tspec = 1
        self.showraw = True
        self.show_scale_bar = 0
        self.itheta = 0


        vbox = QtWidgets.QVBoxLayout()
        hboxT = QtWidgets.QHBoxLayout()
        hboxB = QtWidgets.QHBoxLayout()


        #panel 5
        sizer5 = QtWidgets.QGroupBox('Target Spectra')
        vbox5 = QtWidgets.QVBoxLayout()

        self.tc_speclist = QtWidgets.QListWidget()
        self.tc_speclist.itemClicked.connect(self.OnSpectraListClick)
        #self.tc_speclist.doubleClicked.connect(self.OnEditSpectraListClick)
        vbox5.addWidget(self.tc_speclist)
        sizer5.setLayout(vbox5)



        #panel 1
        vbox1 = QtWidgets.QVBoxLayout()

        self.tc_spmap = QtWidgets.QLabel(self)
        self.tc_spmap.setText("Spectrum composition map")

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.mapfig = Figure((PlotH, PlotH))
        self.MapPanel = FigureCanvas(self.mapfig)
        self.MapPanel.setParent(self)
        fbox.addWidget(self.MapPanel)
        frame.setLayout(fbox)

        vbox1.addWidget(self.tc_spmap)
        vbox1.addWidget(frame)




        #panel 2
        vbox2 = QtWidgets.QVBoxLayout()

        self.tc_tspec = QtWidgets.QLabel(self)
        self.tc_tspec.setText("Target Spectrum: ")
        hbox11 = QtWidgets.QHBoxLayout()

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()
        self.TSpecfig = Figure((PlotW, PlotH))
        self.TSpectrumPanel = FigureCanvas(self.TSpecfig)
        fbox.addWidget(self.TSpectrumPanel)
        frame.setLayout(fbox)

        self.slider_tspec = QtWidgets.QScrollBar(QtCore.Qt.Vertical)
        self.slider_tspec.setFocusPolicy(QtCore.Qt.StrongFocus)
        #self.slider_tspec.setEnabled(False)
        self.slider_tspec.valueChanged[int].connect(self.OnTSScroll)
        self.slider_tspec.setRange(1, 5)


        hbox11.addWidget(frame)
        hbox11.addWidget(self.slider_tspec)

        vbox2.addWidget(self.tc_tspec)
        vbox2.addLayout(hbox11)

        self.slider_theta = QtWidgets.QScrollBar(QtCore.Qt.Horizontal)
        self.slider_theta.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_theta.valueChanged[int].connect(self.OnScrollTheta)
        self.slider_theta.setRange(0, 100)
        self.slider_theta.setMinimumWidth(350)
        self.slider_theta.setVisible(False)
        self.tc_imagetheta = QtWidgets.QLabel(self)
        self.tc_imagetheta.setText("4D Data Angle: ")
        self.tc_imagetheta.setVisible(False)
        hbox21 = QtWidgets.QHBoxLayout()
        hbox21.addWidget(self.tc_imagetheta)
        hbox21.addStretch(1)
        hbox21.addWidget(self.slider_theta)
        hbox21.addStretch(1)
        vbox2.addLayout(hbox21)



        #panel 3
        sizer3 = QtWidgets.QGroupBox('Target spectra')
        vbox3 = QtWidgets.QVBoxLayout()


        self.button_addclspec = QtWidgets.QPushButton('Add cluster spectra')
        self.button_addclspec.setMinimumSize (180 , 0)
        self.button_addclspec.clicked.connect( self.OnAddClusterSpectra)
        self.button_addclspec.setEnabled(False)
        vbox3.addWidget(self.button_addclspec)
        self.button_loadtspec = QtWidgets.QPushButton('Load spectra')
        self.button_loadtspec.clicked.connect(self.OnTSpecFromFile)
        self.button_loadtspec.setEnabled(False)
        vbox3.addWidget(self.button_loadtspec)
        self.button_addflat = QtWidgets.QPushButton('Add flat spectrum')
        self.button_addflat.clicked.connect(self.OnFlatTSpec)
        self.button_addflat.setEnabled(False)
        vbox3.addWidget(self.button_addflat)


        self.button_showrgb = QtWidgets.QPushButton('Composite RGB image...')
        self.button_showrgb.clicked.connect(self.OnCompositeRGB)
        self.button_showrgb.setEnabled(False)
        vbox3.addWidget(self.button_showrgb)

        #ToDo: Repair Histogram Value Cutoff if needed. This is covered by "Spectral ROI"
        #self.button_histogram = QtWidgets.QPushButton('Histogram Value Cutoff...')
        #self.button_histogram.clicked.connect(self.OnHistogram)
        #self.button_histogram.setEnabled(False)
        #vbox3.addWidget(self.button_histogram)

        self.button_save = QtWidgets.QPushButton('Save images...')
        self.button_save.clicked.connect(self.OnSave)
        self.button_save.setEnabled(False)
        vbox3.addWidget(self.button_save)

        self.button_calc4d = QtWidgets.QPushButton('Calculate for all angles')
        self.button_calc4d.clicked.connect(self.OnCalc4D)
        self.button_calc4d.setEnabled(False)
        self.button_calc4d.setVisible(False)
        vbox3.addWidget(self.button_calc4d)

        sizer3.setLayout(vbox3)



        #panel 4
        sizer4 = QtWidgets.QGroupBox('Display')
        vbox4 = QtWidgets.QVBoxLayout()


        sb = QtWidgets.QGroupBox('Spectrum')
        vbox41 = QtWidgets.QVBoxLayout()
        vbox41.setSpacing(10)

        self.textctrl_sp1 =  QtWidgets.QLabel(self)
        vbox41.addWidget(self.textctrl_sp1)
        self.textctrl_sp2 =  QtWidgets.QLabel(self)
        vbox41.addWidget(self.textctrl_sp2)



        self.textctrl_sp1.setText('Common Name: ')
        self.textctrl_sp2.setText('RMS Error: ')

        sb.setLayout(vbox41)
        vbox4.addWidget(sb)




        hbox4b = QtWidgets.QHBoxLayout()

        sb = QtWidgets.QGroupBox('Composition Map')
        vbox42 = QtWidgets.QVBoxLayout()


        self.rb_raw = QtWidgets.QRadioButton( 'Raw', self)
        self.rb_fit = QtWidgets.QRadioButton('Fitted',self)
        self.rb_raw.setChecked(True)
        self.rb_raw.toggled.connect(self.OnRBRawFit)


        vbox42.addWidget(self.rb_raw)
        vbox42.addWidget(self.rb_fit)

        # ToDo: something is wrong here. Scalebar not working. Crucial?
        #self.add_scale_cb = QtWidgets.QCheckBox('Scalebar', self)
        ##self.add_scale_cb.toggle()
        #self.add_scale_cb.stateChanged.connect(self.OnShowScale)
        #vbox42.addWidget(self.add_scale_cb)
        sb.setLayout(vbox42)
        hbox4b.addWidget(sb)




        sb = QtWidgets.QGroupBox('Fit Weights')
        vbox43 = QtWidgets.QVBoxLayout()

        self.tc_spfitlist = QtWidgets.QListWidget()
        vbox43.addWidget(self.tc_spfitlist)

        sb.setLayout(vbox43)
        hbox4b.addWidget(sb)


        vbox44 = QtWidgets.QVBoxLayout()
        self.button_removespec = QtWidgets.QPushButton('Remove Spectrum')
        self.button_removespec.clicked.connect(self.OnRemoveSpectrum)
        self.button_removespec.setEnabled(False)
        vbox44.addWidget(self.button_removespec)
        self.button_movespup = QtWidgets.QPushButton('Move Spectrum Up')
        self.button_movespup.clicked.connect(self.OnMoveSpectrumUp)
        self.button_movespup.setEnabled(False)
        vbox44.addWidget(self.button_movespup)
        self.button_movespdown = QtWidgets.QPushButton('Move Spectrum Down')
        self.button_movespdown.clicked.connect(self.OnMoveSpectrumDown)
        self.button_movespdown.setEnabled(False)
        vbox44.addWidget(self.button_movespdown)
        hbox4b.addLayout(vbox44)

        vbox4.addLayout(hbox4b)
        sizer4.setLayout(vbox4)


        hboxB.addLayout(vbox2)
        hboxB.addStretch(1)
        hboxB.addLayout(vbox1)



        hboxT.addWidget(sizer3)
        hboxT.addWidget(sizer4)
        hboxT.addWidget(sizer5)

        vbox.setContentsMargins(20,20,20,20)

        vbox.addStretch(1)
        vbox.addLayout(hboxT)
        vbox.addStretch(3)
        vbox.addLayout(hboxB)
        vbox.addStretch(1)
        self.setLayout(vbox)



#----------------------------------------------------------------------
    def OnTSpecFromFile(self, event):


        #try:
        if True:
            wildcard = "ASCII files (*.csv *.txt *.xas);;CSV files (*.csv);;TXT files (*.txt);;XAS files (*.xas)"
            directory = self.com.path
            filepath, _filter = QtWidgets.QFileDialog.getOpenFileName(self, 'Choose Spectrum file', directory, wildcard)

            filepath = str(filepath)
            if filepath == '':
                return

            self.filename =  os.path.basename(str(filepath))
            directory =  os.path.dirname(str(filepath))
            self.com.path = directory


            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
            self.anlz.read_target_spectrum(filename=filepath)
            self.com.spec_anl_calculated = 1

            self.i_tspec = self.anlz.n_target_spectra
            self.slider_tspec.setMaximum(self.anlz.n_target_spectra)
            self.slider_tspec.setValue(self.i_tspec)

            if self.anlz.tspec_names[self.i_tspec-1].strip() == '':
                self.anlz.tspec_names[self.i_tspec-1] = 'Component '+str(self.i_tspec)

            self.loadTSpectrum()
            self.loadTargetMap()
            self.ShowSpectraList()

            QtWidgets.QApplication.restoreOverrideCursor()

        #except:
        #    QtWidgets.QApplication.restoreOverrideCursor()
        #    QtWidgets.QMessageBox.warning(self, 'Error', 'Spectra not loaded.')


        self.window().refresh_widgets()


#----------------------------------------------------------------------
    def OnFlatTSpec(self, event):

        try:
            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
            self.anlz.read_target_spectrum(flat=True)
            self.com.spec_anl_calculated = 1

            if self.com.spec_anl4D_calculated == 1:
                self.anlz.calculate_targetmaps_4D()

            self.i_tspec = self.anlz.n_target_spectra
            self.slider_tspec.setMaximum(self.anlz.n_target_spectra)
            self.slider_tspec.setValue(self.i_tspec)

            self.loadTSpectrum()
            self.loadTargetMap()
            self.ShowSpectraList()

            QtWidgets.QApplication.restoreOverrideCursor()

        except:
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, 'Error', 'Flat spectrum not loaded.')


        self.window().refresh_widgets()

#----------------------------------------------------------------------
    def OnAddClusterSpectra(self, event):

        try:

            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
            self.anlz.add_cluster_target_spectra()
            self.com.spec_anl_calculated = 1

            self.i_tspec = self.anlz.n_target_spectra
            self.slider_tspec.setMaximum(self.anlz.n_target_spectra)
            self.slider_tspec.setValue(self.i_tspec)


            self.ShowSpectraList()
            self.loadTSpectrum()
            self.loadTargetMap()


            QtWidgets.QApplication.restoreOverrideCursor()

        except:
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, 'Error', 'Cluster spectra not loaded.')


        self.window().refresh_widgets()

#----------------------------------------------------------------------
    def OnCompositeRGB(self, event):

        compimgwin = ShowCompositeRBGmap(self.window(), self.com, self.anlz)
        compimgwin.show()

#----------------------------------------------------------------------
    def OnHistogram(self, event):

        hwin = ShowMapHistogram(self.window(), self.com, self.anlz)
        hwin.show()

#----------------------------------------------------------------------
    def OnSave(self, event):

        savewin = SaveWinP4(self.window())
        savewin.show()

#----------------------------------------------------------------------
    def OnCalc4D(self, event):


        #if True:
        try:

            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
            self.anlz.calculate_targetmaps_4D()
            self.com.spec_anl4D_calculated = 1

            self.i_tspec = self.anlz.n_target_spectra
            self.slider_tspec.setMaximum(self.anlz.n_target_spectra)
            self.slider_tspec.setValue(self.i_tspec)



            self.anlz.target_svd_maps = self.anlz.target_svd_maps4D[self.itheta]
            self.anlz.original_svd_maps = self.anlz.original_svd_maps4D[self.itheta]
            if len(self.anlz.eigenvecs4D) > 0:
                self.anlz.target_pcafit_maps = self.anlz.target_pcafit_maps4D[self.itheta]
                self.anlz.original_fit_maps = self.anlz.original_fit_maps4D[self.itheta]
                self.anlz.target_pcafit_coeffs = self.anlz.target_pcafit_coeffs4D[self.itheta]
                self.anlz.target_pcafit_spectra = self.anlz.target_pcafit_spectra4D[self.itheta]

            self.slider_theta.setVisible(True)
            self.tc_imagetheta.setVisible(True)
            self.slider_theta.setRange(0, self.stk.n_theta-1)
            self.slider_theta.setValue(self.itheta)
            self.tc_imagetheta.setText("4D Data Angle: "+str(self.stk.theta[self.itheta]))

            self.button_calc4d.setEnabled(True)

            self.ShowSpectraList()
            self.loadTSpectrum()
            self.loadTargetMap()



            QtWidgets.QApplication.restoreOverrideCursor()

        except:
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, 'Error', 'Could not calculate 4D spectra.')


        self.window().refresh_widgets()




#----------------------------------------------------------------------
    def Save(self, filename, path, spec_png = True, spec_pdf = False, spec_svg = False, spec_csv = False,
             img_png = True, img_pdf = False, img_svg = False, img_tif = False):

        self.SaveFileName = os.path.join(path,filename)

        try:
            if img_png:
                self.SaveMaps(imgformat=1)
            if img_pdf:
                self.SaveMaps(imgformat=2)
            if img_svg:
                self.SaveMaps(imgformat=3)
            if img_tif:
                self.SaveMaps(imgformat=0, savetif=True)

            if spec_png:
                self.SaveSpectra(imgformat=1)
            if spec_pdf:
                self.SaveSpectra(imgformat=2)
            if spec_pdf:
                self.SaveSpectra(imgformat=3)
            if spec_csv:
                self.SaveSpectra(imgformat=0, savecsv = True)



        except IOError as e:
            if e.strerror:
                err = e.strerror
            else:
                err = e
            QtWidgets.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err)


#----------------------------------------------------------------------
    def SaveSpectra(self, imgformat=1, savecsv = False):


        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        matplotlib.rcParams['pdf.fonttype'] = 42

        colors=['#FF0000','#000000','#FFFFFF']
        spanclrmap=matplotlib.colors.LinearSegmentedColormap.from_list('spancm',colors)

        if savecsv:
            for i in range (self.anlz.n_target_spectra):
                #Save spectra
                tspectrum = self.anlz.target_spectra[i, :]
                fileName_spec = self.SaveFileName+"_Tspectrum_" +str(i+1)+".csv"
                cname = 'Tspectrum_' +str(i+1)
                self.stk.write_csv(fileName_spec, self.stk.ev, tspectrum, cname = cname)

        if imgformat == 0:
            return

        if imgformat == 1:
            ext = 'png'
        elif imgformat == 2:
            ext = 'pdf'
        elif imgformat == 3:
            ext = 'svg'
        suffix = "." + ext

        for i in range (self.anlz.n_target_spectra):
            #Save spectra
            tspectrum = self.anlz.target_spectra[i, :]


            fig = matplotlib.figure.Figure(figsize =(PlotW*1.21, PlotH*0.48))
            canvas = FigureCanvas(fig)
            fig.clf()
            fig.add_axes((0.15,0.15,0.75,0.75))
            axes = fig.gca()


            line1 = axes.plot(self.stk.ev,tspectrum, color='black', label = 'Raw data')

            if self.com.pca_calculated == 1:
                tspectrumfit = self.anlz.target_pcafit_spectra[i, :]
                diff = np.abs(tspectrum-tspectrumfit)
                line2 = axes.plot(self.stk.ev,tspectrumfit, color='green', label = 'Fit')
                line3 = axes.plot(self.stk.ev,diff, color='grey', label = 'Abs(Raw-Fit)')

            fontP = matplotlib.font_manager.FontProperties()
            fontP.set_size('small')

            axes.legend(loc=4, prop = fontP)

            axes.set_xlabel('Photon Energy [eV]')
            axes.set_ylabel('Optical Density')

            fileName_spec = self.SaveFileName+"_Tspectrum_" +str(i+1)+"."+ext
            fig.savefig(fileName_spec)



        #Save combined:

        fig = matplotlib.figure.Figure(figsize =(PlotW, PlotH))
        canvas = FigureCanvas(fig)
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()

        for i in range (self.anlz.n_target_spectra):
            #Save spectra
            tspectrum = self.anlz.target_spectra[i, :]

            line1 = axes.plot(self.stk.ev,tspectrum, color='black', label = 'Raw data')

            if self.com.pca_calculated == 1:
                tspectrumfit = self.anlz.target_pcafit_spectra[i, :]
                diff = np.abs(tspectrum-tspectrumfit)
                line2 = axes.plot(self.stk.ev,tspectrumfit, color='green', label = 'Fit')
                line3 = axes.plot(self.stk.ev,diff, color='grey', label = 'Abs(Raw-Fit)')

        fontP = matplotlib.font_manager.FontProperties()
        fontP.set_size('small')

        axes.legend(loc=4, prop = fontP)

        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Optical Density')

        fileName_spec = self.SaveFileName+"_Tspectra_composite"+"."+ext
        fig.savefig(fileName_spec)






#----------------------------------------------------------------------
    def SaveMaps(self, imgformat=1, savetif = False):


        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        matplotlib.rcParams['pdf.fonttype'] = 42

        colors=['#FF0000','#000000','#FFFFFF']
        spanclrmap=matplotlib.colors.LinearSegmentedColormap.from_list('spancm',colors)


        if imgformat > 0 :

            if imgformat == 1:
                ext = 'png'
            elif imgformat == 2:
                ext = 'pdf'
            elif imgformat == 3:
                ext = 'svg'
            suffix = "." + ext

            for i in range (self.anlz.n_target_spectra):

                #Save composition maps
                if self.showraw == True:
                    tsmapimage = self.anlz.target_svd_maps[:,:,i]
                else:
                    tsmapimage = self.anlz.target_pcafit_maps[:,:,i]

                fig = matplotlib.figure.Figure(figsize =(PlotH, PlotH))
                canvas = FigureCanvas(fig)
                fig.clf()
                axes = fig.gca()

                divider = make_axes_locatable(axes)
                ax_cb = divider.new_horizontal(size="3%", pad=0.03)

                fig.add_axes(ax_cb)
                axes.set_position([0.03,0.03,0.8,0.94])


                min_val = np.min(tsmapimage)
                max_val = np.max(tsmapimage)
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

                im = axes.imshow(np.rot90(tsmapimage), cmap=spanclrmap, vmin = -bound, vmax = bound)
                cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)

                axes.axis("off")


                fileName_img = self.SaveFileName+"_TSmap_" +str(i+1)+"."+ext
                fig.savefig(fileName_img, pad_inches = 0.0)

        else:
            for i in range (self.anlz.n_target_spectra):

                #Save composition maps
                if self.showraw == True:
                    tsmapimage = self.anlz.target_svd_maps[:,:,i]
                else:
                    tsmapimage = self.anlz.target_pcafit_maps[:,:,i]

                fileName_img = self.SaveFileName+"_TSmap_" +str(i+1)+".tif"

                img1 = Image.fromarray(tsmapimage)
                img1.save(fileName_img)


#----------------------------------------------------------------------
    def OnEditSpectraListClick(self):
        item = self.tc_speclist.currentItem()
        self.tc_speclist.editItem(item)
        self.anlz.tspec_names[self.i_tspec-1] = item.data()
        self.loadTSpectrum()

#----------------------------------------------------------------------
    def OnSpectraListClick(self):
        item = self.tc_speclist.currentRow()

        sel = item
        self.i_tspec = sel+1

        if self.com.spec_anl_calculated == 1:
            self.loadTSpectrum()
            self.loadTargetMap()
            self.slider_tspec.setValue(self.i_tspec)

#----------------------------------------------------------------------
    def OnTSScroll(self, value):

        sel = value
        self.i_tspec = sel

        self.tc_speclist.setCurrentRow(self.i_tspec-1)


        if self.com.spec_anl_calculated == 1:
            self.loadTSpectrum()
            self.loadTargetMap()


#----------------------------------------------------------------------
    def OnScrollTheta(self, value):

        if self.com.spec_anl4D_calculated == 0:
            return

        self.itheta = value

        self.anlz.target_svd_maps = self.anlz.target_svd_maps4D[self.itheta]
        self.anlz.original_svd_maps = self.anlz.original_svd_maps4D[self.itheta]
        if len(self.anlz.eigenvecs4D) > 0:
            self.anlz.target_pcafit_maps = self.anlz.target_pcafit_maps4D[self.itheta]
            self.anlz.original_fit_maps = self.anlz.original_fit_maps4D[self.itheta]
            self.anlz.target_pcafit_coeffs = self.anlz.target_pcafit_coeffs4D[self.itheta]
            self.anlz.target_pcafit_spectra = self.anlz.target_pcafit_spectra4D[self.itheta]

        self.tc_imagetheta.setText("4D Data Angle: "+'{0:5.2f}\t'.format(self.stk.theta[self.itheta]))


        self.loadTSpectrum()
        self.loadTargetMap()


        self.window().page0.itheta = self.itheta
        self.window().page0.slider_theta.setValue(self.itheta)

        self.window().page1.itheta = self.itheta
        self.window().page1.slider_theta.setValue(self.itheta)




#----------------------------------------------------------------------
    def OnRBRawFit(self, enabled):
        state = enabled

        if state:
            self.showraw = True
        else:
            self.showraw = False

        if self.com.spec_anl_calculated == 1:
            self.loadTSpectrum()
            self.loadTargetMap()

#----------------------------------------------------------------------
    def OnShowScale(self, state):

        if state == QtCore.Qt.Checked:
            self.show_scale_bar = 1
        else: self.show_scale_bar = 0

        if self.com.spec_anl_calculated == 1:
            self.loadTSpectrum()
            self.loadTargetMap()

#----------------------------------------------------------------------
    def OnRemoveSpectrum(self, event):
        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
        self.anlz.remove_spectrum(self.i_tspec-1)
        self.com.spec_anl_calculated = 1

        self.i_tspec = self.i_tspec-1
        if self.i_tspec<0:
            self.i_tspec=0
        self.slider_tspec.setMaximum(self.anlz.n_target_spectra)
        self.slider_tspec.setValue(self.i_tspec)

        if self.anlz.tspectrum_loaded == 1:
            if self.com.spec_anl4D_calculated == 1:
                self.anlz.calculate_targetmaps_4D()

            self.loadTSpectrum()
            self.loadTargetMap()
            self.ShowSpectraList()
        else:
            self.com.spec_anl_calculated = 0
            self.com.spec_anl4D_calculated = 0
            self.ClearWidgets()

        QtWidgets.QApplication.restoreOverrideCursor()

#----------------------------------------------------------------------
    def OnMoveSpectrumDown(self, event):

        if self.i_tspec < self.anlz.n_target_spectra:
            self.anlz.move_spectrum(self.i_tspec-1, self.i_tspec)

            self.i_tspec += 1
            self.slider_tspec.setValue(self.i_tspec)

            if self.com.spec_anl4D_calculated == 1:
                self.anlz.calculate_targetmaps_4D()

            self.loadTSpectrum()
            self.loadTargetMap()
            self.ShowSpectraList()


#----------------------------------------------------------------------
    def OnMoveSpectrumUp(self, event):

        if self.i_tspec > 1:
            self.anlz.move_spectrum(self.i_tspec-1, self.i_tspec-2)

            if self.com.spec_anl4D_calculated == 1:
                self.anlz.calculate_targetmaps_4D()

            self.i_tspec -= 1
            self.slider_tspec.setValue(self.i_tspec)
            self.loadTSpectrum()
            self.loadTargetMap()
            self.ShowSpectraList()

#----------------------------------------------------------------------
    def ShowFitWeights(self):

        self.tc_spfitlist.clear()

        norm_factor = 100./np.sum(np.absolute(self.anlz.target_pcafit_coeffs[self.i_tspec-1, :]))

        for i in range(self.anlz.numsigpca):
            textitem = '{0}: {1:5.2f} %'.format(i+1, norm_factor
                                                    *abs(self.anlz.target_pcafit_coeffs[self.i_tspec-1, i]))


            self.tc_spfitlist.addItem(textitem)
        self.tc_spfitlist.setCurrentRow(0)

#----------------------------------------------------------------------
    def ShowSpectraList(self):

        self.tc_speclist.clear()

        for i in range(self.anlz.n_target_spectra):
            self.tc_speclist.addItem(self.anlz.tspec_names[i])


#----------------------------------------------------------------------
    def loadTargetMap(self):

        if self.showraw == True:
            tsmapimage = self.anlz.target_svd_maps[:,:,self.i_tspec-1]
        else:
            tsmapimage = self.anlz.target_pcafit_maps[:,:,self.i_tspec-1]


        colors=['#FF0000','#000000','#FFFFFF']

        spanclrmap=matplotlib.colors.LinearSegmentedColormap.from_list('spancm',colors)

        fig = self.mapfig
        fig.clf()

        axes = fig.gca()

        divider = make_axes_locatable(axes)
        ax_cb = divider.new_horizontal(size="3%", pad=0.03)

        fig.add_axes(ax_cb)

        axes.set_position([0.03,0.03,0.8,0.94])


        min_val = np.min(tsmapimage)
        max_val = np.max(tsmapimage)
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

        im = axes.imshow(np.rot90(tsmapimage), cmap=spanclrmap, vmin = -bound, vmax = bound)
        cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)

        axes.axis("off")
        self.MapPanel.draw()



#----------------------------------------------------------------------
    def loadTSpectrum(self):


        if self.anlz.tspectrum_loaded == 0:
            return

        tspectrum = self.anlz.target_spectra[self.i_tspec-1, :]


        fig = self.TSpecfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()


        line1 = axes.plot(self.stk.ev,tspectrum, color='black', label = 'Raw data')

        if self.com.pca_calculated == 1 and hasattr(self.anlz, 'target_pcafit_spectra'):
            tspectrumfit = self.anlz.target_pcafit_spectra[self.i_tspec-1, :]
            diff = np.abs(tspectrum-tspectrumfit)
            line2 = axes.plot(self.stk.ev,tspectrumfit, color='green', label = 'Fit')

            line3 = axes.plot(self.stk.ev,diff, color='grey', label = 'Abs(Raw-Fit)')
        else:
            QtWidgets.QMessageBox.warning(self, 'Error', 'No PCA conducted. Calculate PCA!')

        fontP = matplotlib.font_manager.FontProperties()
        fontP.set_size('small')

        axes.legend(loc=4, prop = fontP)

        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Optical Density')

        self.TSpectrumPanel.draw()

        self.tc_tspec.setText("Target Spectrum: " +
                               self.anlz.tspec_names[self.i_tspec-1])


        self.textctrl_sp1.setText('Common Name: '+
                                    self.anlz.tspec_names[self.i_tspec-1])
        if self.com.pca_calculated == 1 and hasattr(self.anlz, 'target_pcafit_spectra'):
            self.textctrl_sp2.setText('RMS Error: '+ str('{0:7.5f}').format(self.anlz.target_rms[self.i_tspec-1]))

            self.ShowFitWeights()

#----------------------------------------------------------------------
    def ClearWidgets(self):

        fig = self.mapfig
        fig.clf()
        self.MapPanel.draw()
        fig = self.TSpecfig
        fig.clf()
        self.TSpectrumPanel.draw()

        self.tc_tspec.setText("Target Spectrum: ")
        self.tc_speclist.clear()
        self.tc_spfitlist.clear()

        self.com.spec_anl_calculated = 0
        self.com.spec_anl4D_calculated = 0
        self.i_tspec = 1
        self.showraw = True
        self.rb_raw.setChecked(True)

        self.slider_tspec.setValue(self.i_tspec)

        self.button_calc4d.setEnabled(False)
        self.button_calc4d.setVisible(False)


        self.textctrl_sp1.setText('Common Name: \n')
        self.textctrl_sp2.setText('RMS Error: ')

        self.itheta = 0
        self.slider_theta.setVisible(False)
        self.tc_imagetheta.setVisible(False)

        self.window().refresh_widgets()

