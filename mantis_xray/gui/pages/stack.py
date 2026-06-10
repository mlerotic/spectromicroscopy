
import os
import numpy as np
from PyQt5 import QtWidgets, uic, QtCore, QtGui
from PyQt5.QtCore import Qt
import matplotlib
import matplotlib.cm
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
# from mpl_toolkits.axes_grid1 import make_axes_locatable # Used in commented out code?

from ..widgets import SpecFig, ImgFig
from ..dialogs.save import SaveWin
from ..dialogs.crop import MultiCrop
from ..dialogs.alignment import ImageRegistrationDialog, ImageRegistrationManual, ImageRegistrationFFT
from ..dialogs.dark_signal import DarkSignal
from ..dialogs.artefacts import ShowArtefacts
from ..dialogs import DoseCalculation
from ..dialogs.odmap import ShowODMap
# from ..dialogs.colortable import ColorTableFrame

class PageStack(QtWidgets.QWidget):
    def __init__(self, common, data_struct, stack):
        super(PageStack, self).__init__()
        # Load UI from the same directory
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pagestack.ui'), self)
        self.show()
        self.cmaps = [('Perceptually Uniform Sequential', [
            'viridis', 'plasma', 'inferno', 'magma']),
                      ('Sequential', [
                          'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                          'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                          'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']),
                      ('Sequential (2)', [
                          'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
                          'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
                          'hot', 'afmhot', 'gist_heat', 'copper']),
                      ('Diverging', [
                          'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
                          'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']),
                      ('Qualitative', [
                          'Pastel1', 'Pastel2', 'Paired', 'Accent',
                          'Dark2', 'Set1', 'Set2', 'Set3',
                          'tab10', 'tab20', 'tab20b', 'tab20c']),
                      ('Miscellaneous', [
                          'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
                          'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv',
                          'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar', 'Mantis'])]

        self.initUI(common, data_struct, stack)

#----------------------------------------------------------------------
    def initUI(self, common, data_struct, stack):

        self.spectrum_plotwidget.setBackground("w")

        self.data_struct = data_struct
        self.stk = stack
        self.com = common

        self.filename = " "

        self.ix = 0
        self.iy = 0
        self.iev = 0
        self.itheta = 0
        self.showflux = True
        self.show_colorbar = 0

        self.dispbrightness_min = 0
        self.dispbrightness_max = 100
        self.displaygamma = 10.0
        self.defaultdisplay = 1.0

        self.brightness_min = 0.0
        self.brightness_max = 1.0
        self.gamma = 1.0

        self.colortable = "gray"

        self.addroi = 0
        self.showROImask = 0
        self.line = None
        self.ROIpix = None

        self.show_scale_bar = 1
        self.white_scale_bar = 0

        self.movie_playing = 0
        self.mean_visible = 0


        #Panel Preprocess

        #Align stack...
        self.button_align.clicked.connect(self.OnAlignImgsDialog)
        self.button_align.setEnabled(False)

        #Crop stack 3D/4D...
        self.button_multicrop.clicked.connect( self.OnMultiCrop)
        self.button_multicrop.setEnabled(False)

        #Artefacts && Leveling
        self.button_artefacts.clicked.connect( self.OnArtefacts)
        self.button_artefacts.setEnabled(False)

        #Dark signal subtraction...
        self.button_darksig.clicked.connect(self.OnDarkSignal)
        self.button_darksig.setEnabled(False)

        #Save processed stack
        self.button_savestack.clicked.connect(self.OnSaveStack)
        self.button_savestack.setEnabled(False)

        #Panel Normalize
        #Select I0...
        self.button_i0.clicked.connect( self.OnI0histogram)
        self.button_i0.setEnabled(False)

        #I0 from file...
        self.button_i0ffile.clicked.connect(self.OnI0FFile)
        self.button_i0ffile.setEnabled(False)

        #Show I0...
        self.button_showi0.clicked.connect( self.OnShowI0)
        self.button_showi0.setEnabled(False)

        #Use pre-normalized data
        self.button_prenorm.clicked.connect(self.OnPreNormalizedData)
        self.button_prenorm.setEnabled(False)

        #Load Reference Images
        self.button_refimgs.clicked.connect(self.OnRefImgs)
        self.button_refimgs.setEnabled(False)

        self.button_meanflux.clicked.connect( self.OnShowMean)
        self.button_meanflux.setEnabled(False)
        self.button_meanflux.setChecked(False)
        self.button_slideshow.clicked.connect( self.OnSlideshow)
        self.button_slideshow.setEnabled(False)

        #Save images...
        self.button_save.clicked.connect( self.OnSave)
        self.button_save.setEnabled(False)

        #ROI Dose Calculation...
        self.button_ROIdosecalc.clicked.connect( self.OnROI_DoseCalc)
        self.button_ROIdosecalc.setEnabled(False)

        #Spectral ROI...
        self.button_spectralROI.clicked.connect( self.OnSpectralROI)
        self.button_spectralROI.setEnabled(False)

        self.absimgfig = ImgFig(self,self.canvas)
        self.specfig = SpecFig(self, self.spectrum_plotwidget)

        self.slider_eng.valueChanged[int].connect(self.OnScrollEng)
        self.slider_theta.valueChanged[int].connect(self.OnScrollTheta)
        self.slider_theta.setVisible(False)

        self.pb_copy_img.clicked.connect(self.absimgfig.OnCopy)
        self.pb_copy_img.setEnabled(False)
        self.pb_copy_specimg.clicked.connect(self.specfig.OnCopy)
        self.pb_copy_specimg.setEnabled(False)

        self.SquarePxCheckBox.toggled.connect(lambda: self.absimgfig.OnMetricScale(False, True,self.SquarePxCheckBox.isChecked()))
        self.SquarePxCheckBox.setVisible(True)
        self.ScalebarCheckBox.toggled.connect(lambda: self.absimgfig.OnUpdateScale(self.ScalebarCheckBox.isChecked()))

        self.CMCatBox.addItems([self.cmaps[0][0],self.cmaps[1][0],self.cmaps[2][0],self.cmaps[3][0],self.cmaps[4][0],self.cmaps[5][0]])
        self.CMMapBox.addItems(self.cmaps[2][1])
        self.CMCatBox.setCurrentIndex(2)
        self.CMMapBox.setCurrentIndex(3)
        self.CMCatBox.currentIndexChanged.connect(self.absimgfig.OnCatChanged)
        self.CMMapBox.currentIndexChanged.connect(lambda: self.absimgfig.OnColormapChange(map=self.CMMapBox.currentText(),num_colors=self.StepSpin.value(),fliplut=self.ColorFlipCheckBox.isChecked()))
        self.StepSpin.valueChanged.connect(lambda: self.absimgfig.OnColormapChange(map=self.CMMapBox.currentText(),num_colors=self.StepSpin.value(),fliplut=self.ColorFlipCheckBox.isChecked()))
        self.ColorFlipCheckBox.toggled.connect(lambda: self.absimgfig.OnColormapChange(map=self.CMMapBox.currentText(),num_colors=self.StepSpin.value(),fliplut=self.ColorFlipCheckBox.isChecked()))
        self.ROIShapeBox.addItems(["Lasso", "Rectangle", "Circle", "Ellipse", "Polygon", "Histogram"])
        self.ROIShapeBox.currentTextChanged.connect(self.absimgfig.OnROIShapeChanged)
        self.button_lockspectrum.setEnabled(False)
        self.button_clearlastroi.setEnabled(False)
        self.button_mergeroi.setEnabled(True)
        self.button_subtractroi.setEnabled(True)
        self.button_clearspecfig.setEnabled(False)
        self.ROIShapeBox.setEnabled(False)
        self.ROIvisibleCheckBox.setEnabled(False)
        self.ROIvisibleCheckBox.stateChanged.connect(self.absimgfig.OnROIVisibility)
#----------------------------------------------------------------------

    def OnI0FFile(self, event):


        try:
            directory = self.com.path
            wildcard = "I0 spectra (*.csv *.txt *.xas *.hdr);;CSV I0 files (*.csv);;TXT I0 files (*.txt);;XAS I0 files (*.xas);;SDF I0 files (*.hdr)"
            filepath, _filter = QtWidgets.QFileDialog.getOpenFileName(self, 'Choose I0 spectrum file', directory, wildcard)


            filepath = str(filepath)
            if filepath == '':
                return

            self.filename =  os.path.basename(str(filepath))

            basename, extension = os.path.splitext(self.filename)

            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))

            x = self.stk.n_cols
            y = self.stk.n_rows
            z = self.iev

            self.ix = int(x / 2)
            self.iy = int(y / 2)

            if extension == '.hdr':
                self.stk.read_sdf_i0(filepath)

            elif extension in ['.csv', '.txt', '.xas']:
                self.stk.read_stk_i0(filepath, extension)

        except:
            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
            self.com.i0_loaded = 0
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self,'Error',"I0 file not loaded.")
            import sys; print(sys.exc_info())
        else:
            QtWidgets.QApplication.restoreOverrideCursor()
            self.com.i0_loaded = 1
            self.absimgfig.loadNewImageWithROI()
            self.specfig.ClearandReload()
            self.button_i0.disconnect()
            self.button_i0.setText("Reset I0")
            self.button_i0.clicked.connect(self.OnI0Reset)
            self.window().refresh_widgets()
            self.button_i0ffile.setEnabled(False)
            self.button_prenorm.setEnabled(False)
            self.button_refimgs.setEnabled(False)



#----------------------------------------------------------------------
    def OnI0histogram(self, event):
        self.specfig.OnI0Histogram()

# ----------------------------------------------------------------------
    def OnI0Reset(self, event):
        self.specfig.OnI0Reset()

#----------------------------------------------------------------------
    def I0histogramCalculated(self):
        self.com.i0_loaded = 1
        if self.com.stack_4d == 1:
            self.stk.od3d = self.stk.od4d[:,:,:,self.itheta]
            self.stk.od = self.stk.od3d.copy()
            self.stk.od = np.reshape(self.stk.od, (self.stk.n_cols*self.stk.n_rows, self.stk.n_ev), order='F')

        self.absimgfig.loadNewImageWithROI()
        self.specfig.ClearandReload()

        self.window().refresh_widgets()

#----------------------------------------------------------------------
    def OnShowI0(self, event):
        self.specfig.toggleI0Spectrum()
        if self.button_showi0.isChecked():
            self.button_save.setText('Save I0...')
        else:
            self.button_save.setText('Save...')
        #plot = PlotFrame(self, self.stk.evi0,self.stk.i0data)
        #plot.show()

#-----------------------------------------------------------------------
    def OnArtefacts(self, event):
        artefacts = ShowArtefacts(self.window(), self.com, self.stk)
        artefacts.show()

    # ----------------------------------------------------------------------
    def OnPreNormalizedData(self, event):
        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
        self.stk.UsePreNormalizedData()

        self.com.i0_loaded = 1
        #self.loadSpectrum(self.ix, self.iy)
        self.absimgfig.loadNewImageWithROI()
        self.specfig.ClearandReload()
        #self.parent.I0histogramCalculated()
        self.button_i0.disconnect()
        self.button_i0.setText("Reset I0")
        self.button_i0.clicked.connect(self.OnI0Reset)
        self.window().refresh_widgets()
        self.button_showi0.setEnabled(False)
        self.button_i0ffile.setEnabled(False)
        self.button_prenorm.setEnabled(False)
        self.button_refimgs.setEnabled(False)
        QtWidgets.QApplication.restoreOverrideCursor()
#-----------------------------------------------------------------------
    def OnRefImgs(self, event):

        #Load .xrm reference images
        try:
        #if True:
            wildcard = "Reference images (*.xrm)"

            filepaths, _filter = QtWidgets.QFileDialog.getOpenFileNames(self, 'Select reference files', '', wildcard)


            #Check reference files
            if len(filepaths) != self.stk.n_ev:
                QtWidgets.QMessageBox.warning(self,'Error',"Wrong number of Reference image files.")
                return

            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))


            self.stk.read_xrm_ReferenceImages(filepaths)

            x=self.stk.n_cols
            y=self.stk.n_rows
            self.ix = int(x/2)
            self.iy = int(y/2)

            self.com.i0_loaded = 1
            self.absimgfig.loadNewImageWithROI()
            self.specfig.ClearandReload()
            self.button_i0.disconnect()
            self.button_i0.setText("Reset I0")
            self.button_i0.clicked.connect(self.OnI0Reset)
            self.window().refresh_widgets()
            self.button_i0ffile.setEnabled(False)
            self.button_prenorm.setEnabled(False)
            self.button_refimgs.setEnabled(False)
            QtWidgets.QApplication.restoreOverrideCursor()


        except:
            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
            self.com.i0_loaded = 0
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self,'Error',"Reference image file not loaded.")
            import sys; print(sys.exc_info())

#----------------------------------------------------------------------
    def OnSaveStack(self, event):

        self.window().SaveProcessedStack()

#----------------------------------------------------------------------
    def OnSave(self, event):
        i0_mode = self.button_showi0.isChecked()
        savewin = SaveWin(self, self.com, self.stk, i0_mode=i0_mode)
        if i0_mode:
            savewin.setWindowTitle('Save I0...')
        savewin.show()

# ----------------------------------------------------------------------
    def OnMultiCrop(self, evt):
        multicropwin = MultiCrop(self.window(), self.com, self.stk)
        multicropwin.show()
# ----------------------------------------------------------------------
    def OnAlignImgsDialog(self, event):
        imgregwin = ImageRegistrationDialog(self.window(),self.com)
        imgregwin.show()
#----------------------------------------------------------------------
    def OnAlignImgs(self, event):
        imgregwin = ImageRegistrationManual(self.window(), self.com, self.stk)
        imgregwin.show()

# ----------------------------------------------------------------------
    def OnAlignImgs2(self, event):
        imgreg2 = ImageRegistrationFFT(self.window(), self.com, self.stk)
        imgreg2.show()

#----------------------------------------------------------------------
    def OnDarkSignal(self, event):
        dswin = DarkSignal(self.window(), self.com, self.stk)
        dswin.show()

#----------------------------------------------------------------------
    def OnShowMean(self):
        if (self.com.stack_loaded == 1):
            if (self.mean_visible == 1):
                self.mean_visible = 0
                self.OnScrollEng(self.iev)
                return
            self.mean_visible = 1
            if self.com.i0_loaded == 0:
                # Show flux mean
                image = np.nanmean(self.stk.absdata, axis=2)
            else:
                # Show OD mean
                image = np.nanmean(self.stk.od3d, axis=2)
            self.absimgfig.draw(image)


#----------------------------------------------------------------------
    def OnSlideshow(self, event):

        if (self.com.stack_loaded == 1) and (self.addroi == 0):

            if (self.movie_playing == 1):
                self.movie_playing = 0
                self.button_meanflux.setEnabled(True)
                return
            self.button_meanflux.setEnabled(False)
            self.button_slideshow.setText("Stop stack movie")
            self.movie_playing = 1
            # #self.loadSpectrum(self.ix, self.iy)
            self.updateRate = int(max((2,np.rint(5000/self.stk.n_ev)))) # 5ms min update rate or approximately 5 seconds total duration
            def displayNextImage():
                QtWidgets.qApp.processEvents()
                if self.movie_playing == 0:
                    self.button_slideshow.setText("Play stack movie")
                    return
                self.OnScrollEng(self.iev)
                # at end of stack, start from scratch
                if self.iev == self.stk.n_ev-1:
                    self.iev = 0
                else:
                    self.iev = (self.iev + 1)
                QtCore.QTimer.singleShot(self.updateRate, displayNextImage)
            displayNextImage()
                #self.loadSpectrum(self.ix, self.iy)
#-----------------------------------------------------------------------
    def OnScrollEng(self, value):
        self.iev = value
        self.button_meanflux.setChecked(False)
        self.mean_visible = 0
        if self.com.stack_loaded == 1:
            self.slider_eng.setValue(value)
            if self.defaultdisplay == 1.0:
                # use a pointer to the data not a copy
                if self.com.i0_loaded == 0:
                    # Show flux image
                    image = self.stk.absdata[:, :, self.iev]  # .copy()
                else:
                    # Show OD image
                    image = self.stk.od3d[:, :, self.iev]  # .copy()
            else:
                # Adjustment to the data display setting has been made so make a copy
                #if self.showflux:
                if self.com.i0_loaded == 0:
                    image = self.stk.absdata[:, :, self.iev].copy()
                else:
                    image = self.stk.od3d[:, :, self.iev].copy()
            self.absimgfig.draw(image)
            #self.specfig.setLineIndicator(self.iev)
            self.specfig.LineIndicator.setValue(self.stk.ev[self.iev])

#----------------------------------------------------------------------
    def OnScrollTheta(self, value):
        self.itheta = value
        self.button_meanflux.setChecked(False)
        self.mean_visible = 0
        if self.com.stack_loaded == 1:
            self.slider_theta.setValue(value)
            self.stk.absdata = self.stk.stack4D[:, :, :, self.itheta].copy()
            if self.com.i0_loaded == 0:
                image = self.stk.absdata[:, :, int(self.slider_eng.value())].copy()
            else:
                self.stk.od3d = self.stk.od4d[:, :, :, self.itheta].copy()
                image = self.stk.od3d[:, :, int(self.slider_eng.value())].copy()
            self.stk.calc_histogram()
            #self.tc_imagetheta.setText("4D Data Angle: "+str(self.stk.theta[self.itheta]))
            #self.p1.setTitle("<center>Image at {0:5.2f} eV and {1:5.1f}°</center>".format(float(self.stk.ev[self.iev]),
            #                                                                 float(self.stk.theta[self.itheta])))
            self.absimgfig.draw(image)
            self.specfig.ClearandReload()

#-----------------------------------------------------------------------
    def OnPointSpectrum(self, evt):
        x = evt.xdata
        y = evt.ydata

        if (self.com.stack_loaded == 1) and (self.addroi == 0):
            if x < self.stk.ev[0]:
                sel_ev = 0
            elif x > self.stk.ev[self.stk.n_ev-1]:
                sel_ev = self.stk.n_ev-1
            else:
                indx = np.abs(self.stk.ev - x).argmin()
                sel_ev = indx

            self.iev = sel_ev

            self.loadSpectrum(self.ix, self.iy)
            self.loadImage()

            self.slider_eng.setValue(self.iev)


#----------------------------------------------------------------------
    def OnRb_fluxod(self, enabled):

        state = enabled

        if state:
            self.showflux = True
        else:
            self.showflux = False

        self.ResetDisplaySettings()
        self.loadImage()

#----------------------------------------------------------------------
    def OnShowScale(self,state):

        if state == QtCore.Qt.Checked:
            self.show_scale_bar = 1
        else:
            self.show_scale_bar = 0


        if self.com.stack_loaded == 1:
            self.loadImage()

#----------------------------------------------------------------------
    def OnWhiteScale(self,state):


        if state == QtCore.Qt.Checked:
            self.white_scale_bar = 1
        else:
            self.white_scale_bar = 0

        self.com.white_scale_bar = self.white_scale_bar

        if self.com.stack_loaded == 1:
            self.loadImage()
            #self.window().tab_load.ShowImage()

#----------------------------------------------------------------------
    def OnShowColBar(self, state):

        if state == QtCore.Qt.Checked:
            self.show_colorbar = 1
        else:
            self.show_colorbar = 0

        if self.com.stack_loaded == 1:
            self.loadImage()

#----------------------------------------------------------------------
    def OnResetDisplaySettings(self, event):

        self.ResetDisplaySettings()
        self.loadImage()

#----------------------------------------------------------------------
    def OnScrollBrightnessMin(self, value):

        self.dispbrightness_min = value

        self.brightness_min = float(self.dispbrightness_min)/100.0

        self.defaultdisplay = 0.0

        self.tc_min.setText('Minimum: \t{0:5d}%'.format(int(100*self.brightness_min)))

        if self.com.stack_loaded == 1:
            self.loadImage()

#----------------------------------------------------------------------
    def OnScrollBrightnessMax(self, value):

        self.dispbrightness_max = value

        self.brightness_max = float(self.dispbrightness_max)/100.0

        self.defaultdisplay = 0.0

        self.tc_max.setText('Maximum:{0:5d}%'.format(int(100*self.brightness_max)))

        if self.com.stack_loaded == 1:
            self.loadImage()

#----------------------------------------------------------------------
    def OnScrollGamma(self, value):

        self.displaygamma = value

        self.gamma = float(self.displaygamma)/10.0

        self.defaultdisplay = 0.0

        self.tc_gamma.setText('Gamma:  \t{0:5.2f}'.format(self.gamma))

        if self.com.stack_loaded == 1:
            self.loadImage()
#----------------------------------------------------------------------
    def OnDespike(self, evt):

        image = self.stk.absdata[:,:,self.iev]

        image = self.stk.despike(image)


        self.stk.data_struct.exchange.data = self.stk.absdata

        if self.com.i0_loaded:
            self.stk.calculate_optical_density()

        self.loadImage()

#----------------------------------------------------------------------
    def OnSetColorTable(self, event):

        # colorwin = ColorTableFrame(self.window())
        # colorwin.show()
        pass

# -----------------------------------------------------------------------
#     def loadNewImage(self):
#         fig = self.absimgfig
#         fig.clear()
#         # if self.defaultdisplay == 1.0:
#         #     # use a pointer to the data not a copy
#         #     if self.showflux:
#         #         # Show flux image
#         #         image = self.stk.absdata[:, :, self.iev]  # .copy()
#         #     else:
#         #         # Show OD image
#         #         image = self.stk.od3d[:, :, self.iev]  # .copy()
#         # else:
#         #     # Adjustment to the data display setting has been made so make a copy
#         #     if self.showflux:
#         #         image = self.stk.absdata[:, :, self.iev].copy()
#         #     else:
#         #         image = self.stk.od3d[:, :, self.iev].copy()
#         fig.loadData()
#         #print("setimage")
#         #fig.OnScrollEng(self.iev)
#         #fig.i_item.setImage(image)
#         #im = axes.imshow(np.rot90(image), cmap=matplotlib.cm.get_cmap(self.colortable))
#         # self.tc_imageeng.setText("Image at {0:5.2f} eV".format(float(self.stk.ev[self.iev])))
#-----------------------------------------------------------------------
    def loadImage(self):
        # Placeholder or implemented via Widgets
        self.absimgfig.loadNewImageWithROI()

#-----------------------------------------------------------------------
    def loadSpectrum(self, xpos, ypos):
        if self.com.i0_loaded == 1:
            self.spectrum = self.stk.od3d[int(xpos),int(ypos), :]
        else:
            self.spectrum = self.stk.absdata[int(xpos),int(ypos), :]

        try:
            # Update SpecFig (PyQtGraph)
            # Assuming items[1] is the main PlotCurveItem as managed by SpecFig
            if len(self.specfig.pi.items) > 1:
                self.specfig.pi.items[1].setData(self.stk.ev, self.spectrum)
            
            self.tc_spec.setText('Spectrum at pixel [{0}, {1}] or position [{2:5.2f}, {3:5.2f}]'.format(
                str(xpos),  str(ypos), float(self.stk.x_dist[int(xpos)]), float(self.stk.y_dist[int(ypos)])))
        except Exception as e:
            print(f"Error updating spectrum plot: {e}")

#----------------------------------------------------------------------
    def ResetDisplaySettings(self):

        self.defaultdisplay = 1.0

        self.dispbrightness_min = 0
        self.dispbrightness_max = 100
        self.displaygamma = 10.0

        self.brightness_min = 0.0
        self.brightness_max = 1.0
        self.gamma = 1.0

#----------------------------------------------------------------------
    def CalcROISpectrum(self):

        locked_spectrum = getattr(self.specfig, "last_locked_roi_spectrum", None)
        if locked_spectrum is not None:
            locked_y = np.asarray(locked_spectrum[1])
        else:
            locked_y = None

        self.ROIspectrum = np.zeros((self.stk.n_ev))

        # Get ROI indices from the live ROI if available; otherwise reuse the locked spectrum directly.
        roi_mask = np.sum(self.absimgfig.ROIrgba, axis=2) > 0
        if not np.any(roi_mask):
            if locked_y is not None and locked_y.size == self.stk.n_ev:
                self.ROIspectrum = locked_y.copy()
                return True
            return False
        indices = np.where(roi_mask)
        numroipix = len(indices[0])
        
        if numroipix == 0:
            if locked_y is not None and locked_y.size == self.stk.n_ev:
                self.ROIspectrum = locked_y.copy()
                return True
            # No ROI defined
            return False

        if self.com.i0_loaded == 1:
            for ie in range(self.stk.n_ev):
                thiseng = self.stk.od3d[:, :, ie]
                self.ROIspectrum[ie] = np.sum(thiseng[indices])/numroipix
        else:
            for ie in range(self.stk.n_ev):
                thiseng = self.stk.absdata[:, :, ie]
                self.ROIspectrum[ie] = np.sum(thiseng[indices])/numroipix

        return True

#-----------------------------------------------------------------------
    def ShowROISpectrum(self):

        self.CalcROISpectrum()
        # This calls self.specfig.clf() etc in original code.
        # I need to adapt to SpecFig widget or copy logic.
        pass

#-----------------------------------------------------------------------
    def OnSelectLasso(self,verts):
       pass

#-----------------------------------------------------------------------
    def CalcROI_I0Spectrum(self):
        pass

#----------------------------------------------------------------------
    def OnROI_DoseCalc(self, event):
        roi_spectra = self.specfig.get_roi_spectra_for_dose()
        if not roi_spectra:
            self.CalcROISpectrum()

        if roi_spectra:
            self.ROIspectrum = np.asarray(roi_spectra[0][1]).copy()

        dosewin = DoseCalculation(self, self.stk, self.ROIspectrum, roi_spectra=roi_spectra)
        dosewin.exec_()

#----------------------------------------------------------------------
    def OnSpectralROI(self, evt):
        specroiwin = ShowODMap(self.window(), self.com, self.data_struct, self.stk)
        specroiwin.show()


#----------------------------------------------------------------------
    def point_in_poly(self, x, y, polyx, polyy):
        n = len(polyx)
        inside = False
        p1x, p1y = polyx[0], polyy[0]
        for i in range(1, n + 1):
            p2x, p2y = polyx[i % n], polyy[i % n]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y
        return inside
