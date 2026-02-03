
import os
import numpy as np
from PyQt5 import QtWidgets, uic, QtCore, QtGui
from PyQt5.QtCore import Qt
import matplotlib
import matplotlib.cm
# from mpl_toolkits.axes_grid1 import make_axes_locatable # Unused in view
from ..widgets import SpecFig, ImgFig
from ..dialogs.save import SaveWin

class PagePCACluster(QtWidgets.QWidget):
    def __init__(self, common, data_struct, stack, anlz):
        super(PagePCACluster, self).__init__()
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pagepca.ui'), self)
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
                          'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar'])]

        self.initUI(common, data_struct, stack, anlz)

#----------------------------------------------------------------------
    def initUI(self, common, data_struct, stack, anlz):

        self.spectrum_plotwidget.setBackground("w")

        self.data_struct = data_struct
        self.stk = stack
        self.com = common
        self.anlz = anlz
        self.calcpca = False

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


        self.climgfig = ImgFig(self,self.canvas)
        self.specfig = SpecFig(self, self.spectrum_plotwidget)

        self.slider_cl.valueChanged[int].connect(self.OnScrollCl)
        self.slider_theta.setVisible(False)

        self.pb_copy_img.clicked.connect(self.climgfig.OnCopy)
        self.pb_copy_img.setEnabled(False)
        self.pb_copy_specimg.clicked.connect(self.specfig.OnCopy)
        self.pb_copy_specimg.setEnabled(False)

        #self.MetricCheckBox.toggled.connect(lambda: self.absimgfig.OnMetricScale(self.MetricCheckBox.isChecked(), self.ZeroOriginCheckBox.isChecked(),self.SquarePxCheckBox.isChecked()))
        #self.ZeroOriginCheckBox.toggled.connect(lambda: self.absimgfig.OnMetricScale(self.MetricCheckBox.isChecked(), self.ZeroOriginCheckBox.isChecked(),self.SquarePxCheckBox.isChecked()))
        self.SquarePxCheckBox.toggled.connect(lambda: self.climgfig.OnMetricScale(False, True,self.SquarePxCheckBox.isChecked()))
        self.SquarePxCheckBox.setVisible(True)
        self.ScalebarCheckBox.toggled.connect(lambda: self.climgfig.OnUpdateScale(self.ScalebarCheckBox.isChecked()))

        self.CMCatBox.addItems([self.cmaps[0][0],self.cmaps[1][0],self.cmaps[2][0],self.cmaps[3][0],self.cmaps[4][0],self.cmaps[5][0]])

        self.CMMapBox.addItems(self.cmaps[3][1])
        self.CMCatBox.currentIndexChanged.connect(self.climgfig.OnCatChanged)
        self.CMMapBox.currentIndexChanged.connect(lambda: self.climgfig.OnColormapChange(map=self.CMMapBox.currentText(),num_colors=self.StepSpin.value(),fliplut=self.ColorFlipCheckBox.isChecked()))
        self.CMCatBox.setCurrentIndex(3)
        self.CMMapBox.setCurrentIndex(5)
        self.StepSpin.valueChanged.connect(lambda: self.climgfig.OnColormapChange(map=self.CMMapBox.currentText(),num_colors=self.StepSpin.value(),fliplut=self.ColorFlipCheckBox.isChecked()))
        self.ColorFlipCheckBox.toggled.connect(lambda: self.climgfig.OnColormapChange(map=self.CMMapBox.currentText(),num_colors=self.StepSpin.value(),fliplut=self.ColorFlipCheckBox.isChecked()))
        #self.ROIShapeBox.addItems(["Lasso", "Rectangle", "Circle", "Ellipse", "Polygon", "Histogram"])
        #self.ROIShapeBox.currentTextChanged.connect(self.absimgfig.OnROIShapeChanged)

        self.button_calcpca.clicked.connect( self.OnCalcPCA)
        self.button_calcpca.setEnabled(False)

        self.button_savepca.clicked.connect( self.OnSave)
        self.button_savepca.setEnabled(False)

#----------------------------------------------------------------------
    def OnCalcPCA(self, event):

        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
        self.calcpca = False
        self.selpca = 1
        self.numsigpca = 2
        #self.slidershow.setValue(self.selpca)

        #scrollmax = np.min([self.stk.n_ev, 20])
        #self.slidershow.setMaximum(scrollmax)

        #try:
        self.CalcPCA()
        self.calcpca = True
        #self.loadPCAImage()
        #self.loadPCASpectrum()
        #self.showEvals()
        self.com.pca_calculated = 1
        QtWidgets.QApplication.restoreOverrideCursor()
        self.climgfig.loadPCAImage()
        self.specfig.ClearandReload(self.climgfig)
        #except:
        #    pass
        #    self.com.pca_calculated = 0
        #    QtWidgets.QApplication.restoreOverrideCursor()
         #   QtWidgets.QMessageBox.warning(self, 'Error', 'PCA not calculated.')

        self.window().refresh_widgets()

#----------------------------------------------------------------------
    def CalcPCA(self):

        self.anlz.calculate_pca()

        #Scree plot criterion
        self.numsigpca = self.anlz.numsigpca

        self.npcaspin.setValue(self.numsigpca)

        # cumulative variance
        var = self.anlz.variance[:self.numsigpca].sum()
        self.vartc.setText(str(var.round(decimals=2)*100)+'%')

#----------------------------------------------------------------------
    def CalcPCA4D(self):

        if self.com.stack_4d == 1:
            self.anlz.calculate_pca_4D()
        else:
            return

        #Scree plot criterion
        self.numsigpca = self.anlz.numsigpca
        self.npcaspin.setValue(self.numsigpca)

        # cumulative variance
        var = self.anlz.variance[:self.numsigpca].sum()
        self.vartc.setText(str(var.round(decimals=2)*100)+'%')

        if self.com.spec_anl4D_calculated == 1:
            self.anlz.calculate_targetmaps_4D()

#-----------------------------------------------------------------------
    def OnScrollCl(self, value):
        self.selpca = value
        #self.button_meanflux.setChecked(False)
        #self.mean_visible = 0
        if self.anlz.pca_calculated == 1:
            self.slider_cl.setValue(value)
            image = self.anlz.pcaimages[:, :, self.selpca]  # .copy()
            self.climgfig.draw(image,levels=(-self.anlz.pcaimagebounds[self.selpca],self.anlz.pcaimagebounds[self.selpca]),setpcalabel=True)
            self.specfig.updatePlotDataOnPCA()
            #self.specfig.setLineIndicator(self.iev)
            #self.specfig.LineIndicator.setValue(self.stk.ev[self.iev])


#----------------------------------------------------------------------
    def OnSave(self, event):

        savewin = SaveWin(self.window(), self.com, self.stk)
        savewin.show()
