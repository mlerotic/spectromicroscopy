
import os
import numpy as np
from PyQt5 import QtWidgets, uic, QtCore, QtGui
# import pyqtgraph as pg # Assuming ImgFig uses it internally, but maybe needed for layouts if used directly?
from ..widgets import ImgFig

from ..dialogs.save import SaveWin

class PageLoadData(QtWidgets.QWidget):
    def __init__(self, common, data_struct, stack):
        super(PageLoadData, self).__init__()
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pageloaddata.ui'), self)
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
        #self.scale = 0.000001
        self.data_struct = data_struct
        self.stk = stack
        self.com = common

        self.filename = " "
        self.showflux = True
        self.slider_eng.valueChanged[int].connect(self.OnScrollEng)
        self.slider_theta.valueChanged[int].connect(self.OnScrollTheta)
        self.iev = 0
        self.itheta = 0


        self.absimgfig = ImgFig(self, self.canvas)

        self.button_multiload.clicked.connect( self.OnLoadMulti)
        self.button_multiload.setToolTip('Supported Formats .hdf .hdf5 .ncb .nxs .hdr .stk .tif .tiff .txrm')
        #self.button_multiload.setToolTip('Supported Formats .hdf .hdf5 .ncb .npy .nxs .hdr .stk .tif .tiff .txrm')
        self.button_4d.setToolTip('Supported Formats .hdf5 .ncb')
        self.button_4d.clicked.connect( self.OnLoad4D)

        self.button_sm.setToolTip('Supported Formats .sm, .tif, .xrm')
        self.button_sm.clicked.connect( self.OnBuildStack)

        self.pb_copy_img.clicked.connect(self.absimgfig.OnCopy)
        self.pb_copy_img.setEnabled(False)
        self.pb_copy_specimg.setVisible(False)
        self.MetricCheckBox.toggled.connect(lambda: self.absimgfig.OnMetricScale(self.MetricCheckBox.isChecked(), self.ZeroOriginCheckBox.isChecked(),self.SquarePxCheckBox.isChecked()))
        self.ZeroOriginCheckBox.toggled.connect(lambda: self.absimgfig.OnMetricScale(self.MetricCheckBox.isChecked(), self.ZeroOriginCheckBox.isChecked(),self.SquarePxCheckBox.isChecked()))
        self.SquarePxCheckBox.toggled.connect(lambda: self.absimgfig.OnMetricScale(self.MetricCheckBox.isChecked(), self.ZeroOriginCheckBox.isChecked(),self.SquarePxCheckBox.isChecked()))
        self.SquarePxCheckBox.setVisible(False)
        self.ScalebarCheckBox.toggled.connect(lambda: self.absimgfig.OnUpdateScale(self.ScalebarCheckBox.isChecked()))

        self.CMCatBox.addItems([self.cmaps[0][0],self.cmaps[1][0],self.cmaps[2][0],self.cmaps[3][0],self.cmaps[4][0],self.cmaps[5][0]])
        self.CMMapBox.addItems(self.cmaps[2][1])
        self.CMCatBox.setCurrentIndex(2)
        self.CMMapBox.setCurrentIndex(3)
        self.CMCatBox.currentIndexChanged.connect(self.absimgfig.OnCatChanged)
        self.CMMapBox.currentIndexChanged.connect(lambda: self.absimgfig.OnColormapChange(map=self.CMMapBox.currentText(),num_colors=self.StepSpin.value(),fliplut=self.ColorFlipCheckBox.isChecked()))
        self.StepSpin.valueChanged.connect(lambda: self.absimgfig.OnColormapChange(map=self.CMMapBox.currentText(),num_colors=self.StepSpin.value(),fliplut=self.ColorFlipCheckBox.isChecked()))
        self.ColorFlipCheckBox.toggled.connect(lambda: self.absimgfig.OnColormapChange(map=self.CMMapBox.currentText(),num_colors=self.StepSpin.value(),fliplut=self.ColorFlipCheckBox.isChecked()))
        self.pb_rotate.clicked.connect(self.OnRotate)
        self.pb_mirror.clicked.connect(self.OnMirror)
        self.button_save.clicked.connect( self.OnSave)
        self.button_save.setEnabled(False)
        #self.tc_file.setText('File name')

        #self.tc_path.setText('D:/')
        self.slider_theta.setVisible(False)

#    def keyPressEvent(self, e):
#        if e.key() == 67 and (e.modifiers() & QtCore.Qt.ControlModifier):
#            self.OnCopy()

# ----------------------------------------------------------------------
    def OnSave(self, event):
        savewin = SaveWin(self, self.com, self.stk)
        savewin.show()

# ----------------------------------------------------------------------
    def OnMirror(self):
        if self.com.stack_loaded == 1:
            if self.com.stack_4d == 1:
                self.stk.stack4D = np.flip(self.stk.stack4D, axis=0)
                self.stk.absdata = self.stk.stack4D[:, :, :, self.itheta].copy()
            else:
                self.stk.absdata = np.flip(self.stk.absdata, axis=0)

            if self.com.i0_loaded:
                if self.com.stack_4d:
                    self.stk.od4d = np.flip(self.stk.od4d, axis=0)
                else:
                    self.stk.od3d = np.flip(self.stk.od3d, axis=0)
                    self.stk.od = self.stk.od3d.copy()
                    self.stk.od = np.reshape(self.stk.od, (self.stk.n_rows * self.stk.n_cols, self.stk.n_ev),
                                           order='F')

            self.stk.fill_h5_struct_from_stk()
            if self.com.i0_loaded == 1:
                self.stk.fill_h5_struct_normalization()

            self.OnScrollEng(self.iev)
            # Update/Refresh widgets:
            #if showmaptab:
            #    self.window().tab_map.Clear()
            #    self.window().tab_map.loadData()
            self.window().tab_prep.absimgfig.loadNewImageWithROI()
            self.window().tab_prep.specfig.ClearandReload()
        return

    def OnRotate(self):
        if self.com.stack_loaded == 1:
            if self.com.stack_4d == 1:
                self.stk.stack4D = np.rot90(self.stk.stack4D, 3)
                self.stk.absdata = self.stk.stack4D[:, :, :, self.itheta].copy()
            else:
                self.stk.absdata = np.rot90(self.stk.absdata, 3)

            # Swap x/y constants:
            self.stk.n_cols, self.stk.n_rows = self.stk.n_rows, self.stk.n_cols
            self.stk.x_pxsize, self.stk.y_pxsize = self.stk.y_pxsize, self.stk.x_pxsize
            self.stk.x_start, self.stk.y_start = self.stk.y_start, self.stk.x_start
            self.stk.x_dist, self.stk.y_dist = self.stk.y_dist, self.stk.x_dist

            if self.com.i0_loaded:
                if self.com.stack_4d:
                    self.stk.od4d = np.rot90(self.stk.od4d, 3)
                else:
                    self.stk.od3d =  np.rot90(self.stk.od3d, 3)
                    self.stk.od = self.stk.od3d.copy()
                    self.stk.od = np.reshape(self.stk.od, (self.stk.n_rows * self.stk.n_cols, self.stk.n_ev),
                                           order='F')

            self.stk.fill_h5_struct_from_stk()
            if self.com.i0_loaded == 1:
                self.stk.fill_h5_struct_normalization()

            # Update/Refresh widgets:
            self.OnScrollEng(self.iev)
            self.absimgfig.OnMetricScale(self.MetricCheckBox.isChecked(), self.ZeroOriginCheckBox.isChecked(),self.SquarePxCheckBox.isChecked())
            #if showmaptab:
            #    self.window().tab_map.Clear()
            #    self.window().tab_map.loadData()
            self.window().tab_prep.ix = int(self.stk.n_cols / 2)
            self.window().tab_prep.iy = int(self.stk.n_rows / 2)
            #self.window().tab_prep.loadSpectrum(self.window().tab_prep.ix, self.window().tab_prep.iy)
            self.window().tab_prep.absimgfig.loadNewImageWithROI()
            self.window().tab_prep.specfig.ClearandReload()
        return

#-----------------------------------------------------------------------
    def OnLoadMulti(self, event):

        self.window().LoadStack()


#----------------------------------------------------------------------
    def OnLoad4D(self, event):

        self.window().LoadStack4D()

#----------------------------------------------------------------------
    def OnBuildStack(self, event):

        self.window().BuildStack()

#----------------------------------------------------------------------
    def OnScrollEng(self, value):
        self.iev = value
        if self.com.stack_loaded == 1:
            self.slider_eng.setValue(value)
            image = self.stk.absdata[:, :, int(self.iev)]
            self.absimgfig.draw(image)

#-----------------------------------------------------------------------
    def OnScrollTheta(self, value):
        self.slider_theta.setValue(value)
        self.itheta = value

        self.stk.absdata = self.stk.stack4D[:,:,:,self.itheta].copy()
        image = self.stk.absdata[:, :, int(self.slider_eng.value())].copy()

        self.absimgfig.draw(image)
        #self.window().tab_prep.itheta = self.itheta
        #self.window().tab_prep.slider_theta.setValue(self.itheta)

        #self.window().tab_clus.itheta = self.itheta
        #self.window().tab_clus.slider_theta.setValue(self.itheta)
