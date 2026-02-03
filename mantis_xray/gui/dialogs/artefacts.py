import os
import numpy as np
from PyQt5 import QtWidgets, QtCore, uic
import pyqtgraph as pg
from scipy.interpolate import griddata


class ShowArtefacts(QtWidgets.QDialog):
    def __init__(self, parent, common, stack):
        QtWidgets.QWidget.__init__(self, parent)
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'showartefacts.ui'), self)
        self.parent = parent
        self.com = common
        self.stack = stack
        self.i0_h_warnflag = False
        self.i0_v_warnflag = False
        self.pglayout = pg.GraphicsLayout(border=None)
        self.canvas.setBackground("w") # canvas is a pg.GraphicsView widget
        self.canvas.setCentralWidget(self.pglayout)

        self.vb = self.pglayout.addViewBox()
        self.vb.setAspectLocked()
        self.i_item = pg.ImageItem(border="k",parent= self)


        self.button_ok.clicked.connect(self.OnAccept)
        self.button_cancel.clicked.connect(self.close)
        self.bg_level.toggled.connect(self.ShowImage)
        self.remove_outliers.toggled.connect(self.ShowImage)
        self.cb_h.stateChanged.connect(self.ShowImage)
        self.cb_v.stateChanged.connect(self.ShowImage)
        # I initially wanted to use Qt.QueuedConnection to create a non-blocking Slider. It is good enough.
        self.weight_slider.valueChanged.connect(self.ShowCalcImage)
        self.weight_slider_2.valueChanged.connect(self.ShowCalcImage)
        self.slider_eng.sliderPressed.connect(self.ShowImage)
        self.slider_eng.sliderReleased.connect(self.ShowImage)
        self.slider_eng.valueChanged[int].connect(self.OnScrollEng)

        self.setWindowTitle('Artefacts & Leveling')
        self.slider_eng.setRange(0, self.stack.n_ev - 1)

        self.vb.setMouseEnabled(x=False, y=False)
        self.vb.addItem(self.i_item, ignoreBounds=False)

        self.stack.absdata_level = self.stack.absdata.copy()
        if self.com.i0_loaded:
            self.rb_median_i0.setCheckable(True)
            self.rb_median_i0.setChecked(True)
            self.label_4.setText('I0 mask contribution:')
        else:
            self.sliderWidget.hide()
        self.rb_median.toggled.connect(self.ShowImage)
        self.OnScrollEng(0)

    def ShowCalcImage(self):
        a = self.stack.absdata[:, :, self.slider_eng.value()].astype('float64')
        if self.bg_level.isChecked():
            if any([self.cb_h.isChecked(),self.cb_v.isChecked()]):
                a, wf = self.LevelCalc(a)
                if self.rb_median_i0.isChecked():
                    self.label_3.setText(str('{:d}').format(int(wf))+' %')
        if self.remove_outliers.isChecked():
            a, wf = self.OutlierCalc(a)
            self.label_7.setText(str('<p>exceeding background level &plusmn; {:d} * &sigma;</p>').format(int(wf)))
        self.i_item.setImage(a)

    def CorrectionArray(self,array,axis,wf,mask = None,):
        # calculate median along given axis, ignoring nans if present
        array_i = np.nanmedian(array, axis=axis, keepdims=False)
        if mask is not None and wf != 0.0:
            array = np.where(mask, array, np.nan) # replace masked pixels by nans
            # calculate median along given axis, ignoring nans if present
            array = np.nanmedian(array, axis=axis, keepdims=False)
            bg_median = np.nanmedian(array, axis=0, keepdims=False)
        if mask is None or wf == 0.0:
            # returns the median filtered array neglecting the i0 mask
            return (array_i / np.nanmedian(array_i, axis=0, keepdims=False))
        elif wf == 1:
            return ((array / bg_median)*(wf))
        else:
            return ((array / bg_median)*(wf) + (array_i / bg_median)*(1-wf))

    def LevelCalc(self,a,final=False):
        wf = 0.0
        mask = None
        factor_v = 1
        factor_h = 1
        if self.cb_h.isChecked():
            wf = (self.weight_slider.value()) / 100
            if self.rb_median_i0.isChecked():
                mask = self.stack.i0_mask
                if final:
                    mask = mask[:, :, None]
            factor_h = self.CorrectionArray(a, 0, float(wf), mask)[None, :]
            #a += diff[None, :]
        if self.cb_v.isChecked():
            wf = (self.weight_slider.value()) / 100
            if self.rb_median_i0.isChecked():
                mask = self.stack.i0_mask[:, :]
                if final:
                    mask = mask[:, :, None]
            factor_v = self.CorrectionArray(a, 1, float(wf), mask)[: ,None]
        a = a/factor_h
        a = a/factor_v
        return(a,int(wf*100))

    def CorrectOutliers(self,array,wf,mask = None,):
        # calculate median ignoring nans if present
        if mask is not None:
            array_nan = np.where(mask, array, np.nan) # replace masked pixels by nans
            std = np.nanstd(array_nan, keepdims=False)
            bg_median = np.nanmedian(array_nan, keepdims=False)
        else:
            std = np.nanstd(array, keepdims=False)
            bg_median = np.nanmedian(array, keepdims=False)
        array[(array < (bg_median-wf*int(std))) | (array > (bg_median+wf*int(std)))] = np.nan
        x_idx, y_idx = np.meshgrid(np.arange(0, array.shape[1]),np.arange(0, array.shape[0]))
        array = np.ma.masked_invalid(array) # mask nans
        x = x_idx[~array.mask]
        y = y_idx[~array.mask]
        z = array[~array.mask]
        try:
            returnarray = griddata((x, y), z.ravel(),(x_idx, y_idx), method='nearest') #interpolate array of counts
        except IndexError:
            pass
        return (returnarray)


    def OutlierCalc(self,a,final=False):
        wf = self.weight_slider_2.value()
        mask = None
        #diff_v = 0
        #diff_h = 0
        if self.com.i0_loaded:
            mask = self.stack.i0_mask
            #if final:
            #    mask = mask[:, :, None]
        if final:
            for ev in range(a.shape[-1]):
                a[:, :, ev] = self.CorrectOutliers(a[:,:,ev], wf, mask)
        else:
            a = self.CorrectOutliers(a, wf, mask)
        return(a,wf)
# ----------------------------------------------------------------------
    def OnScrollEng(self, value):
        self.slider_eng.setValue(value)
        self.iev = value

        self.ShowImage()
    def ShowImage(self):
        if (self.slider_eng.isSliderDown()):
            self.i_item.setImage(self.stack.absdata[:, :, int(self.iev)])
        elif any([self.remove_outliers.isChecked(), self.bg_level.isChecked()]):
            self.ShowCalcImage()
        else:
            self.i_item.setImage(self.stack.absdata[:, :, int(self.iev)])

            self.label_3.setText(str(''))
        self.groupBox.setTitle(str('Stack Browser | Image at {0:5.2f} eV').format(float(self.stack.ev[self.iev])))
#----------------------------------------------------------------------
    def OnAccept(self, evt):
        a, wf = self.LevelCalc(self.stack.absdata.astype('float64'),final=True)
        if self.remove_outliers.isChecked():
            a, wf = self.OutlierCalc(a,final=True)
            self.parent.page1.specfig.OnI0Reset()  # Reset I0 on PageStack
            self.parent.window().refresh_widgets()
        if not np.array_equal(self.stack.absdata, a, equal_nan=False):
            self.stack.absdata = a
            self.parent.page0.absimgfig.loadNewImage() # Load new image on PageLoadData
            self.parent.page1.specfig.OnI0Reset() # Reset I0 on PageStack
            self.parent.window().refresh_widgets()
        self.close()
