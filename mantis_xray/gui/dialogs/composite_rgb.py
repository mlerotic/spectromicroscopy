
import os
import numpy as np
from PyQt5 import QtWidgets, QtGui, QtCore
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from ...core.constants import PlotH

class ShowCompositeRBGmap(QtWidgets.QDialog):

    def __init__(self, parent, common, analz):
        QtWidgets.QWidget.__init__(self, parent)

        self.parent = parent


        self.resize(630, 400)
        self.setWindowTitle('Composite RBG Map')

        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)



        self.com = common
        self.anlz = analz


        self.show_info = 0


        self.n_cols = self.anlz.stack.n_cols
        self.n_rows = self.anlz.stack.n_rows

        self.rgbimage = np.zeros((self.n_cols, self.n_rows, 3), dtype=float)

        self.minr = 0
        self.maxr = 100
        self.weightr = 100
        self.ming = 0
        self.maxg = 100
        self.weightg = 100
        self.minb = 0
        self.maxb = 100
        self.weightb = 100

        self.r_spec = 0
        self.g_spec = 1
        self.b_spec = 2


        sizer1 = QtWidgets.QGroupBox('Red spectrum')

        fgs1 = QtWidgets.QGridLayout()

        r = QtWidgets.QLabel(self)
        r.setText('Red')
        rl = QtWidgets.QLabel(self)
        rl.setText('Limits')
        rw = QtWidgets.QLabel(self)
        rw.setText('Weight')


        self.combor = QtWidgets.QComboBox(self)
        self.combor.addItems(self.anlz.tspec_names)
        self.combor.activated[int].connect(self.OnSelectR)

        #self.combor.SetToolTip(wx.ToolTip("select spectrum from dropdown-list"))
        self.combor.setCurrentIndex(self.r_spec)

        hbox12 = QtWidgets.QHBoxLayout()


        self.tcrmin = QtWidgets.QSpinBox()
        self.tcrmin.setRange(0,100)
        self.tcrmin.setValue(0)
        self.tcrmin.valueChanged[int].connect(self.OnLimitMinR)

        self.tcrmax = QtWidgets.QSpinBox()
        self.tcrmax.setRange(0,100)
        self.tcrmax.setValue(100)
        self.tcrmax.valueChanged[int].connect(self.OnLimitMaxR)

        hbox12.addWidget(self.tcrmin)
        hbox12.addWidget(self.tcrmax)

        self.tcrweight = QtWidgets.QSpinBox()
        self.tcrweight.setRange(0,100)
        self.tcrweight.setValue(100)
        self.tcrweight.valueChanged[int].connect(self.OnWeightR)

        fgs1.addWidget(r, 0, 0)
        fgs1.addWidget(self.combor, 0, 1)
        fgs1.addWidget(rl, 1, 0)
        fgs1.addLayout(hbox12, 1, 1)
        fgs1.addWidget(rw, 2, 0)
        fgs1.addWidget(self.tcrweight, 2, 1)

        sizer1.setLayout(fgs1)



        sizer2 = QtWidgets.QGroupBox('Green spectrum')

        fgs2 = QtWidgets.QGridLayout()

        g = QtWidgets.QLabel(self)
        g.setText('Green')
        gl = QtWidgets.QLabel(self)
        gl.setText('Limits')
        gw = QtWidgets.QLabel(self)
        gw.setText('Weight')


        self.combog = QtWidgets.QComboBox(self)
        self.combog.addItems(self.anlz.tspec_names)
        self.combog.activated[int].connect(self.OnSelectG)

        #self.combor.SetToolTip(wx.ToolTip("select spectrum from dropdown-list"))
        self.combog.setCurrentIndex(self.g_spec)

        hbox12 = QtWidgets.QHBoxLayout()


        self.tcgmin = QtWidgets.QSpinBox()
        self.tcgmin.setRange(0,100)
        self.tcgmin.setValue(0)
        self.tcgmin.valueChanged[int].connect(self.OnLimitMinG)

        self.tcgmax = QtWidgets.QSpinBox()
        self.tcgmax.setRange(0,100)
        self.tcgmax.setValue(100)
        self.tcgmax.valueChanged[int].connect(self.OnLimitMaxG)

        hbox12.addWidget(self.tcgmin)
        hbox12.addWidget(self.tcgmax)

        self.tcgweight = QtWidgets.QSpinBox()
        self.tcgweight.setRange(0,100)
        self.tcgweight.setValue(100)
        self.tcgweight.valueChanged[int].connect(self.OnWeightG)


        fgs2.addWidget(g, 0, 0)
        fgs2.addWidget(self.combog, 0, 1)
        fgs2.addWidget(gl, 1, 0)
        fgs2.addLayout(hbox12, 1, 1)
        fgs2.addWidget(gw, 2, 0)
        fgs2.addWidget(self.tcgweight, 2, 1)


        sizer2.setLayout(fgs2)



        sizer3 = QtWidgets.QGroupBox('Blue spectrum')

        fgs3 = QtWidgets.QGridLayout()

        b = QtWidgets.QLabel(self)
        b.setText('Blue')
        bl = QtWidgets.QLabel(self)
        bl.setText('Limits')
        bw = QtWidgets.QLabel(self)
        bw.setText('Weight')


        self.combob = QtWidgets.QComboBox(self)
        self.combob.addItems(self.anlz.tspec_names)
        self.combob.activated[int].connect(self.OnSelectB)

        #self.combor.SetToolTip(wx.ToolTip("select spectrum from dropdown-list"))
        self.combob.setCurrentIndex(self.b_spec)

        hbox12 = QtWidgets.QHBoxLayout()


        self.tcbmin = QtWidgets.QSpinBox()
        self.tcbmin.setRange(0,100)
        self.tcbmin.setValue(0)
        self.tcbmin.valueChanged[int].connect(self.OnLimitMinB)

        self.tcbmax = QtWidgets.QSpinBox()
        self.tcbmax.setRange(0,100)
        self.tcbmax.setValue(100)
        self.tcbmax.valueChanged[int].connect(self.OnLimitMaxB)

        hbox12.addWidget(self.tcbmin)
        hbox12.addWidget(self.tcbmax)

        self.tcbweight = QtWidgets.QSpinBox()
        self.tcbweight.setRange(0,100)
        self.tcbweight.setValue(100)
        self.tcbweight.valueChanged[int].connect(self.OnWeightB)


        fgs3.addWidget(b, 0, 0)
        fgs3.addWidget(self.combob, 0, 1)
        fgs3.addWidget(bl, 1, 0)
        fgs3.addLayout(hbox12, 1, 1)
        fgs3.addWidget(bw, 2, 0)
        fgs3.addWidget(self.tcbweight, 2, 1)


        sizer3.setLayout(fgs3)




        vbox = QtWidgets.QVBoxLayout()
        hbox1 = QtWidgets.QHBoxLayout()
        vbox1 = QtWidgets.QVBoxLayout()

        vbox1.addWidget(sizer1)
        vbox1.addWidget(sizer2)
        vbox1.addWidget(sizer3)



        self.show_info_cb = QtWidgets.QCheckBox( 'Show Info on the Image', self)
        self.show_info_cb.stateChanged.connect(self.OnShowInfo)
        vbox1.addWidget(self.show_info_cb)

        hbox1.addLayout(vbox1)

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()
        self.RGBImagefig = Figure((PlotH, PlotH))
        self.RGBImagePanel = FigureCanvas(self.RGBImagefig)
        fbox.addWidget(self.RGBImagePanel)
        frame.setLayout(fbox)

        hbox1.addWidget(frame)


        vbox.addLayout(hbox1)


        hbox2 = QtWidgets.QHBoxLayout()

        button_save = QtWidgets.QPushButton('Save image')
        button_save.clicked.connect( self.OnSave)
        hbox2.addWidget(button_save)

        button_close = QtWidgets.QPushButton('Dismiss')
        button_close.clicked.connect( self.close)
        hbox2.addWidget(button_close)


        vbox.addLayout(hbox2)


        self.setLayout(vbox)


        self.CalcR()
        self.CalcG()
        self.CalcB()
        self.draw_image()


#----------------------------------------------------------------------
    def OnSelectR(self, value):
        item = value
        self.r_spec = item

        self.CalcR()
        self.draw_image()

#----------------------------------------------------------------------
    def CalcR(self):

        if self.parent.tab_spec.showraw == True:
            tsmap = self.anlz.target_svd_maps[:,:,self.r_spec].copy()
        else:
            tsmap = self.anlz.target_pcafit_maps[:,:,self.r_spec].copy()


        uscale_min = tsmap.min()
        uscale_max = tsmap.max()

        scale_min = uscale_min + (uscale_max-uscale_min)*float(self.minr)/100.
        scale_max = uscale_min + (uscale_max-uscale_min)*float(self.maxr)/100.


        if scale_min >= scale_max:
            tsmap = np.zeros((self.n_cols, self.n_rows), dtype=float)
        else:
            tsmap = tsmap.clip(min=scale_min, max=scale_max)
            tsmap = (tsmap -scale_min) / (scale_max - scale_min)

        indices = np.where(tsmap < 0)
        tsmap[indices] = 0.0
        indices = np.where(tsmap > 1)
        tsmap[indices] = 1.0

        self.rgbimage[:,:,0] = tsmap*float(self.weightr)/100.

#----------------------------------------------------------------------
    def OnLimitMinR(self, value):

        self.minr = value
        #print 'self.minr=', self.minr
        self.CalcR()
        self.draw_image()

#----------------------------------------------------------------------
    def OnLimitMaxR(self, value):

        self.maxr = value
        #print 'self.maxr=', self.maxr
        self.CalcR()
        self.draw_image()


#----------------------------------------------------------------------
    def OnWeightR(self, value):

        self.weightr = value
        #print 'self.weightr=', self.weightr
        self.CalcR()
        self.draw_image()

#----------------------------------------------------------------------
    def OnSelectG(self, value):
        item = value
        self.g_spec = item

        self.CalcG()
        self.draw_image()

#----------------------------------------------------------------------
    def CalcG(self):

        if self.parent.tab_spec.showraw == True:
            tsmap = self.anlz.target_svd_maps[:,:,self.g_spec].copy()
        else:
            tsmap = self.anlz.target_pcafit_maps[:,:,self.g_spec].copy()


        uscale_min = tsmap.min()
        uscale_max = tsmap.max()

        scale_min = uscale_min + (uscale_max-uscale_min)*float(self.ming)/100.
        scale_max = uscale_min + (uscale_max-uscale_min)*float(self.maxg)/100.

        if scale_min >= scale_max:
            tsmap = np.zeros((self.n_cols, self.n_rows), dtype=float)
        else:
            tsmap = tsmap.clip(min=scale_min, max=scale_max)
            tsmap = (tsmap - scale_min) / (scale_max - scale_min)


        indices = np.where(tsmap < 0)
        tsmap[indices] = 0.0
        indices = np.where(tsmap > 1)
        tsmap[indices] = 1.0

        self.rgbimage[:,:,1] = tsmap*float(self.weightg)/100.

#----------------------------------------------------------------------
    def OnLimitMinG(self, value):

        self.ming = value
        #print 'self.ming=', self.ming
        self.CalcG()
        self.draw_image()

#----------------------------------------------------------------------
    def OnLimitMaxG(self, value):

        self.maxg = value
        #print 'self.maxg=', self.maxg
        self.CalcG()
        self.draw_image()

#----------------------------------------------------------------------
    def OnWeightG(self, value):

        self.weightg = value
        #print 'self.weightg=', self.weightg
        self.CalcG()
        self.draw_image()

#----------------------------------------------------------------------
    def OnSelectB(self, value):
        item = value
        self.b_spec = item

        self.CalcB()
        self.draw_image()

#----------------------------------------------------------------------
    def CalcB(self):

        if self.parent.tab_spec.showraw == True:
            tsmap = self.anlz.target_svd_maps[:,:,self.b_spec].copy()
        else:
            tsmap = self.anlz.target_pcafit_maps[:,:,self.b_spec].copy()

        uscale_min = tsmap.min()
        uscale_max = tsmap.max()

        scale_min = uscale_min + (uscale_max-uscale_min)*float(self.minb)/100.
        scale_max = uscale_min + (uscale_max-uscale_min)*float(self.maxb)/100.

        if scale_min >= scale_max:
            tsmap = np.zeros((self.n_cols, self.n_rows), dtype=float)
        else:
            tsmap = tsmap.clip(min=scale_min, max=scale_max)
            tsmap = (tsmap - scale_min) / (scale_max - scale_min)

        indices = np.where(tsmap < 0)
        tsmap[indices] = 0.0
        indices = np.where(tsmap > 1)
        tsmap[indices] = 1.0

        self.rgbimage[:,:,2] = tsmap*float(self.weightb)/100.

#----------------------------------------------------------------------
    def OnLimitMinB(self, value):

        self.minb = value
        #print 'self.minb=', self.minb
        self.CalcB()
        self.draw_image()

#----------------------------------------------------------------------
    def OnLimitMaxB(self, value):

        self.maxb = value
        #print 'self.maxb=', self.maxb
        self.CalcB()
        self.draw_image()

#----------------------------------------------------------------------
    def OnWeightB(self, value):

        self.weightb = value
        #print 'self.weightb=', self.weightb
        self.CalcB()
        self.draw_image()

#----------------------------------------------------------------------
    def OnShowInfo(self, state):

        if state == QtCore.Qt.Checked:
            self.show_info = 1
        else:
            self.show_info = 0

        self.draw_image()

#----------------------------------------------------------------------
    def draw_image(self):


        fig = self.RGBImagefig
        fig.clf()
        fig.add_axes(((0.0,0.0,1.0,1.0)))


        axes = fig.gca()
        fig.patch.set_alpha(1.0)

        im = axes.imshow(np.rot90(self.rgbimage))

        axes.axis("off")

        if self.show_info == 1:
            startx = int(self.n_rows*0.02)
            starty = self.n_cols-int(self.n_cols*0.15)
            info = 'R:%s [%d] \nG:%s [%d] \nB:%s [%d]' % (self.anlz.tspec_names[self.r_spec], self.weightr,
                                                        self.anlz.tspec_names[self.g_spec], self.weightg,
                                                        self.anlz.tspec_names[self.b_spec], self.weightb)
            axes.text(+startx+1,starty+1, info, horizontalalignment='left', verticalalignment='center',
                      color = 'white', fontsize=8)


        self.RGBImagePanel.draw()

#----------------------------------------------------------------------
    def OnSave(self, evt):

        wildcard = "Portable Network Graphics (*.png);;Adobe PDF Files (*.pdf);;"

        SaveFileName, _filter = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Plot', '', wildcard)

        SaveFileName = str(SaveFileName)
        if SaveFileName == '':
            return

        path, ext = os.path.splitext(SaveFileName)
        ext = ext[1:].lower()


        if ext != 'png' and ext != 'pdf':
            error_message = (
                  'Only the PNG and PDF image formats are supported.\n'
                 'A file extension of `png\' or `pdf\' must be used.')

            QtWidgets.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % error_message)
            return


        import matplotlib
        matplotlib.rcParams['pdf.fonttype'] = 42


        fig = self.RGBImagefig
        fig.savefig(SaveFileName, pad_inches = 0.0)
