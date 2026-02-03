import os
import copy
import time
import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtCore import Qt, QCoreApplication, pyqtSignal, pyqtSlot
from queue import SimpleQueue, Empty
import threading

from PIL import Image
from scipy import ndimage
from scipy.stats import linregress
from scipy.interpolate import interp1d
from skimage.registration import phase_cross_correlation
from skimage import filters
from importlib.metadata import version as version_check
from packaging.version import parse as parse_version

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib
import pyqtgraph as pg

class ImageRegistrationDialog(QtWidgets.QDialog):

    def __init__(self, parent, common):
        QtWidgets.QWidget.__init__(self, parent)
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dialogalign.ui'), self)
        self.parent = parent
        self.com = common
        # Currently disabled for Tomo data
        # if self.com.stack_4d == 1:
        #     self.bt_align2.setEnabled(True)
        # else:
        #     self.bt_align2.setEnabled(True)
        self.bt_align.clicked.connect(self.parent.page1.OnAlignImgs)
        self.bt_align2.clicked.connect(self.parent.page1.OnAlignImgs2)
        self.bt_align.clicked.connect(self.done)
        self.bt_align2.clicked.connect(self.done)

#----------------------------------------------------------------------
class ImageRegistrationManual(QtWidgets.QDialog):

    def __init__(self, parent,  common, stack):
        QtWidgets.QWidget.__init__(self, parent)

        self.parent = parent

        self.resize(1050, 750)
        self.setWindowTitle('Stack Alignment')

        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)

        self.stack = stack
        self.com = common

        self.have_ref_image = 0
        self.regist_calculated = 0

        self.iev = 0
        self.itheta = 0
        self.ref_image_index = 0
        self.ref_image_index_theta = 0
        self.ref_image = 0

        self.man_align = 0
        self.man_xref = 0
        self.man_yref = 0

        self.maxshift = 0
        self.auto = True

        self.xleft = 0
        self.xright = self.stack.n_cols
        self.ybottom = 0
        self.ytop = self.stack.n_rows

        if self.com.stack_4d == 0:
            self.aligned_stack = self.stack.absdata.copy()
            self.man_xs = np.zeros((self.stack.n_ev))
            self.man_ys = np.zeros((self.stack.n_ev))
            self.xshifts = np.zeros((self.stack.n_ev))
            self.yshifts = np.zeros((self.stack.n_ev))
        else:
            self.aligned_stack = self.stack.stack4D.copy()
            self.man_xs = np.zeros((self.stack.n_ev,self.stack.n_theta))
            self.man_ys = np.zeros((self.stack.n_ev,self.stack.n_theta))
            self.xshifts = np.zeros((self.stack.n_ev,self.stack.n_theta))
            self.yshifts = np.zeros((self.stack.n_ev,self.stack.n_theta))

        self.minxs = 0
        self.maxxs = 0

        self.minys = 0
        self.maxys = 0

        self.showccorr = 0

        self.edgee = 0

        self.subregion = 0
        self.sr_x1 = 0
        self.sr_x2 = 0
        self.sr_y1 = 0
        self.sr_y2 = 0
        self.patch = None

        #panel 1
        vbox1 = QtWidgets.QVBoxLayout()

        self.tc_imageeng = QtWidgets.QLabel(self)
        self.tc_imageeng.setText("Image at: ")

        gridsizertop = QtWidgets.QGridLayout()

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.absimgfig = Figure((4.0, 4.0))
        self.AbsImagePanel = FigureCanvas(self.absimgfig)
        self.AbsImagePanel.setParent(self)
        self.AbsImagePanel.mpl_connect('button_press_event', self.OnPointCorrimage)
        self.AbsImagePanel.setCursor(Qt.CrossCursor)


        fbox.addWidget(self.AbsImagePanel)
        frame.setLayout(fbox)
        gridsizertop.addWidget(frame, 0, 0, QtCore .Qt. AlignLeft)


        self.slider_eng = QtWidgets.QScrollBar(QtCore.Qt.Vertical)
        self.slider_eng.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_eng.valueChanged[int].connect(self.OnScrollEng)
        self.slider_eng.setRange(0, self.stack.n_ev-1)
        self.slider_eng.setValue(self.iev)

        gridsizertop.addWidget(self.slider_eng, 0, 1, QtCore .Qt. AlignLeft)

        self.slider_theta = QtWidgets.QScrollBar(QtCore.Qt.Horizontal)
        self.slider_theta.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_theta.valueChanged[int].connect(self.OnScrollTheta)
        self.slider_theta.setRange(0, self.stack.n_theta-1)

        self.tc_imagetheta = QtWidgets.QLabel(self)
        self.tc_imagetheta.setText("4D Data Angle: ")
        if self.com.stack_4d == 0 :
            self.tc_imagetheta.setVisible(False)
            self.slider_theta.setVisible(False)
        hbox51 = QtWidgets.QHBoxLayout()
        hbox51.addWidget(self.tc_imagetheta)
        hbox51.addWidget(self.slider_theta)
        gridsizertop.addLayout(hbox51, 1, 0)

        vbox1.addWidget(self.tc_imageeng)
        vbox1.addLayout(gridsizertop)

        self.tc_shift1 = QtWidgets.QLabel(self)
        self.tc_shift2 = QtWidgets.QLabel(self)
        vbox1.addWidget(self.tc_shift1)
        vbox1.addWidget(self.tc_shift2)
        vbox1.addStretch(1)

        if self.com.stack_4d == 0:
            self.tc_shift1.setText('X shift: {0:5.2f} pixels'.format(self.xshifts[self.iev]))
            self.tc_shift2.setText('Y shift: {0:5.2f} pixels'.format(self.yshifts[self.iev]))
        else:
            self.tc_shift1.setText('X shift: {0:5.2f} pixels'.format(self.xshifts[self.iev,self.itheta]))
            self.tc_shift2.setText('Y shift: {0:5.2f} pixels'.format(self.yshifts[self.iev,self.itheta]))


        #panel 2
        vbox2 = QtWidgets.QVBoxLayout()

        tc2 = QtWidgets.QLabel(self)
        tc2.setText('Cross-correlation')

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.cscorrfig = Figure((2.4, 2.4))

        self.CscorrPanel = FigureCanvas(self.cscorrfig)
        self.CscorrPanel.setParent(self)

        fbox.addWidget(self.CscorrPanel)
        frame.setLayout(fbox)

        vbox2.addWidget(tc2)
        vbox2.addWidget(frame)


        #panel 3
        vbox3 = QtWidgets.QVBoxLayout()

        tc3= QtWidgets.QLabel(self)
        tc3.setText('Image shifts')

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.shiftsfig = Figure((4.0, 2.4))
        self.ShiftsPanel = FigureCanvas(self.shiftsfig)
        self.ShiftsPanel.setParent(self)
        self.ShiftsPanel.mpl_connect('button_press_event', self.OnPlotShifts)

        fbox.addWidget(self.ShiftsPanel)
        frame.setLayout(fbox)

        vbox3.addWidget(tc3)
        vbox3.addWidget(frame)


        #panel 9
        vbox9 = QtWidgets.QVBoxLayout()
        vbox9.setSpacing(0)

        groupBox9 = QtWidgets.QGroupBox()
        self.rb_auto = QtWidgets.QRadioButton( 'Automatic Alignment')
        self.rb_man = QtWidgets.QRadioButton('Manual Alignment')
        self.rb_auto.setChecked(True)
        self.rb_auto.toggled.connect(self.Onrb_automanual)

        vbox9.addWidget(self.rb_auto)
        vbox9.addWidget(self.rb_man)
        groupBox9.setLayout(vbox9)

        #panel 8
        sizer8 = QtWidgets.QGroupBox('This Image')
        vbox8 = QtWidgets.QVBoxLayout()
        vbox8.setSpacing(0)

        self.button_refimg = QtWidgets.QPushButton('Set as Reference Image')
        self.button_refimg.clicked.connect(self.SetRefImage)
        vbox8.addWidget(self.button_refimg)

        self.button_refimgsave = QtWidgets.QPushButton('Save Reference Image')
        self.button_refimgsave.clicked.connect(self.SaveRefImage)
        self.button_refimgsave.setEnabled(False)
        vbox8.addWidget(self.button_refimgsave)

        self.button_refimgsload = QtWidgets.QPushButton('Load Reference Image')
        self.button_refimgsload.clicked.connect(self.LoadRefImage)
        vbox8.addWidget(self.button_refimgsload)

        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.HLine)
        line.setFrameShadow(QtWidgets.QFrame.Sunken)

        vbox8.addSpacing(5)
        vbox8.addWidget(line)
        vbox8.addSpacing(5)

        self.button_remove = QtWidgets.QPushButton('Remove energy from stack')
        self.button_remove.clicked.connect(self.OnRemoveImage)
        vbox8.addWidget(self.button_remove)

        sizer8.setLayout(vbox8)


        #panel 4
        sizer4 = QtWidgets.QGroupBox('Automatic Alignment')
        vbox4 = QtWidgets.QVBoxLayout()


        self.button_register = QtWidgets.QPushButton('Calculate image shifts')
        self.button_register.clicked.connect(self.OnCalcRegistration)
        self.button_register.setEnabled(False)
        vbox4.addWidget(self.button_register)
        #vbox4.addStretch(1)

        self.button_subregion = QtWidgets.QPushButton('Select subregion on reference')
        self.button_subregion.clicked.connect(self.OnSelectSubregion)
        self.button_subregion.setEnabled(False)
        vbox4.addWidget(self.button_subregion)

        self.button_delsubregion = QtWidgets.QPushButton('Remove subregion selection')
        self.button_delsubregion.clicked.connect(self.OnDeleteSubregion)
        self.button_delsubregion.setEnabled(False)
        vbox4.addWidget(self.button_delsubregion)
        vbox4.addStretch(1)

        groupBox4 = QtWidgets.QGroupBox()
        hbox43 = QtWidgets.QHBoxLayout()
        self.cb_edgeenh = QtWidgets.QCheckBox('Edge Enhancement', self)
        self.cb_edgeenh.stateChanged.connect(self.OnEdgeE)
        hbox43.addWidget(self.cb_edgeenh)

        self.rb_sobel = QtWidgets.QRadioButton( 'Sobel')
        self.rb_prewitt = QtWidgets.QRadioButton('Prewitt')
        self.rb_prewitt.setChecked(True)
        self.rb_sobel.setEnabled(False)
        self.rb_prewitt.setEnabled(False)

        hbox43.addWidget(self.rb_prewitt)
        hbox43.addWidget(self.rb_sobel)
        groupBox4.setLayout(hbox43)

        vbox4.addWidget(groupBox4)


        hbox42 = QtWidgets.QHBoxLayout()
        text1 = QtWidgets.QLabel(self)
        text1.setText(' Max shift [pixels]: ')

        self.tc_maxshift = QtWidgets.QSpinBox()
        self.tc_maxshift.setMinimum(0)
        self.tc_maxshift.valueChanged[int].connect(self.OnSetMaxShift)


        hbox42.addWidget(text1)
        hbox42.addWidget(self.tc_maxshift)
        vbox4.addLayout(hbox42)
        vbox4.addStretch(1)


        self.showcscor_cb = QtWidgets.QCheckBox('Show Cross-correlation', self)
        self.showcscor_cb.stateChanged.connect(self.OnShowCCorr)
        vbox4.addWidget(self.showcscor_cb)

        vbox4.addStretch()
        sizer4.setLayout(vbox4)


        #panel 5
        vbox5 = QtWidgets.QVBoxLayout()

        self.tc_refimg = QtWidgets.QLabel(self)
        self.tc_refimg.setText('Reference image')

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.refimgfig = Figure((4.0, 4.0))
        self.RefImagePanel = FigureCanvas(self.refimgfig)
        self.RefImagePanel.setParent(self)
        self.RefImagePanel.setCursor(Qt.CrossCursor)

        fbox.addWidget(self.RefImagePanel)
        frame.setLayout(fbox)

        self.RefImagePanel.mpl_connect('button_press_event', self.OnPointRefimage)
        self.RefImagePanel.mpl_connect('button_release_event', self.OnSelection)

        vbox5.addWidget(self.tc_refimg)
        vbox5.addWidget(frame)
        vbox5.addStretch(1)


        #panel 6
        sizer6 = QtWidgets.QGroupBox('Manual Alignment')
        vbox6 = QtWidgets.QVBoxLayout()
        vbox6.setSpacing(0)


        self.button_manalign = QtWidgets.QPushButton('Pick a point on reference image')
        self.button_manalign.clicked.connect(self.OnPickRefPoint)
        self.button_manalign.setEnabled(False)
        vbox6.addWidget(self.button_manalign)

        self.button_pick2ndpoint = QtWidgets.QPushButton('This image: click on same point')
        self.button_pick2ndpoint.clicked.connect(self.OnPickCorrPoint)
        self.button_pick2ndpoint.setEnabled(False)
        vbox6.addWidget(self.button_pick2ndpoint)

        self.button_applyman = QtWidgets.QPushButton('Apply manual shifts')
        self.button_applyman.clicked.connect(self.OnApplyManShifts)
        self.button_applyman.setEnabled(False)
        vbox6.addWidget(self.button_applyman)


        self.textctrl_ms1 = QtWidgets.QLabel(self)
        self.textctrl_ms2 = QtWidgets.QLabel(self)
        vbox6.addWidget(self.textctrl_ms1)
        vbox6.addWidget(self.textctrl_ms2)


        self.textctrl_ms1.setText('X manual shift: ')
        self.textctrl_ms2.setText('Y manual shift: ')


        sizer6.setLayout(vbox6)



        #panel 7
        vbox7 = QtWidgets.QVBoxLayout()
        vbox7.setSpacing(0)

        self.button_saveimg = QtWidgets.QPushButton('Save image shifts plot')
        self.button_saveimg.clicked.connect(self.OnSaveShiftsPlot)
        self.button_saveimg.setEnabled(False)
        vbox7.addWidget(self.button_saveimg)



        self.button_saveshifts = QtWidgets.QPushButton('Save image shifts')
        self.button_saveshifts.clicked.connect(self.OnSaveShifts)
        self.button_saveshifts.setEnabled(False)
        vbox7.addWidget(self.button_saveshifts)

        self.button_loadshifts = QtWidgets.QPushButton('Load image shifts')
        self.button_loadshifts.clicked.connect(self.OnLoadShifts)
        vbox7.addWidget(self.button_loadshifts)

        self.button_crop = QtWidgets.QPushButton('Crop aligned images')
        self.button_crop.clicked.connect(self.OnCropShifts)
        self.button_crop.setEnabled(False)
        vbox7.addWidget(self.button_crop)

        self.button_accept = QtWidgets.QPushButton('Apply Alignment')
        self.button_accept.clicked.connect(self.OnAccept)
        self.button_accept.setEnabled(False)
        vbox7.addWidget(self.button_accept)

        self.button_close = QtWidgets.QPushButton('Dismiss')
        self.button_close.clicked.connect(self.close)
        vbox7.addWidget(self.button_close)





        hboxtop = QtWidgets.QHBoxLayout()

        vboxL = QtWidgets.QVBoxLayout()
        vboxR = QtWidgets.QVBoxLayout()


        vboxL.addWidget(sizer8)
        vboxL.addStretch(1)
        vboxL.addWidget(groupBox9)
        vboxL.addStretch(1)
        vboxL.addWidget(sizer4)
        vboxL.addStretch(1)
        vboxL.addWidget(sizer6)
        vboxL.addStretch(1)
        vboxL.addLayout(vbox7)

        hboxRT = QtWidgets.QHBoxLayout()
        hboxRB = QtWidgets.QHBoxLayout()

        hboxRT.addLayout(vbox1)
        hboxRT.addLayout(vbox5)

        hboxRB.addLayout(vbox3)
        hboxRB.addLayout(vbox2)
        hboxRB.addStretch(1)

        vboxR.addStretch(1)
        vboxR.addLayout(hboxRT)
        vboxR.addStretch(5)
        vboxR.addLayout(hboxRB)
        vboxR.addStretch(5)

        hboxtop.addLayout(vboxL)
        hboxtop.addStretch(1)
        hboxtop.addLayout(vboxR)


        hboxtop.setContentsMargins(20,20,20,20)
        self.setLayout(hboxtop)


        self.ShowImage()



#----------------------------------------------------------------------
    def ShowImage(self):

        if self.com.stack_4d == 0:
            image = self.aligned_stack[:,:,self.iev]
        else:
            image = self.aligned_stack[:,:,self.iev,self.itheta]


        fig = self.absimgfig
        fig.clf()
        fig.add_axes(((0.0,0.0,1.0,1.0)))
        axes = fig.gca()


        _im = axes.imshow(np.rot90(image), cmap=matplotlib.colormaps["gray"])
        axes.autoscale(False)
#         if self.man_align == 2:
#             lx=matplotlib.lines.Line2D([self.man_yref-self.man_ys[self.iev],
#                                     self.man_yref-self.man_ys[self.iev]],
#                                     [0,self.stack.n_cols],color='red')
#             ly=matplotlib.lines.Line2D([0,self.stack.n_rows],
#                                    [self.man_xref-self.man_xs[self.iev],
#                                     self.man_xref-self.man_xs[self.iev]] ,color='red')
#             axes.add_line(lx)
#             axes.add_line(ly)



        axes.axis("off")
        self.AbsImagePanel.draw()


        self.tc_imageeng.setText('Image at: {0:5.2f} eV'.format(float(self.stack.ev[self.iev])))

        if self.com.stack_4d == 0:
            self.tc_shift1.setText('X shift: {0:5.2f} pixels'.format(self.xshifts[self.iev]))
            self.tc_shift2.setText('Y shift: {0:5.2f} pixels'.format(self.yshifts[self.iev]))
        else:
            self.tc_shift1.setText('X shift: {0:5.2f} pixels'.format(self.xshifts[self.iev,self.itheta]))
            self.tc_shift2.setText('Y shift: {0:5.2f} pixels'.format(self.yshifts[self.iev,self.itheta]))


        if (self.man_align == 2):
            if self.com.stack_4d == 0:
                self.textctrl_ms1.setText('X manual shift:  {0:5.2f}  pixels'.format(self.man_xs[self.iev]))
                self.textctrl_ms2.setText('Y manual shift:  {0:5.2f}  pixels'.format(self.man_ys[self.iev]))
            else:
                self.textctrl_ms1.setText('X manual shift:  {0:5.2f}  pixels'.format(self.man_xs[self.iev,self.itheta]))
                self.textctrl_ms2.setText('Y manual shift:  {0:5.2f}  pixels'.format(self.man_ys[self.iev,self.itheta]))



#----------------------------------------------------------------------
    def OnScrollEng(self, value):
        self.iev = value

        self.ShowImage()

#----------------------------------------------------------------------
    def OnScrollTheta(self, value):
        self.itheta = value

        self.tc_imagetheta.setText("4D Data Angle: "+str(self.stack.theta[self.itheta]))

        self.ShowImage()

#----------------------------------------------------------------------
    def SetRefImage(self):

        self.ref_image_index = self.iev
        self.ref_image_index_theta = self.itheta

        if self.com.stack_4d == 0:
            self.ref_image = self.aligned_stack[:,:,self.iev].copy()
        else:
            self.ref_image = self.aligned_stack[:,:,self.iev,self.itheta].copy()


        self.ShowRefImage()
        self.have_ref_image = 1

        self.UpdateWidgets()


#----------------------------------------------------------------------
    def SaveRefImage(self):

        wildcard = "TIFF File (*.tif);;"

        fileName, _filter = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Reference Image', '', wildcard)

        fileName = str(fileName)
        if fileName == '':
            return


        img1 = Image.fromarray(self.ref_image)
        img1.save(fileName)


#----------------------------------------------------------------------
    def LoadRefImage(self):

        wildcard = "TIFF File (*.tif);;"

        fileName, _filter = QtWidgets.QFileDialog.getOpenFileName(self, 'Save Reference Image', '', wildcard)

        fileName = str(fileName)
        if fileName == '':
            return

        from PIL import Image
        img = Image.open(fileName)

        self.ref_image = np.array((img))

        self.ShowRefImage()
        self.have_ref_image = 1

        self.UpdateWidgets()

#----------------------------------------------------------------------
    def ShowRefImage(self):

        fig = self.refimgfig
        fig.clf()
        fig.add_axes(((0.0,0.0,1.0,1.0)))
        axes = fig.gca()


        _im = axes.imshow(np.rot90(self.ref_image), cmap=matplotlib.colormaps["gray"])

        if (self.subregion == 1):

            from matplotlib.path import Path
            import matplotlib.patches as patches

            verts = [
                     (self.sr_x1, self.sr_y1), # left, bottom
                     (self.sr_x1, self.sr_y2), # left, top
                     (self.sr_x2, self.sr_y2), # right, top
                     (self.sr_x2, self.sr_y1), # right, bottom
                     (self.sr_x1, self.sr_y1), # ignored
                     ]

            codes = [Path.MOVETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.CLOSEPOLY,
                     ]

            path = Path(verts, codes)

            self.patch = patches.PathPatch(path, facecolor='red')
            self.patch.set_alpha(0.3)
            axes.add_patch(self.patch)

        if self.man_align > 0:
            lx=matplotlib.lines.Line2D([self.man_xref,self.man_xref], [0,self.stack.n_rows],color='green')
            ly=matplotlib.lines.Line2D([0,self.stack.n_cols], [self.man_yref,self.man_yref] ,color='green')
            axes.add_line(lx)
            axes.add_line(ly)

        axes.axis("off")

        self.RefImagePanel.draw()


#----------------------------------------------------------------------
    def Onrb_automanual(self, enabled):

        state = enabled


        if state:
            self.auto = True
            self.man_align = 0
        else:
            self.auto = False



        if self.com.stack_4d == 0:
            self.aligned_stack = self.stack.absdata.copy()
            self.man_xs = np.zeros((self.stack.n_ev))
            self.man_ys = np.zeros((self.stack.n_ev))
            self.xshifts = np.zeros((self.stack.n_ev))
            self.yshifts = np.zeros((self.stack.n_ev))
        else:
            self.aligned_stack = self.stack.stack4D.copy()
            self.man_xs = np.zeros((self.stack.n_ev,self.stack.n_theta))
            self.man_ys = np.zeros((self.stack.n_ev,self.stack.n_theta))
            self.xshifts = np.zeros((self.stack.n_ev,self.stack.n_theta))
            self.yshifts = np.zeros((self.stack.n_ev,self.stack.n_theta))



        if self.have_ref_image == 1:
            fig = self.shiftsfig
            fig.clf()
            self.ShiftsPanel.draw()

            fig = self.cscorrfig
            fig.clf()
            self.CscorrPanel.draw()

            self.ShowImage()
            self.ShowRefImage()

        self.UpdateWidgets()

#----------------------------------------------------------------------
    def ShowCrossCorrelation(self, ccorr, xshift, yshift):

        fig = self.cscorrfig

        fig.clf()
        fig.add_axes(((0.0,0.0,1.0,1.0)))
        axes = fig.gca()

        im = axes.imshow(np.rot90(ccorr), cmap=matplotlib.colormaps["gray"])

        nx = ccorr.shape[0]
        ny = ccorr.shape[1]
        xcenter = xshift + float(nx)/2.0
        ycenter = yshift + float(ny)/2.0

        xl = xcenter-10
        if xl<0:
            xl=0
        xr = xcenter+10
        if xr>nx-1:
            xr=nx-1
        yl = ycenter-10
        if yl<0:
            yl=0
        yr = ycenter+10
        if yr>ny-1:
            yr=ny-1

        lx=matplotlib.lines.Line2D([xl,xr], [ycenter,ycenter], color='green')
        ly=matplotlib.lines.Line2D([xcenter,xcenter], [yl,yr], color='green')
        axes.add_line(lx)
        axes.add_line(ly)


        axes.axis("off")

        self.CscorrPanel.draw()

#----------------------------------------------------------------------
    def OnShowCCorr(self, state):

        if state == QtCore.Qt.Checked:
            self.showccorr = 1
        else:
            self.showccorr = 0


#----------------------------------------------------------------------
    def OnEdgeE(self, state):

        if state == QtCore.Qt.Checked:
            self.edgee = 1
            self.rb_sobel.setEnabled(True)
            self.rb_prewitt.setEnabled(True)
        else:
            self.edgee = 0
            self.rb_sobel.setEnabled(False)
            self.rb_prewitt.setEnabled(False)

#----------------------------------------------------------------------
    def OnSetMaxShift(self, value):

        self.maxshift = value

#----------------------------------------------------------------------
    def OnRemoveImage(self, event):


        self.stack.absdata = np.delete(self.stack.absdata, self.iev, axis=2)

        self.aligned_stack = np.delete(self.aligned_stack, self.iev, axis=2)

        if self.com.stack_4d == 1:
            self.stack.stack4D = np.delete(self.stack.stack4D, self.iev, axis=2)


        self.stack.n_ev = self.stack.n_ev - 1
        self.stack.ev = np.delete(self.stack.ev, self.iev)




        if self.com.stack_4d == 0:
            self.stack.data_struct.exchange.data = self.stack.absdata
            self.xshifts = np.delete(self.xshifts, self.iev)
            self.yshifts = np.delete(self.yshifts, self.iev)
        else:
            self.stack.data_struct.exchange.data = self.stack.stack4D
            self.xshifts = np.delete(self.xshifts, self.iev, axis=0)
            self.yshifts = np.delete(self.yshifts, self.iev, axis=0)

        self.stack.data_struct.exchange.energy = self.stack.ev

        if self.com.i0_loaded == 1:
            self.stack.calculate_optical_density()

        self.iev = self.iev-1
        if self.iev < 0:
            self.iev = 0


        #self.slider_eng.setRange(0, self.stack.n_ev-1)

        #self.parent.page1.slider_eng.setRange(0, self.stack.n_ev-1)
        #self.parent.page1.iev = self.stack.n_ev/2
        #self.parent.page1.slider_eng.setValue(self.parent.page1.iev)

        #self.parent.page1.loadSpectrum(self.parent.page1.ix, self.parent.page1.iy)
        #self.parent.page1.loadImage()
        self.parent.page1.absimgfig.loadNewImageWithROI()
        self.parent.page0.absimgfig.loadNewImage()
        self.parent.page1.specfig.ClearandReload()

        self.ShowImage()

#----------------------------------------------------------------------
    def OnCalcRegistration(self, event):

        if self.com.stack_4d == 1:
            self.CalcRegistration4D()
            return

        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))

        #Get Edge enhancement info
        edge = 0
        if self.edgee > 0:
            if self.rb_sobel.isChecked():
                edge = 1
            else:
                edge = 2

        #Subregion selection on a reference image
        if self.subregion == 0:
            self.sr_x1 = 0
            self.sr_x2 = 0
            self.sr_y1 = 0
            self.sr_y2 = 0
            referenceimage = self.ref_image
        else:
            referenceimage = self.ref_image[self.sr_x1:self.sr_x2, self.sr_y2:self.sr_y1]

        for i in range(self.stack.n_ev):

            if self.subregion == 0:
                img2 = self.aligned_stack[:,:,i]
            else:
                img2 = self.aligned_stack[self.sr_x1:self.sr_x2, self.sr_y2:self.sr_y1, i]

            if i==0:
                xshift, yshift, ccorr = self.stack.register_images(referenceimage, img2,
                                                          have_ref_img_fft = False,
                                                          edge_enhancement = edge)
            elif i==self.ref_image_index:
                xshift = 0
                yshift = 0
            else:
                xshift, yshift, ccorr = self.stack.register_images(referenceimage, img2,
                                                          have_ref_img_fft = True,
                                                          edge_enhancement =  edge)

            #Limit the shifts to MAXSHIFT chosen by the user
            if (self.maxshift > 0):
                if (abs(xshift) > self.maxshift):
                        xshift = np.sign(xshift)*self.maxshift
                if (abs(yshift) > self.maxshift):
                        yshift = np.sign(yshift)*self.maxshift

            self.xshifts[i] = xshift
            self.yshifts[i] = yshift
            self.PlotShifts()

            if self.showccorr == 1:
                self.ShowCrossCorrelation(ccorr, xshift, yshift)

            QCoreApplication.processEvents()


        #Apply shifts
        for i in range(self.stack.n_ev):
            img = self.aligned_stack[:,:,i]
            if (abs(self.xshifts[i])>0.02) or (abs(self.yshifts[i])>0.02):
                shifted_img = self.stack.apply_image_registration(img,
                                                                  self.xshifts[i],
                                                                  self.yshifts[i])
                self.aligned_stack[:,:,i] = shifted_img


        self.regist_calculated = 1
        self.iev = 0
        self.ShowImage()
        self.slider_eng.setValue(self.iev)

        self.UpdateWidgets()


        min_xshift = np.min(self.xshifts)
        max_xshift = np.max(self.xshifts)

        min_yshift = np.min(self.yshifts)
        max_yshift = np.max(self.yshifts)

        if min_xshift < self.minxs : self.minxs = min_xshift
        if max_xshift > self.maxxs : self.maxxs = max_xshift

        if min_yshift < self.minys : self.minys = min_yshift
        if max_yshift > self.maxys : self.maxys = max_yshift


        QtWidgets.QApplication.restoreOverrideCursor()

#----------------------------------------------------------------------
    def CalcRegistration4D(self):

        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))

        #Get Edge enhancement info
        edge = 0
        if self.edgee > 0:
            if self.rb_sobel.isChecked():
                edge = 1
            else:
                edge = 2

        #Subregion selection on a reference image
        if self.subregion == 0:
            self.sr_x1 = 0
            self.sr_x2 = 0
            self.sr_y1 = 0
            self.sr_y2 = 0
            referenceimage = self.ref_image
        else:
            referenceimage = self.ref_image[self.sr_x1:self.sr_x2, self.sr_y2:self.sr_y1]

        temptheta = self.itheta
        for j in range(self.stack.n_theta):
            self.itheta = j
            for i in range(self.stack.n_ev):


                if self.subregion == 0:
                    img2 = self.aligned_stack[:,:,i,j]
                else:
                    img2 = self.aligned_stack[self.sr_x1:self.sr_x2, self.sr_y2:self.sr_y1, i, j]

                if i==0 and j==0:
                    xshift, yshift, ccorr = self.stack.register_images(referenceimage, img2,
                                                              have_ref_img_fft = False,
                                                              edge_enhancement = edge)
                elif i==self.ref_image_index and j==self.ref_image_index_theta:
                    xshift = 0
                    yshift = 0
                else:
                    xshift, yshift, ccorr = self.stack.register_images(referenceimage, img2,
                                                              have_ref_img_fft = True,
                                                              edge_enhancement =  edge)

                #Limit the shifts to MAXSHIFT chosen by the user
                if (self.maxshift > 0):
                    if (abs(xshift) > self.maxshift):
                            xshift = np.sign(xshift)*self.maxshift
                    if (abs(yshift) > self.maxshift):
                            yshift = np.sign(yshift)*self.maxshift

                self.xshifts[i,j] = xshift
                self.yshifts[i,j] = yshift
                self.PlotShifts()

                if self.showccorr == 1:
                    self.ShowCrossCorrelation(ccorr, xshift, yshift)

                QCoreApplication.processEvents()


        #Apply shifts
        for i in range(self.stack.n_ev):
            for j in range(self.stack.n_theta):
                img = self.aligned_stack[:,:,i,j]
                if (abs(self.xshifts[i,j])>0.02) or (abs(self.yshifts[i,j])>0.02):
                    shifted_img = self.stack.apply_image_registration(img,
                                                                      self.xshifts[i,j],
                                                                      self.yshifts[i,j])
                    self.aligned_stack[:,:,i,j] = shifted_img


        self.itheta = temptheta
        self.regist_calculated = 1
        self.iev = 0
        self.ShowImage()
        self.slider_eng.setValue(self.iev)

        self.UpdateWidgets()


        min_xshift = np.min(self.xshifts)
        max_xshift = np.max(self.xshifts)

        min_yshift = np.min(self.yshifts)
        max_yshift = np.max(self.yshifts)

        if min_xshift < self.minxs : self.minxs = min_xshift
        if max_xshift > self.maxxs : self.maxxs = max_xshift

        if min_yshift < self.minys : self.minys = min_yshift
        if max_yshift > self.maxys : self.maxys = max_yshift


        QtWidgets.QApplication.restoreOverrideCursor()

#----------------------------------------------------------------------
    def OnCropShifts(self, event):

        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))

        self.aligned_stack, self.xleft, self.xright, self.ybottom, self.ytop, = self.stack.crop_registed_images(self.aligned_stack,
                                                             self.minxs, self.maxxs, self.minys, self.maxys)


        self.iev = 0
        self.ShowImage()
        self.slider_eng.setValue(self.iev)

        QtWidgets.QApplication.restoreOverrideCursor()


#----------------------------------------------------------------------
    def OnPickRefPoint(self, event):
        self.man_align = 1


#----------------------------------------------------------------------
    def OnPointRefimage(self, evt):


        x = evt.xdata
        y = evt.ydata

        if (x == None) or (y == None):
            return

        if self.subregion == 1:
            self.sr_x1 = x
            self.sr_y1 = y
            self.button_pressed = True
            self.mousemoveconn = self.RefImagePanel.mpl_connect('motion_notify_event', self.OnSelectionMotion)
            return

        if self.man_align == 0:
            pass

        if (self.man_align == 1):
            self.man_xref = int(np.floor(x))
            self.man_yref = int(np.floor(y))

            if self.man_xref<0 :
                self.man_xref=0
            if self.man_xref>self.stack.n_cols :
                self.man_xref=self.stack.n_cols
            if self.man_yref<0 :
                self.man_yref=0
            if self.man_yref>self.stack.n_rows :
                self.man_yref=self.stack.n_rows


            self.UpdateWidgets()

            self.ShowRefImage()


#----------------------------------------------------------------------
    def OnPickCorrPoint(self):

        self.man_align = 2

#----------------------------------------------------------------------
    def OnPointCorrimage(self, evt):

        x = evt.xdata
        y = evt.ydata

        if (self.man_align == 2):
            xcorr = float(x)
            ycorr = float(y)

            if xcorr<0 :
                xcorr=0
            if xcorr>self.stack.n_cols :
                xcorr=self.stack.n_cols
            if ycorr<0 :
                ycorr=0
            if ycorr>self.stack.n_rows :
                ycorr=self.stack.n_rows


            if self.com.stack_4d == 0:
                self.man_xs[self.iev] = self.man_xref - xcorr
                self.man_ys[self.iev] = -1.0*(self.man_yref - ycorr)


                self.textctrl_ms1.setText('X manual shift:  {0:5.2f}  pixels\n'.format(self.man_xs[self.iev]))
                self.textctrl_ms2.setText('Y manual shift:  {0:5.2f}  pixels'.format(self.man_ys[self.iev]))
            else:
                self.man_xs[self.iev, self.itheta] = self.man_xref - xcorr
                self.man_ys[self.iev, self.itheta] = -1.0*(self.man_yref - ycorr)


                self.textctrl_ms1.setText('X manual shift:  {0:5.2f}  pixels\n'.format(self.man_xs[self.iev, self.itheta]))
                self.textctrl_ms2.setText('Y manual shift:  {0:5.2f}  pixels'.format(self.man_ys[self.iev, self.itheta]))

            self.iev = self.iev + 1
            if self.iev > (self.stack.n_ev-1):
                self.iev = 0

            self.slider_eng.setValue(self.iev)


            self.ShowImage()
            self.UpdateWidgets()

#----------------------------------------------------------------------
    def OnApplyManShifts(self):


        if self.com.stack_4d == 0:
            for i in range(self.stack.n_ev):

                img = self.aligned_stack[:,:,i]
                if (abs(self.man_xs[i])>0.02) or (abs(self.man_ys[i])>0.02):
                    shifted_img = self.stack.apply_image_registration(img,
                                                                      self.man_xs[i],
                                                                      self.man_ys[i])
                    self.aligned_stack[:,:,i] = shifted_img

                    self.xshifts[i] = self.xshifts[i] + self.man_xs[i]
                    self.yshifts[i] = self.yshifts[i] + self.man_ys[i]




                self.man_xs[i] = 0
                self.man_ys[i] = 0

        else:
            for i in range(self.stack.n_ev):
                for j in range(self.stack.n_theta):

                    img = self.aligned_stack[:,:,i,j]
                    if (abs(self.man_xs[i,j])>0.02) or (abs(self.man_ys[i,j])>0.02):
                        shifted_img = self.stack.apply_image_registration(img,
                                                                          self.man_xs[i,j],
                                                                          self.man_ys[i,j])
                        self.aligned_stack[:,:,i,j] = shifted_img

                        self.xshifts[i,j] = self.xshifts[i,j] + self.man_xs[i,j]
                        self.yshifts[i,j] = self.yshifts[i,j] + self.man_ys[i,j]




                    self.man_xs[i,j] = 0
                    self.man_ys[i,j] = 0


        self.regist_calculated = 1
        self.man_align = 0

        self.ShowRefImage()
        self.PlotShifts()

        self.textctrl_ms1.setText('X manual shift: \n')
        self.textctrl_ms2.setText('Y manual shift: ')

        self.UpdateWidgets()

        min_xshift = np.min(self.xshifts)
        max_xshift = np.max(self.xshifts)

        min_yshift = np.min(self.yshifts)
        max_yshift = np.max(self.yshifts)

        if min_xshift < self.minxs : self.minxs = min_xshift
        if max_xshift > self.maxxs : self.maxxs = max_xshift

        if min_yshift < self.minys : self.minys = min_yshift
        if max_yshift > self.maxys : self.maxys = max_yshift


        self.ShowImage()




#----------------------------------------------------------------------
    def OnAccept(self):

        if self.com.stack_4d == 0:
            self.stack.absdata = self.aligned_stack
            self.stack.data_struct.exchange.data = self.stack.absdata
        else:
            self.stack.stack4D = self.aligned_stack
            self.stack.absdata = self.stack.stack4D[:,:,:,self.itheta]
            self.stack.data_struct.exchange.data = self.stack.stack4D


        datadim = np.int32(self.stack.absdata.shape)

        self.stack.n_cols = datadim[0].copy()
        self.stack.n_rows =  datadim[1].copy()

        self.stack.xshifts = self.xshifts
        self.stack.yshifts = self.yshifts

        if self.com.i0_loaded == 1:
            if self.com.stack_4d == 0:
                #Resize optical density
                for i in range(self.stack.n_ev):

                    img = self.stack.od3d[:,:,i]
                    shifted_img = self.stack.apply_image_registration(img, self.xshifts[i], self.yshifts[i])
                    self.stack.od3d[:,:,i] = shifted_img


                self.stack.od3d = self.stack.od3d[self.xleft:self.xright, self.ybottom:self.ytop, :]

                self.stack.od = self.stack.od3d.copy()
                self.stack.od = np.reshape(self.stack.od, (self.stack.n_cols*self.stack.n_rows, self.stack.n_ev), order='F')

            else:
                #Resize optical density for 4D stack

                for i in range(self.stack.n_ev):
                    for j in range(self.stack.n_theta):
                        img = self.stack.od4d[:,:,i,j]
                        shifted_img = self.stack.apply_image_registration(img, self.xshifts[i,j], self.yshifts[i,j])
                        self.stack.od4d[:,:,i,j] = shifted_img

                self.stack.od4d = self.stack.od4d[self.xleft:self.xright, self.ybottom:self.ytop, :, :]

                self.stack.od3d = self.stack.od4d[:,:,:,self.itheta]
                self.stack.od = self.stack.od3d.copy()
                n_pixels = self.stack.n_cols*self.stack.n_rows
                self.stack.od = np.reshape(self.stack.od, (n_pixels, self.stack.n_ev), order='F')

            self.stack.data_struct.spectromicroscopy.optical_density = self.stack.od


        self.stack.data_struct.exchange.energy = self.stack.ev



        self.stack.data_struct.spectromicroscopy.xshifts = self.xshifts
        self.stack.data_struct.spectromicroscopy.yshifts = self.yshifts


        #self.parent.page1.slider_eng.setRange(0,self.stack.n_ev-1)
        #self.parent.page1.iev = int(self.stack.n_ev/2)
        #self.parent.page1.slider_eng.setValue(self.parent.page1.iev)

        #self.parent.page1.ix = int(self.stack.n_cols/2)
        #self.parent.page1.iy = int(self.stack.n_rows/2)

        #self.parent.page1.loadSpectrum(self.parent.page1.ix, self.parent.page1.iy)
        #self.parent.page1.loadImage()
        #self.parent.page9.loadImage()
        self.parent.page1.absimgfig.loadNewImageWithROI()
        self.parent.page0.absimgfig.loadNewImage()
        self.parent.page1.specfig.ClearandReload()

        self.close()





#----------------------------------------------------------------------
    def PlotShifts(self, ):

        fig = self.shiftsfig
        fig.clf()

        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()



        #Matplotlib has inverted axes!
        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Shifts (x-red, y-green) [pixels]')

        if self.com.stack_4d == 0:
            plot = axes.plot(self.stack.ev,self.xshifts, color='green')
            plot = axes.plot(self.stack.ev,self.yshifts, color='red')
        else:
            plot = axes.plot(self.stack.ev,self.xshifts[:,self.itheta], color='green')
            plot = axes.plot(self.stack.ev,self.yshifts[:,self.itheta], color='red')

        self.ShiftsPanel.draw()

#----------------------------------------------------------------------
    def OnPlotShifts(self, evt):
        x = evt.xdata
        y = evt.ydata

        if self.xshifts.any():
            if x < self.stack.ev[0]:
                sel_ev = 0
            elif x > self.stack.ev[self.stack.n_ev-1]:
                sel_ev = self.stack.n_ev-1
            else:
                indx = np.abs(self.stack.ev - x).argmin()
                sel_ev = indx

            self.iev = sel_ev

            self.ShowImage()

            self.slider_eng.setValue(self.iev)

#----------------------------------------------------------------------
    def OnSaveShiftsPlot(self, evt):


        wildcard = "Portable Network Graphics (*.png);;Adobe PDF Files (*.pdf);;"

        self.SaveFileName, _filter = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Image Shifts Plot', '', wildcard)

        self.SaveFileName = str(self.SaveFileName)
        if self.SaveFileName == '':
            return


        path, ext = os.path.splitext(self.SaveFileName)
        ext = ext[1:].lower()

        try:

            matplotlib.rcParams['pdf.fonttype'] = 42
            if ext == 'png':

                fig = self.shiftsfig
                fig.savefig(self.SaveFileName)

            if ext =='pdf':


                fig = self.shiftsfig
                fig.savefig(self.SaveFileName)

        except:
            pass

#----------------------------------------------------------------------
    def OnSelectSubregion(self, event):

        self.subregion = 1
        self.sr_x1 = 0
        self.sr_x2 = 0
        self.sr_y1 = 0
        self.sr_y2 = 0

        self.button_delsubregion.setEnabled(True)


#----------------------------------------------------------------------
    def OnDeleteSubregion(self, event):

        self.subregion = 0
        self.sr_x1 = 0
        self.sr_x2 = 0
        self.sr_y1 = 0
        self.sr_y2 = 0

        self.button_pressed = False
        self.patch = None
        self.RefImagePanel.mpl_disconnect(self.mousemoveconn)

        self.ShowRefImage()


#----------------------------------------------------------------------
    def OnSelection(self, evt):

        if (self.man_align > 0) and (self.subregion == 0):
            return

        x2, y2 = evt.xdata, evt.ydata

        if (x2 == None) or (y2 == None):
            return

        self.sr_x1 = int(self.sr_x1)
        self.sr_x2 = int(x2)

        self.sr_y2 = int(self.sr_y1)
        self.sr_y1 = int(y2)

        if self.sr_x1 > self.sr_x2:
            temp = self.sr_x1
            self.sr_x1 = self.sr_x2
            self.sr_x2 = temp

        if self.sr_y2 > self.sr_y1:
            temp = self.sr_y1
            self.sr_y1 = self.sr_y2
            self.sr_y2 = temp


        self.button_pressed = False

        self.RefImagePanel.mpl_disconnect(self.OnSelectionMotion)

        self.ShowRefImage()

#----------------------------------------------------------------------
    def OnSelectionMotion(self, event):

        x2, y2 = event.xdata, event.ydata

        if (x2 == None) or (y2 == None):
            return

        if self.button_pressed == False:
            return

        self.sr_x2 = int(x2)
        self.sr_y2 = int(y2)

        fig = self.refimgfig

        axes = fig.gca()



        from matplotlib.path import Path
        import matplotlib.patches as patches

        verts = [
                 (self.sr_x1, self.sr_y1), # left, bottom
                 (self.sr_x1, self.sr_y2), # left, top
                 (self.sr_x2, self.sr_y2), # right, top
                 (self.sr_x2, self.sr_y1), # right, bottom
                 (self.sr_x1, self.sr_y1), # ignored
                 ]

        codes = [Path.MOVETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.CLOSEPOLY,
                 ]

        path = Path(verts, codes)

        if self.patch != None:
            self.patch.remove()

        self.patch = patches.PathPatch(path, facecolor='red')
        self.patch.set_alpha(0.3)
        axes.add_patch(self.patch)

        axes.axis("off")

        self.RefImagePanel.draw()

#----------------------------------------------------------------------
    def OnSaveShifts(self, evt):


        wildcard = "CSV files (*.csv)"

        filepath, _filter = QtWidgets.QFileDialog.getSaveFileName(self, 'Please select an alignment file (.csv)', '', wildcard)

        filepath = str(filepath)
        if filepath == '':
            return




        f = open(filepath, 'w')

        if self.com.stack_4d == 0:
            print('*********************  Alignment file  ********************', file=f)
            print('***  for ', self.com.filename, file=f)
            print('***  ev, xshift, yshift', file=f)
            for ie in range(self.stack.n_ev):
                print('%.6f, %.6f, %.6f' %(self.stack.ev[ie], self.xshifts[ie], self.yshifts[ie]), file=f)
        else:
            print('*********************  Alignment file  ********************', file=f)
            print('***  for ', self.com.filename, file=f)
            print('***  ev, theta, xshift, yshift', file=f)
            for i in range(self.stack.n_ev):
                for j in range(self.stack.n_theta):
                    print('%.6f, %.6f, %.6f, %.6f' %(self.stack.ev[i], self.stack.theta[j], self.xshifts[i,j], self.yshifts[i,j]), file=f)


        f.close()

#----------------------------------------------------------------------
    def OnLoadShifts(self, evt):



        wildcard = "I0 CSV files (*.csv)"

        filepath, _filter = QtWidgets.QFileDialog.getOpenFileName(self, 'Please select an alignment file (.csv)', '', wildcard)


        filepath = str(filepath)
        if filepath == '':
            return


        f = open(str(filepath),'r')

        elist = []
        xshiftlist = []
        yshiftlist = []

        for line in f:
            if line.startswith('*'):
                continue
            else:
                e, xs, ys = [float (x) for x in line.split(',')]
                elist.append(e)
                xshiftlist.append(xs)
                yshiftlist.append(ys)

        f.close()


        self.xshifts = np.zeros((self.stack.n_ev))
        self.yshifts = np.zeros((self.stack.n_ev))

        for ie in range(self.stack.n_ev):
            engfl = '%.6f' % ( self.stack.ev[ie])
            eng = float(engfl)
            if eng in elist:
                ind = elist.index(eng)
                self.xshifts[ie] = xshiftlist[ind]
                self.yshifts[ie] = yshiftlist[ind]


        #Apply shifts
        self.PlotShifts()
        for i in range(self.stack.n_ev):
            img = self.aligned_stack[:,:,i]
            if (abs(self.xshifts[i])>0.02) or (abs(self.yshifts[i])>0.02):
                shifted_img = self.stack.apply_image_registration(img,
                                                                  self.xshifts[i],
                                                                  self.yshifts[i])
                self.aligned_stack[:,:,i] = shifted_img


        self.regist_calculated = 1
        self.iev = 0
        self.ShowImage()
        self.slider_eng.setValue(self.iev)

        self.regist_calculated = 1
        self.man_align = 0

        self.textctrl_ms1.setText('X manual shift: \n')
        self.textctrl_ms2.setText('Y manual shift: ')

        self.UpdateWidgets()

        min_xshift = np.min(self.xshifts)
        max_xshift = np.max(self.xshifts)

        min_yshift = np.min(self.yshifts)
        max_yshift = np.max(self.yshifts)

        if min_xshift < self.minxs : self.minxs = min_xshift
        if max_xshift > self.maxxs : self.maxxs = max_xshift

        if min_yshift < self.minys : self.minys = min_yshift
        if max_yshift > self.maxys : self.maxys = max_yshift


#----------------------------------------------------------------------
    def UpdateWidgets(self):

        if self.auto:
            self.button_manalign.setEnabled(False)
            if self.have_ref_image == 1:
                self.button_register.setEnabled(True)
                self.button_subregion.setEnabled(True)
                self.button_refimgsave.setEnabled(True)
            else:
                self.button_register.setEnabled(False)
                self.button_subregion.setEnabled(False)
                self.button_refimgsave.setEnabled(False)
                self.button_delsubregion.setEnabled(False)

            if self.regist_calculated == 1:
                self.button_crop.setEnabled(True)
                self.button_accept.setEnabled(True)
                self.button_saveshifts.setEnabled(True)
                self.button_saveimg.setEnabled(True)
            else:
                self.button_crop.setEnabled(False)
                self.button_accept.setEnabled(False)
                self.button_saveshifts.setEnabled(False)
                self.button_saveimg.setEnabled(False)

        else:
            self.button_register.setEnabled(False)
            self.button_delsubregion.setEnabled(False)
            self.button_subregion.setEnabled(False)

            if self.have_ref_image == 1:
                self.button_manalign.setEnabled(True)
                self.button_refimgsave.setEnabled(True)
            else:
                self.button_manalign.setEnabled(False)
                self.button_refimgsave.setEnabled(False)

            if self.regist_calculated == 1:
                self.button_crop.setEnabled(True)
                self.button_accept.setEnabled(True)
                self.button_saveshifts.setEnabled(True)
            else:
                self.button_crop.setEnabled(False)
                self.button_accept.setEnabled(False)
                self.button_saveshifts.setEnabled(False)

            if self.man_align == 0:
                self.button_pick2ndpoint.setEnabled(False)
                self.button_applyman.setEnabled(False)
            elif self.man_align == 1:
                self.button_pick2ndpoint.setEnabled(True)
            elif self.man_align == 2:
                self.button_applyman.setEnabled(True)

class GeneralPurposeSignals(QtCore.QObject):
    finished = pyqtSignal()
    ithetaprogress = pyqtSignal(int)
# ----------------------------------------------------------------------
class GeneralPurposeProcessor(QtCore.QRunnable):
    def __init__(self, parent, queue):
        super(GeneralPurposeProcessor, self).__init__()
        self.signals = GeneralPurposeSignals()
        self.funcdict = {"ShiftImg": self.ShiftImg, "AlignReferenced": self.AlignReferenced}
        self.parent = parent
        self.queue = queue
        self.current_itheta = 0

    @pyqtSlot()
    def run(self):
        #print('worker', threading.get_ident())
        while True:
            #print(self.parent.pool.pool.activeThreadCount())
            #print(QtCore.QThreadPool.activeThreadCount())
            try:
                workerfunc, *args = self.queue.get(False)
                self.funcdict[workerfunc](*args[0])
                #print("busy with task {}".format(workerfunc))
            except Empty: # if queue empty
                break

    def AlignReferenced(self, data, itheta):
        if self.current_itheta != itheta:
            self.current_itheta = itheta
            self.signals.ithetaprogress.emit(itheta)
        ref_img = self.EdgeDetect(self.Gauss(self.parent.stack.absdata_cropped[:, :, data[0],itheta]))
        mov_img = self.EdgeDetect(self.Gauss(self.parent.stack.absdata_cropped[:, :, data[1],itheta]))
        upsampling = self.UpsamplingFactor()
        drift, error, _ = phase_cross_correlation(ref_img, mov_img,upsample_factor=upsampling,normalization=None)
        self.parent.stack.shiftsdict[itheta]["errors"][data[0]] = round(error,4)
        if self.parent.cb_upsampling.isChecked():
            self.parent.stack.shiftsdict[itheta]["xdots"][data[0]] = round(drift[0],2)
            self.parent.stack.shiftsdict[itheta]["ydots"][data[0]] = round(drift[1],2)
        else:
            self.parent.stack.shiftsdict[itheta]["xdots"][data[0]] = round(drift[0],0)
            self.parent.stack.shiftsdict[itheta]["ydots"][data[0]] = round(drift[1],0)
    def ShiftImg(self, row,x,y,itheta):
        borders, padded = self.PadImg(self.parent.stack.absdata4d[:, :, row,itheta],-x,-y)
        if self.parent.cb_upsampling.isChecked():
            shifted = ndimage.fourier_shift(np.fft.fft2(padded), [float(-x),float(-y)])
        else:
            shifted = ndimage.fourier_shift(np.fft.fft2(padded), [int(round(-x,0)),int(round(-y,0))])
        shifted = np.fft.ifft2(shifted)
        self.parent.stack.absdata4d_shifted[:, :, row, itheta] = shifted.real[borders[0]:padded.shape[0]-borders[1],borders[2]:padded.shape[1]-borders[3]]
        return

    def PadImg(self, img, x, y):
        default = 10 # minimum expansion of image
        borders = 4*[default]
        if x<0:
            borders[1] = abs(int(np.floor(x)))+default
        elif x>0:
            borders[0] = abs(int(np.ceil(x)))+default
        if y<0:
            borders[3] = abs(int(np.floor(y)))+default
        elif y>0:
            borders[2] = abs(int(np.ceil(y)))+default
        padded = np.pad(img,((borders[0],borders[1]),(borders[2],borders[3])),mode = "edge")
        return (borders, padded)
    def Gauss(self, im):
        gaussed = ndimage.gaussian_filter(im, self.parent.spinBoxGauss.value())
        return gaussed
    def EdgeDetect(self, im):
        if self.parent.cb_edgedetect.isChecked():
            im = filters.farid(im)
        return im
    def UpsamplingFactor(self):
        if self.parent.cb_upsampling.isChecked():
            fac = 20 ## equal to 0.05 px precision
        else:
            fac = 5 ## equal to 0.2 px precision
        return fac

# ----------------------------------------------------------------------
class TaskDispatcher(QtCore.QObject):
    def __init__(self,parent):
        print("0 - Task dispatcher called, pool initiated.")
        super(TaskDispatcher, self).__init__()
        try:
            self.worker.signals.ithetaprogress.disconnect()
            self.worker.signals.finished.disconnect()
        except:
            pass
        self.pool = QtCore.QThreadPool.globalInstance()
        try:
            cpus = len(os.sched_getaffinity(0)) # number of cpu threads. not supported on some platforms.
        except:
            cpus = os.cpu_count()
        self.pool.setMaxThreadCount(cpus)
        self.queue = SimpleQueue()
        self.parent = parent

    @pyqtSlot()
    def run(self):
        #print(self.queue.qsize())
        #print('pool', threading.get_ident())
        # if self.parent.com.stack_4d == 0:
        qsize = int(self.queue.qsize())
        maxthreads = self.pool.maxThreadCount()
        preferred_thread_number = min(maxthreads, qsize)
        print("2 - Starting threads: ",preferred_thread_number, "threads needed, number of tasks:",qsize)
        while int(self.queue.qsize()) and self.pool.activeThreadCount() < preferred_thread_number: #start as many threads as needed.
            #print("active threads"+str(self.pool.activeThreadCount())+" qsize "+str(int(self.queue.qsize())))
            worker = GeneralPurposeProcessor(self.parent,self.queue)
            worker.signals.ithetaprogress.connect(self.parent.IThetaProgress)
            self.pool.start(worker)
            time.sleep(0.02) # artifical delay to start threads with a little time separation. Otherwise ShiftImgs freezes in Win10 and MacOS for small stacks!
        self.pool.waitForDone()
        print("3 - All tasks are finished.")
        worker.signals.finished.connect(self.parent.ThreadPoolComplete)
        worker.signals.finished.emit()
    #add a task to the queue
    def enqueuetask(self, func, *args, **kargs):
        self.queue.put((func, args, kargs))
# ----------------------------------------------------------------------

class ImageRegistrationFFT(QtWidgets.QDialog, QtWidgets.QGraphicsScene):
    def __init__(self, parent, common, stack):
        QtWidgets.QWidget.__init__(self, parent)
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'showalign2.ui'), self)
        self.parent = parent
        self.stack = stack
        self.com = common
        self.iev = 0
        self.itheta = 0
        if self.com.stack_loaded == 1:
            if self.com.stack_4d:
                self.stack.absdata4d = self.stack.stack4D.copy()
                #self.button_alignbatch.setVisible(True)
                self.button_align.setText("Align Batch")
                self.thetathread = QtCore.QThread()
            else:
                self.stack.absdata4d = np.expand_dims(self.stack.absdata.copy(), axis=3)
                #self.button_alignbatch.setVisible(False)
                self.button_align.setText("Align")
            self.stack.absdata_cropped = self.stack.absdata4d.copy() # This image stack is used by the alignment/shift routine and cropped to the ROI rectangle
            self.stack.absdata4d_shifted = self.stack.absdata4d.copy() # This is the full-sized output of the alignment/shift routine
            self.stack.absdata4d_shifted_cropped = self.stack.absdata4d.copy() # This is the output cropped to the common region
        self.poolthread = QtCore.QThread()
        self.aligned = False
        self.SetupUI()

    def SetupUI(self):
        self.button_ok.setEnabled(False)
        self.button_rstroi.setEnabled(False)
        self.button_rstalign.setEnabled(False)
        self.button_preview.setEnabled(True)
        self.slider_theta.setVisible(False)
        self.setWindowTitle('FFT Stack Alignment')
        self.pglayout = pg.GraphicsLayout(border=None)
        self.canvas.setBackground("w") # canvas is a pg.GraphicsView widget
        self.canvas.setCentralWidget(self.pglayout)
        self.vb = self.pglayout.addViewBox()
        self.vb.setAspectLocked()
        self.i_item = pg.ImageItem(border="k",parent= self)

        self.vb.setMouseEnabled(x=False, y=False)
        self.vb.addItem(self.i_item, ignoreBounds=False)

        self.DriftsWidget.setBackground("w")

        self.py = self.DriftsWidget.addPlot(row=0, col=0, rowspan=1, colspan=1)
        self.py.setMouseEnabled(x=True, y=True)
        #self.i_item = pg.ImageItem(border="k")
        #self.py.setAspectLocked(lock=True, ratio=1)
        self.py.showAxis("top", show=True)
        self.py.showAxis("bottom", show=True)
        self.py.showAxis("left", show=True)
        self.py.showAxis("right", show=True)
        ay1 = self.py.getAxis("left")
        by1 = self.py.getAxis("right")
        ax1 = self.py.getAxis("bottom")
        bx1 = self.py.getAxis("top")
        ay1.setLabel(text="y-drift", units="px")
        ay1.enableAutoSIPrefix(enable=True)
        ay1.setWidth(w=60)
        ax1.setLabel(text="Photon Energy",units="eV")
        ax1.enableAutoSIPrefix(enable=True)
        ay1.setStyle(tickLength=8)
        ax1.setStyle(tickLength=0)
        #ax1.setHeight(h=46.2)
        ax1.setStyle(tickLength=8)
        by1.setStyle(showValues=False, tickLength=0)
        bx1.setStyle(showValues=False, tickLength=0)

        self.px = self.DriftsWidget.addPlot(row=1, col=0, rowspan=1, colspan=1)
        self.px.setXLink(self.py)
        self.px.setMouseEnabled(x=True, y=True)
        #self.i_item = pg.ImageItem(border="k")
        #self.px.setAspectLocked(lock=True, ratio=1)
        self.px.showAxis("top", show=True)
        self.px.showAxis("bottom", show=True)
        self.px.showAxis("left", show=True)
        self.px.showAxis("right", show=True)
        ay2 = self.px.getAxis("left")
        by2 = self.px.getAxis("right")
        ax2 = self.px.getAxis("bottom")
        bx2 = self.px.getAxis("top")
        ay2.setLabel(text="x-drift",units="px")
        ay2.enableAutoSIPrefix(enable=True)
        ay2.setWidth(w=60)
        ax2.setLabel(text="Photon Energy",units="eV")
        ax2.enableAutoSIPrefix(enable=True)
        ay2.setStyle(tickLength=8)
        ax2.setStyle(tickLength=8)
        by2.setStyle(showValues=False,tickLength=0)
        bx2.setStyle(showValues=False,tickLength=0)

        self.button_ok.clicked.connect(self.OnAccept)
        self.button_cancel.clicked.connect(self.OnCancel)

        if self.com.stack_loaded == 1:
            if self.com.stack_4d:
                self.slider_theta.setVisible(True)
                self.slider_theta.setRange(0, self.stack.n_theta - 1)
                self.slider_theta.valueChanged[int].connect(self.OnScrollTheta)
            #self.maskedvals = [True] * int(self.stack.n_ev)
            self.spinBoxError.setEnabled(False)
            self.slider_eng.sliderPressed.connect(self.ShowImage)
            self.slider_eng.sliderReleased.connect(self.ShowImage)
            self.button_preview.clicked.connect(self.ShowImage)
            self.spinBoxGauss.valueChanged.connect(self.ShowImage)
            self.cb_edgedetect.toggled.connect(self.ShowImage)
            self.slider_eng.valueChanged[int].connect(self.OnScrollEng)
            self.slider_eng.setRange(0, self.stack.n_ev - 1)
            self.refmarkerx = pg.InfiniteLine(angle=90, movable=False, pen=pg.mkPen(color="b", width=2, style=QtCore.Qt.DashLine))
            self.refmarkery = pg.InfiniteLine(angle=90, movable=False, pen=pg.mkPen(color="b", width=2, style=QtCore.Qt.DashLine))
            self.px.addItem(self.refmarkerx, ignoreBounds=True)
            self.py.addItem(self.refmarkery, ignoreBounds=True)

            self.fit_x = pg.PlotCurveItem(pen=pg.mkPen(color="c", width=2))
            self.fit_y = pg.PlotCurveItem(pen=pg.mkPen(color="c", width=2))
            self.fit_x.setZValue(200)
            self.fit_y.setZValue(200)
            self.px.addItem(self.fit_x, ignoreBounds=True)
            self.py.addItem(self.fit_y, ignoreBounds=True)

            self.OnScrollEng(0)
            self.SetupROI()
            self.button_align.clicked.connect(self.ComposeAlignQueue)
            self.xregion = pg.LinearRegionItem(brush=[255, 0, 0, 45], bounds=[self.stack.ev[0], self.stack.ev[-1]])
            self.yregion = pg.LinearRegionItem(brush=[255, 0, 0, 45], bounds=[self.stack.ev[0], self.stack.ev[-1]])
            self.xregion.setRegion([self.stack.ev[0], self.stack.ev[-1]])
            self.yregion.setRegion([self.stack.ev[0], self.stack.ev[-1]])
            self.yregion.setZValue(100)
            self.xregion.setZValue(100)
            self.px.addItem(self.xregion, ignoreBounds=False)
            self.py.addItem(self.yregion, ignoreBounds=False)
            self.xregion.sigRegionChangeFinished.connect(lambda region: self.OnLinRegion(region))
            self.yregion.sigRegionChangeFinished.connect(lambda region: self.OnLinRegion(region))

            self.xscatter = pg.ScatterPlotItem(pxMode=False)
            self.yscatter = pg.ScatterPlotItem(pxMode=False)

            self.InitShiftsDict()

            self.px.addItem(self.xscatter)
            self.py.addItem(self.yscatter)
            self.MakeNewScatterPlots()
            #self.shifts = self.stack.shifts.copy()
            self.xscatter.sigClicked.connect(self.OnPointClicked)
            self.yscatter.sigClicked.connect(self.OnPointClicked)
            self.cb_autocrop.toggled.connect(self.OnAutoCrop)
    def InitShiftsDict(self):
        outer_keys = range(max(self.stack.n_theta,1))
        inner_keys = ["xdots", "ydots", "xshifts", "yshifts", "errors", "errormaskedx","errormaskedy","manualmaskedx","manualmaskedy"]
        self.stack.shiftsdict = {intkey : {key: [False] * int(self.stack.n_ev) for key in inner_keys} for intkey in outer_keys}
        single_keys = ["filter", "method", "autoquality", "extrapolation", "threshold", "regionlimitx", "regionlimity"]
        for intkey in outer_keys:
            self.stack.shiftsdict[intkey].update({key: False for key in single_keys})
    def CreateScatterDots(self,shifts,mask):
        scatterdots = [{'pos': tup[0:2], 'size': 10,
                           #'pen': {'color': 'w', 'width': 2},
                           'brush': QtGui.QColor('red')} if tup[2] else {'pos': tup[0:2], 'size': 10,
                           #'pen': {'color': 'w', 'width': 2},
                           'brush': QtGui.QColor('blue')} for tup in list(zip(self.stack.ev, shifts, mask))]
        return scatterdots

    def OnPointClicked(self, obj, points): # Manually add/remove points to/from fit if in selected region
        selectscatter = {self.xscatter: ["errormaskedx","manualmaskedx",self.xregion], self.yscatter: ["errormaskedy","manualmaskedy",self.yregion]}
        idx = points[0].index()
        mask1 = self.stack.shiftsdict[self.itheta][selectscatter[obj][1]]
        mask2 = self.stack.shiftsdict[self.itheta][selectscatter[obj][0]]
        min_idx, max_idx, *_ = self.getDataClosestToRegion(selectscatter[obj][2], obj)
        if min_idx <= idx <= max_idx :
            mask1[idx] = not np.logical_or(mask1[idx],mask2[idx])
            mask2[idx] = mask1[idx]
            self.ColorizeScatterDots()
            self.OnScrollEng(points[0].index())
            if self.aligned:
                self.ComposeShiftQueue()
    def OnAligned(self):
        self.aligned = True
        self.cb_autoerror.stateChanged.disconnect()
        self.cb_extrapolate.stateChanged.disconnect()
        self.spinBoxFiltersize.valueChanged.disconnect()
        self.spinBoxError.valueChanged.disconnect()
        self.comboBox_approx.currentIndexChanged.disconnect()
        self.MakeNewScatterPlots()
        self.OnLinRegion(self.xregion, update=False)
        self.OnLinRegion(self.yregion, update=False)
        self.ComposeShiftQueue(init=True)

    # ----------------------------------------------------------------------
    def MaskedScatterDotsArray(self,region):
        selection = {self.xregion : ["errormaskedx","manualmaskedx"], self.yregion : ["errormaskedy","manualmaskedy"]}
        array = np.logical_or(self.stack.shiftsdict[self.itheta][selection[region][0]], self.stack.shiftsdict[self.itheta][selection[region][1]])
        return array
    def MakeNewScatterPlots(self):
        errorthreshold = self.stack.shiftsdict[self.itheta]["threshold"]
        errors = self.stack.shiftsdict[self.itheta]["errors"]
        self.MaskScatterDotsAboveErrorThreshold((errors,errorthreshold,self.itheta))
        maskedx = self.MaskedScatterDotsArray(self.xregion)
        maskedy = self.MaskedScatterDotsArray(self.yregion)

        self.cb_autoerror.stateChanged.connect(self.OnAutoError)
        self.cb_extrapolate.stateChanged.connect(self.OnExtrapolate)
        self.spinBoxFiltersize.valueChanged.connect(self.OnFilter)
        self.spinBoxError.valueChanged.connect(lambda value: self.OnSpinBoxError(value))
        self.comboBox_approx.currentIndexChanged.connect(lambda: self.ComposeShiftQueue(init=False))

        self.spinBoxError.blockSignals(True)
        self.spinBoxError.setDecimals(4)
        self.spinBoxError.setMinimum(np.partition(errors, 1)[1])  # makes sure that at least two elements are selected
        self.spinBoxError.setStepType(QtWidgets.QAbstractSpinBox.AdaptiveDecimalStepType)
        self.OnAutoError()
        self.spinBoxError.blockSignals(False)

        xdots = self.stack.shiftsdict[self.itheta]["xdots"]
        ydots = self.stack.shiftsdict[self.itheta]["ydots"]
        if self.rB_consecutive.isChecked(): # pixel shifts relative to neighboring image
            xdots = self.ConvertToCumsum(xdots)
            ydots = self.ConvertToCumsum(ydots)
        self.xscatter.setData(spots=self.CreateScatterDots(xdots,maskedx), pxMode=True)
        self.yscatter.setData(spots=self.CreateScatterDots(ydots,maskedy), pxMode=True)

    def ConvertToCumsum(self,shifts): # Converts list of shifts relative to neighboring images to the cumulative sum of shifts relative to the reference image
        f = [i for i, x in enumerate(shifts) if isinstance(x, bool)][0]  # A little bit hackish. Returns the index of the first boolean value from a list, i.e., the reference image.
        ll = np.flip(shifts[:f]) # shifts on the left, lower energy side of the reference frame
        rl= shifts[f:] # shifts on the right, higher energy side of the reference frame
        lcs = np.flip(np.cumsum(ll))
        rcs = np.cumsum(rl)
        return np.concatenate((lcs, rcs))
    def OnAutoError(self):
        # print("onautoerror",self.stack.shiftsdict[self.itheta]["threshold"])
        self.stack.shiftsdict[self.itheta]["autoquality"] = self.cb_autoerror.isChecked()
        self.stack.shiftsdict[self.itheta]["threshold"] = round(np.mean(self.stack.shiftsdict[self.itheta]["errors"]),
                                                                4)
        if self.cb_autoerror.isChecked() and self.aligned:
            self.spinBoxError.setEnabled(True)
            self.spinBoxError.setValue(self.stack.shiftsdict[self.itheta]["threshold"])
        else:
            self.spinBoxError.setEnabled(False)
            self.spinBoxError.setValue(1)
    def OnSpinBoxError(self,value):
        # print("onspinboxerror")
        self.stack.shiftsdict[self.itheta]["threshold"] = value
        # print(self.stack.shiftsdict[self.itheta]["threshold"])
        self.MaskScatterDotsAboveErrorThreshold((self.stack.shiftsdict[self.itheta]["errors"], value, self.itheta))
        self.ColorizeScatterDots()
        self.ComposeShiftQueue(init=False)
    def OnExtrapolate(self):
        self.ComposeShiftQueue(init=False)
    def OnFilter(self):
        self.ComposeShiftQueue(init=False)

    def ColorizeScatterDots(self):
        brushesx = [QtGui.QColor('red') if bool else QtGui.QColor('blue') for bool in self.MaskedScatterDotsArray(self.xregion)]
        brushesy = [QtGui.QColor('red') if bool else QtGui.QColor('blue') for bool in self.MaskedScatterDotsArray(self.yregion)]
        self.xscatter.setBrush(brushesx,update=True)
        self.yscatter.setBrush(brushesy,update=True)

    def MaskScatterDotsAboveErrorThreshold(self, errorvals=None):
        errors, errorthreshold, itheta = errorvals
        if self.cb_autoerror.isChecked():
            self.stack.shiftsdict[itheta]["errormaskedx"] = (errors > np.float64(errorthreshold))
            self.stack.shiftsdict[itheta]["errormaskedy"] = (errors > np.float64(errorthreshold))
        else:
            self.stack.shiftsdict[itheta]["errormaskedx"] = [False] * len(errors)
            self.stack.shiftsdict[itheta]["errormaskedy"] = [False] * len(errors)

    # ----------------------------------------------------------------------
    def initParams(self, itheta):
        self.stack.shiftsdict[itheta]["regionlimitx"] = self.getDataClosestToRegion(self.xregion,self.xscatter,False)[2:]
        self.stack.shiftsdict[itheta]["regionlimity"] = self.getDataClosestToRegion(self.yregion,self.yscatter,False)[2:]
        self.stack.shiftsdict[itheta]["filter"] = self.spinBoxFiltersize.value()
        self.stack.shiftsdict[itheta]["method"] = self.comboBox_approx.currentIndex()
        #print(self.stack.shiftsdict[itheta]["method"])
        #self.stack.shiftsdict[itheta]["autoquality"] = self.cb_autoerror.isChecked()
        self.stack.shiftsdict[itheta]["extrapolation"] = self.cb_extrapolate.isChecked()
        self.stack.shiftsdict[itheta]["threshold"] = round(np.mean(self.stack.shiftsdict[itheta]["errors"]), 4)
    def restoreParams(self,itheta):
        self.xregion.blockSignals(True)
        self.yregion.blockSignals(True)
        self.cb_autoerror.blockSignals(True)
        self.cb_extrapolate.blockSignals(True)
        self.comboBox_approx.blockSignals(True)
        self.spinBoxFiltersize.blockSignals(True)

        self.xregion.setRegion(self.stack.shiftsdict[itheta]["regionlimitx"])
        self.yregion.setRegion(self.stack.shiftsdict[itheta]["regionlimity"])
        self.cb_autoerror.setChecked(self.stack.shiftsdict[itheta]["autoquality"])
        self.spinBoxError.blockSignals(True)
        #self.spinBoxError.setValue(self.stack.shiftsdict[self.itheta]["threshold"])
        if self.cb_autoerror.isChecked() and self.aligned:
            self.spinBoxError.setEnabled(True)
            self.spinBoxError.setValue(self.stack.shiftsdict[self.itheta]["threshold"])
        else:
            self.spinBoxError.setEnabled(False)
            self.spinBoxError.setValue(1)
        self.spinBoxError.blockSignals(False)
        self.cb_extrapolate.setChecked(self.stack.shiftsdict[itheta]["extrapolation"])
        self.comboBox_approx.setCurrentIndex(self.stack.shiftsdict[itheta]["method"])
        if self.comboBox_approx.currentIndex() == 0:
            self.spinBoxFiltersize.setEnabled(True)
        elif self.comboBox_approx.currentIndex() == 1: # linear regression
            self.spinBoxFiltersize.setEnabled(False)
        self.spinBoxFiltersize.setValue(self.stack.shiftsdict[itheta]["filter"])
        # if limits:
        #     region.setRegion(list(limits))
        # else:
        #     region.setRegion([self.stack.ev[0],self.stack.ev[-1]])
        #self.spinBoxError.blockSignals(False)
        self.xregion.blockSignals(False)
        self.yregion.blockSignals(False)
        self.cb_autoerror.blockSignals(False)
        self.cb_extrapolate.blockSignals(False)
        self.comboBox_approx.blockSignals(False)
        self.spinBoxFiltersize.blockSignals(False)

    def getDataClosestToRegion(self,region,plotitem,snapregion=False):
        selectregion = {self.xregion : "regionlimitx", self.yregion : "regionlimity"}
        limits = selectregion[region]
        minidx, maxidx = region.getRegion()
        data = plotitem.getData()[0]
        index = lambda x: np.argmin(np.abs(data - x))
        minidx = index(minidx)
        maxidx = index(maxidx)
        if minidx == maxidx:
            minidx = 0
            maxidx = np.argmax(data)
        mindata = data[minidx]
        maxdata = data[maxidx]
        if snapregion:
            self.stack.shiftsdict[self.itheta][limits] = (mindata,maxdata)  # snap region to data points
            region.setRegion((mindata,maxdata))
        return minidx, maxidx, mindata, maxdata
    
    def ApplyApproximationFunction(self,region):
        boolarray = self.MaskedScatterDotsArray(region)
        selectscatter = {self.xregion: [self.xscatter, "x"] , self.yregion: [self.yscatter, "y"]}
        selectfit= {self.xregion : self.fit_x, self.yregion : self.fit_y}
        scatter = selectscatter[region][0]
        fit = selectfit[region]
        #min_idx, max_idx, *_ = self.getDataClosestToRegion(region,scatter)
        selected = np.count_nonzero(boolarray == False)
        xdata, ydata = scatter.getData()
        #xdata = xdata[min_idx:max_idx]
        xdata = xdata[~boolarray]
        #ydata = ydata[min_idx:max_idx]
        ydata = ydata[~boolarray]
        #print(selected)
        if selected < 2:
            return []
        if self.comboBox_approx.currentIndex() == 0: # moving average
            self.spinBoxFiltersize.setEnabled(True)
            approximated= ndimage.filters.uniform_filter1d(ydata,self.spinBoxFiltersize.value(),mode = "nearest")
        elif self.comboBox_approx.currentIndex() == 1: # linear regression
            self.spinBoxFiltersize.setEnabled(False)
            reg = linregress([xdata,ydata])
            approximated = [reg.slope * i + reg.intercept for i in xdata]

        if self.cb_extrapolate.isChecked():
            fillval= "extrapolate"
        else:
            fillval=(approximated[0],approximated[-1])

        interpolate_func = interp1d(xdata, approximated,kind="linear",fill_value=fillval,bounds_error=False)
        fitdata = [self.stack.ev,np.around(interpolate_func(self.stack.ev),2)] # round fit to two decimals, i.e., 1/self.upsampling of a px
        fit.setData(x=fitdata[0], y=fitdata[1])
        fit.show()
        return fitdata[1]

    def resetPoolThread(self):
        try:
            self.poolthread.started.disconnect()
            self.poolthread.quit()
            self.poolthread.wait()
        except:
            pass
        self.pool = TaskDispatcher(self)
        self.pool.moveToThread(self.poolthread) # GUI is not blocking during calculation due to this
        self.poolthread.started.connect(self.pool.run)

    def IThetaProgress(self,itheta):
        #print("ithetaprogress")
        # Each thread calls this function. The condition prevents multiple calls.
        if self.slider_theta.value() != itheta:
            #self.slider_theta.blockSignals(True)
            self.slider_theta.setValue(itheta)
            #self.slider_theta.blockSignals(False)
        #print(self.pool.pool.activeThreadCount())

    def ThreadPoolComplete(self):
        #print(str(self.pool.pool.activeThreadCount())+" THREADS REMAINING FROM POOL.")
        #self.slider_theta.setValue(0)
        if not self.aligned and self.pool.pool.activeThreadCount() == 0:
            print("-------------")
            self.OnAligned()

        elif self.aligned and self.pool.pool.activeThreadCount() == 0:
            self.OnAutoCrop()
            print("-------------")
        QtWidgets.QApplication.restoreOverrideCursor()

    def ComposeAlignQueue(self):
        self.button_align.setEnabled(False)
        ref_idx = self.iev
        # Reset reference img:
        self.resetPoolThread()
        # Check scikit-image version:
        if parse_version(version_check('scikit-image')) >= parse_version('0.19.1'):
            AlignFunc = "AlignReferenced"
        else:
            print("scikit-image version <= 0.19.1 is not supported")
            return

        itheta = 0
        #necessary work around for 3d stacks and if 4d stack is loaded with LoadStack()
        if self.com.stack_4d != 1:
            ntheta = 1
        else:
            ntheta = max(self.stack.n_theta, 1)
        while itheta < ntheta:
            idx = copy.copy(self.stack.n_ev)
            while idx: # Generate pairs of indices starting at reference image index.
                next = ref_idx + (self.stack.n_ev - idx)
                prev = ref_idx - (self.stack.n_ev - idx)
                running = 2
                if self.rB_referenced.isChecked(): # compare reference frame to next and previous frame.
                    if next < self.stack.n_ev-1:
                        self.pool.enqueuetask(AlignFunc, (next + 1, ref_idx), itheta)
                    else:
                        running -= 1
                    if prev > 0:
                        self.pool.enqueuetask(AlignFunc, (prev - 1, ref_idx), itheta)
                    else:
                        running -= 1
                    if running:
                        idx -= 1
                    else:
                        break
                elif self.rB_consecutive.isChecked(): # compare frame to next and previous frame.
                    if next < self.stack.n_ev-1:
                        self.pool.enqueuetask(AlignFunc, (next + 1, next), itheta)
                    else:
                        running -= 1
                    if prev > 0:
                        self.pool.enqueuetask(AlignFunc, (prev - 1, prev), itheta)
                    else:
                        running -= 1
                    if running:
                        idx -= 1
                    else:
                        break
            itheta = itheta + 1
        if not self.pool.queue.empty():
            print("1 - Task queue composed. Starting poolthread for displacement estimation.")
            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
            self.poolthread.start()
        
    def ComposeShiftQueue(self, init=False):
        # print("composeshift initial?",str(init))
        self.resetPoolThread()
        #array = np.logical_or(self.MaskedScatterDotsArray(self.xregion),self.MaskedScatterDotsArray(self.yregion))
        if not init:
            ntheta = [self.itheta]
        else:
            if self.com.stack_4d != 1:
                ntheta = range(1)
            else:
                ntheta = range(max(self.stack.n_theta,1)) # necessary work around for 3d stacks and if 4d stack is loaded with LoadStack()
            for itheta in ntheta:
                self.initParams(itheta)
        for itheta in ntheta:
            #print(itheta)
            self.stack.shiftsdict[self.itheta]["filter"] = self.spinBoxFiltersize.value()
            self.stack.shiftsdict[self.itheta]["method"] = self.comboBox_approx.currentIndex()
            self.slider_theta.setValue(itheta)
            #time.sleep(0.2)
            xshifts = self.ApplyApproximationFunction(self.xregion)
            yshifts = self.ApplyApproximationFunction(self.yregion)
            # if empty arrays, i.e., if less than 2 dots selected, do nothing
            if not any(xshifts) and not any(yshifts):
                return
            for ev in range(self.stack.n_ev):
                # Only enqueue after reset or if new shift value is different to previous shift value
                if (isinstance(self.stack.shiftsdict[itheta]["xshifts"][ev], bool)
                        or isinstance(self.stack.shiftsdict[itheta]["yshifts"][ev], bool)
                        or (self.stack.shiftsdict[itheta]["xshifts"][ev],self.stack.shiftsdict[itheta]["yshifts"][ev]) != (xshifts[ev],yshifts[ev])):
                    self.stack.shiftsdict[itheta]["xshifts"][ev] = xshifts[ev]
                    self.stack.shiftsdict[itheta]["yshifts"][ev] = yshifts[ev]
                    self.pool.enqueuetask("ShiftImg", ev, xshifts[ev], yshifts[ev],itheta)
        if not self.pool.queue.empty():
            print("1 - Task queue composed. Starting poolthread for image shifting.")
            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
            self.poolthread.start()

    # ----------------------------------------------------------------------
    def GetIndexPairs(self):
        idxtuplelst = [(i+1,i) for i in range(self.stack.n_ev-1)]
        return idxtuplelst
    # ----------------------------------------------------------------------
    def OnAutoCrop(self):
        if self.aligned:
            self.button_preview.blockSignals(True)
            self.button_preview.setChecked(False)
            self.button_preview.blockSignals(False)
            self.cb_autocrop.blockSignals(True)
            self.CropStack4D()
            self.cb_autocrop.blockSignals(False)
            self.OnScrollEng(self.iev)
            self.button_rstroi.setEnabled(True)
            self.button_rstalign.setEnabled(True)
            self.button_ok.setEnabled(True)

    def CropStack4D(self):
        self.stack.absdata4d_shifted_cropped = self.stack.absdata4d_shifted.copy()
        if self.cb_autocrop.isChecked():
            ntheta = range(max(self.stack.n_theta,1)) # necessary work around for 3d stacks and if 4d stack is loaded with LoadStack()
            globalminx = min([j for i in [self.stack.shiftsdict[theta]["xshifts"] for theta in ntheta] for j in i])
            globalmaxx = max([j for i in [self.stack.shiftsdict[theta]["xshifts"] for theta in ntheta] for j in i])
            globalminy = min([j for i in [self.stack.shiftsdict[theta]["yshifts"] for theta in ntheta] for j in i])
            globalmaxy = max([j for i in [self.stack.shiftsdict[theta]["yshifts"] for theta in ntheta] for j in i])
            if globalminx == globalmaxx == globalminy == globalmaxy == False:
                return
            #if not self.com.stack_4d:
            self.box.hide()
            l = -int(np.floor(globalminx))
            r = -int(np.ceil(globalmaxx))
            cr = r if r < 0 else None
            if l < 0:
                l = 0
                cr = cr - l
            b = -int(np.floor(globalminy))
            t = -int(np.ceil(globalmaxy))
            ct = t if t < 0 else None
            if b < 0:
                b = 0
                ct = ct - b
            if 0 in self.stack.absdata4d_shifted_cropped[l:cr,b:ct,:,self.itheta].shape:
                QtWidgets.QMessageBox.warning(self, 'Error', 'The alignment failed. Cropping would result in a zero-dimensional image. Please check your settings. Auto-crop has been disabled. You can re-enable it manually.')
                self.cb_autocrop.setChecked(False)
                self.OnResetAlignROI()
            else:
                self.stack.absdata4d_shifted_cropped = self.stack.absdata4d_shifted_cropped[l:cr,b:ct,:,:]
                if self.com.i0_loaded == 1:
                    self.stack.i0_mask = self.stack.i0_mask[l: cr, b: ct]
                    print("4 - Stack & I0 mask cropped to common region: ", globalminx, globalmaxx, globalminy, globalmaxy)
                    return
                print("4 - Stack cropped to common region: ", globalminx, globalmaxx, globalminy, globalmaxy)
        else:
            self.box.show()

    def OnScrollEng(self, value):
        self.slider_eng.setValue(value)
        self.iev = value
        self.ShowImage()
        self.refmarkerx.setValue(self.stack.ev[self.iev])
        self.refmarkery.setValue(self.stack.ev[self.iev])
    def OnScrollTheta(self, value):
        #if value != self.itheta:
        #print("Theta slider" + str(value))
        self.slider_theta.setValue(value)
        self.itheta = value
        #self.ClearShifts()
        self.ShowImage()
        if self.aligned:
            #self.MakeNewScatterPlots()
            errorthreshold = self.stack.shiftsdict[self.itheta]["threshold"]
            errors = self.stack.shiftsdict[self.itheta]["errors"]
            if not errorthreshold:
                errorthreshold = round(np.mean(errors), 4)
            self.restoreParams(self.itheta)
            self.MaskScatterDotsAboveErrorThreshold((errors, errorthreshold, self.itheta))
            # if self.cb_autoerror.isChecked():
            #     self.spinBoxError.blockSignals(True)
            #     self.spinBoxError.setValue(errorthreshold)
            #     self.spinBoxError.blockSignals(False)
            xdots = self.stack.shiftsdict[self.itheta]["xdots"]
            ydots = self.stack.shiftsdict[self.itheta]["ydots"]
            maskedx = self.MaskedScatterDotsArray(self.xregion)
            maskedy = self.MaskedScatterDotsArray(self.yregion)
            self.xscatter.setData(spots=self.CreateScatterDots(xdots, maskedx), pxMode=True)
            self.yscatter.setData(spots=self.CreateScatterDots(ydots, maskedy), pxMode=True)
            self.fit_x.setData(x=self.stack.ev, y=self.stack.shiftsdict[self.itheta]["xshifts"])
            self.fit_y.setData(x=self.stack.ev, y=self.stack.shiftsdict[self.itheta]["yshifts"])

            #min_idx, max_idx, min_ev, max_ev = self.getDataClosestToRegion(self.xregion, self.xscatter, True)
            #min_idx, max_idx, min_ev, max_ev = self.getDataClosestToRegion(self.yregion, self.yscatter, True)

            #self.fit_x.show()
            self.OnLinRegion(self.xregion, update=False)
            self.OnLinRegion(self.yregion, update=False)

    def ShowImage(self):
        self.stack.absdata_shifted_cropped = self.stack.absdata4d_shifted_cropped[:, :, :, int(self.itheta)]
        im = self.stack.absdata_shifted_cropped[:, :, int(self.iev)]
        if self.button_preview.isChecked():
            im = ndimage.gaussian_filter(im, self.spinBoxGauss.value())
            if self.cb_edgedetect.isChecked():
                im = filters.farid(im)
        self.i_item.setImage(im)
        if self.com.stack_4d == 1:
            self.groupBox.setTitle(str('Stack Browser | Image at {0:5.2f} eV and {1:5.1f}°').format(float(self.stack.ev[self.iev]),float(self.stack.theta[self.itheta]), ))
        else:
            self.groupBox.setTitle(str('Stack Browser | Image at {0:5.2f} eV').format(float(self.stack.ev[self.iev])))

    ## Setup a ROI for an alignment rectangle. By default the whole image area is used.
    def SetupROI(self):
        self.box = pg.RectROI(self.i_item.boundingRect().topLeft(), self.i_item.boundingRect().bottomRight(),
                              pen=(5, 8), handlePen=QtGui.QPen(QtGui.QColor(255, 0, 128, 255)), centered=False,
                              sideScalers=False, removable=False, scaleSnap=True, translateSnap=True,
                              maxBounds=self.i_item.boundingRect())
        self.vb.addItem(self.box, ignoreBounds=False)
        self.box.sigRegionChangeFinished.connect(self.OnBoxChanged)
        self.box.sigRegionChangeStarted.connect(self.OnBoxChanging)
        self.button_rstroi.clicked.connect(self.OnResetAlignROI)
        self.button_rstalign.clicked.connect(self.OnResetAlign)
    ## The ROI is limited to the visible image area. OnMouseMoveOutside handles the behavior when
    def ClearShifts(self):
        self.aligned = False
        self.cb_autoerror.stateChanged.disconnect()
        self.fit_x.hide()
        self.fit_y.hide()
        self.InitShiftsDict()
        self.button_align.setEnabled(True)
        self.MakeNewScatterPlots()
    def OnResetAlignROI(self):
        self.ClearShifts()
        self.stack.absdata4d_shifted_cropped = self.stack.absdata4d.copy()
        self.OnScrollEng(self.iev)
        self.box.setPos(0, 0, update=False, finish=False)
        self.box.setSize(self.i_item.boundingRect().bottomRight() - self.box.pos(), update=True, snap=True, finish=True)
        self.box.show()
        self.button_rstroi.setEnabled(False)
        self.button_rstalign.setEnabled(False)
        self.button_ok.setEnabled(False)
    def OnResetAlign(self):
        self.ClearShifts()
        self.stack.absdata4d_shifted_cropped = self.stack.absdata4d.copy()
        self.OnScrollEng(self.iev)
        self.box.show()
        self.button_rstalign.setEnabled(False)
        self.button_ok.setEnabled(False)
    def OnBoxChanging(self):
        self.boxsize = self.box.size()
        self.proxy = pg.SignalProxy(self.vb.scene().sigMouseMoved, rateLimit=30, slot=self.OnMouseMoveOutside)
    def OnBoxChanged(self):
        try:
            self.proxy.disconnect()
        except AttributeError:
            pass
        left = int(self.box.pos().x())
        right = left + int(self.box.size().x())
        bottom = int(self.box.pos().y())
        top = bottom + int(self.box.size().y())
        self.stack.absdata_cropped = self.stack.absdata4d[left:right, bottom:top, :,:].copy()

    def OnLinRegion(self, region, update=True):
        #print("OnLinRegion", str(update))
        selectregion= {self.xregion : "manualmaskedx", self.yregion : "manualmaskedy"}
        selectscatter = {self.xregion : self.xscatter, self.yregion : self.yscatter}
        selectregionlimit= {self.xregion : "regionlimitx", self.yregion : "regionlimity"}
        selectplot = {self.xregion : self.px, self.yregion : self.py}
        #print("onlinregion "+ selectregion[region])
        scatter = selectscatter[region]
        min_idx, max_idx, min_ev, max_ev = self.getDataClosestToRegion(region,scatter,update)
        selection = [*range(min_idx, max_idx+1)]
        for idx,val in enumerate(self.stack.shiftsdict[self.itheta]["manualmaskedx"]):
            if idx not in selection:
                self.stack.shiftsdict[self.itheta][selectregion[region]][idx] = True
            else:
                self.stack.shiftsdict[self.itheta][selectregion[region]][idx] = False
        y_vals= selectscatter[region].data["y"][min_idx:max_idx + 1]
        selectplot[region].setRange(yRange=[np.min(y_vals), np.max(y_vals)], disableAutoRange=True, padding=0.1)
        self.ColorizeScatterDots()
        #filter = [True if idx in selection else False for idx,bool in enumerate(self.stack.shiftsdict[self.itheta]["manualmaskedx"])]
        #print(selection,filter)
        #self.UpdateScatterPlots(region, id)
        if not self.aligned:
            for theta in range(self.stack.n_theta):
                self.stack.shiftsdict[theta][selectregionlimit[region]] = (min_ev, max_ev)
        if self.aligned and update:
            self.stack.shiftsdict[self.itheta][selectregionlimit[region]] = (min_ev, max_ev)
            #print("ComposeShiftQueue"+selectregion[region])
            self.ComposeShiftQueue()

    def OnMouseMoveOutside(self, ev):
        mousepos = self.vb.mapSceneToView(ev[0])
        if not self.vb.itemBoundingRect(self.i_item).contains(mousepos):
            maxrect = self.i_item.boundingRect().bottomRight()
            # if bounds exceeded
            out_x = max((mousepos-maxrect).x(),0)
            out_y = max((mousepos-maxrect).y(), 0)
            if self.box.size() != self.boxsize: # prevents taking action when the box is just dragged and not resized
                if out_x and out_y:
                    self.box.setSize(self.i_item.boundingRect().bottomRight()-self.box.pos(), update=True, snap=True, finish=False)
                elif out_x:
                    self.box.setSize([self.i_item.boundingRect().right()-self.box.pos().x(),mousepos.y()-self.box.pos().y()], update=True, snap=True, finish=False)
                elif out_y:
                    self.box.setSize([mousepos.x()-self.box.pos().x(),self.i_item.boundingRect().bottom()-self.box.pos().y()], update=True, snap=True, finish=False)
# ----------------------------------------------------------------------
    def OnCancel(self, evt):
        self.stack.absdata_shifted_cropped = self.stack.absdata
        self.parent.page1.absimgfig.loadNewImageWithROI()
        self.parent.page0.absimgfig.loadNewImage()

        #if showmaptab:
        #    self.parent.page9.Clear()
        #    self.parent.page9.loadData()

        self.close()
    # ----------------------------------------------------------------------
    def OnAccept(self, evt):
        if self.com.stack_4d == 0:
            self.stack.absdata = self.stack.absdata4d_shifted_cropped[:,:,:,0]
            self.stack.data_struct.exchange.data = self.stack.absdata
        else:
            self.stack.stack4D = self.stack.absdata4d_shifted_cropped
            self.stack.absdata = self.stack.stack4D[:, :, :, self.itheta]
            self.stack.data_struct.exchange.data = self.stack.stack4D
            #QtWidgets.QMessageBox.warning(self, 'Error', '4D stack not yet supported.')

        datadim = np.int32(self.stack.absdata.shape)

        self.stack.n_cols = datadim[0].copy()
        self.stack.n_rows =  datadim[1].copy()

        if self.com.i0_loaded == 1:
            self.parent.page1.specfig.I0Update()

        self.stack.data_struct.exchange.energy = self.stack.ev

        #ToDo: Handshake shift data for spectral roi
        #self.stack.data_struct.spectromicroscopy.xshifts = self.x_shiftstemp
        #self.stack.data_struct.spectromicroscopy.yshifts = self.y_shiftstemp

        #self.parent.page1.slider_eng.setRange(0,self.stack.n_ev-1)
        #self.parent.page1.iev = int(self.stack.n_ev/2)
        #self.parent.page1.slider_eng.setValue(self.parent.page1.iev)
        self.parent.page0.absimgfig.loadNewImage()
        self.parent.page1.absimgfig.loadNewImageWithROI()
        self.parent.page1.specfig.ClearandReload()
        self.parent.window().refresh_widgets()
        #self.parent.page1.ix = int(self.stack.n_cols/2)
        #self.parent.page1.iy = int(self.stack.n_rows/2)

        #self.parent.page1.loadSpectrum(self.parent.page1.ix, self.parent.page1.iy)
        #self.parent.page1.loadImage()

        #if showmaptab:
        #    self.parent.page9.Clear()
        #    self.parent.page9.loadData()

        self.close()

