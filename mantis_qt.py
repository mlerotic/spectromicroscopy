# 
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
# 
#   Copyright (C) 2013 Mirna Lerotic, 2nd Look
#   http://2ndlookconsulting.com
#   License: GNU GPL v3
#
#   Mantis is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   any later version.
#
#   Mantis is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details <http://www.gnu.org/licenses/>.


import sys
import os
import numpy as np
import time
import getopt

from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import *
from PyQt4.QtCore import Qt, QCoreApplication

from PIL import Image  

import Tkinter
import FileDialog


import matplotlib 
from numpy import NAN
matplotlib.rcParams['backend.qt4'] = 'PyQt4'
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid import make_axes_locatable
from matplotlib.widgets import LassoSelector
matplotlib.interactive( True )
matplotlib.rcParams['svg.fonttype'] = 'none'

import data_struct
import data_stack
import analyze
import nnma
import henke
import tomo_reconstruction

from helpers import resource_path
import file_plugins
from file_plugins import file_xrm
from file_plugins import file_bim
from file_plugins import file_dataexch_hdf5


version = '2.2.01'

Winsizex = 1000
Winsizey = 700

PlotH = 4.0
PlotW = PlotH*1.61803

ImgDpi = 40

verbose = True

showtomotab = 1



            
#----------------------------------------------------------------------
class common:
    def __init__(self):
        
        self.stack_loaded = 0
        self.stack_4d = 0
        self.i0_loaded = 0
        self.pca_calculated = 0
        self.pca4D_calculated = 0
        self.cluster_calculated = 0
        self.spec_anl_calculated = 0
        self.spec_anl4D_calculated = 0
        self.ica_calculated = 0
        self.xpf_loaded = 0
        
        self.white_scale_bar = 0
        
        self.path = ''
        self.filename = ''

        self.font = ''
        
        
""" ------------------------------------------------------------------------------------------------"""
class PageTomo(QtGui.QWidget):
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
        
        self.fulltomorecdata = []
        self.ncomponents = 0
        
        self.tomo_calculated = 0
        self.full_tomo_calculated = 0
        self.energiesloaded = 0
        
        self.datanames = []
        
        self.haveROI = 0
        self.ROIvol = []
        
        self.select1 = 0        

        self.icomp = 0
        self.islice = 0
        
        self.maxIters = 10
        self.beta = 0.5
        self.samplethick = 0
        

        #panel 1
        sizer1 = QtGui.QGroupBox('Tomo Data')
        vbox1 = QtGui.QVBoxLayout()
        
        self.button_spcomp = QtGui.QPushButton('Load Tomo Data for Spectral Components')
        self.button_spcomp.clicked.connect( self.OnLoadTomoComponents)
        self.button_spcomp.setEnabled(False)
        vbox1.addWidget(self.button_spcomp)
        
        self.button_engdata = QtGui.QPushButton('Load Tomo Data for each Energy')
        self.button_engdata.clicked.connect( self.OnLoadTomoEng)
        self.button_engdata.setEnabled(False)
        vbox1.addWidget(self.button_engdata)
        

        sizer1.setLayout(vbox1)

        #panel 2
        sizer2 = QtGui.QGroupBox('Compressed Sensing Reconstruction')
        vbox2 = QtGui.QVBoxLayout()

        
        self.button_calc1 = QtGui.QPushButton( 'Calculate One Dataset')
        self.button_calc1.clicked.connect( self.OnCalcTomo1)
        self.button_calc1.setEnabled(False)
        vbox2.addWidget(self.button_calc1)
        
        hbox21 = QtGui.QHBoxLayout()
        tc1 = QtGui.QLabel(self)
        hbox21.addWidget(tc1)
        tc1.setText('Choose Tomo Dataset: ')
        self.combonames = QtGui.QComboBox(self)
        self.combonames.activated[int].connect(self.OnSelect1Comp)
        hbox21.addWidget(self.combonames)
        
        vbox2.addLayout(hbox21)
        
        self.button_calcall = QtGui.QPushButton( 'Calculate All Datasets')
        self.button_calcall.clicked.connect( self.OnCalcTomoFull)
        self.button_calcall.setEnabled(False)
        vbox2.addWidget(self.button_calcall)
        
        
        self.button_save = QtGui.QPushButton( 'Save as .mrc')
        self.button_save.clicked.connect( self.OnSave)
        self.button_save.setEnabled(False)
        vbox2.addWidget(self.button_save)
        
        self.button_saveall = QtGui.QPushButton( 'Save All as .mrc')
        self.button_saveall.clicked.connect( self.OnSaveAll)
        self.button_saveall.setEnabled(False)
        vbox2.addWidget(self.button_saveall)
        
        line = QtGui.QFrame()
        line.setFrameShape(QtGui.QFrame.HLine)
        line.setFrameShadow(QtGui.QFrame.Sunken) 
        

        vbox2.addStretch(1)
        vbox2.addWidget(line) 
        vbox2.addStretch(1)  
 
        hbox22 = QtGui.QHBoxLayout()
        text1 = QtGui.QLabel(self)
        text1.setText('Number of iterations')
        hbox22.addWidget(text1)
        hbox22.addStretch(1)
        self.ntc_niterations = QtGui.QLineEdit(self)
        self.ntc_niterations.setFixedWidth(65)
        self.ntc_niterations.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.ntc_niterations.setAlignment(QtCore.Qt.AlignRight)   
        self.ntc_niterations.setText(str(self.maxIters))
        hbox22.addWidget(self.ntc_niterations)  
  
        vbox2.addLayout(hbox22) 
          
                
        hbox23 = QtGui.QHBoxLayout()
        tc1 = QtGui.QLabel(self)
        tc1.setText("CS Parameter Beta")  
        hbox23.addWidget(tc1)      
        self.ntc_beta = QtGui.QLineEdit(self)
        self.ntc_beta.setFixedWidth(65)
        self.ntc_beta.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.ntc_beta.setAlignment(QtCore.Qt.AlignRight)         
        hbox23.addStretch(1)
        self.ntc_beta.setText(str(self.beta))
        hbox23.addWidget(self.ntc_beta) 
        
        vbox2.addLayout(hbox23) 
        
        
        hbox24 = QtGui.QHBoxLayout()
        text1 = QtGui.QLabel(self)
        text1.setText('Sample Thickness')
        hbox24.addWidget(text1)
        hbox24.addStretch(1)
        self.ntc_samplethick = QtGui.QLineEdit(self)
        self.ntc_samplethick.setFixedWidth(65)
        self.ntc_samplethick.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.ntc_samplethick.setAlignment(QtCore.Qt.AlignRight)   
        self.ntc_samplethick.setText(str(self.samplethick))
        hbox24.addWidget(self.ntc_samplethick)  
  
        vbox2.addLayout(hbox24) 
      
        
        sizer2.setLayout(vbox2)
        
        
        sizer1.setLayout(vbox1)

        #panel 3
        sizer3 = QtGui.QGroupBox('ROI')
        vbox3 = QtGui.QVBoxLayout()

        
        self.button_roi = QtGui.QPushButton( 'Select ROI')
        self.button_roi.clicked.connect( self.OnSelectROI)
        self.button_roi.setEnabled(False)
        vbox3.addWidget(self.button_roi)
        
        
        self.button_roispec = QtGui.QPushButton( 'Show ROI Spectrum')
        self.button_roispec.clicked.connect( self.OnShowROISpec)
        self.button_roispec.setEnabled(False)
        vbox3.addWidget(self.button_roispec)
        
        self.button_roidel = QtGui.QPushButton( 'Reset ROI')
        self.button_roidel.clicked.connect( self.OnResetROI)
        self.button_roidel.setEnabled(False)
        vbox3.addWidget(self.button_roidel)
        
        sizer3.setLayout(vbox3)


#         #panel 3
#         sizer3 = QtGui.QGroupBox('File')
#         vbox3 = QtGui.QVBoxLayout()
#  
#   
#         self.tc_file = QtGui.QLabel(self)
#         vbox3.addWidget(self.tc_file)
#         self.tc_file.setText('File name')
#         
#         vbox3.setContentsMargins(20,20,20,30)
#         sizer3.setLayout(vbox3)
#         
# 
#         #panel 4
#         sizer4 = QtGui.QGroupBox('Path')
#         vbox4 = QtGui.QVBoxLayout()
#   
#         self.tc_path = QtGui.QLabel(self)
#         vbox4.addWidget(self.tc_path)
#         self.tc_path.setText('D:/')
#        
#         vbox4.setContentsMargins(20,20,20,30)
#         sizer4.setLayout(vbox4)
                
 
        #panel 5     
        vbox5 = QtGui.QVBoxLayout()
        
        self.tc_imagecomp = QtGui.QLabel(self)
        self.tc_imagecomp.setText("Dataset: ")
        vbox5.addWidget(self.tc_imagecomp)
        

        gridsizertop = QtGui.QGridLayout()
        
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
   
        self.absimgfig = Figure((PlotH*1.25, PlotH*1.25))

        self.AbsImagePanel = FigureCanvas(self.absimgfig)
        self.AbsImagePanel.setParent(self)
        self.cid1 = self.AbsImagePanel.mpl_connect('button_press_event', self.OnPointImage)
        
        
        fbox.addWidget(self.AbsImagePanel)
        frame.setLayout(fbox)
        gridsizertop.addWidget(frame, 0, 0, QtCore .Qt. AlignLeft)
        

        self.slider_slice = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slider_slice.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_slice.valueChanged[int].connect(self.OnScrollSlice)
        self.slider_slice.setRange(0, 100)
        
        gridsizertop.addWidget(self.slider_slice, 0, 1, QtCore .Qt. AlignLeft)
        
        
        self.slider_comp = QtGui.QScrollBar(QtCore.Qt.Horizontal)
        self.slider_comp.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_comp.valueChanged[int].connect(self.OnScrollComp)
        self.slider_comp.setRange(0, 100)    
        self.slider_comp.setEnabled(True) 
        self.tc_comp = QtGui.QLabel(self)
        self.tc_comp.setText("Component: ")
        hbox51 = QtGui.QHBoxLayout()
        hbox51.addWidget(self.tc_comp) 
        hbox51.addWidget(self.slider_comp)
        gridsizertop.addLayout(hbox51, 1, 0)
        

        vbox5.addLayout(gridsizertop)
        vbox5.addStretch(1)
        
       
        vboxtop = QtGui.QVBoxLayout()
        
        hboxtop = QtGui.QHBoxLayout()
        vboxt1 = QtGui.QVBoxLayout()
        vboxt1.addWidget(sizer1)
        vboxt1.addStretch (1)
        vboxt1.addWidget(sizer2)
        vboxt1.addStretch (1)
        vboxt1.addWidget(sizer3)
        
        hboxtop.addStretch (0.5)
        hboxtop.addLayout(vboxt1)
        hboxtop.addStretch (0.5)
        hboxtop.addLayout(vbox5)
        hboxtop.addStretch (0.5)
        
        vboxtop.addStretch (0.5)
        vboxtop.addLayout(hboxtop)
        vboxtop.addStretch (0.9)
#         vboxtop.addWidget(sizer3)
#         vboxtop.addStretch (0.5)
#         vboxtop.addWidget(sizer4)
#         vboxtop.addStretch (0.5)

        vboxtop.setContentsMargins(50,50,50,50)
        self.setLayout(vboxtop)
        
        
        
#----------------------------------------------------------------------          
    def OnLoadTomoEng(self, event):
        
        self.fulltomorecdata = []        
        self.tomo_calculated = 0
        self.full_tomo_calculated = 0

        fig = self.absimgfig
        fig.clf()
        self.AbsImagePanel.draw()   
        
        self.button_save.setEnabled(False)
        self.button_save.setEnabled(False)
        self.tc_comp.setText('Component: ')
        self.slider_comp.setEnabled(False)
        
        
        self.tomodata = self.stack.od4D
        
        for i in range(self.stack.n_ev):
            self.datanames.append(str(self.stack.ev[i]))
            
        self.ncomponents = self.stack.n_ev
        
        self.combonames.clear()
        self.combonames.addItems(self.datanames)
        
        self.tc_imagecomp.setText("Dataset: Energies")
        
        self.button_calcall.setEnabled(True)
        self.button_calc1.setEnabled(True)
        self.energiesloaded = 1
        
        
#----------------------------------------------------------------------          
    def OnLoadTomoComponents(self, event):
        
        self.fulltomorecdata = []        
        self.tomo_calculated = 0
        self.full_tomo_calculated = 0

        fig = self.absimgfig
        fig.clf()
        self.AbsImagePanel.draw()   
        
        self.button_save.setEnabled(False)
        self.button_save.setEnabled(False)
        self.tc_comp.setText('Component: ')
        self.slider_comp.setEnabled(False)
        
        self.ncomponents = self.anlz.n_target_spectra
        
        self.tomodata = np.zeros((self.stack.n_cols, self.stack.n_rows, self.ncomponents, self.stack.n_theta))
        
        
        for i in range (self.stack.n_theta):   
            if self.window().page4.showraw == True:
                self.tomodata[:,:,:,i] = self.anlz.target_svd_maps4D[i]
            else:
                self.tomodata[:,:,:,i] = self.anlz.target_pcafit_maps[i]
         
       
        for i in range(self.ncomponents):
            self.datanames.append(str(self.anlz.tspec_names[i]))
            
        
        self.combonames.clear()
        self.combonames.addItems(self.datanames)
        
        self.tc_imagecomp.setText("Dataset: Spectral Components")
        
        self.button_calcall.setEnabled(True)
        self.button_calc1.setEnabled(True)
        self.energiesloaded = 0
        
        
#----------------------------------------------------------------------          
    def OnCalcTomoFull(self, event):
        
        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        
        
        value = self.ntc_niterations.text()
        self.maxIters = int(value)  
        value = self.ntc_beta.text()      
        self.beta = float(value)        
        value = self.ntc_samplethick.text()
        self.samplethick = int(value) 
        
        self.fulltomorecdata = []

        for i in range(self.ncomponents):
            
            print 'Progress ',i+1,' / ',self.ncomponents
                
            self.tr.calc_tomo1(self.tomodata[:,:,i,:], 
                               self.stack.theta,
                               self.maxIters,
                               self.beta,
                               self.samplethick)
            
            self.fulltomorecdata.append(self.tr.tomorec.copy())
            
        self.tr.tomorec = self.fulltomorecdata[self.icomp]
        
        
        self.full_tomo_calculated = 1
        self.tomo_calculated = 1
        
        dims = self.tr.tomorec.shape
        
        self.nslices = dims[2]
        
        self.islice = int(dims[2]/2)
        
        self.slider_slice.setRange(0, dims[2]-1)
        
        self.ROIvol =  [[]]* dims[2]
        
        self.tc_comp.setText('Component: '+self.datanames[self.icomp])
        
        self.slider_slice.setValue(self.islice)
        self.slider_comp.setRange(0, self.ncomponents-1)
        self.slider_comp.setEnabled(True)
        self.slider_comp.setValue(self.icomp)
        self.button_save.setEnabled(True)
        self.button_saveall.setEnabled(True)
        if self.energiesloaded == 1:
            self.button_roi.setEnabled(True)
        else:
            self.button_roi.setEnabled(False)
          
        
        self.ShowImage()
        
        QtGui.QApplication.restoreOverrideCursor()    
        
#----------------------------------------------------------------------          
    def OnCalcTomo1(self, event):
        
        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        
        
        value = self.ntc_niterations.text()
        self.maxIters = int(value)  
        value = self.ntc_beta.text()      
        self.beta = float(value)        
        value = self.ntc_samplethick.text()
        self.samplethick = int(value) 

                
        self.tr.calc_tomo1(self.tomodata[:,:,self.select1,:], 
                           self.stack.theta,
                           self.maxIters,
                           self.beta,
                           self.samplethick)
        
        self.tomo_calculated = 1
        self.full_tomo_calculated = 0
        
        dims = self.tr.tomorec.shape
        
        self.ROIvol = [[]]* dims[2]
        self.button_roispec.setEnabled(False)
        
        self.islice = int(dims[2]/2)
        self.slider_slice.setValue(self.islice)
        self.slider_slice.setRange(0, dims[2]-1)
        
        self.tc_comp.setText('Component: '+self.datanames[self.select1])
        self.slider_comp.setRange(0, 0)
        self.button_save.setEnabled(True)
        self.button_roi.setEnabled(True)
        
        self.ShowImage()
        
        QtGui.QApplication.restoreOverrideCursor()    
        
       
#----------------------------------------------------------------------           
    def OnSelect1Comp(self, value):
        item = value
        self.select1 = item
        
        
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
        iy = self.stack.n_rows-1-int(np.floor(y))  
                
        if ix<0 :
            ix=0
        if ix>self.stack.n_cols-1 :
            ix=self.stack.n_cols-1
        if iy<0 :
            iy=0
        if iy>self.stack.n_rows-1 :
            iy=self.stack.n_rows-1
            
        spectrum = []

        for i in range(self.ncomponents):
            spectrum.append(self.fulltomorecdata[i][ix,iy,self.islice])
            
        title = 'Point [{0:d}, {0:d}, {0:d}]'.format(ix,iy,self.islice)
                
        plot = PlotFrame(self, self.stack.ev, spectrum, title=title)
        plot.show()       

#----------------------------------------------------------------------    
    def OnSave(self, event): 
        
        
        wildcard = "Mrc files (*.mrc);;"

        SaveFileName = QtGui.QFileDialog.getSaveFileName(self, 'Save Tomo Reconstructions', '', wildcard)

        SaveFileName = str(SaveFileName)
        if SaveFileName == '':
            return
        
        
        data = self.tr.tomorec
                
        self.tr.save_mrc(SaveFileName, data)        
        
        
#----------------------------------------------------------------------    
    def OnSaveAll(self, event): 
        
        
        wildcard = "Mrc files (*.mrc);;"

        SaveFileName = QtGui.QFileDialog.getSaveFileName(self, 'Save Tomo Reconstructions', '', wildcard)

        SaveFileName = str(SaveFileName)
        if SaveFileName == '':
            return
        
        
        basename, extension = os.path.splitext(SaveFileName)  

        for i in range(self.ncomponents):
            data = self.fulltomorecdata[i]        
        
            savefn = basename + '_'+self.datanames[i]+extension
                
            self.tr.save_mrc(savefn, data)        

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

        ROIpix = np.zeros((self.stack.n_cols,self.stack.n_rows))    
        
        for i in range(self.stack.n_cols):
            for j in range(self.stack.n_rows):
                Pinside = path.contains_point((i,j))
                if Pinside == True:
                    ROIpix[i, self.stack.n_rows-1-j] = 255
              


        ROIpix = np.ma.array(ROIpix)
        
        ROIpix_masked =  np.ma.masked_values(ROIpix, 0)
        
        
        self.ROIvol[self.islice] = ROIpix_masked
        
        self.haveROI = 1
        
        self.ShowImage()
        


#----------------------------------------------------------------------           
    def OnSelectROI(self, event):
        
        self.AbsImagePanel.mpl_disconnect(self.cid1)
        self.button_roispec.setEnabled(True)
        self.button_roidel.setEnabled(True)
        
        lineprops = dict(color='red', linestyle='-', linewidth = 1, alpha=1)


        self.lasso = LassoSelector(self.axes, onselect=self.OnSelectLasso, useblit=False, lineprops=lineprops)
        
        
#----------------------------------------------------------------------           
    def OnResetROI(self, event):        
        
        self.button_roispec.setEnabled(False)
        self.button_roidel.setEnabled(False)
        
        self.ROIvol =  [[]]*self.nslices
        
        self.haveROI = 0
               
        self.cid1 = self.AbsImagePanel.mpl_connect('button_press_event', self.OnPointImage)
        
        self.ShowImage()

#----------------------------------------------------------------------    
    def CalcROISpectrum(self):
              
        ROIspectrum = np.zeros((self.stack.n_ev))
               

        for i in range(self.nslices):
            roimask = self.ROIvol[i]
            if roimask != []:
                for ie in range(self.stack.n_ev):
                    indices = np.where(roimask == 255)
                    numroipix = roimask[indices].shape[0]
                
                    roivoxels = self.fulltomorecdata[ie][:,:,i]
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
            return
        
        image = self.tr.tomorec[:,:,self.islice].copy() 

        fig = self.absimgfig
        fig.clf()
        fig.add_axes(((0.0,0.0,1.0,1.0)))
        axes = fig.gca()
        fig.patch.set_alpha(1.0)
         
        im = axes.imshow(np.rot90(image), cmap=matplotlib.cm.get_cmap("gray")) 
        
        if self.haveROI == 1:
            if self.ROIvol[self.islice] != []:
                im_red = axes.imshow(np.rot90(self.ROIvol[self.islice]), cmap=matplotlib.cm.get_cmap("autumn")) 
         
#         if self.window().page1.show_scale_bar == 1:
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
        
        self.xys = np.dstack(np.meshgrid(np.arange(self.stack.n_cols), np.arange(self.stack.n_rows))).reshape(-1,2)
         
        
#----------------------------------------------------------------------        
    def NewStackClear(self):
        
        self.tr = tomo_reconstruction.Ctomo(self.stack)
        
        self.fulltomorecdata = []
        self.ncomponents = 0
        
        self.tomo_calculated = 0
        self.full_tomo_calculated = 0
        self.energiesloaded = 0
        
        self.datanames = []     
        
        self.ROIvol = []
        self.haveROI = 0
        
        
        fig = self.absimgfig
        fig.clf()
        self.AbsImagePanel.draw()   
        
        self.button_spcomp.setEnabled(False)
        self.button_engdata.setEnabled(False)
        self.button_calc1.setEnabled(False)
        self.button_calcall.setEnabled(False)
        self.button_save.setEnabled(False)
        self.button_save.setEnabled(False)
        self.button_roi.setEnabled(False)
        
        self.tc_imagecomp.setText("Dataset: ")
        self.tc_comp.setText('Component: ')
        
        self.slider_comp.setEnabled(False)
        
        self.combonames.clear()
        
        

#----------------------------------------------------------------------
class File_GUI():
    """
    Ask user to choose file and then use an appropriate plugin to read and return a data structure.
    """
    def __init__(self):
        self.last_path = dict([a,dict([t,os.getcwd()] for t in file_plugins.data_types)] for a in file_plugins.actions)
        self.last_filter = dict([a,dict([t,0] for t in file_plugins.data_types)] for a in file_plugins.actions)
        self.supported_filters = file_plugins.supported_filters
        self.filter_list = file_plugins.filter_list
    
    def SelectFile(self,action,data_type):
        dlg=QtGui.QFileDialog(None)
        dlg.setWindowTitle('Choose File')
        dlg.setViewMode(QtGui.QFileDialog.Detail)
        if action == "write":
            dlg.setAcceptMode(QtGui.QFileDialog.AcceptSave)
        dlg.setDirectory(self.last_path[action][data_type])
        dlg.setFilters(self.filter_list[action][data_type])
        dlg.selectFilter(self.filter_list[action][data_type][self.last_filter[action][data_type]])
        if dlg.exec_(): #if not cancelled
            self.last_path[action][data_type] = os.path.split(str(dlg.selectedFiles()[0]))[0]
            chosen_plugin = None
            for i,filt in enumerate(self.filter_list[action][data_type][1:-1]):
                if filt==dlg.selectedFilter():
                    chosen_plugin = file_plugins.plugins[i]
                    break
            if chosen_plugin is not None:
                self.last_filter[action][data_type] = i+1
            return (str(dlg.selectedFiles()[0]),chosen_plugin)
        else:
            print "cancelled"
            return (None,None)




#----------------------------------------------------------------------
    class DataChoiceDialog(QtGui.QDialog):

        def __init__(self,filepath=None,filestruct=None):
            super(File_GUI.DataChoiceDialog, self).__init__()
            self.filepath = filepath
            self.selection = None
            self.setWindowTitle('Choose Dataset')
            
            # A vertical box layout containing rows
            self.MainSizer = QtGui.QVBoxLayout()
            hbox1 = QtGui.QHBoxLayout()
            hbox2 = QtGui.QHBoxLayout()
            hbox3 = QtGui.QHBoxLayout()
            hbox4 = QtGui.QHBoxLayout()
            
            # Add the widgets to the first row
            hbox1.addWidget(QtGui.QLabel("Path:"))
            self.Path_text = QtGui.QLabel("")
            hbox1.addWidget(self.Path_text,stretch=1)
            self.MainSizer.addLayout(hbox1)
            
            # Add the widgets to the second row
            hbox2.addWidget(QtGui.QLabel('File:'))
            self.File_text = QtGui.QLineEdit(self)
            hbox2.addWidget(self.File_text,stretch=1)
            browse_button = QtGui.QPushButton('Browse...')
            browse_button.clicked.connect(self.OnBrowse)
            hbox2.addWidget(browse_button)
            self.MainSizer.addLayout(hbox2)
            
            # Add the widgets to the third row - dynamic set of widgets to display file info
            self.Entry_info = File_GUI.EntryInfoBox()
            self.MainSizer.addWidget(self.Entry_info)
            
            # Add widgets for the fourth row - just OK, cancel buttons
            hbox4.addStretch(1)
            self.button_ok = QtGui.QPushButton('Accept')
            self.button_ok.clicked.connect(self.OnAccept)
            self.button_ok.setEnabled(False)
            hbox4.addWidget(self.button_ok)
            button_cancel = QtGui.QPushButton('Cancel')
            button_cancel.clicked.connect(self.close)
            hbox4.addWidget(button_cancel)
            self.MainSizer.addLayout(hbox4)
            
            # Set self.MainSizer as the layout for the window
            self.setLayout(self.MainSizer)
            self.setModal(True)
            #self.show()
            
            if filepath is not None:
                path, filename = os.path.split(str(self.filepath))
                self.Path_text.setText(path)
                self.File_text.setText(filename)
                if filestruct is None:
                    self.contents = file_plugins.GetFileStructure(str(self.filepath))
                else:
                    self.contents = filestruct
                self.Entry_info.UpdateInfo(self.contents)


        def OnAccept(self):
            self.selection = self.Entry_info.GetSelection()
            if len(self.selection) == 0:
                self.close() #accepting with no regions selected is the same as clicking 'cancel'
            else:
                self.accept()

        def OnBrowse(self):
            filepath, plugin = File_GUI.SelectFile('read','stack')
            if filepath is None:
                return
            path, filename = os.path.split(str(filepath))
            self.filepath = str(filepath)
            self.Path_text.setText(path)
            self.File_text.setText(filename)
            self.contents = file_plugins.GetFileStructure(str(filepath),plugin=plugin)
            if self.contents is None:
                self.selection = [(0,0)]
                self.accept()
            else:
                self.Entry_info.UpdateInfo(self.contents)
            

    #-------------------------------------
    class InfoItem(QtGui.QHBoxLayout):
        '''Row of widgets describing the contents of a file entry that appear within the EntryInfoBox'''
        def __init__(self,name,contents,allowed_application_definitions,index):
            super(File_GUI.InfoItem, self).__init__()
            self.valid_flag = self.checkValidity(contents,allowed_application_definitions)
            self.index = index
            self.populate(name,contents)
            
        def populate(self,name,contents):
            self.checkbox = QtGui.QCheckBox(name+' ('+str(contents.definition)+':'+str(contents.scan_type)+')')
            self.checkbox.clicked.connect(self.setChecked)
            self.addWidget(self.checkbox)
            self.addStretch(1)
            self.addWidget(QtGui.QLabel('Channels:'))
            self.channel_combobox = QtGui.QComboBox()
            for c in contents:
                self.channel_combobox.addItem(c)
            self.addWidget(self.channel_combobox)
            self.addStretch(1)
            Data_Size_Label = QtGui.QLabel('Points: '+str(contents.data_shape))
            if contents.data_axes is not None:
                Data_Size_Label.setToolTip(' | '.join(contents.data_axes))
            self.addWidget(Data_Size_Label)
            if not self.valid_flag:
                self.checkbox.setEnabled(False)
                self.channel_combobox.setEnabled(False)
        
        def setChecked(self,value=True):
            self.checkbox.setChecked(value)

        def isEnabled(self):
            return self.checkbox.isEnabled()

        def isValid(self):
            return self.valid_flag

        def checkValidity(self,contents,allowed_application_definitions):
            '''Check if data file entry appears to have the correct structure.'''
            return True
            #if contents.definition in allowed_application_definitions and contents.scan_type is not None and contents.data_shape is not None and contents.data_axes is not None:
                #return True
            #else:
                #return False
        
        def GetStatus(self):
            return (self.checkbox.isChecked(),self.channel_combobox.currentIndex())

    #---------------------------------------
    class EntryInfoBox(QtGui.QGroupBox):
        '''Widgets giving a summary of the contents of a file via an InfoItem object per file entry.'''
        def __init__(self):
            super(File_GUI.EntryInfoBox, self).__init__('File Summary')
            self.vbox = QtGui.QVBoxLayout()
            self.setLayout(self.vbox)
            self.valid_entry_flag = False
            self.selection = (0,0)
            
        def ClearAll(self):
            # Cycle through children and mark for deletion
            while self.vbox.count():
                row = self.vbox.takeAt(0)
                while row.count():
                    item = row.takeAt(0)
                    if item.widget() is not None:
                        item.widget().deleteLater()
            self.valid_entry_flag = False
            
        def UpdateInfo(self, FileContents):
            self.ClearAll()
            # Now popoulate widgets representing file contents
            for i,entry in enumerate(FileContents):
                entryCheckBox = File_GUI.InfoItem(entry,FileContents[entry],['NXstxm'],i)
                self.vbox.addLayout(entryCheckBox)
                if self.valid_entry_flag is False and entryCheckBox.isValid() is True:
                    entryCheckBox.setChecked(True)
                    self.valid_entry_flag = True
                    self.parent().button_ok.setEnabled(True)
        
        def GetSelection(self):
            selection = []
            for i in range(self.vbox.count()):
                status = self.vbox.itemAt(i).GetStatus()
                if status[0]:
                    selection.append((i,status[1])) # (region,detector)
            return selection

File_GUI = File_GUI() #Create instance so that object can remember things (e.g. last path)

""" ------------------------------------------------------------------------------------------------"""
class PageNNMA(QtGui.QWidget):
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
        sizer1 = QtGui.QGroupBox('NNMA analysis')
        hbox1 = QtGui.QHBoxLayout()
        vbox1 = QtGui.QVBoxLayout()
       
        vbox1.addStretch(1) 
        self.button_calcnnma = QtGui.QPushButton('Calculate NNMA')
        self.button_calcnnma.clicked.connect( self.OnCalcNNMA)   
        self.button_calcnnma.setEnabled(False)
        vbox1.addWidget(self.button_calcnnma)
        self.button_mucluster = QtGui.QPushButton('Load initial cluster spectra')
        self.button_mucluster.clicked.connect( self.OnLoadClusterSpectra)
        self.button_mucluster.setEnabled(False)
        self.button_mufile = QtGui.QPushButton('Load initial standard spectra')
        self.button_mufile.clicked.connect( self.OnLoadStandardSpectra)
        self.button_mufile.setEnabled(False)   
        self.button_murand = QtGui.QPushButton('Load initial random spectra')
        self.button_murand.clicked.connect( self.OnLoadRandomSpectra)         
        self.button_murand.setEnabled(False) 
        self.tc_initspectra = QtGui.QLabel(self)
        self.tc_initspectra.setText('Initial Spectra: ' + self.initMatrices)
        self.button_savennma = QtGui.QPushButton('Save NNMA Results...')
        self.button_savennma.clicked.connect( self.OnSave)
        self.button_savennma.setEnabled(False)
        vbox1.addWidget(self.button_savennma)
        
        vbox1.addStretch(1)
        
        
        sizer11 = QtGui.QGroupBox('NNMA settings')
        vbox11 = QtGui.QVBoxLayout()
                
        hbox11 = QtGui.QHBoxLayout()
        text1 = QtGui.QLabel(self)
        text1.setText('Number of components k')
        hbox11.addWidget(text1)
        hbox11.addStretch(1)
        self.ncompspin = QtGui.QSpinBox()
        self.ncompspin.setRange(2,20)
        self.ncompspin.setValue(self.nComponents)
        self.ncompspin.valueChanged[int].connect(self.OnNNMAspin)
        hbox11.addWidget(self.ncompspin)  
        
        vbox11.addLayout(hbox11) 
        
        line = QtGui.QFrame()
        line.setFrameShape(QtGui.QFrame.HLine)
        line.setFrameShadow(QtGui.QFrame.Sunken) 
        vbox11.addStretch(1)
        vbox11.addWidget(line) 
        vbox11.addStretch(1)    
        
        text1 = QtGui.QLabel(self)
        text1.setText('Regularization parameters')
        vbox11.addWidget(text1)   
        
        hbox15a = QtGui.QHBoxLayout()
        tc1 = QtGui.QLabel(self)
        tc1.setText("\tSpectra similarity")  
        hbox15a.addWidget(tc1)      
        self.ntc_lamsim = QtGui.QLineEdit(self)
        self.ntc_lamsim.setFixedWidth(65)
        self.ntc_lamsim.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.ntc_lamsim.setAlignment(QtCore.Qt.AlignRight)         
        hbox15a.addStretch(1)
        self.ntc_lamsim.setText(str(self.lambdaClusterSim))
        hbox15a.addWidget(self.ntc_lamsim) 
        
        vbox11.addLayout(hbox15a) 
        
        hbox15b = QtGui.QHBoxLayout()
        tc1 = QtGui.QLabel(self)
        tc1.setText("\tSpectra smoothness")  
        hbox15b.addWidget(tc1)      
        self.ntc_lamsmooth = QtGui.QLineEdit(self)
        self.ntc_lamsmooth.setFixedWidth(65)
        self.ntc_lamsmooth.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.ntc_lamsmooth.setAlignment(QtCore.Qt.AlignRight)         
        hbox15b.addStretch(1)
        self.ntc_lamsmooth.setText(str(self.lambdaSmooth))
        hbox15b.addWidget(self.ntc_lamsmooth) 
        
        vbox11.addLayout(hbox15b) 
        
        hbox15c = QtGui.QHBoxLayout()
        tc1 = QtGui.QLabel(self)
        tc1.setText("\tSparseness")  
        hbox15c.addWidget(tc1)      
        self.ntc_lamspar = QtGui.QLineEdit(self)
        self.ntc_lamspar.setFixedWidth(65)
        self.ntc_lamspar.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.ntc_lamspar.setAlignment(QtCore.Qt.AlignRight)         
        hbox15c.addStretch(1)
        self.ntc_lamspar.setText(str(self.lambdaSparse))
        hbox15c.addWidget(self.ntc_lamspar) 
        
        vbox11.addLayout(hbox15c) 
        
        line = QtGui.QFrame()
        line.setFrameShape(QtGui.QFrame.HLine)
        line.setFrameShadow(QtGui.QFrame.Sunken) 
        vbox11.addStretch(1)
        vbox11.addWidget(line) 
        vbox11.addStretch(1)           
         
        hbox12 = QtGui.QHBoxLayout()
        text1 = QtGui.QLabel(self)
        text1.setText('Number of iterations')
        hbox12.addWidget(text1)
        hbox12.addStretch(1)
        self.ntc_niterations = QtGui.QLineEdit(self)
        self.ntc_niterations.setFixedWidth(65)
        self.ntc_niterations.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.ntc_niterations.setAlignment(QtCore.Qt.AlignRight)   
        self.ntc_niterations.setText(str(self.maxIters))
        hbox12.addWidget(self.ntc_niterations)  
  
        vbox11.addLayout(hbox12) 
          
                
        hbox14a = QtGui.QHBoxLayout()
        tc1 = QtGui.QLabel(self)
        tc1.setText("Delta Error Threshold")  
        hbox14a.addWidget(tc1)      
        self.ntc_dthresh = QtGui.QLineEdit(self)
        self.ntc_dthresh.setFixedWidth(65)
        self.ntc_dthresh.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.ntc_dthresh.setAlignment(QtCore.Qt.AlignRight)         
        hbox14a.addStretch(1)
        self.ntc_dthresh.setText(str(self.deltaErrorThreshold))
        hbox14a.addWidget(self.ntc_dthresh) 
        
        vbox11.addLayout(hbox14a) 
        sizer11.setLayout(vbox11)

        line = QtGui.QFrame()
        line.setFrameShape(QtGui.QFrame.HLine)
        line.setFrameShadow(QtGui.QFrame.Sunken) 
        

        vbox1.addStretch(1)
        vbox1.addWidget(line) 
        vbox1.addStretch(1)     
        vbox1.addWidget(self.button_mucluster)
        vbox1.addWidget(self.button_mufile)
        vbox1.addWidget(self.button_murand)
        vbox1.addWidget(self.tc_initspectra)

        
        hbox1.addLayout(vbox1)
        hbox1.addWidget(sizer11)
        
        sizer1.setLayout(hbox1)
       
       
        #panel 2       
        vbox2 = QtGui.QVBoxLayout()
        
        self.tc_NNMAcomp = QtGui.QLabel(self)
        self.tc_NNMAcomp.setText("NNMA thickness map ")
        vbox2.addWidget(self.tc_NNMAcomp)
        
        hbox21 = QtGui.QHBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
        self.nnmaimgfig = Figure((PlotH*1.10, PlotH))
        self.NNMAImagePan = FigureCanvas(self.nnmaimgfig)
        self.NNMAImagePan.setParent(self) 
        
        fbox.addWidget(self.NNMAImagePan)
        frame.setLayout(fbox)
        hbox21.addWidget(frame)                        
            
        self.slidershow = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slidershow.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slidershow.setRange(0,20)
        #self.slidershow.setEnabled(False)          
        self.slidershow.valueChanged[int].connect(self.OnScroll)

        
        hbox21.addWidget(self.slidershow)

        vbox2.addLayout(hbox21)
        
        
        #panel 3
        vbox3 = QtGui.QVBoxLayout()
    
        self.tc_nnmaspec = QtGui.QLabel(self)
        self.tc_nnmaspec.setText("NNMA spectrum ") 
        vbox3.addWidget(self.tc_nnmaspec)
                        
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()

        self.nnmaspecfig = Figure((PlotW, PlotH))
        self.NNMASpecPan = FigureCanvas(self.nnmaspecfig)
        self.NNMASpecPan.setParent(self)
       
        fbox.addWidget(self.NNMASpecPan)
        frame.setLayout(fbox)
        vbox3.addWidget(frame)   
        
        
        #panel 4
        vbox4 = QtGui.QVBoxLayout()  
         
        text4 = QtGui.QLabel(self)
        text4.setText("Cost function")        
        vbox4.addWidget(text4)

        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
       
        self.costffig = Figure((PlotW*0.7, PlotH*0.7))
        self.CostFPan = FigureCanvas(self.costffig)
        self.CostFPan.setParent(self)
        #self.CostFPan.mpl_connect('button_press_event', self.OnPointEvalsImage)
         
        fbox.addWidget(self.CostFPan)
        frame.setLayout(fbox)
        vbox4.addWidget(frame)       
        

        #panel 5
        sizer5 = QtGui.QGroupBox('Display')
        vbox5 = QtGui.QVBoxLayout()         
         
        hbox51 = QtGui.QHBoxLayout()
        self.cb_inputsp = QtGui.QCheckBox('Overlay input spectra', self)
        self.cb_inputsp.stateChanged.connect(self.ShowSpectrum)
        hbox51.addWidget(self.cb_inputsp)
        vbox5.addLayout(hbox51)         
        
        hbox52 = QtGui.QHBoxLayout()
        self.cb_showallsp = QtGui.QCheckBox('Show all spectra', self)
        self.cb_showallsp.stateChanged.connect(self.ShowSpectrum)
        hbox52.addWidget(self.cb_showallsp)
        vbox5.addLayout(hbox52)        
        
        line = QtGui.QFrame()
        line.setFrameShape(QtGui.QFrame.HLine)
        line.setFrameShadow(QtGui.QFrame.Sunken) 
        
        vbox5.addWidget(line)     
        
        hbox52 = QtGui.QHBoxLayout()
        self.cb_totalcost = QtGui.QCheckBox('Total cost', self)
        self.cb_totalcost.setChecked(True)
        self.cb_totalcost.stateChanged.connect(self.ShowCostFunction)
        hbox52.addWidget(self.cb_totalcost)
        vbox5.addLayout(hbox52)  
        
        hbox53 = QtGui.QHBoxLayout()
        self.cb_simlcost = QtGui.QCheckBox('Similarity cost', self)
        self.cb_simlcost.stateChanged.connect(self.ShowCostFunction)
        hbox53.addWidget(self.cb_simlcost)
        vbox5.addLayout(hbox53)  
        
        hbox54 = QtGui.QHBoxLayout()
        self.cb_smoothcost = QtGui.QCheckBox('Smoothness cost', self)
        self.cb_smoothcost.stateChanged.connect(self.ShowCostFunction)
        hbox54.addWidget(self.cb_smoothcost)
        vbox5.addLayout(hbox54)  

        hbox54 = QtGui.QHBoxLayout()
        self.cb_sparcost = QtGui.QCheckBox('Sparsness cost', self)
        self.cb_sparcost.stateChanged.connect(self.ShowCostFunction)
        hbox54.addWidget(self.cb_sparcost)
        vbox5.addLayout(hbox54)  
        
        sizer5.setLayout(vbox5)
         
        
       
        vboxtop = QtGui.QVBoxLayout()
        hboxtopL = QtGui.QHBoxLayout()
        vboxtopL = QtGui.QVBoxLayout()
        vboxtopL.addWidget(sizer1)
        hboxtopL.addLayout(vboxtopL)
        hboxtopL.addWidget(sizer5)
        hboxtopL.addStretch(1)
        
        gridsizertop = QtGui.QGridLayout()
        gridsizertop.setContentsMargins(15,0,0,0)

        
        gridsizertop.addLayout(hboxtopL, 0, 0, QtCore .Qt. AlignLeft)
        gridsizertop.addLayout(vbox4, 0, 1, QtCore .Qt. AlignLeft)
        
        gridsizertop.addLayout(vbox3, 1, 0, QtCore .Qt. AlignCenter)
        gridsizertop.addLayout(vbox2, 1, 1, QtCore .Qt. AlignLeft)
        
        vboxtop.addStretch(1)
        vboxtop.addLayout(gridsizertop)
        vboxtop.addStretch(1)
        self.setLayout(vboxtop)


#----------------------------------------------------------------------          
    def OnLoadClusterSpectra(self, event):
        
        self.nnma.setClusterSpectra(self.anlz.clusterspectra)
                                                        
        self.window().refresh_widgets()     
        
        self.initMatrices = 'Cluster'  
        
        self.tc_initspectra.setText('Initial Spectra: ' + self.initMatrices)   
        
#----------------------------------------------------------------------          
    def OnLoadRandomSpectra(self, event):
        
        self.initMatrices = 'Random'
        
        self.lambdaClusterSim = 0.0
        self.ntc_lamsim.setText(str(self.lambdaClusterSim))
        
        self.tc_initspectra.setText('Initial Spectra: ' + self.initMatrices)
          
#----------------------------------------------------------------------          
    def OnLoadStandardSpectra(self, event):    

        if True:
        #try: 
            
            wildcard = "Spectrum files (*.csv);;Spectrum files (*.xas);;"
            
            #filepath = QtGui.QFileDialog.getOpenFileNames(self, 'Choose Spectra files', '', wildcard, self.DefaultDir)
            filepath = QtGui.QFileDialog.getOpenFileNames(self, 'Choose Spectra files', '', wildcard)

            if filepath == '':
                return
            
            self.filenames = []
            for name in filepath:
                self.filenames.append(str(name))

                
            kNNMA = len(self.filenames)
            muInit = np.zeros((self.stk.n_ev, kNNMA))
                                   
            for i in range(kNNMA):
                thisfn = os.path.basename(str(self.filenames[i]))
                _basename, extension = os.path.splitext(thisfn)      
                if extension == '.csv':
                    spectrum_evdata, spectrum_data, _spectrum_common_name = self.stk.read_csv(self.filenames[i])
                elif extension == '.xas':
                    spectrum_evdata, spectrum_data, _spectrum_common_name = self.stk.read_xas(self.filenames[i])
                                        
                # Map this spectrum onto our energy range - interpolate to ev
                init_spectrum = np.interp(self.stk.ev, spectrum_evdata, spectrum_data)      

                muInit[:, i] = init_spectrum[:]
                
            
            self.nnma.setStandardsSpectra(muInit)   
            self.initMatrices = 'Standards'                     
            self.lambdaClusterSim = 0.0
            self.ntc_lamsim.setText(str(self.lambdaClusterSim))
                                    
                    
            QtGui.QApplication.restoreOverrideCursor()
            
#         except:
#             QtGui.QApplication.restoreOverrideCursor()  
#             QtGui.QMessageBox.warning(self, 'Error', 'Spectra files not loaded.')

        
            self.tc_initspectra.setText('Initial Spectra: ' + self.initMatrices)
                                   
                                 
        self.window().refresh_widgets()  
                      
#----------------------------------------------------------------------          
    def OnCalcNNMA(self, event):
        
        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))

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
        
        QtGui.QApplication.restoreOverrideCursor()    
        
        
        
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
                    fig.savefig(fileName_img, bbox_inches='tight', pad_inches = 0.0)
            
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
                    fig.savefig(fileName_img, bbox_inches='tight', pad_inches = 0.0)
            
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
                    fig.savefig(fileName_img, bbox_inches='tight', pad_inches = 0.0)
            
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
                                
                                
        except IOError, e:
            if e.strerror:
                err = e.strerror 
            else: 
                err = e 
    
            QtGui.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err)
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



#---------------------------------------------------------------------- 
class SaveWinP5(QtGui.QDialog):

    def __init__(self, parent):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        
        self.resize(400, 300)
        self.setWindowTitle('Save')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)


        self.com = self.parent.common          
        
        _path, ext = os.path.splitext(self.com.filename) 
        ext = ext[1:].lower()   
        suffix = "." + ext
        _path, fn = os.path.split(self.com.filename)
        filename = fn[:-len(suffix)]
        
        self.path = self.com.path
        self.filename = filename
                          
        
        vboxtop = QtGui.QVBoxLayout()
        vboxtop.setContentsMargins(20,20,20,20)
        
        gridtop = QtGui.QGridLayout()
        gridtop.setVerticalSpacing(20)
    
        fontb = QtGui.QFont()
        fontb.setBold(True)            
        
        st1 = QtGui.QLabel(self)
        st1.setText('Save')
        st1.setFont(fontb)
        st2 = QtGui.QLabel(self)
        st2.setText('.pdf')
        st2.setFont(fontb)
        st3 = QtGui.QLabel(self)
        st3.setText('.png')
        st3.setFont(fontb)
        st4 = QtGui.QLabel(self)
        st4.setText('.svg')
        st4.setFont(fontb)
        st5 = QtGui.QLabel(self)
        st5.setText('.csv')
        st5.setFont(fontb)        
        st9 = QtGui.QLabel(self)
        st9.setText('.tif (data)')
        st9.setFont(fontb)             
        
        st6 = QtGui.QLabel(self)
        st6.setText('_spectrum')
        
        self.cb11 = QtGui.QCheckBox('', self)
        self.cb11.setChecked(True)
        self.cb12 = QtGui.QCheckBox('', self)
        self.cb13 = QtGui.QCheckBox('', self)   
        self.cb14 = QtGui.QCheckBox('', self)
        
        st7 = QtGui.QLabel(self)
        st7.setText('_map')
        
        self.cb21 = QtGui.QCheckBox('', self)
        self.cb21.setChecked(True)
        self.cb22 = QtGui.QCheckBox('', self)
        self.cb23 = QtGui.QCheckBox('', self)
        self.cb24 = QtGui.QCheckBox('', self)
        
        st8 = QtGui.QLabel(self)
        st8.setText('_costfunction')   
        self.cb31 = QtGui.QCheckBox('', self)
        self.cb32 = QtGui.QCheckBox('', self)
        self.cb33 = QtGui.QCheckBox('', self)


        gridtop.addWidget(st1, 0, 0)
        gridtop.addWidget(st2, 0, 1)
        gridtop.addWidget(st3, 0, 2)
        gridtop.addWidget(st4, 0, 3)
        gridtop.addWidget(st5, 0, 4)
        gridtop.addWidget(st9, 0, 5)
                                
        gridtop.addWidget(st6, 1, 0)
        gridtop.addWidget(self.cb11, 1, 1)
        gridtop.addWidget(self.cb12, 1, 2)  
        gridtop.addWidget(self.cb13, 1, 3)  
        gridtop.addWidget(self.cb14, 1, 4)           
  
        gridtop.addWidget(st7, 2, 0)
        gridtop.addWidget(self.cb21, 2, 1)
        gridtop.addWidget(self.cb22, 2, 2) 
        gridtop.addWidget(self.cb23, 2, 3) 
        gridtop.addWidget(self.cb24, 2, 5)  
        
        gridtop.addWidget(st8, 3, 0)
        gridtop.addWidget(self.cb31, 3, 1)
        gridtop.addWidget(self.cb32, 3, 2) 
        gridtop.addWidget(self.cb33, 3, 3) 
        
        vboxtop.addStretch(0.5)
        vboxtop.addLayout(gridtop)
        vboxtop.addStretch(1)
        
         
        hbox0 = QtGui.QHBoxLayout()
         
        stf = QtGui.QLabel(self)
        stf.setText('Filename:\t')
        self.tc_savefn = QtGui.QLineEdit(self)
        self.tc_savefn.setText(self.filename)

        hbox0.addWidget(stf)
        hbox0.addWidget(self.tc_savefn)         
                  
        hbox1 = QtGui.QHBoxLayout()
                 
        stp = QtGui.QLabel(self)
        stp.setText('Path:  \t')
        self.tc_savepath = QtGui.QLineEdit(self)
        self.tc_savepath.setReadOnly(True)
        self.tc_savepath.setText(self.path)
        self.tc_savepath.setMinimumWidth(100)
        hbox1.addWidget(stp)
        hbox1.addWidget(self.tc_savepath)  
         
        button_path = QtGui.QPushButton('Browse...')
        button_path.clicked.connect(self.OnBrowseDir)
        hbox1.addWidget(button_path)
         
         
        hbox2 = QtGui.QHBoxLayout()
        button_save = QtGui.QPushButton('Save')
        button_save.clicked.connect(self.OnSave)
        hbox2.addWidget(button_save)
         
        button_cancel = QtGui.QPushButton('Cancel')
        button_cancel.clicked.connect(self.close)
        hbox2.addWidget(button_cancel)
        
        vboxtop.addLayout(hbox0)
        vboxtop.addLayout(hbox1)
        vboxtop.addStretch(1.0)
        vboxtop.addLayout(hbox2)

        
        
        self.setLayout(vboxtop)
        
#----------------------------------------------------------------------        
    def OnBrowseDir(self, evt):
        
        directory = QtGui.QFileDialog.getExistingDirectory(self, "Choose a directory", self.path, QtGui.QFileDialog.ShowDirsOnly|QtGui.QFileDialog.ReadOnly)       
                                                        
        
       
        if directory == '':
            return
                 
        directory = str(directory)
        self.com.path = directory
                    
        self.path = directory
        
        self.tc_savepath.setText(self.path)
            
            
                
#----------------------------------------------------------------------        
    def OnSave(self, evt):
        
        self.filename = str(self.tc_savefn.text())
        
        sp_pdf = self.cb11.isChecked()
        sp_png = self.cb12.isChecked()
        sp_svg = self.cb13.isChecked()
        sp_csv = self.cb14.isChecked()
        im_pdf = self.cb21.isChecked()
        im_png = self.cb22.isChecked()
        im_svg = self.cb23.isChecked()
        im_tif = self.cb24.isChecked()
        cf_pdf = self.cb31.isChecked()
        cf_png = self.cb32.isChecked()
        cf_svg = self.cb33.isChecked()
        
        self.close() 
        self.parent.page7.Save(self.filename, self.path,
                               spec_png = sp_png, 
                               spec_pdf = sp_pdf, 
                               spec_svg = sp_svg,
                               spec_csv = sp_csv,
                               map_png = im_png, 
                               map_pdf = im_pdf,
                               map_svg = im_svg,
                               map_tif = im_tif,
                               costf_png = cf_png, 
                               costf_pdf = cf_pdf,
                               costf_svg = cf_svg)        
        

                
""" ------------------------------------------------------------------------------------------------"""
class PageXrayPeakFitting(QtGui.QWidget):
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
        
        vboxL = QtGui.QVBoxLayout()
        vboxR = QtGui.QVBoxLayout()
        hboxT = QtGui.QHBoxLayout()
    

        #panel 
        sizer = QtGui.QGroupBox('Fit Parameters')
        vbox = QtGui.QVBoxLayout()
                
        font = QtGui.QFont()
        font.setPointSize(7)

        fgs1 = QtGui.QGridLayout()
        

        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Position [eV]')  
        textctrl.setFont(font)
        fgs1.addWidget(textctrl, 0, 1)      
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Height [OD]')  
        textctrl.setFont(font)
        fgs1.addWidget(textctrl, 0, 2)   
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Sigma [eV]')  
        textctrl.setFont(font)
        fgs1.addWidget(textctrl, 0, 3)   
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Base')  
        textctrl.setFont(font)
        fgs1.addWidget(textctrl, 16, 0)   


        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Step 1')
        fgs1.addWidget(textctrl, 1, 0)
        
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Step 2')
        fgs1.addWidget(textctrl, 2, 0)

        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Height [OD]')  
        textctrl.setFont(font)
        fgs1.addWidget(textctrl, 3, 2)      
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Position [eV]')  
        textctrl.setFont(font)
        fgs1.addWidget(textctrl, 3, 1)   
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Sigma [eV]')  
        textctrl.setFont(font)
        fgs1.addWidget(textctrl, 3, 3)   
         

        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 1')
        fgs1.addWidget(textctrl, 4, 0)
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 2')
        fgs1.addWidget(textctrl, 5, 0)
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 3')
        fgs1.addWidget(textctrl, 6, 0)
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 4')
        fgs1.addWidget(textctrl, 7, 0)        
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 5')
        fgs1.addWidget(textctrl, 8, 0)
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 6')
        fgs1.addWidget(textctrl, 9, 0)
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 7')
        fgs1.addWidget(textctrl, 10, 0)                
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 8')
        fgs1.addWidget(textctrl, 11, 0)    
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 9')
        fgs1.addWidget(textctrl, 12, 0)
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 10')
        fgs1.addWidget(textctrl, 13, 0)        
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 11')
        fgs1.addWidget(textctrl, 14, 0)                    
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 12')
        fgs1.addWidget(textctrl, 15, 0)
        
        
        self.le_sa1 = QtGui.QLineEdit(self)
        self.le_sa1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sa1.setAlignment(QtCore.Qt.AlignRight)
        self.le_sa1.setMaximumWidth(65)
        self.le_sa1.setText(str(self.stepfitparams[0]))
        fgs1.addWidget(self.le_sa1, 1, 1)

        self.le_sp1 = QtGui.QLineEdit(self)
        self.le_sp1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sp1.setAlignment(QtCore.Qt.AlignRight)
        self.le_sp1.setMaximumWidth(65)
        self.le_sp1.setText(str(self.stepfitparams[1]))
        fgs1.addWidget(self.le_sp1, 1, 2)
        
        self.le_sf1 = QtGui.QLineEdit(self)
        self.le_sf1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sf1.setAlignment(QtCore.Qt.AlignRight)
        self.le_sf1.setMaximumWidth(65)
        self.le_sf1.setText(str(self.stepfitparams[2]))
        fgs1.addWidget(self.le_sf1, 1, 3)
        
        self.le_base = QtGui.QLineEdit(self)
        self.le_base.setValidator(QtGui.QDoubleValidator(-99999, 99999, 2, self))
        self.le_base.setAlignment(QtCore.Qt.AlignRight)
        self.le_base.setMaximumWidth(65)
        self.le_base.setText(str(self.base))
        fgs1.addWidget(self.le_base, 16, 1)

        
        self.le_sa2 = QtGui.QLineEdit(self)
        self.le_sa2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sa2.setAlignment(QtCore.Qt.AlignRight)
        self.le_sa2.setMaximumWidth(65)
        self.le_sa2.setText(str(self.stepfitparams[3]))
        fgs1.addWidget(self.le_sa2, 2, 1)

        self.le_sp2 = QtGui.QLineEdit(self)
        self.le_sp2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sp2.setAlignment(QtCore.Qt.AlignRight)
        self.le_sp2.setMaximumWidth(65)
        self.le_sp2.setText(str(self.stepfitparams[4]))
        fgs1.addWidget(self.le_sp2, 2, 2)
        
        self.le_sf2 = QtGui.QLineEdit(self)
        self.le_sf2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sf2.setAlignment(QtCore.Qt.AlignRight)
        self.le_sf2.setMaximumWidth(65)
        self.le_sf2.setText(str(self.stepfitparams[5]))
        fgs1.addWidget(self.le_sf2, 2, 3)
        
                

        self.le_amplitudes = []
        self.le_mu = []
        self.le_sigma = []

        self.le_ga1 = QtGui.QLineEdit(self)
        self.le_ga1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga1.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga1.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga1, 4, 2)
        self.le_amplitudes.append(self.le_ga1)

        self.le_gm1 = QtGui.QLineEdit(self)
        self.le_gm1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm1.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm1.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm1, 4, 1)
        self.le_mu.append(self.le_gm1)
        
        self.le_gs1 = QtGui.QLineEdit(self)
        self.le_gs1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs1.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs1.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs1, 4, 3)
        self.le_sigma.append(self.le_gs1)
        
        self.le_ga2 = QtGui.QLineEdit(self)
        self.le_ga2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga2.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga2.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga2, 5, 2)
        self.le_amplitudes.append(self.le_ga2)

        self.le_gm2 = QtGui.QLineEdit(self)
        self.le_gm2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm2.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm2.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm2, 5, 1)
        self.le_mu.append(self.le_gm2)
        
        self.le_gs2 = QtGui.QLineEdit(self)
        self.le_gs2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs2.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs2.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs2, 5, 3)
        self.le_sigma.append(self.le_gs2)
        
        self.le_ga3 = QtGui.QLineEdit(self)
        self.le_ga3.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga3.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga3.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga3, 6, 2)
        self.le_amplitudes.append(self.le_ga3)

        self.le_gm3 = QtGui.QLineEdit(self)
        self.le_gm3.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm3.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm3.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm3, 6, 1)
        self.le_mu.append(self.le_gm3)
        
        self.le_gs3 = QtGui.QLineEdit(self)
        self.le_gs3.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs3.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs3.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs3, 6, 3)
        self.le_sigma.append(self.le_gs3)                                                

        self.le_ga4 = QtGui.QLineEdit(self)
        self.le_ga4.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga4.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga4.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga4, 7, 2)
        self.le_amplitudes.append(self.le_ga4)

        self.le_gm4 = QtGui.QLineEdit(self)
        self.le_gm4.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm4.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm4.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm4, 7, 1)
        self.le_mu.append(self.le_gm4)
        
        self.le_gs4 = QtGui.QLineEdit(self)
        self.le_gs4.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs4.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs4.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs4, 7, 3)
        self.le_sigma.append(self.le_gs4)
        
        self.le_ga5 = QtGui.QLineEdit(self)
        self.le_ga5.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga5.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga5.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga5, 8, 2)
        self.le_amplitudes.append(self.le_ga5)

        self.le_gm5 = QtGui.QLineEdit(self)
        self.le_gm5.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm5.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm5.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm5, 8, 1)
        self.le_mu.append(self.le_gm5)
        
        self.le_gs5 = QtGui.QLineEdit(self)
        self.le_gs5.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs5.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs5.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs5, 8, 3)
        self.le_sigma.append(self.le_gs5)
               
        self.le_ga6 = QtGui.QLineEdit(self)
        self.le_ga6.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga6.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga6.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga6, 9, 2)
        self.le_amplitudes.append(self.le_ga6)

        self.le_gm6 = QtGui.QLineEdit(self)
        self.le_gm6.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm6.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm6.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm6, 9, 1)
        self.le_mu.append(self.le_gm6)
        
        self.le_gs6 = QtGui.QLineEdit(self)
        self.le_gs6.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs6.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs6.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs6, 9, 3)
        self.le_sigma.append(self.le_gs6)
 
        self.le_ga7 = QtGui.QLineEdit(self)
        self.le_ga7.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga7.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga7.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga7, 10, 2)
        self.le_amplitudes.append(self.le_ga7)

        self.le_gm7 = QtGui.QLineEdit(self)
        self.le_gm7.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm7.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm7.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm7, 10, 1)
        self.le_mu.append(self.le_gm7)
        
        self.le_gs7 = QtGui.QLineEdit(self)
        self.le_gs7.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs7.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs7.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs7, 10, 3)
        self.le_sigma.append(self.le_gs7) 
       
        self.le_ga8 = QtGui.QLineEdit(self)
        self.le_ga8.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga8.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga8.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga8, 11, 2)
        self.le_amplitudes.append(self.le_ga8)

        self.le_gm8 = QtGui.QLineEdit(self)
        self.le_gm8.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm8.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm8.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm8, 11, 1)
        self.le_mu.append(self.le_gm8)
        
        self.le_gs8 = QtGui.QLineEdit(self)
        self.le_gs8.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs8.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs8.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs8, 11, 3)
        self.le_sigma.append(self.le_gs8) 
        
        self.le_ga9 = QtGui.QLineEdit(self)
        self.le_ga9.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga9.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga9.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga9, 12, 2)
        self.le_amplitudes.append(self.le_ga9)

        self.le_gm9 = QtGui.QLineEdit(self)
        self.le_gm9.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm9.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm9.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm9, 12, 1)
        self.le_mu.append(self.le_gm9)
        
        self.le_gs9 = QtGui.QLineEdit(self)
        self.le_gs9.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs9.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs9.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs9, 12, 3)
        self.le_sigma.append(self.le_gs9)     
        
        self.le_ga10 = QtGui.QLineEdit(self)
        self.le_ga10.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga10.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga10.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga10, 13, 2)
        self.le_amplitudes.append(self.le_ga10)

        self.le_gm10 = QtGui.QLineEdit(self)
        self.le_gm10.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm10.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm10.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm10, 13, 1)
        self.le_mu.append(self.le_gm10)
        
        self.le_gs10 = QtGui.QLineEdit(self)
        self.le_gs10.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs10.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs10.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs10, 13, 3)
        self.le_sigma.append(self.le_gs10)

        self.le_ga11 = QtGui.QLineEdit(self)
        self.le_ga11.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga11.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga11.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga11, 14, 2)
        self.le_amplitudes.append(self.le_ga11)

        self.le_gm11 = QtGui.QLineEdit(self)
        self.le_gm11.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm11.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm11.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm11, 14, 1)
        self.le_mu.append(self.le_gm11)
        
        self.le_gs11 = QtGui.QLineEdit(self)
        self.le_gs11.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gs11.setAlignment(QtCore.Qt.AlignRight)
        self.le_gs11.setMaximumWidth(65)
        fgs1.addWidget(self.le_gs11, 14, 3)
        self.le_sigma.append(self.le_gs11)
               
        self.le_ga12 = QtGui.QLineEdit(self)
        self.le_ga12.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_ga12.setAlignment(QtCore.Qt.AlignRight)
        self.le_ga12.setMaximumWidth(65)
        fgs1.addWidget(self.le_ga12, 15, 2)
        self.le_amplitudes.append(self.le_ga12)

        self.le_gm12 = QtGui.QLineEdit(self)
        self.le_gm12.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_gm12.setAlignment(QtCore.Qt.AlignRight)
        self.le_gm12.setMaximumWidth(65)
        fgs1.addWidget(self.le_gm12, 15, 1)
        self.le_mu.append(self.le_gm12)
        
        self.le_gs12 = QtGui.QLineEdit(self)
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
        vbox2 = QtGui.QVBoxLayout()
         
        self.tc_spec = QtGui.QLabel(self)
        self.tc_spec.setText("X-ray Spectrum: ")
        hbox11 = QtGui.QHBoxLayout()       
         
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
        self.Specfig = Figure((PlotW*1.25, PlotH*1.25))
        self.SpectrumPanel = FigureCanvas(self.Specfig)
        fbox.addWidget(self.SpectrumPanel)
        frame.setLayout(fbox)
        
        self.slider_spec = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slider_spec.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_spec.setEnabled(False)
        self.slider_spec.valueChanged[int].connect(self.OnSScroll)
        self.slider_spec.setRange(1, 5)
           
        hbox11.addWidget(frame)
        hbox11.addWidget(self.slider_spec)
           
        vbox2.addWidget(self.tc_spec)       
        vbox2.addLayout(hbox11)
         
         
         
        #panel 3
        sizer3 = QtGui.QGroupBox('Load Spectrum')
        vbox3 = QtGui.QVBoxLayout()    
        
        self.button_addclspec = QtGui.QPushButton('Add Cluster Spectra')
        self.button_addclspec.clicked.connect( self.OnAddClusterSpectra)   
        self.button_addclspec.setEnabled(False)
        vbox3.addWidget(self.button_addclspec)    
         
        self.button_loadtspec = QtGui.QPushButton('Load Spectrum')
        self.button_loadtspec.clicked.connect(self.OnSpecFromFile)
        vbox3.addWidget(self.button_loadtspec)            
 
        self.button_save = QtGui.QPushButton('Save Images...')
        self.button_save.clicked.connect( self.OnSave)
        self.button_save.setEnabled(False)
        vbox3.addWidget(self.button_save)

        sizer3.setLayout(vbox3)
         
 
         
        #panel 4
        sizer4 = QtGui.QGroupBox('Fit Settings')
        vbox4 = QtGui.QVBoxLayout()
                
           
        hbox41 = QtGui.QHBoxLayout()
        text1 = QtGui.QLabel(self)
        text1.setText('Number of steps')
        hbox41.addWidget(text1)
        self.nstepsspin = QtGui.QSpinBox()
        self.nstepsspin.setMinimumWidth(85)
        self.nstepsspin.setRange(0,2)
        self.nstepsspin.setValue(0)
        self.nstepsspin.valueChanged[int].connect(self.OnNstepsspin)
        hbox41.addWidget(text1)
        hbox41.addWidget(self.nstepsspin)  
        vbox4.addLayout(hbox41)
        
        hbox42 = QtGui.QHBoxLayout()
        text2 = QtGui.QLabel(self)
        text2.setText('Number of peaks')
        hbox42.addWidget(text2)
        self.ngaussspin = QtGui.QSpinBox()
        self.ngaussspin.setMinimumWidth(85)
        self.ngaussspin.setRange(0,12)
        self.ngaussspin.setValue(0)
        self.ngaussspin.valueChanged[int].connect(self.OnNgaussspin)
        hbox42.addWidget(text2)
        hbox42.addWidget(self.ngaussspin)  
        vbox4.addLayout(hbox42)

        
#         self.button_initfitparams = QtGui.QPushButton('Initialize Fit Parameters')
#         self.button_initfitparams.clicked.connect( self.OnInitFitParams)   
#         self.button_initfitparams.setEnabled(False)
#         vbox4.addWidget(self.button_initfitparams)
        
        self.button_fitspec = QtGui.QPushButton('Fit Spectrum')
        self.button_fitspec.clicked.connect( self.OnFitSpectrum)   
        self.button_fitspec.setEnabled(False)
        vbox4.addWidget(self.button_fitspec)
        
        sizer4.setLayout(vbox4)


        #panel 5
        
        sizer5 = QtGui.QGroupBox('Peak Positions [eV]')
        vbox5 = QtGui.QVBoxLayout()
                   
        self.tc_peakid = []
        for i in range(12):
            tc = QtGui.QLabel(self)
            self.tc_peakid.append(tc)
            vbox5.addWidget(tc)          
            
        vbox5.addStretch(1)
        
        self.cb_showengs = QtGui.QCheckBox('Show eV lines', self) 
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
        
        try: 
            
            wildcard = "Spectrum files (*.csv)"
            
            filepath = QtGui.QFileDialog.getOpenFileName(self, 'Choose Spectrum file', '', wildcard)
            

            filepath = str(filepath)
            if filepath == '':
                return
            
            self.filename =  os.path.basename(str(filepath))
            directory =  os.path.dirname(str(filepath))
            self.DefaultDir = directory            
                                                        
            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))    
                                            
            self.anlz.load_xraypeakfit_spectrum(filename=filepath)

            self.com.xpf_loaded = 1
            self.spectrumfitted.append(0)
            self.firstinit.append(1)
            self.initdone.append(0)
            self.nsteps.append(1)
            self.npeaks.append(4)
            self.fits.append(0)
            self.fits_sep.append(0)
            
            self.i_spec = self.anlz.n_xrayfitsp      
            self.slider_spec.setMaximum(self.anlz.n_xrayfitsp)
            self.slider_spec.setValue(self.i_spec)
            self.nstepsspin.setValue(self.nsteps[self.i_spec-1])
            self.ngaussspin.setValue(self.npeaks[self.i_spec-1])
            
            self.loadSpectrum()
            self.updatewidgets()
                    
            QtGui.QApplication.restoreOverrideCursor()
            
        except:
            QtGui.QApplication.restoreOverrideCursor()  
            QtGui.QMessageBox.warning(self, 'Error', 'Spectrum file not loaded.')
                                   
                                 
        
        self.window().refresh_widgets()


#----------------------------------------------------------------------
    def OnAddClusterSpectra(self, event):
        
        try:

            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor)) 
            

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
                    
            QtGui.QApplication.restoreOverrideCursor()
            
        except:
            QtGui.QApplication.restoreOverrideCursor()  
            QtGui.QMessageBox.warning(self, 'Error', 'Cluster spectra not loaded.')
                                   
                                 
        
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
        
        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        
        self.SetFitParams()
        
        fitted_spectrum, separate_peaks = self.anlz.fit_spectrum(self.i_spec-1, self.nsteps[self.i_spec-1], self.npeaks[self.i_spec-1])
        
        self.fits[self.i_spec-1] = fitted_spectrum
        self.fits_sep[self.i_spec-1] = separate_peaks
        
        self.spectrumfitted[self.i_spec-1] = 1
        
        self.loadSpectrum()
        
        QtGui.QApplication.restoreOverrideCursor() 
        self.updatewidgets()
        
       

#----------------------------------------------------------------------
    def OnSave(self, event):
        
        
        #Save images      
        wildcard = "Portable Network Graphics (*.png);;Adobe PDF Files (*.pdf);; SVG (*.svg)"

        fileName = QtGui.QFileDialog.getSaveFileName(self, 'Save Fit Plot', '', wildcard)

        fileName = str(fileName)
        if fileName == '':
            return
        
        path, ext = os.path.splitext(fileName) 
        ext = ext[1:].lower() 
        
       
        if ext != 'png' and ext != 'pdf' and ext != 'svg': 
            error_message = ( 
                  'Only the PNG, PDF and SVG image formats are supported.\n' 
                 'A file extension of `png\' or `pdf\' must be used.') 

            QtGui.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % error_message)
            return 
   
   
                          
        try: 

            matplotlib.rcParams['pdf.fonttype'] = 42
            
            fig = self.Specfig
            fig.savefig(fileName)

            
        except IOError, e:
            if e.strerror:
                err = e.strerror 
            else: 
                err = e 
   
            
            QtGui.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err)
            
            
        #Save text file with fit info
        textfilepath = path+'_'+self.anlz.xfspec_names[self.i_spec-1]+'_fitinfo.txt'
        f = open(textfilepath, 'w')
        print>>f, '*********************  Fit Results  ********************'

        print>>f, '\n'
        print>>f, 'Base:\t\t'+'{0:04.3f}'.format(self.anlz.xfitpars[self.i_spec-1].base)

        print>>f, '\n'
        text = 'Step Inflection Point [eV]:\t' 
        for i in range(self.nsteps[self.i_spec-1]):
            text = text + '{0:04.3f}'.format(self.anlz.xfitpars[self.i_spec-1].stepfitparams[i*3]) + ',  '
        print>>f, text
        
        text = 'Step Height:\t' 
        for i in range(self.nsteps[self.i_spec-1]):
            text = text + '{0:04.3f}'.format(self.anlz.xfitpars[self.i_spec-1].stepfitparams[i*3+1]) + ',  '
        print>>f, text
 
        text = 'Step FWHM:\t' 
        for i in range(self.nsteps[self.i_spec-1]):
            text = text + '{0:04.3f}'.format(self.anlz.xfitpars[self.i_spec-1].stepfitparams[i*3+2]) + ',  '
        print>>f, text       
        
        print>>f, '\n'
        text = 'Peak Positions:\t'
        for i in range(self.npeaks[self.i_spec-1]):
            text = text + '{0:04.3f}'.format(self.anlz.xfitpars[self.i_spec-1].gauss_fp_m[i]) + ',  '
        print>>f, text

        text = 'Peak Sigma:\t'
        for i in range(self.npeaks[self.i_spec-1]):
            text = text + '{0:04.3f}'.format(self.anlz.xfitpars[self.i_spec-1].gauss_fp_s[i]) + ',  '
        print>>f, text
        
        text = 'Peak Amplitude:\t'
        for i in range(self.npeaks[self.i_spec-1]):
            text = text + '{0:04.3f}'.format(self.anlz.xfitpars[self.i_spec-1].gauss_fp_a[i]) + ',  '
        print>>f, text        
        
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
            peak_engs = self.window().page5.peak_engs
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

              
            
                            

#---------------------------------------------------------------------- 
class FitParams(QtGui.QDialog):

    def __init__(self, parent, title, fitparams, readonly, initialization):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent
        
        self.base = fitparams[0]
        self.stepfitparams = fitparams[1]
        self.gauss_fp_a = fitparams[2]
        self.gauss_fp_m = fitparams[3]
        self.gauss_fp_s = fitparams[4]
        
        self.init = initialization        


        self.resize(600, 450)
        self.setWindowTitle(title)
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal) 
                
             
        vboxtop = QtGui.QVBoxLayout()
        
        #panel 
        sizer = QtGui.QGroupBox('Fit Parameters')
        vbox = QtGui.QVBoxLayout()
                

        fgs1 = QtGui.QGridLayout()
        

        
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Amplitude')  
        fgs1.addWidget(textctrl, 0, 1)      
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Position')  
        fgs1.addWidget(textctrl, 0, 2)   
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('FWHM')  
        fgs1.addWidget(textctrl, 0, 3)   
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Base')  
        fgs1.addWidget(textctrl, 0, 4)   


        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Step 1')
        fgs1.addWidget(textctrl, 1, 0)
        
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Step 2')
        fgs1.addWidget(textctrl, 2, 0)

        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Amplitude')  
        fgs1.addWidget(textctrl, 3, 1)      
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Mu')  
        fgs1.addWidget(textctrl, 3, 2)   
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Sigma')  
        fgs1.addWidget(textctrl, 3, 3)   
         

        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 1')
        fgs1.addWidget(textctrl, 4, 0)
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 2')
        fgs1.addWidget(textctrl, 5, 0)
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 3')
        fgs1.addWidget(textctrl, 6, 0)
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 4')
        fgs1.addWidget(textctrl, 7, 0)        
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 5')
        fgs1.addWidget(textctrl, 8, 0)
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 6')
        fgs1.addWidget(textctrl, 9, 0)
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 7')
        fgs1.addWidget(textctrl, 10, 0)                
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 8')
        fgs1.addWidget(textctrl, 11, 0)    
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 9')
        fgs1.addWidget(textctrl, 12, 0)
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 10')
        fgs1.addWidget(textctrl, 13, 0)        
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 11')
        fgs1.addWidget(textctrl, 14, 0)                    
        textctrl =  QtGui.QLabel(self)
        textctrl.setText('Gauss 12')
        fgs1.addWidget(textctrl, 15, 0)
        
        
        self.le_sa1 = QtGui.QLineEdit(self)
        self.le_sa1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sa1.setReadOnly(readonly)
        self.le_sa1.setText(str(self.stepfitparams[0]))
        fgs1.addWidget(self.le_sa1, 1, 1)

        self.le_sp1 = QtGui.QLineEdit(self)
        self.le_sp1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sp1.setReadOnly(readonly)
        self.le_sp1.setText(str(self.stepfitparams[1]))
        fgs1.addWidget(self.le_sp1, 1, 2)
        
        self.le_sf1 = QtGui.QLineEdit(self)
        self.le_sf1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sf1.setReadOnly(readonly)
        self.le_sf1.setText(str(self.stepfitparams[2]))
        fgs1.addWidget(self.le_sf1, 1, 3)
        
        self.le_base = QtGui.QLineEdit(self)
        self.le_base.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_base.setReadOnly(readonly)
        self.le_base.setText(str(self.base))
        fgs1.addWidget(self.le_base, 1, 4)

        
        self.le_sa2 = QtGui.QLineEdit(self)
        self.le_sa2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sa2.setReadOnly(readonly)
        self.le_sa2.setText(str(self.stepfitparams[3]))
        fgs1.addWidget(self.le_sa2, 2, 1)

        self.le_sp2 = QtGui.QLineEdit(self)
        self.le_sp2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sp2.setReadOnly(readonly)
        self.le_sp2.setText(str(self.stepfitparams[4]))
        fgs1.addWidget(self.le_sp2, 2, 2)
        
        self.le_sf2 = QtGui.QLineEdit(self)
        self.le_sf2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sf2.setReadOnly(readonly)
        self.le_sf2.setText(str(self.stepfitparams[5]))
        fgs1.addWidget(self.le_sf2, 2, 3)
        
                

        self.le_amplitudes = []
        self.le_mu = []
        self.le_sigma = []

        self.le_ga1 = QtGui.QLineEdit(self)
        self.le_ga1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga1, 4, 1)
        self.le_amplitudes.append(self.le_ga1)

        self.le_gm1 = QtGui.QLineEdit(self)
        self.le_gm1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm1, 4, 2)
        self.le_mu.append(self.le_gm1)
        
        self.le_gs1 = QtGui.QLineEdit(self)
        self.le_gs1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs1, 4, 3)
        self.le_sigma.append(self.le_gs1)
        
        self.le_ga2 = QtGui.QLineEdit(self)
        self.le_ga2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga2, 5, 1)
        self.le_amplitudes.append(self.le_ga2)

        self.le_gm2 = QtGui.QLineEdit(self)
        self.le_gm2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm2, 5, 2)
        self.le_mu.append(self.le_gm2)
        
        self.le_gs2 = QtGui.QLineEdit(self)
        self.le_gs2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs2, 5, 3)
        self.le_sigma.append(self.le_gs2)
        
        self.le_ga3 = QtGui.QLineEdit(self)
        self.le_ga3.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga3, 6, 1)
        self.le_amplitudes.append(self.le_ga3)

        self.le_gm3 = QtGui.QLineEdit(self)
        self.le_gm3.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm3, 6, 2)
        self.le_mu.append(self.le_gm3)
        
        self.le_gs3 = QtGui.QLineEdit(self)
        self.le_gs3.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs3, 6, 3)
        self.le_sigma.append(self.le_gs3)                                                

        self.le_ga4 = QtGui.QLineEdit(self)
        self.le_ga4.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga4, 7, 1)
        self.le_amplitudes.append(self.le_ga4)

        self.le_gm4 = QtGui.QLineEdit(self)
        self.le_gm4.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm4, 7, 2)
        self.le_mu.append(self.le_gm4)
        
        self.le_gs4 = QtGui.QLineEdit(self)
        self.le_gs4.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs4, 7, 3)
        self.le_sigma.append(self.le_gs4)
        
        self.le_ga5 = QtGui.QLineEdit(self)
        self.le_ga5.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga5, 8, 1)
        self.le_amplitudes.append(self.le_ga5)

        self.le_gm5 = QtGui.QLineEdit(self)
        self.le_gm5.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm5, 8, 2)
        self.le_mu.append(self.le_gm5)
        
        self.le_gs5 = QtGui.QLineEdit(self)
        self.le_gs5.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs5, 8, 3)
        self.le_sigma.append(self.le_gs5)
               
        self.le_ga6 = QtGui.QLineEdit(self)
        self.le_ga6.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga6, 9, 1)
        self.le_amplitudes.append(self.le_ga6)

        self.le_gm6 = QtGui.QLineEdit(self)
        self.le_gm6.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm6, 9, 2)
        self.le_mu.append(self.le_gm6)
        
        self.le_gs6 = QtGui.QLineEdit(self)
        self.le_gs6.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs6, 9, 3)
        self.le_sigma.append(self.le_gs6)
 
        self.le_ga7 = QtGui.QLineEdit(self)
        self.le_ga7.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga7, 10, 1)
        self.le_amplitudes.append(self.le_ga7)

        self.le_gm7 = QtGui.QLineEdit(self)
        self.le_gm7.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm7, 10, 2)
        self.le_mu.append(self.le_gm7)
        
        self.le_gs7 = QtGui.QLineEdit(self)
        self.le_gs7.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs7, 10, 3)
        self.le_sigma.append(self.le_gs7) 
       
        self.le_ga8 = QtGui.QLineEdit(self)
        self.le_ga8.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga8, 11, 1)
        self.le_amplitudes.append(self.le_ga8)

        self.le_gm8 = QtGui.QLineEdit(self)
        self.le_gm8.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm8, 11, 2)
        self.le_mu.append(self.le_gm8)
        
        self.le_gs8 = QtGui.QLineEdit(self)
        self.le_gs8.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs8, 11, 3)
        self.le_sigma.append(self.le_gs8) 
        
        self.le_ga9 = QtGui.QLineEdit(self)
        self.le_ga9.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga9, 12, 1)
        self.le_amplitudes.append(self.le_ga9)

        self.le_gm9 = QtGui.QLineEdit(self)
        self.le_gm9.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm9, 12, 2)
        self.le_mu.append(self.le_gm9)
        
        self.le_gs9 = QtGui.QLineEdit(self)
        self.le_gs9.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs9, 12, 3)
        self.le_sigma.append(self.le_gs9)     
        
        self.le_ga10 = QtGui.QLineEdit(self)
        self.le_ga10.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga10, 13, 1)
        self.le_amplitudes.append(self.le_ga10)

        self.le_gm10 = QtGui.QLineEdit(self)
        self.le_gm10.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm10, 13, 2)
        self.le_mu.append(self.le_gm10)
        
        self.le_gs10 = QtGui.QLineEdit(self)
        self.le_gs10.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs10, 13, 3)
        self.le_sigma.append(self.le_gs10)

        self.le_ga11 = QtGui.QLineEdit(self)
        self.le_ga11.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga11, 14, 1)
        self.le_amplitudes.append(self.le_ga11)

        self.le_gm11 = QtGui.QLineEdit(self)
        self.le_gm11.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm11, 14, 2)
        self.le_mu.append(self.le_gm11)
        
        self.le_gs11 = QtGui.QLineEdit(self)
        self.le_gs11.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs11, 14, 3)
        self.le_sigma.append(self.le_gs11)
               
        self.le_ga12 = QtGui.QLineEdit(self)
        self.le_ga12.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga12, 15, 1)
        self.le_amplitudes.append(self.le_ga12)

        self.le_gm12 = QtGui.QLineEdit(self)
        self.le_gm12.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm12, 15, 2)
        self.le_mu.append(self.le_gm12)
        
        self.le_gs12 = QtGui.QLineEdit(self)
        self.le_gs12.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs12, 15, 3)
        self.le_sigma.append(self.le_gs12)
        
        for i in range(12):
            self.le_amplitudes[i].setReadOnly(readonly)
            self.le_mu[i].setReadOnly(readonly)
            self.le_sigma[i].setReadOnly(readonly)
            
            self.le_amplitudes[i].setText(str(self.gauss_fp_a[i]))
            self.le_mu[i].setText(str(self.gauss_fp_m[i]))
            self.le_sigma[i].setText(str(self.gauss_fp_s[i]))
                    
        vbox.addLayout(fgs1)
        sizer.setLayout(vbox)
        
        vboxtop.addWidget(sizer)
        
        if initialization:
            bname = 'Confirm'
        else:
            bname = 'Close'
        self.button_close = QtGui.QPushButton(bname)
        self.button_close.clicked.connect( self.OnClose)   
        vboxtop.addWidget(self.button_close)
        
        
        self.setLayout(vboxtop)
                                 

#----------------------------------------------------------------------
    def OnClose(self, event):
        
        
        if self.init:
            
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
            
            
            self.parent.page6.SetInitFitParams(self.base, self.stepfitparams, self.gauss_fp_a, self.gauss_fp_m, self.gauss_fp_s)
        
        self.close()             
        
                            
                
""" ------------------------------------------------------------------------------------------------"""
class PagePeakID(QtGui.QWidget):
    def __init__(self, common, data_struct, stack, anlz):
        super(PagePeakID, self).__init__()

        self.initUI(common, data_struct, stack, anlz)
        
#----------------------------------------------------------------------          
    def initUI(self, common, data_struct, stack, anlz): 
        
        
        self.com = common 
        self.data_struct = data_struct
        self.stk = stack       
        self.anlz = anlz
        
       
        
        pw = PlotW*0.8
        ph = PlotH*0.8
        
                
        self.peak_engs = [285.0, 286.5, 288.5, 289.5, 290.5]
        self.peak_names = ['aromatic', 'phenolic', 'carboxyl', 'alkyl', 'carbonyl']
        
        self.i_spectrum = 1
        self.i_eng = 0
        
        self.n_spectra = 0
        
        self.spectrum_loaded = 0
        
          
        #panel 1        
        vbox1 = QtGui.QVBoxLayout()
               
        self.tc_1 = QtGui.QLabel(self)
        self.tc_1.setText("X-ray Spectrum")    
        hbox11 = QtGui.QHBoxLayout() 
        
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
   
        self.kespecfig = Figure((pw, ph))
        self.KESpecPan = FigureCanvas(self.kespecfig)
        self.KESpecPan.mpl_connect('button_press_event', self.OnPointSpectrum)
        fbox.addWidget(self.KESpecPan)
        frame.setLayout(fbox)

        self.slider_spec = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slider_spec.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_spec.setEnabled(False)
        self.slider_spec.valueChanged[int].connect(self.OnSScroll)
        self.slider_spec.setRange(1, 5)
                    
        hbox11.addWidget(frame)
        hbox11.addWidget(self.slider_spec)

                            
        vbox1.addStretch(1)
        vbox1.addWidget(self.tc_1)        
        vbox1.addLayout(hbox11)
        vbox1.addStretch(1)
 
       
        
                    
        #panel 2
        sizer2 = QtGui.QGroupBox('Load Spectrum')
        vbox2 = QtGui.QVBoxLayout()     
        
        self.button_addclspec = QtGui.QPushButton('Add Cluster Spectra')
        self.button_addclspec.clicked.connect( self.OnAddClusterSpectra)   
        self.button_addclspec.setEnabled(False)
        vbox2.addWidget(self.button_addclspec)
           
         
        self.button_loadtspec = QtGui.QPushButton('Load Spectrum')
        self.button_loadtspec.clicked.connect(self.OnSpecFromFile)
        vbox2.addWidget(self.button_loadtspec)
            
 
        self.button_save = QtGui.QPushButton('Save Images...')
        self.button_save.clicked.connect( self.OnSave)
        self.button_save.setEnabled(False)
        vbox2.addWidget(self.button_save)

        sizer2.setLayout(vbox2)


        #panel 3
        vbox3 = QtGui.QVBoxLayout()
    
        t1 = QtGui.QLabel(self)
        t1.setText("Peak ID Energies")       
         
        self.lc_1 = QListWidget()   
        self.lc_1.itemClicked.connect(self.OnEngListClick)
        self.lc_1.setMinimumSize(200, 400)
         
        vbox3.addWidget(t1)
        vbox3.addWidget(self.lc_1)
                 
      
        #panel 4     
        vbox4 = QtGui.QVBoxLayout()
         
        self.tc_imageeng = QtGui.QLabel(self)
        self.tc_imageeng.setText("Image at peak energy: ")
        vbox4.addWidget(self.tc_imageeng)
         
        hbox41 = QtGui.QHBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
  
        self.absimgfig = Figure((ph, ph))
        self.AbsImagePanel = FigureCanvas(self.absimgfig)
   
        fbox.addWidget(self.AbsImagePanel)
        frame.setLayout(fbox)
        hbox41.addWidget(frame)       

        self.slider_eng = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slider_eng.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_eng.valueChanged[int].connect(self.OnScrollEng)
        self.slider_eng.setRange(0, len(self.peak_names)-1)
        hbox41.addWidget(self.slider_eng)
        hbox41.addStretch(1)
        vbox4.addLayout(hbox41)
         
 

        

        vboxtop1 = QtGui.QVBoxLayout()
                     
        vboxtop1.addStretch(1) 
        vboxtop1.addWidget(sizer2)
        vboxtop1.addStretch(1) 
        vboxtop1.addLayout(vbox3)
         
        vboxtop2 = QtGui.QVBoxLayout()  
        vboxtop2.addStretch(1)       
        vboxtop2.addLayout(vbox1)
        vboxtop2.addStretch(1)
        vboxtop2.addLayout(vbox4)

                
        hboxtop = QtGui.QHBoxLayout()  
        hboxtop.addStretch(1)     
        hboxtop.addLayout(vboxtop1)
        hboxtop.addStretch(1) 
        hboxtop.addLayout(vboxtop2)
        hboxtop.addStretch(1) 

        
        hboxtop.setContentsMargins(20,20,20,20)
        self.setLayout(hboxtop) 


        
        self.ShowPeakEngs()


#----------------------------------------------------------------------
    def OnSpecFromFile(self, event):
        

        try:
            
            wildcard = "Spectrum files (*.csv)"
            
            filepath = QtGui.QFileDialog.getOpenFileName(self, 'Choose Spectrum file', '', wildcard)
            

            filepath = str(filepath)
            if filepath == '':
                return
            
            self.filename =  os.path.basename(str(filepath))
            directory =  os.path.dirname(str(filepath))
            self.DefaultDir = directory            
                                                        
            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))    
                                            
            self.anlz.load_xraypeakfit_spectrum(filename=filepath)

            self.com.xpf_loaded = 1
            if self.window().page6 != None:
                self.window().page6.spectrumfitted.append(0)
                self.window().page6.firstinit.append(1)
                self.window().page6.initdone.append(0)
                self.window().page6.nsteps.append(1)
                self.window().page6.npeaks.append(4)
                self.window().page6.fits.append(0)
                self.window().page6.fits_sep.append(0)
                
                self.window().page6.i_spec = self.anlz.n_xrayfitsp      
                self.window().page6.slider_spec.setMaximum(self.anlz.n_xrayfitsp)
                self.window().page6.slider_spec.setValue(self.window().page6.i_spec)
                self.window().page6.nstepsspin.setValue(self.window().page6.nsteps[self.window().page6.i_spec-1])
                self.window().page6.ngaussspin.setValue(self.window().page6.npeaks[self.window().page6.i_spec-1])                
                self.window().page6.loadSpectrum()
                self.window().page6.updatewidgets()

            self.spectrum_loaded = 1
            self.updatewidgets()
            
            
            QtGui.QApplication.restoreOverrideCursor()
            
        except:
            QtGui.QApplication.restoreOverrideCursor()  
            QtGui.QMessageBox.warning(self, 'Error', 'Spectrum file not loaded.')
                                   
                                 
        
        self.window().refresh_widgets()
        
        self.slider_spec.setMaximum(self.n_spectra)
        self.slider_spec.setEnabled(True)
        
        self.ShowODPlot()
        

#----------------------------------------------------------------------
    def OnAddClusterSpectra(self, event):
        
        try:

            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor)) 
            

            for i in range(self.anlz.nclusters):
                self.anlz.load_xraypeakfit_clusterspectrum(i)

                if self.window().page6 != None:
                    self.window().page6.spectrumfitted.append(0)
                    self.window().page6.firstinit.append(1)
                    self.window().page6.initdone.append(0)
                    self.window().page6.nsteps.append(1)
                    self.window().page6.npeaks.append(4)
                    self.window().page6.fits.append(0)
                    self.window().page6.fits_sep.append(0)
                
                
            self.spectrum_loaded = 1
            
            
            self.com.xpf_loaded = 1
            if self.window().page6 != None:
                self.window().page6.i_spec = self.anlz.n_xrayfitsp      
                self.window().page6.slider_spec.setMaximum(self.anlz.n_xrayfitsp)
                self.window().page6.slider_spec.setValue(self.window().page6.i_spec)
                self.window().page6.nstepsspin.setValue(self.window().page6.nsteps[self.window().page6.i_spec-1])
                self.window().page6.ngaussspin.setValue(self.window().page6.npeaks[self.window().page6.i_spec-1])                
                self.window().page6.loadSpectrum()
                self.window().page6.updatewidgets()

                    
            QtGui.QApplication.restoreOverrideCursor()
            
        except:
            QtGui.QApplication.restoreOverrideCursor()  
            QtGui.QMessageBox.warning(self, 'Error', 'Cluster spectra not loaded.')
                                   
                                 
        
        
        self.window().refresh_widgets()
        
        self.updatewidgets()
        
        self.ShowODPlot()
        self.ShowImage()

#----------------------------------------------------------------------        
    def OnSScroll(self, value):
        
        
        sel = value
        self.i_spectrum = sel
                

        self.ShowODPlot()
        self.ShowImage()
                           
#----------------------------------------------------------------------            
    def OnScrollEng(self, value):
        self.i_eng = value


        self.ShowImage()
        self.ShowODPlot()



#----------------------------------------------------------------------  
    def OnPointSpectrum(self, evt):
        x = evt.xdata
        
        if (self.com.i0_loaded != 0) or (self.spectrum_loaded == 1):      
            if x < self.stk.ev[0]:
                sel_ev = 0
            elif x > self.stk.ev[self.stk.n_ev-1]:
                sel_ev = self.stk.n_ev-1
            else:
                indx = np.abs(self.stk.ev - x).argmin()
                sel_ev = indx
                
            
            self.i_eng=(np.abs(self.peak_engs-self.stk.ev[sel_ev])).argmin()                 

            self.ShowImage()
            self.ShowODPlot()
            
            self.slider_eng.setValue(self.i_eng)
                              
#----------------------------------------------------------------------     
    def ShowODPlot(self):

        if (self.spectrum_loaded == 0) and (self.com.i0_loaded == 0):
            return
        
        if (self.com.i0_loaded == 0):
            spectrum = self.anlz.xrayfitspectra[self.i_spectrum-1, :]
            self.tc_1.setText("X-ray Spectrum: " + 
                           self.anlz.xfspec_names[self.i_spectrum-1])            
        else:
            if (self.i_spectrum == 1):
                spectrum = self.stk.od3d.sum(axis=0)   
                spectrum = spectrum.sum(axis=0)/(self.stk.n_rows*self.stk.n_cols)      
            
                self.tc_1.setText("X-ray Spectrum: Average OD")
            else:
                spectrum = self.anlz.xrayfitspectra[self.i_spectrum-2, :]
                self.tc_1.setText("X-ray Spectrum: " + 
                               self.anlz.xfspec_names[self.i_spectrum-2])                

                
        fig = self.kespecfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        
        
        specplot = axes.plot(self.stk.ev,spectrum)

        for i in range(len(self.peak_engs)):
            if (self.peak_engs[i]>self.stk.ev[0]) and (self.peak_engs[i]<self.stk.ev[-1]):
                axes.axvline(x=self.peak_engs[i], color = 'g', alpha=0.5)

        if (self.peak_engs[self.i_eng]>self.stk.ev[0]) and (self.peak_engs[self.i_eng]<self.stk.ev[-1]):
            axes.axvline(x=self.peak_engs[self.i_eng], color = 'g', alpha=1.0)                
                        
        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Optical Density')
        
        self.KESpecPan.draw()

#----------------------------------------------------------------------        
    def OnEngListClick(self):
        item = self.lc_1.currentRow()
                
        self.i_eng = item
        
        #self.ShowPeakEngs()     
        self.ShowImage()
        self.ShowODPlot()
        
#----------------------------------------------------------------------           
    def ShowPeakEngs(self):    
        
        self.lc_1.clear()
        
        for i in range(len(self.peak_engs)):
            self.lc_1.addItem('{0:8.2f}'.format(self.peak_engs[i])+'\t'+self.peak_names[i])
            

#----------------------------------------------------------------------        
    def ShowImage(self):
        
        if self.com.i0_loaded == 0:
            return
            
        iev=(np.abs(self.stk.ev-self.peak_engs[self.i_eng])).argmin() 
        image = self.stk.absdata[:,:,int(iev)].copy() 

        fig = self.absimgfig
        fig.clf()
        fig.add_axes(((0.0,0.0,1.0,1.0)))
        axes = fig.gca()
        fig.patch.set_alpha(1.0)
        
        im = axes.imshow(np.rot90(image), cmap=matplotlib.cm.get_cmap("gray")) 
        
        #Show Scale Bar
        startx = int(self.stk.n_rows*0.05)
        starty = int(self.stk.n_cols*0.05)+self.stk.scale_bar_pixels_y
        um_string = ' $\mathrm{\mu m}$'
        microns = '$'+self.stk.scale_bar_string+' $'+um_string
        axes.text(self.stk.scale_bar_pixels_x+startx+1,starty+1, microns, horizontalalignment='left', verticalalignment='center',
                  color = 'black', fontsize=14)
        #Matplotlib has flipped scales so I'm using rows instead of cols!
        p = matplotlib.patches.Rectangle((startx,starty), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                               color = 'black', fill = True)
        axes.add_patch(p)
            
       
        axes.axis("off")      
        self.AbsImagePanel.draw()
        
        self.tc_imageeng.setText('Image at peak energy: {0:5.2f} eV'.format(float(self.stk.ev[iev])))
 
        self.lc_1.setCurrentRow(self.i_eng) 

#----------------------------------------------------------------------
    def updatewidgets(self):

        if self.spectrum_loaded:
            self.n_spectra = self.anlz.n_xrayfitsp
        else:
            self.n_spectra = 0
        
        if self.com.i0_loaded == 1:
            self.n_spectra += 1
            
        self.i_spectrum = self.n_spectra
            
        self.slider_spec.setMaximum(self.n_spectra)
        self.slider_spec.setEnabled(True)
        self.slider_spec.setValue(self.i_spectrum)
        self.ShowODPlot()
        self.ShowImage()        
                
        
#----------------------------------------------------------------------
    def OnSave(self, event):


        #Save images      
        wildcard = "Portable Network Graphics (*.png);;Adobe PDF Files (*.pdf);;"

        fileName = QtGui.QFileDialog.getSaveFileName(self, 'Save OD Plot with Key Energies', '', wildcard)

        fileName = str(fileName)
        if fileName == '':
            return
        
        path, ext = os.path.splitext(fileName) 
        ext = ext[1:].lower() 
        
       
        if ext != 'png' and ext != 'pdf': 
            error_message = ( 
                  'Only the PNG and PDF image formats are supported.\n' 
                 'A file extension of `png\' or `pdf\' must be used.') 

            QtGui.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % error_message)
            return 
   
   
                          
        try: 

            matplotlib.rcParams['pdf.fonttype'] = 42
            
            fig = self.kespecfig
            fig.savefig(fileName)

            
        except IOError, e:
            if e.strerror:
                err = e.strerror 
            else: 
                err = e 
   
            
            QtGui.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err)
            
              
            
        return
          
 
    

""" ------------------------------------------------------------------------------------------------"""
class PageSpectral(QtGui.QWidget):
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
        

        vbox = QtGui.QVBoxLayout()
        hboxT = QtGui.QHBoxLayout()
        hboxB = QtGui.QHBoxLayout()
    
    
        #panel 5
        sizer5 = QtGui.QGroupBox('Target Spectra')
        vbox5 = QtGui.QVBoxLayout()
 
        self.tc_speclist =  QListWidget()  
        self.tc_speclist.itemClicked.connect(self.OnSpectraListClick)
        #self.tc_speclist.doubleClicked.connect(self.OnEditSpectraListClick)      
        vbox5.addWidget(self.tc_speclist)
        sizer5.setLayout(vbox5)
        
        
        
        #panel 1        
        vbox1 = QtGui.QVBoxLayout()     
  
        self.tc_spmap = QtGui.QLabel(self)
        self.tc_spmap.setText("Spectrum composition map")

        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()

        self.mapfig = Figure((PlotH, PlotH))
        self.MapPanel = FigureCanvas(self.mapfig)
        self.MapPanel.setParent(self)
        fbox.addWidget(self.MapPanel)
        frame.setLayout(fbox)
        
        vbox1.addWidget(self.tc_spmap)        
        vbox1.addWidget(frame)
  
                  
     
     
        #panel 2
        vbox2 = QtGui.QVBoxLayout()
         
        self.tc_tspec = QtGui.QLabel(self)
        self.tc_tspec.setText("Target Spectrum: ")
        hbox11 = QtGui.QHBoxLayout()       
         
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
        self.TSpecfig = Figure((PlotW, PlotH))
        self.TSpectrumPanel = FigureCanvas(self.TSpecfig)
        fbox.addWidget(self.TSpectrumPanel)
        frame.setLayout(fbox)
        
        self.slider_tspec = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slider_tspec.setFocusPolicy(QtCore.Qt.StrongFocus)
        #self.slider_tspec.setEnabled(False)
        self.slider_tspec.valueChanged[int].connect(self.OnTSScroll)
        self.slider_tspec.setRange(1, 5)
           
 
        hbox11.addWidget(frame)
        hbox11.addWidget(self.slider_tspec)
           
        vbox2.addWidget(self.tc_tspec)       
        vbox2.addLayout(hbox11)
        
        self.slider_theta = QtGui.QScrollBar(QtCore.Qt.Horizontal)
        self.slider_theta.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_theta.valueChanged[int].connect(self.OnScrollTheta)
        self.slider_theta.setRange(0, 100)  
        self.slider_theta.setMinimumWidth(350)   
        self.slider_theta.setVisible(False)  
        self.tc_imagetheta = QtGui.QLabel(self)
        self.tc_imagetheta.setText("4D Data Angle: ")
        self.tc_imagetheta.setVisible(False)
        hbox21 = QtGui.QHBoxLayout()
        hbox21.addWidget(self.tc_imagetheta) 
        hbox21.addStretch(1)
        hbox21.addWidget(self.slider_theta)
        hbox21.addStretch(1)
        vbox2.addLayout(hbox21)
         
         
         
        #panel 3
        sizer3 = QtGui.QGroupBox('Target Spectra')
        vbox3 = QtGui.QVBoxLayout()
        

        self.button_addclspec = QtGui.QPushButton('Add Cluster Spectra')
        self.button_addclspec.setMinimumSize (180 , 0)
        self.button_addclspec.clicked.connect( self.OnAddClusterSpectra)   
        self.button_addclspec.setEnabled(False)
        vbox3.addWidget(self.button_addclspec)         
        self.button_loadtspec = QtGui.QPushButton('Load Spectrum')
        self.button_loadtspec.clicked.connect(self.OnTSpecFromFile)
        self.button_loadtspec.setEnabled(False)
        vbox3.addWidget(self.button_loadtspec)
        self.button_addflat = QtGui.QPushButton('Add Flat Spectrum')
        self.button_addflat.clicked.connect(self.OnFlatTSpec)
        self.button_addflat.setEnabled(False)
        vbox3.addWidget(self.button_addflat)

         
        self.button_showrgb = QtGui.QPushButton('Composite RGB image...')
        self.button_showrgb.clicked.connect(self.OnCompositeRGB)   
        self.button_showrgb.setEnabled(False)
        vbox3.addWidget(self.button_showrgb)      
        
        self.button_histogram = QtGui.QPushButton('Histogram Value Cutoff...')
        self.button_histogram.clicked.connect(self.OnHistogram)   
        self.button_histogram.setEnabled(False)
        vbox3.addWidget(self.button_histogram)     
 
        self.button_save = QtGui.QPushButton('Save Images...')
        self.button_save.clicked.connect(self.OnSave)
        self.button_save.setEnabled(False)
        vbox3.addWidget(self.button_save)
        
        self.button_calc4d = QtGui.QPushButton('Calculate for all angles')
        self.button_calc4d.clicked.connect(self.OnCalc4D)
        self.button_calc4d.setEnabled(False)
        self.button_calc4d.setVisible(False)
        vbox3.addWidget(self.button_calc4d)

        sizer3.setLayout(vbox3)
         
 
         
        #panel 4
        sizer4 = QtGui.QGroupBox('Display')
        vbox4 = QtGui.QVBoxLayout()
                

        sb = QtGui.QGroupBox('Spectrum')
        vbox41 = QtGui.QVBoxLayout()
        vbox41.setSpacing(10)

        self.textctrl_sp1 =  QtGui.QLabel(self)
        vbox41.addWidget(self.textctrl_sp1)
        self.textctrl_sp2 =  QtGui.QLabel(self)
        vbox41.addWidget(self.textctrl_sp2)
               
        
         
        self.textctrl_sp1.setText('Common Name: ')
        self.textctrl_sp2.setText('RMS Error: ')
        
        sb.setLayout(vbox41)
        vbox4.addWidget(sb)


         
        
        hbox4b = QtGui.QHBoxLayout() 
          
        sb = QtGui.QGroupBox('Composition Map')
        vbox42 = QtGui.QVBoxLayout()
        
        
        self.rb_raw = QtGui.QRadioButton( 'Raw', self)
        self.rb_fit = QtGui.QRadioButton('Fitted',self)
        self.rb_raw.setChecked(True)
        self.rb_raw.toggled.connect(self.OnRBRawFit)
        

        vbox42.addWidget(self.rb_raw)
        vbox42.addWidget(self.rb_fit)


        self.add_scale_cb = QtGui.QCheckBox('Scalebar', self) 
        #self.add_scale_cb.toggle()
        self.add_scale_cb.stateChanged.connect(self.OnShowScale)
        vbox42.addWidget(self.add_scale_cb)
        sb.setLayout(vbox42)
        hbox4b.addWidget(sb)
        
        
        
                  
        sb = QtGui.QGroupBox('Fit Weights')
        vbox43 = QtGui.QVBoxLayout()

        self.tc_spfitlist = QListWidget()
        vbox43.addWidget(self.tc_spfitlist)
          
        sb.setLayout(vbox43)
        hbox4b.addWidget(sb)
          
          
        vbox44 = QtGui.QVBoxLayout()
        self.button_removespec = QtGui.QPushButton('Remove Spectrum')
        self.button_removespec.clicked.connect(self.OnRemoveSpectrum)   
        self.button_removespec.setEnabled(False)    
        vbox44.addWidget(self.button_removespec)
        self.button_movespup = QtGui.QPushButton('Move Spectrum Up')
        self.button_movespup.clicked.connect(self.OnMoveSpectrumUp)   
        self.button_movespup.setEnabled(False)  
        vbox44.addWidget(self.button_movespup)
        self.button_movespdown = QtGui.QPushButton('Move Spectrum Down')
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
        

        try: 
        
            wildcard = "Spectrum files (*.csv);;Spectrum files (*.xas);;"
            
            #filepath = QtGui.QFileDialog.getOpenFileName(self, 'Choose Spectrum file', '', wildcard, self.DefaultDir)
            filepath = QtGui.QFileDialog.getOpenFileName(self, 'Choose Spectrum file', '', wildcard)
            

            filepath = str(filepath)
            if filepath == '':
                return
            
            self.filename =  os.path.basename(str(filepath))
            directory =  os.path.dirname(str(filepath))
            self.DefaultDir = directory            
                                                        
            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))    
                                            
            self.anlz.read_target_spectrum(filename=filepath)
            self.com.spec_anl_calculated = 1
            
            self.i_tspec = self.anlz.n_target_spectra      
            self.slider_tspec.setMaximum(self.anlz.n_target_spectra)
            self.slider_tspec.setValue(self.i_tspec)
            
            
            self.loadTSpectrum()
            self.loadTargetMap()    
            self.ShowSpectraList()
                    
            QtGui.QApplication.restoreOverrideCursor()
            
        except:
            QtGui.QApplication.restoreOverrideCursor()  
            QtGui.QMessageBox.warning(self, 'Error', 'Spectrum file not loaded.')
                                   
                                 
        self.window().refresh_widgets()
        

#----------------------------------------------------------------------
    def OnFlatTSpec(self, event):

        try: 
            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor)) 
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
        
            QtGui.QApplication.restoreOverrideCursor()
            
        except:
            QtGui.QApplication.restoreOverrideCursor()  
            QtGui.QMessageBox.warning(self, 'Error', 'Flat spectrum not loaded.')
            
                                                      
        self.window().refresh_widgets()
        
#----------------------------------------------------------------------
    def OnAddClusterSpectra(self, event):

        try:

            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor)) 
            self.anlz.add_cluster_target_spectra()
            self.com.spec_anl_calculated = 1
            
            self.i_tspec = self.anlz.n_target_spectra      
            self.slider_tspec.setMaximum(self.anlz.n_target_spectra)
            self.slider_tspec.setValue(self.i_tspec)
            
            
            self.ShowSpectraList() 
            self.loadTSpectrum()
            self.loadTargetMap()  
             
        
            QtGui.QApplication.restoreOverrideCursor()
            
        except:
            QtGui.QApplication.restoreOverrideCursor()  
            QtGui.QMessageBox.warning(self, 'Error', 'Cluster spectra not loaded.')
 
                                                        
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

            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor)) 
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
            
                     
        
            QtGui.QApplication.restoreOverrideCursor()
            
        except:
            QtGui.QApplication.restoreOverrideCursor()  
            QtGui.QMessageBox.warning(self, 'Error', 'Could not calculate 4D spectra.')
 
                                                        
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
                
                
            
        except IOError, e:
            if e.strerror:
                err = e.strerror 
            else: 
                err = e 
            QtGui.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err)

            
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
                fig.savefig(fileName_img, bbox_inches='tight', pad_inches = 0.0)
                
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
        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor)) 
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
        
        QtGui.QApplication.restoreOverrideCursor()
        
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

        if self.com.pca_calculated == 1:       
            tspectrumfit = self.anlz.target_pcafit_spectra[self.i_tspec-1, :]  
            diff = np.abs(tspectrum-tspectrumfit)      
            line2 = axes.plot(self.stk.ev,tspectrumfit, color='green', label = 'Fit')
        
            line3 = axes.plot(self.stk.ev,diff, color='grey', label = 'Abs(Raw-Fit)')
        

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
        if self.com.pca_calculated == 1:       
            self.textctrl_sp2.setText('RMS Error: '+ str('{0:7.5f}').format(self.anlz.target_rms[self.i_tspec-1]))
            
            self.ShowFitWeights()
        


#---------------------------------------------------------------------- 
class ShowCompositeRBGmap(QtGui.QDialog):

    def __init__(self, parent, common, analz):    
        QtGui.QWidget.__init__(self, parent)
        
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
        

        sizer1 = QtGui.QGroupBox('Red spectrum')

        fgs1 = QtGui.QGridLayout()
        
        r = QtGui.QLabel(self) 
        r.setText('Red')
        rl = QtGui.QLabel(self) 
        rl.setText('Limits')
        rw = QtGui.QLabel(self) 
        rw.setText('Weight')     
        
        
        self.combor = QtGui.QComboBox(self)
        self.combor.addItems(self.anlz.tspec_names)
        self.combor.activated[int].connect(self.OnSelectR) 

        #self.combor.SetToolTip(wx.ToolTip("select spectrum from dropdown-list"))
        self.combor.setCurrentIndex(self.r_spec)
        
        hbox12 = QtGui.QHBoxLayout() 
        
        
        self.tcrmin = QtGui.QSpinBox()
        self.tcrmin.setRange(0,100)
        self.tcrmin.setValue(0)
        self.tcrmin.valueChanged[int].connect(self.OnLimitMinR)
        
        self.tcrmax = QtGui.QSpinBox()
        self.tcrmax.setRange(0,100)
        self.tcrmax.setValue(100)
        self.tcrmax.valueChanged[int].connect(self.OnLimitMaxR)        
                        
        hbox12.addWidget(self.tcrmin)
        hbox12.addWidget(self.tcrmax)
        
        self.tcrweight = QtGui.QSpinBox()
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
        
        
  
        sizer2 = QtGui.QGroupBox('Green spectrum')

        fgs2 = QtGui.QGridLayout()
        
        g = QtGui.QLabel(self) 
        g.setText('Green')
        gl = QtGui.QLabel(self) 
        gl.setText('Limits')
        gw = QtGui.QLabel(self) 
        gw.setText('Weight')     
        
        
        self.combog = QtGui.QComboBox(self)
        self.combog.addItems(self.anlz.tspec_names)
        self.combog.activated[int].connect(self.OnSelectG) 

        #self.combor.SetToolTip(wx.ToolTip("select spectrum from dropdown-list"))
        self.combog.setCurrentIndex(self.g_spec)
        
        hbox12 = QtGui.QHBoxLayout() 
        
        
        self.tcgmin = QtGui.QSpinBox()
        self.tcgmin.setRange(0,100)
        self.tcgmin.setValue(0)
        self.tcgmin.valueChanged[int].connect(self.OnLimitMinG)
        
        self.tcgmax = QtGui.QSpinBox()
        self.tcgmax.setRange(0,100)
        self.tcgmax.setValue(100)
        self.tcgmax.valueChanged[int].connect(self.OnLimitMaxG)        
                        
        hbox12.addWidget(self.tcgmin)
        hbox12.addWidget(self.tcgmax)
        
        self.tcgweight = QtGui.QSpinBox()
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
                 


        sizer3 = QtGui.QGroupBox('Blue spectrum')

        fgs3 = QtGui.QGridLayout()
        
        b = QtGui.QLabel(self) 
        b.setText('Blue')
        bl = QtGui.QLabel(self) 
        bl.setText('Limits')
        bw = QtGui.QLabel(self) 
        bw.setText('Weight')     
        
        
        self.combob = QtGui.QComboBox(self)
        self.combob.addItems(self.anlz.tspec_names)
        self.combob.activated[int].connect(self.OnSelectB) 

        #self.combor.SetToolTip(wx.ToolTip("select spectrum from dropdown-list"))
        self.combob.setCurrentIndex(self.b_spec)
        
        hbox12 = QtGui.QHBoxLayout() 
        
        
        self.tcbmin = QtGui.QSpinBox()
        self.tcbmin.setRange(0,100)
        self.tcbmin.setValue(0)
        self.tcbmin.valueChanged[int].connect(self.OnLimitMinB)
        
        self.tcbmax = QtGui.QSpinBox()
        self.tcbmax.setRange(0,100)
        self.tcbmax.setValue(100)
        self.tcbmax.valueChanged[int].connect(self.OnLimitMaxB)        
                        
        hbox12.addWidget(self.tcbmin)
        hbox12.addWidget(self.tcbmax)
        
        self.tcbweight = QtGui.QSpinBox()
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
        
          

        
        vbox = QtGui.QVBoxLayout()
        hbox1 = QtGui.QHBoxLayout()
        vbox1 = QtGui.QVBoxLayout()
                       
        vbox1.addWidget(sizer1)
        vbox1.addWidget(sizer2)
        vbox1.addWidget(sizer3)
        
        
        
        self.show_info_cb = QtGui.QCheckBox( 'Show Info on the Image', self)
        self.show_info_cb.stateChanged.connect(self.OnShowInfo)
        vbox1.addWidget(self.show_info_cb)
        
        hbox1.addLayout(vbox1)
        
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
        self.RGBImagefig = Figure((PlotH, PlotH))
        self.RGBImagePanel = FigureCanvas(self.RGBImagefig)
        fbox.addWidget(self.RGBImagePanel)
        frame.setLayout(fbox)

        hbox1.addWidget(frame)
        
        
        vbox.addLayout(hbox1) 
        
             
        hbox2 = QtGui.QHBoxLayout()
        
        button_save = QtGui.QPushButton('Save image')
        button_save.clicked.connect( self.OnSave)   
        hbox2.addWidget(button_save)
        
        button_close = QtGui.QPushButton('Dismiss')
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
        
        if self.parent.page4.showraw == True:
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

        if self.parent.page4.showraw == True:
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
        
        if self.parent.page4.showraw == True:
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

        SaveFileName = QtGui.QFileDialog.getSaveFileName(self, 'Save Plot', '', wildcard)

        SaveFileName = str(SaveFileName)
        if SaveFileName == '':
            return
        
        path, ext = os.path.splitext(SaveFileName) 
        ext = ext[1:].lower() 
        
       
        if ext != 'png' and ext != 'pdf': 
            error_message = ( 
                  'Only the PNG and PDF image formats are supported.\n' 
                 'A file extension of `png\' or `pdf\' must be used.') 

            QtGui.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % error_message)
            return 
   

        matplotlib.rcParams['pdf.fonttype'] = 42

                            
        fig = self.RGBImagefig
        fig.savefig(SaveFileName, bbox_inches='tight', pad_inches = 0.0)
                
   

#---------------------------------------------------------------------- 
class ShowMapHistogram(QtGui.QDialog):

    def __init__(self, parent, common, analz):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent
        
        self.limit = 1
        
        self.histmax = None

        
        self.resize(600, 500)
        self.setWindowTitle('Histogram')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)
                
        self.com = common 
        self.anlz = analz
      
        vbox = QtGui.QVBoxLayout()
               
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
   
        self.histfig = Figure((6.0, 4.2))
        self.HistogramPanel = FigureCanvas(self.histfig)
        self.HistogramPanel.setParent(self)
        self.HistogramPanel.mpl_connect('button_press_event', self.OnClick)
        
        
        fbox.addWidget(self.HistogramPanel)
        frame.setLayout(fbox)
        vbox.addWidget(frame)
        
        
                        
        vbox1 = QtGui.QVBoxLayout()
        sizer1 = QtGui.QGroupBox('Histogram Cutoff')
        
        st = QtGui.QLabel(self) 
        st.setText('Select a cutoff values on the histogram. All the values outside the defined limits will be set to cutoff limit value.')
     
        vbox1.addWidget(st)
           
        hbox1 = QtGui.QHBoxLayout()
        self.rb_min = QtGui.QRadioButton( 'Lower Limit', self)
        self.rb_max = QtGui.QRadioButton('Upper Limit',self)
        self.rb_min.setChecked(True)
        self.rb_min.toggled.connect(self.OnRb_limit)
        

        hbox1.addWidget(self.rb_min)
        hbox1.addWidget(self.rb_max)
        hbox1.addStretch (1)
        vbox1.addLayout(hbox1)
        
        self.tl_cutmin = QtGui.QLabel(self) 
        self.tl_cutmin.setText('Lower Cutoff Value: ')
        vbox1.addWidget(self.tl_cutmin)
        
        self.tl_cutmax = QtGui.QLabel(self) 
        self.tl_cutmax.setText('Upper Cutoff Value: ')
        vbox1.addWidget(self.tl_cutmax)

        sizer1.setLayout(vbox1)
        vbox.addWidget(sizer1)
                
        hbox2 = QtGui.QHBoxLayout()
        self.button_ok = QtGui.QPushButton('Accept')
        self.button_ok.clicked.connect(self.OnAccept)
        self.button_ok.setEnabled(False)
        hbox2.addWidget(self.button_ok)
                
        button_cancel = QtGui.QPushButton('Cancel')
        button_cancel.clicked.connect(self.close)
        hbox2.addWidget(button_cancel)
        
        vbox.addLayout(hbox2)
        
        self.setLayout(vbox)
        
        if len(self.anlz.original_svd_maps4D) == 0:
            self.draw_histogram()
        else:
            self.draw_histogram4D()



        
#----------------------------------------------------------------------        
    def draw_histogram(self):
        
     
        fig = self.histfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        self.axes = fig.gca()

        #target_svd_maps;target_pcafit_maps
        if self.parent.page4.showraw == True:
            self.histogram = self.anlz.original_svd_maps
        else:
            self.histogram = self.anlz.original_fit_maps
            
        
        histdata = np.reshape(self.histogram, (self.anlz.stack.n_cols*self.anlz.stack.n_rows*self.anlz.n_target_spectra), order='F')
        
        self.n, self.bins, patches = self.axes.hist(histdata, 200, normed=1, facecolor='green', alpha=0.75)
        
        self.axes.set_xlabel('Thickness per Pixel in Spectral Maps')
        self.axes.set_ylabel('Percentage of Pixels')
        
        self.HistogramPanel.draw()
        

#----------------------------------------------------------------------        
    def draw_histogram4D(self):
        
     
        fig = self.histfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        self.axes = fig.gca()

        #target_svd_maps;target_pcafit_maps
        if self.parent.page4.showraw == True:
            self.histogram = self.anlz.original_svd_maps4D
        else:
            self.histogram = self.anlz.original_fit_maps4D
            
        
        histdata = np.reshape(self.histogram, (self.anlz.stack.n_cols*self.anlz.stack.n_rows*self.anlz.n_target_spectra*self.anlz.stack.n_theta), order='F')
        
        self.n, self.bins, patches = self.axes.hist(histdata, 200, normed=1, facecolor='green', alpha=0.75)
        
        self.axes.set_xlabel('Thickness per Pixel in Spectral Maps')
        self.axes.set_ylabel('Percentage of Pixels')
        
        self.HistogramPanel.draw()        

#----------------------------------------------------------------------        
    def OnClick(self, evt):
        
        x1 = evt.xdata
        
        if x1 == None:
            return
        
        if self.limit == 1:
            self.tl_cutmin.setText('Lower Cutoff Value: '+str(x1))
        
            self.histmin = x1
            
        else:
            self.tl_cutmax.setText('Upper Cutoff Value: '+str(x1))
        
            self.histmax = x1
            
        self.button_ok.setEnabled(True)


#----------------------------------------------------------------------          
    def OnRb_limit(self, enabled):
        
        state = enabled      
      
        if state:
            self.limit = 1
        else:        
            self.limit = 2
        
        
#----------------------------------------------------------------------        
    def OnAccept(self, evt):
        
        if self.parent.page4.showraw == True:
            self.anlz.svd_map_threshold(self.histmin, self.histmax, svd=True)
        else:
            self.anlz.svd_map_threshold(self.histmin, self.histmax, pca=True)
        self.parent.page4.loadTargetMap()
        self.close()


#---------------------------------------------------------------------- 
class SaveWinP4(QtGui.QDialog):

    def __init__(self, parent):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        
        self.resize(400, 300)
        self.setWindowTitle('Save')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)


        self.com = self.parent.common          
        
        path, ext = os.path.splitext(self.com.filename) 
        ext = ext[1:].lower()   
        suffix = "." + ext
        path, fn = os.path.split(self.com.filename)
        filename = fn[:-len(suffix)]
        
        self.path = self.com.path
        self.filename = filename
                          
        
        vboxtop = QtGui.QVBoxLayout()
        vboxtop.setContentsMargins(20,20,20,20)
        
        gridtop = QtGui.QGridLayout()
        gridtop.setVerticalSpacing(20)
    
        fontb = QtGui.QFont()
        fontb.setBold(True)            
        
        st1 = QtGui.QLabel(self)
        st1.setText('Save')
        st1.setFont(fontb)
        st2 = QtGui.QLabel(self)
        st2.setText('.pdf')
        st2.setFont(fontb)
        st3 = QtGui.QLabel(self)
        st3.setText('.png')
        st3.setFont(fontb)
        st4 = QtGui.QLabel(self)
        st4.setText('.svg')
        st4.setFont(fontb)
        st5 = QtGui.QLabel(self)
        st5.setText('.csv')
        st5.setFont(fontb)        
        st8 = QtGui.QLabel(self)
        st8.setText('.tif (data)')
        st8.setFont(fontb)          
        
        st6 = QtGui.QLabel(self)
        st6.setText('_spectrum')
        
        self.cb11 = QtGui.QCheckBox('', self)
        self.cb11.setChecked(True)
        self.cb12 = QtGui.QCheckBox('', self)
        self.cb13 = QtGui.QCheckBox('', self)   
        self.cb14 = QtGui.QCheckBox('', self)
        
        st7 = QtGui.QLabel(self)
        st7.setText('_images')
        
        self.cb21 = QtGui.QCheckBox('', self)
        self.cb21.setChecked(True)
        self.cb22 = QtGui.QCheckBox('', self)
        self.cb23 = QtGui.QCheckBox('', self)
        self.cb24 = QtGui.QCheckBox('', self)
        

        
        gridtop.addWidget(st1, 0, 0)
        gridtop.addWidget(st2, 0, 1)
        gridtop.addWidget(st3, 0, 2)
        gridtop.addWidget(st4, 0, 3)
        gridtop.addWidget(st5, 0, 4)
        gridtop.addWidget(st8, 0, 5)
                                
        gridtop.addWidget(st6, 1, 0)
        gridtop.addWidget(self.cb11, 1, 1)
        gridtop.addWidget(self.cb12, 1, 2)  
        gridtop.addWidget(self.cb13, 1, 3)  
        gridtop.addWidget(self.cb14, 1, 4)           
  
        gridtop.addWidget(st7, 2, 0)
        gridtop.addWidget(self.cb21, 2, 1)
        gridtop.addWidget(self.cb22, 2, 2) 
        gridtop.addWidget(self.cb23, 2, 3)  
        gridtop.addWidget(self.cb24, 2, 5) 
        

                
        vboxtop.addStretch(0.5)
        vboxtop.addLayout(gridtop)
        vboxtop.addStretch(1)
        
         
        hbox0 = QtGui.QHBoxLayout()
         
        stf = QtGui.QLabel(self)
        stf.setText('Filename:\t')
        self.tc_savefn = QtGui.QLineEdit(self)
        self.tc_savefn.setText(self.filename)

        hbox0.addWidget(stf)
        hbox0.addWidget(self.tc_savefn)         
                  
        hbox1 = QtGui.QHBoxLayout()
                 
        stp = QtGui.QLabel(self)
        stp.setText('Path:  \t')
        self.tc_savepath = QtGui.QLineEdit(self)
        self.tc_savepath.setReadOnly(True)
        self.tc_savepath.setText(self.path)
        self.tc_savepath.setMinimumWidth(100)
        hbox1.addWidget(stp)
        hbox1.addWidget(self.tc_savepath)  
         
        button_path = QtGui.QPushButton('Browse...')
        button_path.clicked.connect(self.OnBrowseDir)
        hbox1.addWidget(button_path)
         
         
        hbox2 = QtGui.QHBoxLayout()
        button_save = QtGui.QPushButton('Save')
        button_save.clicked.connect(self.OnSave)
        hbox2.addWidget(button_save)
         
        button_cancel = QtGui.QPushButton('Cancel')
        button_cancel.clicked.connect(self.close)
        hbox2.addWidget(button_cancel)
        
        vboxtop.addLayout(hbox0)
        vboxtop.addLayout(hbox1)
        vboxtop.addStretch(1.0)
        vboxtop.addLayout(hbox2)

        
        
        self.setLayout(vboxtop)
        
#----------------------------------------------------------------------        
    def OnBrowseDir(self, evt):
        
        directory = QtGui.QFileDialog.getExistingDirectory(self, "Choose a directory", self.path, QtGui.QFileDialog.ShowDirsOnly|QtGui.QFileDialog.ReadOnly)       
                                                        
        
       
        if directory == '':
            return
                 
        directory = str(directory)
        self.com.path = directory
                    
        self.path = directory
        
        self.tc_savepath.setText(self.path)
            
            
                
#----------------------------------------------------------------------        
    def OnSave(self, evt):
        
        self.filename = str(self.tc_savefn.text())
        
        sp_pdf = self.cb11.isChecked()
        sp_png = self.cb12.isChecked()
        sp_svg = self.cb13.isChecked()
        sp_csv = self.cb14.isChecked()
        im_pdf = self.cb21.isChecked()
        im_png = self.cb22.isChecked()
        im_svg = self.cb23.isChecked()
        im_tif = self.cb24.isChecked()

        
        self.close() 
        self.parent.page4.Save(self.filename, self.path,
                                         spec_png = sp_png, 
                                         spec_pdf = sp_pdf, 
                                         spec_svg = sp_svg,
                                         spec_csv = sp_csv,
                                         img_png = im_png, 
                                         img_pdf = im_pdf,
                                         img_svg = im_svg,
                                         img_tif = im_tif)


    
    
""" ------------------------------------------------------------------------------------------------"""
class PageCluster(QtGui.QWidget):
    def __init__(self, common, data_struct, stack, anlz):
        super(PageCluster, self).__init__()

        self.initUI(common, data_struct, stack, anlz)
        
#----------------------------------------------------------------------          
    def initUI(self, common, data_struct, stack, anlz): 
        
        self.data_struct = data_struct
        self.stk = stack
        self.com = common
        self.anlz = anlz
              
        
        self.selcluster = 1
        self.numclusters = 0
        self.init_nclusters = 5
        self.wo_1st_pca = 0
        self.sigma_split = 0
        self.showallspectra = 0
        self.pcscalingfactor = 0.0
                
        self.MakeColorTable()             
         
        
        #panel 1
        sizer1 = QtGui.QGroupBox('Cluster analysis')
        vbox1 = QtGui.QVBoxLayout()
       
        
        self.button_calcca = QtGui.QPushButton('Calculate Clusters')
        self.button_calcca.clicked.connect( self.OnCalcClusters)   
        self.button_calcca.setEnabled(False)
        vbox1.addWidget(self.button_calcca)
        self.button_scatterplots = QtGui.QPushButton('Show scatter plots...')
        self.button_scatterplots.clicked.connect( self.OnShowScatterplots)
        self.button_scatterplots.setEnabled(False)
        self.button_savecluster = QtGui.QPushButton('Save CA Results...')
        self.button_savecluster.clicked.connect( self.OnSave)
        self.button_savecluster.setEnabled(False)
        
        vbox1.addStretch(1)
        
                
        hbox11 = QtGui.QHBoxLayout()
        text1 = QtGui.QLabel(self)
        text1.setText('Number of clusters')
        hbox11.addWidget(text1)
        self.nclusterspin = QtGui.QSpinBox()
        self.nclusterspin.setRange(2,20)
        self.nclusterspin.setValue(self.init_nclusters)
        self.nclusterspin.valueChanged[int].connect(self.OnNClusterspin)
        hbox11.addWidget(text1)
        hbox11.addWidget(self.nclusterspin)  
        
        vbox1.addLayout(hbox11) 
         
        hbox12 = QtGui.QHBoxLayout()
        text1a = QtGui.QLabel(self)
        text1a.setText("Number of clusters found")
        hbox12.addWidget(text1a)
        
        self.ntc_clusters_found = QtGui.QLabel(self)
        self.ntc_clusters_found.setText(str(self.numclusters))
        hbox12.addWidget(self.ntc_clusters_found)
  
        vbox1.addLayout(hbox12) 
          
         
        hbox13 = QtGui.QHBoxLayout()
        self.remove1stpcacb = QtGui.QCheckBox('Reduce thickness effects', self)
        self.remove1stpcacb.stateChanged.connect(self.OnRemove1stpca)
        hbox13.addWidget(self.remove1stpcacb)
        
        vbox1.addLayout(hbox13) 
 
        hbox14 = QtGui.QHBoxLayout()
        self.cb_splitclusters = QtGui.QCheckBox('Divide clusters with large Sigma', self)
        self.cb_splitclusters.stateChanged.connect(self.OnSplitClusters)
        hbox14.addWidget(self.cb_splitclusters)
        
        vbox1.addLayout(hbox14) 
                
        hbox14a = QtGui.QHBoxLayout()
        tc1 = QtGui.QLabel(self)
        tc1.setText("PC scaling factor")  
        hbox14a.addWidget(tc1)      
        self.ntc_pcscaling = QtGui.QLineEdit(self)
        self.ntc_pcscaling.setFixedWidth(65)
        self.ntc_pcscaling.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.ntc_pcscaling.setAlignment(QtCore.Qt.AlignRight)         
        
        self.ntc_pcscaling.setText(str(self.pcscalingfactor))
        hbox14a.addWidget(self.ntc_pcscaling) 
        hbox14a.addStretch(1)
        vbox1.addLayout(hbox14a) 

        line = QtGui.QFrame()
        line.setFrameShape(QtGui.QFrame.HLine)
        line.setFrameShadow(QtGui.QFrame.Sunken) 
        

        vbox1.addStretch(1)
        vbox1.addWidget(line) 
        vbox1.addStretch(1)     
        vbox1.addWidget(self.button_scatterplots)
        vbox1.addWidget(self.button_savecluster)
        
        
        sizer1.setLayout(vbox1)
       
                
                 
        #panel 2        
        vbox2 = QtGui.QVBoxLayout()
         
        tc_clustercomp = QtGui.QLabel(self)
        tc_clustercomp.setText("Composite cluster image")  
        vbox2.addWidget(tc_clustercomp)      
          
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
   
        self.clusterimgfig = Figure((PlotH, PlotH))
        self.ClusterImagePan = FigureCanvas(self.clusterimgfig)
        self.ClusterImagePan.mpl_connect('button_press_event', self.OnPointClusterImage)
        self.ClusterImagePan.setParent(self)
        fbox.addWidget(self.ClusterImagePan)
        frame.setLayout(fbox)
        vbox2.addWidget(frame)         
         
        #panel 3 
        vbox3 = QtGui.QVBoxLayout()
        fgs = QtGui.QGridLayout()
         
        self.tc_cluster = QtGui.QLabel(self)
        self.tc_cluster.setText("Cluster ")
        fgs.addWidget(self.tc_cluster, 0, 0, QtCore .Qt. AlignLeft)

        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()

        self.clusterindvimgfig = Figure((PlotH*0.73, PlotH*0.73))
        self.ClusterIndvImagePan = FigureCanvas(self.clusterindvimgfig)
        self.ClusterIndvImagePan.setParent(self)
        fbox.addWidget(self.ClusterIndvImagePan)
        frame.setLayout(fbox)   
        fgs.addWidget(frame, 1, 0, QtCore .Qt. AlignLeft)
             
        self.slidershow = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slidershow.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slidershow.setEnabled(False)
        self.slidershow.valueChanged[int].connect(self.OnClusterScroll)
        self.slidershow.setRange(1, 20)
        fgs.addWidget(self.slidershow, 1, 1, QtCore .Qt. AlignLeft)
            
         
        text3 = QtGui.QLabel(self)
        text3.setText('Cluster Error Map')
        fgs.addWidget(text3, 0, 2, QtCore .Qt. AlignLeft)
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
        self.clusterdistmapfig = Figure((PlotH*0.73, PlotH*0.73))
        self.ClusterDistMapPan = FigureCanvas(self.clusterdistmapfig)
        self.ClusterDistMapPan.setParent(self)
        fbox.addWidget(self.ClusterDistMapPan)
        frame.setLayout(fbox)
        fgs.addWidget(frame, 1, 2, QtCore .Qt. AlignLeft)          
 
        vbox3.addLayout(fgs)
 
         
        

        #panel 4 
        vbox4 = QtGui.QVBoxLayout()
         
        self.tc_clustersp = QtGui.QLabel(self)
        self.tc_clustersp.setText("Cluster spectrum")   
        vbox4.addWidget(self.tc_clustersp)      
 
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()

        self.clusterspecfig = Figure((PlotW, PlotH))
        self.ClusterSpecPan = FigureCanvas(self.clusterspecfig)
        self.ClusterSpecPan.setParent(self)
 
        fbox.addWidget(self.ClusterSpecPan)
        frame.setLayout(fbox)
        vbox4.addWidget(frame) 
         

        #panel 5
        sizer5 = QtGui.QGroupBox('Display')
        vbox5 = QtGui.QVBoxLayout()         
         
        hbox51 = QtGui.QHBoxLayout()
        self.showallspectracb = QtGui.QCheckBox('Show all spectra', self)
        self.showallspectracb.stateChanged.connect(self.OnShowallspectra)
        hbox51.addWidget(self.showallspectracb)
        
        vbox5.addLayout(hbox51)         
        
        sizer5.setLayout(vbox5)
        
        
        
        vboxtop = QtGui.QVBoxLayout()
        hboxtopL = QtGui.QHBoxLayout()
        vboxtopL = QtGui.QVBoxLayout()
        vboxtopL.addWidget(sizer1)
        vboxtopL.addWidget(sizer5)
        hboxtopL.addLayout(vboxtopL)
        hboxtopL.addStretch(1)
        
        gridsizertop = QtGui.QGridLayout()
        gridsizertop.setContentsMargins(15,0,0,0)

        
        gridsizertop.addLayout(hboxtopL, 0, 0, QtCore .Qt. AlignLeft)
        gridsizertop.addLayout(vbox2, 1, 0, QtCore .Qt. AlignLeft)
        
        gridsizertop.addLayout(vbox3, 0, 1, QtCore .Qt. AlignLeft)
        gridsizertop.addLayout(vbox4, 1, 1, QtCore .Qt. AlignLeft)
        
        vboxtop.addStretch(1)
        vboxtop.addLayout(gridsizertop)
        vboxtop.addStretch(1)
        self.setLayout(vboxtop)
                
        
#----------------------------------------------------------------------     
    def MakeColorTable(self):
        self.maxclcolors = 11
        colors_i = np.linspace(0,self.maxclcolors,self.maxclcolors+1)

       
#         self.colors=['#0000FF','#FF0000','#DFE32D','#36F200','#B366FF',
#                 '#FF470A','#33FFFF','#006600','#CCCC99','#993300',
#                 '#000000']

#         self.colors=['#D98619','#ED2024','#98CC31','#861F78','#007FFF',
#                 '#6FDBDB','#5C3F32','#FF6EC7','#CCCC99','#993300',
#                 '#000000']
        
        self.colors=['#007FFF','#ED2024','#98CC31','#861F78','#D98619',
                '#6FDBDB','#5C3F32','#FF6EC7','#CCCC99','#993300',
                '#000000']
        
        
        self.clusterclrmap1=matplotlib.colors.LinearSegmentedColormap.from_list('clustercm',self.colors)
     
        self.bnorm1 = matplotlib.colors.BoundaryNorm(colors_i, self.clusterclrmap1.N)
        
        colors_i = np.linspace(0,self.maxclcolors+2,self.maxclcolors+3)
        
        #use black color for clusters > maxclcolors, the other 2 colors are for background      
#         colors2=['#0000FF','#FF0000','#DFE32D','#36F200','#B366FF',
#                 '#FF470A','#33FFFF','#006600','#CCCC99','#993300',
#                 '#000000','#FFFFFF','#EEEEEE']

#         colors2=['#D98619','#ED2024','#98CC31','#861F78','#007FFF',
#                 '#6FDBDB','#5C3F32','#FF6EC7','#CCCC99','#993300',
#                 '#000000','#FFFFFF','#EEEEEE']
        
        colors2=['#007FFF','#ED2024','#98CC31','#861F78','#D98619',
                '#6FDBDB','#5C3F32','#FF6EC7','#CCCC99','#993300',
                '#000000','#FFFFFF','#EEEEEE']
                
        self.clusterclrmap2=matplotlib.colors.LinearSegmentedColormap.from_list('clustercm2',colors2)
     
        self.bnorm2 = matplotlib.colors.BoundaryNorm(colors_i, self.clusterclrmap2.N)
        


        
        
#----------------------------------------------------------------------
    def OnCalcClusters(self, event):
       
        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        self.calcclusters = False  
        

        try: 
            
            value = self.ntc_pcscaling.text()
            #try:
            self.pcscalingfactor = float(value)
#             except:
#                 self.pcscalingfactor = 0.0
#                 self.ntc_pcscaling.setText('0.0')
                
            
            self.CalcClusters()
            
            self.calcclusters = True
            
            self.selcluster = 1
            self.slidershow.setValue(self.selcluster)
            self.slidershow.setMaximum(self.numclusters)
            
            self.showClusterImage()
            self.showClusterSpectrum()
            self.showIndvClusterImage()     
            self.showClusterDistanceMap()
            self.com.cluster_calculated = 1       
            QtGui.QApplication.restoreOverrideCursor()
            
        except:
            self.com.cluster_calculated = 0
            QtGui.QApplication.restoreOverrideCursor()     
            
        self.window().refresh_widgets()
            
#----------------------------------------------------------------------        
    def OnNClusterspin(self, value):
        num = value
        self.init_nclusters = num
                       
 
#----------------------------------------------------------------------        
    def OnClusterScroll(self, value):
        sel = value
        self.selcluster = sel
        if self.com.cluster_calculated == 1:
            self.showClusterSpectrum()
            self.showIndvClusterImage()
            
#----------------------------------------------------------------------            
    def OnClusterSpinUp(self, event):
        if (self.com.cluster_calculated == 1) and (self.selcluster > 1):
            self.selcluster = self.selcluster - 1
            self.slidershow.setValue(self.selcluster)

            self.showClusterSpectrum()
            self.showIndvClusterImage()
            
#----------------------------------------------------------------------            
    def OnClusterSpinDown(self, event):
        if (self.com.cluster_calculated == 1) and (self.selcluster < self.numclusters):
            self.selcluster = self.selcluster + 1
            self.slidershow.setValue(self.selcluster) 
            
            self.showClusterSpectrum()
            self.showIndvClusterImage()
                       
#----------------------------------------------------------------------  
    def OnPointClusterImage(self, evt):
        x = evt.xdata
        y = evt.ydata
        

        if self.com.cluster_calculated == 1:   
            try:  
                self.ix = int(np.floor(y))           
                self.iy = int(np.floor(x))  
                        
                if self.ix<0 :
                    self.ix=0
                if self.ix>self.stk.n_cols :
                    self.ix=self.stk.n_cols
                if self.iy<0 :
                    self.iy=0
                if self.iy>self.stk.n_rows :
                    self.iy=self.stk.n_rows 
                    
                    
                self.selcluster = self.anlz.cluster_indices[self.ix,self.iy] + 1
                
                self.slidershow.setValue(self.selcluster)
                
                self.showClusterSpectrum()
                self.showIndvClusterImage()
            except:
                pass
            
                       
#----------------------------------------------------------------------
    def CalcClusters(self):
        
        nclusters = self.anlz.calculate_clusters(self.init_nclusters, 
                                                 self.wo_1st_pca, 
                                                 self.sigma_split, 
                                                 pcscalingfactor = self.pcscalingfactor)
        #nclusters = self.anlz.calculate_clusters_kmeansangle(self.init_nclusters, self.wo_1st_pca, self.sigma_split)
        self.numclusters = nclusters
        self.ntc_clusters_found.setText(str(self.numclusters))
    
        
        
        
#----------------------------------------------------------------------     
#Show composite cluster image  
    def showClusterImage(self):       
        
        self.clusterimage = self.anlz.cluster_indices  
        
        #print self.selpca

        fig = self.clusterimgfig
        fig.clf()
        fig.add_axes(((0.0,0.0,1.0,1.0)))
        axes = fig.gca()
        
        
        
        im = axes.imshow(np.rot90(self.clusterimage), cmap=self.clusterclrmap1, norm=self.bnorm1)
        axes.axis("off")
        #cbar = axes.figure.colorbar(im)         
        self.ClusterImagePan.draw()


      
#----------------------------------------------------------------------     
#Show composite cluster image  
    def showIndvClusterImage(self):
        
        indvclusterimage = np.zeros((self.anlz.stack.n_cols, self.anlz.stack.n_rows))+20.      
        ind = np.where(self.anlz.cluster_indices == self.selcluster-1)    
        colorcl = min(self.selcluster-1,self.maxclcolors-1)
        indvclusterimage[ind] = colorcl

        fig = self.clusterindvimgfig
        fig.clf()
        fig.add_axes(((0.0,0.0,1.0,1.0)))
        axes = fig.gca()
            
        
        im = axes.imshow(np.rot90(indvclusterimage), cmap=self.clusterclrmap2, norm=self.bnorm2)
        axes.axis("off")

        self.ClusterIndvImagePan.draw()
        
        self.tc_cluster.setText("Cluster " + str(self.selcluster))
         
   
#----------------------------------------------------------------------     
    def showClusterSpectrum(self):
        
         
        clusterspectrum = self.anlz.clusterspectra[self.selcluster-1, ]

       
        fig = self.clusterspecfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        
        
        if self.showallspectra == 0:
            if self.selcluster >= self.maxclcolors:
                clcolor = self.colors[self.maxclcolors-1]
            else:
                clcolor = self.colors[self.selcluster-1]
            specplot = axes.plot(self.anlz.stack.ev,clusterspectrum, color = clcolor)
            self.tc_clustersp.setText("Cluster " + str(self.selcluster)+ " spectrum"  ) 

        else:
            #Show all spectra
            for i in range(1, self.numclusters+1):
                if i >= self.maxclcolors:
                    clcolor = self.colors[self.maxclcolors-1]
                else:
                    clcolor = self.colors[i-1]  
                clusterspectrum = self.anlz.clusterspectra[i-1, ]/np.amax(self.anlz.clusterspectra[i-1, ])           
                specplot = axes.plot(self.anlz.stack.ev,clusterspectrum, color = clcolor)
            
            self.tc_clustersp.setText(" Normalized Cluster spectra"  ) 
            
        
        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Optical Density')
        
        
        self.ClusterSpecPan.draw()
        
        
        
        
#----------------------------------------------------------------------     
    def showClusterDistanceMap(self):       
        
        mapimage = self.anlz.cluster_distances
        
        #print self.selpca

        fig = self.clusterdistmapfig
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        axes = fig.gca()
        
#         divider = make_axes_locatable(axes)
#         axcb = divider.new_horizontal(size="3%", pad=0.03)  
#         fig.add_axes(axcb)  
#         axes.set_position([0.03,0.03,0.8,0.94])
               
        
        
        im = axes.imshow(np.rot90(mapimage), cmap=matplotlib.cm.get_cmap('gray'))
        
        #cbar = axes.figure.colorbar(im, orientation='vertical',cax=axcb) 
        
        axes.axis("off")
       
        self.ClusterDistMapPan.draw()
        
#----------------------------------------------------------------------           
    def OnRemove1stpca(self, state):
        if state == QtCore.Qt.Checked:
            self.wo_1st_pca = 1
        else: self.wo_1st_pca = 0

#----------------------------------------------------------------------           
    def OnSplitClusters(self, state):
        if state == QtCore.Qt.Checked:
            self.sigma_split = 1
        else: self.sigma_split = 0        

#----------------------------------------------------------------------           
    def OnShowallspectra(self, state):
        if state == QtCore.Qt.Checked:
            self.showallspectra = 1
        else: self.showallspectra = 0     
        
        if self.com.cluster_calculated == 1:   
            self.showClusterSpectrum()
                
#----------------------------------------------------------------------    
    def OnSave(self, event):     
               
        savewin = SaveWinP3(self.window())
        savewin.show()
        
        
#----------------------------------------------------------------------    
    def Save(self, filename, path, spec_png = True, spec_pdf = False, spec_svg = False, spec_csv = False,
             img_png = True, img_pdf = False, img_svg = False, img_tif = False,
             indimgs_png = True, indimgs_pdf = False, indimgs_svg = False, indimgs_tif = False,
             scatt_png = True, scatt_pdf = False, scatt_svg = False): 
        
        self.SaveFileName = os.path.join(path,filename)
        
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas    
        matplotlib.rcParams['pdf.fonttype'] = 42
   
        try: 
            
            if img_png:
                ext = 'png'
            
                fig = matplotlib.figure.Figure(figsize = (float(self.stk.n_rows)/10, float(self.stk.n_cols)/10))
                canvas = FigureCanvas(fig)
                fig.clf()
                fig.add_axes((0.0,0.0,1.0,1.0))
                axes = fig.gca()         
        
                im = axes.imshow(np.rot90(self.clusterimage), cmap=self.clusterclrmap1, norm=self.bnorm1)
                axes.axis("off")
                
                fileName_caimg = self.SaveFileName+"_CAcimg."+ext       
                fig.savefig(fileName_caimg, dpi=ImgDpi, pad_inches = 0.0)
                
            
            if img_pdf:
                ext = 'pdf'
                suffix = "." + ext
            
                fig = matplotlib.figure.Figure(figsize = (float(self.stk.n_rows)/30, float(self.stk.n_cols)/30))
                canvas = FigureCanvas(fig)
                fig.clf()
                fig.add_axes((0.0,0.0,1.0,1.0))
                axes = fig.gca() 
        
                im = axes.imshow(np.rot90(self.clusterimage), cmap=self.clusterclrmap1, norm=self.bnorm1)
                axes.axis("off")
                
                fileName_caimg = self.SaveFileName+"_CAcimg."+ext       
                fig.savefig(fileName_caimg, dpi=300, pad_inches = 0.0)
                            
            if img_svg:
                ext = 'svg'
                suffix = "." + ext
            
                fig = matplotlib.figure.Figure(figsize = (float(self.stk.n_rows)/30, float(self.stk.n_cols)/30))
                canvas = FigureCanvas(fig)
                fig.clf()
                fig.add_axes((0.0,0.0,1.0,1.0))
                axes = fig.gca() 
        
                im = axes.imshow(np.rot90(self.clusterimage), cmap=self.clusterclrmap1, norm=self.bnorm1)
                axes.axis("off")
                
                fileName_caimg = self.SaveFileName+"_CAcimg."+ext       
                fig.savefig(fileName_caimg, dpi=300, pad_inches = 0.0)
                            
            if img_tif:
                fileName_caimg = self.SaveFileName+"_CAcimg.tif"   
                img1 = Image.fromarray(self.clusterimage)
                img1.save(fileName_caimg)        
                
                  
            ext = 'png'
            suffix = "." + ext
                
            if indimgs_png:
                for i in range (self.numclusters):
              
                    indvclusterimage = np.zeros((self.anlz.stack.n_cols, self.anlz.stack.n_rows))+20.      
                    ind = np.where(self.anlz.cluster_indices == i)    
                    colorcl = min(i,9)
                    indvclusterimage[ind] = colorcl

                    fig = matplotlib.figure.Figure(figsize =(float(self.stk.n_rows)/10, float(self.stk.n_cols)/10))
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.0,0.0,1.0,1.0))
                    axes = fig.gca()            
                    im = axes.imshow(np.rot90(indvclusterimage), cmap=self.clusterclrmap2, norm=self.bnorm2)
                    axes.axis("off")
                   
                    fileName_img = self.SaveFileName+"_CAimg_" +str(i+1)+"."+ext               
                    fig.savefig(fileName_img, dpi=ImgDpi, pad_inches = 0.0)
                
            if spec_png:
                for i in range (self.numclusters):
                   
                    clusterspectrum = self.anlz.clusterspectra[i, ]
                    fig = matplotlib.figure.Figure(figsize =(PlotW, PlotH))
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.15,0.15,0.75,0.75))
                    axes = fig.gca()
                    if i >= self.maxclcolors:
                        clcolor = self.colors[self.maxclcolors-1]
                    else:
                        clcolor = self.colors[i]
        
                    specplot = axes.plot(self.anlz.stack.ev,clusterspectrum, color = clcolor)
        
                    axes.set_xlabel('Photon Energy [eV]')
                    axes.set_ylabel('Optical Density')

                    fileName_spec = self.SaveFileName+"_CAspectrum_" +str(i+1)+"."+ext
                    fig.savefig(fileName_spec)   
                    
                #Save all spectra in one plot
                fig = matplotlib.figure.Figure(figsize =(PlotW, PlotH))
                canvas = FigureCanvas(fig)
                fig.add_axes((0.15,0.15,0.75,0.75))
                axes = fig.gca()                
                               
                
                for i in range(1, self.numclusters+1):
                   
                    clusterspectrum = self.anlz.clusterspectra[i-1, ]/np.amax(self.anlz.clusterspectra[i-1, ])

                    if i >= self.maxclcolors:
                        clcolor = self.colors[self.maxclcolors-1]
                    else:
                        clcolor = self.colors[i-1]
        
                    specplot = axes.plot(self.anlz.stack.ev,clusterspectrum, color = clcolor)
        
                    axes.set_xlabel('Photon Energy [eV]')
                    axes.set_ylabel('Optical Density')

                fileName_spec = self.SaveFileName+"_CAspectra"+"."+ext
                fig.savefig(fileName_spec)  
                    
            if spec_csv:
                for i in range (self.numclusters):
                    clusterspectrum = self.anlz.clusterspectra[i, ]
                    fileName_spec = self.SaveFileName+"_CAspectrum_" +str(i+1)+".csv"
                    cname = 'CAspectrum_' +str(i+1)
                    self.stk.write_csv(fileName_spec, self.anlz.stack.ev, clusterspectrum, cname=cname)
                                                     
                
            ext = 'pdf'
            suffix = "." + ext
                
            if indimgs_pdf:
                for i in range (self.numclusters):
              
                    indvclusterimage = np.zeros((self.anlz.stack.n_cols, self.anlz.stack.n_rows))+20.      
                    ind = np.where(self.anlz.cluster_indices == i)    
                    colorcl = min(i,9)
                    indvclusterimage[ind] = colorcl

                    fig = matplotlib.figure.Figure(figsize =(float(self.stk.n_rows)/30, float(self.stk.n_cols)/30))
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.0,0.0,1.0,1.0))
                    axes = fig.gca()       
                    im = axes.imshow(np.rot90(indvclusterimage), cmap=self.clusterclrmap2, norm=self.bnorm2)
                    axes.axis("off")
                   
                    fileName_img = self.SaveFileName+"_CAimg_" +str(i+1)+"."+ext               
                    fig.savefig(fileName_img, dpi=300, pad_inches = 0.0)
                
            if spec_pdf:
                for i in range (self.numclusters):
                   
                    clusterspectrum = self.anlz.clusterspectra[i, ]
                    fig = matplotlib.figure.Figure(figsize =(PlotW, PlotH))
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.15,0.15,0.75,0.75))
                    axes = fig.gca()
                    if i >= self.maxclcolors:
                        clcolor = self.colors[self.maxclcolors-1]
                    else:
                        clcolor = self.colors[i]
        
                    specplot = axes.plot(self.anlz.stack.ev,clusterspectrum, color = clcolor)
        
                    axes.set_xlabel('Photon Energy [eV]')
                    axes.set_ylabel('Optical Density')

                    fileName_spec = self.SaveFileName+"_CAspectrum_" +str(i+1)+"."+ext
                    fig.savefig(fileName_spec) 
                    
                #Save all spectra in one plot
                fig = matplotlib.figure.Figure(figsize =(PlotW, PlotH))
                canvas = FigureCanvas(fig)
                fig.add_axes((0.15,0.15,0.75,0.75))
                axes = fig.gca()                
                
                for i in range(1, self.numclusters+1):
                   
                    clusterspectrum = self.anlz.clusterspectra[i-1, ]/np.amax(self.anlz.clusterspectra[i-1, ])

                    if i >= self.maxclcolors:
                        clcolor = self.colors[self.maxclcolors-1]
                    else:
                        clcolor = self.colors[i-1]
        
                    specplot = axes.plot(self.anlz.stack.ev,clusterspectrum, color = clcolor)
        
                    axes.set_xlabel('Photon Energy [eV]')
                    axes.set_ylabel('Optical Density')

                fileName_spec = self.SaveFileName+"_CAspectra"+"."+ext
                fig.savefig(fileName_spec)  
 
            ext = 'svg'
            suffix = "." + ext
                
            if indimgs_svg:
                for i in range (self.numclusters):
              
                    indvclusterimage = np.zeros((self.anlz.stack.n_cols, self.anlz.stack.n_rows))+20.      
                    ind = np.where(self.anlz.cluster_indices == i)    
                    colorcl = min(i,9)
                    indvclusterimage[ind] = colorcl

                    fig = matplotlib.figure.Figure(figsize =(float(self.stk.n_rows)/30, float(self.stk.n_cols)/30))
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.0,0.0,1.0,1.0))
                    axes = fig.gca()       
                    im = axes.imshow(np.rot90(indvclusterimage), cmap=self.clusterclrmap2, norm=self.bnorm2)
                    axes.axis("off")
                   
                    fileName_img = self.SaveFileName+"_CAimg_" +str(i+1)+"."+ext               
                    fig.savefig(fileName_img, dpi=300, pad_inches = 0.0)
                
            if spec_svg:
                for i in range (self.numclusters):
                   
                    clusterspectrum = self.anlz.clusterspectra[i, ]
                    fig = matplotlib.figure.Figure(figsize =(PlotW, PlotH))
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.15,0.15,0.75,0.75))
                    axes = fig.gca()
                    if i >= self.maxclcolors:
                        clcolor = self.colors[self.maxclcolors-1]
                    else:
                        clcolor = self.colors[i]
        
                    specplot = axes.plot(self.anlz.stack.ev,clusterspectrum, color = clcolor)
        
                    axes.set_xlabel('Photon Energy [eV]')
                    axes.set_ylabel('Optical Density')

                    fileName_spec = self.SaveFileName+"_CAspectrum_" +str(i+1)+"."+ext
                    fig.savefig(fileName_spec) 
                    
                #Save all spectra in one plot
                fig = matplotlib.figure.Figure(figsize =(PlotW, PlotH))
                canvas = FigureCanvas(fig)
                fig.add_axes((0.15,0.15,0.75,0.75))
                axes = fig.gca()                
                
                for i in range(1, self.numclusters+1):
                   
                    clusterspectrum = self.anlz.clusterspectra[i-1, ]/np.amax(self.anlz.clusterspectra[i-1, ])

                    if i >= self.maxclcolors:
                        clcolor = self.colors[self.maxclcolors-1]
                    else:
                        clcolor = self.colors[i-1]
        
                    specplot = axes.plot(self.anlz.stack.ev,clusterspectrum, color = clcolor)
        
                    axes.set_xlabel('Photon Energy [eV]')
                    axes.set_ylabel('Optical Density')

                fileName_spec = self.SaveFileName+"_CAspectra"+"."+ext
                fig.savefig(fileName_spec)  


            if indimgs_tif:
                for i in range (self.numclusters):
              
                    indvclusterimage = np.zeros((self.anlz.stack.n_cols, self.anlz.stack.n_rows))+20.      
                    ind = np.where(self.anlz.cluster_indices == i)    
                    colorcl = min(i,9)
                    indvclusterimage[ind] = colorcl

                    fileName_img = self.SaveFileName+"_CAimg_" +str(i+1)+".tif" 
                    img1 = Image.fromarray(indvclusterimage)
                    img1.save(fileName_img) 
                        
                    
            if scatt_png:
                self.SaveScatt(png_pdf = 1)
            if scatt_pdf:
                self.SaveScatt(png_pdf = 2)      
            if scatt_svg:
                self.SaveScatt(png_pdf = 3)              
            
        except IOError, e:
            if e.strerror:
                err = e.strerror 
            else: 
                err = e 
   
            QtGui.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err)

  
  
#----------------------------------------------------------------------    
#If png_pdg = 1 save png, if =2 save pdf, if =3 save svg
    def SaveScatt(self, png_pdf = 1): 
          
        od_reduced = self.anlz.pcaimages[:,:,0:self.anlz.numsigpca]        
        od_reduced = np.reshape(od_reduced, (self.stk.n_cols*self.stk.n_rows,self.anlz.numsigpca), order='F')

        clindices = self.anlz.cluster_indices
        clindices = np.reshape(clindices, (self.stk.n_cols*self.stk.n_rows), order='F')
        
        path, ext = os.path.splitext(self.SaveFileName) 
        ext = ext[1:].lower() 
        
        if png_pdf == 1:
            ext = 'png'
        elif png_pdf == 2:
            ext = 'pdf'
        elif png_pdf == 3:
            ext = 'svg'
        
        
   
        try: 
            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
            
            suffix = "." + ext
            
            nplots = 0
            for ip in range(self.anlz.numsigpca):
                for jp in range(self.anlz.numsigpca):
                    if jp >= (ip+1):
                        nplots = nplots+1    
            nplotsrows = np.ceil(nplots/2)    
            
            plotsize = 2.5    
            
            from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas  
            
            if nplots > 1 :
                fig = matplotlib.figure.Figure(figsize =(6.0,plotsize*nplotsrows))
                fig.subplots_adjust(wspace = 0.4, hspace = 0.4)
            else:
                fig = matplotlib.figure.Figure(figsize =(3.0,2.5))
                fig.subplots_adjust(bottom = 0.2, left = 0.2)
                
            canvas = FigureCanvas(fig)
            #axes = fig.gca()
            matplotlib.rcParams['font.size'] = 6   
            
                                   
            pplot = 1
            for ip in range(self.anlz.numsigpca):
                for jp in range(self.anlz.numsigpca):
                    if jp >= (ip+1):
                        

                        x_comp = od_reduced[:,ip]
                        y_comp = od_reduced[:,jp]
                        if nplots > 1 :
                            axes = fig.add_subplot(nplotsrows,2, pplot)
                        else:
                            axes = fig.add_subplot(1,1,1)
                            
                        pplot = pplot+1
                        
                        for i in range(self.numclusters):
                            thiscluster = np.where(clindices == i)
                            axes.plot(x_comp[thiscluster], y_comp[thiscluster],'.',color=self.colors[i],alpha=0.5)
                        axes.set_xlabel('Component '+str(ip+1))
                        axes.set_ylabel('Component '+str(jp+1))
                            
    
            fileName_sct = self.SaveFileName+"_CAscatterplot_" +str(i+1)+"."+ext
            matplotlib.rcParams['pdf.fonttype'] = 42
            fig.savefig(fileName_sct)
            
            QtGui.QApplication.restoreOverrideCursor()
     
            
        except IOError, e:
            QtGui.QApplication.restoreOverrideCursor()
            if e.strerror:
                err = e.strerror 
            else: 
                err = e 
   
            QtGui.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err)

            
                
#----------------------------------------------------------------------     
    def OnShowScatterplots(self, evt):   
        
        scattplwin = Scatterplots(self.window(), self.com, self.anlz)
        scattplwin.show()


#---------------------------------------------------------------------- 
class Scatterplots(QtGui.QDialog):

    def __init__(self, parent,  common, analz):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        
        self.resize(600, 470)
        self.setWindowTitle('Scatter plots')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)
                
         
         
        self.colors = self.parent.page3.colors
         
        self.com = common       
         
        self.anlz = analz
        self.numsigpca = self.anlz.numsigpca  
        self.ncols = self.anlz.stack.n_cols
        self.nrows = self.anlz.stack.n_rows
         
        self.od_reduced = self.anlz.pcaimages[:,:,0:self.numsigpca]        
        self.od_reduced = np.reshape(self.od_reduced, (self.ncols*self.nrows,self.numsigpca), order='F')
         
        self.clindices = self.anlz.cluster_indices
        self.clindices = np.reshape(self.clindices, (self.ncols*self.nrows), order='F')
        self.numclusters = self.parent.page3.numclusters
              
                      
        self.pca_y = 1
        self.pca_x = 1
         
 
                 
        vbox = QtGui.QVBoxLayout()
           
        grid1 = QtGui.QGridLayout()

        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
        self.scattplfig = Figure((5.0, 4.8))
        self.ScatterPPanel = FigureCanvas(self.scattplfig)
        self.ScatterPPanel.setParent(self)
        fbox.addWidget(self.ScatterPPanel)
        frame.setLayout(fbox)
        
        self.slidershow_y = QtGui.QSlider(QtCore.Qt.Vertical)
        self.slidershow_y.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slidershow_y.setRange(1, self.numsigpca)
        self.slidershow_y.setValue(self.pca_y)          
        self.slidershow_y.valueChanged[int].connect(self.OnSliderScroll_y)
               
        grid1.addWidget(self.slidershow_y, 0, 0)
        grid1.addWidget(frame, 0, 1)
                
         
        self.slidershow_x = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.slidershow_x.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slidershow_x.setRange(1, self.numsigpca)
        self.slidershow_x.setValue(self.pca_x)          
        self.slidershow_x.valueChanged[int].connect(self.OnSliderScroll_x)
         
        #grid1.addWidget(wx.StaticText(panel, -1, ''))
        grid1.addWidget(self.slidershow_x, 1,  1)
         
        hbox = QtGui.QVBoxLayout()
         
        button_close = QtGui.QPushButton('Close')
        button_close.clicked.connect( self.close)
                
        hbox.addStretch(1)
        hbox.addWidget(button_close)
         
        vbox.addLayout(grid1)
        vbox.addLayout(hbox)
                 
         
        self.setLayout(vbox)
         
        self.draw_scatterplot()
        
        
#----------------------------------------------------------------------        
    def OnSliderScroll_x(self, value):
        self.pca_x = value
        
        self.draw_scatterplot()

#----------------------------------------------------------------------        
    def OnSliderScroll_y(self, value):
        self.pca_y = value
        
        self.draw_scatterplot()        

      
#----------------------------------------------------------------------        
    def draw_scatterplot(self):
                
        x_comp = self.od_reduced[:,self.pca_x-1]
        y_comp = self.od_reduced[:,self.pca_y-1]
        
        fig = self.scattplfig
        fig.clf()
        axes = fig.gca()
             
        
        for i in range(self.numclusters):
            thiscluster = np.where(self.clindices == i)
            axes.plot(x_comp[thiscluster], y_comp[thiscluster],'.',color=self.colors[i],alpha=0.5)
            
        axes.set_xlabel('Component '+str(self.pca_x))
        axes.set_ylabel('Component '+str(self.pca_y))
      
        self.ScatterPPanel.draw()
       
  
     
#---------------------------------------------------------------------- 
class SaveWinP3(QtGui.QDialog):

    def __init__(self, parent):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        
        self.resize(400, 300)
        self.setWindowTitle('Save')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)


        self.com = self.parent.common          
        
        path, ext = os.path.splitext(self.com.filename) 
        ext = ext[1:].lower()   
        suffix = "." + ext
        path, fn = os.path.split(self.com.filename)
        filename = fn[:-len(suffix)]
        
        self.path = self.com.path
        self.filename = filename
                          
        
        vboxtop = QtGui.QVBoxLayout()
        vboxtop.setContentsMargins(20,20,20,20)
        
        gridtop = QtGui.QGridLayout()
        gridtop.setVerticalSpacing(20)
    
        fontb = QtGui.QFont()
        fontb.setBold(True)            
        
        st1 = QtGui.QLabel(self)
        st1.setText('Save')
        st1.setFont(fontb)
        st2 = QtGui.QLabel(self)
        st2.setText('.pdf')
        st2.setFont(fontb)
        st3 = QtGui.QLabel(self)
        st3.setText('.png')
        st3.setFont(fontb)
        st4 = QtGui.QLabel(self)
        st4.setText('.svg')
        st4.setFont(fontb)
        st5 = QtGui.QLabel(self)
        st5.setText('.csv')
        st5.setFont(fontb)        
        st10 = QtGui.QLabel(self)
        st10.setText('.tif (data)')
        st10.setFont(fontb)                
        
        st6 = QtGui.QLabel(self)
        st6.setText('_spectrum')
        
        self.cb11 = QtGui.QCheckBox('', self)
        self.cb11.setChecked(True)
        self.cb12 = QtGui.QCheckBox('', self)
        self.cb13 = QtGui.QCheckBox('', self)   
        self.cb14 = QtGui.QCheckBox('', self)
        
        st7 = QtGui.QLabel(self)
        st7.setText('_composite_images')
        
        self.cb21 = QtGui.QCheckBox('', self)
        self.cb21.setChecked(True)
        self.cb22 = QtGui.QCheckBox('', self)
        self.cb23 = QtGui.QCheckBox('', self)
        self.cb24 = QtGui.QCheckBox('', self)
        
        st8 = QtGui.QLabel(self)
        st8.setText('_individual_images')  
        
        self.cb31 = QtGui.QCheckBox('', self)
        self.cb31.setChecked(True) 
        self.cb32 = QtGui.QCheckBox('', self)
        self.cb33 = QtGui.QCheckBox('', self)
        self.cb34 = QtGui.QCheckBox('', self)

        st9 = QtGui.QLabel(self)
        st9.setText('_scatter_plots')  
        
        self.cb41 = QtGui.QCheckBox('', self)
        self.cb41.setChecked(True) 
        self.cb42 = QtGui.QCheckBox('', self)
        self.cb43 = QtGui.QCheckBox('', self)
        
        
        gridtop.addWidget(st1, 0, 0)
        gridtop.addWidget(st2, 0, 1)
        gridtop.addWidget(st3, 0, 2)
        gridtop.addWidget(st4, 0, 3)
        gridtop.addWidget(st5, 0, 4)
        gridtop.addWidget(st10, 0, 5)
                                
        gridtop.addWidget(st6, 1, 0)
        gridtop.addWidget(self.cb11, 1, 1)
        gridtop.addWidget(self.cb12, 1, 2)  
        gridtop.addWidget(self.cb13, 1, 3)  
        gridtop.addWidget(self.cb14, 1, 4)           
  
        gridtop.addWidget(st7, 2, 0)
        gridtop.addWidget(self.cb21, 2, 1)
        gridtop.addWidget(self.cb22, 2, 2) 
        gridtop.addWidget(self.cb23, 2, 3)  
        gridtop.addWidget(self.cb24, 2, 5) 
        
        gridtop.addWidget(st8, 3, 0)
        gridtop.addWidget(self.cb31, 3, 1)
        gridtop.addWidget(self.cb32, 3, 2)  
        gridtop.addWidget(self.cb33, 3, 3)
        gridtop.addWidget(self.cb34, 3, 5)
        
        gridtop.addWidget(st9, 4, 0)
        gridtop.addWidget(self.cb41, 4, 1)
        gridtop.addWidget(self.cb42, 4, 2)  
        gridtop.addWidget(self.cb43, 4, 3)
        
                
        vboxtop.addStretch(0.5)
        vboxtop.addLayout(gridtop)
        vboxtop.addStretch(1)
        
         
        hbox0 = QtGui.QHBoxLayout()
         
        stf = QtGui.QLabel(self)
        stf.setText('Filename:\t')
        self.tc_savefn = QtGui.QLineEdit(self)
        self.tc_savefn.setText(self.filename)

        hbox0.addWidget(stf)
        hbox0.addWidget(self.tc_savefn)         
                  
        hbox1 = QtGui.QHBoxLayout()
                 
        stp = QtGui.QLabel(self)
        stp.setText('Path:  \t')
        self.tc_savepath = QtGui.QLineEdit(self)
        self.tc_savepath.setReadOnly(True)
        self.tc_savepath.setText(self.path)
        self.tc_savepath.setMinimumWidth(100)
        hbox1.addWidget(stp)
        hbox1.addWidget(self.tc_savepath)  
         
        button_path = QtGui.QPushButton('Browse...')
        button_path.clicked.connect(self.OnBrowseDir)
        hbox1.addWidget(button_path)
         
         
        hbox2 = QtGui.QHBoxLayout()
        button_save = QtGui.QPushButton('Save')
        button_save.clicked.connect(self.OnSave)
        hbox2.addWidget(button_save)
         
        button_cancel = QtGui.QPushButton('Cancel')
        button_cancel.clicked.connect(self.close)
        hbox2.addWidget(button_cancel)
        
        vboxtop.addLayout(hbox0)
        vboxtop.addLayout(hbox1)
        vboxtop.addStretch(1.0)
        vboxtop.addLayout(hbox2)

        
        
        self.setLayout(vboxtop)
        
#----------------------------------------------------------------------        
    def OnBrowseDir(self, evt):
        
        directory = QtGui.QFileDialog.getExistingDirectory(self, "Choose a directory", self.path, QtGui.QFileDialog.ShowDirsOnly|QtGui.QFileDialog.ReadOnly)       
                                                        
        
       
        if directory == '':
            return
                 
        directory = str(directory)
        self.com.path = directory
                    
        self.path = directory
        
        self.tc_savepath.setText(self.path)
            
            
                
#----------------------------------------------------------------------        
    def OnSave(self, evt):
        
        self.filename = str(self.tc_savefn.text())
        
        sp_pdf = self.cb11.isChecked()
        sp_png = self.cb12.isChecked()
        sp_svg = self.cb13.isChecked()
        sp_csv = self.cb14.isChecked()
        im_pdf = self.cb21.isChecked()
        im_png = self.cb22.isChecked()
        im_svg = self.cb23.isChecked()
        im_tif = self.cb24.isChecked()
        indim_pdf = self.cb31.isChecked()
        indim_png = self.cb32.isChecked()
        indim_svg = self.cb33.isChecked()
        indim_tif = self.cb34.isChecked()
        scatt_pdf = self.cb41.isChecked()
        scatt_png = self.cb42.isChecked()
        scatt_svg = self.cb43.isChecked()
        
        self.close() 
        self.parent.page3.Save(self.filename, self.path,
                                         spec_png = sp_png, 
                                         spec_pdf = sp_pdf, 
                                         spec_svg = sp_svg,
                                         spec_csv = sp_csv,
                                         img_png = im_png, 
                                         img_pdf = im_pdf,
                                         img_svg = im_svg,
                                         img_tif = im_tif,
                                         indimgs_png = indim_png, 
                                         indimgs_pdf = indim_pdf,
                                         indimgs_svg = indim_svg,
                                         indimgs_tif = indim_tif,
                                         scatt_png = scatt_png,
                                         scatt_pdf = scatt_pdf,
                                         scatt_svg = scatt_svg)   
        
                 
    
""" ------------------------------------------------------------------------------------------------"""
class PagePCA(QtGui.QWidget):
    def __init__(self, common, data_struct, stack, anlz):
        super(PagePCA, self).__init__()

        self.initUI(common, data_struct, stack, anlz)
        
#----------------------------------------------------------------------          
    def initUI(self, common, data_struct, stack, anlz): 

        self.com = common 
        self.data_struct = data_struct
        self.stk = stack       
        self.anlz = anlz
        
        
        self.selpca = 1       
        self.numsigpca = 2
        self.itheta = 0
        
        #panel 1        
        vbox1 = QtGui.QVBoxLayout()
        
        self.tc_PCAcomp = QtGui.QLabel(self)
        self.tc_PCAcomp.setText("PCA component ")
        vbox1.addWidget(self.tc_PCAcomp)
        
        hbox11 = QtGui.QHBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
        self.pcaimgfig = Figure((PlotH*1.10, PlotH))
        self.PCAImagePan = FigureCanvas(self.pcaimgfig)
        self.PCAImagePan.setParent(self) 
        
        fbox.addWidget(self.PCAImagePan)
        frame.setLayout(fbox)
        hbox11.addWidget(frame)                        
            
        self.slidershow = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slidershow.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slidershow.setRange(1,20)
        self.slidershow.setEnabled(False)          
        self.slidershow.valueChanged[int].connect(self.OnPCAScroll)

        
        hbox11.addWidget(self.slidershow)
        vbox1.addLayout(hbox11)        
        
        
        self.slider_theta = QtGui.QScrollBar(QtCore.Qt.Horizontal)
        self.slider_theta.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_theta.valueChanged[int].connect(self.OnScrollTheta)
        self.slider_theta.setRange(0, 100)  
        self.slider_theta.setMinimumWidth(250)   
        self.slider_theta.setVisible(False)  
        self.tc_imagetheta = QtGui.QLabel(self)
        self.tc_imagetheta.setText("4D Data Angle: ")
        self.tc_imagetheta.setVisible(False)
        hbox51 = QtGui.QHBoxLayout()
        hbox51.addWidget(self.tc_imagetheta) 
        hbox51.addStretch(1)
        hbox51.addWidget(self.slider_theta)
        hbox51.addStretch(1)
        vbox1.addLayout(hbox51)

   
                
        #panel 2
        vbox2 = QtGui.QVBoxLayout()
        vbox2.setContentsMargins(20,20,20,20)
        sizer2 = QtGui.QGroupBox('PCA')
        vbox21 = QtGui.QVBoxLayout()        
        
        self.button_calcpca = QtGui.QPushButton('Calculate PCA')
        self.button_calcpca.clicked.connect( self.OnCalcPCA)     
        self.button_calcpca.setEnabled(False)   
        vbox21.addWidget(self.button_calcpca)
        self.button_savepca = QtGui.QPushButton('Save PCA Results...')
        self.button_savepca.clicked.connect( self.OnSave)
        self.button_savepca.setEnabled(False)
        vbox21.addWidget(self.button_savepca)
        

        hbox21 = QtGui.QHBoxLayout()
        text1 = QtGui.QLabel(self)
        text1.setText('Number of significant components')
    
        
        self.npcaspin = QtGui.QSpinBox()
        self.npcaspin.setRange(1,20)      
        self.npcaspin.valueChanged[int].connect(self.OnNPCAspin)
        
        hbox21.addWidget(text1)
        hbox21.addWidget(self.npcaspin)
        vbox21.addLayout(hbox21)
              
        hbox22 = QtGui.QHBoxLayout()
        text2 = QtGui.QLabel(self)
        text2.setText( 'Cumulative variance')
        self.vartc = QtGui.QLabel(self)
        self.vartc.setText('0%')
        hbox22.addWidget(text2)
        hbox22.addWidget(self.vartc)

        vbox21.addLayout(hbox22)      

        self.button_movepcup = QtGui.QPushButton('Move PC up')
        self.button_movepcup.clicked.connect( self.OnMovePCUP)     
        self.button_movepcup.setEnabled(False)   
        vbox21.addWidget(self.button_movepcup)
        

#         line = QtGui.QFrame()
#         line.setFrameShape(QtGui.QFrame.HLine)
#         line.setFrameShadow(QtGui.QFrame.Sunken)         
#         vbox21.addWidget(line) 

        
        self.button_calcpca4D = QtGui.QPushButton('Calculate PCA for all angles')
        self.button_calcpca4D.clicked.connect( self.OnCalcPCA4D)     
        self.button_calcpca4D.setEnabled(False)   
        self.button_calcpca4D.setVisible(False) 
        vbox21.addWidget(self.button_calcpca4D)       
                 
        
        sizer2.setLayout(vbox21)
        vbox2.addStretch(1)
        vbox2.addWidget(sizer2)
        vbox2.addStretch(3)

        
        #panel 3
        vbox3 = QtGui.QVBoxLayout()
   

 
        self.text_pcaspec = QtGui.QLabel(self)
        self.text_pcaspec.setText("PCA spectrum ") 
        vbox3.addWidget(self.text_pcaspec)
                        
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()

        self.pcaspecfig = Figure((PlotW, PlotH))
        self.PCASpecPan = FigureCanvas(self.pcaspecfig)
        self.PCASpecPan.setParent(self)
       
        fbox.addWidget(self.PCASpecPan)
        frame.setLayout(fbox)
        vbox3.addWidget(frame)       
    
        
        #panel 4
        vbox4 = QtGui.QVBoxLayout()  
         
        text4 = QtGui.QLabel(self)
        text4.setText("PCA eigenvalues ")        
        vbox4.addWidget(text4)

        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
       
        self.pcaevalsfig = Figure((PlotW, PlotH*0.75))
        self.PCAEvalsPan = FigureCanvas(self.pcaevalsfig)
        self.PCAEvalsPan.setParent(self)
        self.PCAEvalsPan.mpl_connect('button_press_event', self.OnPointEvalsImage)
         
        fbox.addWidget(self.PCAEvalsPan)
        frame.setLayout(fbox)
        vbox4.addWidget(frame)                 

    
        
        vboxtop = QtGui.QVBoxLayout()
        gridsizertop = QtGui.QGridLayout()
        
        gridsizertop.addLayout(vbox2, 0, 0, QtCore .Qt. AlignLeft)
        gridsizertop.addLayout(vbox4, 0, 1)
        gridsizertop.addLayout(vbox1, 1, 0, QtCore .Qt. AlignLeft)
        gridsizertop.addLayout(vbox3, 1, 1)
        
        vboxtop.addStretch(1)
        vboxtop.addLayout(gridsizertop)
        vboxtop.addStretch(1)
        self.setLayout(vboxtop)

        
#----------------------------------------------------------------------
    def OnCalcPCA(self, event):
       
        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        self.calcpca = False  
        self.selpca = 1       
        self.numsigpca = 2
        self.slidershow.setValue(self.selpca)
        
        scrollmax = np.min([self.stk.n_ev, 20])
        self.slidershow.setMaximum(scrollmax)

        try: 
            self.CalcPCA()
            self.calcpca = True
            self.loadPCAImage()
            self.loadPCASpectrum()
            self.showEvals()
            self.com.pca_calculated = 1
            QtGui.QApplication.restoreOverrideCursor()
        except:
            self.com.pca_calculated = 0
            QtGui.QApplication.restoreOverrideCursor()
            QtGui.QMessageBox.warning(self, 'Error', 'PCA not calculated.')
        
        self.window().refresh_widgets()
        
        
#----------------------------------------------------------------------
    def OnCalcPCA4D(self, event):
       
        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        self.calcpca = False  
        self.selpca = 1       
        self.numsigpca = 2
        self.slidershow.setValue(self.selpca)
        
        scrollmax = np.min([self.stk.n_ev, 20])
        self.slidershow.setMaximum(scrollmax)



        try: 
            self.CalcPCA4D()
            self.calcpca = True
            self.loadPCAImage()
            self.loadPCASpectrum()
            self.showEvals()
            self.com.pca_calculated = 1
            self.com.pca4D_calculated = 1
            
            self.slider_theta.setVisible(True)  
            self.tc_imagetheta.setVisible(True)
            self.slider_theta.setRange(0, self.stk.n_theta-1)
            self.slider_theta.setValue(self.itheta)
            self.tc_imagetheta.setText("4D Data Angle: "+str(self.stk.theta[self.itheta]))  
            
            self.anlz.pcaimages = self.anlz.pcaimages4D[self.itheta]
            self.anlz.eigenvals = self.anlz.eigenvals4D[self.itheta]
            self.anlz.eigenvecs = self.anlz.eigenvecs4D[self.itheta]
            self.anlz.variance = self.anlz.variance4D[self.itheta]  
            self.anlz.pcaimagebounds = self.anlz.pcaimagebounds4D[self.itheta]        

            QtGui.QApplication.restoreOverrideCursor()
        except:
            self.com.pca_calculated = 0
            QtGui.QApplication.restoreOverrideCursor()
            QtGui.QMessageBox.warning(self, 'Error', 'PCA not calculated.')
            
            

        
        self.window().refresh_widgets()

#----------------------------------------------------------------------        
    def OnNPCAspin(self, value):
        num = value
        self.numsigpca = num
        
        if self.com.pca_calculated == 1:      
            self.anlz.numsigpca = self.numsigpca 
        
            # cumulative variance
            var = self.anlz.variance[:self.numsigpca].sum()
            self.vartc.setText(str(var.round(decimals=2)*100)+'%')
                 

#----------------------------------------------------------------------        
    def OnMovePCUP(self):       
        
        thiscomponent = self.selpca-1
        
        if thiscomponent == 0:
            return
        
        self.anlz.move_pc_up(thiscomponent)
 
        self.selpca = self.selpca-1
        self.loadPCAImage()
        self.loadPCASpectrum()
        self.showEvals()
        
        self.slidershow.setValue(self.selpca)
                   
        
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
               

 
#----------------------------------------------------------------------        
    def OnPCAScroll(self, value):
        self.sel = value
        self.selpca = self.sel
        if self.calcpca == True:
            self.loadPCAImage()
            self.loadPCASpectrum()


#----------------------------------------------------------------------            
    def OnScrollTheta(self, value):
        
        if self.com.pca4D_calculated == 0:
            return
        
        self.itheta = value
        
        self.anlz.pcaimages = self.anlz.pcaimages4D[self.itheta]
        self.anlz.eigenvals = self.anlz.eigenvals4D[self.itheta]
        self.anlz.eigenvecs = self.anlz.eigenvecs4D[self.itheta]
        self.anlz.variance = self.anlz.variance4D[self.itheta]
        self.anlz.pcaimagebounds = self.anlz.pcaimagebounds4D[self.itheta]
        
        self.tc_imagetheta.setText("4D Data Angle: "+'{0:5.2f}\t'.format(self.stk.theta[self.itheta]))

        var = self.anlz.variance[:self.numsigpca].sum()
        self.vartc.setText(str(var.round(decimals=2)*100)+'%')         
        
        self.loadPCAImage()
        self.loadPCASpectrum()
        self.showEvals()
        
        
        self.window().page0.itheta = self.itheta
        self.window().page0.slider_theta.setValue(self.itheta)
        
        self.window().page1.itheta = self.itheta
        self.window().page1.slider_theta.setValue(self.itheta)
            
#----------------------------------------------------------------------  
    def OnPointEvalsImage(self, evt):
        x = evt.xdata
        y = evt.ydata
                
        if self.com.pca_calculated == 1:     
            #Find the closest point to the point clicked on the plot
            self.selpca = int(np.round(x))
            
            if self.selpca < 1:
                self.selpca = 1
                       
            self.loadPCAImage()
            self.loadPCASpectrum()
            


#----------------------------------------------------------------------    
    def OnSave(self, event):     

        savewin = SaveWinP2(self.window())
        savewin.show()


            
#----------------------------------------------------------------------    
    def Save(self, filename, path, spec_png = True, spec_pdf = False, spec_svg = False, spec_csv = False,
             img_png = True, img_pdf = False, img_svg = False, img_tif = False,
             evals_png = True, evals_pdf = False, evals_svg = False): 
        
        
        self.SaveFileName = os.path.join(path,filename)
   
        try: 

            matplotlib.rcParams['pdf.fonttype'] = 42
            if evals_png:
                ext = 'png'
                suffix = "." + ext
                fileName_evals = self.SaveFileName+"_PCAevals."+ext
                            
                fig = self.pcaevalsfig
                fig.savefig(fileName_evals)
                
            if evals_pdf:
                ext = 'pdf'
                suffix = "." + ext
                fileName_evals = self.SaveFileName+"_PCAevals."+ext
                            
                fig = self.pcaevalsfig
                fig.savefig(fileName_evals)     
                
            if evals_svg:
                ext = 'svg'
                suffix = "." + ext
                fileName_evals = self.SaveFileName+"_PCAevals."+ext
                            
                fig = self.pcaevalsfig
                fig.savefig(fileName_evals)             

                        
            from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
            matplotlib.rcParams['pdf.fonttype'] = 42
            
            ext = 'png'
            suffix = "." + ext        
                
            if img_png:
                for i in range (self.numsigpca):
              
                    self.pcaimage = self.anlz.pcaimages[:,:,i]
              
                    fig = matplotlib.figure.Figure(figsize =(PlotH*1.15, PlotH))
                    canvas = FigureCanvas(fig)
                    axes = fig.gca()
                    divider = make_axes_locatable(axes)
                    ax_cb = divider.new_horizontal(size="3%", pad=0.03)  
                    fig.add_axes(ax_cb)
                    axes.set_position([0.03,0.03,0.8,0.94])
                    bound = self.anlz.pcaimagebounds[i]       
        
                    im = axes.imshow(np.rot90(self.pcaimage), cmap=matplotlib.cm.get_cmap("seismic_r"), vmin = -bound, vmax = bound)
                    cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)  
                    axes.axis("off") 
                                
                    fileName_img = self.SaveFileName+"_PCA_" +str(i+1)+"."+ext
                    fig.savefig(fileName_img, bbox_inches='tight', pad_inches = 0.0)
            
            if spec_png:
                for i in range (self.numsigpca):
                
                    pcaspectrum = self.anlz.eigenvecs[:,i]
                    fig = matplotlib.figure.Figure(figsize =(PlotW, PlotH))
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.15,0.15,0.75,0.75))
                    axes = fig.gca()
                    specplot = axes.plot(self.stk.ev, pcaspectrum)    
                    axes.set_xlabel('Photon Energy [eV]')
                    axes.set_ylabel('Optical Density')
                
                    fileName_spec = self.SaveFileName+"_PCAspectrum_" +str(i+1)+"."+ext
                    fig.savefig(fileName_spec)
                    
            if spec_csv:
                for i in range (self.numsigpca):
                    pcaspectrum = self.anlz.eigenvecs[:,i]
                    fileName_spec = self.SaveFileName+"_PCAspectrum_" +str(i+1)+".csv"
                    cname = "PCAspectrum_" +str(i+1)
                    self.stk.write_csv(fileName_spec, self.stk.ev, pcaspectrum, cname = cname)
                    
                
            ext = 'pdf'
            suffix = "." + ext        
                
            if img_pdf:
                for i in range (self.numsigpca):
              
                    self.pcaimage = self.anlz.pcaimages[:,:,i]
              
                    fig = matplotlib.figure.Figure(figsize =(PlotH*1.15, PlotH))
                    canvas = FigureCanvas(fig)
                    axes = fig.gca()
                    divider = make_axes_locatable(axes)
                    ax_cb = divider.new_horizontal(size="3%", pad=0.03)  
                    fig.add_axes(ax_cb)
                    axes.set_position([0.03,0.03,0.8,0.94])
                    bound = self.anlz.pcaimagebounds[i]            
        
                    im = axes.imshow(np.rot90(self.pcaimage), cmap=matplotlib.cm.get_cmap("seismic_r"), vmin = -bound, vmax = bound)
                    cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)  
                    axes.axis("off") 
                                
                    fileName_img = self.SaveFileName+"_PCA_" +str(i+1)+"."+ext
                    fig.savefig(fileName_img, bbox_inches='tight', pad_inches = 0.0)
            
            if spec_pdf:
                for i in range (self.numsigpca):
                
                    self.pcaspectrum = self.anlz.eigenvecs[:,i]
                    fig = matplotlib.figure.Figure(figsize =(PlotW, PlotH))
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.15,0.15,0.75,0.75))
                    axes = fig.gca()
                    specplot = axes.plot(self.stk.ev,self.pcaspectrum)    
                    axes.set_xlabel('Photon Energy [eV]')
                    axes.set_ylabel('Optical Density')
                
                    fileName_spec = self.SaveFileName+"_PCAspectrum_" +str(i+1)+"."+ext
                    fig.savefig(fileName_spec)                

            ext = 'svg'
            suffix = "." + ext        
                
            if img_svg:
                for i in range (self.numsigpca):
              
                    self.pcaimage = self.anlz.pcaimages[:,:,i]
              
                    fig = matplotlib.figure.Figure(figsize =(PlotH*1.15, PlotH))
                    canvas = FigureCanvas(fig)
                    axes = fig.gca()
                    divider = make_axes_locatable(axes)
                    ax_cb = divider.new_horizontal(size="3%", pad=0.03)  
                    fig.add_axes(ax_cb)
                    axes.set_position([0.03,0.03,0.8,0.94])
                    bound = self.anlz.pcaimagebounds[i]            
        
                    im = axes.imshow(np.rot90(self.pcaimage), cmap=matplotlib.cm.get_cmap("seismic_r"), vmin = -bound, vmax = bound)
                    cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)  
                    axes.axis("off") 
                                
                    fileName_img = self.SaveFileName+"_PCA_" +str(i+1)+"."+ext
                    fig.savefig(fileName_img, bbox_inches='tight', pad_inches = 0.0)
            
            if spec_svg:
                for i in range (self.numsigpca):
                
                    self.pcaspectrum = self.anlz.eigenvecs[:,i]
                    fig = matplotlib.figure.Figure(figsize =(PlotW, PlotH))
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.15,0.15,0.75,0.75))
                    axes = fig.gca()
                    specplot = axes.plot(self.stk.ev,self.pcaspectrum)    
                    axes.set_xlabel('Photon Energy [eV]')
                    axes.set_ylabel('Optical Density')
                
                    fileName_spec = self.SaveFileName+"_PCAspectrum_" +str(i+1)+"."+ext
                    fig.savefig(fileName_spec)  
                    
            if img_tif:
                for i in range (self.numsigpca):
              
                    self.pcaimage = self.anlz.pcaimages[:,:,i]
                              
                    fileName_img = self.SaveFileName+"_PCA_" +str(i+1)+".tif"

                    img1 = Image.fromarray(self.pcaimage)
                    img1.save(fileName_img)
                                
        except IOError, e:
            if e.strerror:
                err = e.strerror 
            else: 
                err = e 
   
            QtGui.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err)

            
        
     
#----------------------------------------------------------------------      
    def showEvals(self):
        
        evalmax = np.min([self.stk.n_ev, 40])
        self.pcaevals = self.anlz.eigenvals[0:evalmax]
        

        fig = self.pcaevalsfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        
       
        evalsplot = axes.semilogy(np.arange(1,evalmax+1), self.pcaevals,'b.')    
        
        axes.set_xlabel('Principal Component')
        axes.set_ylabel('Log(Eigenvalue)')
         

        self.PCAEvalsPan.draw()


#----------------------------------------------------------------------      
    def loadPCAImage(self):
        
        self.tc_PCAcomp.setText("PCA component " + str(self.selpca))
        self.text_pcaspec.setText("PCA spectrum "+ str(self.selpca))          
        
        self.pcaimage = self.anlz.pcaimages[:,:,self.selpca-1]
        
        
        fig = self.pcaimgfig
        fig.clf()
     
        axes = fig.gca()
    
        divider = make_axes_locatable(axes)
        ax_cb = divider.new_horizontal(size="3%", pad=0.03)  

        fig.add_axes(ax_cb)
        
        axes.set_position([0.03,0.03,0.8,0.94])
        
        bound = self.anlz.pcaimagebounds[self.selpca-1]
        
     
        im = axes.imshow(np.rot90(self.pcaimage), cmap=matplotlib.cm.get_cmap("seismic_r"), vmin = -bound, vmax = bound)
        cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)  
    
        axes.axis("off") 
        self.PCAImagePan.draw()

        
        
#----------------------------------------------------------------------     
    def loadPCASpectrum(self):

        self.pcaspectrum = self.anlz.eigenvecs[:,self.selpca-1]
            
        
        fig = self.pcaspecfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        

        specplot = axes.plot(self.stk.ev,self.pcaspectrum)
                        
        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Optical Density')
        
        self.PCASpecPan.draw()
        

#---------------------------------------------------------------------- 
class SaveWinP2(QtGui.QDialog):

    def __init__(self, parent):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        
        self.resize(400, 300)
        self.setWindowTitle('Save')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)


        self.com = self.parent.common          
        
        path, ext = os.path.splitext(self.com.filename) 
        ext = ext[1:].lower()   
        suffix = "." + ext
        path, fn = os.path.split(self.com.filename)
        filename = fn[:-len(suffix)]
        
        self.path = self.com.path
        self.filename = filename
                          
        
        vboxtop = QtGui.QVBoxLayout()
        vboxtop.setContentsMargins(20,20,20,20)
        
        gridtop = QtGui.QGridLayout()
        gridtop.setVerticalSpacing(20)
    
        fontb = QtGui.QFont()
        fontb.setBold(True)            
        
        st1 = QtGui.QLabel(self)
        st1.setText('Save')
        st1.setFont(fontb)
        st2 = QtGui.QLabel(self)
        st2.setText('.pdf')
        st2.setFont(fontb)
        st3 = QtGui.QLabel(self)
        st3.setText('.png')
        st3.setFont(fontb)
        st4 = QtGui.QLabel(self)
        st4.setText('.svg')
        st4.setFont(fontb)
        st9 = QtGui.QLabel(self)
        st9.setText('.tif (data)')
        st9.setFont(fontb)    
        st5 = QtGui.QLabel(self)
        st5.setText('.csv')
        st5.setFont(fontb)        
        
        
        st6 = QtGui.QLabel(self)
        st6.setText('_spectrum')
        
        self.cb11 = QtGui.QCheckBox('', self)
        self.cb11.setChecked(True)
        self.cb12 = QtGui.QCheckBox('', self)
        self.cb13 = QtGui.QCheckBox('', self)   
        self.cb14 = QtGui.QCheckBox('', self)
        
        st7 = QtGui.QLabel(self)
        st7.setText('_image')
        
        self.cb21 = QtGui.QCheckBox('', self)
        self.cb21.setChecked(True)
        self.cb22 = QtGui.QCheckBox('', self)
        self.cb23 = QtGui.QCheckBox('', self)
        self.cb24 = QtGui.QCheckBox('', self)
        
        st8 = QtGui.QLabel(self)
        st8.setText('_eigenvals')   
        self.cb31 = QtGui.QCheckBox('', self)
        self.cb32 = QtGui.QCheckBox('', self)
        self.cb33 = QtGui.QCheckBox('', self)


        gridtop.addWidget(st1, 0, 0)
        gridtop.addWidget(st2, 0, 1)
        gridtop.addWidget(st3, 0, 2)
        gridtop.addWidget(st4, 0, 3)
        gridtop.addWidget(st9, 0, 5)
        gridtop.addWidget(st5, 0, 4)
                                
        gridtop.addWidget(st6, 1, 0)
        gridtop.addWidget(self.cb11, 1, 1)
        gridtop.addWidget(self.cb12, 1, 2)  
        gridtop.addWidget(self.cb13, 1, 3)  
        gridtop.addWidget(self.cb14, 1, 4)           
  
        gridtop.addWidget(st7, 2, 0)
        gridtop.addWidget(self.cb21, 2, 1)
        gridtop.addWidget(self.cb22, 2, 2) 
        gridtop.addWidget(self.cb23, 2, 3)  
        gridtop.addWidget(self.cb24, 2, 5)  
        
        gridtop.addWidget(st8, 3, 0)
        gridtop.addWidget(self.cb31, 3, 1)
        gridtop.addWidget(self.cb32, 3, 2) 
        gridtop.addWidget(self.cb33, 3, 3) 
        
        vboxtop.addStretch(0.5)
        vboxtop.addLayout(gridtop)
        vboxtop.addStretch(1)
        
         
        hbox0 = QtGui.QHBoxLayout()
         
        stf = QtGui.QLabel(self)
        stf.setText('Filename:\t')
        self.tc_savefn = QtGui.QLineEdit(self)
        self.tc_savefn.setText(self.filename)

        hbox0.addWidget(stf)
        hbox0.addWidget(self.tc_savefn)         
                  
        hbox1 = QtGui.QHBoxLayout()
                 
        stp = QtGui.QLabel(self)
        stp.setText('Path:  \t')
        self.tc_savepath = QtGui.QLineEdit(self)
        self.tc_savepath.setReadOnly(True)
        self.tc_savepath.setText(self.path)
        self.tc_savepath.setMinimumWidth(100)
        hbox1.addWidget(stp)
        hbox1.addWidget(self.tc_savepath)  
         
        button_path = QtGui.QPushButton('Browse...')
        button_path.clicked.connect(self.OnBrowseDir)
        hbox1.addWidget(button_path)
         
         
        hbox2 = QtGui.QHBoxLayout()
        button_save = QtGui.QPushButton('Save')
        button_save.clicked.connect(self.OnSave)
        hbox2.addWidget(button_save)
         
        button_cancel = QtGui.QPushButton('Cancel')
        button_cancel.clicked.connect(self.close)
        hbox2.addWidget(button_cancel)
        
        vboxtop.addLayout(hbox0)
        vboxtop.addLayout(hbox1)
        vboxtop.addStretch(1.0)
        vboxtop.addLayout(hbox2)

        
        
        self.setLayout(vboxtop)
        
#----------------------------------------------------------------------        
    def OnBrowseDir(self, evt):
        
        directory = QtGui.QFileDialog.getExistingDirectory(self, "Choose a directory", self.path, QtGui.QFileDialog.ShowDirsOnly|QtGui.QFileDialog.ReadOnly)       
                                                        
        
       
        if directory == '':
            return
                 
        directory = str(directory)
        self.com.path = directory
                    
        self.path = directory
        
        self.tc_savepath.setText(self.path)
            
            
                
#----------------------------------------------------------------------        
    def OnSave(self, evt):
        
        self.filename = str(self.tc_savefn.text())
        
        sp_pdf = self.cb11.isChecked()
        sp_png = self.cb12.isChecked()
        sp_svg = self.cb13.isChecked()
        sp_csv = self.cb14.isChecked()
        im_pdf = self.cb21.isChecked()
        im_png = self.cb22.isChecked()
        im_svg = self.cb23.isChecked()
        im_tif = self.cb24.isChecked()
        ev_pdf = self.cb31.isChecked()
        ev_png = self.cb32.isChecked()
        ev_svg = self.cb33.isChecked()
        
        self.close() 
        self.parent.page2.Save(self.filename, self.path,
                               spec_png = sp_png, 
                               spec_pdf = sp_pdf, 
                               spec_svg = sp_svg,
                               spec_csv = sp_csv,
                               img_png = im_png, 
                               img_pdf = im_pdf,
                               img_svg = im_svg,
                               img_tif = im_tif,
                               evals_png = ev_png, 
                               evals_pdf = ev_pdf,
                               evals_svg = ev_svg)        
        

""" ------------------------------------------------------------------------------------------------"""
class PageStack(QtGui.QWidget):
    def __init__(self, common, data_struct, stack):
        super(PageStack, self).__init__()

        self.initUI(common, data_struct, stack)
        
#----------------------------------------------------------------------          
    def initUI(self, common, data_struct, stack): 
        
        self.data_struct = data_struct
        self.stk = stack
        self.com = common                  
        
        self.filename = " "
       
        self.ix = 0
        self.iy = 0
        self.iev = 50  
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
        

        #panel 1
        sizer1 = QtGui.QGroupBox('Preprocess')
        vbox1 = QtGui.QVBoxLayout()
        vbox1.setSpacing(0)
        
        self.button_align = QtGui.QPushButton('Align stack...')
        self.button_align.clicked.connect(self.OnAlignImgs)
        self.button_align.setEnabled(False)
        vbox1.addWidget(self.button_align) 
        
  
        
        self.button_limitev = QtGui.QPushButton('Limit energy range...')
        self.button_limitev.clicked.connect( self.OnLimitEv)
        self.button_limitev.setEnabled(False)
        vbox1.addWidget(self.button_limitev)
        
        self.button_subregion = QtGui.QPushButton('Clip to subregion...')
        self.button_subregion.clicked.connect(self.OnCliptoSubregion)
        self.button_subregion.setEnabled(False)
        vbox1.addWidget(self.button_subregion)    
        
        self.button_darksig = QtGui.QPushButton('Dark signal subtraction...')
        self.button_darksig.clicked.connect(self.OnDarkSignal)
        self.button_darksig.setEnabled(False)
        vbox1.addWidget(self.button_darksig)     
        
            
        self.button_savestack = QtGui.QPushButton('Save processed stack')
        self.button_savestack.clicked.connect(self.OnSaveStack)
        self.button_savestack.setEnabled(False)          
        vbox1.addWidget(self.button_savestack)
        
        sizer1.setLayout(vbox1)
        

        #panel 1B
        sizer1b = QtGui.QGroupBox('Normalize')
        vbox1b = QtGui.QVBoxLayout()
        vbox1b.setSpacing(0)
        
        self.button_i0ffile = QtGui.QPushButton('I0 from file...')
        self.button_i0ffile.clicked.connect(self.OnI0FFile)
        self.button_i0ffile.setEnabled(False)
        vbox1b.addWidget(self.button_i0ffile)
        self.button_i0histogram = QtGui.QPushButton('I0 from histogram...')
        self.button_i0histogram.clicked.connect( self.OnI0histogram)   
        self.button_i0histogram.setEnabled(False)     
        vbox1b.addWidget(self.button_i0histogram)
        self.button_showi0 = QtGui.QPushButton('Show I0...')
        self.button_showi0.clicked.connect( self.OnShowI0)   
        self.button_showi0.setEnabled(False)
        vbox1b.addWidget(self.button_showi0)
        
        self.button_prenorm = QtGui.QPushButton('Use pre-normalized data')
        self.button_prenorm.clicked.connect(self.OnPreNormalizedData)   
        self.button_prenorm.setEnabled(False)
        vbox1b.addWidget(self.button_prenorm)
        
        self.button_refimgs = QtGui.QPushButton('Load Reference Images')
        self.button_refimgs.clicked.connect(self.OnRefImgs)
        self.button_refimgs.setEnabled(False)
        vbox1b.addWidget(self.button_refimgs)    
        
        self.button_reseti0 = QtGui.QPushButton('Reset I0')
        self.button_reseti0.clicked.connect(self.OnResetI0)   
        self.button_reseti0.setEnabled(False)
        vbox1b.addWidget(self.button_reseti0)
        
        self.button_saveod = QtGui.QPushButton('Save OD data')
        self.button_saveod.clicked.connect(self.OnSaveOD)   
        self.button_saveod.setEnabled(False)
        vbox1b.addWidget(self.button_saveod)        
        
        sizer1b.setLayout(vbox1b)  


        #panel 2
        sizer2 = QtGui.QGroupBox('Display')
        vbox2 = QtGui.QVBoxLayout()
        vbox2.setSpacing(0)
        
        hbox20 = QtGui.QHBoxLayout()
        
        sizer21 = QtGui.QGroupBox('File')
        vbox21 = QtGui.QVBoxLayout()

        self.textctrl = QtGui.QLabel(self)
        vbox21.addWidget(self.textctrl)
        self.textctrl.setText('File name')
        
        sizer21.setLayout(vbox21)
        
        hbox20.addWidget(sizer21)
        hbox20.addSpacing(15)
        
        vbox22 = QtGui.QVBoxLayout()
        self.button_slideshow = QtGui.QPushButton('Play stack movie')
        self.button_slideshow.setMaximumSize (150 , 150)
        self.button_slideshow.clicked.connect( self.OnSlideshow)
        self.button_slideshow.setEnabled(False)
        vbox22.addWidget(self.button_slideshow)
        
        self.button_save = QtGui.QPushButton( 'Save images...')
        self.button_save.setMaximumSize (150 , 150)
        self.button_save.clicked.connect( self.OnSave)
        self.button_save.setEnabled(False)          
        vbox22.addWidget(self.button_save)
              
        hbox20.addLayout(vbox22)
        vbox2.addLayout(hbox20)
        vbox2.addSpacing(5)
        

        hbox21 = QtGui.QHBoxLayout()
        
        sizer22 = QtGui.QGroupBox('Image')
        vbox23 = QtGui.QVBoxLayout()

        self.rb_flux = QtGui.QRadioButton( 'Flux', self)
        self.rb_od = QtGui.QRadioButton('Optical Density',self)
        self.rb_flux.setChecked(True)
        self.rb_flux.toggled.connect(self.OnRb_fluxod)
        

        vbox23.addWidget(self.rb_flux)
        vbox23.addWidget(self.rb_od)
        vbox23.addStretch (1)
         
        self.rb_flux.setEnabled(False)
        self.rb_od.setEnabled(False)
 
        self.add_scale_cb = QtGui.QCheckBox('Scalebar', self) 
        self.add_scale_cb.setChecked(True)
        self.add_scale_cb.stateChanged.connect(self.OnShowScale)
        vbox23.addWidget(self.add_scale_cb)
        
        hbox211 = QtGui.QHBoxLayout()
        hbox211.addSpacing(20)
        self.cb_white_scale_bar = QtGui.QCheckBox('White', self) 
        self.cb_white_scale_bar.setChecked(False)
        self.cb_white_scale_bar.stateChanged.connect(self.OnWhiteScale)
        hbox211.addWidget(self.cb_white_scale_bar)
        vbox23.addLayout(hbox211)
         
        self.add_colbar_cb = QtGui.QCheckBox('Colorbar', self)
        self.add_colbar_cb.stateChanged.connect(self.OnShowColBar)
        vbox23.addWidget(self.add_colbar_cb)
        sizer22.setLayout(vbox23)
         
        hbox21.addWidget(sizer22)
        hbox21.addSpacing(5)
        vbox2.addLayout(hbox21)
                 
        sizer23 = QtGui.QGroupBox('Display settings')
        hbox23 = QtGui.QHBoxLayout() 
 
        fgs21 = QtGui.QGridLayout()
         
        self.tc_min = QtGui.QLabel(self)
        self.tc_min.setText('Minimum: \t{0:5d}%'.format(int(100*self.brightness_min)))
         
        self.tc_max = QtGui.QLabel(self)
        self.tc_max.setText('Maximum:{0:5d}%'.format(int(100*self.brightness_max)))
 
        self.slider_brightness_min = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        self.slider_brightness_min.setRange(0,49)
        self.slider_brightness_min.setValue(self.dispbrightness_min)
        self.slider_brightness_min.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_brightness_min.valueChanged[int].connect(self.OnScrollBrightnessMin)

                 
        self.slider_brightness_max = QtGui.QSlider(QtCore.Qt.Horizontal, self)      
        self.slider_brightness_max.setRange(50,120)
        self.slider_brightness_max.setValue(self.dispbrightness_max)
        self.slider_brightness_max.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_brightness_max.valueChanged[int].connect(self.OnScrollBrightnessMax)  
         
        self.tc_gamma = QtGui.QLabel(self)
        self.tc_gamma.setText('Gamma:  \t{0:5.2f}'.format(self.gamma))
         
        self.slider_gamma = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        self.slider_gamma.setRange(1,20)
        self.slider_gamma.setValue(self.displaygamma)
        self.slider_gamma.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_gamma.valueChanged[int].connect(self.OnScrollGamma)  
               
         
        fgs21.addWidget(self.tc_min, 0, 0)
        fgs21.addWidget(self.slider_brightness_min, 0, 1)
        fgs21.addWidget(self.tc_max, 1, 0) 
        fgs21.addWidget(self.slider_brightness_max, 1, 1)
        fgs21.addWidget(self.tc_gamma, 2, 0)
        fgs21.addWidget(self.slider_gamma, 2, 1)
         
        hbox23.addLayout(fgs21)
                
    
        vbox24 = QtGui.QVBoxLayout()
        self.button_despike = QtGui.QPushButton('Despike')
        self.button_despike.clicked.connect( self.OnDespike)   
        self.button_despike.setEnabled(False)     
        vbox24.addWidget(self.button_despike)
        self.button_resetdisplay = QtGui.QPushButton( 'Reset')
        self.button_resetdisplay.clicked.connect( self.OnResetDisplaySettings)   
        self.button_resetdisplay.setEnabled(False)     
        vbox24.addWidget(self.button_resetdisplay)
        self.button_displaycolor = QtGui.QPushButton('Color Table...   ')
        self.button_displaycolor.clicked.connect( self.OnSetColorTable)   
        self.button_displaycolor.setEnabled(False)     
        vbox24.addWidget(self.button_displaycolor)
        
        hbox23.addSpacing(20)
        hbox23.addLayout(vbox24)
        
        sizer23.setLayout(hbox23)
        hbox21.addWidget(sizer23)       
        sizer2.setLayout(vbox2)
        
        
        
        #panel 3
        sizer3 = QtGui.QGroupBox('Region of Interest')
        vbox3 = QtGui.QVBoxLayout()
        vbox3.setSpacing(0)

        
        self.button_addROI = QtGui.QPushButton('Add ROI')
        self.button_addROI.clicked.connect( self.OnAddROI)
        self.button_addROI.setEnabled(False)
        vbox3.addWidget(self.button_addROI)
        
        self.button_acceptROI = QtGui.QPushButton('Accept ROI')
        self.button_acceptROI.clicked.connect( self.OnAcceptROI)   
        self.button_acceptROI.setEnabled(False)     
        vbox3.addWidget(self.button_acceptROI)
        
        self.button_resetROI = QtGui.QPushButton('Reset ROI')
        self.button_resetROI.clicked.connect( self.OnResetROI)
        self.button_resetROI.setEnabled(False)
        vbox3.addWidget(self.button_resetROI) 

        self.button_setROII0 = QtGui.QPushButton('Set ROI As I0')
        self.button_setROII0.clicked.connect( self.OnSetROII0)
        self.button_setROII0.setEnabled(False)
        vbox3.addWidget(self.button_setROII0)
        
        self.button_saveROIspectr = QtGui.QPushButton( 'Save ROI Spectrum...')
        self.button_saveROIspectr.clicked.connect( self.OnSaveROISpectrum)   
        self.button_saveROIspectr.setEnabled(False)     
        vbox3.addWidget(self.button_saveROIspectr)
        
        self.button_ROIdosecalc = QtGui.QPushButton('ROI Dose Calculation...')
        self.button_ROIdosecalc.clicked.connect( self.OnROI_DoseCalc)   
        self.button_ROIdosecalc.setEnabled(False)     
        vbox3.addWidget(self.button_ROIdosecalc)        
        
        self.button_spectralROI = QtGui.QPushButton('Spectral ROI...')
        self.button_spectralROI.clicked.connect( self.OnSpectralROI)   
        self.button_spectralROI.setEnabled(False)     
        vbox3.addWidget(self.button_spectralROI)
        sizer3.setLayout(vbox3)
        

        

        #panel 4     
        #vbox4 = QtGui.QVBoxLayout()
        gridsizer4 = QtGui.QGridLayout() 
        
        self.tc_imageeng = QtGui.QLabel(self)
        self.tc_imageeng.setText("Image at energy: ")
        gridsizer4.addWidget(self.tc_imageeng, 0, 0, QtCore .Qt. AlignLeft)
        

        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
   
        self.absimgfig = Figure((PlotH, PlotH))
        self.AbsImagePanel = FigureCanvas(self.absimgfig)
        self.AbsImagePanel.setParent(self)
        self.AbsImagePanel.mpl_connect('button_press_event', self.OnPointAbsimage)
        fbox.addWidget(self.AbsImagePanel)
        frame.setLayout(fbox)
        gridsizer4.addWidget(frame, 1, 0, QtCore .Qt. AlignLeft)
       

        self.slider_eng = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slider_eng.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_eng.valueChanged[int].connect(self.OnScrollEng)
        self.slider_eng.setRange(0, 100)
        gridsizer4.addWidget(self.slider_eng, 1, 1, QtCore .Qt. AlignLeft)

        self.slider_theta = QtGui.QScrollBar(QtCore.Qt.Horizontal)
        self.slider_theta.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_theta.valueChanged[int].connect(self.OnScrollTheta)
        self.slider_theta.setRange(0, 100)     
        self.slider_theta.setVisible(False)  
        self.tc_imagetheta = QtGui.QLabel(self)
        self.tc_imagetheta.setText("4D Data Angle: ")
        self.tc_imagetheta.setVisible(False)
        hbox41 = QtGui.QHBoxLayout()
        hbox41.addWidget(self.tc_imagetheta) 
        hbox41.addWidget(self.slider_theta)
        gridsizer4.addLayout(hbox41, 2, 0)
        

        #panel 5             
        self.tc_spec = QtGui.QLabel(self)
        self.tc_spec.setText("Spectrum ")
        gridsizer4.addWidget(self.tc_spec, 0, 2, QtCore.Qt.AlignLeft)
        
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
      
        self.specfig = Figure((PlotW, PlotH))
        self.SpectrumPanel = FigureCanvas(self.specfig)
        self.SpectrumPanel.setParent(self)
        self.SpectrumPanel.mpl_connect('button_press_event', self.OnPointSpectrum)

        fbox.addWidget(self.SpectrumPanel)
        frame.setLayout(fbox)
        gridsizer4.addWidget(frame, 1, 2,  QtCore.Qt.AlignLeft)             
    
    
        vboxtop = QtGui.QVBoxLayout()
        
        hboxtop = QtGui.QHBoxLayout()
        hboxtop.addWidget(sizer1)
        hboxtop.addWidget(sizer1b)
        hboxtop.addWidget(sizer2)
        hboxtop.addWidget(sizer3)
        
        hboxbott = QtGui.QHBoxLayout()
        hboxbott2 = QtGui.QHBoxLayout()
        hboxbott2.addLayout(gridsizer4)
        hboxbott.addLayout(hboxbott2)
              
        vboxtop.addStretch (1)
        vboxtop.addLayout(hboxtop)
        vboxtop.addStretch (1)
        vboxtop.addLayout(hboxbott)
        vboxtop.addStretch (1)


        vboxtop.setContentsMargins(20,20,20,20)
        self.setLayout(vboxtop)


#----------------------------------------------------------------------
        
    def OnI0FFile(self, event):

            
        try:  
            wildcard = "I0 CSV files (*.csv);; I0 files (*.xas);;SDF I0 files (*.hdr)"
            
            filepath = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '', wildcard)
            

            filepath = str(filepath)
            if filepath == '':
                return
            
            filepath_i0 =  os.path.dirname(str(filepath))
            self.filename =  os.path.basename(str(filepath))
            

                                                        
            basename, extension = os.path.splitext(self.filename)      
            
            
            if extension == '.hdr':
                QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))                           

                x=self.stk.n_cols
                y=self.stk.n_rows
                z=self.iev               

                self.ix = x/2
                self.iy = y/2
                self.stk.read_sdf_i0(filepath)
                self.com.i0_loaded = 1
                self.loadSpectrum(self.ix, self.iy)
                self.loadImage()
                
                QtGui.QApplication.restoreOverrideCursor()
                
                
            elif extension == '.xas':
                QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))                          

                x=self.stk.n_cols
                y=self.stk.n_rows
                z=self.iev               

                self.ix = x/2
                self.iy = y/2

                self.stk.read_stk_i0(filepath, extension)
                
                self.com.i0_loaded = 1
                self.loadSpectrum(self.ix, self.iy)
                self.loadImage()
                
                QtGui.QApplication.restoreOverrideCursor()

            elif extension == '.csv':
                QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))                          

                x=self.stk.n_cols
                y=self.stk.n_rows
                z=self.iev               

                self.ix = x/2
                self.iy = y/2

                self.stk.read_stk_i0(filepath, extension)

                self.com.i0_loaded = 1
                self.loadSpectrum(self.ix, self.iy)
                self.loadImage()
                
                QtGui.QApplication.restoreOverrideCursor()
            
        except:
            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))  
            self.com.i0_loaded = 0        
            QtGui.QApplication.restoreOverrideCursor()
            QtGui.QMessageBox.warning(self,'Error',"I0 file not loaded.")
            import sys; print sys.exc_info()
                       
                          
        self.window().refresh_widgets()

        
        
        
#----------------------------------------------------------------------       
    def OnI0histogram(self, event):    
        #self.window().Hide()
        histogram = ShowHistogram(self, self.stk)
        histogram.show()
         

#----------------------------------------------------------------------       
    def I0histogramCalculated(self):    
        
        plot = PlotFrame(self, self.stk.evi0hist,self.stk.i0datahist)
        plot.show()
        
        self.com.i0_loaded = 1
        
        if self.com.stack_4d == 1:
            self.stk.od3d = self.stk.od4D[:,:,:,self.itheta]
            self.stk.od = self.stk.od3d.copy()
            self.stk.od = np.reshape(self.stk.od, (self.stk.n_cols*self.stk.n_rows, self.stk.n_ev), order='F')      
         
        self.loadSpectrum(self.ix, self.iy)
        self.loadImage()
        
        self.window().refresh_widgets()
        
#----------------------------------------------------------------------    
    def OnShowI0(self, event):     

        plot = PlotFrame(self, self.stk.evi0,self.stk.i0data)
        plot.show()
        
#----------------------------------------------------------------------    
    def OnPreNormalizedData(self, event):     

        self.stk.UsePreNormalizedData()
        
        self.com.i0_loaded = 1
        self.loadSpectrum(self.ix, self.iy)
        self.loadImage()
        
        self.window().refresh_widgets()
        
#----------------------------------------------------------------------    
    def OnRefImgs(self, event):     

        #Load .xrm reference images
        try:
        #if True:        
            wildcard = "Reference images (*.xrm)"
            
            filepaths = QtGui.QFileDialog.getOpenFileNames(self, 'Select reference files', '', wildcard)
            
                        
            #Check reference files
            if len(filepaths) != self.stk.n_ev:
                QtGui.QMessageBox.warning(self,'Error',"Wrong number of Reference image files.")
                return

            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))                           
            

            self.stk.read_xrm_ReferenceImages(filepaths)
            
            x=self.stk.n_cols
            y=self.stk.n_rows           
            self.ix = x/2
            self.iy = y/2
                        
            self.loadSpectrum(self.ix, self.iy)
            self.loadImage()
            self.com.i0_loaded = 1
            QtGui.QApplication.restoreOverrideCursor()
                          
            
        except:
            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))  
            self.com.i0_loaded = 0        
            QtGui.QApplication.restoreOverrideCursor()
            QtGui.QMessageBox.warning(self,'Error',"Reference image file not loaded.")
            import sys; print sys.exc_info()
                       
                          
        self.window().refresh_widgets()

        
#----------------------------------------------------------------------    
    def OnResetI0(self, event):     

        self.stk.reset_i0()
        
        self.com.i0_loaded = 0
        
        self.showflux = True
        self.rb_flux.setChecked(True)
        
        self.loadSpectrum(self.ix, self.iy)
        self.loadImage()
        self.window().refresh_widgets()
        
        
#----------------------------------------------------------------------    
    def OnSaveStack(self, event): 
        
        self.window().SaveProcessedStack()
               
#----------------------------------------------------------------------    
    def OnSave(self, event):     
        
        savewin = SaveWinP1(self.window())
        savewin.show()

#----------------------------------------------------------------------        
    def OnLimitEv(self, evt):   
        
        limitevwin = LimitEv(self.window(), self.com, self.stk)
        limitevwin.show()

#----------------------------------------------------------------------        
    def OnCliptoSubregion(self, evt):    
        clipwin = CliptoSubregion(self.window(), self.com, self.stk)
        clipwin.show() 
   
#----------------------------------------------------------------------    
    def Save(self, filename, path, spec_png = True, spec_pdf = False, spec_svg = False, sp_csv = False, 
             img_png = True, img_pdf = False, img_svg = False, img_tif = False, img_all = False, img_all_tif = False): 

        self.SaveFileName = os.path.join(path,filename)
      
        
        try: 
            ext = 'png'
            suffix = "." + ext
            matplotlib.rcParams['pdf.fonttype'] = 42
            
            if spec_png:
                fileName_spec = self.SaveFileName+"_spectrum."+ext
                
                fig = self.specfig
                fig.savefig(fileName_spec)

            if img_png:
                
                fileName_img = self.SaveFileName+"_" +str(self.stk.ev[self.iev])+"eV."+ext
                fig = self.absimgfig
                fig.savefig(fileName_img, pad_inches = 0.0)
                
            #Save all images in the stack
            if img_all:
                QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
                matplotlib.rcParams['pdf.fonttype'] = 42

                for i in range (self.stk.n_ev):
                    if self.showflux:
                        #Show flux image      
                        image = self.stk.absdata[:,:,i] 
                    else:
                        #Show OD image
                        image = self.stk.od3d[:,:,i]

                    fig = matplotlib.figure.Figure(figsize =(float(self.stk.n_rows)/10, float(self.stk.n_cols)/10))
                    fig.clf()
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.0,0.0,1.0,1.0))
                    axes = fig.gca()
                    im = axes.imshow(np.rot90(image), cmap=matplotlib.cm.get_cmap(self.colortable)) 
                    axes.axis("off") 
                                
                    fileName_img = self.SaveFileName+"_imnum_" +str(i+1)+"."+ext
                    fig.savefig(fileName_img,  dpi=ImgDpi, pad_inches = 0.0)
                QtGui.QApplication.restoreOverrideCursor()
                
            #Save all images in the stack
            if img_all_tif:
                QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
     
                for i in range (self.stk.n_ev):
                    if self.showflux:
                        #Show flux image      
                        image = self.stk.absdata[:,:,i] 
                    else:
                        #Show OD image
                        image = self.stk.od3d[:,:,i]
                        
                    fileName_img = self.SaveFileName+"_imnum_" +str(i+1)+".tif"
                    img1 = Image.fromarray(image)
                    img1.save(fileName_img)
                                
                QtGui.QApplication.restoreOverrideCursor()
                    
            ext = 'pdf'
            suffix = "." + ext
            
            if spec_pdf:
                fileName_spec = self.SaveFileName+"_spectrum."+ext
                fig = self.specfig
                fig.savefig(fileName_spec)

            if img_pdf:
                fileName_img = self.SaveFileName+"_" +str(self.stk.ev[self.iev])+"eV."+ext
                fig = self.absimgfig
                fig.savefig(fileName_img, bbox_inches='tight', pad_inches = 0.0)
                
            if sp_csv:
                fileName_spec = self.SaveFileName+"_spectrum.csv"
                self.stk.write_csv(fileName_spec, self.stk.ev, self.spectrum)
                
            ext = 'svg'
            suffix = "." + ext
            
            if spec_svg:
                fileName_spec = self.SaveFileName+"_spectrum."+ext
                fig = self.specfig
                fig.savefig(fileName_spec)

            if img_svg:
                fileName_img = self.SaveFileName+"_" +str(self.stk.ev[self.iev])+"eV."+ext
                fig = self.absimgfig
                fig.savefig(fileName_img, bbox_inches='tight', pad_inches = 0.0)
                
            if img_tif:
                fileName_img = self.SaveFileName+"_" +str(self.stk.ev[self.iev])+"eV.tif"
                if self.showflux:     
                    image = self.stk.absdata[:,:,self.iev] 
                else:
                    image = self.stk.od3d[:,:,self.iev]           
                img1 = Image.fromarray(image)
                img1.save(fileName_img)

                
        except IOError, e:
            if e.strerror:
                err = e.strerror 
            else: 
                err = e 
   
            QtGui.QMessageBox.warning(self,'Error','Could not save file: %s' % err)


#----------------------------------------------------------------------    
    def OnSaveOD(self, event): 
        
        """
        Browse for tiff 
        """
        
        try:
        #if True:
            wildcard = "TIFF files (.tif)"

            filepath = QtGui.QFileDialog.getSaveFileName(self, 'Save OD', '', wildcard)

            filepath = str(filepath)
            if filepath == '':
                return
            
            
            directory =  os.path.dirname(str(filepath))
        
            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
            
            
            self.stk.write_tif(filepath, self.stk.od3d) 
                
 
         
            QtGui.QApplication.restoreOverrideCursor()

        except:
     
            QtGui.QApplication.restoreOverrideCursor()
                
            QtGui.QMessageBox.warning(self, 'Error', 'Could not save OD stack file.')
                   

        
        return
        
#----------------------------------------------------------------------    
    def OnAlignImgs(self, event):  
        
        #self.window().Hide()
        imgregwin = ImageRegistration(self.window(), self.com, self.stk)
        imgregwin.show()
        
#----------------------------------------------------------------------    
    def OnDarkSignal(self, event):  
        

        dswin = DarkSignal(self.window(), self.com, self.stk)
        dswin.show()

#----------------------------------------------------------------------    
    def OnSlideshow(self, event):  

        if (self.com.stack_loaded == 1) and (self.addroi == 0):    
            
            if (self.movie_playing == 1):
                self.movie_playing = 0
                return
            
            self.button_slideshow.setText("Stop stack movie")
            old_iev =  self.iev
            self.movie_playing = 1
            
            for i in range(self.stk.n_ev):   
                QtGui.qApp.processEvents()
                if self.movie_playing == 0:
                    break
                self.iev = i                   
                self.loadImage()
                self.slider_eng.setValue(self.iev)
                self.loadSpectrum(self.ix, self.iy)
                
                
            self.iev = old_iev 
            self.loadImage()
            self.slider_eng.setValue(self.iev)
            self.loadSpectrum(self.ix, self.iy)
            self.button_slideshow.setText("Play stack movie")
            self.movie_playing = 0
            

                        
#----------------------------------------------------------------------            
    def OnScrollEng(self, value):
        self.iev = value
        

        if self.com.stack_loaded == 1:
            self.loadImage()
            self.loadSpectrum(self.ix, self.iy)      
            
#----------------------------------------------------------------------            
    def OnScrollTheta(self, value):
        self.itheta = value
        

        self.stk.absdata = self.stk.stack4D[:,:,:,self.itheta]
        if self.com.i0_loaded:
            self.stk.od3d = self.stk.od4D[:,:,:,self.itheta]
            self.stk.od = self.stk.od3d.copy()
            n_pixels = self.stk.n_cols*self.stk.n_rows
            self.stk.od = np.reshape(self.stk.od, (n_pixels, self.stk.n_ev), order='F')        
            
        self.tc_imagetheta.setText("4D Data Angle: "+str(self.stk.theta[self.itheta]))
        
        
        if self.com.stack_loaded == 1:
            self.loadImage()
            self.loadSpectrum(self.ix, self.iy)  
            
            
        self.window().page0.itheta = self.itheta
        self.window().page0.slider_theta.setValue(self.itheta)
        
        self.window().page2.itheta = self.itheta
        self.window().page2.slider_theta.setValue(self.itheta)

#----------------------------------------------------------------------  
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
    def OnPointAbsimage(self, evt):

        
        x = evt.xdata
        y = evt.ydata
        

        if (x == None) or (y == None):
            return
        
        if (self.com.stack_loaded == 1) and (self.addroi == 0):      
            self.ix = int(np.floor(x))           
            self.iy = self.stk.n_rows-1-int(np.floor(y))  
                    
            if self.ix<0 :
                self.ix=0
            if self.ix>self.stk.n_cols-1 :
                self.ix=self.stk.n_cols-1
            if self.iy<0 :
                self.iy=0
            if self.iy>self.stk.n_rows-1 :
                self.iy=self.stk.n_rows-1
            

            self.loadSpectrum(self.ix, self.iy)
            self.loadImage()

            
        if (self.com.stack_loaded == 1) and (self.addroi == 1):
            if self.line == None: # if there is no line, create a line
                self.line = matplotlib.lines.Line2D([x,  x], [y, y], marker = '.', color = 'red')
                self.start_point = [x,y]
                self.previous_point =  self.start_point
                self.roixdata.append(x)
                self.roiydata.append(y)
                self.loadImage()
            # add a segment
            else: # if there is a line, create a segment
                self.roixdata.append(x)
                self.roiydata.append(y)
                self.line.set_data(self.roixdata,self.roiydata)
                self.previous_point = [x,y]
                if len(self.roixdata) == 3:
                    self.button_acceptROI.setEnabled(True)
                self.loadImage()

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
            self.window().page0.ShowImage()
                    
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
        
        colorwin = ColorTableFrame(self.window())
        colorwin.show()
                                    
#----------------------------------------------------------------------        
    def loadImage(self):
        
        
        if (self.defaultdisplay == 1.0):
            #use a pointer to the data not a copy
            if self.showflux:
                #Show flux image      
                image = self.stk.absdata[:,:,self.iev]#.copy() 
            else:
                #Show OD image
                image = self.stk.od3d[:,:,self.iev]#.copy()
        else:   
            #Adjustment to the data display setting has been made so make a copy
            if self.showflux:
                image = self.stk.absdata[:,:,self.iev].copy() 
            else:
                image = self.stk.od3d[:,:,self.iev].copy() 
                 
 
                       
        fig = self.absimgfig
        fig.clf()
         
        if self.show_colorbar == 0:
            fig.add_axes(((0.0,0.0,1.0,1.0)))
            axes = fig.gca()
 
        else:
            axes = fig.gca()
            divider = make_axes_locatable(axes)
            axcb = divider.new_horizontal(size="3%", pad=0.03)  
 
            fig.add_axes(axcb)
         
            axes.set_position([0.03,0.03,0.8,0.94])
             
 
        fig.patch.set_alpha(1.0)
         
        if (self.line != None) and (self.addroi == 1):
            axes.add_line(self.line)
                     
 
        if self.defaultdisplay == 1.0:
            im = axes.imshow(np.rot90( image), cmap=matplotlib.cm.get_cmap(self.colortable))
        else:
            imgmax = np.amax(image)
            imgmin = np.amin(image)
            if (self.gamma != 1.0) or (imgmin < 0.0):
                image = (image-imgmin)/(imgmax-imgmin)
                imgmax = 1.0
                imgmin = 0.0
                if (self.gamma != 1.0):
                    image = np.power(image, self.gamma)
            vmin=(imgmin+imgmax*self.brightness_min)
            vmax=imgmax*self.brightness_max
            if vmin > vmax : vmax = vmin + 0.1
            im = axes.imshow(np.rot90( image), cmap=matplotlib.cm.get_cmap(self.colortable), 
                             vmin=vmin,vmax=vmax)
             

        if (self.showROImask == 1) and (self.addroi == 1):
            im_red = axes.imshow(np.rot90( self.ROIpix_masked), cmap=matplotlib.cm.get_cmap("autumn")) 
        
          
        axes.axis("off") 
         
        if self.show_colorbar == 1:
            cbar = axes.figure.colorbar(im, orientation='vertical',cax=axcb) 
         
        if self.show_scale_bar == 1:
            if self.white_scale_bar == 1:
                sbcolor = 'white'
            else:
                sbcolor = 'black'
            startx = int(self.stk.n_cols*0.05)
            starty = self.stk.n_rows-int(self.stk.n_rows*0.05)-self.stk.scale_bar_pixels_y
            um_string = ' $\mathrm{\mu m}$'
            microns = '$'+self.stk.scale_bar_string+' $'+um_string
            axes.text(self.stk.scale_bar_pixels_x+startx+1,starty+1, microns, horizontalalignment='left', verticalalignment='center',
                      color = sbcolor, fontsize=14)
            #Matplotlib has flipped scales so I'm using rows instead of cols!
            p = matplotlib.patches.Rectangle((startx,starty), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                                   color = sbcolor, fill = True)
            axes.add_patch(p)
         
        self.AbsImagePanel.draw()
         
        self.tc_imageeng.setText("Image at energy: {0:5.2f} eV".format(float(self.stk.ev[self.iev])))


#----------------------------------------------------------------------          
    def loadSpectrum(self, xpos, ypos):
        
        
        fig = self.specfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        
         
        if self.com.i0_loaded == 1:
            self.spectrum = self.stk.od3d[xpos,ypos, :]

            axes.set_xlabel('Photon Energy [eV]')
            axes.set_ylabel('Optical Density')
        else:

            self.spectrum = self.stk.absdata[xpos,ypos, :]
            axes.set_xlabel('Photon Energy [eV]')
            axes.set_ylabel('Flux')
             
 
        specplot = axes.plot(self.stk.ev,self.spectrum)
         
     
        axes.axvline(x=self.stk.ev[self.iev], color = 'g', alpha=0.5)
 
         
        self.SpectrumPanel.draw()
         
#         self.tc_spec.setText("Spectrum at pixel [" +str(ypos)+", " + str(xpos)+"] or position ["+
#                               str(self.stk.x_dist[xpos])+", "+ str(self.stk.y_dist[ypos])+ "]")

        #self.tc_spec.setText('Spectrum at pixel [{0}, {1}] or position [{2:5.2f}, {3:5.2f}]'.format(str(ypos),  str(xpos), self.stk.x_dist[xpos], self.stk.y_dist[ypos]))
#         print("self.stk.x_dist[xpos] = ", type(ypos))
#         print("self.stk.y_dist[ypos] = ", type(xpos))

        self.tc_spec.setText('Spectrum at pixel [{0}, {1}] or position [{2:5.2f}, {3:5.2f}]'.format(str(ypos),  str(xpos), np.float(self.stk.x_dist[xpos]), np.float(self.stk.y_dist[ypos])))

#----------------------------------------------------------------------
    def ResetDisplaySettings(self):
        


        self.defaultdisplay = 1.0
         
        self.dispbrightness_min = 0
        self.dispbrightness_max = 100
        self.displaygamma = 10.0
         
        self.brightness_min = 0.0
        self.brightness_max = 1.0
        self.gamma = 1.0
         
        self.slider_brightness_max.setValue(self.dispbrightness_max)
        self.slider_brightness_min.setValue(self.dispbrightness_min) 
        self.slider_gamma.setValue(self.displaygamma)      
 
        self.tc_min.setText('Minimum: \t{0:5d}%'.format(int(100*self.brightness_min)))
        self.tc_max.setText('Maximum:{0:5d}%'.format(int(100*self.brightness_max)))        
        self.tc_gamma.setText('Gamma:  \t{0:5.2f}'.format(self.gamma))      


                
#----------------------------------------------------------------------         
# Determine if a point is inside a given polygon or not. The algorithm is called
# "Ray Casting Method".
    def point_in_poly(self, x, y, polyx, polyy):

        n = len(polyx)
        inside = False

        p1x = polyx[0]
        p1y = polyy[0]
        for i in range(n+1):
            p2x = polyx[i % n]
            p2y = polyy[i % n]
            if y > min(p1y,p2y):
                if y <= max(p1y,p2y):
                    if x <= max(p1x,p2x):
                        if p1y != p2y:
                            xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                            if p1x == p2x or x <= xinters:
                                inside = not inside
            p1x,p1y = p2x,p2y

        return inside
        
#----------------------------------------------------------------------    
    def OnAddROI(self, evt):    
        self.addroi = 1
        self.previous_point = []
        self.start_point = []
        self.end_point = []
        self.line = None
        self.roixdata = []
        self.roiydata = []
                
        fig = self.specfig
        fig.clf()
        self.SpectrumPanel.draw()
        self.tc_spec.setText("Average ROI Spectrum: ")
        
        self.button_acceptROI.setEnabled(False)
        self.button_resetROI.setEnabled(True)
        self.button_ROIdosecalc.setEnabled(False) 
        self.window().refresh_widgets()

        return
        
        
#----------------------------------------------------------------------    
    def CalcROISpectrum(self):
              
        self.ROIspectrum = np.zeros((self.stk.n_ev))
               
        indices = np.where(self.ROIpix == 255)
        numroipix = self.ROIpix[indices].shape[0]
            
        for ie in range(self.stk.n_ev):
            thiseng_od = self.stk.od3d[:,:,ie]
            self.ROIspectrum[ie] = np.sum(thiseng_od[indices])/numroipix
                 
#----------------------------------------------------------------------          
    def ShowROISpectrum(self):
        
        self.CalcROISpectrum()

        fig = self.specfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        

        specplot = axes.plot(self.stk.ev,self.ROIspectrum)
        
        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Optical Density')
        
        self.SpectrumPanel.draw()
        
        self.tc_spec.setText("Average ROI Spectrum: ")
        
#----------------------------------------------------------------------    
    def OnAcceptROI(self, evt):    
        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        self.roixdata.append(self.start_point[0])
        self.roiydata.append(self.start_point[1])
        self.line.set_data(self.roixdata,self.roiydata)
        #self.loadImage()
        
        #find pixels inside the polygon 
        if self.ROIpix == None:
            self.ROIpix = np.zeros((self.stk.n_cols,self.stk.n_rows))    
        
        for i in range(self.stk.n_cols):
            for j in range(self.stk.n_rows):
                Pinside = self.point_in_poly(i, j, self.roixdata, self.stk.n_rows-1-self.roiydata)
                if Pinside == True:
                    self.ROIpix[i, j] = 255
              


        self.ROIpix = np.ma.array(self.ROIpix)
        
        
        self.ROIpix_masked =  np.ma.masked_values(self.ROIpix, 0)
        
        
        self.showROImask = 1
        self.line = None
        self.previous_point = []
        self.start_point = []
        self.end_point = []
        self.roixdata = []
        self.roiydata = []

        self.button_saveROIspectr.setEnabled(True)
        self.button_setROII0.setEnabled(True)
        self.button_ROIdosecalc.setEnabled(True) 
        self.window().refresh_widgets()
                
        self.loadImage()
        if (self.com.i0_loaded == 1):
            self.ShowROISpectrum()
            
        QtGui.QApplication.restoreOverrideCursor()
        
    
#----------------------------------------------------------------------    
    def OnResetROI(self, evt): 
        self.addroi = 0   
        self.showROImask = 0
        self.ROIpix = None
        
        self.button_acceptROI.setEnabled(False)
        self.button_setROII0.setEnabled(False)
        self.button_resetROI.setEnabled(False)
        self.button_saveROIspectr.setEnabled(False)
        self.button_ROIdosecalc.setEnabled(False) 
        self.window().refresh_widgets()
        
        self.loadImage()
        if (self.com.i0_loaded == 1):
            self.loadSpectrum(self.ix, self.iy)
        pass
    
#----------------------------------------------------------------------    
    def CalcROI_I0Spectrum(self):
   
        self.ROIspectrum = np.zeros((self.stk.n_ev))
               
        indices = np.where(self.ROIpix == 255)
        numroipix = self.ROIpix[indices].shape[0]
            
        for ie in range(self.stk.n_ev):
            thiseng_abs = self.stk.absdata[:,:,ie]
            self.ROIspectrum[ie] = np.sum(thiseng_abs[indices])/numroipix
                
#----------------------------------------------------------------------    
    def OnSetROII0(self, evt):    
        self.CalcROI_I0Spectrum()   
        
        self.stk.set_i0(self.ROIspectrum, self.stk.ev)  
        
        plot = PlotFrame(self, self.stk.evi0,self.stk.i0data)
        plot.show()              
         
        x=self.stk.n_cols
        y=self.stk.n_rows
             
        self.ix = int(x/2)
        self.iy = int(y/2)
        
        self.com.i0_loaded = 1
        
        self.addroi = 0   
        self.showROImask = 0
        self.ROIpix = None
        
        self.loadSpectrum(self.ix, self.iy)
        self.loadImage()
        
        self.button_acceptROI.setEnabled(False)
        self.button_setROII0.setEnabled(False)
        self.button_resetROI.setEnabled(False)
        self.button_saveROIspectr.setEnabled(False)
        self.button_ROIdosecalc.setEnabled(False) 
        
        self.window().refresh_widgets()
        
#----------------------------------------------------------------------    
    def OnROI_DoseCalc(self, event):  
        
        self.CalcROISpectrum()
        
        dosewin = DoseCalculation(self, self.stk, self.ROIspectrum)
        dosewin.show()
        
        
#----------------------------------------------------------------------    
    def OnSaveROISpectrum(self, event):  
        



        wildcard = "CSV files (*.csv)"

        fileName = QtGui.QFileDialog.getSaveFileName(self, 'Save ROI Spectrum (.csv)', '', wildcard)

        fileName = str(fileName)
        if fileName == '':
            return

        path, ext = os.path.splitext(fileName) 
        ext = ext[1:].lower() 
   
        try:
            if (self.com.i0_loaded == 1):
                self.stk.write_csv(fileName, self.stk.ev, self.ROIspectrum, cname='from ROI')
            else:
                self.CalcROI_I0Spectrum()
                self.stk.write_csv(fileName, self.stk.ev, self.ROIspectrum, cname='from ROI')
                     
                
        except IOError, e:
            if e.strerror:
                err = e.strerror 
            else: 
                err = e 
   
            QtGui.QMessageBox.warning(self,'Error','Could not save file: %s' % err)
            

#----------------------------------------------------------------------        
    def OnSpectralROI(self, evt):    
        specroiwin = SpectralROI(self, self.com, self.stk)
        specroiwin.show()     
     
#---------------------------------------------------------------------- 
class SaveWinP1(QtGui.QDialog):

    def __init__(self, parent):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        
        self.resize(400, 300)
        self.setWindowTitle('Save')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)


        self.com = self.parent.common          
        
        path, ext = os.path.splitext(self.com.filename) 
        ext = ext[1:].lower()   
        suffix = "." + ext
        path, fn = os.path.split(self.com.filename)
        filename = fn[:-len(suffix)]
        
        self.path = self.com.path
        self.filename = filename
                          
        
        vboxtop = QtGui.QVBoxLayout()
        vboxtop.setContentsMargins(20,20,20,20)
        
        gridtop = QtGui.QGridLayout()
        gridtop.setVerticalSpacing(20)
    
        fontb = QtGui.QFont()
        fontb.setBold(True)            
        
        st1 = QtGui.QLabel(self)
        st1.setText('Save')
        st1.setFont(fontb)
        st2 = QtGui.QLabel(self)
        st2.setText('.pdf')
        st2.setFont(fontb)
        st3 = QtGui.QLabel(self)
        st3.setText('.png')
        st3.setFont(fontb)
        st4 = QtGui.QLabel(self)
        st4.setText('.svg')
        st4.setFont(fontb)
        st9 = QtGui.QLabel(self)
        st9.setText('.tif (data)')
        st9.setFont(fontb)
        st5 = QtGui.QLabel(self)
        st5.setText('.csv')
        st5.setFont(fontb)        
        
        
        st6 = QtGui.QLabel(self)
        st6.setText('_spectrum')
        
        self.cb11 = QtGui.QCheckBox('', self)
        self.cb11.setChecked(True)
        
        self.cb12 = QtGui.QCheckBox('', self)
        
        self.cb13 = QtGui.QCheckBox('', self)   
        
        self.cb14 = QtGui.QCheckBox('', self)
        
        st7 = QtGui.QLabel(self)
        st7.setText('_image')
        
        self.cb21 = QtGui.QCheckBox('', self)
        self.cb21.setChecked(True)
        
        self.cb22 = QtGui.QCheckBox('', self)
        
        self.cb23 = QtGui.QCheckBox('', self)
        
        self.cb24 = QtGui.QCheckBox('', self)
        
        st8 = QtGui.QLabel(self)
        st8.setText('all images')   
        
        self.cb32 = QtGui.QCheckBox('', self)
        
        self.cb34 = QtGui.QCheckBox('', self)


        gridtop.addWidget(st1, 0, 0)
        gridtop.addWidget(st2, 0, 1)
        gridtop.addWidget(st3, 0, 2)
        gridtop.addWidget(st4, 0, 3)
        gridtop.addWidget(st9, 0, 5)
        gridtop.addWidget(st5, 0, 4)
                                
        gridtop.addWidget(st6, 1, 0)
        gridtop.addWidget(self.cb11, 1, 1)
        gridtop.addWidget(self.cb12, 1, 2)  
        gridtop.addWidget(self.cb13, 1, 3)  
        gridtop.addWidget(self.cb14, 1, 4)           
  
        gridtop.addWidget(st7, 2, 0)
        gridtop.addWidget(self.cb21, 2, 1)
        gridtop.addWidget(self.cb22, 2, 2) 
        gridtop.addWidget(self.cb23, 2, 3) 
        gridtop.addWidget(self.cb24, 2, 5)  
        
        gridtop.addWidget(st8, 3, 0)
        gridtop.addWidget(self.cb32, 3, 2)  
        gridtop.addWidget(self.cb34, 3, 5)  
        
        vboxtop.addStretch(0.5)
        vboxtop.addLayout(gridtop)
        vboxtop.addStretch(1)
        
         
        hbox0 = QtGui.QHBoxLayout()
         
        stf = QtGui.QLabel(self)
        stf.setText('Filename:\t')
        self.tc_savefn = QtGui.QLineEdit(self)
        self.tc_savefn.setText(self.filename)

        hbox0.addWidget(stf)
        hbox0.addWidget(self.tc_savefn)         
                  
        hbox1 = QtGui.QHBoxLayout()
                 
        stp = QtGui.QLabel(self)
        stp.setText('Path:  \t')
        self.tc_savepath = QtGui.QLineEdit(self)
        self.tc_savepath.setReadOnly(True)
        self.tc_savepath.setText(self.path)
        self.tc_savepath.setMinimumWidth(100)
        hbox1.addWidget(stp)
        hbox1.addWidget(self.tc_savepath)  
         
        button_path = QtGui.QPushButton('Browse...')
        button_path.clicked.connect(self.OnBrowseDir)
        hbox1.addWidget(button_path)
         
         
        hbox2 = QtGui.QHBoxLayout()
        button_save = QtGui.QPushButton('Save')
        button_save.clicked.connect(self.OnSave)
        hbox2.addWidget(button_save)
         
        button_cancel = QtGui.QPushButton('Cancel')
        button_cancel.clicked.connect(self.close)
        hbox2.addWidget(button_cancel)
        
        vboxtop.addLayout(hbox0)
        vboxtop.addLayout(hbox1)
        vboxtop.addStretch(1.0)
        vboxtop.addLayout(hbox2)

        
        
        self.setLayout(vboxtop)
        
#----------------------------------------------------------------------        
    def OnBrowseDir(self, evt):
        
        directory = QtGui.QFileDialog.getExistingDirectory(self, "Choose a directory", self.path, QtGui.QFileDialog.ShowDirsOnly|QtGui.QFileDialog.ReadOnly)       
                                                        
        
       
        if directory == '':
            return
                 
        directory = str(directory)
        self.com.path = directory
                    
        self.path = directory
        
        self.tc_savepath.setText(self.path)
            
            
                
#----------------------------------------------------------------------        
    def OnSave(self, evt):
        
        self.filename = str(self.tc_savefn.text())
        
        sp_pdf = self.cb11.isChecked()
        sp_png = self.cb12.isChecked()
        sp_svg = self.cb13.isChecked()
        sp_csv = self.cb14.isChecked()
        im_pdf = self.cb21.isChecked()
        im_png = self.cb22.isChecked()
        im_svg = self.cb23.isChecked()
        im_tif = self.cb24.isChecked()
        im_all = self.cb32.isChecked()
        im_all_tif = self.cb34.isChecked()
        
        self.close() 
        self.parent.page1.Save(self.filename, self.path,
                                         spec_png = sp_png, 
                                         spec_pdf = sp_pdf, 
                                         spec_svg = sp_svg,
                                         sp_csv = sp_csv,
                                         img_png = im_png, 
                                         img_pdf = im_pdf,
                                         img_svg = im_svg,
                                         img_tif = im_tif,
                                         img_all = im_all,
                                         img_all_tif = im_all_tif)




#---------------------------------------------------------------------- 
class ShowHistogram(QtGui.QDialog):

    def __init__(self, parent, stack):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        
        self.resize(600, 600)
        self.setWindowTitle('Histogram')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)
                
        self.stack = stack
      
        
        self.stack.calc_histogram()
        averagefluxmax = np.max(self.stack.histogram)
        self.histmin = 0.98*averagefluxmax
        self.histmax = averagefluxmax
        

        vbox = QtGui.QVBoxLayout()
               
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
   
        self.histfig = Figure((6.0, 4.2))
        self.HistogramPanel = FigureCanvas(self.histfig)
        self.HistogramPanel.setParent(self)
        self.HistogramPanel.mpl_connect('button_press_event', self.OnSelection1)
        self.HistogramPanel.mpl_connect('button_release_event', self.OnSelection2)
        
        
        fbox.addWidget(self.HistogramPanel)
        frame.setLayout(fbox)
        vbox.addWidget(frame)
                        
        vbox1 = QtGui.QVBoxLayout()
        sizer1 = QtGui.QGroupBox('I0 pixels')

        self.textctrl = QtGui.QLabel(self)
        self.textctrl.setText('Selection: [ {0:5.2f} kHz, {1:5.2f} kHz ]'.format(float(self.histmin), float(self.histmax)))
        vbox1.addWidget(self.textctrl)

 
        self.absimgfig = Figure((2.5,2.5))
        self.AbsImagePanel = FigureCanvas(self.absimgfig)
        self.AbsImagePanel.setParent(self)

        vbox1.addWidget(self.AbsImagePanel,0,QtCore .Qt. AlignLeft)
        sizer1.setLayout(vbox1)
        vbox.addWidget(sizer1)
                
        hbox2 = QtGui.QHBoxLayout()
        button_ok = QtGui.QPushButton('Accept')
        button_ok.clicked.connect(self.OnAccept)
        hbox2.addWidget(button_ok)
                
        button_cancel = QtGui.QPushButton('Cancel')
        button_cancel.clicked.connect(self.close)
        hbox2.addWidget(button_cancel)
        
        vbox.addLayout(hbox2)
        
        self.setLayout(vbox)
        
        self.draw_histogram()
        self.draw_image()


        
#----------------------------------------------------------------------        
    def draw_histogram(self):
        
     
        fig = self.histfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        self.axes = fig.gca()
        
        
        histdata =  np.reshape(self.stack.histogram, (self.stack.n_cols*self.stack.n_rows), order='F')
        
        self.n, self.bins, patches = self.axes.hist(histdata, 200, normed=1, facecolor='green', alpha=0.75)
        
        self.axes.set_xlabel('Average Flux [kHz]')
        self.axes.set_ylabel('Percentage of Pixels')

        self.patch = self.axes.axvspan(self.histmin, self.histmax, facecolor='r', alpha=0.3)

        
        self.HistogramPanel.draw()
        
        
    
    
#----------------------------------------------------------------------        
    def draw_image(self):
        
   
        image = self.stack.absdata[:,:,self.stack.n_ev/2].copy() 
       
        fluxmin = self.histmin
        fluxmax = self.histmax
                
        hist_indices = np.where((fluxmin<self.stack.averageflux)&(self.stack.averageflux<fluxmax))

        redpix = np.zeros((self.stack.n_cols,self.stack.n_rows))    
        redpix[hist_indices] = 255
              
        redpix = np.ma.array(redpix)
        
        redpix_masked =  np.ma.masked_values(redpix, 0)
        
               
        fig = self.absimgfig
        fig.clf()
        fig.add_axes(((0.0,0.0,1.0,1.0)))
        
        
        axes = fig.gca()
        fig.patch.set_alpha(1.0) 
      
        im = axes.imshow(np.rot90(image), cmap=matplotlib.cm.get_cmap("gray")) 

        im_red = axes.imshow(np.rot90(redpix_masked),cmap=matplotlib.cm.get_cmap("autumn"))  
         
        axes.axis("off")  
        self.AbsImagePanel.draw()

#----------------------------------------------------------------------        
    def OnSelection1(self, evt):
        
        x1 = evt.xdata
                
        self.button_pressed = True
        self.conn = self.HistogramPanel.mpl_connect('motion_notify_event', self.OnSelectionMotion) 
        
        if x1 == None:
            return
        
        self.histmin = x1

    

        
#----------------------------------------------------------------------        
    def OnSelection2(self, evt):
        
        x2 = evt.xdata
            

        self.button_pressed = False
        self.HistogramPanel.mpl_disconnect(self.conn)    
        
        if x2 == None:
            return        
        
        self.histmax = x2

        self.textctrl.setText('Selection: [ {0:5.2f} kHz, {1:5.2f} kHz ]'.format(float(self.histmin), float(self.histmax)))
    
        
        self.draw_histogram()
        self.draw_image()       

        
#----------------------------------------------------------------------        
    def OnSelectionMotion(self, event):        

        x2 = event.xdata
        
        if x2 == None:
            return  
        self.histmax = x2
        
        fig = self.histfig

        axes = fig.gca()

        if self.patch != None:
            self.patch.remove()
        self.patch = axes.axvspan(self.histmin, self.histmax, facecolor='white', alpha=0.3)

        self.HistogramPanel.draw()        
        
        
#----------------------------------------------------------------------        
    def OnAccept(self, evt):
        
        self.stack.i0_from_histogram(self.histmin, self.histmax)
        self.parent.I0histogramCalculated()
        self.close()




#---------------------------------------------------------------------- 
class LimitEv(QtGui.QDialog):

    def __init__(self, parent,  common, stack):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        self.stack = stack
        self.com = common  
        
        self.resize(630, 450)
        self.setWindowTitle('Limit energy range')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)

        
        self.evlimited = 0
        self.limitevmin = 0
        self.limitevmax = self.stack.n_ev-1
        
        self.patch = None
        
        vbox = QtGui.QVBoxLayout()
        
        
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
   
        self.specfig = Figure((6.0, 4.2))
        self.SpectrumPanel = FigureCanvas(self.specfig)
        self.SpectrumPanel.setParent(self)
        self.SpectrumPanel.mpl_connect('button_press_event', self.OnSelection1)
        self.SpectrumPanel.mpl_connect('button_release_event', self.OnSelection2)
        
        
        fbox.addWidget(self.SpectrumPanel)
        frame.setLayout(fbox)
        vbox.addWidget(frame)
                
       
        hbox2 = QtGui.QHBoxLayout()
        sizer2 = QtGui.QGroupBox('Energy')
        self.textctrl = QtGui.QLabel(self)
        self.textctrl.setText(' ')
        hbox2.addWidget(self.textctrl, 0)
        sizer2.setLayout(hbox2)
        vbox.addWidget(sizer2)
        
        hbox = QtGui.QHBoxLayout()
        
        
        button_ok = QtGui.QPushButton('Accept')
        button_ok.clicked.connect(self.OnAccept)
        hbox.addWidget(button_ok)
        
        button_cancel = QtGui.QPushButton('Cancel')
        button_cancel.clicked.connect(self.close)
        hbox.addWidget(button_cancel)
        
        vbox.addLayout(hbox)
        
        self.setLayout(vbox)
        
        self.draw_limitev_plot()
        
      
#----------------------------------------------------------------------        
    def draw_limitev_plot(self):
        
        
        if self.com.i0_loaded == 1: 
            odtotal = self.stack.od3d.sum(axis=0)   
        else:
            odtotal = self.stack.absdata.sum(axis=0)   
        
        odtotal = odtotal.sum(axis=0)/(self.stack.n_rows*self.stack.n_cols)             
        
        fig = self.specfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        self.axes = fig.gca()
        

        specplot = self.axes.plot(self.stack.ev,odtotal)
        
        self.axes.set_xlabel('Photon Energy [eV]')
        self.axes.set_ylabel('Optical Density')
        
        if self.evlimited == 1:
            self.patch = self.axes.axvspan(self.stack.ev[self.limitevmin], self.stack.ev[self.limitevmax], facecolor='g', alpha=0.5)
        
        
        self.SpectrumPanel.draw()
        
        self.textctrl.setText('Min energy {0:5.2f} eV\n'.format(float(self.stack.ev[self.limitevmin]))
                              + 'Max energy {0:5.2f} eV'.format(float(self.stack.ev[self.limitevmax])))

#----------------------------------------------------------------------        
    def OnSelection1(self, evt):


        x1 = evt.xdata
                
        self.button_pressed = True
        self.conn = self.SpectrumPanel.mpl_connect('motion_notify_event', self.OnSelectionMotion) 
        
        if x1 == None:
            return
                
        
        self.limitevmin = np.abs(self.stack.ev-x1).argmin()

             
        
#----------------------------------------------------------------------        
    def OnSelection2(self, evt):

        x2 = evt.xdata
            

        self.button_pressed = False
        self.SpectrumPanel.mpl_disconnect(self.conn)    
        
        if x2 == None:
            return    
        
                
        #self.limitevmin = np.abs(self.stack.ev-x1).argmin()
        self.limitevmax = np.abs(self.stack.ev-x2).argmin()
        
        self.evlimited = 1
             
        self.draw_limitev_plot()
        
        
#----------------------------------------------------------------------        
    def OnSelectionMotion(self, event):        

        x2 = event.xdata
        
        if x2 == None:
            return  
        
        self.limitevmax = np.abs(self.stack.ev-x2).argmin()
        
        fig = self.specfig

        axes = fig.gca()

        if self.patch != None:
            self.patch.remove()
        self.patch = self.axes.axvspan(self.stack.ev[self.limitevmin], self.stack.ev[self.limitevmax], facecolor='w', alpha=0.5)

        self.SpectrumPanel.draw()      
        
        
#----------------------------------------------------------------------        
    def OnAccept(self, evt):
        #change the energy range to limitevmin-limitev-max  
        #print self.stack.n_ev, self.stack.ev.shape
        self.stack.n_ev = self.limitevmax+1-self.limitevmin
        self.stack.ev = self.stack.ev[self.limitevmin:self.limitevmax+1]
        self.stack.data_dwell = self.stack.data_dwell[self.limitevmin:self.limitevmax+1]
        
        #print self.stack.n_ev, self.stack.ev.shape
        
        
        self.stack.absdata = self.stack.absdata[:,:,self.limitevmin:self.limitevmax+1]
        
        if self.com.i0_loaded == 1:   
            self.stack.od3d = self.stack.od3d[:,:,self.limitevmin:self.limitevmax+1]
        
            self.stack.od = self.stack.od3d.copy()
        
            self.stack.od = np.reshape(self.stack.od, (self.stack.n_rows*self.stack.n_cols, self.stack.n_ev), order='F')
        
        self.stack.fill_h5_struct_from_stk()
        if self.com.i0_loaded == 1: 
            self.stack.fill_h5_struct_normalization()
            
            
        if self.com.stack_4d == 1:
            self.stack.stack4D = self.stack.stack4D[:,:,self.limitevmin:self.limitevmax+1,:]
            if self.com.i0_loaded == 1:   
                self.stack.od4D = self.stack.od4D[:,:,self.limitevmin:self.limitevmax+1,:]            

        
        #Fix the slider on Page 1! 
        self.parent.page1.slider_eng.setRange(0,self.stack.n_ev-1)
        self.parent.page1.iev = self.stack.n_ev/2
        self.parent.page1.slider_eng.setValue(self.parent.page1.iev)
        
        self.parent.page0.slider_eng.setRange(0,self.stack.n_ev-1)
        self.parent.page0.iev = self.stack.n_ev/2
        self.parent.page0.slider_eng.setValue(self.parent.page1.iev)
        
        self.parent.page1.loadSpectrum(self.parent.page1.ix, self.parent.page1.iy)
        self.parent.page1.loadImage()
        
        self.close()
        
        
#---------------------------------------------------------------------- 
class CliptoSubregion(QtGui.QDialog):

    def __init__(self, parent,  common, stack):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        self.stack = stack
        self.com = common
        
        self.resize(500, 470)
        self.setWindowTitle('Clip to Subregion')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)
        
        self.new_x1 = int(self.stack.n_cols*0.10)
        self.new_x2 = int(self.stack.n_cols*0.90)
        self.new_y2 = self.stack.n_rows-1-int(self.stack.n_rows*0.10)
        self.new_y1 = self.stack.n_rows-1-int(self.stack.n_rows*0.90)
        
        self.new_ncols = self.new_x2 - self.new_x1
        self.new_nrows = self.new_y2 - self.new_y1
        

        
        vbox = QtGui.QVBoxLayout()             
        
        sizer = QtGui.QGroupBox('Select new stack size')
        vbox1 = QtGui.QVBoxLayout()   
        self.textctrl1 = QtGui.QLabel(self)
        self.textctrl1.setText('Original stack size:\t{0:5d}   x{1:5d} '.format(self.stack.n_cols, self.stack.n_rows))
        vbox1.addWidget(self.textctrl1)
 
        self.textctrl2 = QtGui.QLabel(self)
        self.textctrl2.setText('New stack size:\t{0:5d}   x{1:5d} '.format(self.new_ncols, self.new_nrows))
        vbox1.addWidget(self.textctrl2)
        
        self.textctrl3 = QtGui.QLabel(self)
        self.textctrl3.setText('Clip coordinates [[x1, x2], [y1, y2]] : [[{0:5d},{1:5d}], [{2:5d},{3:5d}]]'.format(
                                    self.new_x1, self.new_x2, self.new_y1, self.new_y2))
        vbox1.addWidget(self.textctrl3)  
        
        self.absimgfig = Figure((PlotH,PlotH))
        self.AbsImagePanel = FigureCanvas(self.absimgfig)
        self.AbsImagePanel.setParent(self)
        self.AbsImagePanel.mpl_connect('button_press_event', self.OnSelection1)
        self.AbsImagePanel.mpl_connect('button_release_event', self.OnSelection2)                       
        
        vbox1.addWidget(self.AbsImagePanel)
        sizer.setLayout(vbox1)
            
        vbox.addWidget(sizer)
         
         
        hbox = QtGui.QHBoxLayout() 
              
        button_ok = QtGui.QPushButton('Accept')
        button_ok.clicked.connect(self.OnAccept)
        hbox.addWidget(button_ok)
        
        button_cancel = QtGui.QPushButton('Cancel')
        button_cancel.clicked.connect(self.close)
        hbox.addWidget(button_cancel)
        
        vbox.addLayout(hbox)
        
        self.setLayout(vbox)

        self.draw_image()


    
    
#----------------------------------------------------------------------        
    def draw_image(self):
        
   
        image = self.stack.absdata[:,:,self.stack.n_ev/2].copy() 
       
       
               
        fig = self.absimgfig
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        
        
        axes = fig.gca()
        fig.patch.set_alpha(1.0) 
      
        im = axes.imshow(np.rot90(image), cmap=matplotlib.cm.get_cmap("gray")) 
        
        # Draw the rectangle
        line1=matplotlib.lines.Line2D([self.new_x1,self.new_x2], [self.new_y1,self.new_y1] ,color="red")
        line1.set_clip_on(False)
        self.l1 = axes.add_line(line1)   

        line2=matplotlib.lines.Line2D([self.new_x1,self.new_x2], [self.new_y2,self.new_y2] ,color="red")
        line2.set_clip_on(False)
        self.l2 = axes.add_line(line2)   
        
        line3=matplotlib.lines.Line2D([self.new_x1,self.new_x1], [self.new_y1,self.new_y2] ,color="red")
        line3.set_clip_on(False)
        self.l3 = axes.add_line(line3)   
        
        line4=matplotlib.lines.Line2D([self.new_x2,self.new_x2], [self.new_y1,self.new_y2] ,color="red")
        line4.set_clip_on(False)
        self.l4 = axes.add_line(line4)   
         
        axes.axis("off")  
        self.AbsImagePanel.draw()


#----------------------------------------------------------------------        
    def OnSelection1(self, evt):
        
        x1, y1 = evt.xdata, evt.ydata

        self.button_pressed = True
        self.conn = self.AbsImagePanel.mpl_connect('motion_notify_event', self.OnSelectionMotion) 
        
        if (x1 == None) or (y1 == None):
            return
                
        self.new_y1 = int(y1)
        self.new_x1 = int(x1)

        
        self.new_ncols = self.new_x2 - self.new_x1 + 1
        self.new_nrows = self.new_y1 - self.new_y2 + 1

#         self.textctrl2.SetValue('New stack size:\t{0:5d}   x{1:5d} '.format(self.new_ncols, self.new_nrows))
#         self.textctrl3.SetValue('Clip coordinates [[x1, x2], [y1, y2]]:[[{0:5d},{1:5d}], [{2:5d},{3:5d}]]'.format(
#                                     self.new_x1, self.new_x2, self.new_y1, self.new_y2))
#     
#         self.draw_image()

#----------------------------------------------------------------------        
    def OnSelection2(self, evt):
        
        x2, y2 = evt.xdata, evt.ydata
        
        self.button_pressed = False
        self.AbsImagePanel.mpl_disconnect(self.conn) 
        
        if (x2 == None) or (y2 == None):
            return
        
        self.new_y2 = int(y2)
        self.new_x2 = int(x2)
        
        if self.new_x1 > self.new_x2:
            temp = self.new_x1
            self.new_x1 = self.new_x2
            self.new_x2 = temp
        
        if self.new_y1 > self.new_y2:
            temp = self.new_y1
            self.new_y1 = self.new_y2
            self.new_y2 = temp
        
        self.new_ncols = self.new_x2 - self.new_x1 + 1
        self.new_nrows = self.new_y2 - self.new_y1 + 1

        self.textctrl2.setText('New stack size:\t{0:5d}   x{1:5d} '.format(self.new_ncols, self.new_nrows))
        self.textctrl3.setText('Clip coordinates [[x1, x2], [y1, y2]]:[[{0:5d},{1:5d}], [{2:5d},{3:5d}]]'.format(
                                    self.new_x1, self.new_x2, self.new_y1, self.new_y2))
    
        self.draw_image()

#----------------------------------------------------------------------        
    def OnSelectionMotion(self, event):        

        x2, y2 = event.xdata, event.ydata
        
        if x2 == None:
            return  
        
        self.new_y2 = int(y2)
        self.new_x2 = int(x2)
        
        fig = self.absimgfig

        axes = fig.gca()
        
        self.l1.remove()
        self.l2.remove()
        self.l3.remove()
        self.l4.remove()

        line1=matplotlib.lines.Line2D([self.new_x1,self.new_x2], [self.new_y1,self.new_y1] ,color="red")
        line1.set_clip_on(False)
        self.l1 = axes.add_line(line1)   

        line2=matplotlib.lines.Line2D([self.new_x1,self.new_x2], [self.new_y2,self.new_y2] ,color="red")
        line2.set_clip_on(False)
        self.l2 = axes.add_line(line2)   
        
        line3=matplotlib.lines.Line2D([self.new_x1,self.new_x1], [self.new_y1,self.new_y2] ,color="red")
        line3.set_clip_on(False)
        self.l3 = axes.add_line(line3)   
        
        line4=matplotlib.lines.Line2D([self.new_x2,self.new_x2], [self.new_y1,self.new_y2] ,color="red")
        line4.set_clip_on(False)
        self.l4 = axes.add_line(line4)   
        

        self.AbsImagePanel.draw()     
        
              
#----------------------------------------------------------------------        
    def OnAccept(self, evt):
        
        #change the stack size to [x1,x2], [y1,y2] 
        
        self.stack.absdata = self.stack.absdata[ self.new_x1:self.new_x2+1, self.new_y1:self.new_y2+1, :]
        
              
        self.stack.n_cols = self.stack.absdata.shape[0]
        self.stack.n_rows = self.stack.absdata.shape[1]

        
        if self.com.i0_loaded == 1:        
            self.stack.od3d = self.stack.od3d[ self.new_x1:self.new_x2+1, self.new_y1:self.new_y2+1, :]
        
            self.stack.od = self.stack.od3d.copy()
        
            self.stack.od = np.reshape(self.stack.od, (self.stack.n_rows*self.stack.n_cols, self.stack.n_ev), order='F')

        
        #Fix the slider on Page 1! 
        self.parent.page1.ix = self.stack.n_cols/2
        self.parent.page1.iy = self.stack.n_rows/2
        
        self.stack.fill_h5_struct_from_stk()
        
        self.parent.page1.loadSpectrum(self.parent.page1.ix, self.parent.page1.iy)
        self.parent.page1.loadImage()
        self.parent.page0.ShowImage()

        self.close()


#---------------------------------------------------------------------- 
class ImageRegistration(QtGui.QDialog):

    def __init__(self, parent,  common, stack):    
        QtGui.QWidget.__init__(self, parent)
        
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
        vbox1 = QtGui.QVBoxLayout()
        
        self.tc_imageeng = QtGui.QLabel(self)
        self.tc_imageeng.setText("Image at energy: ")
       
        gridsizertop = QtGui.QGridLayout()
        
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
   
        self.absimgfig = Figure((4.0, 4.0))
        self.AbsImagePanel = FigureCanvas(self.absimgfig)
        self.AbsImagePanel.setParent(self)
        self.AbsImagePanel.mpl_connect('button_press_event', self.OnPointCorrimage)
        self.AbsImagePanel.setCursor(Qt.CrossCursor)

        
        fbox.addWidget(self.AbsImagePanel)
        frame.setLayout(fbox)
        gridsizertop.addWidget(frame, 0, 0, QtCore .Qt. AlignLeft)
        

        self.slider_eng = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slider_eng.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_eng.valueChanged[int].connect(self.OnScrollEng)
        self.slider_eng.setRange(0, self.stack.n_ev-1)
        self.slider_eng.setValue(self.iev)
        
        gridsizertop.addWidget(self.slider_eng, 0, 1, QtCore .Qt. AlignLeft)      
        
        self.slider_theta = QtGui.QScrollBar(QtCore.Qt.Horizontal)
        self.slider_theta.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_theta.valueChanged[int].connect(self.OnScrollTheta)
        self.slider_theta.setRange(0, self.stack.n_theta-1)     
          
        self.tc_imagetheta = QtGui.QLabel(self)
        self.tc_imagetheta.setText("4D Data Angle: ")
        if self.com.stack_4d == 0 : 
            self.tc_imagetheta.setVisible(False)
            self.slider_theta.setVisible(False)
        hbox51 = QtGui.QHBoxLayout()
        hbox51.addWidget(self.tc_imagetheta) 
        hbox51.addWidget(self.slider_theta)
        gridsizertop.addLayout(hbox51, 1, 0) 
      
        vbox1.addWidget(self.tc_imageeng)        
        vbox1.addLayout(gridsizertop)

        self.tc_shift1 = QtGui.QLabel(self)
        self.tc_shift2 = QtGui.QLabel(self)
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
        vbox2 = QtGui.QVBoxLayout()
        
        tc2 = QtGui.QLabel(self)
        tc2.setText('Cross-correlation')

        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
   
        self.cscorrfig = Figure((2.4, 2.4))
        self.CscorrPanel = FigureCanvas(self.cscorrfig)
        self.CscorrPanel.setParent(self)
        
        fbox.addWidget(self.CscorrPanel)
        frame.setLayout(fbox)
        
        vbox2.addWidget(tc2)        
        vbox2.addWidget(frame)

        
        
        #panel 3
        vbox3 = QtGui.QVBoxLayout()
         
        tc3= QtGui.QLabel(self)
        tc3.setText('Image shifts')
         
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
    
        self.shiftsfig = Figure((4.0, 2.4))
        self.ShiftsPanel = FigureCanvas(self.shiftsfig)
        self.ShiftsPanel.setParent(self)
        self.ShiftsPanel.mpl_connect('button_press_event', self.OnPlotShifts)
         
        fbox.addWidget(self.ShiftsPanel)
        frame.setLayout(fbox)
           
        vbox3.addWidget(tc3)       
        vbox3.addWidget(frame)
         
         
         
        #panel 9
        vbox9 = QtGui.QVBoxLayout()
        vbox9.setSpacing(0)
         
        groupBox9 = QtGui.QGroupBox()
        self.rb_auto = QtGui.QRadioButton( 'Automatic Alignment')
        self.rb_man = QtGui.QRadioButton('Manual Alignment')
        self.rb_auto.setChecked(True)
        self.rb_auto.toggled.connect(self.Onrb_automanual)
         
        vbox9.addWidget(self.rb_auto)
        vbox9.addWidget(self.rb_man)
        groupBox9.setLayout(vbox9)
 
                
         
        #panel 8
        sizer8 = QtGui.QGroupBox('This Image')
        vbox8 = QtGui.QVBoxLayout()
        vbox8.setSpacing(0)
 
         
        self.button_refimg = QtGui.QPushButton('Set as Reference Image')
        self.button_refimg.clicked.connect(self.SetRefImage)
        vbox8.addWidget(self.button_refimg)
        
        self.button_refimgsave = QtGui.QPushButton('Save Reference Image')
        self.button_refimgsave.clicked.connect(self.SaveRefImage)
        self.button_refimgsave.setEnabled(False)
        vbox8.addWidget(self.button_refimgsave)
        
        self.button_refimgsload = QtGui.QPushButton('Load Reference Image')
        self.button_refimgsload.clicked.connect(self.LoadRefImage)
        vbox8.addWidget(self.button_refimgsload)
        
        line = QtGui.QFrame()
        line.setFrameShape(QtGui.QFrame.HLine)
        line.setFrameShadow(QtGui.QFrame.Sunken) 
        
        
        vbox8.addSpacing(5)
        vbox8.addWidget(line) 
        vbox8.addSpacing(5)  
        
         
        self.button_remove = QtGui.QPushButton('Remove image from stack')
        self.button_remove.clicked.connect(self.OnRemoveImage)    
        vbox8.addWidget(self.button_remove)
               
        sizer8.setLayout(vbox8)
 
            
         
        #panel 4
        sizer4 = QtGui.QGroupBox('Automatic Alignment')
        vbox4 = QtGui.QVBoxLayout()
   
         
        self.button_register = QtGui.QPushButton('Calculate image shifts')
        self.button_register.clicked.connect(self.OnCalcRegistration)   
        self.button_register.setEnabled(False)     
        vbox4.addWidget(self.button_register)
        #vbox4.addStretch(1)
         
        self.button_subregion = QtGui.QPushButton('Select subregion on reference')
        self.button_subregion.clicked.connect(self.OnSelectSubregion)   
        self.button_subregion.setEnabled(False)     
        vbox4.addWidget(self.button_subregion)
          
        self.button_delsubregion = QtGui.QPushButton('Remove subregion selection')
        self.button_delsubregion.clicked.connect(self.OnDeleteSubregion)   
        self.button_delsubregion.setEnabled(False)     
        vbox4.addWidget(self.button_delsubregion)
        vbox4.addStretch(1)

        groupBox4 = QtGui.QGroupBox()      
        hbox43 = QtGui.QHBoxLayout()    
        self.cb_edgeenh = QtGui.QCheckBox('Edge Enhancement', self) 
        self.cb_edgeenh.stateChanged.connect(self.OnEdgeE)
        hbox43.addWidget(self.cb_edgeenh)
        
        self.rb_sobel = QtGui.QRadioButton( 'Sobel')
        self.rb_prewitt = QtGui.QRadioButton('Prewitt')
        self.rb_prewitt.setChecked(True)
        self.rb_sobel.setEnabled(False)
        self.rb_prewitt.setEnabled(False)
         
        hbox43.addWidget(self.rb_prewitt) 
        hbox43.addWidget(self.rb_sobel)
        groupBox4.setLayout(hbox43)
        
        vbox4.addWidget(groupBox4)
        
          
        hbox42 = QtGui.QHBoxLayout()
        text1 = QtGui.QLabel(self)
        text1.setText(' Max shift [pixels]: ')
          
        self.tc_maxshift = QtGui.QSpinBox()
        self.tc_maxshift.setMinimum(0)    
        self.tc_maxshift.valueChanged[int].connect(self.OnSetMaxShift)
          
      
        hbox42.addWidget(text1)
        hbox42.addWidget(self.tc_maxshift)
        vbox4.addLayout(hbox42)
        vbox4.addStretch(1)
          
  
        self.showcscor_cb = QtGui.QCheckBox('Show Cross-correlation', self) 
        self.showcscor_cb.stateChanged.connect(self.OnShowCCorr)
        vbox4.addWidget(self.showcscor_cb)
 
        vbox4.addStretch()
        sizer4.setLayout(vbox4)
         
         
        #panel 5        
        vbox5 = QtGui.QVBoxLayout()
         
        self.tc_refimg = QtGui.QLabel(self)
        self.tc_refimg.setText('Reference image')
         
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
    
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
        sizer6 = QtGui.QGroupBox('Manual Alignment')
        vbox6 = QtGui.QVBoxLayout()
        vbox6.setSpacing(0)
             
                 
        self.button_manalign = QtGui.QPushButton('Pick a point on reference image')
        self.button_manalign.clicked.connect(self.OnPickRefPoint)
        self.button_manalign.setEnabled(False)
        vbox6.addWidget(self.button_manalign)
         
        self.button_pick2ndpoint = QtGui.QPushButton('This image: click on same point')
        self.button_pick2ndpoint.clicked.connect(self.OnPickCorrPoint)
        self.button_pick2ndpoint.setEnabled(False)
        vbox6.addWidget(self.button_pick2ndpoint)
         
        self.button_applyman = QtGui.QPushButton('Apply manual shifts')
        self.button_applyman.clicked.connect(self.OnApplyManShifts)
        self.button_applyman.setEnabled(False)
        vbox6.addWidget(self.button_applyman)
         
         
        self.textctrl_ms1 = QtGui.QLabel(self)
        self.textctrl_ms2 = QtGui.QLabel(self)
        vbox6.addWidget(self.textctrl_ms1)
        vbox6.addWidget(self.textctrl_ms2)
         
         
        self.textctrl_ms1.setText('X manual shift: ')
        self.textctrl_ms2.setText('Y manual shift: ')
               
         
        sizer6.setLayout(vbox6)
         
         
         
        #panel 7
        vbox7 = QtGui.QVBoxLayout()
        vbox7.setSpacing(0)
         
        self.button_saveimg = QtGui.QPushButton('Save image shifts plot')
        self.button_saveimg.clicked.connect(self.OnSaveShiftsPlot)   
        self.button_saveimg.setEnabled(False)
        vbox7.addWidget(self.button_saveimg)
         

 
        self.button_saveshifts = QtGui.QPushButton('Save image shifts')
        self.button_saveshifts.clicked.connect(self.OnSaveShifts)
        self.button_saveshifts.setEnabled(False)
        vbox7.addWidget(self.button_saveshifts)
         
        self.button_loadshifts = QtGui.QPushButton('Load image shifts')
        self.button_loadshifts.clicked.connect(self.OnLoadShifts)
        vbox7.addWidget(self.button_loadshifts)
        
        self.button_crop = QtGui.QPushButton('Crop aligned images')
        self.button_crop.clicked.connect(self.OnCropShifts)   
        self.button_crop.setEnabled(False)
        vbox7.addWidget(self.button_crop)    
                 
        self.button_accept = QtGui.QPushButton('Apply Alignment')
        self.button_accept.clicked.connect(self.OnAccept)
        self.button_accept.setEnabled(False)
        vbox7.addWidget(self.button_accept)
         
        self.button_close = QtGui.QPushButton('Dismiss')
        self.button_close.clicked.connect(self.close)
        vbox7.addWidget(self.button_close)
         
         
         
               
        
        hboxtop = QtGui.QHBoxLayout()
        
        vboxL = QtGui.QVBoxLayout()
        vboxR = QtGui.QVBoxLayout()
        

        vboxL.addWidget(sizer8)
        vboxL.addStretch(1)
        vboxL.addWidget(groupBox9)
        vboxL.addStretch(1)
        vboxL.addWidget(sizer4)
        vboxL.addStretch(1)
        vboxL.addWidget(sizer6)
        vboxL.addStretch(1)
        vboxL.addLayout(vbox7)
        
        hboxRT = QtGui.QHBoxLayout()
        hboxRB = QtGui.QHBoxLayout()
        
        hboxRT.addLayout(vbox1)
        hboxRT.addLayout(vbox5)
         
        hboxRB.addLayout(vbox3)
        hboxRB.addLayout(vbox2)
        hboxRB.addStretch(1)
        
        vboxR.addStretch(0.2)
        vboxR.addLayout(hboxRT)
        vboxR.addStretch(1)
        vboxR.addLayout(hboxRB)
        vboxR.addStretch(1)
        
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
             

        _im = axes.imshow(np.rot90(image), cmap=matplotlib.cm.get_cmap('gray')) 
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
        
        
        self.tc_imageeng.setText('Image at energy: {0:5.2f} eV'.format(float(self.stack.ev[self.iev])))
        
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

        fileName = QtGui.QFileDialog.getSaveFileName(self, 'Save Reference Image', '', wildcard)

        fileName = str(fileName)
        if fileName == '':
            return
    
        
        img1 = Image.fromarray(self.ref_image)
        img1.save(fileName)
        
        
#----------------------------------------------------------------------        
    def LoadRefImage(self):
        
        wildcard = "TIFF File (*.tif);;"

        fileName = QtGui.QFileDialog.getOpenFileName(self, 'Save Reference Image', '', wildcard)

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
        
    
        _im = axes.imshow(np.rot90(self.ref_image), cmap=matplotlib.cm.get_cmap('gray')) 
        
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
    
        im = axes.imshow(np.rot90(ccorr), cmap=matplotlib.cm.get_cmap('gray')) 
        
        nx = ccorr.shape[0]
        ny = ccorr.shape[1]
        xcenter = xshift + np.float(nx)/2.0
        ycenter = yshift + np.float(ny)/2.0
        
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
        
        
        self.xshifts = np.delete(self.xshifts, self.iev) 
        self.yshifts = np.delete(self.yshifts, self.iev) 
        
        if self.com.stack_4d == 0:
            self.stack.data_struct.exchange.data = self.stack.absdata
        else:
            self.stack.data_struct.exchange.data = self.stack.stack4D
        self.stack.data_struct.exchange.energy = self.stack.ev
        
        if self.com.i0_loaded == 1:
            self.stack.calculate_optical_density()
            
        self.iev = self.iev-1
        if self.iev < 0:
            self.iev = 0

        
        self.slider_eng.setRange(0, self.stack.n_ev-1)

        self.parent.page1.slider_eng.setRange(0, self.stack.n_ev-1)
        self.parent.page1.iev = self.stack.n_ev/2
        self.parent.page1.slider_eng.setValue(self.parent.page1.iev)
        
        self.parent.page1.loadSpectrum(self.parent.page1.ix, self.parent.page1.iy)
        self.parent.page1.loadImage() 
        
        self.ShowImage() 

#----------------------------------------------------------------------            
    def OnCalcRegistration(self, event):
        
        if self.com.stack_4d == 1:
            self.CalcRegistration4D()
            return
        
        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        
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
        
                    
        QtGui.QApplication.restoreOverrideCursor()
        
#----------------------------------------------------------------------            
    def CalcRegistration4D(self):
        
        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        
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
        
                    
        QtGui.QApplication.restoreOverrideCursor()
        
#----------------------------------------------------------------------            
    def OnCropShifts(self, event):
        
        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        
        self.aligned_stack, self.xleft, self.xright, self.ybottom, self.ytop, = self.stack.crop_registed_images(self.aligned_stack, 
                                                             self.minxs, self.maxxs, self.minys, self.maxys)
        
        
        self.iev = 0
        self.ShowImage()
        self.slider_eng.setValue(self.iev)
        
        QtGui.QApplication.restoreOverrideCursor()
        
        
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
        else:
            self.stack.stack4D = self.aligned_stack  
            self.stack.absdata = self.stack.stack4D[:,:,:,self.itheta]
            
        
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
                        img = self.stack.od4D[:,:,i,j]
                        shifted_img = self.stack.apply_image_registration(img, self.xshifts[i,j], self.yshifts[i,j])         
                        self.stack.od4D[:,:,i,j] = shifted_img       
                        
                self.stack.od4D = self.stack.od4D[self.xleft:self.xright, self.ybottom:self.ytop, :, :]                  

                self.stack.od3d = self.stack.od4D[:,:,:,self.itheta]
                self.stack.od = self.stack.od3d.copy()
                n_pixels = self.stack.n_cols*self.stack.n_rows
                self.stack.od = np.reshape(self.stack.od, (n_pixels, self.stack.n_ev), order='F')      
                
            self.stack.data_struct.spectromicroscopy.optical_density = self.stack.od   

        self.stack.data_struct.exchange.data = self.stack.absdata
        self.stack.data_struct.exchange.energy = self.stack.ev
        

        
        self.stack.data_struct.spectromicroscopy.xshifts = self.xshifts
        self.stack.data_struct.spectromicroscopy.yshifts = self.yshifts
        
        
        self.parent.page1.slider_eng.setRange(0,self.stack.n_ev-1)
        self.parent.page1.iev = self.stack.n_ev/2
        self.parent.page1.slider_eng.setValue(self.parent.page1.iev)
        
        self.parent.page1.ix = int(self.stack.n_cols/2)
        self.parent.page1.iy = int(self.stack.n_rows/2)
        
        self.parent.page1.loadSpectrum(self.parent.page1.ix, self.parent.page1.iy)
        self.parent.page1.loadImage()
        
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

        self.SaveFileName = QtGui.QFileDialog.getSaveFileName(self, 'Save Image Shifts Plot', '', wildcard)

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

        filepath = QtGui.QFileDialog.getSaveFileName(self, 'Please select an alignment file (.csv)', '', wildcard)

        filepath = str(filepath)
        if filepath == '':
            return
        
        f = open(filepath, 'w')
        print>>f, '*********************  Alignment file  ********************'
        print>>f, '***  for ', self.com.filename
        print>>f, '***  ev, xshift, yshift'           
        for ie in range(self.stack.n_ev):
            print>>f, '%.6f, %.6f, %.6f' %(self.stack.ev[ie], self.xshifts[ie], self.yshifts[ie])
    
        f.close()            
        
#----------------------------------------------------------------------  
    def OnLoadShifts(self, evt):
        


        wildcard = "I0 CSV files (*.csv)"
        
        filepath = QtGui.QFileDialog.getOpenFileName(self, 'Please select an alignment file (.csv)', '', wildcard)
        

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
        
        self.UpdateWidgets()                    
           
            
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
               
        
                
        
#---------------------------------------------------------------------- 
class SpectralROI(QtGui.QDialog):

    def __init__(self, parent,  common, stack):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        self.stack = stack
        self.com = common  
        
        self.resize(630, 700)
        self.setWindowTitle('Spectral Regions of Interest')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal) 

        self.stack = stack
        self.com = common
               
        
        self.imin = 0
        self.imax = 0
        self.i0min = 0
        self.i0max = 0
        self.iselected = 0
        self.i0selected = 0        
        
                
        self.odtotal = self.stack.od3d.sum(axis=0)   
        self.odtotal = self.odtotal.sum(axis=0)/(self.stack.n_rows*self.stack.n_cols) 
        
        self.image_i0 = np.zeros((self.stack.n_cols, self.stack.n_rows))
        self.image_i = np.zeros((self.stack.n_cols, self.stack.n_rows))
        self.odthickmap = np.zeros((self.stack.n_cols, self.stack.n_rows))
        

        
        vbox = QtGui.QVBoxLayout()
        text = QtGui.QLabel(self)
        text.setText('First select I0 region below the edge, then select I region above the edge:')
        vbox.addWidget(text)
        
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
   
        self.specfig = Figure((6.0, 4.2))
        self.SpectrumPanel = FigureCanvas(self.specfig)
        self.SpectrumPanel.setParent(self)
        self.SpectrumPanel.mpl_connect('button_press_event', self.OnSelection1)
        self.SpectrumPanel.mpl_connect('button_release_event', self.OnSelection2)
               
        fbox.addWidget(self.SpectrumPanel)
        frame.setLayout(fbox)
        vbox.addWidget(frame)
       
       
        hbox2 = QtGui.QHBoxLayout()
               
        
        sizer2 =QtGui.QGroupBox('Selected Spectral Regions')
        sizer2.setMinimumWidth(350)
        vbox2 = QtGui.QVBoxLayout()
        text = QtGui.QLabel(self)
        text.setText('I Selection (red): ')
        self.textctrl1 = QtGui.QLabel(self)
        self.textctrl1.setText('[  ]' )
        vbox2.addWidget(text)
        vbox2.addWidget(self.textctrl1)
        
        text = QtGui.QLabel(self)
        text.setText('I0 Selection (green): ')
        self.textctrl2 = QtGui.QLabel(self)
        self.textctrl2.setText('[  ]' )
        vbox2.addWidget(text)
        vbox2.addWidget(self.textctrl2)
        sizer2.setLayout(vbox2)
        
        text = QtGui.QLabel(self)
        text.setText('Optical density map')
        vbox.addWidget(text)

        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
   
        self.odmfig = Figure((2.4,2.4))
        self.ODMImagePanel = FigureCanvas(self.odmfig)
        self.ODMImagePanel.setParent(self)
               
        fbox.addWidget(self.ODMImagePanel, stretch=0)
        frame.setLayout(fbox)
        hbox2.addWidget(frame, stretch=0)   
        hbox2.addStretch(0.5)
        hbox2.addWidget(sizer2)
        
        vbox.addLayout(hbox2)    
          
        hbox = QtGui.QHBoxLayout()
        
        button_save = QtGui.QPushButton('Save')
        button_save.clicked.connect(self.OnSave)
        hbox.addWidget(button_save)
        
        button_cancel = QtGui.QPushButton('Dismiss')
        button_cancel.clicked.connect(self.close)
        hbox.addWidget(button_cancel)
        
        vbox.addLayout(hbox)
    
        
        self.setLayout(vbox)
        
        self.draw_spectrum()
        
        
#----------------------------------------------------------------------        
    def draw_spectrum(self):

        
        fig = self.specfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        self.axes = fig.gca()
        

        specplot = self.axes.plot(self.stack.ev,self.odtotal)
        
        if self.i0selected == 1:
            self.axes.axvspan(self.stack.ev[self.i0min], self.stack.ev[self.i0max], facecolor='g', alpha=0.5)
            
        if self.iselected == 1:
            self.axes.axvspan(self.stack.ev[self.imin], self.stack.ev[self.imax], facecolor='r', alpha=0.5)

        
        self.axes.set_xlabel('Photon Energy [eV]')
        self.axes.set_ylabel('Optical Density')
        

        self.SpectrumPanel.draw()
        
    
#----------------------------------------------------------------------        
    def draw_image(self):

               
        fig = self.odmfig
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        
        axes = fig.gca()
        divider = make_axes_locatable(axes)
        axcb = divider.new_horizontal(size="3%", pad=0.03)  

        fig.add_axes(axcb)
        
        axes.set_position([0.03,0.03,0.8,0.94])
        
        
        im = axes.imshow(np.rot90(self.odthickmap), cmap=matplotlib.cm.get_cmap("gray")) 

        cbar = axes.figure.colorbar(im, orientation='vertical',cax=axcb) 

        #Show Scale Bar
        startx = int(self.stack.n_rows*0.05)
        starty = self.stack.n_cols-int(self.stack.n_cols*0.05)-self.stack.scale_bar_pixels_y
        um_string = ' $\mathrm{\mu m}$'
        microns = '$'+self.stack.scale_bar_string+' $'+um_string
        axes.text(self.stack.scale_bar_pixels_x+startx+1,starty+1, microns, horizontalalignment='left', verticalalignment='center',
                  color = 'white', fontsize=14)
        #Matplotlib has flipped scales so I'm using rows instead of cols!
        p = matplotlib.patches.Rectangle((startx,starty), self.stack.scale_bar_pixels_x, self.stack.scale_bar_pixels_y,
                               color = 'white', fill = True)
        axes.add_patch(p)

         
        axes.axis("off")  
        self.ODMImagePanel.draw()
        

#----------------------------------------------------------------------        
    def OnSelection1(self, evt):
        
        x1 = evt.xdata
                
        self.button_pressed = True
        self.patch = None
        self.conn = self.SpectrumPanel.mpl_connect('motion_notify_event', self.OnSelectionMotion) 
        
        if x1 == None:
            return
        
        self.x1 = x1

#----------------------------------------------------------------------        
    def OnSelection2(self, evt):
        
        x2 = evt.xdata
            

        self.button_pressed = False
        self.SpectrumPanel.mpl_disconnect(self.conn)    
        
        if x2 == None:
            return    
        
        x1 = self.x1

        
        if (self.i0selected == 1) and (self.iselected ==1):
            self.i0selected = 0
            self.iselected = 0
        
        if self.i0selected == 0:       
            self.i0min = np.abs(self.stack.ev - x1).argmin()
            self.i0max = np.abs(self.stack.ev - x2).argmin()
            
            self.image_i0 = np.sum(self.stack.absdata[:, :, self.i0min:self.i0max+1], axis=2)/(self.i0max+1-self.i0min)

            self.textctrl1.setText('Selection: [ '+str(self.stack.ev[self.i0min]) + ' eV, '+ str(self.stack.ev[self.i0max])+' eV ]' )
            self.i0selected = 1
            
        elif self.iselected == 0:
            self.imin = np.abs(self.stack.ev - x1).argmin()
            self.imax = np.abs(self.stack.ev - x2).argmin()
                       
            self.image_i = np.sum(self.stack.absdata[:, :, self.imin:self.imax+1], axis=2)/(self.imax+1-self.imin)
            
            self.textctrl2.setText('Selection: [ '+str(self.stack.ev[self.imin]) + ' eV, '+ str(self.stack.ev[self.imax])+' eV ]' )
            self.iselected = 1        
            
        if (self.i0selected == 1) and (self.iselected ==1):              
            nonzeroind = self.image_i0.nonzero()
            self.odthickmap = np.zeros((self.stack.n_cols, self.stack.n_rows))
            self.odthickmap[nonzeroind] = - np.log(self.image_i[nonzeroind]/self.image_i0[nonzeroind])
            self.draw_image()
            
        
        self.draw_spectrum()


#----------------------------------------------------------------------        
    def OnSelectionMotion(self, event):        

        x2 = event.xdata
        
        if x2 == None:
            return  
        
        x1 = self.x1
        
                
        fig = self.specfig

        axes = fig.gca()

        if self.patch != None:
            self.patch.remove()
        self.patch = self.axes.axvspan(x1, x2, facecolor='w', alpha=0.5)

        self.SpectrumPanel.draw()    
        
        
#----------------------------------------------------------------------        
    def OnSave(self, evt):
        #Save images

        wildcard = "Portable Network Graphics (*.png);;Adobe PDF Files (*.pdf);;TIFF File (*.tif);;"

        fileName = QtGui.QFileDialog.getSaveFileName(self, 'Save OD Map', '', wildcard)

        fileName = str(fileName)
        if fileName == '':
            return
        
        path, ext = os.path.splitext(fileName) 
        ext = ext[1:].lower() 
                               

       
        if ext != 'png' and ext != 'pdf' and ext != 'tif': 
            error_message = ( 
                  'Only the PNG and PDF image formats are supported.\n' 
                 'A file extension of `png\' or `pdf\' must be used.') 
            QtGui.QMessageBox.warning(self, 'Error', 'Error - Could not save file.') 
            return 
        
        
        if ext == 'tif':
            from PIL import Image  
            img1 = Image.fromarray(self.odthickmap)
            img1.save(fileName)
        else:
            
            try: 
                matplotlib.rcParams['pdf.fonttype'] = 42
                
                fig = self.odmfig
                fig.savefig(fileName)
    
                
            except IOError, e:
                if e.strerror:
                    err = e.strerror 
                else: 
                    err = e 
       
                QtGui.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err) 
            



#---------------------------------------------------------------------- 
class DoseCalculation(QtGui.QDialog):

    def __init__(self, parent,  stack, ROIspectrum):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent
        
        self.resize(300, 170)
        self.setWindowTitle('Dose Calculation')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal) 
                
        self.stack = stack
        self.ROIspectrum = ROIspectrum
               

        
        vboxtop = QtGui.QVBoxLayout()
        
        
        gridtop = QtGui.QGridLayout()

        
        #fontb = wx.SystemSettings_GetFont(wx.SYS_DEFAULT_GUI_FONT)
        #fontb.SetWeight(wx.BOLD)
        
        
        st1 = QtGui.QLabel(self)
        st1.setText('Detector efficiency [%]:')
        #st1.SetFont(fontb)
        st2 = QtGui.QLabel(self)
        st2.setText('I region composition:')
        #st2.SetFont(fontb)
#        st3 = QtGui.QLabel(self,'Xray absorption length:')
#        st3.SetFont(fontb)
        st4 = QtGui.QLabel(self)
        st4.setText('Dose [Gray]:')
        #st4.SetFont(fontb)

        
        self.tc_1 = QtGui.QLineEdit(self)
        self.tc_1.setText('30')
        
        self.tc_2 = QtGui.QLineEdit(self)
        
#        self.tc_3 = wx.TextCtrl(panel1, -1, size=((200,-1)), style=wx.TE_RICH|wx.VSCROLL|wx.TE_READONLY, 
#                                         value=' ')   
        
        self.tc_4 = QtGui.QLabel(self)
        

        gridtop.addWidget(st1, 0,0)
        gridtop.addWidget( self.tc_1, 0,1)
        gridtop.addWidget(st2, 1,0)
        gridtop.addWidget( self.tc_2, 1,1)
#        gridtop.addWidget(st3, 0)
#        gridtop.addWidget( self.tc_3, 0)
        gridtop.addWidget(st4, 2,0)        
        gridtop.addWidget( self.tc_4, 2,1)             
          


        button_calcdose = QtGui.QPushButton('Calculate Dose')
        button_calcdose.clicked.connect(self.OnCalcDose)
                
        button_cancel = QtGui.QPushButton('Dismiss')
        button_cancel.clicked.connect(self.close)


        vboxtop.addLayout(gridtop)
        vboxtop.addWidget(button_calcdose) 
        vboxtop.addWidget(button_cancel) 
        
        self.setLayout(vboxtop)
        
              
#---------------------------------------------------------------------- 
    def CalcDose(self):
                      
        
        try:
            detector_eff = 0.01*float(self.tc_1.text())
        except:
            QtGui.QMessageBox.warning(self, 'Error', 'Please enter numeric number for detector efficiency.')
            print 'Please enter numeric number for detector efficiency.'
            return
            
        
        i_composition = str(self.tc_2.text())
        
        dose = 0.
        
        Chenke = henke.henke()
        
        #Check if composition array is recognizable
        try:
            z_array, atwt = Chenke.compound(i_composition,1.0)
        except:
            QtGui.QMessageBox.warning(self, 'Error', "Please enter new compound.")
            return
        
        try:
            dose = Chenke.dose_calc(self.stack, i_composition, self.ROIspectrum, self.stack.i0data, detector_eff)
        except:
            QtGui.QMessageBox.warning(self, 'Error', "Could not calculate dose. Please enter new compound.")
            return            
        
        self.tc_4.setText(str(dose))
        
        return
        
        
#---------------------------------------------------------------------- 
    def OnCalcDose(self, evt):

        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        self.CalcDose()
        QtGui.QApplication.restoreOverrideCursor()
        
        
#---------------------------------------------------------------------- 
class DarkSignal(QtGui.QDialog):

    def __init__(self, parent, common, stack):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        self.stack = stack
        self.com = common
        
        self.resize(300, 170)
        self.setWindowTitle('Dark Signal Correction')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal) 
                
        
        vboxtop = QtGui.QVBoxLayout()
        
        
        gridtop = QtGui.QGridLayout()

        
        #fontb = wx.SystemSettings_GetFont(wx.SYS_DEFAULT_GUI_FONT)
        #fontb.SetWeight(wx.BOLD)
        
        
        st1 = QtGui.QLabel(self)
        st1.setText('Dark Signal Value:')

        
        self.ntc_ds = QtGui.QLineEdit(self)
        self.ntc_ds.setFixedWidth(150)
        self.ntc_ds.setValidator(QtGui.QDoubleValidator(-99999, 99999, 2, self))
        self.ntc_ds.setAlignment(QtCore.Qt.AlignRight)         
        
        self.ntc_ds.setText(str(0.0))
 

        gridtop.addWidget(st1, 0,0)
        gridtop.addWidget(self.ntc_ds, 0,1)
        

        button_dscalc = QtGui.QPushButton('Subtract Dark Signal')
        button_dscalc.clicked.connect(self.OnDSCalc)
                
        button_cancel = QtGui.QPushButton('Dismiss')
        button_cancel.clicked.connect(self.close)


        vboxtop.addLayout(gridtop)
        vboxtop.addWidget(button_dscalc) 
        vboxtop.addWidget(button_cancel) 
        
        self.setLayout(vboxtop)
        
              
#---------------------------------------------------------------------- 
    def OnDSCalc(self):
                      
        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        
        darksig = 0.0
        
        try:
            value = self.ntc_ds.text()
            darksig = float(value)
        except:
            QtGui.QApplication.restoreOverrideCursor()
            QtGui.QMessageBox.warning(self, 'Error', 'Please enter numeric number for dark signal.')
            return
            

        self.stack.absdata = self.stack.absdata - darksig
        
        if self.com.i0_loaded == 1:   
            self.stack.calculate_optical_density()
            
        self.stack.fill_h5_struct_from_stk()
        if self.com.i0_loaded == 1: 
            self.stack.fill_h5_struct_normalization()
        
        
        self.parent.page1.loadSpectrum(self.parent.page1.ix, self.parent.page1.iy)
        self.parent.page1.loadImage()
        
        
        
        QtGui.QApplication.restoreOverrideCursor()
        
        self.close()
        
        return
        
        
        
                
#---------------------------------------------------------------------- 
class PlotFrame(QtGui.QDialog):

    def __init__(self, parent, datax, datay, title = "I0 data"):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent
        
        self.title = title

        self.datax = datax
        self.datay = datay
        
        self.resize(630, 500)
        self.setWindowTitle(title)
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal) 
                
        

        vbox = QtGui.QVBoxLayout()
        
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
   
        self.plotfig = Figure((6.0, 4.2))
        self.PlotPanel = FigureCanvas(self.plotfig)
        self.PlotPanel.setParent(self)
        
        fbox.addWidget(self.PlotPanel)
        frame.setLayout(fbox)
        vbox.addWidget(frame)
        

        hbox = QtGui.QHBoxLayout()
        
        button_save = QtGui.QPushButton('Save Spectrum')
        button_save.clicked.connect(self.OnSave)
        hbox.addWidget(button_save)        
                
        button_close = QtGui.QPushButton('Close')
        button_close.clicked.connect(self.close)
        hbox.addWidget(button_close)
        
        vbox.addLayout(hbox)

        self.setLayout(vbox)

        
        self.draw_plot(datax,datay)
        
      
#----------------------------------------------------------------------        
    def draw_plot(self, datax, datay):
        
        
        fig = self.plotfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        self.axes = fig.gca()
        
        
        plot = self.axes.plot(datax,datay)
        
        self.axes.set_xlabel('Photon Energy [eV]')
        self.axes.set_ylabel('I0 Flux')
        

        
        self.PlotPanel.draw()
        
#----------------------------------------------------------------------
    def OnSave(self, event):


        try: 
            wildcard = "CSV files (*.csv)"

            filepath = QtGui.QFileDialog.getSaveFileName(self, 'Save plot as .csv (.csv)', '', wildcard)

            filepath = str(filepath)
            if filepath == '':
                return
        
                            
            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))                 
            self.Save(filepath)    
            QtGui.QApplication.restoreOverrideCursor()    

        except:
 
            QtGui.QApplication.restoreOverrideCursor()
            QtGui.QMessageBox.warning(self, 'Error', "Could not save .csv file.")
                   
        self.close()

        
        return

#----------------------------------------------------------------------
    def Save(self, filename):
            
        f = open(filename, 'w')
        print>>f, '*********************  X-ray Absorption Data  ********************'
        print>>f, '*'
        print>>f, '* Formula: '
        print>>f, '* Common name: ', self.title
        print>>f, '* Edge: '
        print>>f, '* Acquisition mode: '
        print>>f, '* Source and purity: ' 
        print>>f, '* Comments: Stack list ROI ""'
        print>>f, '* Delta eV: '
        print>>f, '* Min eV: '
        print>>f, '* Max eV: '
        print>>f, '* Y axis: '
        print>>f, '* Contact person: '
        print>>f, '* Write date: '
        print>>f, '* Journal: '
        print>>f, '* Authors: '
        print>>f, '* Title: '
        print>>f, '* Volume: '
        print>>f, '* Issue number: '
        print>>f, '* Year: '
        print>>f, '* Pages: '
        print>>f, '* Booktitle: '
        print>>f, '* Editors: '
        print>>f, '* Publisher: '
        print>>f, '* Address: '
        print>>f, '*--------------------------------------------------------------'
        dim = self.datax.shape
        n=dim[0]
        for ie in range(n):
            print>>f, '{0:06.2f}, {1:06f}'.format(self.datax[ie], self.datay[ie])
        
        f.close()
    
        

#---------------------------------------------------------------------- 
class ColorTableFrame(QtGui.QDialog):

    def __init__(self, parent):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        
        self.resize(200, 430)
        self.setWindowTitle('Pick Color Table')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)    
        
        self.colors= ["gray","jet","autumn","bone", "cool","copper", "flag","hot","hsv","pink",
                      "prism","spring","summer","winter", "spectral"]
        
        vboxtop = QtGui.QVBoxLayout()
        
        
               
        sizer1 = QtGui.QGroupBox('Color Tables')
        vbox = QtGui.QVBoxLayout()
        
        self.rb_grey = QtGui.QRadioButton(self.colors[0], self)
        self.rb_jet  = QtGui.QRadioButton(self.colors[1], self)
        self.rb_autumn = QtGui.QRadioButton(self.colors[2], self)
        self.rb_bone = QtGui.QRadioButton(self.colors[3], self)
        self.rb_cool = QtGui.QRadioButton(self.colors[4], self)
        self.rb_copper = QtGui.QRadioButton(self.colors[5], self)
        self.rb_flag = QtGui.QRadioButton(self.colors[6], self)
        self.rb_hot = QtGui.QRadioButton(self.colors[7], self)
        self.rb_hsv = QtGui.QRadioButton(self.colors[8], self)
        self.rb_pink = QtGui.QRadioButton(self.colors[9], self)
        self.rb_prism = QtGui.QRadioButton(self.colors[10], self)
        self.rb_spring = QtGui.QRadioButton(self.colors[11], self)
        self.rb_summer = QtGui.QRadioButton(self.colors[12], self)
        self.rb_winter = QtGui.QRadioButton(self.colors[13], self)
        self.rb_spectral = QtGui.QRadioButton(self.colors[14], self)
        
        
        self.radios = [] 
        
        self.radios.append(self.rb_grey)
        self.radios.append(self.rb_jet)
        self.radios.append(self.rb_autumn)
        self.radios.append(self.rb_bone)
        self.radios.append(self.rb_cool)
        self.radios.append(self.rb_copper)
        self.radios.append(self.rb_flag)
        self.radios.append(self.rb_hot)
        self.radios.append(self.rb_hsv)
        self.radios.append(self.rb_pink)
        self.radios.append(self.rb_prism)
        self.radios.append(self.rb_spring)
        self.radios.append(self.rb_summer)
        self.radios.append(self.rb_winter)
        self.radios.append(self.rb_spectral)
        
        self.ct_dict = dict([(self.colors[x], self.radios[x]) for x in range(len(self.colors))])
        

        self.rb_grey.clicked.connect(self.OnColorTable)
        self.rb_jet.clicked.connect(self.OnColorTable)
        self.rb_autumn.clicked.connect(self.OnColorTable)
        self.rb_bone.clicked.connect(self.OnColorTable)
        self.rb_cool.clicked.connect(self.OnColorTable)
        self.rb_copper.clicked.connect(self.OnColorTable)
        self.rb_flag.clicked.connect(self.OnColorTable)
        self.rb_hot.clicked.connect(self.OnColorTable)
        self.rb_hsv.clicked.connect(self.OnColorTable)
        self.rb_pink.clicked.connect(self.OnColorTable)
        self.rb_prism.clicked.connect(self.OnColorTable)
        self.rb_spring.clicked.connect(self.OnColorTable)
        self.rb_summer.clicked.connect(self.OnColorTable)
        self.rb_winter.clicked.connect(self.OnColorTable)
        self.rb_spectral.clicked.connect(self.OnColorTable)
        
        vbox.addWidget(self.rb_grey)
        vbox.addWidget(self.rb_jet)
        vbox.addWidget(self.rb_autumn)
        vbox.addWidget(self.rb_bone)
        vbox.addWidget(self.rb_cool)
        vbox.addWidget(self.rb_copper)
        vbox.addWidget(self.rb_flag)
        vbox.addWidget(self.rb_hot)
        vbox.addWidget(self.rb_hsv)
        vbox.addWidget(self.rb_pink)
        vbox.addWidget(self.rb_prism)
        vbox.addWidget(self.rb_spring)
        vbox.addWidget(self.rb_summer)
        vbox.addWidget(self.rb_winter)
        vbox.addWidget(self.rb_spectral)

        sizer1.setLayout(vbox)
        vboxtop.addWidget(sizer1)

        
        button_close = QtGui.QPushButton('Close')
        button_close.clicked.connect(self.close)
        vboxtop.addWidget(button_close)
        
        self.setLayout(vboxtop)
        
        self.ct_dict[self.parent.page1.colortable].setChecked(True)
        
        
#----------------------------------------------------------------------          
    def OnColorTable(self, event):
    
        for radioButton in self.findChildren(QtGui.QRadioButton):
            if radioButton.isChecked():
                radioButtonText = radioButton.text()
                break
            

        self.parent.page1.colortable = str(radioButtonText)
        
        self.parent.page1.loadImage()
                
                                                        

""" ------------------------------------------------------------------------------------------------"""
class PageLoadData(QtGui.QWidget):
    def __init__(self, common, data_struct, stack):
        super(PageLoadData, self).__init__()

        self.initUI(common, data_struct, stack)
        
#----------------------------------------------------------------------          
    def initUI(self, common, data_struct, stack): 
        
        self.data_struct = data_struct
        self.stk = stack
        self.com = common                  
        
        self.filename = " "
        
        self.iev = 0
        self.itheta = 0
        
        self.initMatplotlib()

        #panel 1
        sizer1 = QtGui.QGroupBox('Load Data Stack')
        vbox1 = QtGui.QVBoxLayout()
        
        self.button_multiload = QtGui.QPushButton('  Load Stack  ')
        self.button_multiload.clicked.connect( self.OnLoadMulti)
        vbox1.addWidget(self.button_multiload)
        
        
        self.button_4d = QtGui.QPushButton( 'Load 4D stack TOMO-XANES')
        self.button_4d.clicked.connect( self.OnLoad4D)
        vbox1.addWidget(self.button_4d) 

        sizer1.setLayout(vbox1)

        #panel 2
        sizer2 = QtGui.QGroupBox('Build a stack from a set of files')
        vbox2 = QtGui.QVBoxLayout()

        button_sm = QtGui.QPushButton( '  Select a directory with stack files  ')
        button_sm.clicked.connect( self.OnBuildStack)
        vbox2.addWidget(button_sm)
        
        sizer2.setLayout(vbox2)


        #panel 3
        sizer3 = QtGui.QGroupBox('File')
        vbox3 = QtGui.QVBoxLayout()
 
  
        self.tc_file = QtGui.QLabel(self)
        vbox3.addWidget(self.tc_file)
        self.tc_file.setText('File name')
        
        vbox3.setContentsMargins(20,20,20,30)
        sizer3.setLayout(vbox3)
        

        #panel 4
        sizer4 = QtGui.QGroupBox('Path')
        vbox4 = QtGui.QVBoxLayout()
  
        self.tc_path = QtGui.QLabel(self)
        vbox4.addWidget(self.tc_path)
        self.tc_path.setText('D:/')
       
        vbox4.setContentsMargins(20,20,20,30)
        sizer4.setLayout(vbox4)
                
 
        #panel 5     
        vbox5 = QtGui.QVBoxLayout()
        
        self.tc_imageeng = QtGui.QLabel(self)
        self.tc_imageeng.setText("Image at energy: ")
        vbox5.addWidget(self.tc_imageeng)
        
        

        gridsizertop = QtGui.QGridLayout()
        
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
   
        self.absimgfig = Figure((PlotH, PlotH))

        self.AbsImagePanel = FigureCanvas(self.absimgfig)
        self.AbsImagePanel.setParent(self)
        
        
        fbox.addWidget(self.AbsImagePanel)
        frame.setLayout(fbox)
        gridsizertop.addWidget(frame, 0, 0, QtCore .Qt. AlignLeft)
        

        self.slider_eng = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slider_eng.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_eng.valueChanged[int].connect(self.OnScrollEng)
        self.slider_eng.setRange(0, 100)
        
        gridsizertop.addWidget(self.slider_eng, 0, 1, QtCore .Qt. AlignLeft)
        
        
        self.slider_theta = QtGui.QScrollBar(QtCore.Qt.Horizontal)
        self.slider_theta.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_theta.valueChanged[int].connect(self.OnScrollTheta)
        self.slider_theta.setRange(0, 100)     
        self.slider_theta.setVisible(False)  
        self.tc_imagetheta = QtGui.QLabel(self)
        self.tc_imagetheta.setText("4D Data Angle: ")
        self.tc_imagetheta.setVisible(False)
        hbox51 = QtGui.QHBoxLayout()
        hbox51.addWidget(self.tc_imagetheta) 
        hbox51.addWidget(self.slider_theta)
        gridsizertop.addLayout(hbox51, 1, 0)
        

        vbox5.addLayout(gridsizertop)
        vbox5.addStretch(1)
        
       
        vboxtop = QtGui.QVBoxLayout()
        
        hboxtop = QtGui.QHBoxLayout()
        vboxt1 = QtGui.QVBoxLayout()
        vboxt1.addStretch (1)
        vboxt1.addWidget(sizer1)
        vboxt1.addStretch (1)
        vboxt1.addWidget(sizer2)
        vboxt1.addStretch (1)
        
        hboxtop.addStretch (0.5)
        hboxtop.addLayout(vboxt1)
        hboxtop.addStretch (0.5)
        hboxtop.addLayout(vbox5)
        hboxtop.addStretch (0.5)
        
        vboxtop.addStretch (0.5)
        vboxtop.addLayout(hboxtop)
        vboxtop.addStretch (0.5)
        vboxtop.addWidget(sizer3)
        vboxtop.addStretch (0.5)
        vboxtop.addWidget(sizer4)
        vboxtop.addStretch (0.5)

        vboxtop.setContentsMargins(50,50,50,50)
        self.setLayout(vboxtop)


    
#----------------------------------------------------------------------   
    def initMatplotlib(self):  
               
       
        matplotlib.rcParams['figure.facecolor'] = 'white'

        matplotlib.rcParams['font.size'] = 10.0
        
        
#----------------------------------------------------------------------          
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
            self.ShowImage()
            
#----------------------------------------------------------------------            
    def OnScrollTheta(self, value):
        self.itheta = value

        self.stk.absdata = self.stk.stack4D[:,:,:,self.itheta]
        self.tc_imagetheta.setText("4D Data Angle: "+str(self.stk.theta[self.itheta]))
        
        if self.com.stack_loaded == 1:
            self.ShowImage()
            
        self.window().page1.itheta = self.itheta
        self.window().page1.slider_theta.setValue(self.itheta)
        
        self.window().page2.itheta = self.itheta
        self.window().page2.slider_theta.setValue(self.itheta)
            
                    
#----------------------------------------------------------------------        
    def ShowImage(self):
        
        
        image = self.stk.absdata[:,:,int(self.iev)].copy() 

        fig = self.absimgfig
        fig.clf()
        fig.add_axes(((0.0,0.0,1.0,1.0)))
        axes = fig.gca()
        fig.patch.set_alpha(1.0)
         
        im = axes.imshow(np.rot90( image ), cmap=matplotlib.cm.get_cmap("gray")) 
         
        if self.window().page1.show_scale_bar == 1:
            #Show Scale Bar
            if self.com.white_scale_bar == 1:
                sbcolor = 'white'
            else:
                sbcolor = 'black'
            startx = int(self.stk.n_cols*0.05)
            starty = self.stk.n_rows-int(self.stk.n_rows*0.05)-self.stk.scale_bar_pixels_y
            um_string = ' $\mathrm{\mu m}$'
            microns = '$'+self.stk.scale_bar_string+' $'+um_string
            axes.text(self.stk.scale_bar_pixels_x+startx+1,starty+1, microns, horizontalalignment='left', verticalalignment='center',
                      color = sbcolor, fontsize=14)
            #Matplotlib has flipped scales so I'm using rows instead of cols!
            p = matplotlib.patches.Rectangle((startx,starty), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                                   color = sbcolor, fill = True)
            axes.add_patch(p)
             
        
        axes.axis("off")      
        self.AbsImagePanel.draw()
         
        self.tc_imageeng.setText('Image at energy: {0:5.2f} eV'.format(float(self.stk.ev[self.iev])))


        
#----------------------------------------------------------------------        
    def ShowInfo(self, filename, filepath):
        
        self.ShowImage()
        
        self.tc_file.setText(filename)
        self.tc_path.setText(filepath)
        
#---------------------------------------------------------------------- 
class StackListFrame(QtGui.QDialog):

    def __init__(self, parent, filepath, com, stack, data_struct):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        self.data_struct = data_struct
        self.stk = stack
        self.common = com
        
        self.resize(600, 500)
        self.setWindowTitle('Stack File List')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)
        
        self.filepath = str(filepath)
        
        self.have1st = 0
        self.havelast = 0
        self.file1st = ' '
        self.filelast = ' '
        
        self.filetype = ''
        
        vbox = QtGui.QVBoxLayout()
                
        self.textt = QtGui.QLabel(self)
        self.textt.setText('Select first stack file')    
      
        vbox.addStretch(1)
        vbox.addWidget(self.textt)  
        

        self.filelist = QtGui.QTableWidget()
        self.filelist.setMinimumHeight(450)
        self.filelist.setColumnCount(4)
        self.filelist.setHorizontalHeaderLabels(('File list', 'X', 'Y', 'eV'))
        self.filelist.setShowGrid(False)
        self.filelist.verticalHeader().setVisible(False)
        
        self.filelist.setColumnWidth(0,400)
        self.filelist.setColumnWidth(1,50)
        self.filelist.setColumnWidth(2,50)
        self.filelist.setColumnWidth(3,50)
        
        self.filelist.setRowCount(0)
        
        self.filelist.cellClicked.connect(self.OnFileList)



        vbox.addWidget(self.filelist)
        vbox.addStretch(1) 
        
        self.tc_first = QtGui.QLabel(self)
        self.tc_first.setText('First stack file: ')
        self.tc_last = QtGui.QLabel(self)
        self.tc_last.setText('Last stack file: ')

        
        vbox.addWidget(self.tc_first)
        vbox.addWidget(self.tc_last)
        vbox.addStretch(1)
        
        
        hbox = QtGui.QHBoxLayout()
        
        
        self.button_accept = QtGui.QPushButton('Accept')
        #self.button_accept.setEnabled(False)
        self.button_accept.clicked.connect( self.OnAccept)
        hbox.addWidget(self.button_accept)
        
        button_cancel = QtGui.QPushButton('Cancel')
        button_cancel.clicked.connect( self.close)
        hbox.addWidget(button_cancel)
        
        vbox.addLayout(hbox)
        
        vbox.addStretch(0.5)
                        
        
        self.setLayout(vbox)
        
        self.ShowFileList()

#----------------------------------------------------------------------          
    def OnFileList(self, row, column):
    
    
        if (self.have1st == 1) and (self.havelast==1):
            self.have1st = 0
            self.havelast = 0
            #self.button_accept.setEnabled(False)
            self.tc_first.setText('First stack file: ')
            self.tc_last.setText('Last stack file: ')
            
        item = self.filelist.item(row, 0)
        fn =  item.text()
        if self.have1st == 0:
            self.tc_first.setText('First stack file: ' + fn) 
            self.textt.setText('Select last stack file')
            self.file1st = fn
            self.have1st = 1
        elif self.havelast == 0:
            self.tc_last.setText('Last stack file: ' + fn) 
            self.textt.setText('Select first stack file')         
            self.filelast = fn
            #self.button_accept.setEnabled(True)
            self.havelast = 1
                    
#----------------------------------------------------------------------           
    def ShowFileList(self):    
        
        filepath = str(self.filepath)
        self.sm_files = [x for x in os.listdir(filepath) if x.endswith('.sm')]
        
        if self.sm_files:
            
            self.filetype = 'sm'
            
            try: 
                from netCDF4 import Dataset
                from file_plugins import file_sm_netcdf
            
            except:
                QtGui.QMessageBox.warning(self, 'Error', "Could not import netCDF4 library.")
                return

            count = 0
        
            for i in range(len(self.sm_files)):

                filename = self.sm_files[i]
                thisfile = os.path.join(filepath, filename)
            
                filever, ncols, nrows, iev = file_sm_netcdf.read_sm_header(thisfile)            
                
                if filever > 0:  
                    self.filelist.insertRow(count)
                    self.filelist.setRowHeight(count,20)

                    self.filelist.setItem(count, 0, QtGui.QTableWidgetItem(filename))
                    self.filelist.setItem(count, 1, QtGui.QTableWidgetItem(str(ncols)))
                    self.filelist.setItem(count, 2, QtGui.QTableWidgetItem(str(nrows)))
                    self.filelist.setItem(count, 3, QtGui.QTableWidgetItem('{0:5.2f}'.format(iev)))
                                     
                    count += 1
                else:
                    continue
                
            
         
            
        self.xrm_files = [x for x in os.listdir(filepath) if x.endswith('.xrm')] 
        

        if self.xrm_files:        
            
            self.filetype = 'xrm'

            count = 0
        
            for i in range(len(self.xrm_files)):

                filename = self.xrm_files[i]
                thisfile = os.path.join(filepath, filename)

                ncols, nrows, iev = file_xrm.read_xrm_fileinfo(thisfile)
  
                if ncols > 0:                       
                    self.filelist.insertRow(count)
                    self.filelist.setRowHeight(count,20)

                    self.filelist.setItem(count, 0, QtGui.QTableWidgetItem(filename))
                    self.filelist.setItem(count, 1, QtGui.QTableWidgetItem(str(ncols)))
                    self.filelist.setItem(count, 2, QtGui.QTableWidgetItem(str(nrows)))
                    self.filelist.setItem(count, 3, QtGui.QTableWidgetItem('{0:5.2f}'.format(iev)))
                            
                    count += 1

            
            self.sm_files = self.xrm_files
        
        
        self.bim_files = [x for x in os.listdir(filepath) if x.endswith('.bim')]
        
        
        if self.bim_files:
            
            self.filetype = 'bim'
            
            

            count = 0
        
            for i in range(len(self.bim_files)):
                #print sm_files
                filename = self.bim_files[i]
                thisfile = os.path.join(filepath, filename)
            
                ncols, nrows, iev = file_bim.read_bim_info(thisfile)

                if ncols >0 :
                    self.filelist.insertRow(count)
                    self.filelist.setRowHeight(count,20)
    
                    self.filelist.setItem(count, 0, QtGui.QTableWidgetItem(filename))
                    self.filelist.setItem(count, 1, QtGui.QTableWidgetItem(str(ncols)))
                    self.filelist.setItem(count, 2, QtGui.QTableWidgetItem(str(nrows)))
                    self.filelist.setItem(count, 3, QtGui.QTableWidgetItem('{0:5.2f}'.format(iev)))
                                     
                    count += 1

                
            self.sm_files = self.bim_files
        
        return
        

        
#----------------------------------------------------------------------        
    def OnAccept(self, evt):
        
        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        
        self.parent.new_stack_refresh()
        self.stk.new_data()
        
        ind1st = self.sm_files.index(self.file1st) 
        indlast = self.sm_files.index(self.filelast)
        
        filelist = self.sm_files[ind1st:indlast+1]
        

        if self.filetype == 'sm':
            from file_plugins import file_sm_netcdf
            file_sm_netcdf.read_sm_list(self, filelist, self.filepath, self.data_struct)
        elif self.filetype == 'xrm':
            file_xrm.read_xrm_list(self, filelist, self.filepath, self.data_struct)
        elif self.filetype == 'bim':
            file_bim.read_bim_list(self, filelist, self.filepath, self.data_struct)
        else:
            print 'Wrong file type'
            return
        
        
        #fill the gui structure data
        self.stk.absdata = self.data_struct.exchange.data
        
        datadim = np.int32(self.stk.absdata.shape)
        self.stk.n_cols = datadim[0].copy()
        self.stk.n_rows =  datadim[1].copy()
        self.stk.ev = self.data_struct.exchange.energy 
        self.stk.n_ev = np.int32(self.stk.ev.shape[0]).copy()
        
        npixels = self.stk.n_cols*self.stk.n_rows*self.stk.n_ev
        
        
        self.stk.x_dist = self.data_struct.exchange.x 
        self.stk.y_dist = self.data_struct.exchange.y
        
        self.stk.data_dwell = self.data_struct.spectromicroscopy.data_dwell
        
        
        self.stk.fill_h5_struct_from_stk()
        
        self.stk.scale_bar()
                   

        self.parent.page1.iev = int(self.stk.n_ev/3)
   
        self.parent.ix = int(self.stk.n_cols/2)
        self.parent.iy = int(self.stk.n_rows/2)
                        
        self.common.stack_loaded = 1
        
        self.parent.refresh_widgets()
        self.parent.page1.ResetDisplaySettings()
        self.parent.page1.filename = filelist[0]
        self.parent.page1.textctrl.setText(filelist[0])
 
        self.parent.page0.slider_eng.setRange(0,self.stk.n_ev-1)
        self.parent.page0.iev = self.stk.n_ev/2
        self.parent.page0.slider_eng.setValue(self.parent.page1.iev)       
        
        self.parent.page1.slider_eng.setRange(0,self.stk.n_ev-1)
        self.parent.page1.iev = self.stk.n_ev/2
        self.parent.page1.slider_eng.setValue(self.parent.page1.iev)
        
        self.parent.page1.loadSpectrum(self.parent.page1.ix, self.parent.page1.iy)
        self.parent.page1.loadImage()
        
        self.parent.page0.ShowInfo(filelist[0], self.filepath)
        
        QtGui.QApplication.restoreOverrideCursor()
        self.close()
        
                
#----------------------------------------------------------------------  
class InputRegionDialog(QtGui.QDialog):

    def __init__(self, parent, nregions, title='Multi Region Stack'):
    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent
    
        mainLayout = QtGui.QVBoxLayout()
    
        layout = QtGui.QHBoxLayout()
        label = QtGui.QLabel()
        label.setText('Select Region to load:')
        layout.addWidget(label)
        combo = QtGui.QComboBox(self)
        combo.addItem("All")
        for i in range(nregions):
            combo.addItem(str(i+1))
        combo.setCurrentIndex(1)
        self.parent.loadregion = 1
        combo.activated[int].connect(self.OnSelectRegion) 

        layout.addWidget(combo)
    
        mainLayout.addLayout(layout)
    

        layout = QtGui.QHBoxLayout()
        button = QtGui.QPushButton("Submit") #string or icon
        self.connect(button, QtCore.SIGNAL("clicked()"), self.close)
        layout.addWidget(button)
    
        mainLayout.addLayout(layout)
        self.setLayout(mainLayout)
    
        self.resize(250, 60)
        self.setWindowTitle(title)
                               
#----------------------------------------------------------------------           
    def OnSelectRegion(self, value):
        item = value
        selregion = item
        self.parent.loadregion = selregion
        
        
        
#---------------------------------------------------------------------- 
class AboutFrame(QtGui.QDialog):

    def __init__(self, parent = None, title='About'):    
        QtGui.QWidget.__init__(self, parent)

        self.resize(360, 660)
        self.setWindowTitle('About Mantis')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)
        
        vbox = QtGui.QVBoxLayout()
        
       

        self.image = QtGui.QImage(resource_path(os.path.join('images','Mantis_logo_about.png')))
        
        self.imageLabel = QtGui.QLabel()
        self.imageLabel.setBackgroundRole(QtGui.QPalette.Base)
        
        self.imageLabel.setPixmap(QtGui.QPixmap.fromImage(self.image))
        vbox.addWidget(self.imageLabel)
        


        text1 = QtGui.QLabel(self)
        text1.setText("www.2ndlookconsulting.com")
        text1.setStyleSheet('color: rgb(53,159,217);font-size: 14pt; font-family: SansSerif;')
        #text1.setFont(QtGui.QFont('SansSerif', 14))
        
        #font2 = wx.Font(12, wx.SWISS, wx.NORMAL, wx.NORMAL)
        text2 = QtGui.QLabel(self)
        text2.setText('Mantis '+version)  
        text2.setStyleSheet('font-size: 12pt')
        #text2.SetFont(font2)        

        #font3 = wx.Font(8, wx.SWISS, wx.NORMAL, wx.NORMAL)
        text3 = QtGui.QLabel(self)
        text3.setText( '''
Developed by Mirna Lerotic, based on earlier programs by Mirna 
Lerotic and Chris Jacobsen. Initial development supported by 
Argonne National Laboratory LDRD 2010-193-R1 9113. ''')  
        #text3.SetFont(font3)          

        
        #font4 = wx.Font(8, wx.SWISS, wx.NORMAL, wx.NORMAL)
        text4 = QtGui.QLabel(self)
        text4.setText( '''        
Mantis is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published 
by the Free Software Foundation, either version 3 of the License, 
or any later version.

Mantis is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty 
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU General Public License for more details 
http://www.gnu.org/licenses/.''')  
        #text4.SetFont(font4)          

        vbox.addStretch(1)
        hbox = QtGui.QHBoxLayout()
        hbox.addStretch(1)
        vbox2 = QtGui.QVBoxLayout()
        vbox2.addWidget(text1)  
        vbox2.addStretch(0.5)
        vbox2.addWidget(text2)   
        vbox2.addWidget(text3)   
        vbox2.addWidget(text4)  
        hbox.addLayout(vbox2)
        hbox.addStretch(1)
        vbox.addLayout(hbox)
        vbox.addStretch(1) 
        
        button_close = QtGui.QPushButton('Close')
        button_close.clicked.connect( self.close)
        vbox.addWidget(button_close)
        vbox.addStretch(0.5)
                        
        
        
        self.setLayout(vbox)
        
""" ------------------------------------------------------------------------------------------------"""
class MainFrame(QtGui.QMainWindow):
    
    def __init__(self):
        super(MainFrame, self).__init__()
        
        self.initUI()
        


#----------------------------------------------------------------------          
    def initUI(self):   
        
        self.data_struct = data_struct.h5()
        self.stk = data_stack.data(self.data_struct)
        self.anlz = analyze.analyze(self.stk)
        self.nnma = nnma.nnma(self.stk)
        self.common = common()
               

        self.resize(Winsizex, Winsizey)
        self.setWindowTitle('Mantis')
        
        self.initToolbar()

                            
        ico = QtGui.QIcon(resource_path(os.path.join('images','logo-2l-32.ico')))
        self.setWindowIcon(ico)  
        
        tabs = QtGui.QTabWidget()
        


        
        # create the page windows as tabs
        self.page0 = PageLoadData(self.common, self.data_struct, self.stk)
        self.page1 = PageStack(self.common, self.data_struct, self.stk)
        self.page2 = PagePCA(self.common, self.data_struct, self.stk, self.anlz)
        self.page3 = PageCluster(self.common, self.data_struct, self.stk, self.anlz)
        self.page4 = PageSpectral(self.common, self.data_struct, self.stk, self.anlz)
        self.page5 = PagePeakID(self.common, self.data_struct, self.stk, self.anlz)
        self.page6 = PageXrayPeakFitting(self.common, self.data_struct, self.stk, self.anlz)
        self.page7 = PageNNMA(self.common, self.data_struct, self.stk, self.anlz, self.nnma)
        
        
        tabs.addTab(self.page0,"Load Data")
        tabs.addTab(self.page1,"Preprocess Data")
        tabs.addTab(self.page2,"PCA")
        tabs.addTab(self.page3,"Cluster Analysis")
        tabs.addTab(self.page4,"Spectral Maps")
        tabs.addTab(self.page7, "NNMA Analysis")
        tabs.addTab(self.page5,"Peak ID")
        tabs.addTab(self.page6, "XrayPeakFitting")  
         
        if showtomotab:
            self.page8 = PageTomo(self.common, self.data_struct, self.stk, self.anlz)
            tabs.addTab(self.page8, "Tomography")  
            
                    
        tabs.tabBar().setTabTextColor(0, QtGui.QColor('green'))
        tabs.tabBar().setTabTextColor(1, QtGui.QColor('green'))
        tabs.tabBar().setTabTextColor(2, QtGui.QColor('darkRed'))
        tabs.tabBar().setTabTextColor(3, QtGui.QColor('darkRed'))
        tabs.tabBar().setTabTextColor(4, QtGui.QColor('darkRed'))  
        tabs.tabBar().setTabTextColor(5, QtGui.QColor('darkRed'))      
        tabs.tabBar().setTabTextColor(6, QtGui.QColor('purple'))
        tabs.tabBar().setTabTextColor(7, QtGui.QColor('purple')) 
        if showtomotab:
            tabs.tabBar().setTabTextColor(8, QtGui.QColor('darkblue')) 

       
        
        
        # Only add "expert" pages if option "--key" is given in command line
        try:
            options, extraParams = getopt.getopt(sys.argv[1:], '', ['wx', 'batch', 'nnma', 'ica', 'keyeng'])
        except:
            print 'Error - wrong command line option used. Available options are --wx, --batch and --nnma'
            return
        
#         for opt, arg in options:        
#             if opt in '--nnma':
#                 if verbose: print "Running with NNMA."
#                 self.page7 = PageNNMA(self.common, self.data_struct, self.stk, self.anlz, self.nnma)
#                 tabs.addTab(self.page7, "NNMA Analysis")

        
        

    
        layout = QVBoxLayout()

        layout.addWidget(tabs)
        self.setCentralWidget(tabs)
       
                              
        self.show()
        if sys.platform == "darwin":
            self.raise_()

        
        

#----------------------------------------------------------------------   
    def initToolbar(self):   
        
        self.actionOpen = QtGui.QAction(self)
        self.actionOpen.setObjectName('actionOpen')
        self.actionOpen.setIcon(QtGui.QIcon(resource_path(os.path.join('images','document-open.png'))))
        self.toolbar = self.addToolBar('actionOpen') 
        self.toolbar.addAction(self.actionOpen)
        self.actionOpen.triggered.connect(self.LoadStack)
        
        self.actionOpenSL = QtGui.QAction(self)
        self.actionOpenSL.setObjectName('actionOpenSL')
        self.actionOpenSL.setIcon(QtGui.QIcon(resource_path(os.path.join('images','open-sl.png'))))
        self.toolbar.addAction(self.actionOpenSL)
        self.actionOpenSL.triggered.connect(self.BuildStack)
        
        self.actionSave = QtGui.QAction(self)
        self.actionSave.setObjectName('actionSave')
        self.actionSave.setIcon(QtGui.QIcon(resource_path(os.path.join('images','media-floppy.png'))))
        self.toolbar.addAction(self.actionSave)
        self.actionSave.triggered.connect(self.onSaveResultsToH5)
        self.actionSave.setEnabled(False)

        self.actionInfo = QtGui.QAction(self)
        self.actionInfo.setObjectName('actionInfo')
        self.actionInfo.setIcon(QtGui.QIcon(resource_path(os.path.join('images','help-browser.png'))))
        self.toolbar.addAction(self.actionInfo)
        self.actionInfo.triggered.connect(self.onAbout)
        
#----------------------------------------------------------------------
    def LoadStack(self):
        """
        Browse for a stack file:
        """

        filepath, plugin = File_GUI.SelectFile('read','stack')
        if filepath is not None:
            if plugin is None:
                plugin = file_plugins.identify(filepath)
            FileStruct = file_plugins.GetFileStructure(filepath, plugin=plugin)
            FileInternalSelection = [(0,0)]
            if FileStruct is not None:
                dlg = File_GUI.DataChoiceDialog(filepath=filepath, filestruct=FileStruct)
                if not dlg.exec_():
                    return # do nothing if GUI is cancelled
                FileInternalSelection = dlg.selection
                if dlg.filepath != filepath:
                    filepath = dlg.filepath
                    plugin = file_plugins.identify(dlg.filepath)
                    
                    
            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
            
            if self.common.stack_loaded == 1:
                self.new_stack_refresh()
                self.stk.new_data()
                self.anlz.delete_data()
            file_plugins.load(filepath, stack_object=self.stk, plugin=plugin,selection=FileInternalSelection)
            directory =  os.path.dirname(str(filepath))
            self.page1.filename =  os.path.basename(str(filepath))
                
            
        
            #Update widgets 
            x=self.stk.n_cols
            y=self.stk.n_rows
            self.page1.imgrgb = np.zeros(x*y*3,dtype = "uint8")        
            self.page1.maxval = np.amax(self.stk.absdata)
            
            
            self.ix = int(x/2)
            self.iy = int(y/2)
            
            self.page1.ix = self.ix
            self.page1.iy = self.iy
            
            self.iev = int(self.stk.n_ev/2)
            self.page0.slider_eng.setRange(0,self.stk.n_ev-1)
            self.page0.iev = self.iev
            self.page0.slider_eng.setValue(self.iev)
                     
            self.page1.slider_eng.setRange(0,self.stk.n_ev-1)
            self.page1.iev = self.iev
            self.page1.slider_eng.setValue(self.iev)
            
            self.stk.scale_bar()
            self.common.stack_loaded = 1
            
            if self.stk.data_struct.spectromicroscopy.normalization.white_spectrum is not None:
                self.common.i0_loaded = 1
                self.stk.calculate_optical_density()
                self.stk.fill_h5_struct_normalization()
            

            self.page1.ResetDisplaySettings()
            self.page1.loadImage()
            #print (x,y), (self.ix,self.iy), self.stk.absdata.shape
            self.page1.loadSpectrum(self.ix, self.iy)
            self.page1.textctrl.setText(self.page1.filename)
            
            self.page0.ShowInfo(self.page1.filename, directory)
            
            self.page5.updatewidgets()

            QtGui.QApplication.restoreOverrideCursor()
                 
        self.refresh_widgets()


#----------------------------------------------------------------------
    def BuildStack(self):
        """
        Browse for .sm files
        """
        
        #try:
        if True:
            directory = QtGui.QFileDialog.getExistingDirectory(self, "Choose a directory", '', QtGui.QFileDialog.ShowDirsOnly|QtGui.QFileDialog.ReadOnly )       
                                                        
        
       
            if directory == '':
                return
                 
            self.common.path = directory
 
            stackframe = StackListFrame(self, directory, self.common, self.stk, self.data_struct)
            stackframe.show()
             
#         except:
#             print 'Error could not build stack list.'
#             self.common.stack_loaded = 0 
#             self.common.i0_loaded = 0
#             self.new_stack_refresh()
#             self.refresh_widgets()
#                                    
#             QtGui.QMessageBox.warning(self,'Error',"Error could not build stack list")
#             import sys; print sys.exc_info()
            
            
#----------------------------------------------------------------------
    def LoadStack4D(self):

        try:
        #if True:

            wildcard =  "Supported 4D formats (*.hdf5 *.ncb);;HDF5 files (*.hdf5);;NCB files (*.ncb);;"  

            filepath = QtGui.QFileDialog.getOpenFileNames(self, 'Open files', '', wildcard)

            
            if filepath == '':
                return
            
            filenames = []
            for name in filepath:
                filenames.append(str(name))
                       
            
            directory =  os.path.dirname(str(filenames[0]))
            self.page4.DefaultDir = directory
            self.page1.filename =  os.path.basename(str(filenames[0]))
        
            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
            basename, extension = os.path.splitext(self.page1.filename)      
            
            self.common.path = directory
            self.common.filename = self.page1.filename
            

            if extension == '.hdf5':
                if self.common.stack_loaded == 1:
                    self.new_stack_refresh()  
                    self.stk.new_data()
                    self.anlz.delete_data()                
 
 
                self.stk.read_h54D(filenames[0])
            
                            
            elif extension == '.ncb':              
                if self.common.stack_loaded == 1:
                    self.new_stack_refresh()  
                    self.stk.new_data()
                    #self.stk.data_struct.delete_data()
                    self.anlz.delete_data()       
                self.stk.read_ncb4D(filenames)    
                #Get energy list
                wildcard =  "Text file (*.txt);;"  
                engfilepath = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '', wildcard)
                self.stk.read_ncb4Denergy(str(engfilepath))
                
            
            #Update widgets 
            
            self.common.stack_4d = 1
        
            x=self.stk.n_cols
            y=self.stk.n_rows  
            self.page1.imgrgb = np.zeros(x*y*3,dtype = "uint8")        
            self.page1.maxval = np.amax(self.stk.absdata)
            
            
            self.ix = int(x/2)
            self.iy = int(y/2)
            
            self.page1.ix = self.ix
            self.page1.iy = self.iy
            
            self.iev = int(self.stk.n_ev/2)
            self.page0.slider_eng.setRange(0,self.stk.n_ev-1)
            self.page0.iev = self.iev
            self.page0.slider_eng.setValue(self.iev)
                     
            self.page1.slider_eng.setRange(0,self.stk.n_ev-1)
            self.page1.iev = self.iev
            self.page1.slider_eng.setValue(self.iev)
            
            self.page0.slider_theta.setVisible(True)  
            self.page0.tc_imagetheta.setVisible(True)
            self.itheta = 0
            self.page0.slider_theta.setRange(0,self.stk.n_theta-1)
            self.page0.itheta = self.itheta
            self.page0.slider_theta.setValue(self.itheta)
            self.page0.tc_imagetheta.setText("4D Data Angle: "+str(self.stk.theta[self.itheta]))


            self.page1.slider_theta.setVisible(True)  
            self.page1.tc_imagetheta.setVisible(True)
            self.page1.slider_theta.setRange(0,self.stk.n_theta-1)
            self.page1.itheta = self.itheta
            self.page1.slider_theta.setValue(self.itheta)
            self.page1.tc_imagetheta.setText("4D Data Angle: "+str(self.stk.theta[self.itheta]))
            
            
            self.page2.button_calcpca4D.setVisible(True) 
            self.page4.button_calc4d.setVisible(True)
                        

            self.common.stack_loaded = 1
            
            if self.stk.data_struct.spectromicroscopy.normalization.white_spectrum is not None:
                self.common.i0_loaded = 1
                
            
            self.page1.ResetDisplaySettings()
            self.page1.loadImage()
            self.page1.loadSpectrum(self.ix, self.iy)
            self.page1.textctrl.setText(self.page1.filename)
            
            self.page0.ShowInfo(self.page1.filename, directory)
            
            self.page5.updatewidgets()

            QtGui.QApplication.restoreOverrideCursor()
                 
        except:
         
            self.common.stack_loaded = 0 
            self.common.i0_loaded = 0
            self.new_stack_refresh()
                                        
            QtGui.QApplication.restoreOverrideCursor()
            QtGui.QMessageBox.warning(self, 'Error', 'Image stack not loaded.')
        
            import sys
            print sys.exc_info()
                   
        
        self.refresh_widgets()
            
#----------------------------------------------------------------------
    def OnSaveProcessedStack(self, event):
        self.SaveProcessedStack()
       
        
#----------------------------------------------------------------------
    def SaveProcessedStack(self):

        """
        Browse for .hdf5 file or .ncb or tiff or .stk
        """
        
        #try:
        if True:
            #wildcard = "HDF5 file (*.hdf5);;aXis2000 NCB file (*.ncb);;TIFF file (.tif);;STK file (*.stk);;"
            wildcard = "HDF5 file (*.hdf5);;aXis2000 NCB file (*.ncb);;TIFF file (.tif);;STK file (*.stk);;"

            filepath = QtGui.QFileDialog.getSaveFileName(self, 'Save processed stack', '', wildcard)

            filepath = str(filepath)
            if filepath == '':
                return
               
            directory =  os.path.dirname(str(filepath))
        
            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))

            
            basename, extension = os.path.splitext(filepath)      
                       
            if extension == '.hdf5':            
                file_dataexch_hdf5.write_h5(filepath, self.data_struct) 
                           
#             elif extension == '.ncb':    
#                 self.stk.write_ncb(filepath, self.data_struct) 
#            
#             elif extension == '.tif':    
#                 self.stk.write_tif(filepath, self.stk.absdata, energies=self.stk.ev) 
#                 
#             elif extension == '.stk':    
#                 self.stk.write_stk(filepath) 
            
         
            QtGui.QApplication.restoreOverrideCursor()

#         except:
#       
#             QtGui.QApplication.restoreOverrideCursor()
#                  
#             QtGui.QMessageBox.warning(self, 'Error', 'Could not save processed stack file.')
                   

        self.refresh_widgets()
        
        return

#----------------------------------------------------------------------
    def onSaveResultsToH5(self, event):
        self.SaveResultsToH5()     
        
#----------------------------------------------------------------------
    def SaveResultsToH5(self):

        """
        Browse for .hdf5 file
        """

        try:

            wildcard = "HDF5 files (*.hdf5)"

            filepath = QtGui.QFileDialog.getSaveFileName(self, 'Save as .hdf5', '', wildcard)

            filepath = str(filepath)
            if filepath == '':
                return
            
            
            directory =  os.path.dirname(str(filepath))
            self.page1.filename =  os.path.basename(str(filepath))
        
            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
            
            self.common.path = directory
            self.common.filename = self.page1.filename
            

            file_dataexch_hdf5.write_results_h5(filepath, self.data_struct, self.anlz)   

            QtGui.QApplication.restoreOverrideCursor()

        except:
   
            QtGui.QApplication.restoreOverrideCursor()
            QtGui.QMessageBox.warning(self, 'Error', 'Could not save HDF5 file.')
                   

        self.refresh_widgets()
        
        return
                  
#----------------------------------------------------------------------
    def onAbout(self):

        self.popup = AboutFrame(self)
        self.popup.show()


        
#----------------------------------------------------------------------        
    def refresh_widgets(self):
        
        

        if self.common.stack_loaded == 0:
            self.page1.button_i0ffile.setEnabled(False)
            self.page1.button_i0histogram.setEnabled(False) 
            self.page1.button_prenorm.setEnabled(False)
            self.page1.button_refimgs.setEnabled(False)
            self.page1.button_limitev.setEnabled(False)
            self.page1.button_subregion.setEnabled(False)
            self.page1.button_darksig.setEnabled(False)
            self.page1.button_save.setEnabled(False) 
            self.page1.button_savestack.setEnabled(False)
            self.page1.button_align.setEnabled(False)
            self.page1.button_slideshow.setEnabled(False)
            self.page1.button_addROI.setEnabled(False)
            self.page1.button_spectralROI.setEnabled(False)
            self.page1.button_resetdisplay.setEnabled(False) 
            self.page1.button_despike.setEnabled(False)   
            self.page1.button_displaycolor.setEnabled(False)
            self.actionSave.setEnabled(False)      
                 
        else:
            self.page1.button_i0ffile.setEnabled(True)
            self.page1.button_i0histogram.setEnabled(True) 
            self.page1.button_prenorm.setEnabled(True)
            self.page1.button_limitev.setEnabled(True)
            if self.common.stack_4d == 0: 
                self.page1.button_refimgs.setEnabled(True)
                self.page1.button_subregion.setEnabled(True)
                self.page1.button_darksig.setEnabled(True)
            else:
                self.page1.button_refimgs.setEnabled(False)
                self.page1.button_subregion.setEnabled(False)
                self.page1.button_darksig.setEnabled(False)                
            self.page1.button_save.setEnabled(True) 
            self.page1.button_savestack.setEnabled(True)  
            self.page1.button_align.setEnabled(True)  
            self.page1.button_slideshow.setEnabled(True)
            self.page1.button_addROI.setEnabled(True)  
            self.page1.button_spectralROI.setEnabled(True)
            self.page1.button_resetdisplay.setEnabled(True) 
            self.page1.button_despike.setEnabled(True) 
            self.page1.button_displaycolor.setEnabled(True)
            self.actionSave.setEnabled(True)  
             
             
        if self.common.i0_loaded == 0:
            self.page1.button_showi0.setEnabled(False) 
            self.page1.rb_flux.setEnabled(False)
            self.page1.rb_od.setEnabled(False)
            self.page1.button_reseti0.setEnabled(False)
            self.page1.button_saveod.setEnabled(False)
            self.page2.button_calcpca.setEnabled(False)
            self.page2.button_calcpca4D.setEnabled(False)
            self.page4.button_loadtspec.setEnabled(False)
            self.page4.button_addflat.setEnabled(False)
            self.page5.button_save.setEnabled(False)
            if self.page7:
                self.page7.button_calcnnma.setEnabled(False)
                self.page7.button_mufile.setEnabled(False)
                self.page7.button_murand.setEnabled(False)
            if showtomotab:
                self.page8.button_engdata.setEnabled(False)
        else:
            self.page1.button_showi0.setEnabled(True)
            self.page1.rb_flux.setEnabled(True)
            self.page1.rb_od.setEnabled(True)   
            self.page1.button_reseti0.setEnabled(True)
            self.page1.button_saveod.setEnabled(True)
            self.page2.button_calcpca.setEnabled(True) 
            self.page2.button_calcpca4D.setEnabled(True) 
            self.page4.button_loadtspec.setEnabled(True)
            self.page4.button_addflat.setEnabled(True)   
            self.page5.button_save.setEnabled(True)
            if self.page7:
                self.page7.button_calcnnma.setEnabled(True)
                self.page7.button_mufile.setEnabled(True)
                self.page7.button_murand.setEnabled(True)             
            if showtomotab:
                self.page8.button_engdata.setEnabled(True)    
             
        if self.common.pca_calculated == 0:      
            self.page2.button_savepca.setEnabled(False)
            self.page2.slidershow.setEnabled(False) 
            self.page2.button_movepcup.setEnabled(False)
            self.page3.button_calcca.setEnabled(False)
            self.page4.rb_fit.setEnabled(False)
        else:
            self.page2.button_savepca.setEnabled(True)
            self.page2.slidershow.setEnabled(True)
            self.page2.button_movepcup.setEnabled(True)
            self.page3.button_calcca.setEnabled(True)  
            self.page4.rb_fit.setEnabled(True)       
             
        if self.common.cluster_calculated == 0:   
            self.page3.button_scatterplots.setEnabled(False)
            self.page3.button_savecluster.setEnabled(False)
            self.page3.slidershow.setEnabled(False)
            self.page4.button_addclspec.setEnabled(False)
            self.page5.button_addclspec.setEnabled(False)
            if self.page7:
                self.page7.button_mucluster.setEnabled(False)
        else:
            self.page3.button_scatterplots.setEnabled(True)
            self.page3.button_savecluster.setEnabled(True)  
            self.page3.slidershow.setEnabled(True)
            self.page4.button_addclspec.setEnabled(True)
            self.page5.button_addclspec.setEnabled(True)
            if self.page7:
                self.page7.button_mucluster.setEnabled(True)            
            
             
        if self.common.spec_anl_calculated == 0:
            self.page4.button_removespec.setEnabled(False)
            self.page4.button_movespdown.setEnabled(False)
            self.page4.button_movespup.setEnabled(False)
            self.page4.button_save.setEnabled(False)
            self.page4.button_showrgb.setEnabled(False) 
            self.page4.button_histogram.setEnabled(False) 
            self.page4.button_calc4d.setEnabled(False)

        else:
            self.page4.button_removespec.setEnabled(True)
            self.page4.button_movespdown.setEnabled(True)
            self.page4.button_movespup.setEnabled(True)
            self.page4.button_save.setEnabled(True) 
            self.page4.button_showrgb.setEnabled(True)    
            self.page4.button_histogram.setEnabled(True) 
            self.page4.button_calc4d.setEnabled(True)
                
        
        if showtomotab:
            if self.common.spec_anl4D_calculated == 0:
                self.page8.button_spcomp.setEnabled(False)
            else:
                self.page8.button_spcomp.setEnabled(True)
            
            
        if self.page6 != None:
            if self.common.cluster_calculated == 0:   
                self.page6.button_addclspec.setEnabled(False)
            else:
                self.page6.button_addclspec.setEnabled(True)                
            
            if self.common.xpf_loaded == 1:
                self.page6.slider_spec.setEnabled(True)
            else:
                self.page6.slider_spec.setEnabled(False)
        
                                                      
             
        self.page1.ResetDisplaySettings()
             
             
            
#----------------------------------------------------------------------        
    def new_stack_refresh(self):
        
        
        self.common.i0_loaded = 0
        self.common.pca_calculated = 0
        self.common.cluster_calculated = 0
        self.common.spec_anl_calculated = 0
        self.common.xpf_loaded = 0
        self.common.stack_4d = 0
        self.common.pca4D_calculated = 0
        self.common.spec_anl4D_calculated = 0
        
        self.refresh_widgets()
        
        #page 0
        self.page0.slider_theta.setVisible(False)  
        self.page0.tc_imagetheta.setVisible(False)               
        
        #page 1
        self.page1.rb_flux.setChecked(True)
        self.page1.rb_od.setChecked(False)
        self.page1.showflux = True
         
        fig = self.page1.specfig
        fig.clf()
        self.page1.SpectrumPanel.draw()
        self.page1.tc_spec.setText("Spectrum at point: ")       
         
        fig = self.page1.absimgfig
        fig.clf()
        self.page1.AbsImagePanel.draw()        
        self.page1.tc_imageeng.setText("Image at energy: ")
         
        self.page1.textctrl.setText(' ')
         
        self.page1.ResetDisplaySettings()
        #page 0
        self.page1.slider_theta.setVisible(False)  
        self.page1.tc_imagetheta.setVisible(False)            
         
        #page 2
        fig = self.page2.pcaevalsfig
        fig.clf()
        self.page2.PCAEvalsPan.draw()
         
        fig = self.page2.pcaimgfig
        fig.clf()
        self.page2.PCAImagePan.draw()
         
        fig = self.page2.pcaspecfig
        fig.clf()
        self.page2.PCASpecPan.draw()
         
        self.page2.vartc.setText('0%')
        self.page2.npcaspin.setValue(1)
        self.page2.tc_PCAcomp.setText("PCA component ")
        self.page2.text_pcaspec.setText("PCA spectrum ")
         
        self.page2.selpca = 1       
        self.page2.numsigpca = 2
        self.page2.slidershow.setValue(self.page2.selpca)
        
        self.page2.slider_theta.setVisible(False)  
        self.page2.tc_imagetheta.setVisible(False)  
        self.page2.button_calcpca4D.setVisible(False) 
         
        #page 3
        fig = self.page3.clusterimgfig
        fig.clf()
        self.page3.ClusterImagePan.draw()
         
        fig = self.page3.clusterindvimgfig
        fig.clf()
        self.page3.ClusterIndvImagePan.draw()
         
        fig = self.page3.clusterspecfig
        fig.clf()
        self.page3.ClusterSpecPan.draw()
         
        fig = self.page3.clusterdistmapfig
        fig.clf()
        self.page3.ClusterDistMapPan.draw()       
        
        self.page3.selcluster = 1
        self.page3.slidershow.setValue(self.page3.selcluster)
        self.page3.numclusters = 5
        self.page3.nclusterspin.setValue(self.page3.numclusters)
        self.page3.tc_cluster.setText("Cluster ")
        self.page3.tc_clustersp.setText("Cluster spectrum")
        self.page3.wo_1st_pca = 0
        self.page3.remove1stpcacb.setChecked(False)
 
         
        #page 4
        self.page4.ClearWidgets()
         
        #page 5
        fig = self.page5.kespecfig
        fig.clf()
        self.page5.KESpecPan.draw()         
        fig = self.page5.absimgfig
        fig.clf()
        self.page5.AbsImagePanel.draw()   
 
        
        #page 6
        if self.page6 != None:
            fig = self.page6.Specfig
            fig.clf()
            self.page6.SpectrumPanel.draw()
            self.page6.slider_spec.setEnabled(False)
            
        #page8
        if showtomotab:
            self.page8.NewStackClear()
        
                
""" ------------------------------------------------------------------------------------------------"""
                        
def main():
    
    
    app = QtGui.QApplication(sys.argv)
    frame = MainFrame()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
