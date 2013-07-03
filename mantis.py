'''
Created on Jun 11, 2013

@author: Mirna Lerotic
'''

import sys
import os
import numpy as npy

from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import *
from PyQt4.QtCore import Qt


import matplotlib 
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid import make_axes_locatable

import data_struct
import data_stack
import analyze
import logos
import nnma
import henke


Winsizex = 1000
Winsizey = 700

PlotH = 4.0
PlotW = PlotH*1.61803

#----------------------------------------------------------------------
class common:
    def __init__(self):
        self.fontsize = 8
        
        self.stack_loaded = 0
        self.i0_loaded = 0
        self.pca_calculated = 0
        self.cluster_calculated = 0
        self.spec_anl_calculated = 0
        self.ica_calculated = 0
        
        self.path = ''
        self.filename = ''

        self.font = ''
        
""" ------------------------------------------------------------------------------------------------"""
class PageKeyEng(QtGui.QWidget):
    def __init__(self):
        super(PageKeyEng, self).__init__()

        self.initUI()
        
#----------------------------------------------------------------------          
    def initUI(self): 
        pass
    

""" ------------------------------------------------------------------------------------------------"""
class PageSpectral(QtGui.QWidget):
    def __init__(self):
        super(PageSpectral, self).__init__()

        self.initUI()
        
#----------------------------------------------------------------------          
    def initUI(self): 
        pass
    
    
""" ------------------------------------------------------------------------------------------------"""
class PageCluster(QtGui.QWidget):
    def __init__(self):
        super(PageCluster, self).__init__()

        self.initUI()
        
#----------------------------------------------------------------------          
    def initUI(self): 
        pass
    
    
""" ------------------------------------------------------------------------------------------------"""
class PagePCA(QtGui.QWidget):
    def __init__(self):
        super(PagePCA, self).__init__()

        self.initUI()
        
#----------------------------------------------------------------------          
    def initUI(self): 


        
        #panel 1        
        vbox1 = QtGui.QVBoxLayout()
        
        self.tc_PCAcomp = QtGui.QLabel(self)
        #self.tc_PCAcomp.SetFont(self.com.font)
        self.tc_PCAcomp.setText("PCA component ")
        vbox1.addWidget(self.tc_PCAcomp)
        
        hbox11 = QtGui.QHBoxLayout()
   
        #i1panel = wx.Panel(panel1, -1, style = wx.SUNKEN_BORDER)
        #self.PCAImagePan = wxmpl.PlotPanel(i1panel, -1, size =(ph*1.10, ph), cursor=False, crosshairs=False, location=False, zoom=False)
        self.image = QtGui.QImage(os.path.join('images','Mantis_img.jpg'))
        
        self.imageLabel = QtGui.QLabel()
        self.imageLabel.setBackgroundRole(QtGui.QPalette.Base)
        #self.imageLabel.setSizePolicy()
        #self.imageLabel.setScaledContents(True)
        
        self.imageLabel.setPixmap(QtGui.QPixmap.fromImage(self.image))
        hbox11.addWidget(self.imageLabel)
                       
            
        self.slidershow = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slidershow.setFocusPolicy(QtCore.Qt.StrongFocus)
        #self.slidershow.Disable()          
        #self.slidershow.SetFocus()
        #self.Bind(wx.EVT_SCROLL, self.OnPCAScroll, self.slidershow)
        
        hbox11.addWidget(self.slidershow)

        vbox1.addLayout(hbox11)
        
       
                
        #panel 2
        vbox2 = QtGui.QVBoxLayout()
        sizer2 = QtGui.QGroupBox('PCA')
        vbox21 = QtGui.QVBoxLayout()        
        
        self.button_calcpca = QtGui.QPushButton('Calculate PCA')
        #self.button_calcpca.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnCalcPCA, id=self.button_calcpca.GetId())     
        #self.button_calcpca.Disable()   
        vbox21.addWidget(self.button_calcpca)
        self.button_savepca = QtGui.QPushButton('Save PCA Results...')
        #self.button_savepca.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnSave, id=self.button_savepca.GetId())
        #self.button_savepca.Disable()
        vbox21.addWidget(self.button_savepca)
        
        hbox21 = QtGui.QHBoxLayout()
        text1 = QtGui.QLabel(self)
        text1.setText('Number of significant components')
        #text1.SetFont(self.com.font)
        #self.npcaspin = wx.SpinCtrl(panel2, -1, '',  size= (60, -1), style=wx.ALIGN_LEFT)
        #self.npcaspin.SetRange(1,20)
        #self.Bind(wx.EVT_SPINCTRL, self.OnNPCAspin, self.npcaspin)
        vbox21.addWidget(text1)
              
        hbox22 = QtGui.QHBoxLayout()
        text2 = QtGui.QLabel(self)
        text2.setText( 'Cumulative variance')
        #text2.SetFont(self.com.font)
        self.vartc = QtGui.QLabel(self)
        self.vartc.setText('0%')
        hbox22.addWidget(text2)
        hbox22.addWidget(self.vartc)

        vbox21.addLayout(hbox22)      
        sizer2.setLayout(vbox21)
        vbox2.addStretch(1)
        vbox2.addWidget(sizer2)
        vbox2.addStretch(3)

        
        #panel 3
        vbox3 = QtGui.QVBoxLayout()
   
        #i3panel = wx.Panel(panel3, -1, style = wx.SUNKEN_BORDER)
        #self.PCASpecPan = wxmpl.PlotPanel(i3panel, -1, size =(pw, ph), cursor=False, crosshairs=False, location=False, zoom=False)
 
        self.text_pcaspec = QtGui.QLabel(self)
        #self.text_pcaspec.SetFont(self.com.font)
        self.text_pcaspec.setText("PCA spectrum ") 
        vbox3.addWidget(self.text_pcaspec)
                        
        image = QtGui.QImage(os.path.join('images','spectrum.jpg'))       
        self.imageLabel = QtGui.QLabel()
        self.imageLabel.setBackgroundRole(QtGui.QPalette.Base)
        #self.imageLabel.setSizePolicy()
        #self.imageLabel.setScaledContents(True)
        self.imageLabel.setPixmap(QtGui.QPixmap.fromImage(image))
        vbox3.addWidget(self.imageLabel)   

 
    
        
        #panel 4
        vbox4 = QtGui.QVBoxLayout()
 
        #i4panel = wx.Panel(panel4, -1, style = wx.SUNKEN_BORDER)
        #self.PCAEvalsPan = wxmpl.PlotPanel(i4panel, -1, size =(pw, ph*0.75), cursor=False, crosshairs=False, location=False, zoom=False)
        #wxmpl.EVT_POINT(i4panel, self.PCAEvalsPan.GetId(), self.OnPointEvalsImage)   
         
        text4 = QtGui.QLabel(self)
        #text4.SetFont(self.com.font)
        text4.setText("PCA eigenvalues ")        
        vbox4.addWidget(text4)
        
        image2 = QtGui.QImage(os.path.join('images','spectrum.jpg'))       
        self.imageLabel2 = QtGui.QLabel()
        self.imageLabel2.setBackgroundRole(QtGui.QPalette.Base)
        #self.imageLabel.setSizePolicy()
        #self.imageLabel.setScaledContents(True)
        self.imageLabel2.setPixmap(QtGui.QPixmap.fromImage(image2))
        vbox4.addWidget(self.imageLabel2)           


    
        
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
        self.showflux = True
        self.fontsize = self.com.fontsize
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
        
        self.movie_playing = 0
        

        #panel 1
        sizer1 = QtGui.QGroupBox('Preprocess')
        vbox1 = QtGui.QVBoxLayout()
        vbox1.setSpacing(0)
        
        self.button_align = QtGui.QPushButton('Align images...')
        #self.Bind(wx.EVT_BUTTON, self.OnAlignImgs, id=self.button_align.GetId())
        #self.button_align.SetFont(self.com.font)
        #self.button_align.Disable()
        vbox1.addWidget(self.button_align) 
        
        self.button_i0ffile = QtGui.QPushButton('I0 from file...')
        #self.button_i0ffile.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnI0FFile, id=self.button_i0ffile.GetId())
        #self.button_i0ffile.Disable()
        vbox1.addWidget(self.button_i0ffile)
        self.button_i0histogram = QtGui.QPushButton('I0 from histogram...')
        #self.button_i0histogram.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnI0histogram, id=self.button_i0histogram.GetId())   
        #self.button_i0histogram.Disable()     
        vbox1.addWidget(self.button_i0histogram)
        self.button_showi0 = QtGui.QPushButton('Show I0...')
        #self.button_showi0.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnShowI0, id=self.button_showi0.GetId())   
        #self.button_showi0.Disable()
        vbox1.addWidget(self.button_showi0)
        
        
        self.button_limitev = QtGui.QPushButton('Limit energy range...')
        #self.Bind(wx.EVT_BUTTON, self.OnLimitEv, id=self.button_limitev.GetId())
        #self.button_limitev.SetFont(self.com.font)
        #self.button_limitev.Disable()
        vbox1.addWidget(self.button_limitev)
        
        self.button_subregion = QtGui.QPushButton('Clip to subregion...')
        #self.Bind(wx.EVT_BUTTON, self.OnCliptoSubregion, id=self.button_subregion.GetId())
        #self.button_subregion.SetFont(self.com.font)
        #self.button_subregion.Disable()
        vbox1.addWidget(self.button_subregion)       
            
        self.button_savestack = QtGui.QPushButton('Save preprocessed stack')
        #self.button_savestack.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnSaveStack, id=self.button_savestack.GetId())
        #self.button_savestack.Disable()          
        vbox1.addWidget(self.button_savestack)
        
        sizer1.setLayout(vbox1)


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
        #self.Bind(wx.EVT_BUTTON, self.OnSlideshow, id=self.button_slideshow.GetId())
        #self.button_slideshow.SetFont(self.com.font)
        #self.button_slideshow.Disable()
        vbox22.addWidget(self.button_slideshow)
        
        self.button_save = QtGui.QPushButton( 'Save images...')
        self.button_save.setMaximumSize (150 , 150)
        #self.button_save.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnSave, id=self.button_save.GetId())
        #self.button_save.Disable()          
        vbox22.addWidget(self.button_save)
              
        hbox20.addLayout(vbox22)
        vbox2.addLayout(hbox20)
        vbox2.addSpacing(5)
        

        hbox21 = QtGui.QHBoxLayout()
        
        sizer22 = QtGui.QGroupBox('Image')
        vbox23 = QtGui.QVBoxLayout()

        self.rb_flux = QtGui.QRadioButton( 'Flux', self)
        self.rb_od = QtGui.QRadioButton('Optical Density',self)
        #self.rb_flux.SetFont(self.com.font)
        #self.rb_od.SetFont(self.com.font)
        #self.Bind(wx.EVT_RADIOBUTTON, self.onrb_fluxod, id=self.rb_flux.GetId())
        #self.Bind(wx.EVT_RADIOBUTTON, self.onrb_fluxod, id=self.rb_od.GetId())
        vbox23.addWidget(self.rb_flux)
        vbox23.addWidget(self.rb_od)
        vbox23.addStretch (1)
         
        #self.rb_flux.Disable()
        #self.rb_od.Disable()
 
        self.add_scale_cb = QtGui.QCheckBox('Scalebar', self) 
        #self.add_scale_cb.SetFont(self.com.font)
        #self.Bind(wx.EVT_CHECKBOX, self.OnShowScale, self.add_scale_cb)
        #self.add_scale_cb.SetValue(True)
        vbox23.addWidget(self.add_scale_cb)
         
        self.add_colbar_cb = QtGui.QCheckBox('Colorbar', self)
        #self.add_colbar_cb.SetFont(self.com.font)
        #self.Bind(wx.EVT_CHECKBOX, self.OnShowColBar, self.add_colbar_cb)
        vbox23.addWidget(self.add_colbar_cb)
        sizer22.setLayout(vbox23)
         
        hbox21.addWidget(sizer22)
        hbox21.addSpacing(5)
        vbox2.addLayout(hbox21)
                 
        sizer23 = QtGui.QGroupBox('Display settings')
        hbox23 = QtGui.QHBoxLayout() 
 
        fgs21 = QtGui.QGridLayout()
         
        self.tc_min = QtGui.QLabel(self)
        #self.tc_min.SetFont(self.com.font)
        self.tc_min.setText('Minimum: \t{0:5d}%'.format(int(100*self.brightness_min)))
         
        self.tc_max = QtGui.QLabel(self)
        #self.tc_max.SetFont(self.com.font)
        self.tc_max.setText('Maximum:{0:5d}%'.format(int(100*self.brightness_max)))
 
        self.slider_brightness_min = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        #self.slider_brightness_min = wx.Slider(panel4, -1, self.dispbrightness_min, 0, 49, size= (50,20), style=wx.SL_HORIZONTAL)        
        #self.slider_brightness_min.SetFocus()
        #self.Bind(wx.EVT_SCROLL, self.OnScrollBrightnessMin, self.slider_brightness_min)
         
        self.slider_brightness_max = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        #self.slider_brightness_max = wx.Slider(panel4, -1, self.dispbrightness_max, 50, 120, style=wx.SL_HORIZONTAL)        
        #self.slider_brightness_max.SetFocus()
        #self.Bind(wx.EVT_SCROLL, self.OnScrollBrightnessMax, self.slider_brightness_max)        
         
        self.tc_gamma = QtGui.QLabel(self)
        #self.tc_gamma.SetFont(self.com.font)
        self.tc_gamma.setText('Gamma:  \t{0:5.2f}'.format(self.gamma))
         
        self.slider_gamma = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        #self.slider_gamma = wx.Slider(panel4, -1, self.displaygamma, 1, 20, style=wx.SL_HORIZONTAL)        
        #self.slider_gamma.SetFocus()
        #self.Bind(wx.EVT_SCROLL, self.OnScrollGamma, self.slider_gamma)
       
         
        fgs21.addWidget(self.tc_min, 0, 0)
        fgs21.addWidget(self.slider_brightness_min, 0, 1)
        fgs21.addWidget(self.tc_max, 1, 0) 
        fgs21.addWidget(self.slider_brightness_max, 1, 1)
        fgs21.addWidget(self.tc_gamma, 2, 0)
        fgs21.addWidget(self.slider_gamma, 2, 1)
         
        hbox23.addLayout(fgs21)
                
    
        vbox24 = QtGui.QVBoxLayout()
        self.button_despike = QtGui.QPushButton('Despike')
        #self.button_despike.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnDespike, id=self.button_despike.GetId())   
        #self.button_despike.Disable()     
        vbox24.addWidget(self.button_despike)
        self.button_resetdisplay = QtGui.QPushButton( 'Reset')
        #self.button_resetdisplay.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.onResetDisplaySettings, id=self.button_resetdisplay.GetId())   
        #self.button_resetdisplay.Disable()     
        vbox24.addWidget(self.button_resetdisplay)
        self.button_displaycolor = QtGui.QPushButton('Color Table...   ')
        #self.button_displaycolor.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.onSetColorTable, id=self.button_displaycolor.GetId())   
        #self.button_displaycolor.Disable()     
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
        #self.button_addROI.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnAddROI, id=self.button_addROI.GetId())
        #self.button_addROI.Disable()
        vbox3.addWidget(self.button_addROI)
        
        self.button_acceptROI = QtGui.QPushButton('Accept ROI')
        #self.button_acceptROI.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnAcceptROI, id=self.button_acceptROI.GetId())   
        #self.button_acceptROI.Disable()     
        vbox3.addWidget(self.button_acceptROI)
        
        self.button_resetROI = QtGui.QPushButton('Reset ROI')
        #self.button_resetROI.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnResetROI, id=self.button_resetROI.GetId())
        #self.button_resetROI.Disable()
        vbox3.addWidget(self.button_resetROI) 

        self.button_setROII0 = QtGui.QPushButton('Set ROI As I0')
        #self.button_setROII0.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnSetROII0, id=self.button_setROII0.GetId())
        #self.button_setROII0.Disable()
        vbox3.addWidget(self.button_setROII0)
        
        self.button_saveROIspectr = QtGui.QPushButton( 'Save ROI Spectrum...')
        #self.button_saveROIspectr.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnSaveROISpectrum, id=self.button_saveROIspectr.GetId())   
        #self.button_saveROIspectr.Disable()     
        vbox3.addWidget(self.button_saveROIspectr)
        
        self.button_ROIdosecalc = QtGui.QPushButton('ROI Dose Calculation...')
        #self.button_ROIdosecalc.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnROI_DoseCalc, id=self.button_ROIdosecalc.GetId())   
        #self.button_ROIdosecalc.Disable()     
        vbox3.addWidget(self.button_ROIdosecalc)        
        
        self.button_spectralROI = QtGui.QPushButton('Spectral ROI...')
        #self.button_spectralROI.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnSpectralROI, id=self.button_spectralROI.GetId())   
        #self.button_spectralROI.Disable()     
        vbox3.addWidget(self.button_spectralROI)
        sizer3.setLayout(vbox3)
        

       
        
        

        #panel 4     
        vbox4 = QtGui.QVBoxLayout()
        
        self.tc_imageeng = QtGui.QLabel(self)
        #self.tc_imageeng.SetFont(self.com.font)
        self.tc_imageeng.setText("Image at energy: ")
        vbox4.addWidget(self.tc_imageeng)
        
        
        hbox41 = QtGui.QHBoxLayout()
        self.absimgfig = Figure((PlotH, PlotH))
        self.AbsImagePanel = FigureCanvas(self.absimgfig)
        self.AbsImagePanel.setParent(self)
        self.AbsImagePanel.mpl_connect('button_press_event', self.OnPointAbsimage)
        
        hbox41.addWidget(self.AbsImagePanel)   

       

        self.slider_eng = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slider_eng.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_eng.valueChanged[int].connect(self.OnScrollEng)
        self.slider_eng.setRange(0, 100)
        hbox41.addWidget(self.slider_eng)

        
        vbox4.addLayout(hbox41)
        


        #panel 5     
        vbox5 = QtGui.QVBoxLayout()
        
        self.tc_spec = QtGui.QLabel(self)
        #self.tc_spec.SetFont(self.com.font)
        self.tc_spec.setText("Spectrum ")
        vbox5.addWidget(self.tc_spec)
        
        

        self.specfig = Figure((PlotW, PlotH))
        self.SpectrumPanel = FigureCanvas(self.specfig)
        self.SpectrumPanel.setParent(self)
        self.SpectrumPanel.mpl_connect('button_press_event', self.OnPointSpectrum)
        vbox5.addWidget(self.SpectrumPanel)
               
    
        vboxtop = QtGui.QVBoxLayout()
        
        hboxtop = QtGui.QHBoxLayout()
        hboxtop.addWidget(sizer1)
        hboxtop.addWidget(sizer2)
        hboxtop.addWidget(sizer3)
        
        hboxbott = QtGui.QHBoxLayout()
        hboxbott2 = QtGui.QHBoxLayout()
        hboxbott2.addLayout(vbox4)
        hboxbott2.addStretch(1)
        hboxbott2.addLayout(vbox5)
        hboxbott.addLayout(hboxbott2)
              
        vboxtop.addStretch (1)
        vboxtop.addLayout(hboxtop)
        vboxtop.addStretch (1)
        vboxtop.addLayout(hboxbott)
        vboxtop.addStretch (1)


        vboxtop.setContentsMargins(20,20,20,20)
        self.setLayout(vboxtop)


#----------------------------------------------------------------------            
    def OnScrollEng(self, value):
        self.iev = value

        if self.com.stack_loaded == 1:
            self.loadImage()
            self.loadSpectrum(self.ix, self.iy)      

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
                indx = npy.abs(self.stk.ev - x).argmin()
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
            self.ix = int(npy.floor(y))           
            self.iy = int(npy.floor(x))  
                    
            if self.ix<0 :
                self.ix=0
            if self.ix>self.stk.n_cols :
                self.ix=self.stk.n_cols
            if self.iy<0 :
                self.iy=0
            if self.iy>self.stk.n_rows :
                self.iy=self.stk.n_rows 
            

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
                    self.button_acceptROI.Enable()
                self.loadImage()


#----------------------------------------------------------------------        
    def loadImage(self):
        
                  
        
        if self.defaultdisplay == 1.0:
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
            fig.add_axes((0.02,0.02,0.96,0.96))
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
            im = axes.imshow(image, cmap=matplotlib.cm.get_cmap(self.colortable)) 
        else:
            imgmax = npy.amax(image)
            imgmin = npy.amin(image)
            if (self.gamma != 1.0) or (imgmin < 0.0):
                image = (image-imgmin)/(imgmax-imgmin)
                imgmax = 1.0
                imgmin = 0.0
                if (self.gamma != 1.0):
                    image = npy.power(image, self.gamma)
            im = axes.imshow(image, cmap=matplotlib.cm.get_cmap(self.colortable), 
                             vmin=(imgmin+imgmax*self.brightness_min),vmax=imgmax*self.brightness_max)
             
        if (self.showROImask == 1) and (self.addroi == 1):
            im_red = axes.imshow(self.ROIpix_masked,cmap=matplotlib.cm.get_cmap("autumn")) 
          
        axes.axis("off") 
         
        if self.show_colorbar == 1:
            cbar = axes.figure.colorbar(im, orientation='vertical',cax=axcb) 
         
        if self.show_scale_bar == 1:
            startx = int(self.stk.n_rows*0.05)
            starty = self.stk.n_cols-int(self.stk.n_cols*0.05)-self.stk.scale_bar_pixels_y
            um_string = '$\mu m$'
            microns = '$'+self.stk.scale_bar_string+' $'+um_string
            axes.text(self.stk.scale_bar_pixels_x+startx+1,starty+1, microns, horizontalalignment='left', verticalalignment='center',
                      color = 'black', fontsize=14)
            #Matplotlib has flipped scales so I'm using rows instead of cols!
            p = matplotlib.patches.Rectangle((startx,starty), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                                   color = 'black', fill = True)
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
 
         
        matplotlib.rcParams['font.size'] = self.fontsize
 
        specplot = axes.plot(self.stk.ev,self.spectrum)
         
 
         
        axes.axvline(x=self.stk.ev[self.iev], color = 'g', alpha=0.5)
 
         
        self.SpectrumPanel.draw()
         
        self.tc_spec.setText("Spectrum at pixel [" +str(ypos)+", " + str(xpos)+"] or position ["+
                              str(self.stk.x_dist[xpos])+", "+ str(self.stk.y_dist[ypos])+ "]")


#----------------------------------------------------------------------
    def ResetDisplaySettings(self):
        
        pass

#         self.defaultdisplay = 1.0
#         
#         self.dispbrightness_min = 0
#         self.dispbrightness_max = 100
#         self.displaygamma = 10.0
#         
#         self.brightness_min = 0.0
#         self.brightness_max = 1.0
#         self.gamma = 1.0
#         
#         self.slider_brightness_max.SetValue(self.dispbrightness_max)
#         self.slider_brightness_min.SetValue(self.dispbrightness_min) 
#         self.slider_gamma.SetValue(self.displaygamma)      
# 
#         self.tc_min.Clear()
#         self.tc_min.AppendText('Minimum: \t{0:5d}%'.format(int(100*self.brightness_min)))
#         self.tc_max.Clear()
#         self.tc_max.AppendText('Maximum:{0:5d}%'.format(int(100*self.brightness_max)))        
#         self.tc_gamma.Clear()
#         self.tc_gamma.AppendText('Gamma:  \t{0:5.2f}'.format(self.gamma))      
                

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
       
        self.fontsize = self.com.fontsize
        
        self.iev = 0
        

        #panel 1
        sizer1 = QtGui.QGroupBox('Load Data Stack')
        vbox1 = QtGui.QVBoxLayout()
        
        self.button_hdf5 = QtGui.QPushButton('Load HDF5 Stack (*.hdf5)')
        self.button_hdf5.clicked.connect( self.OnLoadHDF5)
        vbox1.addWidget(self.button_hdf5)
        
        self.button_sdf = QtGui.QPushButton('Load SDF Stack (*.hrd)')
        self.button_sdf.clicked.connect( self.OnLoadSDF)   
        vbox1.addWidget(self.button_sdf)
        
        self.button_stk = QtGui.QPushButton('Load STK Stack (*.stk)')
        self.button_stk.clicked.connect( self.OnLoadSTK)   
        vbox1.addWidget(self.button_stk)
        
        self.button_xrm = QtGui.QPushButton('Load XRM Image (*.xrm)')
        self.button_xrm.clicked.connect( self.OnLoadXRM)
        vbox1.addWidget(self.button_xrm)
        
        self.button_txrm = QtGui.QPushButton('Load TXRM Stack (*.txrm)')
        self.button_txrm.clicked.connect( self.OnLoadTXRM)
        vbox1.addWidget(self.button_txrm)    
        
        self.button_tif = QtGui.QPushButton( 'Load TIF Stack (*.tif)')
        self.button_tif.clicked.connect( self.OnLoadTIF)
        vbox1.addWidget(self.button_tif)      

        sizer1.setLayout(vbox1)

        #panel 2
        sizer2 = QtGui.QGroupBox('Build a stack from a set of files')
        vbox2 = QtGui.QVBoxLayout()

        self.button_sm = QtGui.QPushButton( 'Build a stack from a set of NetCDF (*.sm) files')
        self.button_sm.clicked.connect( self.OnBuildStack)
        vbox2.addWidget(self.button_sm)
        
        self.button_xrm_list = QtGui.QPushButton( 'Build a stack from a set of XRM (*.xrm) files')
        self.button_xrm_list.clicked.connect( self.OnBuildStack)
        vbox2.addWidget(self.button_xrm_list)
        
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
        
        
        hbox51 = QtGui.QHBoxLayout()
   
        self.absimgfig = Figure((PlotH*.9, PlotH*.9))
        self.AbsImagePanel = FigureCanvas(self.absimgfig)
        self.AbsImagePanel.setParent(self)
        hbox51.addWidget(self.AbsImagePanel)
        

        self.slider_eng = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slider_eng.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_eng.valueChanged[int].connect(self.OnScrollEng)
        self.slider_eng.setRange(0, 100)
        
        hbox51.addWidget(self.slider_eng)        

        vbox5.addLayout(hbox51)
        
       
        vboxtop = QtGui.QVBoxLayout()
        
        hboxtop = QtGui.QHBoxLayout()
        vboxt1 = QtGui.QVBoxLayout()
        vboxt1.addWidget(sizer1)
        vboxt1.addStretch (1)
        vboxt1.addWidget(sizer2)
        
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
    def OnLoadHDF5(self, event):

        wildcard =  'HDF5 files (*.hdf5)'
        self.window().LoadStack(wildcard)
        
#----------------------------------------------------------------------          
    def OnLoadSDF(self, event):

        wildcard =  "SDF files (*.hdr)" 
        self.window().LoadStack(wildcard)
                
#----------------------------------------------------------------------          
    def OnLoadSTK(self, event):

        wildcard =  "STK files (*.stk)" 
        self.window().LoadStack(wildcard)
        
#----------------------------------------------------------------------          
    def OnLoadTXRM(self, event):

        wildcard =  "TXRM (*.txrm)" 
        self.window().LoadStack(wildcard)

#----------------------------------------------------------------------          
    def OnLoadXRM(self, event):

        wildcard =  "XRM (*.xrm)" 
        self.window().LoadStack(wildcard)

#----------------------------------------------------------------------          
    def OnLoadTIF(self, event):

        wildcard =  "TIF (*.tif)" 
        self.window().LoadStack(wildcard)
                
#----------------------------------------------------------------------          
    def OnBuildStack(self, event):

        self.window().BuildStack()
        
#----------------------------------------------------------------------            
    def OnScrollEng(self, value):
        self.iev = value

        if self.com.stack_loaded == 1:
            self.ShowImage()
            
                    
#----------------------------------------------------------------------        
    def ShowImage(self):
        
        image = self.stk.absdata[:,:,int(self.iev)].copy() 

        fig = self.absimgfig
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        axes = fig.gca()
        fig.patch.set_alpha(1.0)
         
        im = axes.imshow(image, cmap=matplotlib.cm.get_cmap("gray")) 
         
        if self.window().page1.show_scale_bar == 1:
            #Show Scale Bar
            startx = int(self.stk.n_rows*0.05)
            starty = self.stk.n_cols-int(self.stk.n_cols*0.05)-self.stk.scale_bar_pixels_y
            um_string = '$\mu m$'
            microns = '$'+self.stk.scale_bar_string+' $'+um_string
            axes.text(self.stk.scale_bar_pixels_x+startx+1,starty+1, microns, horizontalalignment='left', verticalalignment='center',
                      color = 'black', fontsize=14)
            #Matplotlib has flipped scales so I'm using rows instead of cols!
            p = matplotlib.patches.Rectangle((startx,starty), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                                   color = 'black', fill = True)
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

        self.data_struct = data_struct
        self.stk = stack
        self.common = com
        
        self.resize(400, 500)
        self.setWindowTitle('Stack File List')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)
        
        self.filepath = filepath
        
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
        self.filelist.setColumnCount(4)
        self.filelist.setHorizontalHeaderLabels(('File list', 'X', 'Y', 'eV'))
        self.filelist.setShowGrid(False)
        
        self.filelist.setColumnWidth(0,150)
        self.filelist.setColumnWidth(1,50)
        self.filelist.setColumnWidth(2,50)
        self.filelist.setColumnWidth(3,50)
        
        self.filelist.setRowCount(3)

        self.filelist.setItem(1, 0, QtGui.QTableWidgetItem('3'))


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
        #self.button_accept.Disable()
        #self.button_accept.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnAccept, id=self.button_accept.GetId())
        hbox.addWidget(self.button_accept)
        
        button_cancel = QtGui.QPushButton('Cancel')
        #self.Bind(wx.EVT_BUTTON, self.OnCancel, id=button_cancel.GetId())
        #button_cancel.SetFont(self.com.font)
        hbox.addWidget(button_cancel)
        
        vbox.addLayout(hbox)
        
        vbox.addStretch(0.5)
                        
        
        
        self.setLayout(vbox)

        
#----------------------------------------------------------------------           
    def ShowFileList(self):    
        
                
        self.sm_files = [x for x in os.listdir(self.filepath) if x.endswith('.sm')]
        
        if self.sm_files:
            
            self.filetype = 'sm'
            
            try: 
                from netCDF4 import Dataset
                import sm_netcdf
                self.sm = sm_netcdf.sm(data_struct)
            
            except:
                QtGui.QMessageBox.warning(self, 'Error', "Could not import netCDF4 library.")
                return

            count = 0
        
            for i in range(len(self.sm_files)):
                #print sm_files
                filename = self.sm_files[i]
                file = os.path.join(self.filepath, filename)
            
                filever, ncols, nrows, iev = self.sm.read_sm_header(file)
            
                if filever > 0:           
                    self.filelist.InsertStringItem(count,filename)
                    self.filelist.SetStringItem(count,1,str(ncols))
                    self.filelist.SetStringItem(count,2,str(nrows))
                    self.filelist.SetStringItem(count,3,'{0:5.2f}'.format(iev))
                    count += 1
                else:
                    continue
                
            return
         
            
        self.xrm_files = [x for x in os.listdir(self.filepath) if x.endswith('.xrm')] 
        

        if self.xrm_files:        
            
            self.filetype = 'xrm'
            
            import xradia_xrm
            self.xrm = xradia_xrm.xrm()

            count = 0
        
            for i in range(len(self.xrm_files)):

                filename = self.xrm_files[i]
                file = os.path.join(self.filepath, filename)

                ncols, nrows, iev = self.xrm.read_xrm_fileinfo(file)
  
                if ncols > 0:           
                    self.filelist.InsertStringItem(count,filename)
                    self.filelist.SetStringItem(count,1,str(ncols))
                    self.filelist.SetStringItem(count,2,str(nrows))
                    self.filelist.SetStringItem(count,3,'{0:5.2f}'.format(iev))
                    count += 1

            
        self.sm_files = self.xrm_files
        return
        
        
                        
        
#---------------------------------------------------------------------- 
class AboutFrame(QtGui.QDialog):

    def __init__(self, parent = None):    
        QtGui.QWidget.__init__(self, parent)

        self.resize(360, 660)
        self.setWindowTitle('About Mantis')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)
        
        vbox = QtGui.QVBoxLayout()
        
       

        self.image = QtGui.QImage(os.path.join('images','Mantis_logo_about.png'))
        
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
        text2.setText('''Mantis 2.0''')  
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
        self.nnma = nnma.nnma(self.stk, self.data_struct, self.anlz)
        self.common = common()
               

        self.resize(Winsizex, Winsizey)
        self.setWindowTitle('Mantis')
        
        self.initToolbar()

                            
        ico = QtGui.QIcon(os.path.join('images','logo-2l-32.ico'))
        self.setWindowIcon(ico)  
        
        tabs    = QtGui.QTabWidget()
        

  
        self.page2 = PagePCA()
        self.page3 = PageCluster()
        self.page4 = PageSpectral()
        self.page5 = PageKeyEng()
        self.page7 = None
        
        # create the page windows as tabs
        self.page0 = PageLoadData(self.common, self.data_struct, self.stk)
        self.page1 = PageStack(self.common, self.data_struct, self.stk)
#         self.page2 = PagePCA(nb, self.common, self.data_struct, self.stk, self.anlz)
#         self.page3 = PageCluster(nb, self.common, self.data_struct, self.stk, self.anlz)
#         self.page4 = PageSpectral(nb, self.common, self.data_struct, self.stk, self.anlz)
#         self.page7 = None
        
        
        tabs.addTab(self.page0,"Load Data")
        tabs.addTab(self.page1,"Preprocess Data")
        tabs.addTab(self.page2,"PCA")
        tabs.addTab(self.page3,"Cluster Analysis")
        tabs.addTab(self.page4,"Spectral Analysis")
        tabs.addTab(self.page5,"Key Energies")
    
        layout = QVBoxLayout()
#         tempLayout = QHBoxLayout()
        layout.addWidget(tabs)
        self.setCentralWidget(tabs)
#         rightLayout.addLayout(tempLayout)
#         rightLayout.addWidget(QListView())
#         
#         layout.addLayout(leftLayout)
#         layout.addLayout(rightLayout)
        #self.setLayout(layout)
                              
        self.show()
        

#----------------------------------------------------------------------   
    def initToolbar(self):   
        
        self.actionOpen = QtGui.QAction(self)
        self.actionOpen.setObjectName('actionOpen')
        self.actionOpen.setIcon(QtGui.QIcon(os.path.join('images','document-open.png')))
        self.toolbar = self.addToolBar('actionOpen') 
        self.toolbar.addAction(self.actionOpen)
        self.actionOpen.triggered.connect(self.LoadStack)
        
        self.actionOpenSL = QtGui.QAction(self)
        self.actionOpenSL.setObjectName('actionOpenSL')
        self.actionOpenSL.setIcon(QtGui.QIcon(os.path.join('images','open-sl.png')))
        #toolbar = self.addToolBar('actionOpenSL') 
        self.toolbar.addAction(self.actionOpenSL)
        self.actionOpenSL.triggered.connect(self.BuildStack)
        
        self.actionSave = QtGui.QAction(self)
        self.actionSave.setObjectName('actionSave')
        self.actionSave.setIcon(QtGui.QIcon(os.path.join('images','media-floppy.png')))
        #toolbar = self.addToolBar('actionSave') 
        self.toolbar.addAction(self.actionSave)
        self.actionSave.triggered.connect(self.onSaveAsH5)

        self.actionInfo = QtGui.QAction(self)
        self.actionInfo.setObjectName('actionInfo')
        self.actionInfo.setIcon(QtGui.QIcon(os.path.join('images','help-browser.png')))
        #toolbar = self.addToolBar('actionInfo') 
        self.toolbar.addAction(self.actionInfo)
        self.actionInfo.triggered.connect(self.onAbout)
        
#----------------------------------------------------------------------
    def LoadStack(self, wildcard = ''):
        """
        Browse for a stack file:
        """

        if True:
        #try:
            if wildcard == False:
                wildcard =  "HDF5 files (*.hdf5);;SDF files (*.hdr);;STK files (*.stk);;TXRM (*.txrm);;XRM (*.xrm);;TIF (*.tif)" 

            filepath = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '', wildcard)
            

            filepath = str(filepath)
            if filepath == '':
                return
            
            
            directory =  os.path.dirname(str(filepath))
            self.page1.filename =  os.path.basename(str(filepath))
        
            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
            basename, extension = os.path.splitext(self.page1.filename)      
            
            self.common.path = directory
            self.common.filename = self.page1.filename
                       
            
            if extension == '.hdr':            
                if self.common.stack_loaded == 1:
                    self.new_stack_refresh()  
                    self.stk.new_data()
                    #self.stk.data_struct.delete_data()
                    self.anlz.delete_data()    
                self.stk.read_sdf(filepath)        
                           
            
            if extension == '.stk':            
                if self.common.stack_loaded == 1:
                    self.new_stack_refresh()  
                    self.stk.new_data()
                    #self.stk.data_struct.delete_data()
                    self.anlz.delete_data()       
                self.stk.read_stk(filepath)     
                

                
            if extension == '.hdf5':
                if self.common.stack_loaded == 1:
                    self.new_stack_refresh()  
                    self.stk.new_data()
                    #self.stk.data_struct.delete_data()
                    self.anlz.delete_data()                

                self.stk.read_h5(filepath)
                         

                
            if extension == '.txrm':            
                if self.common.stack_loaded == 1:
                    self.new_stack_refresh()  
                    self.stk.new_data()
                    #self.stk.data_struct.delete_data()
                    self.anlz.delete_data()  
                         
                self.stk.read_txrm(filepath)        

                
                              
            if extension == '.xrm':              
                if self.common.stack_loaded == 1:
                    self.new_stack_refresh()  
                    self.stk.new_data()
                    #self.stk.data_struct.delete_data()
                    self.anlz.delete_data()  
                         
                self.stk.read_xrm(filepath)        
                
            if extension == '.tif':              
                if self.common.stack_loaded == 1:
                    self.new_stack_refresh()  
                    self.stk.new_data()
                    #self.stk.data_struct.delete_data()
                    self.anlz.delete_data()  
                         
                self.stk.read_tiff(filepath)    
                self.page1.show_scale_bar = 0
                self.page1.add_scale_cb.SetValue(False)


            #Update widgets 
            self.iev = int(self.stk.n_ev/2)
            self.page0.slider_eng.setRange(0,self.stk.n_ev-1)
            self.page0.iev = self.iev
            self.page0.slider_eng.setValue(self.iev)
                     
            self.page1.slider_eng.setRange(0,self.stk.n_ev-1)
            self.page1.iev = self.iev
            self.page1.slider_eng.setValue(self.iev)
            
        
            x=self.stk.n_cols
            y=self.stk.n_rows
            z=self.iev               
            self.page1.imgrgb = npy.zeros(x*y*3,dtype = "uint8")        
            self.page1.maxval = npy.amax(self.stk.absdata)
            

            self.ix = int(x/2)
            self.iy = int(y/2)
                    
            self.common.stack_loaded = 1
            
            if self.stk.data_struct.spectromicroscopy.normalization.white_spectrum is not None:
                self.common.i0_loaded = 1
            
            self.page1.ix = self.ix
            self.page1.iy = self.iy
            
            self.page1.ResetDisplaySettings()
            self.page1.loadImage()
            self.page1.loadSpectrum(self.ix, self.iy)
            self.page1.textctrl.setText(self.page1.filename)
            
            self.page0.ShowInfo(self.page1.filename, directory)
            

            QtGui.QApplication.restoreOverrideCursor()
                
#         except:
#  
#             self.common.stack_loaded = 0 
#             self.common.i0_loaded = 0
#             self.new_stack_refresh()
#                                 
#             QtGui.QApplication.restoreOverrideCursor()
#             QtGui.QMessageBox.warning(self, 'Error', 'Image stack not loaded.')
# 
#             import sys; print sys.exc_info()
                   

        #self.refresh_widgets()


#----------------------------------------------------------------------
    def BuildStack(self):
        """
        Browse for .sm files
        """
        
        if True:
        #try:
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
#             wx.MessageBox(".sm files not loaded.")
#             import sys; print sys.exc_info()
            
#----------------------------------------------------------------------
    def onSaveAsH5(self, event):
        self.SaveStackH5()
        
        
#----------------------------------------------------------------------
    def SaveStackH5(self):

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
            
                
            self.stk.write_h5(filepath, self.data_struct)    
            QtGui.QApplication.restoreOverrideCursor()

        except:
 
            QtGui.QApplication.restoreOverrideCursor()
            
            QtGui.QMessageBox.warning(self, 'Error', 'Could not save HDF5 file.')
                   

        #self.refresh_widgets()
        
        return
       
#----------------------------------------------------------------------
    def onAbout(self):

        self.popup = AboutFrame(self)
        self.popup.show()


#----------------------------------------------------------------------        
    def new_stack_refresh(self):
        
        
        self.common.i0_loaded = 0
        self.common.pca_calculated = 0
        self.common.cluster_calculated = 0
        
        #self.refresh_widgets()
        
#         #page 1
#         self.page1.rb_flux.SetValue(True)
#         self.page1.rb_od.SetValue(False)
#         self.page1.showflux = True
#         
#         fig = self.page1.SpectrumPanel.get_figure()
#         fig.clf()
#         self.page1.SpectrumPanel.draw()
#         self.page1.tc_spec.SetValue("Spectrum at point: ")       
#         
#         fig = self.page1.AbsImagePanel.get_figure()
#         fig.clf()
#         self.page1.AbsImagePanel.draw()        
#         self.page1.tc_imageeng.SetValue("Image at energy: ")
#         
#         self.page1.textctrl.SetValue(' ')
#         
#         self.page1.ResetDisplaySettings()
#         
#         
#         #page 2
#         fig = self.page2.PCAEvalsPan.get_figure()
#         fig.clf()
#         self.page2.PCAEvalsPan.draw()
#         
#         fig = self.page2.PCAImagePan.get_figure()
#         fig.clf()
#         self.page2.PCAImagePan.draw()
#         
#         fig = self.page2.PCASpecPan.get_figure()
#         fig.clf()
#         self.page2.PCASpecPan.draw()
#         
#         self.page2.vartc.SetLabel('0%')
#         self.page2.npcaspin.SetValue(1)
#         self.page2.tc_PCAcomp.SetValue("PCA component ")
#         self.page2.text_pcaspec.SetValue("PCA spectrum ")
#         
#         self.page2.selpca = 1       
#         self.page2.numsigpca = 2
#         self.page2.slidershow.SetValue(self.page2.selpca)
#         
#         #page 3
#         fig = self.page3.ClusterImagePan.get_figure()
#         fig.clf()
#         self.page3.ClusterImagePan.draw()
#         
#         fig = self.page3.ClusterIndvImagePan.get_figure()
#         fig.clf()
#         self.page3.ClusterIndvImagePan.draw()
#         
#         fig = self.page3.ClusterSpecPan.get_figure()
#         fig.clf()
#         self.page3.ClusterSpecPan.draw()
#         
#         fig = self.page3.ClusterDistMapPan.get_figure()
#         fig.clf()
#         self.page3.ClusterDistMapPan.draw()       
#        
#         self.page3.selcluster = 1
#         self.page3.slidershow.SetValue(self.page3.selcluster)
#         self.page3.numclusters = 5
#         self.page3.nclusterspin.SetValue(self.page3.numclusters)
#         self.page3.tc_cluster.SetValue("Cluster ")
#         self.page3.tc_clustersp.SetValue("Cluster spectrum")
#         self.page3.wo_1st_pca = 0
#         self.page3.remove1stpcacb.SetValue(False)
# 
#         
#         #page 4
#         self.page4.ClearWidgets()
#         
#         #page 7
#         if self.page7:
#             self.page7.button_calckeng.Disable()
#             fig = self.page7.KESpecPan.get_figure()
#             fig.clf()
#             self.page7.KESpecPan.draw()         
#             fig = self.page7.AbsImagePanel.get_figure()
#             fig.clf()
#             self.page7.AbsImagePanel.draw()   
#             self.page7.lc_1.DeleteAllItems()   
#             self.page7.keyenergies = []
#             self.page7.keyengs_calculated = 0    
#                
""" ------------------------------------------------------------------------------------------------"""
                        
def main():
    
    app = QtGui.QApplication(sys.argv)
    frame = MainFrame()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()