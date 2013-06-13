'''
Created on Jun 11, 2013

@author: Mirna Lerotic
'''

import sys
from PyQt4 import QtCore, QtGui
import os
from PyQt4.QtGui import *

import logos

Winsizex = 1000
Winsizey = 740

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
        pass


""" ------------------------------------------------------------------------------------------------"""
class PageStack(QtGui.QWidget):
    def __init__(self):
        super(PageStack, self).__init__()

        self.initUI()
        
#----------------------------------------------------------------------          
    def initUI(self): 
        
        self.brightness_min = 0.0
        self.brightness_max = 1.0
        self.gamma = 1.0
        

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
                
    
        vbox44 = QtGui.QVBoxLayout()
        self.button_despike = QtGui.QPushButton('Despike')
        #self.button_despike.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnDespike, id=self.button_despike.GetId())   
        #self.button_despike.Disable()     
        vbox44.addWidget(self.button_despike)
        self.button_resetdisplay = QtGui.QPushButton( 'Reset')
        #self.button_resetdisplay.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.onResetDisplaySettings, id=self.button_resetdisplay.GetId())   
        #self.button_resetdisplay.Disable()     
        vbox44.addWidget(self.button_resetdisplay)
        self.button_displaycolor = QtGui.QPushButton('Color Table...   ')
        #self.button_displaycolor.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.onSetColorTable, id=self.button_displaycolor.GetId())   
        #self.button_displaycolor.Disable()     
        vbox44.addWidget(self.button_displaycolor)
        
        hbox23.addSpacing(20)
        hbox23.addLayout(vbox44)
        
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
        
        
        vboxtop = QtGui.QVBoxLayout()
        
        hboxtop = QtGui.QHBoxLayout()
        hboxtop.addWidget(sizer1)
        hboxtop.addWidget(sizer2)
        hboxtop.addWidget(sizer3)
        
#         hboxtop.addStretch (0.5)
#         hboxtop.addLayout(vboxt1)
#         hboxtop.addStretch (0.5)
#         hboxtop.addLayout(vbox5)
#         hboxtop.addStretch (0.5)
#         
        vboxtop.addStretch (1)
        vboxtop.addLayout(hboxtop)
        vboxtop.addStretch (1)
#         vboxtop.addWidget(sizer3)
#         vboxtop.addStretch (0.5)
#         vboxtop.addWidget(sizer4)
#         vboxtop.addStretch (0.5)

        vboxtop.setContentsMargins(20,20,20,20)
        self.setLayout(vboxtop)
        

""" ------------------------------------------------------------------------------------------------"""
class PageLoadData(QtGui.QWidget):
    def __init__(self):
        super(PageLoadData, self).__init__()

        self.initUI()
        
#----------------------------------------------------------------------          
    def initUI(self): 
        

        #panel 1
        sizer1 = QtGui.QGroupBox('Load Data Stack')
        vbox1 = QtGui.QVBoxLayout()
        
        self.button_hdf5 = QtGui.QPushButton('Load HDF5 Stack (*.hdf5)')
        #self.button_hdf5.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnLoadHDF5, id=self.button_hdf5.GetId())
        vbox1.addWidget(self.button_hdf5)
        
        self.button_sdf = QtGui.QPushButton('Load SDF Stack (*.hrd)')
        #self.button_sdf.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnLoadSDF, id=self.button_sdf.GetId())   
        vbox1.addWidget(self.button_sdf)
        
        self.button_stk = QtGui.QPushButton('Load STK Stack (*.stk)')
        #self.button_stk.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnLoadSTK, id=self.button_stk.GetId())   
        vbox1.addWidget(self.button_stk)
        
        self.button_xrm = QtGui.QPushButton('Load XRM Image (*.xrm)')
        #self.Bind(wx.EVT_BUTTON, self.OnLoadXRM, id=self.button_xrm.GetId())
        #self.button_xrm.SetFont(self.com.font)
        vbox1.addWidget(self.button_xrm)
        
        self.button_txrm = QtGui.QPushButton('Load TXRM Stack (*.txrm)')
        #self.Bind(wx.EVT_BUTTON, self.OnLoadTXRM, id=self.button_txrm.GetId())
        #self.button_txrm.SetFont(self.com.font)
        vbox1.addWidget(self.button_txrm)    
        
        self.button_tif = QtGui.QPushButton( 'Load TIF Stack (*.tif)')
        #self.Bind(wx.EVT_BUTTON, self.OnLoadTIF, id=self.button_tif.GetId())
        #self.button_tif.SetFont(self.com.font)
        vbox1.addWidget(self.button_tif)      

        #vbox1.setContentsMargins(15,15,15,15)
        sizer1.setLayout(vbox1)
                
        
        

        #panel 2
        sizer2 = QtGui.QGroupBox('Build a stack from a set of files')
        vbox2 = QtGui.QVBoxLayout()

        self.button_sm = QtGui.QPushButton( 'Build a stack from a set of NetCDF (*.sm) files')
        #self.button_sm.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnBuildStack, id=self.button_sm.GetId())
        vbox2.addWidget(self.button_sm)
        
        self.button_xrm_list = QtGui.QPushButton( 'Build a stack from a set of XRM (*.xrm) files')
        #self.button_xrm_list.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.OnBuildStack, id=self.button_xrm_list.GetId())
        vbox2.addWidget(self.button_xrm_list)
        
        #vbox2.setContentsMargins(15,15,15,15)
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
        #self.tc_imageeng.SetFont(self.com.font)
        self.tc_imageeng.setText("Image at energy: ")
        vbox5.addWidget(self.tc_imageeng)
        
        
        hbox51 = QtGui.QHBoxLayout()
   
        #i1panel = wx.Panel(panel5, -1, style = wx.SUNKEN_BORDER)
        #self.AbsImagePanel = wxmpl.PlotPanel(i1panel, -1, size =(PlotH*.8, PlotH*.8), cursor=False, crosshairs=False, location=False, zoom=False)
        self.image = QtGui.QImage(os.path.join('images','Mantis_img.jpg'))
        
        self.imageLabel = QtGui.QLabel()
        self.imageLabel.setBackgroundRole(QtGui.QPalette.Base)
        #self.imageLabel.setSizePolicy()
        #self.imageLabel.setScaledContents(True)
        
        self.imageLabel.setPixmap(QtGui.QPixmap.fromImage(self.image))
        hbox51.addWidget(self.imageLabel)
        

        self.slider_eng = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slider_eng.setFocusPolicy(QtCore.Qt.StrongFocus)
#         vbox51 = wx.BoxSizer(wx.VERTICAL)
#         self.slider_eng = wx.Slider(panel5, 0, self.iev, 0, 100, style=wx.SL_LEFT|wx.SL_VERTICAL)        
#         self.slider_eng.SetFocus()
#         self.Bind(wx.EVT_SCROLL, self.OnScrollEng, self.slider_eng)
        hbox51.addWidget(self.slider_eng)
# 
#         self.engspin = wx.SpinButton(panel5, -1, size = ((8,-1)), style=wx.SP_ARROW_KEYS)
#         self.Bind(wx.EVT_SPIN_UP, self.OnEngspinUp, self.engspin)
#         self.Bind(wx.EVT_SPIN_DOWN, self.OnEngspinDown, self.engspin)
        

        vbox5.addLayout(hbox51)
        
#         hbox51.Add(i1panel, 0)
#         vbox51.Add((0,3))
#         vbox51.Add(self.slider_eng, 1,  wx.EXPAND) 
#         vbox51.Add(self.engspin, 0,  wx.EXPAND)         
#         hbox51.Add(vbox51, 0,  wx.EXPAND) 
#         
#         vbox5.Add(self.tc_imageeng, 0)        
#         vbox5.Add(hbox51, 0)

 
       
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

       
""" ------------------------------------------------------------------------------------------------"""
class MainFrame(QtGui.QMainWindow):
    
    def __init__(self):
        super(MainFrame, self).__init__()
        
        self.initUI()

#----------------------------------------------------------------------          
    def initUI(self):    
               
        #self.setGeometry(300, 300, 250, 150)
        self.resize(Winsizex, Winsizey)
        self.setWindowTitle('Mantis')
        
        self.initToolbar()

                            
        ico = QtGui.QIcon(os.path.join('images','logo-2l-32.ico'))
        self.setWindowIcon(ico)  
        
        tabs    = QtGui.QTabWidget()
        pushButton1 = QtGui.QPushButton("QPushButton 1")
        pushButton2 = QtGui.QPushButton("QPushButton 2")
        

        tab1 = PageLoadData()   
        tab2 = PageStack()   
        tab3 = PagePCA()
        tab4 = PageCluster()
        tab5 = PageSpectral()
        tab6 = PageKeyEng()
        

        
        tabs.addTab(tab1,"Load Data")
        tabs.addTab(tab2,"Preprocess Data")
        tabs.addTab(tab3,"PCA")
        tabs.addTab(tab4,"Cluster Analysis")
        tabs.addTab(tab5,"Spectral Analysis")
        tabs.addTab(tab6,"Key Energies")
    
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
        
        self.actionOpenSL = QtGui.QAction(self)
        self.actionOpenSL.setObjectName('actionOpenSL')
        self.actionOpenSL.setIcon(QtGui.QIcon(os.path.join('images','open-sl.png')))
        #toolbar = self.addToolBar('actionOpenSL') 
        self.toolbar.addAction(self.actionOpenSL)
        
        self.actionSave = QtGui.QAction(self)
        self.actionSave.setObjectName('actionSave')
        self.actionSave.setIcon(QtGui.QIcon(os.path.join('images','media-floppy.png')))
        #toolbar = self.addToolBar('actionSave') 
        self.toolbar.addAction(self.actionSave)

        self.actionInfo = QtGui.QAction(self)
        self.actionInfo.setObjectName('actionInfo')
        self.actionInfo.setIcon(QtGui.QIcon(os.path.join('images','help-browser.png')))
        #toolbar = self.addToolBar('actionInfo') 
        self.toolbar.addAction(self.actionInfo)
        
                
def main():
    
    app = QtGui.QApplication(sys.argv)
    frame = MainFrame()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()