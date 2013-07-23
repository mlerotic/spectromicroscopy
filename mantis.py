'''
Created on Jun 11, 2013

@author: Mirna Lerotic
'''

import sys
import os
import numpy as npy
import time

from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import *
from PyQt4.QtCore import Qt, QCoreApplication


import matplotlib 
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid import make_axes_locatable
matplotlib.interactive( True )

import data_struct
import data_stack
import analyze
import nnma
import henke


Winsizex = 1000
Winsizey = 700

PlotH = 4.0
PlotW = PlotH*1.61803

ImgDpi = 40

#----------------------------------------------------------------------
class common:
    def __init__(self):
        
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
    def __init__(self, common, data_struct, stack, anlz):
        super(PageKeyEng, self).__init__()

        self.initUI(common, data_struct, stack, anlz)
        
#----------------------------------------------------------------------          
    def initUI(self, common, data_struct, stack, anlz): 
        
        
        self.com = common 
        self.data_struct = data_struct
        self.stk = stack       
        self.anlz = anlz
        
        self.selica = 1       
        self.numica = 2
        
        
        pw = PlotW*0.8
        ph = PlotH*0.8
        
        
        self.i_eng = 0
        self.keyengs_calculated = 0
       
         
          
        #panel 1        
        vbox1 = QtGui.QVBoxLayout()
               
        self.tc_1 = QtGui.QLabel(self)
        self.tc_1.setText("Average Optical Density")    
        
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
   
        self.kespecfig = Figure((pw, ph))
        self.KESpecPan = FigureCanvas(self.kespecfig)
        self.KESpecPan.mpl_connect('button_press_event', self.OnPointSpectrum)
        fbox.addWidget(self.KESpecPan)
        frame.setLayout(fbox)
            
        vbox1.addStretch(1)
        vbox1.addWidget(self.tc_1)        
        vbox1.addWidget(frame)
        vbox1.addStretch(1)
 
                    
        #panel 2
        sizer2 = QtGui.QGroupBox('Key Energies Analysis')
        vbox2 = QtGui.QVBoxLayout()
                
        self.button_calckeng = QtGui.QPushButton('Find Key Energies')
        self.button_calckeng.clicked.connect(self.OnCalcKeyEng)     
        self.button_calckeng.setEnabled(False)  
        vbox2.addWidget( self.button_calckeng) 
        self.button_save = QtGui.QPushButton('Save Results...')
        self.button_save.clicked.connect(self.OnSave)
        self.button_save.setEnabled(False)
        
        
        hbox21 = QtGui.QHBoxLayout()
        text1 = QtGui.QLabel(self)
        text1.setText("Threshold")
        
    
        self.tc_keyengthresh = QtGui.QDoubleSpinBox()
        self.tc_keyengthresh.setRange(0,5)
        self.tc_keyengthresh.setValue(0.1) 
        self.tc_keyengthresh.setSingleStep(0.1)   


        hbox21.addWidget(text1)
        hbox21.addWidget(self.tc_keyengthresh)    
        
        line = QtGui.QFrame()
        line.setFrameShape(QtGui.QFrame.HLine)
        line.setFrameShadow(QtGui.QFrame.Sunken) 
              
        vbox2.addLayout(hbox21)    
        vbox2.addWidget(line)  
        vbox2.addWidget( self.button_save)        
        

        sizer2.setLayout(vbox2)


        #panel 3
        vbox3 = QtGui.QVBoxLayout()
    
        t1 = QtGui.QLabel(self)
        t1.setText("Key Energies")       
         
        self.lc_1 = QListWidget()   
        self.lc_1.itemClicked.connect(self.OnEngListClick)
        self.lc_1.setMinimumSize(200, 400)
         
        vbox3.addWidget(t1)
        vbox3.addWidget(self.lc_1)
                 
      
        #panel 4     
        vbox4 = QtGui.QVBoxLayout()
         
        self.tc_imageeng = QtGui.QLabel(self)
        self.tc_imageeng.setText("Image at key energy: ")
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
        self.slider_eng.setRange(0, 100)
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


        
        
#----------------------------------------------------------------------
    def OnCalcKeyEng(self, event):
        

        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        
        threshold = self.tc_keyengthresh.value()
  
        self.keyenergies = []
        
        self.keyenergies = self.anlz.calc_key_engs(threshold)
        
        if len(self.keyenergies) > 0:        
            self.keyengs_calculated = 1
        
            self.i_eng = 0
            self.slider_eng.setRange(0,len(self.keyenergies)-1)
            self.slider_eng.setValue(self.i_eng)

            self.ShowPlots()
            self.ShowImage()
            self.ShowKEngs()
            
            self.button_save.setEnabled(True)
            
        else:
            self.ShowPlots()
            fig = self.absimgfig
            fig.clf()
            self.AbsImagePanel.draw()               
            self.lc_1.clear()
            
            self.button_save.setEnabled(False)
            
        QtGui.QApplication.restoreOverrideCursor() 

        
#----------------------------------------------------------------------            
    def OnScrollEng(self, value):
        self.i_eng = value

        if self.keyengs_calculated == 1:
            self.ShowImage()
            self.ShowKEngs()
            
#----------------------------------------------------------------------            
    def OnEngspinUp(self, event):
        if (self.keyengs_calculated == 1) and (self.i_eng > 0):
            self.i_eng = self.i_eng - 1
            self.slider_eng.setValue(self.i_eng)

            self.ShowImage()
            self.ShowKEngs()

            
#----------------------------------------------------------------------            
    def OnEngspinDown(self, event):
        if (self.keyengs_calculated == 1) and (self.i_eng < len(self.keyenergies)-1):
            self.i_eng = self.i_eng + 1
            self.slider_eng.setValue(self.i_eng)

            self.ShowImage()
            self.ShowKEngs()

#----------------------------------------------------------------------  
    def OnPointSpectrum(self, evt):
        x = evt.xdata
        
        if (self.keyengs_calculated == 1):      
            if x < self.stk.ev[0]:
                sel_ev = 0
            elif x > self.stk.ev[self.stk.n_ev-1]:
                sel_ev = self.stk.n_ev-1
            else:
                indx = npy.abs(self.stk.ev - x).argmin()
                sel_ev = indx
                
            
            self.i_eng=(npy.abs(self.keyenergies-self.stk.ev[sel_ev])).argmin()                 

            self.ShowImage()
            self.ShowKEngs()
            
            self.slider_eng.setValue(self.i_eng)
                              
#----------------------------------------------------------------------     
    def ShowPlots(self):


        odtotal = self.stk.od3d.sum(axis=0)   
        odtotal = odtotal.sum(axis=0)/(self.stk.n_rows*self.stk.n_cols)      
        odtotal /= odtotal.max()/0.7

        
        fig = self.kespecfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        

        
        specplot = axes.plot(self.stk.ev,odtotal)
#        for i in range(self.anlz.numsigpca):
#            pcaspectrum = self.anlz.eigenvecs[:,i]
#            specplot = axes.plot(self.stk.ev,pcaspectrum)


        for i in range(len(self.keyenergies)):
            axes.axvline(x=self.keyenergies[i], color = 'g', alpha=0.5)
                        
        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Optical Density')
        
        self.KESpecPan.draw()

#----------------------------------------------------------------------        
    def OnEngListClick(self):
        item = self.lc_1.currentRow()
                
        self.i_eng = item
        
        self.ShowKEngs()     
        self.ShowImage()
         
#----------------------------------------------------------------------           
    def ShowKEngs(self):    
        
        self.lc_1.clear()
        
        for i in range(len(self.keyenergies)):
            self.lc_1.addItem('{0:08.2f}'.format(self.keyenergies[i]))
            

#----------------------------------------------------------------------        
    def ShowImage(self):
        
        iev=(npy.abs(self.stk.ev-self.keyenergies[self.i_eng])).argmin() 
        image = self.stk.absdata[:,:,int(iev)].copy() 

        fig = self.absimgfig
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        axes = fig.gca()
        fig.patch.set_alpha(1.0)
        
        im = axes.imshow(image, cmap=matplotlib.cm.get_cmap("gray")) 
        
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
        
        self.tc_imageeng.setText('Image at key energy: {0:5.2f} eV'.format(float(self.stk.ev[iev])))
 
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
            
            
        #Save text file with list of energies
        textfilepath = path+'_keyenergies.csv'
        f = open(textfilepath, 'w')
        print>>f, '*********************  Key Energies  ********************'
        for i in range(len(self.keyenergies)):
            print>>f, '%.6f' %(self.keyenergies[i])
        
        f.close()        
            
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
         
         
         
        #panel 3
        sizer3 = QtGui.QGroupBox('Target Spectrum')
        vbox3 = QtGui.QVBoxLayout()
        
         
        self.button_loadtspec = QtGui.QPushButton('Load Spectrum')
        self.button_loadtspec.clicked.connect(self.OnTSpecFromFile)
        self.button_loadtspec.setEnabled(False)
        vbox3.addWidget(self.button_loadtspec)
        self.button_addflat = QtGui.QPushButton('Add Flat Spectrum')
        self.button_addflat.clicked.connect( self.OnFlatTSpec)
        self.button_addflat.setEnabled(False)
        vbox3.addWidget(self.button_addflat)
        self.button_addclspec = QtGui.QPushButton('Add Cluster Spectra')
        self.button_addclspec.clicked.connect( self.OnAddClusterSpectra)   
        self.button_addclspec.setEnabled(False)
        vbox3.addWidget(self.button_addclspec)
         
        self.button_showrgb = QtGui.QPushButton('Composite RGB image...')
        self.button_showrgb.clicked.connect( self.OnCompositeRGB)   
        self.button_showrgb.setEnabled(False)
        vbox3.addWidget(self.button_showrgb)        
 
        self.button_save = QtGui.QPushButton('Save Images...')
        self.button_save.clicked.connect( self.OnSave)
        self.button_save.setEnabled(False)
        vbox3.addWidget(self.button_save)

        sizer3.setLayout(vbox3)
         
 
         
        #panel 4
        sizer4 = QtGui.QGroupBox('Display')
        vbox4 = QtGui.QVBoxLayout()
                

        sb = QtGui.QGroupBox('Spectrum')
        vbox41 = QtGui.QVBoxLayout()

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


        self.add_scale_cb = QtGui.QCheckBox('Scale', self) 
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
        

        #try: 
        if True:
            
            wildcard = "Spectrum files (*.csv)"
            
            filepath = QtGui.QFileDialog.getOpenFileName(self, 'Choose Spectrum file', '', wildcard)
            

            filepath = str(filepath)
            if filepath == '':
                return
            
            self.filename =  os.path.basename(str(filepath))
            
                                                        
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
            
#         except:
#             QtGui.QApplication.restoreOverrideCursor()  
#             QtGui.QMessageBox.warning(self, 'Error', 'Spectrum file not loaded.')
                                   
                                 
        self.window().refresh_widgets()
        

#----------------------------------------------------------------------
    def OnFlatTSpec(self, event):

        try: 
            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor)) 
            self.anlz.read_target_spectrum(flat=True)
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
            QtGui.QMessageBox.warning(self, 'Error', 'Flat spectrum not loaded.')
            
                                                      
        self.window().refresh_widgets()
        
#----------------------------------------------------------------------
    def OnAddClusterSpectra(self, event):

        #try:
        if True: 
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
            
#        except:
#            QtGui.QApplication.restoreOverrideCursor()  
#            wx.MessageBox("Cluster spectra not loaded.")
 
                                                        
        self.window().refresh_widgets()
        
#----------------------------------------------------------------------
    def OnCompositeRGB(self, event):

        compimgwin = ShowCompositeRBGmap(self.window(), self.com, self.anlz)
        compimgwin.show()
                
#----------------------------------------------------------------------
    def OnSave(self, event):
        
        pass
        #SaveWinP4().Show()
        
        
#----------------------------------------------------------------------
    def Save(self, filename, path, spec_png = True, spec_pdf = False, spec_csv = False, img_png = True, img_pdf = False):

        self.SaveFileName = os.path.join(path,filename)
   
        try: 
            if img_png:
                self.SaveMaps(png_pdf=1)
            if img_pdf:
                self.SaveMaps(png_pdf=2)
                
            if spec_png:    
                self.SaveSpectra(png_pdf=1)
            if spec_pdf:
                self.SaveSpectra(png_pdf=2)
            if spec_csv:
                self.SaveSpectra(savecsv = True)
                
                
            
        except IOError, e:
            if e.strerror:
                err = e.strerror 
            else: 
                err = e 
            QtGui.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err)

            
#----------------------------------------------------------------------
    def SaveSpectra(self, png_pdf=1, savecsv = False):
        
        
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas   
        matplotlib.rcParams['pdf.fonttype'] = 42
            
        colors=['#FF0000','#000000','#FFFFFF']
        spanclrmap=matplotlib.colors.LinearSegmentedColormap.from_list('spancm',colors)
        
        if png_pdf == 1:   
            ext = 'png'
        else:
            ext = 'pdf'
        suffix = "." + ext
        
        
        for i in range (self.anlz.n_target_spectra):
            #Save spectra images
            tspectrum = self.anlz.target_spectra[i, :]
                        
        
            fig = matplotlib.figure.Figure(figsize =(PlotW, PlotH))
            canvas = FigureCanvas(fig)
            fig.clf()
            fig.add_axes((0.15,0.15,0.75,0.75))
            axes = fig.gca()
        
            matplotlib.rcParams['font.size'] = self.fontsize

            line1 = axes.plot(self.stk.ev,tspectrum, color='black', label = 'Raw data')
            
            if self.com.pca_calculated == 1: 
                tspectrumfit = self.anlz.target_pcafit_spectra[i, :]
                diff = npy.abs(tspectrum-tspectrumfit)
                line2 = axes.plot(self.stk.ev,tspectrumfit, color='green', label = 'Fit')
                line3 = axes.plot(self.stk.ev,diff, color='grey', label = 'Abs(Raw-Fit)')
            
            fontP = matplotlib.font_manager.FontProperties()
            fontP.set_size('small')
       
            axes.legend(loc=4, prop = fontP)
                        
            axes.set_xlabel('Photon Energy [eV]')
            axes.set_ylabel('Optical Density')

            fileName_spec = self.SaveFileName+"_Tspectrum_" +str(i+1)+"."+ext
            fig.savefig(fileName_spec)    
            
            if savecsv:
                fileName_spec = self.SaveFileName+"_Tspectrum_" +str(i+1)+".csv"
                self.stk.write_csv(fileName_spec, self.stk.ev, tspectrum)
#----------------------------------------------------------------------
    def SaveMaps(self, png_pdf=1):            
            
            
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas  
        matplotlib.rcParams['pdf.fonttype'] = 42 
            
        colors=['#FF0000','#000000','#FFFFFF']
        spanclrmap=matplotlib.colors.LinearSegmentedColormap.from_list('spancm',colors)
            
        if png_pdf == 1:   
            ext = 'png'
        else:
            ext = 'pdf'
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
        
        
            min_val = npy.min(tsmapimage)
            max_val = npy.max(tsmapimage)
            bound = npy.max((npy.abs(min_val), npy.abs(max_val)))
        
            if self.show_scale_bar == 1:
                um_string = '$\mu m$'
                microns = '$'+self.stk.scale_bar_string+' $'+um_string
                axes.text(self.stk.scale_bar_pixels_x+10,self.stk.n_cols-9, microns, horizontalalignment='left', verticalalignment='center',
                              color = 'white', fontsize=14)
                #Matplotlib has flipped scales so I'm using rows instead of cols!
                p = matplotlib.patches.Rectangle((5,self.stk.n_cols-10), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                                            color = 'white', fill = True)
                axes.add_patch(p)     
     
            im = axes.imshow(tsmapimage, cmap=spanclrmap, vmin = -bound, vmax = bound)
            cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)  
    
            axes.axis("off") 
                
                   
            fileName_img = self.SaveFileName+"_TSmap_" +str(i+1)+"."+ext               
            fig.savefig(fileName_img, bbox_inches='tight', pad_inches = 0.0)
            
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
    def OnTspecSpinUp(self, event):
        if (self.com.spec_anl_calculated == 1) and (self.i_tspec > 1):
            self.i_tspec = self.i_tspec - 1
            self.slider_tspec.setValue(self.i_tspec)

            self.loadTSpectrum()
            self.loadTargetMap()
            
#----------------------------------------------------------------------            
    def OnTspecSpinDown(self, event):
        if (self.com.spec_anl_calculated == 1) and (self.i_tspec < self.anlz.n_target_spectra):
            self.i_tspec = self.i_tspec + 1
            self.slider_tspec.setValue(self.i_tspec) 
            
            self.loadTSpectrum()
            self.loadTargetMap()


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
            self.loadTSpectrum()
            self.loadTargetMap()  
            self.ShowSpectraList()  
        else:
            self.com.spec_anl_calculated = 0
            self.ClearWidgets()
        
        QtGui.QApplication.restoreOverrideCursor()
        
#----------------------------------------------------------------------           
    def OnMoveSpectrumDown(self, event):
        
        if self.i_tspec < self.anlz.n_target_spectra:
            self.anlz.move_spectrum(self.i_tspec-1, self.i_tspec)
            
            self.i_tspec += 1
            self.slider_tspec.setValue(self.i_tspec)
            self.loadTSpectrum()
            self.loadTargetMap()  
            self.ShowSpectraList()
        
        
#----------------------------------------------------------------------           
    def OnMoveSpectrumUp(self, event):        
        
        if self.i_tspec > 1:
            self.anlz.move_spectrum(self.i_tspec-1, self.i_tspec-2)      
            
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
        self.i_tspec = 1
        self.showraw = True
        self.rb_raw.setChecked(True)
        
        self.slider_tspec.setValue(self.i_tspec)
        
        
        self.textctrl_sp1.setText('Common Name: \n')
        self.textctrl_sp2.setText('RMS Error: ')
        
        self.window().refresh_widgets()
            
#----------------------------------------------------------------------           
    def ShowFitWeights(self):    
        
        self.tc_spfitlist.clear()     
        
        norm_factor = 100./npy.sum(npy.absolute(self.anlz.target_pcafit_coeffs[self.i_tspec-1, :]))
        
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
        
        
        min_val = npy.min(tsmapimage)
        max_val = npy.max(tsmapimage)
        bound = npy.max((npy.abs(min_val), npy.abs(max_val)))
        
        if self.show_scale_bar == 1:
            um_string = '$\mu m$'
            microns = '$'+self.stk.scale_bar_string+' $'+um_string
            axes.text(self.stk.scale_bar_pixels_x+10,self.stk.n_cols-9, microns, horizontalalignment='left', verticalalignment='center',
                      color = 'white', fontsize=14)
            #Matplotlib has flipped scales so I'm using rows instead of cols!
            p = matplotlib.patches.Rectangle((5,self.stk.n_cols-10), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                                   color = 'white', fill = True)
            axes.add_patch(p)     
     
        im = axes.imshow(tsmapimage, cmap=spanclrmap, vmin = -bound, vmax = bound)
        cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)  
    
        axes.axis("off") 
        self.MapPanel.draw()

        
        
#----------------------------------------------------------------------     
    def loadTSpectrum(self):

        tspectrum = self.anlz.target_spectra[self.i_tspec-1, :]
                 
        
        fig = self.TSpecfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        

        line1 = axes.plot(self.stk.ev,tspectrum, color='black', label = 'Raw data')

        if self.com.pca_calculated == 1:       
            tspectrumfit = self.anlz.target_pcafit_spectra[self.i_tspec-1, :]  
            diff = npy.abs(tspectrum-tspectrumfit)      
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
                                    self.anlz.tspec_names[self.i_tspec-1]+'\n')
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
         
        self.rgbimage = npy.zeros((self.n_cols, self.n_rows, 3), dtype=float)
        
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
            tsmap = npy.zeros((self.n_cols, self.n_rows), dtype=float)
        else:
            tsmap = tsmap.clip(min=scale_min, max=scale_max)
            tsmap = (tsmap -scale_min) / (scale_max - scale_min)
            
        indices = npy.where(tsmap < 0)
        tsmap[indices] = 0.0
        indices = npy.where(tsmap > 1)
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
            tsmap = npy.zeros((self.n_cols, self.n_rows), dtype=float)
        else:
            tsmap = tsmap.clip(min=scale_min, max=scale_max)
            tsmap = (tsmap - scale_min) / (scale_max - scale_min)


        indices = npy.where(tsmap < 0)
        tsmap[indices] = 0.0
        indices = npy.where(tsmap > 1)
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
            tsmap = npy.zeros((self.n_cols, self.n_rows), dtype=float)
        else:
            tsmap = tsmap.clip(min=scale_min, max=scale_max)
            tsmap = (tsmap - scale_min) / (scale_max - scale_min)

        indices = npy.where(tsmap < 0)
        tsmap[indices] = 0.0
        indices = npy.where(tsmap > 1)
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
        fig.add_axes((0.02,0.02,0.96,0.96))
        
        
        axes = fig.gca()
        fig.patch.set_alpha(1.0) 
      
        im = axes.imshow(self.rgbimage) 
        
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
    def OnClose(self, evt):
        self.Destroy()             
                
    
    
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
        
        
        line = QtGui.QFrame()
        line.setFrameShape(QtGui.QFrame.HLine)
        line.setFrameShadow(QtGui.QFrame.Sunken) 
        
        
        vbox1.addLayout(hbox14) 
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
        text3.setText('Cluster Distance Map')
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
        
        vboxtopL = QtGui.QVBoxLayout()
        vboxtopL.addWidget(sizer1)
        vboxtopL.addWidget(sizer5)
        
        gridsizertop = QtGui.QGridLayout()
        
        gridsizertop.addLayout(vboxtopL, 0, 0, QtCore .Qt. AlignLeft)
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
        colors_i = npy.linspace(0,self.maxclcolors,self.maxclcolors+1)

       
        self.colors=['#0000FF','#FF0000','#DFE32D','#36F200','#B366FF',
                '#FF470A','#33FFFF','#006600','#CCCC99','#993300',
                '#000000']


        
        self.clusterclrmap1=matplotlib.colors.LinearSegmentedColormap.from_list('clustercm',self.colors)
     
        self.bnorm1 = matplotlib.colors.BoundaryNorm(colors_i, self.clusterclrmap1.N)
        
        colors_i = npy.linspace(0,self.maxclcolors+2,self.maxclcolors+3)
        
        #use black color for clusters > maxclcolors, the other 2 colors are for background      
        colors2=['#0000FF','#FF0000','#DFE32D','#36F200','#B366FF',
                '#FF470A','#33FFFF','#006600','#CCCC99','#993300',
                '#000000','#FFFFFF','#EEEEEE']
        
        self.clusterclrmap2=matplotlib.colors.LinearSegmentedColormap.from_list('clustercm2',colors2)
     
        self.bnorm2 = matplotlib.colors.BoundaryNorm(colors_i, self.clusterclrmap2.N)
        


        
        
#----------------------------------------------------------------------
    def OnCalcClusters(self, event):
       
        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        self.calcclusters = False  
        
        if True:
        #try: 
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
            
#        except:
#            self.com.cluster_calculated = 0
#            QtGui.QApplication.restoreOverrideCursor()     
            
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
                    
                    
                self.selcluster = self.anlz.cluster_indices[self.ix,self.iy] + 1
                
                self.slidershow.setValue(self.selcluster)
                
                self.showClusterSpectrum()
                self.showIndvClusterImage()
            except:
                pass
            
                       
#----------------------------------------------------------------------
    def CalcClusters(self):
        nclusters = self.anlz.calculate_clusters(self.init_nclusters, self.wo_1st_pca, self.sigma_split)
        self.numclusters = nclusters
        self.ntc_clusters_found.setText(str(self.numclusters))
    
        
        
        
#----------------------------------------------------------------------     
#Show composite cluster image  
    def showClusterImage(self):       
        
        self.clusterimage = self.anlz.cluster_indices  
        
        #print self.selpca

        fig = self.clusterimgfig
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        axes = fig.gca()
        
        
        
        im = axes.imshow(self.clusterimage, cmap=self.clusterclrmap1, norm=self.bnorm1)
        axes.axis("off")
        #cbar = axes.figure.colorbar(im)         
        self.ClusterImagePan.draw()


      
#----------------------------------------------------------------------     
#Show composite cluster image  
    def showIndvClusterImage(self):
        
        indvclusterimage = npy.zeros((self.anlz.stack.n_cols, self.anlz.stack.n_rows))+20.      
        ind = npy.where(self.anlz.cluster_indices == self.selcluster-1)    
        colorcl = min(self.selcluster-1,self.maxclcolors-1)
        indvclusterimage[ind] = colorcl

        fig = self.clusterindvimgfig
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        axes = fig.gca()
            
        
        im = axes.imshow(indvclusterimage, cmap=self.clusterclrmap2, norm=self.bnorm2)
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
                clusterspectrum = self.anlz.clusterspectra[i-1, ]/npy.amax(self.anlz.clusterspectra[i-1, ])           
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
               
        
        
        im = axes.imshow(mapimage, cmap=matplotlib.cm.get_cmap('gray'))
        
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
               
        pass
        #SaveWinP3().Show()
        
        
#----------------------------------------------------------------------    
    def Save(self, filename, path, spec_png = True, spec_pdf = False, spec_csv = False,
             img_png = True, img_pdf = False, 
             indimgs_png = True, indimgs_pdf = False,
             scatt_png = True, scatt_pdf = False): 
        
        self.SaveFileName = os.path.join(path,filename)
        
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas    
        matplotlib.rcParams['pdf.fonttype'] = 42
   
        try: 
            
            if img_png:
                ext = 'png'
                suffix = "." + ext
            
                fig = matplotlib.figure.Figure(figsize = (float(self.stk.n_rows)/10, float(self.stk.n_cols)/10))
                canvas = FigureCanvas(fig)
                fig.clf()
                fig.add_axes((0.0,0.0,1.0,1.0))
                axes = fig.gca()         
        
                im = axes.imshow(self.clusterimage, cmap=self.clusterclrmap1, norm=self.bnorm1)
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

        
                im = axes.imshow(self.clusterimage, cmap=self.clusterclrmap1, norm=self.bnorm1)
                axes.axis("off")
                
                fileName_caimg = self.SaveFileName+"_CAcimg."+ext       
                fig.savefig(fileName_caimg, dpi=300, pad_inches = 0.0)
                            
            

                  
            ext = 'png'
            suffix = "." + ext
                
            if indimgs_png:
                for i in range (self.numclusters):
              
                    indvclusterimage = npy.zeros((self.anlz.stack.n_cols, self.anlz.stack.n_rows))+20.      
                    ind = npy.where(self.anlz.cluster_indices == i)    
                    colorcl = min(i,9)
                    indvclusterimage[ind] = colorcl

                    fig = matplotlib.figure.Figure(figsize =(float(self.stk.n_rows)/10, float(self.stk.n_cols)/10))
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.0,0.0,1.0,1.0))
                    axes = fig.gca()            
                    im = axes.imshow(indvclusterimage, cmap=self.clusterclrmap2, norm=self.bnorm2)
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
                    
            if spec_csv:
                for i in range (self.numclusters):
                    clusterspectrum = self.anlz.clusterspectra[i, ]
                    fileName_spec = self.SaveFileName+"_CAspectrum_" +str(i+1)+".csv"
                    self.stk.write_csv(fileName_spec, self.anlz.stack.ev, clusterspectrum)
                                                     
                
            ext = 'pdf'
            suffix = "." + ext
                
            if indimgs_pdf:
                for i in range (self.numclusters):
              
                    indvclusterimage = npy.zeros((self.anlz.stack.n_cols, self.anlz.stack.n_rows))+20.      
                    ind = npy.where(self.anlz.cluster_indices == i)    
                    colorcl = min(i,9)
                    indvclusterimage[ind] = colorcl

                    fig = matplotlib.figure.Figure(figsize =(float(self.stk.n_rows)/30, float(self.stk.n_cols)/30))
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.0,0.0,1.0,1.0))
                    axes = fig.gca()       
                    im = axes.imshow(indvclusterimage, cmap=self.clusterclrmap2, norm=self.bnorm2)
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
                    
            if scatt_png:
                self.SaveScatt(png_pdf = 1)
            if scatt_pdf:
                self.SaveScatt(png_pdf = 2)                   
            
        except IOError, e:
            if e.strerror:
                err = e.strerror 
            else: 
                err = e 
   
            QtGui.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err)

  
  
#----------------------------------------------------------------------    
#If png_pdg = 1 save png, if =2 save pdf
    def SaveScatt(self, png_pdf = 1): 
          
        od_reduced = self.anlz.pcaimages[:,:,0:self.anlz.numsigpca]        
        od_reduced = npy.reshape(od_reduced, (self.stk.n_cols*self.stk.n_rows,self.anlz.numsigpca), order='F')

        clindices = self.anlz.cluster_indices
        clindices = npy.reshape(clindices, (self.stk.n_cols*self.stk.n_rows), order='F')
        
        path, ext = os.path.splitext(self.SaveFileName) 
        ext = ext[1:].lower() 
        
        if png_pdf == 1:
            ext = 'png'
        else:
            ext = 'pdf'
        
   
        try: 
            QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
            
            suffix = "." + ext
            
            nplots = 0
            for ip in range(self.anlz.numsigpca):
                for jp in range(self.anlz.numsigpca):
                    if jp >= (ip+1):
                        nplots = nplots+1    
            nplotsrows = npy.ceil(nplots/2)    
            
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
                            thiscluster = npy.where(clindices == i)
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
        self.od_reduced = npy.reshape(self.od_reduced, (self.ncols*self.nrows,self.numsigpca), order='F')
         
        self.clindices = self.anlz.cluster_indices
        self.clindices = npy.reshape(self.clindices, (self.ncols*self.nrows), order='F')
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
            thiscluster = npy.where(self.clindices == i)
            axes.plot(x_comp[thiscluster], y_comp[thiscluster],'.',color=self.colors[i],alpha=0.5)
            
        axes.set_xlabel('Component '+str(self.pca_x))
        axes.set_ylabel('Component '+str(self.pca_y))
      
        self.ScatterPPanel.draw()
       
  
        
        
            
    
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
        
       
                
        #panel 2
        vbox2 = QtGui.QVBoxLayout()
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
        
        scrollmax = npy.min([self.stk.n_ev, 20])
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
    def OnNPCAspin(self, value):
        num = value
        self.numsigpca = num
        
        if self.com.pca_calculated == 1:      
            self.anlz.numsigpca = self.numsigpca 
        
            # cumulative variance
            var = self.anlz.variance[:self.numsigpca].sum()
            self.vartc.setText(str(var.round(decimals=2)*100)+'%')
                 
       
        
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
    def OnPCAScroll(self, value):
        self.sel = value
        self.selpca = self.sel
        if self.calcpca == True:
            self.loadPCAImage()
            self.loadPCASpectrum()


 
            
#----------------------------------------------------------------------  
    def OnPointEvalsImage(self, evt):
        x = evt.xdata
        y = evt.ydata
                
        if self.com.pca_calculated == 1:     
            #Find the closest point to the point clicked on the plot
            self.selpca = int(npy.round(x))
                       
            self.loadPCAImage()
            self.loadPCASpectrum()
            


#----------------------------------------------------------------------    
    def OnSave(self, event):     

        pass
        #SaveWinP2().Show()


            
#----------------------------------------------------------------------    
    def Save(self, filename, path, spec_png = True, spec_pdf = False, spec_csv = False,
             img_png = True, img_pdf = False, evals_png = True, evals_pdf = False): 
        
        
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
        
                    im = axes.imshow(self.pcaimage, cmap=matplotlib.cm.get_cmap("seismic_r"), vmin = -bound, vmax = bound)
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
                    self.stk.write_csv(fileName_spec, self.stk.ev, pcaspectrum)
                    
                
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
        
                    im = axes.imshow(self.pcaimage, cmap=matplotlib.cm.get_cmap("seismic_r"), vmin = -bound, vmax = bound)
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
            
        except IOError, e:
            if e.strerror:
                err = e.strerror 
            else: 
                err = e 
   
            QtGui.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err)

            
        
     
#----------------------------------------------------------------------      
    def showEvals(self):
        
        evalmax = npy.min([self.stk.n_ev, 40])
        self.pcaevals = self.anlz.eigenvals[0:evalmax]
        

        fig = self.pcaevalsfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        
       
        evalsplot = axes.semilogy(npy.arange(1,evalmax+1), self.pcaevals,'b.')    
        
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
        
     
        im = axes.imshow(self.pcaimage, cmap=matplotlib.cm.get_cmap("seismic_r"), vmin = -bound, vmax = bound)
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
        self.button_align.clicked.connect(self.OnAlignImgs)
        self.button_align.setEnabled(False)
        vbox1.addWidget(self.button_align) 
        
        self.button_i0ffile = QtGui.QPushButton('I0 from file...')
        self.button_i0ffile.clicked.connect(self.OnI0FFile)
        self.button_i0ffile.setEnabled(False)
        vbox1.addWidget(self.button_i0ffile)
        self.button_i0histogram = QtGui.QPushButton('I0 from histogram...')
        self.button_i0histogram.clicked.connect( self.OnI0histogram)   
        self.button_i0histogram.setEnabled(False)     
        vbox1.addWidget(self.button_i0histogram)
        self.button_showi0 = QtGui.QPushButton('Show I0...')
        self.button_showi0.clicked.connect( self.OnShowI0)   
        self.button_showi0.setEnabled(False)
        vbox1.addWidget(self.button_showi0)
        
        
        self.button_limitev = QtGui.QPushButton('Limit energy range...')
        self.button_limitev.clicked.connect( self.OnLimitEv)
        self.button_limitev.setEnabled(False)
        vbox1.addWidget(self.button_limitev)
        
        self.button_subregion = QtGui.QPushButton('Clip to subregion...')
        self.button_subregion.clicked.connect(self.OnCliptoSubregion)
        self.button_subregion.setEnabled(False)
        vbox1.addWidget(self.button_subregion)       
            
        self.button_savestack = QtGui.QPushButton('Save preprocessed stack')
        self.button_savestack.clicked.connect(self.OnSaveStack)
        self.button_savestack.setEnabled(False)          
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
        vbox4 = QtGui.QVBoxLayout()
        
        self.tc_imageeng = QtGui.QLabel(self)
        self.tc_imageeng.setText("Image at energy: ")
        vbox4.addWidget(self.tc_imageeng)
        
        
        hbox41 = QtGui.QHBoxLayout()
        
        
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
   
        self.absimgfig = Figure((PlotH, PlotH))
        self.AbsImagePanel = FigureCanvas(self.absimgfig)
        self.AbsImagePanel.setParent(self)
        self.AbsImagePanel.mpl_connect('button_press_event', self.OnPointAbsimage)
        fbox.addWidget(self.AbsImagePanel)
        frame.setLayout(fbox)
        hbox41.addWidget(frame)        
       

        self.slider_eng = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slider_eng.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_eng.valueChanged[int].connect(self.OnScrollEng)
        self.slider_eng.setRange(0, 100)
        hbox41.addWidget(self.slider_eng)

        
        vbox4.addLayout(hbox41)
        

        #panel 5     
        vbox5 = QtGui.QVBoxLayout()
        
        self.tc_spec = QtGui.QLabel(self)
        self.tc_spec.setText("Spectrum ")
        vbox5.addWidget(self.tc_spec)
        
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
      
        self.specfig = Figure((PlotW, PlotH))
        self.SpectrumPanel = FigureCanvas(self.specfig)
        self.SpectrumPanel.setParent(self)
        self.SpectrumPanel.mpl_connect('button_press_event', self.OnPointSpectrum)

        fbox.addWidget(self.SpectrumPanel)
        frame.setLayout(fbox)
        vbox5.addWidget(frame)                 
    
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
        
    def OnI0FFile(self, event):


        if True:
            
        #try:
                       
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
                self.loadSpectrum(self.ix, self.iy)
                self.loadImage()
                self.com.i0_loaded = 1
                QtGui.QApplication.restoreOverrideCursor()
                
                
            elif extension == '.xas':
                QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))                          

                x=self.stk.n_cols
                y=self.stk.n_rows
                z=self.iev               

                self.ix = x/2
                self.iy = y/2

                self.stk.read_stk_i0(filepath, extension)

                self.loadSpectrum(self.ix, self.iy)
                self.loadImage()
                self.com.i0_loaded = 1
                QtGui.QApplication.restoreOverrideCursor()

            elif extension == '.csv':
                QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))                          

                x=self.stk.n_cols
                y=self.stk.n_rows
                z=self.iev               

                self.ix = x/2
                self.iy = y/2

                self.stk.read_stk_i0(filepath, extension)

                self.loadSpectrum(self.ix, self.iy)
                self.loadImage()
                self.com.i0_loaded = 1
                QtGui.QApplication.restoreOverrideCursor()
            
#         except:
#             QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))  
#             self.com.i0_loaded = 0        
#             QtGui.QApplication.restoreOverrideCursor()
#             QtGui.QMessageBox.warning(self,'Error',"I0 file not loaded.")
#             import sys; print sys.exc_info()
                       
                          
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
        
         
        self.loadSpectrum(self.ix, self.iy)
        self.loadImage()
        
        self.com.i0_loaded = 1
        self.window().refresh_widgets()
        
#----------------------------------------------------------------------    
    def OnShowI0(self, event):     

        plot = PlotFrame(self, self.stk.evi0,self.stk.i0data)
        plot.show()
 
#----------------------------------------------------------------------    
    def OnSaveStack(self, event): 
        
        self.window().SaveStackH5()
               
#----------------------------------------------------------------------    
    def OnSave(self, event):     
        
        savewin = SaveWinP1(self)
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
    def Save(self, filename, path, spec_png = True, spec_pdf = False, sp_csv = False, img_png = True, img_pdf = False, img_all = False): 

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
                    im = axes.imshow(image, cmap=matplotlib.cm.get_cmap(self.colortable)) 
                    axes.axis("off") 
                                
                    fileName_img = self.SaveFileName+"_imnum_" +str(i+1)+"."+ext
                    fig.savefig(fileName_img,  dpi=ImgDpi, pad_inches = 0.0)
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
            
        except IOError, e:
            if e.strerror:
                err = e.strerror 
            else: 
                err = e 
   
            QtGui.QMessageBox.warning(self,'Error','Could not save file: %s' % err)


        
#----------------------------------------------------------------------    
    def OnAlignImgs(self, event):  
        
        #self.window().Hide()
        imgregwin = ImageRegistration(self.window(), self.com, self.stk)
        imgregwin.show()

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
                
                time.sleep(0.01)
                
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
        
        colorwin = ColorTableFrame(self)
        colorwin.show()
                                    
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
            vmin=(imgmin+imgmax*self.brightness_min)
            vmax=imgmax*self.brightness_max
            if vmin > vmax : vmax = vmin + 0.1
            im = axes.imshow(image, cmap=matplotlib.cm.get_cmap(self.colortable), 
                             vmin=vmin,vmax=vmax)
             
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
 
    
 
        specplot = axes.plot(self.stk.ev,self.spectrum)
         
 
         
        axes.axvline(x=self.stk.ev[self.iev], color = 'g', alpha=0.5)
 
         
        self.SpectrumPanel.draw()
         
        self.tc_spec.setText("Spectrum at pixel [" +str(ypos)+", " + str(xpos)+"] or position ["+
                              str(self.stk.x_dist[xpos])+", "+ str(self.stk.y_dist[ypos])+ "]")


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
              
        self.ROIspectrum = npy.zeros((self.stk.n_ev))
               
        indices = npy.where(self.ROIpix == 255)
        numroipix = self.ROIpix[indices].shape[0]
            
        for ie in range(self.stk.n_ev):
            thiseng_od = self.stk.od3d[:,:,ie]
            self.ROIspectrum[ie] = npy.sum(thiseng_od[indices])/numroipix
                 
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
        self.loadImage()
        
        #find pixels inside the polygon 
        if self.ROIpix == None:
            self.ROIpix = npy.zeros((self.stk.n_cols,self.stk.n_rows))    
        
        for i in range(self.stk.n_cols):
            for j in range(self.stk.n_rows):
                Pinside = self.point_in_poly(i, j, self.roixdata, self.roiydata)
                if Pinside == True:
                    self.ROIpix[j,i] = 255
              
        self.ROIpix = npy.ma.array(self.ROIpix)
        
        self.ROIpix_masked =  npy.ma.masked_values(self.ROIpix, 0)
        
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
   
        self.ROIspectrum = npy.zeros((self.stk.n_ev))
               
        indices = npy.where(self.ROIpix == 255)
        numroipix = self.ROIpix[indices].shape[0]
            
        for ie in range(self.stk.n_ev):
            thiseng_abs = self.stk.absdata[:,:,ie]
            self.ROIspectrum[ie] = npy.sum(thiseng_abs[indices])/numroipix
                
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
        self.window().refresh_widgets()
        
#----------------------------------------------------------------------    
    def OnROI_DoseCalc(self, event):  
        
        self.CalcROISpectrum()
        
        DoseCalculation(self.com, self.stk, self.ROIspectrum).Show()
        
        
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
                self.stk.write_csv(fileName, self.stk.ev, self.ROIspectrum)
            else:
                self.CalcROI_I0Spectrum()
                self.stk.write_csv(fileName, self.stk.ev, self.ROIspectrum)
                     
                
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

        
        self.resize(400, 500)
        self.setWindowTitle('Stack File List')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)


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
        averagefluxmax = npy.max(self.stack.histogram)
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
        
        
        histdata =  npy.reshape(self.stack.histogram, (self.stack.n_cols*self.stack.n_rows), order='F')
        
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
                
        hist_indices = npy.where((fluxmin<self.stack.averageflux)&(self.stack.averageflux<fluxmax))

        redpix = npy.zeros((self.stack.n_cols,self.stack.n_rows))    
        redpix[hist_indices] = 255
              
        redpix = npy.ma.array(redpix)
        
        redpix_masked =  npy.ma.masked_values(redpix, 0)
        
               
        fig = self.absimgfig
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        
        
        axes = fig.gca()
        fig.patch.set_alpha(1.0) 
      
        im = axes.imshow(image, cmap=matplotlib.cm.get_cmap("gray")) 

        im_red = axes.imshow(redpix_masked,cmap=matplotlib.cm.get_cmap("autumn"))  
         
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
        self.close()
        self.parent.I0histogramCalculated()




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
        
        odtotal = self.stack.od3d.sum(axis=0)   
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
                
        
        self.limitevmin = npy.abs(self.stack.ev-x1).argmin()

             
        
#----------------------------------------------------------------------        
    def OnSelection2(self, evt):

        x2 = evt.xdata
            

        self.button_pressed = False
        self.SpectrumPanel.mpl_disconnect(self.conn)    
        
        if x2 == None:
            return    
        
                
        #self.limitevmin = npy.abs(self.stack.ev-x1).argmin()
        self.limitevmax = npy.abs(self.stack.ev-x2).argmin()
        
        self.evlimited = 1
             
        self.draw_limitev_plot()
        
        
#----------------------------------------------------------------------        
    def OnSelectionMotion(self, event):        

        x2 = event.xdata
        
        if x2 == None:
            return  
        
        self.limitevmax = npy.abs(self.stack.ev-x2).argmin()
        
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
        
        #print self.stack.n_ev, self.stack.ev.shape
        
        
        self.stack.absdata = self.stack.absdata[:,:,self.limitevmin:self.limitevmax+1]
        
        if self.com.i0_loaded == 1:   
            self.stack.od3d = self.stack.od3d[:,:,self.limitevmin:self.limitevmax+1]
        
            self.stack.od = self.stack.od3d.copy()
        
            self.stack.od = npy.reshape(self.stack.od, (self.stack.n_rows*self.stack.n_cols, self.stack.n_ev), order='F')
        
        self.stack.fill_h5_struct_from_stk()
        
        #Fix the slider on Page 1! 
        self.parent.page1.slider_eng.setRange(0,self.stack.n_ev-1)
        self.parent.page1.iev = self.stack.n_ev/2
        self.parent.page1.slider_eng.setValue(self.parent.page1.iev)
        
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
 
        
        self.new_y1 = int(self.stack.n_cols*0.10)
        self.new_y2 = int(self.stack.n_cols*0.90)
        self.new_x1 = int(self.stack.n_rows*0.10)
        self.new_x2 = int(self.stack.n_rows*0.90)
        
        self.new_ncols = self.new_y2 - self.new_y1
        self.new_nrows = self.new_x2 - self.new_x1
    
    
        
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
      
        im = axes.imshow(image, cmap=matplotlib.cm.get_cmap("gray")) 
        
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

        
        self.new_ncols = self.new_y1 - self.new_y2 + 1
        self.new_nrows = self.new_x2 - self.new_x1 + 1

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
        
        self.new_ncols = self.new_y2 - self.new_y1 + 1
        self.new_nrows = self.new_x2 - self.new_x1 + 1

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
        #Matlab axis are inverted
        self.stack.absdata = self.stack.absdata[self.new_y1:self.new_y2+1, self.new_x1:self.new_x2+1, :]
              
        self.stack.n_cols = self.stack.absdata.shape[0]
        self.stack.n_rows = self.stack.absdata.shape[1]

        
        if self.com.i0_loaded == 1:        
            self.stack.od3d = self.stack.od3d[self.new_y1:self.new_y2+1, self.new_x1:self.new_x2+1, :]
        
            self.stack.od = self.stack.od3d.copy()
        
            self.stack.od = npy.reshape(self.stack.od, (self.stack.n_rows*self.stack.n_cols, self.stack.n_ev), order='F')
        
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
        
        self.resize(850, 700)
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
        self.ref_image_index = 0
        self.ref_image = 0
        
        self.man_align = 0
        self.man_xref = 0
        self.man_yref = 0
        
        self.maxshift = 0
        self.auto = True
        
        self.man_xs = npy.zeros((self.stack.n_ev))
        self.man_ys = npy.zeros((self.stack.n_ev))
        
        self.aligned_stack = self.stack.absdata.copy()
        
        
        self.xshifts = npy.zeros((self.stack.n_ev))
        self.yshifts = npy.zeros((self.stack.n_ev))
        
        self.showccorr = 0
        
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
       
        hbox11 = QtGui.QHBoxLayout()
        
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
   
        self.absimgfig = Figure((3.0, 3.0))
        self.AbsImagePanel = FigureCanvas(self.absimgfig)
        self.AbsImagePanel.setParent(self)
        self.AbsImagePanel.mpl_connect('button_press_event', self.OnPointCorrimage)
        
        fbox.addWidget(self.AbsImagePanel)
        frame.setLayout(fbox)
        hbox11.addWidget(frame)
        

        self.slider_eng = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slider_eng.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_eng.valueChanged[int].connect(self.OnScrollEng)
        self.slider_eng.setRange(0, self.stack.n_ev-1)
        self.slider_eng.setValue(self.iev)
        
        hbox11.addWidget(self.slider_eng)        
      
        vbox1.addWidget(self.tc_imageeng)        
        vbox1.addLayout(hbox11)

             
        
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
         
        self.rb_auto = QtGui.QRadioButton( 'Automatic Alignment', self)
        self.rb_man = QtGui.QRadioButton('Manual Alignment',self)
        self.rb_auto.setChecked(True)
        self.rb_auto.toggled.connect(self.Onrb_automanual)
         
        vbox9.addWidget(self.rb_auto)
        vbox9.addWidget(self.rb_man)
 
                
         
        #panel 8
        sizer8 = QtGui.QGroupBox('This Image')
        vbox8 = QtGui.QVBoxLayout()
        vbox8.setSpacing(0)
 
         
        self.button_refimg = QtGui.QPushButton('Set as Reference Image')
        self.button_refimg.clicked.connect(self.SetRefImage)
        vbox8.addWidget(self.button_refimg)
         
        self.button_remove = QtGui.QPushButton('Remove image')
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
        vbox4.addStretch(1)
         
        self.button_subregion = QtGui.QPushButton('Select subregion on reference')
        self.button_subregion.clicked.connect(self.OnSelectSubregion)   
        self.button_subregion.setEnabled(False)     
        vbox4.addWidget(self.button_subregion)
          
        self.button_delsubregion = QtGui.QPushButton('Remove subregion selection')
        self.button_delsubregion.clicked.connect(self.OnDeleteSubregion)   
        self.button_delsubregion.setEnabled(False)     
        vbox4.addWidget(self.button_delsubregion)
        vbox4.addStretch(1)
          
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
          
        self.tc_shift1 = QtGui.QLabel(self)
        self.tc_shift2 = QtGui.QLabel(self)
        vbox4.addWidget(self.tc_shift1)
        vbox4.addWidget(self.tc_shift2)
        vbox4.addStretch(1)
          
        self.tc_shift1.setText('X shift: {0:5.2f} pixels'.format(self.yshifts[self.iev]))
        self.tc_shift2.setText('Y shift: {0:5.2f} pixels'.format(self.xshifts[self.iev]))
          
  
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
    
        self.refimgfig = Figure((3.0, 3.0))
        self.RefImagePanel = FigureCanvas(self.refimgfig)
        self.RefImagePanel.setParent(self)
         
        fbox.addWidget(self.RefImagePanel)
        frame.setLayout(fbox)
                 
        self.RefImagePanel.mpl_connect('button_press_event', self.OnPointRefimage)
        self.RefImagePanel.mpl_connect('button_release_event', self.OnSelection)
        
        
         
        vbox5.addWidget(self.tc_refimg)        
        vbox5.addWidget(frame)
 
         
         
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
         
        self.button_crop = QtGui.QPushButton('Crop aligned images')
        self.button_crop.clicked.connect(self.OnCropShifts)   
        self.button_crop.setEnabled(False)
        vbox7.addWidget(self.button_crop)
 
        self.button_saveshifts = QtGui.QPushButton('Save image shifts')
        self.button_saveshifts.clicked.connect(self.OnSaveShifts)
        self.button_saveshifts.setEnabled(False)
        vbox7.addWidget(self.button_saveshifts)
         
        self.button_loadshifts = QtGui.QPushButton('Load image shifts')
        self.button_loadshifts.clicked.connect(self.OnLoadShifts)
        vbox7.addWidget(self.button_loadshifts)
                 
        self.button_accept = QtGui.QPushButton('Accept changes')
        self.button_accept.clicked.connect(self.OnAccept)
        self.button_accept.setEnabled(False)
        vbox7.addWidget(self.button_accept)
         
        self.button_close = QtGui.QPushButton('Dismiss changes')
        self.button_close.clicked.connect(self.close)
        vbox7.addWidget(self.button_close)
         
         
         
               
        
        hboxtop = QtGui.QHBoxLayout()
        
        vboxL = QtGui.QVBoxLayout()
        vboxR = QtGui.QVBoxLayout()
        

        vboxL.addWidget(sizer8)
        vboxL.addStretch(1)
        vboxL.addLayout(vbox9)
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
               
        image = self.aligned_stack[:,:,self.iev]
            
            
        fig = self.absimgfig
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        axes = fig.gca()  
             

        im = axes.imshow(image, cmap=matplotlib.cm.get_cmap('gray')) 
       
        if self.man_align == 2:           
            lx=matplotlib.lines.Line2D([self.man_yref-self.man_ys[self.iev],
                                    self.man_yref-self.man_ys[self.iev]], 
                                    [0,self.stack.n_cols],color='red')
            ly=matplotlib.lines.Line2D([0,self.stack.n_rows], 
                                   [self.man_xref-self.man_xs[self.iev],
                                    self.man_xref-self.man_xs[self.iev]] ,color='red')
            axes.add_line(lx)
            axes.add_line(ly)
            
        axes.axis("off") 
        self.AbsImagePanel.draw()
        
        
        self.tc_imageeng.setText('Image at energy: {0:5.2f} eV'.format(float(self.stack.ev[self.iev])))
        
        self.tc_shift1.setText('X shift: {0:5.2f} pixels\n'.format(self.yshifts[self.iev]))
        self.tc_shift2.setText('Y shift: {0:5.2f} pixels'.format(self.xshifts[self.iev]))

        
        if (self.man_align == 2):                  
            self.textctrl_ms1.setText('X manual shift:  {0:5.2f}  pixels\n'.format(self.man_ys[self.iev]))
            self.textctrl_ms2.setText('Y manual shift:  {0:5.2f}  pixels'.format(self.man_xs[self.iev]))
            

        
        
#----------------------------------------------------------------------            
    def OnScrollEng(self, value):
        self.iev = value
        
        self.ShowImage()
            

#----------------------------------------------------------------------        
    def SetRefImage(self):
        
        self.ref_image_index = self.iev
               
        self.ref_image = self.aligned_stack[:,:,self.iev].copy()
        
        self.ShowRefImage()
        self.have_ref_image = 1
        
        self.UpdateWidgets()
        

#----------------------------------------------------------------------        
    def ShowRefImage(self):

        fig = self.refimgfig
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        axes = fig.gca()
        
    
        im = axes.imshow(self.ref_image, cmap=matplotlib.cm.get_cmap('gray')) 
        
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
            lx=matplotlib.lines.Line2D([self.man_yref,self.man_yref], [0,self.stack.n_cols],color='green')
            ly=matplotlib.lines.Line2D([0,self.stack.n_rows], [self.man_xref,self.man_xref] ,color='green')
            axes.add_line(lx)
            axes.add_line(ly)
    
        axes.axis("off") 

        self.RefImagePanel.draw()
        
        self.tc_refimg.setText('Reference image at energy: {0:5.2f} eV'.format(float(self.stack.ev[self.ref_image_index])))

#----------------------------------------------------------------------          
    def Onrb_automanual(self, enabled):
        
        state = enabled      
        
      
        if state:
            self.auto = True
            self.man_align = 0
        else:        
            self.auto = False
            
        self.man_xs = npy.zeros((self.stack.n_ev))
        self.man_ys = npy.zeros((self.stack.n_ev))
        
        self.aligned_stack = self.stack.absdata.copy()      
        self.xshifts = npy.zeros((self.stack.n_ev))
        self.yshifts = npy.zeros((self.stack.n_ev))
        
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
        fig.add_axes((0.02,0.02,0.96,0.96))
        axes = fig.gca()
    
        im = axes.imshow(ccorr, cmap=matplotlib.cm.get_cmap('gray')) 
        
        nx = ccorr.shape[0]
        ny = ccorr.shape[1]
        xcenter = xshift + npy.float(nx)/2.0
        ycenter = yshift + npy.float(ny)/2.0
        
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
                   
        lx=matplotlib.lines.Line2D([ycenter,ycenter], [xl,xr],color='green')
        ly=matplotlib.lines.Line2D([yl,yr], [xcenter,xcenter] ,color='green')
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
    def OnSetMaxShift(self, value):

        self.maxshift = value      

#----------------------------------------------------------------------        
    def OnRemoveImage(self, event):
        
        self.stack.absdata = npy.delete(self.stack.absdata, self.iev, axis=2)  
        
        self.aligned_stack = npy.delete(self.aligned_stack, self.iev, axis=2)
        
        self.stack.n_ev = self.stack.n_ev - 1
        self.stack.ev = npy.delete(self.stack.ev, self.iev) 
        
        
        self.xshifts = npy.delete(self.xshifts, self.iev) 
        self.yshifts = npy.delete(self.yshifts, self.iev) 
        
        self.stack.data_struct.exchange.data = self.stack.absdata
        self.stack.data_struct.exchange.energy = self.stack.ev
        
        if self.com.i0_loaded == 1:
            self.stack.calculate_optical_density()
            
        self.iev = self.iev-1
        if self.iev < 0:
            self.iev = 0

        self.ShowImage() 
        
        
        self.parent.page1.slider_eng.setRange(0,self.stack.n_ev-1)
        self.parent.page1.iev = self.stack.n_ev/2
        self.parent.page1.slider_eng.setValue(self.parent.page1.iev)
        
        self.parent.page1.loadSpectrum(self.parent.page1.ix, self.parent.page1.iy)
        self.parent.page1.loadImage() 

#----------------------------------------------------------------------            
    def OnCalcRegistration(self, event):
        
        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        
        #Subregion selection on a reference image
        if self.subregion == 0:
            self.sr_x1 = 0
            self.sr_x2 = 0
            self.sr_y1 = 0
            self.sr_y2 = 0            
            referenceimage = self.ref_image            
        else:    
            referenceimage = self.ref_image[self.sr_y2:self.sr_y1, self.sr_x1:self.sr_x2]  

        for i in range(self.stack.n_ev):

            if self.subregion == 0:
                img2 = self.aligned_stack[:,:,i]  
            else:
                img2 = self.aligned_stack[self.sr_y2:self.sr_y1, self.sr_x1:self.sr_x2, i]  
               
            if i==0:     
                xshift, yshift, ccorr = self.stack.register_images(referenceimage, img2, 
                                                          have_ref_img_fft = False)   
            elif i==self.ref_image_index:
                xshift = 0
                yshift = 0       
            else:
                xshift, yshift, ccorr = self.stack.register_images(referenceimage, img2, 
                                                          have_ref_img_fft = True)
            
            #Limit the shifts to MAXSHIFT chosen by the user
            if (self.maxshift > 0):
                if (abs(xshift) > self.maxshift):
                        xshift = npy.sign(xshift)*self.maxshift
                if (abs(yshift) > self.maxshift):
                        yshift = npy.sign(yshift)*self.maxshift
            
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
        
            
        QtGui.QApplication.restoreOverrideCursor()
        
#----------------------------------------------------------------------            
    def OnCropShifts(self, event):
        
        QtGui.QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        
        self.aligned_stack = self.stack.crop_registed_images(self.aligned_stack, 
                                                             self.xshifts,
                                                             self.yshifts)
        
        
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
            self.man_xref = int(npy.floor(y))           
            self.man_yref = int(npy.floor(x))  
                    
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
            xcorr = int(npy.floor(y))           
            ycorr = int(npy.floor(x))  
                    
            if xcorr<0 :
                xcorr=0
            if xcorr>self.stack.n_cols :
                xcorr=self.stack.n_cols
            if ycorr<0 :
                ycorr=0
            if ycorr>self.stack.n_rows :
                ycorr=self.stack.n_rows 
                
        
            self.man_xs[self.iev] = self.man_xref - xcorr
            self.man_ys[self.iev] = self.man_yref - ycorr
                
            
            self.textctrl_ms1.setText('X manual shift:  {0:5.2f}  pixels\n'.format(self.man_ys[self.iev]))
            self.textctrl_ms2.setText('Y manual shift:  {0:5.2f}  pixels'.format(self.man_xs[self.iev]))
            
            self.iev = self.iev + 1
            if self.iev > (self.stack.n_ev-1):
                self.iev = 0
            
            self.slider_eng.setValue(self.iev)
                    
           
            self.ShowImage()
            self.UpdateWidgets()
            
#----------------------------------------------------------------------            
    def OnApplyManShifts(self):
        
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
                
        self.regist_calculated = 1
        self.man_align = 0
        
        self.ShowRefImage()
        self.PlotShifts()
        
        self.textctrl_ms1.setText('X manual shift: \n')
        self.textctrl_ms2.setText('Y manual shift: ')
        
        self.UpdateWidgets()
        

        self.ShowImage()
        
        
      
                    
#----------------------------------------------------------------------            
    def OnAccept(self):
        
        self.stack.absdata = self.aligned_stack  
        
        datadim = npy.int32(self.stack.absdata.shape)
        
        
        self.stack.n_cols = datadim[0].copy()
        self.stack.n_rows =  datadim[1].copy()
        
        self.stack.xshifts = self.xshifts
        self.stack.yshifts = self.yshifts
                      
        if self.com.i0_loaded == 1:
            self.stack.calculate_optical_density()

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
    def PlotShifts(self):
        
        fig = self.shiftsfig
        fig.clf()
        
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        
        
        
        #Matplotlib has inverted axes! 
        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Shifts (x-red, y-green) [pixels]')

        plot = axes.plot(self.stack.ev,self.xshifts, color='green')
        plot = axes.plot(self.stack.ev,self.yshifts, color='red')
        
        
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
                indx = npy.abs(self.stack.ev - x).argmin()
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
        

        self.xshifts = npy.zeros((self.stack.n_ev))
        self.yshifts = npy.zeros((self.stack.n_ev)) 
       
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
            else:
                self.button_register.setEnabled(False)
                self.button_subregion.setEnabled(False)
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
            else:
                self.button_manalign.setEnabled(False)
            
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
        
        self.resize(400, 500)
        self.setWindowTitle('Spectral Regions of Interest')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal) 


#---------------------------------------------------------------------- 
class DoseCalculation(QtGui.QDialog):

    def __init__(self, parent,  common, stack):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        self.stack = stack
        self.com = common  
        
        self.resize(400, 500)
        self.setWindowTitle('Dose Calculation')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal) 
                

#---------------------------------------------------------------------- 
class PlotFrame(QtGui.QDialog):

    def __init__(self, parent, datax, datay, title = "I0 data"):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        self.datax = datax
        self.datay = datay
        
        self.resize(400, 500)
        self.setWindowTitle(title)
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal) 

#---------------------------------------------------------------------- 
class ColorTableFrame(QtGui.QDialog):

    def __init__(self, parent):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        
        self.resize(400, 500)
        self.setWindowTitle('Pick Color Table')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)                                                    

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
        
        self.initMatplotlib()

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
        
        frame = QtGui.QFrame()
        frame.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        fbox = QtGui.QHBoxLayout()
   
        self.absimgfig = Figure((PlotH*.9, PlotH*.9))
        self.AbsImagePanel = FigureCanvas(self.absimgfig)
        self.AbsImagePanel.setParent(self)
        
        fbox.addWidget(self.AbsImagePanel)
        frame.setLayout(fbox)
        hbox51.addWidget(frame)
        

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
    def initMatplotlib(self):  
               
       
        matplotlib.rcParams['figure.facecolor'] = 'white'

        matplotlib.rcParams['font.size'] = 10.0
        
        
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
        
        self.parent = parent

        self.data_struct = data_struct
        self.stk = stack
        self.common = com
        
        self.resize(400, 500)
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
        
        self.filelist.setColumnWidth(0,200)
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
                import sm_netcdf
                self.sm = sm_netcdf.sm(data_struct)
            
            except:
                QtGui.QMessageBox.warning(self, 'Error', "Could not import netCDF4 library.")
                return

            count = 0
        
            for i in range(len(self.sm_files)):
                #print sm_files
                filename = self.sm_files[i]
                file = os.path.join(filepath, filename)
            
                filever, ncols, nrows, iev = self.sm.read_sm_header(file)
            
                
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
                
            return
         
            
        self.xrm_files = [x for x in os.listdir(filepath) if x.endswith('.xrm')] 
        

        if self.xrm_files:        
            
            self.filetype = 'xrm'
            
            import xradia_xrm
            self.xrm = xradia_xrm.xrm()

            count = 0
        
            for i in range(len(self.xrm_files)):

                filename = self.xrm_files[i]
                file = os.path.join(filepath, filename)

                ncols, nrows, iev = self.xrm.read_xrm_fileinfo(file)
  
                if ncols > 0:                       
                    self.filelist.insertRow(count)
                    self.filelist.setRowHeight(count,20)

                    self.filelist.setItem(count, 0, QtGui.QTableWidgetItem(filename))
                    self.filelist.setItem(count, 1, QtGui.QTableWidgetItem(str(ncols)))
                    self.filelist.setItem(count, 2, QtGui.QTableWidgetItem(str(nrows)))
                    self.filelist.setItem(count, 3, QtGui.QTableWidgetItem('{0:5.2f}'.format(iev)))
                            
                    count += 1

            
        self.sm_files = self.xrm_files
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
            self.sm.read_sm_list(filelist, self.filepath, self.data_struct)
        elif self.filetype == 'xrm':
            self.xrm.read_xrm_list(filelist, self.filepath, self.data_struct)
        else:
            print 'Wrong file type'
            return
        
        
        #fill the gui structure data
        self.stk.absdata = self.data_struct.exchange.data
        
        datadim = npy.int32(self.stk.absdata.shape)
        self.stk.n_cols = datadim[0].copy()
        self.stk.n_rows =  datadim[1].copy()
        self.stk.ev = self.data_struct.exchange.energy 
        self.stk.n_ev = npy.int32(self.stk.ev.shape[0]).copy()
        
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
        


        
        # create the page windows as tabs
        self.page0 = PageLoadData(self.common, self.data_struct, self.stk)
        self.page1 = PageStack(self.common, self.data_struct, self.stk)
        self.page2 = PagePCA(self.common, self.data_struct, self.stk, self.anlz)
        self.page3 = PageCluster(self.common, self.data_struct, self.stk, self.anlz)
        self.page4 = PageSpectral(self.common, self.data_struct, self.stk, self.anlz)
        self.page5 = PageKeyEng(self.common, self.data_struct, self.stk, self.anlz)
        
        
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
        self.toolbar.addAction(self.actionOpenSL)
        self.actionOpenSL.triggered.connect(self.BuildStack)
        
        self.actionSave = QtGui.QAction(self)
        self.actionSave.setObjectName('actionSave')
        self.actionSave.setIcon(QtGui.QIcon(os.path.join('images','media-floppy.png')))
        self.toolbar.addAction(self.actionSave)
        self.actionSave.triggered.connect(self.onSaveAsH5)
        self.actionSave.setEnabled(False)

        self.actionInfo = QtGui.QAction(self)
        self.actionInfo.setObjectName('actionInfo')
        self.actionInfo.setIcon(QtGui.QIcon(os.path.join('images','help-browser.png')))
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
                   

        self.refresh_widgets()


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
#             QtGui.QMessageBox.warning(self,'Error',".sm files not loaded.")
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
            self.page1.button_limitev.setEnabled(False)
            self.page1.button_subregion.setEnabled(False)
            self.page1.button_showi0.setEnabled(False) 
            self.page1.rb_flux.setEnabled(False)
            self.page1.rb_od.setEnabled(False)
            self.page2.button_calcpca.setEnabled(False)
            self.page4.button_loadtspec.setEnabled(False)
            self.page4.button_addflat.setEnabled(False)
        else:
            self.page1.button_limitev.setEnabled(True)
            self.page1.button_subregion.setEnabled(True)
            self.page1.button_showi0.setEnabled(True)
            self.page1.rb_flux.setEnabled(True)
            self.page1.rb_od.setEnabled(True)   
            self.page2.button_calcpca.setEnabled(True) 
            self.page4.button_loadtspec.setEnabled(True)
            self.page4.button_addflat.setEnabled(True)   
             
             
             
        if self.common.pca_calculated == 0:      
            self.page2.button_savepca.setEnabled(False)
            self.page2.slidershow.setEnabled(False) 
            self.page3.button_calcca.setEnabled(False)
            self.page4.rb_fit.setEnabled(False)
            self.page5.button_calckeng.setEnabled(False)
        else:
            self.page2.button_savepca.setEnabled(True)
            self.page2.slidershow.setEnabled(True)
            self.page3.button_calcca.setEnabled(True)  
            self.page4.rb_fit.setEnabled(True)  
            self.page5.button_calckeng.setEnabled(True)       
             
        if self.common.cluster_calculated == 0:   
            self.page3.button_scatterplots.setEnabled(False)
            self.page3.button_savecluster.setEnabled(False)
            self.page3.slidershow.setEnabled(False)
            self.page4.button_addclspec.setEnabled(False)
        else:
            self.page3.button_scatterplots.setEnabled(True)
            self.page3.button_savecluster.setEnabled(True)  
            self.page3.slidershow.setEnabled(True)
            self.page4.button_addclspec.setEnabled(True)
             
        if self.common.spec_anl_calculated == 0:
            self.page4.button_removespec.setEnabled(False)
            self.page4.button_movespdown.setEnabled(False)
            self.page4.button_movespup.setEnabled(False)
            self.page4.button_save.setEnabled(False)
            self.page4.button_showrgb.setEnabled(False) 
        else:
            self.page4.button_removespec.setEnabled(True)
            self.page4.button_movespdown.setEnabled(True)
            self.page4.button_movespup.setEnabled(True)
            self.page4.button_save.setEnabled(True) 
            self.page4.button_showrgb.setEnabled(True)          
                   
             
        self.page1.ResetDisplaySettings()
             
             
            
#----------------------------------------------------------------------        
    def new_stack_refresh(self):
        
        
        self.common.i0_loaded = 0
        self.common.pca_calculated = 0
        self.common.cluster_calculated = 0
        
        self.refresh_widgets()
               
        
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
        self.page5.button_calckeng.setEnabled(False)
        fig = self.page5.kespecfig
        fig.clf()
        self.page5.KESpecPan.draw()         
        fig = self.page5.absimgfig
        fig.clf()
        self.page5.AbsImagePanel.draw()   
        self.page5.lc_1.clear()  
        self.page5.keyenergies = []
        self.page5.keyengs_calculated = 0    
                
""" ------------------------------------------------------------------------------------------------"""
                        
def main():
    
    app = QtGui.QApplication(sys.argv)
    frame = MainFrame()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()