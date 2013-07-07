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
from PyQt4.QtCore import Qt


import matplotlib 
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid import make_axes_locatable

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
    def __init__(self, common, data_struct, stack, anlz):
        super(PageKeyEng, self).__init__()

        self.initUI(common, data_struct, stack, anlz)
        
#----------------------------------------------------------------------          
    def initUI(self, common, data_struct, stack, anlz): 
        pass
    

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
        
        self.SetBackgroundColour("White")
           
        self.fontsize = self.com.fontsize        
        
        vbox = QtGui.QVBoxLayout()
        hboxT = QtGui.QHBoxLayout()
        hboxB = QtGui.QHBoxLayout()
    
        #panel 1        
        vbox1 = QtGui.QVBoxLayout()     
  
        self.tc_spmap = wx.TextCtrl(panel1, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_spmap.SetFont(self.com.font)
        self.tc_spmap.SetValue("Spectrum composition map")

        i1panel = wx.Panel(panel1, -1, style = wx.SUNKEN_BORDER)
        self.MapPanel = wxmpl.PlotPanel(i1panel, -1, size =(PlotH, PlotH), cursor=False, crosshairs=False, location=False, zoom=False)                            
  
        vbox1.Add((0,10))
        vbox1.Add(self.tc_spmap,1, wx.LEFT | wx.EXPAND, 20)        
        vbox1.Add(i1panel, 0,  wx.LEFT, 20)


     
     
#         #panel 2
#         panel2 = wx.Panel(self, -1)
#         vbox2 = wx.BoxSizer(wx.VERTICAL)
#         panel2.SetBackgroundColour("White")
#         
#         self.tc_tspec = wx.TextCtrl(panel2, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
#         self.tc_tspec.SetFont(self.com.font)
#         self.tc_tspec.SetValue("Target Spectrum: ")
#         hbox11 = wx.BoxSizer(wx.HORIZONTAL)         
#         
#         i2panel = wx.Panel(panel2, -1, style = wx.SUNKEN_BORDER)
#         self.TSpectrumPanel = wxmpl.PlotPanel(i2panel, -1, size=(PlotW, PlotH), cursor=False, crosshairs=False, location=False, zoom=False)
# 
#         self.slider_tspec = wx.Slider(panel2, -1, 1, 1, 5, style=wx.SL_LEFT|wx.SL_VERTICAL )        
#         self.slider_tspec.SetFocus()
#         self.Bind(wx.EVT_SCROLL, self.OnTSScroll, self.slider_tspec)
#         
#         vbox21 = wx.BoxSizer(wx.VERTICAL)               
#         self.tspecspin = wx.SpinButton(panel2, -1, size = ((8,-1)), style=wx.SP_ARROW_KEYS)
#         self.Bind(wx.EVT_SPIN_UP, self.OnTspecSpinUp, self.tspecspin)
#         self.Bind(wx.EVT_SPIN_DOWN, self.OnTspecSpinDown, self.tspecspin)
#         
#         vbox21.Add((0,3))
#         vbox21.Add(self.slider_tspec, 1,  wx.EXPAND) 
#         vbox21.Add(self.tspecspin, 0,  wx.EXPAND)      
# 
#         hbox11.Add(i2panel, 0)
#         hbox11.Add(vbox21, 0,  wx.EXPAND)
#           
#         vbox2.Add((0,10))
#         vbox2.Add(self.tc_tspec, 1, wx.LEFT | wx.EXPAND, 20)       
#         vbox2.Add(hbox11, 0, wx.LEFT , 20)
#         
#         panel2.SetSizer(vbox2)
#         
#         
#         #panel 3
#         panel3 = wx.Panel(self, -1)
#         sb = wx.StaticBox(panel3, -1, 'Target Spectrum')
#         sb.SetBackgroundColour("white")
#         sizer1 = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)
#         vbox31 = wx.BoxSizer(wx.VERTICAL)
#         vbox31.Add((0,10)) 
#         
#         self.button_loadtspec = wx.Button(panel3, -1, 'Load Spectrum')
#         self.button_loadtspec.SetFont(self.com.font)
#         self.Bind(wx.EVT_BUTTON, self.OnTSpecFromFile, id=self.button_loadtspec.GetId())
#         self.button_loadtspec.Disable()
#         vbox31.Add(self.button_loadtspec, 0, wx.EXPAND)
#         self.button_addflat = wx.Button(panel3, -1, 'Add Flat Spectrum')
#         self.button_addflat.SetFont(self.com.font)
#         self.Bind(wx.EVT_BUTTON, self.OnFlatTSpec, id=self.button_addflat.GetId())
#         self.button_addflat.Disable()
#         vbox31.Add(self.button_addflat, 0, wx.EXPAND)
#         self.button_addclspec = wx.Button(panel3, -1, 'Add Cluster Spectra')
#         self.button_addclspec.SetFont(self.com.font)
#         self.Bind(wx.EVT_BUTTON, self.OnAddClusterSpectra, id=self.button_addclspec.GetId())   
#         self.button_addclspec.Disable()     
#         vbox31.Add(self.button_addclspec, 0, wx.EXPAND)
#         
#         self.button_showrgb = wx.Button(panel3, -1, 'Composite RGB image...')
#         self.button_showrgb.SetFont(self.com.font)
#         self.Bind(wx.EVT_BUTTON, self.OnCompositeRGB, id=self.button_showrgb.GetId())   
#         self.button_showrgb.Disable()     
#         vbox31.Add(self.button_showrgb, 0, wx.EXPAND)        
# 
#         self.button_save = wx.Button(panel3, -1, 'Save Images...', (10,10))
#         self.button_save.SetFont(self.com.font)
#         self.Bind(wx.EVT_BUTTON, self.OnSave, id=self.button_save.GetId())
#         self.button_save.Disable()          
#         vbox31.Add(self.button_save, 0, wx.EXPAND)
#         sizer1.Add(vbox31,1, wx.LEFT|wx.RIGHT|wx.EXPAND,2)
#         panel3.SetSizer(sizer1)
#         
# 
#         
#         #panel 4
#         panel4 = wx.Panel(self, -1)
#         vbox4 = wx.BoxSizer(wx.VERTICAL)
#         sb = wx.StaticBox(panel4, -1, 'Display')
#         sb.SetBackgroundColour("white")
#         sizer4 = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)
# 
#         sb = wx.StaticBox(panel4, -1, 'Spectrum')
#         sb.SetBackgroundColour("white")
#         sizer41 = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)
#         self.textctrl_sp = wx.TextCtrl(panel4, -1, style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
#         self.textctrl_sp.SetFont(self.com.font)
#         sizer41.Add(self.textctrl_sp, 1, wx.EXPAND|wx.TOP|wx.LEFT, 5)
#       
#         hbox40 = wx.BoxSizer(wx.HORIZONTAL)    
#         hbox40.Add(sizer41, 1, wx.EXPAND)
#         vbox4.Add(hbox40, 0, wx.EXPAND)
#         
#         self.textctrl_sp.AppendText('Common Name: \n')
#         self.textctrl_sp.AppendText('RMS Error: ')
#         
#         hbox41 = wx.BoxSizer(wx.HORIZONTAL)
#         
#         sb = wx.StaticBox(panel4, -1, 'Composition Map')
#         sb.SetBackgroundColour("white")
#         sizer42 = wx.StaticBoxSizer(sb,  orient=wx.VERTICAL)
#         self.rb_raw = wx.RadioButton(panel4, -1, 'Raw', style=wx.RB_GROUP)
#         self.rb_fit = wx.RadioButton(panel4, -1, 'Fitted')
#         self.rb_raw.SetFont(self.com.font)
#         self.rb_fit.SetFont(self.com.font)
#         self.Bind(wx.EVT_RADIOBUTTON, self.OnRBRawFit, id=self.rb_raw.GetId())
#         self.Bind(wx.EVT_RADIOBUTTON, self.OnRBRawFit, id=self.rb_fit.GetId())
#         self.rb_raw.SetValue(True)
# 
#         self.add_scale_cb = wx.CheckBox(panel4, -1, '  Scale')
#         self.add_scale_cb.SetFont(self.com.font)
#         self.Bind(wx.EVT_CHECKBOX, self.OnShowScale, self.add_scale_cb)
#         
#         sizer42.Add((0,3))
#         sizer42.Add(self.rb_raw)
#         sizer42.Add((0,5))
#         sizer42.Add(self.rb_fit)
#         sizer42.Add((0,10))
#         sizer42.Add(self.add_scale_cb)
#         
#         hbox41.Add(sizer42, 1, wx.EXPAND)
#                 
#         sb = wx.StaticBox(panel4, -1, 'Fit Weights')
#         sb.SetBackgroundColour("white")
#         sizer43 = wx.StaticBoxSizer(sb,  orient=wx.VERTICAL)
#         hbox42 = wx.BoxSizer(wx.HORIZONTAL)
#         hbox42.Add((3,0))
#         
#         self.tc_spfitlist = wx.TextCtrl(panel4, -1, style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
#         self.tc_spfitlist.SetFont(self.com.font)
#         
#         hbox42.Add(self.tc_spfitlist,1,wx.EXPAND)
#         
#         sizer43.Add(hbox42,1,wx.EXPAND)       
# 
#         
#         vbox43 = wx.BoxSizer(wx.VERTICAL)
#         self.button_removespec = wx.Button(panel4, -1, 'Remove Spectrum')
#         self.button_removespec.SetFont(self.com.font)
#         self.Bind(wx.EVT_BUTTON, self.OnRemoveSpectrum, id=self.button_removespec.GetId())   
#         self.button_removespec.Disable()     
#         vbox43.Add(self.button_removespec, 0, wx.EXPAND)
#         self.button_movespup = wx.Button(panel4, -1, 'Move Spectrum Up')
#         self.button_movespup.SetFont(self.com.font)
#         self.Bind(wx.EVT_BUTTON, self.OnMoveSpectrumUp, id=self.button_movespup.GetId())   
#         self.button_movespup.Disable()     
#         vbox43.Add(self.button_movespup, 0, wx.EXPAND)
#         self.button_movespdown = wx.Button(panel4, -1, 'Move Spectrum Down')
#         self.button_movespdown.SetFont(self.com.font)
#         self.Bind(wx.EVT_BUTTON, self.OnMoveSpectrumDown, id=self.button_movespdown.GetId())   
#         self.button_movespdown.Disable()     
#         vbox43.Add(self.button_movespdown, 0, wx.EXPAND)
#                        
# 
#                
#             
#         hbox41.Add((10,0))
#         hbox41.Add(sizer43, 1, wx.EXPAND)
#         hbox41.Add(vbox43, 1, wx.EXPAND|wx.ALL, 5)
#         vbox4.Add(hbox41, 1, wx.EXPAND)
# 
#         sizer4.Add(vbox4,1, wx.EXPAND)
#         
#         panel4.SetSizer(sizer4)
#         
#         
#         #panel 5
#         panel5 = wx.Panel(self, -1)
#         sb = wx.StaticBox(panel5, -1, 'Target Spectra')
#         sb.SetBackgroundColour("white")
#         sizer5 = wx.StaticBoxSizer(sb, orient= wx.VERTICAL)
#         
#         hbox51 = wx.BoxSizer(wx.HORIZONTAL)
#         hbox51.Add((0,2))
# 
#         self.tc_speclist =  wx.ListCtrl(panel5, -1, 
#                                         style=wx.LC_REPORT|wx.LC_NO_HEADER|wx.NO_BORDER|wx.LC_EDIT_LABELS|wx.LC_SINGLE_SEL)
#         self.tc_speclist.InsertColumn(0, 'Spectra')
#         self.Bind(wx.EVT_LIST_ITEM_FOCUSED , self.OnSpectraListClick, self.tc_speclist)
#         self.Bind(wx.EVT_LIST_END_LABEL_EDIT, self.OnEditSpectraListClick, self.tc_speclist)
#         self.tc_speclist.SetBackgroundColour('white')
#         self.tc_speclist.SetFont(self.com.font)
#         hbox51.Add(self.tc_speclist, 1, wx.EXPAND)
#         sizer5.Add(hbox51,1, wx.ALL|wx.EXPAND,2)        
#         panel5.SetSizer(sizer5)
#         
# #        self.tc_speclist = wx.TextCtrl(panel5, -1, style=wx.TE_MULTILINE|wx.TE_RICH|wx.BORDER_NONE)
# #        self.tc_speclist.SetFont(self.com.font)
# #        hbox51.Add(self.tc_speclist, 1, wx.EXPAND)
# #        sizer5.Add(hbox51,1, wx.ALL|wx.EXPAND,2)        
# #        panel5.SetSizer(sizer5)
# 
#         
#         hboxB.Add(panel2, 0, wx.BOTTOM | wx.TOP, 9)
#         hboxB.Add(panel1, 0, wx.BOTTOM | wx.TOP, 9)
#         hboxT.Add((10,0)) 
#                
#         hboxT.Add(panel3, 0, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 9)
#         hboxT.Add(panel4, 2.5, wx.LEFT | wx.RIGHT |wx.TOP | wx.EXPAND, 9)
#         hboxT.Add(panel5, 1, wx.LEFT | wx.RIGHT |wx.TOP | wx.EXPAND, 9)
# 
#         vbox.Add(hboxT, 0, wx.ALL, 5)
#         
#         vbox.Add((0, 5))
#         
#         vbox.Add(hboxB, 0, wx.LEFT | wx.RIGHT, 5)
#   
#         self.SetSizer(vbox) 
        
        
    
    
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
        
        self.MakeColorTable()             
        

    
        
        #panel 1
        sizer1 = QtGui.QGroupBox('Cluster analysis')
        vbox1 = QtGui.QVBoxLayout()
        #vbox1.setSpacing(0)
        
        self.button_calcca = QtGui.QPushButton('Calculate Clusters')
        #self.button_calcca.clicked.connect( self.OnCalcClusters)   
        self.button_calcca.setEnabled(False)
        vbox1.addWidget(self.button_calcca)
        self.button_scatterplots = QtGui.QPushButton('Show scatter plots...')
        #self.button_scatterplots.clicked.connect( self.OnShowScatterplots)
        self.button_scatterplots.setEnabled(False)
        vbox1.addWidget(self.button_scatterplots)
        self.button_savecluster = QtGui.QPushButton('Save CA Results...')
        #self.button_savecluster.clicked.connect( self.OnSave)
        self.button_savecluster.setEnabled(False)
        vbox1.addWidget(self.button_savecluster)
        

                
        hbox11 = QtGui.QHBoxLayout()
        text1 = QtGui.QLabel(self)
        text1.setText('Number of clusters')
        hbox11.addWidget(text1)
        self.nclusterspin = QtGui.QSpinBox()
        self.nclusterspin.setRange(2,20)
        self.nclusterspin.setValue(self.init_nclusters)
        #self.Bind(wx.EVT_SPINCTRL, self.OnNClusterspin, self.nclusterspin)
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
        #self.remove1stpcacb.stateChanged.connect(self.OnRemove1stpca)
        hbox13.addWidget(self.remove1stpcacb)
        
        vbox1.addLayout(hbox13) 
 
        hbox14 = QtGui.QHBoxLayout()
        self.cb_splitclusters = QtGui.QCheckBox('Divide clusters with large Sigma', self)
        #self.cb_splitclusters.stateChanged.connect(self.OnSplitClusters)
        hbox14.addWidget(self.cb_splitclusters)
        
        vbox1.addLayout(hbox14) 
        
        sizer1.setLayout(vbox1)
        
                
                 
        #panel 2        
        vbox2 = QtGui.QVBoxLayout()
         
        tc_clustercomp = QtGui.QLabel(self)
        tc_clustercomp.setText("Composite cluster image")  
        vbox2.addWidget(tc_clustercomp)      
          
        self.clusterimgfig = Figure((PlotH, PlotH))
        self.ClusterImagePan = FigureCanvas(self.clusterimgfig)
        #wxmpl.EVT_POINT(i2panel, self.ClusterImagePan.GetId(), self.OnPointClusterImage) 
        self.ClusterImagePan.setParent(self)
        vbox2.addWidget(self.ClusterImagePan)
         
         
        #panel 3 
        vbox3 = QtGui.QVBoxLayout()
        fgs = QtGui.QGridLayout()
         
        self.tc_cluster = QtGui.QLabel(self)
        self.tc_cluster.setText("Cluster ")
        fgs.addWidget(self.tc_cluster, 0, 0, QtCore .Qt. AlignLeft)
 
        self.clusterindvimgfig = Figure((PlotH*0.73, PlotH*0.73))
        self.ClusterIndvImagePan = FigureCanvas(self.clusterindvimgfig)
        self.ClusterIndvImagePan.setParent(self)
        fgs.addWidget(self.ClusterIndvImagePan, 1, 0, QtCore .Qt. AlignLeft)
         
        self.slidershow = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slidershow.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slidershow.setEnabled(False)
        #self.slidershow.valueChanged[int].connect(self.OnClusterScroll)
        self.slidershow.setRange(1, 20)
        fgs.addWidget(self.slidershow, 1, 1, QtCore .Qt. AlignLeft)
            
         
        text3 = QtGui.QLabel(self)
        text3.setText('Cluster Distance Map')
        fgs.addWidget(text3, 0, 2, QtCore .Qt. AlignLeft)
        self.clusterdistmapfig = Figure((PlotH*0.73, PlotH*0.73))
        self.ClusterDistMapPan = FigureCanvas(self.clusterdistmapfig)
        self.ClusterDistMapPan.setParent(self)
        fgs.addWidget(self.ClusterDistMapPan, 1, 2, QtCore .Qt. AlignLeft)          
 
        vbox3.addLayout(fgs)
 
         
        

        #panel 4 
        vbox4 = QtGui.QVBoxLayout()
         
        self.tc_clustersp = QtGui.QLabel(self)
        self.tc_clustersp.setText("Cluster spectrum")   
        vbox4.addWidget(self.tc_clustersp)      
 
        self.clusterspecfig = Figure((PlotW, PlotH))
        self.ClusterSpecPan = FigureCanvas(self.clusterspecfig)
        self.ClusterSpecPan.setParent(self)
        vbox4.addWidget(self.ClusterSpecPan)
 
         
 
        
        
        vboxtop = QtGui.QVBoxLayout()
        gridsizertop = QtGui.QGridLayout()
        
        gridsizertop.addWidget(sizer1, 0, 0, QtCore .Qt. AlignLeft)
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
        
    
    
""" ------------------------------------------------------------------------------------------------"""
class PagePCA(QtGui.QWidget):
    def __init__(self, common, data_struct, stack, anlz):
        super(PagePCA, self).__init__()

        self.initUI(common, data_struct, stack, anlz)
        
#----------------------------------------------------------------------          
    def initUI(self, common, data_struct, stack, anlz): 


        
        #panel 1        
        vbox1 = QtGui.QVBoxLayout()
        
        self.tc_PCAcomp = QtGui.QLabel(self)
        self.tc_PCAcomp.setText("PCA component ")
        vbox1.addWidget(self.tc_PCAcomp)
        
        hbox11 = QtGui.QHBoxLayout()
   
        self.pcaimgfig = Figure((PlotH*1.10, PlotH))
        self.PCAImagePan = FigureCanvas(self.pcaimgfig)
        self.PCAImagePan.setParent(self)
        hbox11.addWidget(self.PCAImagePan)
                       
            
        self.slidershow = QtGui.QScrollBar(QtCore.Qt.Vertical)
        self.slidershow.setFocusPolicy(QtCore.Qt.StrongFocus)
        #self.slidershow.setEnabled(False)          
        #self.slidershow.SetFocus()
        #self.Bind(wx.EVT_SCROLL, self.OnPCAScroll, self.slidershow)
        
        hbox11.addWidget(self.slidershow)

        vbox1.addLayout(hbox11)
        
       
                
        #panel 2
        vbox2 = QtGui.QVBoxLayout()
        sizer2 = QtGui.QGroupBox('PCA')
        vbox21 = QtGui.QVBoxLayout()        
        
        self.button_calcpca = QtGui.QPushButton('Calculate PCA')
        #.clicked.connect( self.OnCalcPCA, id=self.button_calcpca.GetId())     
        #self.button_calcpca.setEnabled(False)   
        vbox21.addWidget(self.button_calcpca)
        self.button_savepca = QtGui.QPushButton('Save PCA Results...')
        #.clicked.connect( self.OnSave, id=self.button_savepca.GetId())
        #self.button_savepca.setEnabled(False)
        vbox21.addWidget(self.button_savepca)
        
        hbox21 = QtGui.QHBoxLayout()
        text1 = QtGui.QLabel(self)
        text1.setText('Number of significant components')
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
   

 
        self.text_pcaspec = QtGui.QLabel(self)
        #self.text_pcaspec.SetFont(self.com.font)
        self.text_pcaspec.setText("PCA spectrum ") 
        vbox3.addWidget(self.text_pcaspec)
                        

        self.pcaspecfig = Figure((PlotW, PlotH))
        self.PCASpecPan = FigureCanvas(self.pcaspecfig)
        self.PCASpecPan.setParent(self)
        vbox3.addWidget(self.PCASpecPan)
        
    
        
        #panel 4
        vbox4 = QtGui.QVBoxLayout()
  
         
        text4 = QtGui.QLabel(self)
        text4.setText("PCA eigenvalues ")        
        vbox4.addWidget(text4)
        
        self.pcaevalsfig = Figure((PlotW, PlotH*0.75))
        self.PCAEvalsPan = FigureCanvas(self.pcaevalsfig)
        self.PCAEvalsPan.setParent(self)
        vbox4.addWidget(self.PCAEvalsPan)
                 

    
        
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
        self.add_scale_cb.toggle()
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
        
        limitevwin = LimitEv(self, self.com, self.stk)
        limitevwin.show()

#----------------------------------------------------------------------        
    def OnCliptoSubregion(self, evt):    
        clipwin = CliptoSubregion(self, self.com, self.stk)
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
                
                            
                fig = self.SpectrumPanel.get_figure()
                fig.savefig(fileName_spec)

            if img_png:
                
                fileName_img = self.SaveFileName+"_" +str(self.stk.ev[self.iev])+"eV."+ext
                fig = self.AbsImagePanel.get_figure()
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
            
                fig = self.SpectrumPanel.get_figure()
                fig.savefig(fileName_spec)

            if img_pdf:
                fileName_img = self.SaveFileName+"_" +str(self.stk.ev[self.iev])+"eV."+ext
                fig = self.AbsImagePanel.get_figure()
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
        imgregwin = ImageRegistration(self, self.com, self.stk)
        imgregwin.Show()

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
 
         
        matplotlib.rcParams['font.size'] = self.fontsize
 
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
        
        matplotlib.rcParams['font.size'] = self.fontsize

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

        self.stk = stack
        
        self.resize(400, 500)
        self.setWindowTitle('Stack File List')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)
                


#---------------------------------------------------------------------- 
class LimitEv(QtGui.QDialog):

    def __init__(self, parent,  common, stack):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        self.stack = stack
        self.com = common  
        
        self.resize(400, 500)
        self.setWindowTitle('Stack File List')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)


#---------------------------------------------------------------------- 
class CliptoSubregion(QtGui.QDialog):

    def __init__(self, parent,  common, stack):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        self.stack = stack
        self.com = common  
        
        self.resize(400, 500)
        self.setWindowTitle('Stack File List')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)
 
 
#---------------------------------------------------------------------- 
class ImageRegistration(QtGui.QDialog):

    def __init__(self, parent,  common, stack):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        self.stack = stack
        self.com = common  
        
        self.resize(400, 500)
        self.setWindowTitle('Stack File List')
        
        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)    
        
        
#---------------------------------------------------------------------- 
class SpectralROI(QtGui.QDialog):

    def __init__(self, parent,  common, stack):    
        QtGui.QWidget.__init__(self, parent)
        
        self.parent = parent

        self.stack = stack
        self.com = common  
        
        self.resize(400, 500)
        self.setWindowTitle('Stack File List')
        
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
        self.setWindowTitle('Stack File List')
        
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
        self.setWindowTitle('Stack File List')
        
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
        self.setWindowTitle('Stack File List')
        
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
                   

        #self.refresh_widgets()
        
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
            #self.page2.button_calcpca.setEnabled(False)
            #self.page4.button_loadtspec.setEnabled(False)
            #self.page4.button_addflat.setEnabled(False)
        else:
            self.page1.button_limitev.setEnabled(True)
            self.page1.button_subregion.setEnabled(True)
            self.page1.button_showi0.setEnabled(True)
            self.page1.rb_flux.setEnabled(True)
            self.page1.rb_od.setEnabled(True)   
            #self.page2.button_calcpca.setEnabled(True) 
            #self.page4.button_loadtspec.setEnabled(True)
            #self.page4.button_addflat.setEnabled(True)   
             
             
             
#         if self.common.pca_calculated == 0:      
#             self.page2.button_savepca.setEnabled(False)
#             self.page2.slidershow.setEnabled(False) 
#             self.page3.button_calcca.setEnabled(False)
#             self.page4.rb_fit.setEnabled(False)
#             if self.page7: self.page7.button_calckeng.setEnabled(False)
#         else:
#             self.page2.button_savepca.setEnabled(True)
#             self.page2.slidershow.setEnabled(True)
#             self.page3.button_calcca.setEnabled(True)  
#             self.page4.rb_fit.setEnabled(True)  
#             if self.page7: self.page7.button_calckeng.setEnabled(True)       
#             
#         if self.common.cluster_calculated == 0:   
#             self.page3.button_scatterplots.setEnabled(False)
#             self.page3.button_savecluster.setEnabled(False)
#             self.page3.slidershow.setEnabled(False)
#             self.page4.button_addclspec.setEnabled(False)
#         else:
#             self.page3.button_scatterplots.setEnabled(True)
#             self.page3.button_savecluster.setEnabled(True)  
#             self.page3.slidershow.setEnabled(True)
#             self.page4.button_addclspec.setEnabled(True)
#             
#         if self.common.spec_anl_calculated == 0:
#             self.page4.button_removespec.setEnabled(False)
#             self.page4.button_movespdown.setEnabled(False)
#             self.page4.button_movespup.setEnabled(False)
#             self.page4.button_save.setEnabled(False)
#             self.page4.button_showrgb.setEnabled(False) 
#         else:
#             self.page4.button_removespec.setEnabled(True)
#             self.page4.button_movespdown.setEnabled(True)
#             self.page4.button_movespup.setEnabled(True)
#             self.page4.button_save.setEnabled(True) 
#             self.page4.button_showrgb.setEnabled(True)          
#                   
#             
#         self.page1.ResetDisplaySettings()
#             
#             
            
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
#             self.page7.button_calckeng.setEnabled(False)
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