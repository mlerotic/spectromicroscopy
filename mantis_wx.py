# 
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
# 
#   Copyright (C) 2011 Mirna Lerotic, 2nd Look
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


from __future__ import division
import wx
import wx.lib.intctrl
import wx.lib.masked.numctrl
import matplotlib as mtplot 
mtplot.interactive( True )
mtplot.use( 'WXAgg', warn=False )
mtplot.rcParams['svg.fonttype'] = 'none'

import wxmpl
import numpy as npy
import os.path
from mpl_toolkits.axes_grid import make_axes_locatable
import time
import sys
import getopt

import data_struct
import data_stack
import analyze
import logos
#import nnma
import henke

from helpers import resource_path

Winsizex = 1000
Winsizey = 740

ImgDpi = 40

PlotH = 3.46
PlotW = PlotH*1.61803

verbose = False

defaultDir = ''

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

        self.font = wx.SystemSettings_GetFont(wx.SYS_DEFAULT_GUI_FONT)


""" ------------------------------------------------------------------------------------------------"""
class PageKeyEng(wx.Panel):
    def __init__(self, parent, common, data_struct, stack, anlz):
        wx.Panel.__init__(self, parent)
        
        self.SetBackgroundColour("White")
        
        self.com = common 
        self.data_struct = data_struct
        self.stk = stack       
        self.anlz = anlz
        
        self.selica = 1       
        self.numica = 2
        
        
        pw = PlotW*0.8
        ph = PlotH*0.8
        
          
        self.fontsize = self.com.fontsize
        
        self.i_eng = 0
        self.keyengs_calculated = 0
       

          
        #panel 1        
        panel1 = wx.Panel(self, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_1 = wx.TextCtrl(panel1, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_1.SetFont(self.com.font)
        self.tc_1.SetValue("Average Optical Density")        
   
        i1panel = wx.Panel(panel1, -1, style = wx.SUNKEN_BORDER)
        self.KESpecPan = wxmpl.PlotPanel(i1panel, -1, size =(pw, ph), cursor=False, crosshairs=False, location=False, zoom=False)
        wxmpl.EVT_POINT(i1panel, self.KESpecPan.GetId(), self.OnPointSpectrum)
        
        vbox1.Add(self.tc_1, 1, wx.EXPAND )        
        vbox1.Add(i1panel, 0)

        panel1.SetSizer(vbox1)
        
        
                
        #panel 2
        panel2 = wx.Panel(self, -1)
        vbox2 = wx.BoxSizer(wx.VERTICAL)
                
        sizer1 = wx.StaticBoxSizer(wx.StaticBox(panel2, -1, 'Key Energies Analysis'), wx.VERTICAL)
        self.button_calckeng = wx.Button(panel2, -1, 'Find Key Energies', (90, 10))
        self.button_calckeng.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnCalcKeyEng, id=self.button_calckeng.GetId())     
        self.button_calckeng.Disable()   
        sizer1.Add(self.button_calckeng, 0, wx.EXPAND|wx.LEFT|wx.RIGHT, 5)
        self.button_save = wx.Button(panel2, -1, 'Save Results...', (90, 10))
        self.button_save.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnSave, id=self.button_save.GetId())
        self.button_save.Disable()
        
        
        hbox21 = wx.BoxSizer(wx.HORIZONTAL)
        text1 = wx.StaticText(panel2, label="Threshold")
        self.tc_keyengthresh = wx.lib.masked.numctrl.NumCtrl(panel2, -1, 
                                              value = float(0.1),  
                                              integerWidth = 5,
                                              fractionWidth = 2, min = 0, max = 5,
                                              limited = False)

        hbox21.Add(text1, 0, wx.TOP, 20)
        hbox21.Add((10,0))
        hbox21.Add(self.tc_keyengthresh, 0, wx.TOP, 15)            
        sizer1.Add(hbox21, 0, wx.EXPAND|wx.LEFT|wx.RIGHT, 5)    
              
        sizer1.Add((0,10))        
        sizer1.Add(wx.StaticLine(panel2), 0, wx.ALL|wx.EXPAND, 5)        
        sizer1.Add((0,10))  
        
        sizer1.Add(self.button_save, 0, wx.EXPAND|wx.LEFT|wx.RIGHT, 5)        
        
        vbox2.Add(sizer1,0)

        panel2.SetSizer(vbox2)


        #panel 3
        panel3 = wx.Panel(self, -1)
        vbox3 = wx.BoxSizer(wx.VERTICAL)
   
        t1 = wx.StaticText(panel3, label="Key Energies")

        
        
        self.lc_1 = wx.ListCtrl(panel3, -1, size = (200, 390), style = wx.LC_REPORT|wx.LC_NO_HEADER)
        self.lc_1.InsertColumn(0, 'KEng')
        self.lc_1.SetColumnWidth(0, 160)  
        
        vbox3.Add(t1, 0)
        vbox3.Add(self.lc_1, 0,  wx.EXPAND)
        panel3.SetSizer(vbox3)
            
        self.Bind(wx.EVT_LIST_ITEM_SELECTED , self.OnEngListClick, self.lc_1)   
        
        #panel 4     
        panel4 = wx.Panel(self, -1)
        vbox4 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_imageeng = wx.TextCtrl(panel4, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE, size=(160,-1))
        self.tc_imageeng.SetFont(self.com.font)
        self.tc_imageeng.SetValue("Image at key energy: ")
        
        hbox41 = wx.BoxSizer(wx.HORIZONTAL)
   
        i1panel = wx.Panel(panel4, -1, style = wx.SUNKEN_BORDER)
        self.AbsImagePanel = wxmpl.PlotPanel(i1panel, -1, size =(PlotH*.8, PlotH*.8), cursor=False, crosshairs=False, location=False, zoom=False)
        
        vbox41 = wx.BoxSizer(wx.VERTICAL)
        self.slider_eng = wx.Slider(panel4, -1, self.i_eng, 0, 100, style=wx.SL_LEFT|wx.SL_VERTICAL)        
        self.slider_eng.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnScrollEng, self.slider_eng)

        self.engspin = wx.SpinButton(panel4, -1, size = ((8,-1)), style=wx.SP_ARROW_KEYS)
        self.Bind(wx.EVT_SPIN_UP, self.OnEngspinUp, self.engspin)
        self.Bind(wx.EVT_SPIN_DOWN, self.OnEngspinDown, self.engspin)
        

        hbox41.Add(i1panel, 0)
        vbox41.Add((0,3))
        vbox41.Add(self.slider_eng, 1,  wx.EXPAND) 
        vbox41.Add(self.engspin, 0,  wx.EXPAND)         
        hbox41.Add(vbox41, 0,  wx.EXPAND) 
        
        vbox4.Add(self.tc_imageeng, 0)        
        vbox4.Add(hbox41, 0)

        panel4.SetSizer(vbox4)
        

        vboxtop1 = wx.BoxSizer(wx.VERTICAL)
                     
        vboxtop1.Add((0,40))
        vboxtop1.Add(panel2, 0, wx.LEFT, 40)
        vboxtop1.Add((0,40))
        vboxtop1.Add(panel3, 0, wx.LEFT, 40)
         
        vboxtop2 = wx.BoxSizer(wx.VERTICAL)          
        vboxtop2.Add((0,20))
        vboxtop2.Add(panel1, 0, wx.LEFT, 40)
        vboxtop2.Add((0,20))
        vboxtop2.Add(panel4, 0, wx.LEFT, 40)
                
        hboxtop = wx.BoxSizer(wx.HORIZONTAL)        
        hboxtop.Add(vboxtop1)
        hboxtop.Add((20,0))
        hboxtop.Add(vboxtop2)
        
        
        self.SetSizer(hboxtop) 
        
        
#----------------------------------------------------------------------
    def OnCalcKeyEng(self, event):
        

        wx.BeginBusyCursor()
        
        threshold = self.tc_keyengthresh.GetValue()
  
        self.keyenergies = []
        
        self.keyenergies = self.anlz.calc_key_engs(threshold)
        
        if len(self.keyenergies) > 0:        
            self.keyengs_calculated = 1
        
            self.i_eng = 0
            self.slider_eng.SetRange(0,len(self.keyenergies)-1)
            self.slider_eng.SetValue(self.i_eng)

            self.ShowPlots()
            self.ShowImage()
            self.ShowKEngs()
            
            self.button_save.Enable()
            
        else:
            self.ShowPlots()
            fig = self.AbsImagePanel.get_figure()
            fig.clf()
            self.AbsImagePanel.draw()               
            self.lc_1.DeleteAllItems()
            
            self.button_save.Disable()
            
        wx.EndBusyCursor() 

        
#----------------------------------------------------------------------            
    def OnScrollEng(self, event):
        self.i_eng = event.GetInt()

        if self.keyengs_calculated == 1:
            self.ShowImage()
            self.ShowKEngs()
            
#----------------------------------------------------------------------            
    def OnEngspinUp(self, event):
        if (self.keyengs_calculated == 1) and (self.i_eng > 0):
            self.i_eng = self.i_eng - 1
            self.slider_eng.SetValue(self.i_eng)

            self.ShowImage()
            self.ShowKEngs()

            
#----------------------------------------------------------------------            
    def OnEngspinDown(self, event):
        if (self.keyengs_calculated == 1) and (self.i_eng < len(self.keyenergies)-1):
            self.i_eng = self.i_eng + 1
            self.slider_eng.SetValue(self.i_eng)

            self.ShowImage()
            self.ShowKEngs()

#----------------------------------------------------------------------  
    def OnPointSpectrum(self, evt):
        x = evt.xdata
        y = evt.ydata
        
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
            
            self.slider_eng.SetValue(self.i_eng)
                              
#----------------------------------------------------------------------     
    def ShowPlots(self):


        odtotal = self.stk.od3d.sum(axis=0)   
        odtotal = odtotal.sum(axis=0)/(self.stk.n_rows*self.stk.n_cols)      
        odtotal /= odtotal.max()/0.7

        
        fig = self.KESpecPan.get_figure()
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        
        mtplot.rcParams['font.size'] = self.fontsize
        
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
    def OnEngListClick(self, event):
        
                
        self.i_eng = event.m_itemIndex
        
        event.Skip()
        
        self.ShowKEngs()     
        self.ShowImage()
         
#----------------------------------------------------------------------           
    def ShowKEngs(self):    
        
        self.lc_1.DeleteAllItems()
        
        for i in range(len(self.keyenergies)):
            self.lc_1.InsertStringItem(i,'{0:08.2f}'.format(self.keyenergies[i]))
            
        self.lc_1.SetItemBackgroundColour(self.i_eng, 'light blue')

#----------------------------------------------------------------------        
    def ShowImage(self):
        
        iev=(npy.abs(self.stk.ev-self.keyenergies[self.i_eng])).argmin() 
        image = self.stk.absdata[:,:,int(iev)].copy() 

        fig = self.AbsImagePanel.get_figure()
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        axes = fig.gca()
        fig.patch.set_alpha(1.0)
        
        im = axes.imshow(image, cmap=mtplot.cm.get_cmap("gray")) 
        
        #Show Scale Bar
        startx = int(self.stk.n_rows*0.05)
        starty = self.stk.n_cols-int(self.stk.n_cols*0.05)-self.stk.scale_bar_pixels_y
        um_string = ' $\mathrm{\mu m}$'
        microns = '$'+self.stk.scale_bar_string+' $'+um_string
        axes.text(self.stk.scale_bar_pixels_x+startx+1,starty+1, microns, horizontalalignment='left', verticalalignment='center',
                  color = 'black', fontsize=14)
        #Matplotlib has flipped scales so I'm using rows instead of cols!
        p = mtplot.patches.Rectangle((startx,starty), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                               color = 'black', fill = True)
        axes.add_patch(p)
            
       
        axes.axis("off")      
        self.AbsImagePanel.draw()
        
        self.tc_imageeng.SetValue('Image at key energy: {0:5.2f} eV'.format(float(self.stk.ev[iev])))
 
#----------------------------------------------------------------------
    def OnSave(self, event):


        #Save images
                       
        fileName = wx.FileSelector('Save OD Plot with Key Energies', default_extension='png', 
                                   wildcard=('Portable Network Graphics (*.png)|*.png|' 
                                             + 'Adobe PDF Files (*.pdf)|*.pdf|All files (*.*)|*.*'), 
                                              parent=self, flags=wx.SAVE|wx.OVERWRITE_PROMPT) 
   
        if not fileName: 
            return 

        path, ext = os.path.splitext(fileName) 
        ext = ext[1:].lower() 
       
        if ext != 'png' and ext != 'pdf': 
            error_message = ( 
                  'Only the PNG and PDF image formats are supported.\n' 
                 'A file extension of `png\' or `pdf\' must be used.') 
            wx.MessageBox(error_message, 'Error - Could not save file.', 
                  parent=self, style=wx.OK|wx.ICON_ERROR) 
            return 
   
        try: 

            mtplot.rcParams['pdf.fonttype'] = 42
            
            fig = self.KESpecPan.get_figure()
            fig.savefig(fileName)

            
        except IOError, e:
            if e.strerror:
                err = e.strerror 
            else: 
                err = e 
   
            wx.MessageBox('Could not save file: %s' % err, 'Error', 
                          parent=self, style=wx.OK|wx.ICON_ERROR) 
            
            
        #Save text file with list of energies
        textfilepath = path+'_keyenergies.csv'
        file = open(textfilepath, 'w')
        print>>file, '*********************  Key Energies  ********************'
        for i in range(len(self.keyenergies)):
            print>>file, '%.6f' %(self.keyenergies[i])
        
        file.close()        
            
        return
          
 
""" ------------------------------------------------------------------------------------------------"""
class PageICA(wx.Panel):
    def __init__(self, parent, common, data_struct, stack, anlz):
        wx.Panel.__init__(self, parent)
        
        self.SetBackgroundColour("White")
        
        self.com = common 
        self.data_struct = data_struct
        self.stk = stack       
        self.anlz = anlz
        
        self.selica = 1       
        self.numica = 2
        
        
        pw = PlotW*0.97
        ph = PlotH*0.97
        
          
        self.fontsize = self.com.fontsize
        

          
        #panel 1        
        panel1 = wx.Panel(self, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_ICAsp = wx.TextCtrl(panel1, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_ICAsp.SetFont(self.com.font)
        self.tc_ICAsp.SetValue("ICA spectrum ")
        
        hbox11 = wx.BoxSizer(wx.HORIZONTAL)
   
        i1panel = wx.Panel(panel1, -1, style = wx.SUNKEN_BORDER)
        self.ICASpecPan = wxmpl.PlotPanel(i1panel, -1, size =(pw, ph), cursor=False, crosshairs=False, location=False, zoom=False)
               
        vbox11 = wx.BoxSizer(wx.VERTICAL)               
        self.slidershow = wx.Slider(panel1, -1, self.selica, 1, 20, style=wx.SL_LEFT|wx.SL_VERTICAL)
        #self.slidershow.Disable()          
        self.slidershow.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnICAScroll, self.slidershow)
        

        self.icaspin = wx.SpinButton(panel1, -1, size = ((8,-1)), style=wx.SP_ARROW_KEYS)
        self.Bind(wx.EVT_SPIN_UP, self.OnICASpinUp, self.icaspin)
        self.Bind(wx.EVT_SPIN_DOWN, self.OnICASpinDown, self.icaspin)
        
        hbox11.Add(i1panel, 0)
        vbox11.Add((0,3))
        vbox11.Add(self.slidershow, 1,  wx.EXPAND) 
        vbox11.Add(self.icaspin, 0,  wx.EXPAND)         
        hbox11.Add(vbox11, 0,  wx.EXPAND) 

        
        vbox1.Add(self.tc_ICAsp, 1, wx.EXPAND )        
        vbox1.Add(hbox11, 0)

        panel1.SetSizer(vbox1)
        
        
                
        #panel 2
        panel2 = wx.Panel(self, -1)
        vbox2 = wx.BoxSizer(wx.VERTICAL)
                
        sizer1 = wx.StaticBoxSizer(wx.StaticBox(panel2, -1, 'Independent Component Analysis'), wx.VERTICAL)
        self.button_calcica = wx.Button(panel2, -1, 'Calculate ICA', (80, 10))
        self.button_calcica.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnCalcICA, id=self.button_calcica.GetId())     
        #self.button_calcica.Disable()   
        sizer1.Add(self.button_calcica, 0, wx.EXPAND|wx.LEFT|wx.RIGHT, 20)
        self.button_saveica = wx.Button(panel2, -1, 'Save ICA Results...', (80, 10))
        self.button_saveica.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnSave, id=self.button_saveica.GetId())
        self.button_saveica.Disable()
        sizer1.Add(self.button_saveica, 0, wx.EXPAND|wx.LEFT|wx.RIGHT, 20)
        
#        hbox21 = wx.BoxSizer(wx.HORIZONTAL)
#        text1 = wx.StaticText(panel2, -1, 'Number of independent components',  style=wx.ALIGN_LEFT)
#        text1.SetFont(self.com.font)
#        self.nicaspin = wx.SpinCtrl(panel2, -1, '',  size= (60, -1), style=wx.ALIGN_LEFT)
#        self.nicaspin.SetRange(1,20)
#        #self.Bind(wx.EVT_SPINCTRL, self.OnNPCAspin, self.nicaspin)
#        hbox21.Add(text1, 0, wx.TOP, 20)
#        hbox21.Add((10,0))
#        hbox21.Add(self.nicaspin, 0, wx.TOP, 15)            
#        sizer1.Add(hbox21, 0, wx.EXPAND)    
              
        
        vbox2.Add(sizer1,0)

        panel2.SetSizer(vbox2)

      
        
#        #panel 4
#        panel4 = wx.Panel(self, -1)
#             
#        i4panel = wx.Panel(panel4, -1, style = wx.SUNKEN_BORDER)
#        self.ICAEvalsPan = wxmpl.PlotPanel(i4panel, -1, size =(pw, ph*0.75), cursor=False, crosshairs=False, location=False, zoom=False)
#        #wxmpl.EVT_POINT(i4panel, self.ICAEvalsPan.GetId(), self.OnPointEvalsImage)   
#        
#        vbox4 = wx.BoxSizer(wx.VERTICAL)
#        text4 = wx.TextCtrl(panel4, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
#        text4.SetFont(self.com.font)
#        text4.SetValue("ICA eigenvalues ??")        
#        vbox4.Add(text4, 0)
#        vbox4.Add(i4panel, 0)
#        
#        panel4.SetSizer(vbox4)      
        

        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
#        gridtop = wx.FlexGridSizer(1, 2, vgap=10, hgap=20)
#        gridtop.Add(panel2, 0, wx.LEFT|wx.TOP, 20)
#        gridtop.Add(panel4, 0, wx.ALIGN_LEFT)
              
        vboxtop.Add((0,40))
        vboxtop.Add(panel2, 0, wx.LEFT, 40)
        vboxtop.Add((0,40))
        vboxtop.Add(panel1, 0, wx.LEFT, 40)
         
        
        self.SetSizer(vboxtop) 
        
        
#----------------------------------------------------------------------
    def OnCalcICA(self, event):
        
        #self.anlz.calculate_fastica(self.stk.od, 4)
        
        try:
            import mdp
        except:
            print 'ERROR: Could not find MDP library.'
            return
            
        if self.com.cluster_calculated == 0:
            print 'ERROR: Calculate cluster spectra before ICA.'
            return

        wx.BeginBusyCursor()
        self.calcica = 0   
        self.selpca = 1   
        self.slidershow.SetValue(self.selpca)
 
        #Use cluster spectra
        X = self.anlz.clusterspectra.T
        
        scrollmax = npy.min([self.anlz.clusterspectra.shape[0], 20])
        self.slidershow.SetMax(scrollmax)        
        
        try:

            #ica = mdp.nodes.FastICANode( verbose=True)
            ica = mdp.nodes.CuBICANode(limit=0.0001, verbose=True)
            ica.train(X)
            comp = ica.execute(X)
        
            self.icasig = comp        
            self.recica = npy.dot(comp,ica.get_recmatrix())        
              
            self.com.ica_calculated = 1
            self.showICASpectrum()
            wx.EndBusyCursor() 
            self.button_saveica.Enable()
        except:
            self.com.ica_calculated = 0
            wx.EndBusyCursor()
            wx.MessageBox("ICA not calculated.")
        
        wx.GetApp().TopWindow.refresh_widgets()

                 
#----------------------------------------------------------------------        
    def OnICAScroll(self, event):
        self.sel = event.GetInt()
        self.selica = self.sel
        if self.com.ica_calculated == 1:
            self.showICASpectrum()       
 
#----------------------------------------------------------------------            
    def OnICASpinUp(self, event):
        if (self.com.ica_calculated == 1) and (self.selica > 1):
            self.selica = self.selica - 1
            self.slidershow.SetValue(self.selica)

            self.showICASpectrum() 

            
#----------------------------------------------------------------------            
    def OnICASpinDown(self, event):
        if (self.com.ica_calculated == 1) and (self.selica < self.icasig.shape[1]-1):
            self.selica = self.selica + 1
            self.slidershow.SetValue(self.selica) 
            
            self.showICASpectrum() 
            
#----------------------------------------------------------------------
    def OnSave(self, event):


        try: 
            wildcard = "PNG files (*.png)|*.png"
            dialog = wx.FileDialog(None, "Save as .png", wildcard=wildcard,
                                    style=wx.SAVE|wx.OVERWRITE_PROMPT)

            if dialog.ShowModal() == wx.ID_OK:
                filepath = dialog.GetPath()
                
                dialog.Destroy()
                            
            wx.BeginBusyCursor()   
                           
            from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
            mtplot.rcParams['pdf.fonttype'] = 42
            
            ext = 'png'
            suffix = "." + ext    
            
            for i in range(self.icasig.shape[1]):
            
                icaspectrum = self.icasig[:, i]
                fig = mtplot.figure.Figure(figsize =(PlotW, PlotH))
                canvas = FigureCanvas(fig)
                fig.add_axes((0.15,0.15,0.75,0.75))
                axes = fig.gca()
                mtplot.rcParams['font.size'] = self.fontsize
                specplot = axes.plot(self.stk.ev, icaspectrum)    
                axes.set_xlabel('Photon Energy [eV]')
                axes.set_ylabel('Optical Density')
            
                fileName_spec = filepath+"_ICAspectrum_" +str(i+1)+"."+ext
                fig.savefig(fileName_spec)  
             
            wx.EndBusyCursor()      

        except:

            wx.EndBusyCursor()
            wx.MessageBox("Could not save .png file.")
                   
        
        return
                   
#----------------------------------------------------------------------
    def CalcICA(self):
        pass

        

#----------------------------------------------------------------------     
    def showICASpectrum(self):

        if self.com.ica_calculated == 0:
            return
        
        self.icaspectrum = self.icasig[:, self.selica-1]
        #self.icaspectrum = self.recica[:, self.selica-1]
        
            
        
        fig = self.ICASpecPan.get_figure()
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        
        mtplot.rcParams['font.size'] = self.fontsize

        specplot = axes.plot(self.stk.ev,self.icaspectrum)
                        
        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Optical Density')
        
        self.ICASpecPan.draw()
 
 

""" ------------------------------------------------------------------------------------------------"""
class PageNNMA(wx.Panel):

    def __init__(self, parent, common, data_struct, stack, nnma):

        wx.Panel.__init__(self, parent)

        self.data_struct = data_struct
        self.stk = stack
        self.com = common
        self.nnma = nnma
        self.kNNMA = self.nnma.kNNMA	# number of chemical components to look for in NNMA
        self.fontsize = self.com.fontsize        
        self.iev = 0
        #self.calcnnma = False
        

# ------------------------------------------------------------------------------------------------
# Subclass of PageNNMA to display optical density results
# ------------------------------------------------------------------------------------------------
class PageNNMAOptDensity(PageNNMA):

    def __init__(self, parent, common, data_struct, stack, nnma):

        PageNNMA.__init__(self, parent, common, data_struct, stack, nnma)

        self.data_struct = data_struct
        self.stk = stack
        self.com = common
        self.nnma = nnma
        self.kNNMA = self.nnma.kNNMA	# number of chemical components to look for in NNMA
        self.fontsize = self.com.fontsize        
        self.iev = 0
        
        # Panel 1: For users to enter NNMA parameters/settings -------------------------------------
        
        panel1 = wx.Panel(self, -1)
        sizer1 = wx.StaticBoxSizer(wx.StaticBox(panel1, -1, 'NNMA parameters'), wx.VERTICAL)
        
        # A spinner to allow users to set number of NNMA components
        hbox11 = wx.BoxSizer(wx.HORIZONTAL)
        text11 = wx.StaticText(panel1, -1, 'Number of chemical components: ', style=wx.ALIGN_LEFT)
        text11.SetFont(self.com.font)
        self.kNNMASpin = wx.SpinCtrl(panel1, -1, '',  size=(60, -1), style=wx.ALIGN_LEFT)
        self.kNNMASpin.SetRange(1,100)
        self.kNNMASpin.SetValue(5)
        self.Bind(wx.EVT_SPINCTRL, self.onkNNMASpin, self.kNNMASpin)
        if verbose: print("self.kNNMASpin.GetId() = ", self.kNNMASpin.GetId())
        hbox11.Add(text11, 0, wx.TOP, 20)
        hbox11.Add((10,0))
        hbox11.Add(self.kNNMASpin, 0, wx.TOP, 15)            
        sizer1.Add(hbox11, 0, wx.EXPAND)

        # A text field to allow users to set number of NNMA iterations
        hbox12 = wx.BoxSizer(wx.HORIZONTAL)
        text12 = wx.StaticText(panel1, -1, 'Number of iterations: ', style=wx.ALIGN_LEFT)
        text12.SetFont(self.com.font)
        self.itersCtrl = wx.TextCtrl(panel1, -1, style=wx.ALIGN_LEFT)
        self.itersCtrl.SetValue('100')
        hbox12.Add(text12, 1, wx.TOP, 15)
        hbox12.Add((10, 0))
        hbox12.Add(self.itersCtrl, 1, wx.TOP, 15)
        sizer1.Add((0, 10))
        sizer1.Add(hbox12, 0, wx.EXPAND)

        # ComboBox to display available NNMA algorithms
        algoNNMA = ['Basic', 'Smoothness', 'Sparsity']
        self.comboBoxAlgoNNMA = wx.ComboBox(panel1, -1, choices=algoNNMA, style=wx.CB_READONLY)
        self.comboBoxAlgoNNMA.SetValue(algoNNMA[0])
        self.Bind(wx.EVT_COMBOBOX, self.onSelectComboBoxAlgoNNMA, id=self.comboBoxAlgoNNMA.GetId())
        text12 = wx.StaticText(panel1, -1, 'NNMA algorithm: ', style=wx.ALIGN_LEFT)
        text12.SetFont(self.com.font)
        hbox12 = wx.BoxSizer(wx.HORIZONTAL)
        hbox12.Add(text12, 0, wx.TOP, 20)
        hbox12.Add((10, 0))
        sizer1.Add(hbox12, 0, wx.EXPAND)
        sizer1.Add(self.comboBoxAlgoNNMA, 0, wx.EXPAND)

        # A text field to allow users to set sparseness of t matrix
        hbox13 = wx.BoxSizer(wx.HORIZONTAL)
        text13 = wx.StaticText(panel1, -1, 'Sparseness (0 to 1): ', style=wx.ALIGN_LEFT)
        text13.SetFont(self.com.font)
        self.sparsenessCtrl = wx.TextCtrl(panel1, -1, style=wx.ALIGN_LEFT)
        self.sparsenessCtrl.SetValue('0.5')
        hbox13.Add(text13, 1, wx.TOP, 15)
        hbox13.Add((10, 0))
        hbox13.Add(self.sparsenessCtrl, 1, wx.TOP, 15)
        sizer1.Add((0, 10))
        sizer1.Add(hbox13, 0, wx.EXPAND)

        # ComboBox for choosing initial matrices
        initMatricesNNMA = ['Random', 'Cluster']
        self.comboBoxInitMatricesNNMA = wx.ComboBox(panel1, -1, choices=initMatricesNNMA, style=wx.CB_READONLY)
        self.comboBoxInitMatricesNNMA.SetValue(initMatricesNNMA[0])
        self.Bind(wx.EVT_COMBOBOX, self.onSelectComboBoxInitMatricesNNMA, id=self.comboBoxInitMatricesNNMA.GetId())
        text14 = wx.StaticText(panel1, -1, 'Initial matrices: ', style=wx.ALIGN_LEFT)
        text14.SetFont(self.com.font)
        hbox14 = wx.BoxSizer(wx.HORIZONTAL)
        hbox14.Add(text14, 0, wx.TOP, 20)
        hbox14.Add((10, 0))
        sizer1.Add(hbox14, 0, wx.EXPAND)
        sizer1.Add(self.comboBoxInitMatricesNNMA, 0, wx.EXPAND)
        
        # Button to calculate NNMA
        self.button_calcNNMA = wx.Button(panel1, -1, 'Calculate NNMA', (10, 10))
        self.button_calcNNMA.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.onCalcNNMA, id=self.button_calcNNMA.GetId())
        #self.button_calcNNMA.Disable()
        sizer1.Add((0, 10))
        sizer1.Add(self.button_calcNNMA, 0, wx.EXPAND)

        # Button for testing
        self.button_test = wx.Button(panel1, -1, 'Print test', (10, 10))
        self.button_test.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnTestButton, id=self.button_test.GetId())
        sizer1.Add(self.button_test, 0, wx.EXPAND)
        
        panel1.SetSizer(sizer1)
        
        # Panel 2: Display optical density reconstruction results 
        panel2 = wx.Panel(self, -1)
        sizer21 = wx.StaticBoxSizer(wx.StaticBox(panel2, -1, 'Optical density'), wx.VERTICAL)
        
        scaleFactorW = 0.6		# factor to scale plot width
        scaleFactorH = 0.7		# factor to scale plot height
        
        # hbox21 contains the original and reconstructed optical density
        hbox21 = wx.BoxSizer(wx.HORIZONTAL)
        vbox21 = wx.BoxSizer(wx.VERTICAL)
        self.ODPanel = wxmpl.PlotPanel(panel2, -1, size=(PlotW*scaleFactorW, PlotH*scaleFactorH), cursor=False, crosshairs=False, location=False, zoom=False)
        textOD = wx.TextCtrl(panel2, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE, size=(100, 20))
        textOD.SetFont(self.com.font)
        textOD.SetValue("Original")
        vbox21.Add(textOD, 0, wx.EXPAND)
        vbox21.Add(self.ODPanel, 0, wx.TOP)
        vbox22 = wx.BoxSizer(wx.VERTICAL)
        self.ODReconPanel = wxmpl.PlotPanel(panel2, -1, size=(PlotW*scaleFactorW, PlotH*scaleFactorH), cursor=False, crosshairs=False, location=False, zoom=False)
        textODRecon = wx.TextCtrl(panel2, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        textODRecon.SetFont(self.com.font)
        textODRecon.SetValue('Reconstructed')
        vbox22.Add(textODRecon, 0, wx.EXPAND)
        vbox22.Add(self.ODReconPanel, 0, wx.TOP)
        hbox21.Add(vbox21)
        hbox21.Add(vbox22)
        
        # vbox23 contains the difference between original and reconstructed optical density
        vbox23 = wx.BoxSizer(wx.VERTICAL)
        self.ODErrorPanel = wxmpl.PlotPanel(panel2, -1, size=(PlotW*scaleFactorW, PlotH*scaleFactorH), cursor=False, crosshairs=False, location=False, zoom=False)
        textODError = wx.TextCtrl(panel2, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        textODError.SetFont(self.com.font)
        textODError.SetValue("% Difference")
        vbox23.Add(textODError, 0, wx.EXPAND)
        vbox23.Add(self.ODErrorPanel, 0, wx.TOP)
        
        flexGridSizer = wx.FlexGridSizer(2, 1, vgap=10, hgap=20)
        flexGridSizer.Add(hbox21)
        #flexGridSizer.Add((0, 20))
        flexGridSizer.Add(vbox23)
        panel2.SetSizer(flexGridSizer)
        
        sizer21.Add(panel2)
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(panel1, 0)
        hbox.Add(panel2, 0)
        self.SetSizer(hbox)

#----------------------------------------------------------------------
    def OnTestButton(self, event):

        self.nnma.printTest()

#----------------------------------------------------------------------
    def onkNNMASpin(self, event):
        
        self.nnma.kNNMA = event.GetInt()
        if verbose: print("nnma.kNNMA = ", self.nnma.kNNMA)

#----------------------------------------------------------------------
    def onSelectComboBoxAlgoNNMA(self, event):

        self.nnma.algoNNMA = event.GetString()
        if verbose: print("nnma.algoNNMA = ", self.nnma.algoNNMA)

#----------------------------------------------------------------------
    def onSelectComboBoxInitMatricesNNMA(self, event):

        self.nnma.initMatrices = event.GetString()
        if verbose: print("nnma.initMatrices = ", self.nnma.initMatrices)

#----------------------------------------------------------------------
    def onScrollEnergy(self, event):

        self.iev = event.GetInt()
        if self.com.stack_loaded == 1:
            self.loadODImage()

#----------------------------------------------------------------------
    def calcNNMA(self):

        self.kNNMA = self.nnma.kNNMA
        self.algoNNMA = self.nnma.algoNNMA
        try:
            self.nnma.maxIters = int(self.itersCtrl.GetValue())
        except:
            wx.MessageBox("Please enter valid number of iterations.")
        try:
            self.nnma.sparsenessT = float(self.sparsenessCtrl.GetValue())
        except:
            wx.MessageBox("Please enter sparseness parameter between 0 and 1.")
        self.nnma.calcNNMA(self.algoNNMA, kComponents=self.kNNMA)

#----------------------------------------------------------------------
    def onCalcNNMA(self, event):

        wx.BeginBusyCursor()
        PageNNMA.calcnnma = False
        #self.calcnnma = False		# boolean for whether NNMA has been calculated

        try:
            self.calcNNMA()
            #self.calcnnma = True
            PageNNMA.calcnnma = True
            self.loadODImage()
            self.loadODReconImage()
            self.loadODErrorImage()
            wx.EndBusyCursor()

        except:
            wx.EndBusyCursor()
            wx.MessageBox("NNMA not calculated.")

        #wx.GetApp().TopWindow.refresh_widgets() 	# <-- Is this necessary??

#----------------------------------------------------------------------
    def loadODImage(self):

        self.show_colorbar = 1
        self.show_scale_bar = 0
        
        image = self.nnma.OD[self.iev, :, :] 

        fig = self.ODPanel.get_figure()
        fig.clf()
        
        if self.show_colorbar == 0:
            fig.add_axes([0.15, 0.15, 0.8, 0.8])
            axes = fig.gca()
        else:
            axes = fig.gca()
            divider = make_axes_locatable(axes)
            axcb = divider.new_horizontal(size="3%", pad=0.03)  
            fig.add_axes(axcb)
            axes.set_position([0.15, 0.15, 0.8, 0.8]) 	# <-- Is this necessary?
        
        fig.patch.set_alpha(1.0) 	# set figure transparency
        
        self.defaultdisplay = 0
        self.colortable = "gray"
        if self.defaultdisplay == 1.0:
            im = axes.imshow(image, cmap=mtplot.cm.get_cmap(self.colortable)) 
        else:
            imgmax = npy.amax(image)
            imgmin = npy.amin(image)
            if (imgmin < 0.0):
                image = (image-imgmin)/(imgmax-imgmin)
                imgmax = 1.0
                imgmin = 0.0
            self.brightness_min = 0.0
            self.brightness_max = 1.0
            im = axes.imshow(image, cmap=mtplot.cm.get_cmap(self.colortable), vmin=(imgmin+imgmax*self.brightness_min),vmax=imgmax*self.brightness_max)
        
        if self.show_colorbar == 1:
            cbar = axes.figure.colorbar(im, orientation='vertical', cax=axcb) 
        
        if self.show_scale_bar == 1:
            startx = int(self.stk.n_rows*0.05)
            starty = self.stk.n_cols-int(self.stk.n_cols*0.05)-self.stk.scale_bar_pixels_y
            um_string = ' $\mathrm{\mu m}$'
            microns = '$'+self.stk.scale_bar_string+' $'+um_string
            axes.text(self.stk.scale_bar_pixels_x+startx+1,starty+1, microns, horizontalalignment='left', verticalalignment='center', color='black', fontsize=14)
            #Matplotlib has flipped scales so I'm using rows instead of cols!
            p = mtplot.patches.Rectangle((startx,starty), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y, color = 'black', fill = True)
            axes.add_patch(p)
        
        self.ODPanel.draw()
        
        #self.tc_imageeng.SetValue("Image at energy: {0:5.2f} eV".format(float(self.stk.ev[self.iev])))
 
#----------------------------------------------------------------------
    def loadODReconImage(self):

        self.show_colorbar = 1
        self.show_scale_bar = 0
        
        image = self.nnma.ODRecon[self.iev, :, :] 	# <--- C-style row/col major issue: transpose()?

        fig = self.ODReconPanel.get_figure()
        fig.clf()
        
        if self.show_colorbar == 0:
            fig.add_axes([0.15, 0.15, 0.8, 0.8])
            axes = fig.gca()
        else:
            axes = fig.gca()
            divider = make_axes_locatable(axes)
            axcb = divider.new_horizontal(size="3%", pad=0.03)  
            fig.add_axes(axcb)
            axes.set_position([0.15, 0.15, 0.8, 0.8])
        
        fig.patch.set_alpha(1.0) 	# set figure transparency
        
        self.defaultdisplay = 0
        self.colortable = "gray"
        if self.defaultdisplay == 1.0:
            im = axes.imshow(image, cmap=mtplot.cm.get_cmap(self.colortable)) 
        else:
            imgmax = npy.amax(image)
            imgmin = npy.amin(image)
            if (imgmin < 0.0):
                image = (image-imgmin)/(imgmax-imgmin)
                imgmax = 1.0
                imgmin = 0.0
            self.brightness_min = 0.0
            self.brightness_max = 1.0
            im = axes.imshow(image, cmap=mtplot.cm.get_cmap(self.colortable), vmin=(imgmin+imgmax*self.brightness_min),vmax=imgmax*self.brightness_max)
        
        if self.show_colorbar == 1:
            cbar = axes.figure.colorbar(im, orientation='vertical', cax=axcb) 
        
        if self.show_scale_bar == 1:
            startx = int(self.stk.n_rows*0.05)
            starty = self.stk.n_cols-int(self.stk.n_cols*0.05)-self.stk.scale_bar_pixels_y
            um_string = ' $\mathrm{\mu m}$'
            microns = '$'+self.stk.scale_bar_string+' $'+um_string
            axes.text(self.stk.scale_bar_pixels_x+startx+1,starty+1, microns, horizontalalignment='left', verticalalignment='center', color='black', fontsize=14)
            #Matplotlib has flipped scales so I'm using rows instead of cols!
            p = mtplot.patches.Rectangle((startx,starty), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y, color = 'black', fill = True)
            axes.add_patch(p)
        
        self.ODReconPanel.draw()
        
        #self.tc_imageeng.SetValue("Image at energy: {0:5.2f} eV".format(float(self.stk.ev[self.iev])))

#----------------------------------------------------------------------
    def loadODErrorImage(self):

        self.show_colorbar = 1
        self.show_scale_bar = 0
        
        image = self.nnma.ODError[self.iev, :, :] 

        fig = self.ODErrorPanel.get_figure()
        fig.clf()
        
        if self.show_colorbar == 0:
            fig.add_axes([0.15, 0.15, 0.8, 0.8])
            axes = fig.gca()
        else:
            axes = fig.gca()
            divider = make_axes_locatable(axes)
            axcb = divider.new_horizontal(size="3%", pad=0.03)  
            fig.add_axes(axcb)
            axes.set_position([0.15, 0.15, 0.8, 0.8])
        
        fig.patch.set_alpha(1.0) 	# set figure transparency
        
        self.defaultdisplay = 0
        self.colortable = "gray"
        if self.defaultdisplay == 1.0:
            im = axes.imshow(image, cmap=mtplot.cm.get_cmap(self.colortable)) 
        else:
            imgmax = npy.amax(image)
            imgmin = npy.amin(image)
            if (imgmin < 0.0):
                image = (image-imgmin)/(imgmax-imgmin)
                imgmax = 1.0
                imgmin = 0.0
            self.brightness_min = 0.0
            self.brightness_max = 1.0
            im = axes.imshow(image, cmap=mtplot.cm.get_cmap(self.colortable), vmin=(imgmin+imgmax*self.brightness_min),vmax=imgmax*self.brightness_max)
        
        if self.show_colorbar == 1:
            cbar = axes.figure.colorbar(im, orientation='vertical', cax=axcb) 
        
        if self.show_scale_bar == 1:
            startx = int(self.stk.n_rows*0.05)
            starty = self.stk.n_cols-int(self.stk.n_cols*0.05)-self.stk.scale_bar_pixels_y
            um_string = ' $\mathrm{\mu m}$'
            microns = '$'+self.stk.scale_bar_string+' $'+um_string
            axes.text(self.stk.scale_bar_pixels_x+startx+1,starty+1, microns, horizontalalignment='left', verticalalignment='center', color='black', fontsize=14)
            #Matplotlib has flipped scales so I'm using rows instead of cols!
            p = mtplot.patches.Rectangle((startx,starty), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y, color = 'black', fill = True)
            axes.add_patch(p)
        
        self.ODErrorPanel.draw()
        
        #self.tc_imageeng.SetValue("Image at energy: {0:5.2f} eV".format(float(self.stk.ev[self.iev])))


# ------------------------------------------------------------------------------------------------
# Subclass of PageNNMA to display absorption spectra results
# ------------------------------------------------------------------------------------------------
class PageNNMASpectra(PageNNMA):

    def __init__(self, parent, common, data_struct, stack, nnma):

        #wx.Panel.__init__(self, parent)
        PageNNMA.__init__(self, parent, common, data_struct, stack, nnma)
            
        self.data_struct = data_struct
        self.stk = stack
        self.com = common
        self.nnma = nnma
        self.kNNMA = self.nnma.kNNMA	# number of chemical components to look for in NNMA
        self.fontsize = self.com.fontsize        
        self.iev = 0 

        # Panel 2: Display results for reconstructed absorption spectra
        panel2 = wx.Panel(self, -1)
        sizer21 = wx.StaticBoxSizer(wx.StaticBox(panel2, -1, 'Absorption spectra'), wx.VERTICAL)

        vbox21 = wx.BoxSizer(wx.VERTICAL)
        scaleFactorW = 0.9
        scaleFactorH = 0.8
        self.NNMASpectraImagePanel = wxmpl.PlotPanel(panel2, -1, size=(PlotW*scaleFactorW, PlotH*scaleFactorH), cursor=False, crosshairs=False, location=False, zoom=False)
        wxmpl.EVT_POINT(panel2, self.NNMASpectraImagePanel.GetId(), self.onPointNNMASpectraImage)
        vbox21.Add(self.NNMASpectraImagePanel, 0, wx.TOP)
        sizer21.Add(vbox21)

        vbox22 = wx.BoxSizer(wx.VERTICAL)
        scaleFactorW = 0.9
        scaleFactorH = 0.8
        self.NNMASpectraPanel = wxmpl.PlotPanel(panel2, -1, size=(PlotW*scaleFactorW, PlotH*scaleFactorH), cursor=False, crosshairs=False, location=False, zoom=False)
        vbox22.Add(self.NNMASpectraPanel, 0, wx.TOP)
        vbox22.Add((0, 10))
        sizer21.Add(vbox22)

        # Button to display NNMA spectra
        self.button_showNNMASpectra = wx.Button(panel2, -1, 'Show spectra', (10, 10))
        self.button_showNNMASpectra.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.onShowNNMASpectra, id=self.button_showNNMASpectra.GetId())
        #self.button_calcNNMA.Disable()
        sizer21.Add(self.button_showNNMASpectra, 0, wx.EXPAND)
        
        # Button to clear displayed NNMA spectra
        self.button_clearNNMASpectra = wx.Button(panel2, -1, 'Clear spectra', (10, 10))
        self.button_clearNNMASpectra.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.onClearNNMASpectra, id=self.button_clearNNMASpectra.GetId())
        #self.button_calcNNMA.Disable()
        sizer21.Add(self.button_clearNNMASpectra, 0, wx.EXPAND)
        
        panel2.SetSizer(sizer21)
        
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        gridtop = wx.FlexGridSizer(1, 0, vgap=10, hgap=20)
        gridtop.Add(panel2, 0)
        vboxtop.Add((0, 10))
        vboxtop.Add(gridtop, 0, wx.LEFT, 20)
        self.SetSizer(gridtop) 


#----------------------------------------------------------------------
    def onShowNNMASpectra(self, event):

        wx.BeginBusyCursor()
        
        try:
            if PageNNMA.calcnnma == True:	# if NNMA already calculated, then just load the spectra
                if verbose: print("PageNNMA.calcnnma == True")
                self.loadNNMASpectraImage()
                wx.EndBusyCursor()
        except:
            try:
                PageNNMA.calcnnma = False
                self.calcNNMA()
                if verbose: print("Successfully calculated self.calcNNMA()")
                PageNNMA.calcnnma = True
                self.loadNNMASpectraImage()
                wx.EndBusyCursor()
            except:
                wx.EndBusyCursor()
                wx.MessageBox("NNMA not calculated.")
          
        #wx.GetApp().TopWindow.refresh_widgets() 	# <-- Is this necessary??

#----------------------------------------------------------------------
    def calcNNMA(self):

        self.kNNMA = self.nnma.kNNMA
        self.algoNNMA = self.nnma.algoNNMA
        try:
            self.nnma.maxIters = int(self.itersCtrl.GetValue())
        except:
            wx.MessageBox("Please enter valid number of iterations.")
        self.nnma.calcNNMA(self.algoNNMA, kComponents=self.kNNMA)

#----------------------------------------------------------------------
    def onPointNNMASpectraImage(self, evt):
        if verbose: print("Inside onPointNNMASpectraImage()")
        x = evt.xdata
        y = evt.ydata
        if verbose: print("x = ", x)
        if verbose: print("y = ", y)
        energyIndex = npy.floor(x)
        muIndex = npy.floor(y)
        if verbose: print("energyIndex = ", energyIndex)
        if verbose: print("muIndex = ", muIndex)
        if self.calcnnma == True:
            self.loadNNMASpectra(muIndex)

#----------------------------------------------------------------------
    def loadNNMASpectraImage(self):
          
        self.defaultdisplay = 1.0
        self.show_colorbar = 1
        self.show_scale_bar = 0
        self.colortable = 'gist_rainbow'
 
        if self.defaultdisplay == 1.0:
            image = self.nnma.mu.T		# pointer to data (not a copy); transposed so (normalized) energy is along horizontal axis
        else:   
            # Adjustment to the data display setting has been made so make a copy
            image = self.nnma.mu.T.copy()

        fig = self.NNMASpectraImagePanel.get_figure()
        fig.clf()

        if self.show_colorbar == 0:
            fig.add_axes((0.15, 0.15, 0.75, 0.75))
            axes = fig.gca()
        else:
            axes = fig.gca()
            divider = make_axes_locatable(axes)
            axcb = divider.new_horizontal(size="3%", pad=0.03)  
            fig.add_axes(axcb)
            axes.set_position([0.15, 0.15, 0.75, 0.75])

        axes.set_xlabel("Energy index")
        axes.set_ylabel("NNMA component number") 
        fig.patch.set_alpha(1.0)

        maxEnergyIndex = self.stk.n_ev
        maxMuIndex = self.nnma.kNNMA 
        if self.defaultdisplay == 1.0:
            im = axes.imshow(image, cmap=mtplot.cm.get_cmap(self.colortable), interpolation='nearest', extent=[0, maxEnergyIndex, maxMuIndex, 0], aspect='auto')
        else:
            imgmax = npy.amax(image)
            imgmin = npy.amin(image)
            image = (image-imgmin)/(imgmax-imgmin)
            imgmax = 1.0
            imgmin = 0.0
            self.brightness_min = 0.0
            self.brightness_max = 1.0
            im = axes.imshow(image, cmap=mtplot.cm.get_cmap(self.colortable), 
                	 vmin=(imgmin+imgmax*self.brightness_min),vmax=imgmax*self.brightness_max, interpolation='nearest', extent=[0, maxEnergyIndex, maxMuIndex, 0], aspect='auto')
              
        if self.show_colorbar == 1:
            cbar = axes.figure.colorbar(im, orientation='vertical', cax=axcb) 

        self.NNMASpectraImagePanel.draw()

        if verbose: print("Done loadNNMASpectraImage()")
  
#----------------------------------------------------------------------         
    def loadNNMASpectra(self, muIndex):

        fig = self.NNMASpectraPanel.get_figure()
        axes = fig.add_axes([0.15, 0.15, 0.75, 0.75])
        axes.plot(self.stk.ev, self.nnma.mu[:, muIndex])
        axes.set_xlabel("Energy (eV)")
        axes.set_ylabel("Absorption")
        self.NNMASpectraPanel.draw()

#----------------------------------------------------------------------         
    def onClearNNMASpectra(self, evt):

        fig = self.NNMASpectraPanel.get_figure()
        axes = fig.gca()
        axes.cla()
        self.NNMASpectraPanel.draw()

       
 

# ------------------------------------------------------------------------------------------------
# Subclass of PageNNMA to display thickness map results
# ------------------------------------------------------------------------------------------------
class PageNNMAThickness(PageNNMA):

    def __init__(self, parent, common, data_struct, stack, nnma):

        #wx.Panel.__init__(self, parent)
        PageNNMA.__init__(self, parent, common, data_struct, stack, nnma)
        


""" ------------------------------------------------------------------------------------------------"""
class PageSpectral(wx.Panel):
    def __init__(self, parent, common, data_struct, stack, anlz):
        wx.Panel.__init__(self, parent)
        
        self.data_struct = data_struct
        self.stk = stack
        self.com = common
        self.anlz = anlz
        
        self.i_tspec = 1
        self.showraw = True
        self.show_scale_bar = 0
        
        self.SetBackgroundColour("White")
           
        self.fontsize = self.com.fontsize        
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        hboxT = wx.BoxSizer(wx.HORIZONTAL)
        hboxB = wx.BoxSizer(wx.HORIZONTAL)
    
        #panel 1        
        panel1 = wx.Panel(self, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        panel1.SetBackgroundColour("White")      
  
        self.tc_spmap = wx.TextCtrl(panel1, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_spmap.SetFont(self.com.font)
        self.tc_spmap.SetValue("Spectrum composition map")

        i1panel = wx.Panel(panel1, -1, style = wx.SUNKEN_BORDER)
        self.MapPanel = wxmpl.PlotPanel(i1panel, -1, size =(PlotH, PlotH), cursor=False, crosshairs=False, location=False, zoom=False)                            
  
        vbox1.Add((0,10))
        vbox1.Add(self.tc_spmap,1, wx.LEFT | wx.EXPAND, 20)        
        vbox1.Add(i1panel, 0,  wx.LEFT, 20)

        panel1.SetSizer(vbox1)
     
     
        #panel 2
        panel2 = wx.Panel(self, -1)
        vbox2 = wx.BoxSizer(wx.VERTICAL)
        panel2.SetBackgroundColour("White")
        
        self.tc_tspec = wx.TextCtrl(panel2, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_tspec.SetFont(self.com.font)
        self.tc_tspec.SetValue("Target Spectrum: ")
        hbox11 = wx.BoxSizer(wx.HORIZONTAL)         
        
        i2panel = wx.Panel(panel2, -1, style = wx.SUNKEN_BORDER)
        self.TSpectrumPanel = wxmpl.PlotPanel(i2panel, -1, size=(PlotW, PlotH), cursor=False, crosshairs=False, location=False, zoom=False)

        self.slider_tspec = wx.Slider(panel2, -1, 1, 1, 5, style=wx.SL_LEFT|wx.SL_VERTICAL )        
        self.slider_tspec.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnTSScroll, self.slider_tspec)
        
        vbox21 = wx.BoxSizer(wx.VERTICAL)               
        self.tspecspin = wx.SpinButton(panel2, -1, size = ((8,-1)), style=wx.SP_ARROW_KEYS)
        self.Bind(wx.EVT_SPIN_UP, self.OnTspecSpinUp, self.tspecspin)
        self.Bind(wx.EVT_SPIN_DOWN, self.OnTspecSpinDown, self.tspecspin)
        
        vbox21.Add((0,3))
        vbox21.Add(self.slider_tspec, 1,  wx.EXPAND) 
        vbox21.Add(self.tspecspin, 0,  wx.EXPAND)      

        hbox11.Add(i2panel, 0)
        hbox11.Add(vbox21, 0,  wx.EXPAND)
          
        vbox2.Add((0,10))
        vbox2.Add(self.tc_tspec, 1, wx.LEFT | wx.EXPAND, 20)       
        vbox2.Add(hbox11, 0, wx.LEFT , 20)
        
        panel2.SetSizer(vbox2)
        
        
        #panel 3
        panel3 = wx.Panel(self, -1)
        sb = wx.StaticBox(panel3, -1, 'Target Spectrum')
        sb.SetBackgroundColour("white")
        sizer1 = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)
        vbox31 = wx.BoxSizer(wx.VERTICAL)
        vbox31.Add((0,10)) 
        
        self.button_loadtspec = wx.Button(panel3, -1, 'Load Spectrum')
        self.button_loadtspec.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnTSpecFromFile, id=self.button_loadtspec.GetId())
        self.button_loadtspec.Disable()
        vbox31.Add(self.button_loadtspec, 0, wx.EXPAND)
        self.button_addflat = wx.Button(panel3, -1, 'Add Flat Spectrum')
        self.button_addflat.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnFlatTSpec, id=self.button_addflat.GetId())
        self.button_addflat.Disable()
        vbox31.Add(self.button_addflat, 0, wx.EXPAND)
        self.button_addclspec = wx.Button(panel3, -1, 'Add Cluster Spectra')
        self.button_addclspec.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnAddClusterSpectra, id=self.button_addclspec.GetId())   
        self.button_addclspec.Disable()     
        vbox31.Add(self.button_addclspec, 0, wx.EXPAND)
        
        self.button_showrgb = wx.Button(panel3, -1, 'Composite RGB image...')
        self.button_showrgb.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnCompositeRGB, id=self.button_showrgb.GetId())   
        self.button_showrgb.Disable()     
        vbox31.Add(self.button_showrgb, 0, wx.EXPAND)        

        self.button_save = wx.Button(panel3, -1, 'Save Images...', (10,10))
        self.button_save.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnSave, id=self.button_save.GetId())
        self.button_save.Disable()          
        vbox31.Add(self.button_save, 0, wx.EXPAND)
        sizer1.Add(vbox31,1, wx.LEFT|wx.RIGHT|wx.EXPAND,2)
        panel3.SetSizer(sizer1)
        

        
        #panel 4
        panel4 = wx.Panel(self, -1)
        vbox4 = wx.BoxSizer(wx.VERTICAL)
        sb = wx.StaticBox(panel4, -1, 'Display')
        sb.SetBackgroundColour("white")
        sizer4 = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)

        sb = wx.StaticBox(panel4, -1, 'Spectrum')
        sb.SetBackgroundColour("white")
        sizer41 = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)
        self.textctrl_sp = wx.TextCtrl(panel4, -1, style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.textctrl_sp.SetFont(self.com.font)
        sizer41.Add(self.textctrl_sp, 1, wx.EXPAND|wx.TOP|wx.LEFT, 5)
      
        hbox40 = wx.BoxSizer(wx.HORIZONTAL)    
        hbox40.Add(sizer41, 1, wx.EXPAND)
        vbox4.Add(hbox40, 0, wx.EXPAND)
        
        self.textctrl_sp.AppendText('Common Name: \n')
        self.textctrl_sp.AppendText('RMS Error: ')
        
        hbox41 = wx.BoxSizer(wx.HORIZONTAL)
        
        sb = wx.StaticBox(panel4, -1, 'Composition Map')
        sb.SetBackgroundColour("white")
        sizer42 = wx.StaticBoxSizer(sb,  orient=wx.VERTICAL)
        self.rb_raw = wx.RadioButton(panel4, -1, 'Raw', style=wx.RB_GROUP)
        self.rb_fit = wx.RadioButton(panel4, -1, 'Fitted')
        self.rb_raw.SetFont(self.com.font)
        self.rb_fit.SetFont(self.com.font)
        self.Bind(wx.EVT_RADIOBUTTON, self.OnRBRawFit, id=self.rb_raw.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.OnRBRawFit, id=self.rb_fit.GetId())
        self.rb_raw.SetValue(True)

        self.add_scale_cb = wx.CheckBox(panel4, -1, '  Scale')
        self.add_scale_cb.SetFont(self.com.font)
        self.Bind(wx.EVT_CHECKBOX, self.OnShowScale, self.add_scale_cb)
        
        sizer42.Add((0,3))
        sizer42.Add(self.rb_raw)
        sizer42.Add((0,5))
        sizer42.Add(self.rb_fit)
        sizer42.Add((0,10))
        sizer42.Add(self.add_scale_cb)
        
        hbox41.Add(sizer42, 1, wx.EXPAND)
                
        sb = wx.StaticBox(panel4, -1, 'Fit Weights')
        sb.SetBackgroundColour("white")
        sizer43 = wx.StaticBoxSizer(sb,  orient=wx.VERTICAL)
        hbox42 = wx.BoxSizer(wx.HORIZONTAL)
        hbox42.Add((3,0))
        
        self.tc_spfitlist = wx.TextCtrl(panel4, -1, style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_spfitlist.SetFont(self.com.font)
        
        hbox42.Add(self.tc_spfitlist,1,wx.EXPAND)
        
        sizer43.Add(hbox42,1,wx.EXPAND)       

        
        vbox43 = wx.BoxSizer(wx.VERTICAL)
        self.button_removespec = wx.Button(panel4, -1, 'Remove Spectrum')
        self.button_removespec.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnRemoveSpectrum, id=self.button_removespec.GetId())   
        self.button_removespec.Disable()     
        vbox43.Add(self.button_removespec, 0, wx.EXPAND)
        self.button_movespup = wx.Button(panel4, -1, 'Move Spectrum Up')
        self.button_movespup.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnMoveSpectrumUp, id=self.button_movespup.GetId())   
        self.button_movespup.Disable()     
        vbox43.Add(self.button_movespup, 0, wx.EXPAND)
        self.button_movespdown = wx.Button(panel4, -1, 'Move Spectrum Down')
        self.button_movespdown.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnMoveSpectrumDown, id=self.button_movespdown.GetId())   
        self.button_movespdown.Disable()     
        vbox43.Add(self.button_movespdown, 0, wx.EXPAND)
                       

               
            
        hbox41.Add((10,0))
        hbox41.Add(sizer43, 1, wx.EXPAND)
        hbox41.Add(vbox43, 1, wx.EXPAND|wx.ALL, 5)
        vbox4.Add(hbox41, 1, wx.EXPAND)

        sizer4.Add(vbox4,1, wx.EXPAND)
        
        panel4.SetSizer(sizer4)
        
        
        #panel 5
        panel5 = wx.Panel(self, -1)
        sb = wx.StaticBox(panel5, -1, 'Target Spectra')
        sb.SetBackgroundColour("white")
        sizer5 = wx.StaticBoxSizer(sb, orient= wx.VERTICAL)
        
        hbox51 = wx.BoxSizer(wx.HORIZONTAL)
        hbox51.Add((0,2))

        self.tc_speclist =  wx.ListCtrl(panel5, -1, 
                                        style=wx.LC_REPORT|wx.LC_NO_HEADER|wx.NO_BORDER|wx.LC_EDIT_LABELS|wx.LC_SINGLE_SEL)
        self.tc_speclist.InsertColumn(0, 'Spectra')
        self.Bind(wx.EVT_LIST_ITEM_FOCUSED , self.OnSpectraListClick, self.tc_speclist)
        self.Bind(wx.EVT_LIST_END_LABEL_EDIT, self.OnEditSpectraListClick, self.tc_speclist)
        self.tc_speclist.SetBackgroundColour('white')
        self.tc_speclist.SetFont(self.com.font)
        hbox51.Add(self.tc_speclist, 1, wx.EXPAND)
        sizer5.Add(hbox51,1, wx.ALL|wx.EXPAND,2)        
        panel5.SetSizer(sizer5)
        
#        self.tc_speclist = wx.TextCtrl(panel5, -1, style=wx.TE_MULTILINE|wx.TE_RICH|wx.BORDER_NONE)
#        self.tc_speclist.SetFont(self.com.font)
#        hbox51.Add(self.tc_speclist, 1, wx.EXPAND)
#        sizer5.Add(hbox51,1, wx.ALL|wx.EXPAND,2)        
#        panel5.SetSizer(sizer5)

        
        hboxB.Add(panel2, 0, wx.BOTTOM | wx.TOP, 9)
        hboxB.Add(panel1, 0, wx.BOTTOM | wx.TOP, 9)
        hboxT.Add((10,0)) 
               
        hboxT.Add(panel3, 0, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 9)
        hboxT.Add(panel4, 2.5, wx.LEFT | wx.RIGHT |wx.TOP | wx.EXPAND, 9)
        hboxT.Add(panel5, 1, wx.LEFT | wx.RIGHT |wx.TOP | wx.EXPAND, 9)

        vbox.Add(hboxT, 0, wx.ALL, 5)
        
        vbox.Add((0, 5))
        
        vbox.Add(hboxB, 0, wx.LEFT | wx.RIGHT, 5)
  
        self.SetSizer(vbox) 
        
        
        
#----------------------------------------------------------------------
    def OnTSpecFromFile(self, event):
        

        #try: 
        if True:
            wildcard = "Spectrum files (*.csv)|*.csv"
            dialog = wx.FileDialog(None, "Choose Spectrum file",
                                   defaultDir = defaultDir,
                                   style=wx.OPEN)
            dialog.SetWildcard(wildcard)
            if dialog.ShowModal() == wx.ID_OK:
                filepath = dialog.GetPath()
                self.filename = dialog.GetFilename()
                                                        
            wx.BeginBusyCursor()    
                                            
            self.anlz.read_target_spectrum(filename=filepath)
            self.com.spec_anl_calculated = 1
            
            self.i_tspec = self.anlz.n_target_spectra      
            self.slider_tspec.SetMax(self.anlz.n_target_spectra)
            self.slider_tspec.SetValue(self.i_tspec)
            
            self.loadTSpectrum()
            self.loadTargetMap()    
            self.ShowSpectraList()
                    
            wx.EndBusyCursor()
            
#         except:
#             wx.EndBusyCursor()  
#             wx.MessageBox("Spectrum file not loaded.")
                                   
        dialog.Destroy()
                                 
        wx.GetApp().TopWindow.refresh_widgets()
        

#----------------------------------------------------------------------
    def OnFlatTSpec(self, event):

        try: 
            wx.BeginBusyCursor() 
            self.anlz.read_target_spectrum(flat=True)
            self.com.spec_anl_calculated = 1
            
            self.i_tspec = self.anlz.n_target_spectra      
            self.slider_tspec.SetMax(self.anlz.n_target_spectra)
            self.slider_tspec.SetValue(self.i_tspec)
            
            self.loadTSpectrum()
            self.loadTargetMap()
            self.ShowSpectraList()
        
            wx.EndBusyCursor()
            
        except:
            wx.EndBusyCursor()  
            wx.MessageBox("Flat spectrum not loaded.")
            
                                                      
        wx.GetApp().TopWindow.refresh_widgets()
        
#----------------------------------------------------------------------
    def OnAddClusterSpectra(self, event):

        #try:
        if True: 
            wx.BeginBusyCursor() 
            self.anlz.add_cluster_target_spectra()
            self.com.spec_anl_calculated = 1
            
            self.i_tspec = self.anlz.n_target_spectra      
            self.slider_tspec.SetMax(self.anlz.n_target_spectra)
            self.slider_tspec.SetValue(self.i_tspec)
            
            self.ShowSpectraList() 
            self.loadTSpectrum()
            self.loadTargetMap()  
             
        
            wx.EndBusyCursor()
            
#        except:
#            wx.EndBusyCursor()  
#            wx.MessageBox("Cluster spectra not loaded.")
 
                                                        
        wx.GetApp().TopWindow.refresh_widgets()
        
#----------------------------------------------------------------------
    def OnCompositeRGB(self, event):

        ShowCompositeRBGmap(self.com, self.anlz).Show()
                
#----------------------------------------------------------------------
    def OnSave(self, event):
        
        SaveWinP4().Show()
        
        
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
   
            wx.MessageBox('Could not save file: %s' % err, 'Error', 
                          parent=self, style=wx.OK|wx.ICON_ERROR) 
            
#----------------------------------------------------------------------
    def SaveSpectra(self, png_pdf=1, savecsv = False):
        
        
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas   
        mtplot.rcParams['pdf.fonttype'] = 42
            
        colors=['#FF0000','#000000','#FFFFFF']
        spanclrmap=mtplot.colors.LinearSegmentedColormap.from_list('spancm',colors)
        
        if png_pdf == 1:   
            ext = 'png'
        else:
            ext = 'pdf'
        suffix = "." + ext
        
        
        for i in range (self.anlz.n_target_spectra):
            #Save spectra images
            tspectrum = self.anlz.target_spectra[i, :]
                        
        
            fig = mtplot.figure.Figure(figsize =(PlotW, PlotH))
            canvas = FigureCanvas(fig)
            fig.clf()
            fig.add_axes((0.15,0.15,0.75,0.75))
            axes = fig.gca()
        
            mtplot.rcParams['font.size'] = self.fontsize

            line1 = axes.plot(self.stk.ev,tspectrum, color='black', label = 'Raw data')
            
            if self.com.pca_calculated == 1: 
                tspectrumfit = self.anlz.target_pcafit_spectra[i, :]
                diff = npy.abs(tspectrum-tspectrumfit)
                line2 = axes.plot(self.stk.ev,tspectrumfit, color='green', label = 'Fit')
                line3 = axes.plot(self.stk.ev,diff, color='grey', label = 'Abs(Raw-Fit)')
            
            fontP = mtplot.font_manager.FontProperties()
            fontP.set_size('small')
       
            axes.legend(loc=4, prop = fontP)
                        
            axes.set_xlabel('Photon Energy [eV]')
            axes.set_ylabel('Optical Density')

            fileName_spec = self.SaveFileName+"_Tspectrum_" +str(i+1)+"."+ext
            fig.savefig(fileName_spec)    
            
            if savecsv:
                fileName_spec = self.SaveFileName+"_Tspectrum_" +str(i+1)+".csv"
                cname = 'Tspectrum_' +str(i+1)
                self.stk.write_csv(fileName_spec, self.stk.ev, tspectrum, cname = cname)
#----------------------------------------------------------------------
    def SaveMaps(self, png_pdf=1):            
            
            
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas  
        mtplot.rcParams['pdf.fonttype'] = 42 
            
        colors=['#FF0000','#000000','#FFFFFF']
        spanclrmap=mtplot.colors.LinearSegmentedColormap.from_list('spancm',colors)
            
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
  
            fig = mtplot.figure.Figure(figsize =(PlotH, PlotH))
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
                um_string = ' $\mathrm{\mu m}$'
                microns = '$'+self.stk.scale_bar_string+' $'+um_string
                axes.text(self.stk.scale_bar_pixels_x+10,self.stk.n_cols-9, microns, horizontalalignment='left', verticalalignment='center',
                              color = 'white', fontsize=14)
                #Matplotlib has flipped scales so I'm using rows instead of cols!
                p = mtplot.patches.Rectangle((5,self.stk.n_cols-10), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                                            color = 'white', fill = True)
                axes.add_patch(p)     
     
            im = axes.imshow(tsmapimage, cmap=spanclrmap, vmin = -bound, vmax = bound)
            cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)  
    
            axes.axis("off") 
                
                   
            fileName_img = self.SaveFileName+"_TSmap_" +str(i+1)+"."+ext               
            fig.savefig(fileName_img, bbox_inches='tight', pad_inches = 0.0)
            
#----------------------------------------------------------------------        
    def OnEditSpectraListClick(self, event):
        self.anlz.tspec_names[self.i_tspec-1] = event.GetText()
        self.loadTSpectrum()

#----------------------------------------------------------------------        
    def OnSpectraListClick(self, event):
        sel = event.m_itemIndex
        self.i_tspec = sel+1
        
        if self.com.spec_anl_calculated == 1:
            self.loadTSpectrum()
            self.loadTargetMap()
            self.slider_tspec.SetValue(self.i_tspec)
            
#----------------------------------------------------------------------        
    def OnTSScroll(self, event):
        self.tc_speclist.SetItemState(self.i_tspec-1, 0, wx.LIST_STATE_SELECTED) 
        
        sel = event.GetInt()
        self.i_tspec = sel
        if self.com.spec_anl_calculated == 1:
            self.loadTSpectrum()
            self.loadTargetMap()
            
#----------------------------------------------------------------------            
    def OnTspecSpinUp(self, event):
        if (self.com.spec_anl_calculated == 1) and (self.i_tspec > 1):
            self.i_tspec = self.i_tspec - 1
            self.slider_tspec.SetValue(self.i_tspec)

            self.loadTSpectrum()
            self.loadTargetMap()
            
#----------------------------------------------------------------------            
    def OnTspecSpinDown(self, event):
        if (self.com.spec_anl_calculated == 1) and (self.i_tspec < self.anlz.n_target_spectra):
            self.i_tspec = self.i_tspec + 1
            self.slider_tspec.SetValue(self.i_tspec) 
            
            self.loadTSpectrum()
            self.loadTargetMap()


#----------------------------------------------------------------------          
    def OnRBRawFit(self, evt):
        state = self.rb_raw.GetValue()
           
        if state:
            self.showraw = True
        else:        
            self.showraw = False
            
        if self.com.spec_anl_calculated == 1:
            self.loadTSpectrum()
            self.loadTargetMap()
            
#----------------------------------------------------------------------           
    def OnShowScale(self, event):
        if self.add_scale_cb.GetValue():
            self.show_scale_bar = 1
        else: self.show_scale_bar = 0
        
        if self.com.spec_anl_calculated == 1:
            self.loadTSpectrum()
            self.loadTargetMap()
            
#----------------------------------------------------------------------           
    def OnRemoveSpectrum(self, event):
        wx.BeginBusyCursor() 
        self.anlz.remove_spectrum(self.i_tspec-1)
        self.com.spec_anl_calculated = 1
            
        self.i_tspec = self.i_tspec-1
        if self.i_tspec<0:
            self.i_tspec=0
        self.slider_tspec.SetMax(self.anlz.n_target_spectra)
        self.slider_tspec.SetValue(self.i_tspec)
            
        if self.anlz.tspectrum_loaded == 1:
            self.loadTSpectrum()
            self.loadTargetMap()  
            self.ShowSpectraList()  
        else:
            self.com.spec_anl_calculated = 0
            self.ClearWidgets()
        
        wx.EndBusyCursor()
        
#----------------------------------------------------------------------           
    def OnMoveSpectrumDown(self, event):
        
        if self.i_tspec < self.anlz.n_target_spectra:
            self.anlz.move_spectrum(self.i_tspec-1, self.i_tspec)
            
            self.i_tspec += 1
            self.slider_tspec.SetValue(self.i_tspec)
            self.loadTSpectrum()
            self.loadTargetMap()  
            self.ShowSpectraList()
        
        
#----------------------------------------------------------------------           
    def OnMoveSpectrumUp(self, event):        
        
        if self.i_tspec > 1:
            self.anlz.move_spectrum(self.i_tspec-1, self.i_tspec-2)      
            
            self.i_tspec -= 1
            self.slider_tspec.SetValue(self.i_tspec)
            self.loadTSpectrum()
            self.loadTargetMap() 
            self.ShowSpectraList()
        
#----------------------------------------------------------------------           
    def ClearWidgets(self):
        
        fig = self.MapPanel.get_figure()
        fig.clf()
        self.MapPanel.draw()
        fig = self.TSpectrumPanel.get_figure()
        fig.clf()
        self.TSpectrumPanel.draw()
        
        self.tc_tspec.SetValue("Target Spectrum: ")
        self.tc_speclist.DeleteAllItems()
        self.tc_spfitlist.Clear()
        
        self.com.spec_anl_calculated = 0
        self.i_tspec = 1
        self.showraw = True
        self.rb_raw.SetValue(True)
        
        self.slider_tspec.SetValue(self.i_tspec)
        
        self.textctrl_sp.Clear()
        
        self.textctrl_sp.AppendText('Common Name: \n')
        self.textctrl_sp.AppendText('RMS Error: ')
        
        wx.GetApp().TopWindow.refresh_widgets()
            
#----------------------------------------------------------------------           
    def ShowFitWeights(self):    
        
        self.tc_spfitlist.Clear()     
        
        norm_factor = 100./npy.sum(npy.absolute(self.anlz.target_pcafit_coeffs[self.i_tspec-1, :]))
        
        for i in range(self.anlz.numsigpca):
            self.tc_spfitlist.AppendText('{0}: {1:5.2f} %'.format(i+1, norm_factor
                                                    *abs(self.anlz.target_pcafit_coeffs[self.i_tspec-1, i])))
            if i < self.anlz.numsigpca-1:
                self.tc_spfitlist.AppendText('\n')
        self.tc_spfitlist.SetInsertionPoint(0)

#----------------------------------------------------------------------           
    def ShowSpectraList(self):    
        
        self.tc_speclist.DeleteAllItems()   
        
        for i in range(self.anlz.n_target_spectra):
            self.tc_speclist.InsertStringItem(i, self.anlz.tspec_names[i])
            
        
#----------------------------------------------------------------------      
    def loadTargetMap(self):

        if self.showraw == True:
            tsmapimage = self.anlz.target_svd_maps[:,:,self.i_tspec-1]
        else:
            tsmapimage = self.anlz.target_pcafit_maps[:,:,self.i_tspec-1] 
        
        
        colors=['#FF0000','#000000','#FFFFFF']
        
        spanclrmap=mtplot.colors.LinearSegmentedColormap.from_list('spancm',colors)
                
        fig = self.MapPanel.get_figure()
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
            um_string = ' $\mathrm{\mu m}$'
            microns = '$'+self.stk.scale_bar_string+' $'+um_string
            axes.text(self.stk.scale_bar_pixels_x+10,self.stk.n_cols-9, microns, horizontalalignment='left', verticalalignment='center',
                      color = 'white', fontsize=14)
            #Matplotlib has flipped scales so I'm using rows instead of cols!
            p = mtplot.patches.Rectangle((5,self.stk.n_cols-10), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                                   color = 'white', fill = True)
            axes.add_patch(p)     
     
        im = axes.imshow(tsmapimage, cmap=spanclrmap, vmin = -bound, vmax = bound)
        cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)  
    
        axes.axis("off") 
        self.MapPanel.draw()

        
        
#----------------------------------------------------------------------     
    def loadTSpectrum(self):

        tspectrum = self.anlz.target_spectra[self.i_tspec-1, :]
                 
        
        fig = self.TSpectrumPanel.get_figure()
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        
        mtplot.rcParams['font.size'] = self.fontsize

        line1 = axes.plot(self.stk.ev,tspectrum, color='black', label = 'Raw data')

        if self.com.pca_calculated == 1:       
            tspectrumfit = self.anlz.target_pcafit_spectra[self.i_tspec-1, :]  
            diff = npy.abs(tspectrum-tspectrumfit)      
            line2 = axes.plot(self.stk.ev,tspectrumfit, color='green', label = 'Fit')
        
            line3 = axes.plot(self.stk.ev,diff, color='grey', label = 'Abs(Raw-Fit)')
        

        fontP = mtplot.font_manager.FontProperties()
        fontP.set_size('small')
       
        axes.legend(loc=4, prop = fontP)
                        
        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Optical Density')
        
        self.TSpectrumPanel.draw()
        
        self.tc_tspec.SetValue("Target Spectrum: " + 
                               self.anlz.tspec_names[self.i_tspec-1])
        
        self.textctrl_sp.Clear()
        
        self.textctrl_sp.AppendText('Common Name: '+ 
                                    self.anlz.tspec_names[self.i_tspec-1]+'\n')
        if self.com.pca_calculated == 1:       
            self.textctrl_sp.AppendText('RMS Error: '+ str('{0:7.5f}').format(self.anlz.target_rms[self.i_tspec-1]))
            
            self.ShowFitWeights()
        
        
#---------------------------------------------------------------------- 
class ShowCompositeRBGmap(wx.Frame):

    title = "Composite RBG Map"

    def __init__(self, common, anlz):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=self.title, size=(630, 500))
               
        ico = logos.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White") 
        
        self.com = common 
        self.anlz = anlz
        
        self.com = wx.GetApp().TopWindow.common         
        self.fontsize = self.com.fontsize
        
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
        
    
        vboxtop = wx.BoxSizer(wx.HORIZONTAL)
        
        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)   
        
        vbox1 = wx.BoxSizer(wx.VERTICAL)    
       
        sizer1 = wx.StaticBoxSizer(wx.StaticBox(panel, -1, 'Red spectrum'), orient=wx.VERTICAL)

        fgs1 = wx.FlexGridSizer(3, 2, 2, 5)
        r = wx.StaticText(panel, label="Red")
        r.SetFont(self.com.font)
        rl = wx.StaticText(panel, label="Limits")
        rl.SetFont(self.com.font)
        rw = wx.StaticText(panel, label="Weight")
        rw.SetFont(self.com.font)        
        
        
        self.combor = wx.ComboBox(panel, size=(150, -1), choices=self.anlz.tspec_names, style=wx.CB_READONLY)       
        self.Bind(wx.EVT_COMBOBOX, self.OnSelectR, self.combor)
        self.combor.SetToolTip(wx.ToolTip("select spectrum from dropdown-list"))
        self.combor.SetValue(self.anlz.tspec_names[self.r_spec])
        
        hbox12 = wx.BoxSizer(wx.HORIZONTAL)   
        
        bsize = 75
        
        self.tcrmin = wx.lib.intctrl.IntCtrl( panel, size=( bsize, -1 ), value = 0 , limited = True )
        self.tcrmax = wx.lib.intctrl.IntCtrl( panel, size=( bsize, -1 ), value = 100, limited = True )
        
        self.tcrmin.SetMin(0)
        self.tcrmin.SetMax(100)
        self.tcrmax.SetMin(0)
        self.tcrmax.SetMax(100)
        
                        
        hbox12.Add(self.tcrmin)
        hbox12.Add(self.tcrmax)
        
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnLimitMinR, self.tcrmin)
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnLimitMaxR, self.tcrmax)
        
        self.tcrweight = wx.lib.intctrl.IntCtrl( panel, value = 100, limited = True )
        
        self.tcrweight.SetMin(0)
        self.tcrweight.SetMax(100)
        
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnWeightR, self.tcrweight)
        
        fgs1.AddMany([(r, 0, wx.ALIGN_CENTER_VERTICAL), (self.combor, 0, wx.EXPAND), (rl, 0, wx.ALIGN_CENTER_VERTICAL), 
            (hbox12, 0, wx.EXPAND),(rw, 0, wx.ALIGN_CENTER_VERTICAL), (self.tcrweight, 0, wx.EXPAND)])
        
        sizer1.Add(fgs1, 0, wx.EXPAND|wx.ALL, 10)
        
        
        sizer2 = wx.StaticBoxSizer(wx.StaticBox(panel, -1, 'Green spectrum'), orient=wx.VERTICAL)
        
        fgs2 = wx.FlexGridSizer(3, 2, 2, 5)
        g = wx.StaticText(panel, label="Green")
        g.SetFont(self.com.font)
        gl = wx.StaticText(panel, label="Limits")
        gl.SetFont(self.com.font)
        gw = wx.StaticText(panel, label="Weight")
        gw.SetFont(self.com.font)        
        
        
        self.combog = wx.ComboBox(panel, size=(150, -1), choices=self.anlz.tspec_names, style=wx.CB_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.OnSelectG, self.combog)        
        self.combog.SetToolTip(wx.ToolTip("select spectrum from dropdown-list"))
        self.combog.SetValue(self.anlz.tspec_names[self.g_spec])
        
        hbox22 = wx.BoxSizer(wx.HORIZONTAL)   
        
        self.tcgmin = wx.lib.intctrl.IntCtrl( panel, size=( bsize, -1 ), value = 0, limited = True   )
        self.tcgmax = wx.lib.intctrl.IntCtrl( panel, size=( bsize, -1 ), value = 100, limited = True  )
        
        self.tcgmin.SetMin(0)
        self.tcgmin.SetMax(100)
        self.tcgmax.SetMin(0)
        self.tcgmax.SetMax(100)
        
        hbox22.Add(self.tcgmin)
        hbox22.Add(self.tcgmax)
        
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnLimitMinG, self.tcgmin)
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnLimitMaxG, self.tcgmax)
        
        self.tcgweight = wx.lib.intctrl.IntCtrl( panel, value = 100, limited = True )
        
        self.tcgweight.SetMin(0)
        self.tcgweight.SetMax(100)
        
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnWeightG, self.tcgweight)
        
        fgs2.AddMany([(g, 0, wx.ALIGN_CENTER_VERTICAL), (self.combog, 0, wx.EXPAND), (gl, 0, wx.ALIGN_CENTER_VERTICAL), 
            (hbox22, 0, wx.EXPAND),(gw, 0, wx.ALIGN_CENTER_VERTICAL), (self.tcgweight, 0, wx.EXPAND)])
        
        sizer2.Add(fgs2, 0, wx.EXPAND|wx.ALL, 10)
        
        
        sizer3 = wx.StaticBoxSizer(wx.StaticBox(panel, -1, 'Blue spectrum'), orient=wx.VERTICAL)
        
        fgs3 = wx.FlexGridSizer(3, 2, 2, 5)
        b = wx.StaticText(panel, label="Blue")
        b.SetFont(self.com.font)
        bl = wx.StaticText(panel, label="Limits")
        bl.SetFont(self.com.font)
        bw = wx.StaticText(panel, label="Weight")
        bw.SetFont(self.com.font)        
        
               
        self.combob = wx.ComboBox(panel, size=(150, -1), choices=self.anlz.tspec_names, style=wx.CB_READONLY)        
        self.Bind(wx.EVT_COMBOBOX, self.OnSelectB, self.combob)
        self.combob.SetToolTip(wx.ToolTip("select spectrum from dropdown-list"))
        self.combob.SetValue(self.anlz.tspec_names[self.b_spec])
        
        hbox32 = wx.BoxSizer(wx.HORIZONTAL)   
        
        self.tcbmin = wx.lib.intctrl.IntCtrl( panel, size=( bsize, -1 ), value = 0, limited = True   )
        self.tcbmax = wx.lib.intctrl.IntCtrl( panel, size=( bsize, -1 ), value = 100, limited = True  )
        
        self.tcbmin.SetMin(0)
        self.tcbmin.SetMax(100)
        self.tcbmax.SetMin(0)
        self.tcbmax.SetMax(100)
        
        hbox32.Add(self.tcbmin)
        hbox32.Add(self.tcbmax)
        
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnLimitMinB, self.tcbmin)
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnLimitMaxB, self.tcbmax)
        
        self.tcbweight = wx.lib.intctrl.IntCtrl( panel, value = 100, limited = True  )
        
        self.tcbweight.SetMin(0)
        self.tcbweight.SetMax(100)
        
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnWeightB, self.tcbweight)
        
        fgs3.AddMany([(b, 0, wx.ALIGN_CENTER_VERTICAL), (self.combob, 0, wx.EXPAND), (bl, 0, wx.ALIGN_CENTER_VERTICAL), 
            (hbox32, 0, wx.EXPAND),(bw, 0, wx.ALIGN_CENTER_VERTICAL), (self.tcbweight, 0, wx.EXPAND)])
        
        sizer3.Add(fgs3, 0, wx.EXPAND|wx.ALL, 10)
        
                
        vbox1.Add(sizer1, 0, wx.EXPAND)
        vbox1.Add(sizer2, 0, wx.EXPAND)
        vbox1.Add(sizer3, 0, wx.EXPAND)
        
        self.show_info_cb = wx.CheckBox(panel, -1, '  Show Info on the Image')
        self.show_info_cb.SetFont(self.com.font)
        self.Bind(wx.EVT_CHECKBOX, self.OnShowInfo, self.show_info_cb)
        vbox1.Add((0,10))
        vbox1.Add(self.show_info_cb, 0, wx.EXPAND)
        
        hbox1.Add(vbox1, 0,  wx.LEFT|wx.RIGHT ,20)

        i1panel = wx.Panel(panel, -1, style = wx.SUNKEN_BORDER)
        self.RGBImagePanel = wxmpl.PlotPanel(i1panel, -1, size =(PlotH,PlotH), cursor=False, crosshairs=False, location=False, zoom=False)
        hbox1.Add(i1panel, 0)
        
        vbox.Add(hbox1, 0, wx.EXPAND| wx.TOP, 10) 
        
             
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
               
        button_save = wx.Button(panel, -1, 'Save image')
        self.Bind(wx.EVT_BUTTON, self.OnSave, id=button_save.GetId())
        hbox2.Add(button_save, 1, wx.ALL,20)
        
        button_close = wx.Button(panel, -1, 'Dismiss')
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=button_close.GetId())
        hbox2.Add(button_close, 1, wx.ALL ,20)
        
        vbox.Add(hbox2, 0, wx.EXPAND|wx.TOP, 10 )
        
        panel.SetSizer(vbox)
        
        vboxtop.Add(panel,1, wx.EXPAND )
        
        self.SetSizer(vboxtop)
        
        self.SetPosition((220, 150))
        
        self.CalcR()
        self.CalcG()
        self.CalcB()        
        self.draw_image()        
        
        
#----------------------------------------------------------------------           
    def OnSelectR(self, event):
        item = event.GetSelection()
        self.r_spec = item
        
        self.CalcR()
        self.draw_image()
        
#----------------------------------------------------------------------           
    def CalcR(self):
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
    def OnLimitMinR(self, event):
        ctl = event.GetEventObject()
        value = ctl.GetValue()
        self.minr = value
        #print 'self.minr=', self.minr
        self.CalcR()
        self.draw_image()
    
#----------------------------------------------------------------------           
    def OnLimitMaxR(self, event):
        ctl = event.GetEventObject()
        value = ctl.GetValue()
        self.maxr = value
        #print 'self.maxr=', self.maxr
        self.CalcR()
        self.draw_image()
        
            
#----------------------------------------------------------------------           
    def OnWeightR(self, event):
        ctl = event.GetEventObject()
        value = ctl.GetValue()
        self.weightr = value
        #print 'self.weightr=', self.weightr
        self.CalcR()
        self.draw_image()
        
#----------------------------------------------------------------------           
    def OnSelectG(self, event):
        item = event.GetSelection()
        self.g_spec = item

        self.CalcG()
        self.draw_image()
        
#----------------------------------------------------------------------           
    def CalcG(self):
        
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
    def OnLimitMinG(self, event):
        ctl = event.GetEventObject()
        value = ctl.GetValue()
        self.ming = value
        #print 'self.ming=', self.ming
        self.CalcG()
        self.draw_image()
        
#----------------------------------------------------------------------           
    def OnLimitMaxG(self, event):
        ctl = event.GetEventObject()
        value = ctl.GetValue()
        self.maxg = value
        #print 'self.maxg=', self.maxg
        self.CalcG()
        self.draw_image()
            
#----------------------------------------------------------------------           
    def OnWeightG(self, event):
        ctl = event.GetEventObject()
        value = ctl.GetValue()
        self.weightg = value
        #print 'self.weightg=', self.weightg
        self.CalcG()
        self.draw_image()
        
#----------------------------------------------------------------------           
    def OnSelectB(self, event):
        item = event.GetSelection()
        self.b_spec = item
        
        self.CalcB()
        self.draw_image()
        
#----------------------------------------------------------------------           
    def CalcB(self):
        
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
    def OnLimitMinB(self, event):
        ctl = event.GetEventObject()
        value = ctl.GetValue()
        self.minb = value
        #print 'self.minb=', self.minb
        self.CalcB()
        self.draw_image()
        
#----------------------------------------------------------------------           
    def OnLimitMaxB(self, event):
        ctl = event.GetEventObject()
        value = ctl.GetValue()
        self.maxb = value
        #print 'self.maxb=', self.maxb
        self.CalcB()
        self.draw_image()
            
#----------------------------------------------------------------------           
    def OnWeightB(self, event):
        ctl = event.GetEventObject()
        value = ctl.GetValue()
        self.weightb = value
        #print 'self.weightb=', self.weightb
        self.CalcB()
        self.draw_image()
        
#----------------------------------------------------------------------           
    def OnShowInfo(self, event):
        if self.show_info_cb.GetValue():
            self.show_info = 1
        else: 
            self.show_info = 0
        
        self.draw_image()
        
#----------------------------------------------------------------------        
    def draw_image(self):
               
               
        fig = self.RGBImagePanel.get_figure()
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

        wildcard = 'Portable Network Graphics (*.png)|*.png|Adobe PDF Files (*.pdf)|*.pdf'               
        SaveFileName = wx.FileSelector('Save Plot', default_extension='png', 
                                   wildcard=wildcard, parent=self, flags=wx.SAVE|wx.OVERWRITE_PROMPT) 
        
        if not SaveFileName: 
            return      
        path, ext = os.path.splitext(SaveFileName) 
        ext = ext[1:].lower() 
        
       
        if ext != 'png' and ext != 'pdf': 
            error_message = ( 
                  'Only the PNG and PDF image formats are supported.\n' 
                 'A file extension of `png\' or `pdf\' must be used.') 
            wx.MessageBox(error_message, 'Error - Could not save file.'+error_message, 
                  parent=self, style=wx.OK|wx.ICON_ERROR) 
            return 
   

        mtplot.rcParams['pdf.fonttype'] = 42

                            
        fig = self.RGBImagePanel.get_figure()
        fig.savefig(SaveFileName, bbox_inches='tight', pad_inches = 0.0)
                
   

#----------------------------------------------------------------------              
    def OnClose(self, evt):
        self.Destroy()             
        
#---------------------------------------------------------------------- 
class SaveWinP4(wx.Frame):
    
    title = "Save"

    def __init__(self):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=self.title, size=(400, 330))
               
        ico = logos.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White") 

        
        self.com = wx.GetApp().TopWindow.common         
        self.fontsize = self.com.fontsize   
        
        
        path, ext = os.path.splitext(self.com.filename) 
        ext = ext[1:].lower()   
        suffix = "." + ext
        path, fn = os.path.split(self.com.filename)
        filename = fn[:-len(suffix)]
        
        self.path = self.com.path
        self.filename = filename
            
        
        
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        
        panel1 = wx.Panel(self, -1)
        
        gridtop = wx.FlexGridSizer(3, 4, vgap=20, hgap=20)
    
        fontb = wx.SystemSettings_GetFont(wx.SYS_DEFAULT_GUI_FONT)
        fontb.SetWeight(wx.BOLD)
        
        
        st1 = wx.StaticText(panel1, -1, 'Save',  style=wx.ALIGN_LEFT)
        st1.SetFont(fontb)
        st2 = wx.StaticText(panel1, -1, '.pdf',  style=wx.ALIGN_LEFT)
        st2.SetFont(fontb)
        st3 = wx.StaticText(panel1, -1, '.png',  style=wx.ALIGN_LEFT)
        st3.SetFont(fontb)
        st3a = wx.StaticText(panel1, -1, '.csv',  style=wx.ALIGN_LEFT)
        st3a.SetFont(fontb)
                
        st4 = wx.StaticText(panel1, -1, '_spectrum',  style=wx.ALIGN_LEFT)
        st4.SetFont(self.com.font)
        
        self.cb1 = wx.CheckBox(panel1, -1, '')
        self.cb1.SetFont(self.com.font)
        self.cb1.Set3StateValue(wx.CHK_CHECKED)
        
        self.cb2 = wx.CheckBox(panel1, -1, '')
        self.cb2.SetFont(self.com.font)
        self.cb2a = wx.CheckBox(panel1, -1, '')
        self.cb2a.SetFont(self.com.font)
        
        st5 = wx.StaticText(panel1, -1, '_image',  style=wx.ALIGN_LEFT)
        st5.SetFont(self.com.font)
        
        self.cb3 = wx.CheckBox(panel1, -1, '')
        self.cb3.SetFont(self.com.font)
        self.cb3.Set3StateValue(wx.CHK_CHECKED)
        
        self.cb4 = wx.CheckBox(panel1, -1, '')
        self.cb4.SetFont(self.com.font)

        
        gridtop.Add(st1, 0)
        gridtop.Add(st2, 0)
        gridtop.Add(st3, 0)
        gridtop.Add(st3a, 0)
                
        gridtop.Add(st4, 0)
        gridtop.Add(self.cb1, 0)
        gridtop.Add(self.cb2, 0)  
        gridtop.Add(self.cb2a, 0)            
  
        gridtop.Add(st5, 0)
        gridtop.Add(self.cb3, 0)
        gridtop.Add(self.cb4, 0)  
        gridtop.Add(wx.StaticText(panel1, -1, ' '))
        
        panel1.SetSizer(gridtop)
        
        panel2 = wx.Panel(self, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        
        hbox0 = wx.BoxSizer(wx.HORIZONTAL)
        
        st0 = wx.StaticText(panel2, -1, 'Filename: ')
        st0.SetFont(self.com.font)
        self.tc_savefn = wx.TextCtrl(panel2, -1,  style=wx.TE_RICH, size =((100, -1)),
                                         value=self.filename)
        self.tc_savefn.SetFont(self.com.font)
        hbox0.Add(st0, 0, flag = wx.ALIGN_CENTER_VERTICAL)
        hbox0.Add(self.tc_savefn,1, wx.EXPAND|wx.LEFT, 2)         
        
        
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
                
        st1 = wx.StaticText(panel2, label='Path: ')
        st1.SetFont(self.com.font)
        self.tc_savepath = wx.TextCtrl(panel2, -1,  style=wx.TE_RICH|wx.TE_READONLY, size =((100, -1)),
                                         value=self.path)
        self.tc_savepath.SetFont(self.com.font)
        hbox1.Add(st1, 0, wx.ALIGN_CENTER_VERTICAL|wx.RIGHT, 20)
        hbox1.Add(self.tc_savepath, 1, wx.EXPAND|wx.LEFT, 2)  
        
        button_path = wx.Button(panel2, -1, 'Browse...')
        self.Bind(wx.EVT_BUTTON, self.OnBrowseDir, id=button_path.GetId())
        hbox1.Add(button_path, 0)
        
        
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        button_save = wx.Button(panel2, -1, 'Save')
        self.Bind(wx.EVT_BUTTON, self.OnSave, id=button_save.GetId())
        hbox2.Add(button_save, 0, wx.ALIGN_RIGHT)
        
        button_cancel = wx.Button(panel2, -1, 'Cancel')
        self.Bind(wx.EVT_BUTTON, self.OnCancel, id=button_cancel.GetId())
        hbox2.Add(button_cancel, 0, wx.ALIGN_RIGHT|wx.LEFT,10)
        
        vbox1.Add(hbox0, 1, wx.EXPAND)
        vbox1.Add((0,5))
        vbox1.Add(hbox1, 1, wx.EXPAND)
        vbox1.Add((0,20))
        vbox1.Add(hbox2, 0, wx.ALIGN_RIGHT)
        panel2.SetSizer(vbox1)
        
        
        vboxtop.Add(panel1, 1, wx.ALL | wx.EXPAND, 20)
        vboxtop.Add(panel2, 0, wx.ALL | wx.EXPAND, 20) 
        
        self.SetSizer(vboxtop)    
        
        
#----------------------------------------------------------------------        
    def OnBrowseDir(self, evt):
        
        dialog = wx.DirDialog(None, "Choose a directory",
                               style=wx.DD_DIR_MUST_EXIST,
                               defaultPath=self.path)

        if dialog.ShowModal() == wx.ID_OK:
            directory = dialog.GetPath()
            
            self.path = directory
            
            self.tc_savepath.SetValue(self.path)
            
            
                
#----------------------------------------------------------------------        
    def OnSave(self, evt):
        
        self.filename = self.tc_savefn.GetValue()
        
        sp_pdf = self.cb1.GetValue()
        sp_png = self.cb2.GetValue()
        sp_csv = self.cb2a.GetValue()
        im_pdf = self.cb3.GetValue()
        im_png = self.cb4.GetValue()
        
        self.Destroy() 
        wx.GetApp().TopWindow.page4.Save(self.filename, self.path,
                                         spec_png = sp_png, 
                                         spec_pdf = sp_pdf, 
                                         spec_csv = sp_csv,
                                         img_png = im_png, 
                                         img_pdf = im_pdf)
        

#---------------------------------------------------------------------- 
    def OnCancel(self, evt):
        self.Destroy()  
                
                     

        
""" ------------------------------------------------------------------------------------------------"""
class PageCluster(wx.Panel):
    def __init__(self, parent, common, data_struct, stack, anlz):
        wx.Panel.__init__(self, parent)
        
        self.data_struct = data_struct
        self.stk = stack
        self.com = common
        self.anlz = anlz
        
        self.SetBackgroundColour("White")
        
        
        self.selcluster = 1
        self.numclusters = 0
        self.init_nclusters = 5
        self.wo_1st_pca = 0
        self.sigma_split = 0
        
        self.MakeColorTable()             
        self.fontsize = self.com.fontsize
        
        
        #panel 1
        panel1 = wx.Panel(self, -1)
        sizer1 = wx.StaticBoxSizer(wx.StaticBox(panel1, -1, 'Cluster analysis'), wx.VERTICAL)
        self.button_calcca = wx.Button(panel1, -1, 'Calculate Clusters', (10, 10))
        self.button_calcca.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnCalcClusters, id=self.button_calcca.GetId())   
        self.button_calcca.Disable()     
        self.button_scatterplots = wx.Button(panel1, -1, 'Show scatter plots...', (10, 10))
        self.button_scatterplots.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnShowScatterplots, id=self.button_scatterplots.GetId())
        self.button_scatterplots.Disable()
        self.button_savecluster = wx.Button(panel1, -1, 'Save CA Results...', (10, 10))
        self.button_savecluster.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnSave, id=self.button_savecluster.GetId())
        self.button_savecluster.Disable()
        
        
        hbox11 = wx.BoxSizer(wx.HORIZONTAL)
        text1 = wx.StaticText(panel1, -1, 'Number of clusters',  style=wx.ALIGN_LEFT)
        text1.SetFont(self.com.font)
        self.nclusterspin = wx.SpinCtrl(panel1, -1, '',  size= (70, -1), style=wx.ALIGN_LEFT)
        self.nclusterspin.SetRange(2,20)
        self.nclusterspin.SetValue(self.init_nclusters)
        self.Bind(wx.EVT_SPINCTRL, self.OnNClusterspin, self.nclusterspin)
        hbox11.Add((20,0))
        hbox11.Add(text1, 0, wx.TOP, 20)
        hbox11.Add((10,0))
        hbox11.Add(self.nclusterspin, 0, wx.TOP, 15)  
        hbox11.Add((20,0))    
        
        hbox11a = wx.BoxSizer(wx.HORIZONTAL)
        hbox11a.Add((20,0))
        text1a = wx.StaticText(panel1, label="Number of clusters found")
        self.ntc_clusters_found = wx.lib.intctrl.IntCtrl( panel1, size=( 45, -1 ), 
                                                     style = wx.TE_CENTRE|wx.TE_READONLY, 
                                                     value = self.numclusters)

        hbox11a.Add(text1a, 0, wx.TOP, 10)
        hbox11a.Add((5,0))
        hbox11a.Add(self.ntc_clusters_found, 0, wx.TOP, 5)   
            
        
        hbox12 = wx.BoxSizer(wx.HORIZONTAL)
        hbox12.Add((20,0))
        self.remove1stpcacb = wx.CheckBox(panel1, -1, 'Reduce thickness effects')
        self.remove1stpcacb.SetFont(self.com.font)
        self.Bind(wx.EVT_CHECKBOX, self.OnRemove1stpca, self.remove1stpcacb)
        hbox12.Add(self.remove1stpcacb, 0, wx.EXPAND|wx.TOP, 10)
        hbox12.Add((20,0))

        hbox13 = wx.BoxSizer(wx.HORIZONTAL)
        hbox13.Add((20,0))
        self.cb_splitclusters = wx.CheckBox(panel1, -1, 'Divide clusters with large Sigma')
        self.cb_splitclusters.SetFont(self.com.font)
        self.Bind(wx.EVT_CHECKBOX, self.OnSplitClusters, self.cb_splitclusters)
        hbox13.Add(self.cb_splitclusters, 0, wx.EXPAND, 15)
        hbox13.Add((20,0))        
        
        sizer1.Add((0,10)) 
        sizer1.Add(self.button_calcca, 0, wx.EXPAND|wx.LEFT|wx.RIGHT, 20)
        sizer1.Add(hbox11, 0, wx.EXPAND)
        sizer1.Add(hbox11a, 0, wx.EXPAND)
        sizer1.Add(hbox12, 0, wx.EXPAND)
        sizer1.Add(hbox13, 0, wx.EXPAND|wx.TOP, 3)
                
        sizer1.Add((0,5))        
        sizer1.Add(wx.StaticLine(panel1), 0, wx.ALL|wx.EXPAND, 5)        
        sizer1.Add((0,5)) 
      
        sizer1.Add(self.button_scatterplots, 0, wx.EXPAND|wx.LEFT|wx.RIGHT, 20)
        sizer1.Add(self.button_savecluster, 0, wx.EXPAND|wx.LEFT|wx.RIGHT, 20)       
        sizer1.Add((0,5))
         
        panel1.SetSizer(sizer1)
        
        
                
        #panel 2        
        panel2 = wx.Panel(self, -1)
        vbox2 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_clustercomp = wx.TextCtrl(panel2, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_clustercomp.SetValue("Composite cluster image")        
        self.tc_clustercomp.SetFont(self.com.font)
        
        i2panel = wx.Panel(panel2, -1, style = wx.SUNKEN_BORDER)
        self.ClusterImagePan = wxmpl.PlotPanel(i2panel, -1, size =(PlotH, PlotH), cursor=False, crosshairs=True, location=False, zoom=False)                              
        wxmpl.EVT_POINT(i2panel, self.ClusterImagePan.GetId(), self.OnPointClusterImage)   
        vbox2.Add(self.tc_clustercomp, 0, wx.EXPAND) 
        vbox2.Add(i2panel, 0)   

        panel2.SetSizer(vbox2)
        
        
        #panel 3 
        panel3 = wx.Panel(self, -1)
        vbox3 = wx.BoxSizer(wx.VERTICAL)
        panel3.SetBackgroundColour("White")  
        fgs = wx.FlexGridSizer(2, 3, 0, 0)
        
        self.tc_cluster = wx.TextCtrl(panel3, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE, size=(250,-1))
        self.tc_cluster.SetValue("Cluster ")
        self.tc_cluster.SetFont(self.com.font)

        i3panel = wx.Panel(panel3, -1, style = wx.SUNKEN_BORDER)
        self.ClusterIndvImagePan = wxmpl.PlotPanel(i3panel, -1, size =(PlotH*0.73, PlotH*0.73), cursor=False, crosshairs=False, location=False, zoom=False)
    
        self.slidershow = wx.Slider(panel3, -1, self.selcluster, 1, 20, style=wx.SL_LEFT|wx.SL_VERTICAL)   
        self.slidershow.Disable()    
        self.slidershow.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnClusterScroll, self.slidershow)
        
        vbox31 = wx.BoxSizer(wx.VERTICAL)               
        self.clusterspin = wx.SpinButton(panel3, -1, size = ((8,-1)), style=wx.SP_ARROW_KEYS)
        self.Bind(wx.EVT_SPIN_UP, self.OnClusterSpinUp, self.clusterspin)
        self.Bind(wx.EVT_SPIN_DOWN, self.OnClusterSpinDown, self.clusterspin)
        
        vbox31.Add((0,3))
        vbox31.Add(self.slidershow, 1,  wx.EXPAND) 
        vbox31.Add(self.clusterspin, 0,  wx.EXPAND)         
        
        
        text3 = wx.StaticText(panel3, -1, 'Cluster Distance Map',  style=wx.ALIGN_LEFT)
        text3.SetFont(self.com.font)
        i4panel = wx.Panel(panel3, -1, style = wx.SUNKEN_BORDER)
        self.ClusterDistMapPan = wxmpl.PlotPanel(i4panel, -1, size =(PlotH*0.73, PlotH*0.73), cursor=False, crosshairs=False, location=False, zoom=False)
        
          
        fgs.AddMany([(self.tc_cluster), (wx.StaticText(panel3, -1, ' ')), (text3, 0, wx.LEFT, 15), 
                     (i3panel), (vbox31, 0,  wx.EXPAND), (i4panel, 0, wx.LEFT, 20)])

        vbox3.Add(fgs)

        panel3.SetSizer(vbox3)
        
        
        #panel 4 
        panel4 = wx.Panel(self, -1)
        vbox4 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_clustersp = wx.TextCtrl(panel4, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_clustersp.SetValue("Cluster spectrum")
        self.tc_clustersp.SetFont(self.com.font)        

        i5panel = wx.Panel(panel4, -1, style = wx.SUNKEN_BORDER)
        self.ClusterSpecPan = wxmpl.PlotPanel(i5panel, -1, size =(PlotW, PlotH), cursor=False, crosshairs=False, location=False, zoom=False)
        
        vbox4.Add(self.tc_clustersp, 0, wx.EXPAND)        
        vbox4.Add(i5panel, 0)

        panel4.SetSizer(vbox4)
        
 

        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        gridtop = wx.FlexGridSizer(2, 2, vgap=10, hgap=50)
        gridtop.Add(panel1, 0, wx.LEFT|wx.TOP, 20)
        gridtop.Add(panel3, 0)
        
        gridtop.Add(panel2, 0)
        gridtop.Add(panel4, 0)
              
        vboxtop.Add((0,10))
        vboxtop.Add(gridtop, 0, wx.LEFT, 20)
         
        
        self.SetSizer(vboxtop) 
        
        
        
        
#----------------------------------------------------------------------
    def OnCalcClusters(self, event):
       
        wx.BeginBusyCursor()
        self.calcclusters = False  
        
        if True:
        #try: 
            self.CalcClusters()
            
            self.calcclusters = True
            
            self.selcluster = 1
            self.slidershow.SetValue(self.selcluster)
            self.slidershow.SetMax(self.numclusters)
            
            self.showClusterImage()
            self.showClusterSpectrum()
            self.showIndvClusterImage()     
            self.showClusterDistanceMap()
            self.com.cluster_calculated = 1       
            wx.EndBusyCursor() 
            
#        except:
#            self.com.cluster_calculated = 0
#            wx.EndBusyCursor()      
            
        wx.GetApp().TopWindow.refresh_widgets()
            
#----------------------------------------------------------------------        
    def OnNClusterspin(self, event):
        num = event.GetInt()
        self.init_nclusters = num
                       
 
#----------------------------------------------------------------------        
    def OnClusterScroll(self, event):
        sel = event.GetInt()
        self.selcluster = sel
        if self.com.cluster_calculated == 1:
            self.showClusterSpectrum()
            self.showIndvClusterImage()
            
#----------------------------------------------------------------------            
    def OnClusterSpinUp(self, event):
        if (self.com.cluster_calculated == 1) and (self.selcluster > 1):
            self.selcluster = self.selcluster - 1
            self.slidershow.SetValue(self.selcluster)

            self.showClusterSpectrum()
            self.showIndvClusterImage()
            
#----------------------------------------------------------------------            
    def OnClusterSpinDown(self, event):
        if (self.com.cluster_calculated == 1) and (self.selcluster < self.numclusters):
            self.selcluster = self.selcluster + 1
            self.slidershow.SetValue(self.selcluster) 
            
            self.showClusterSpectrum()
            self.showIndvClusterImage()
                       
#----------------------------------------------------------------------  
    def OnPointClusterImage(self, evt):
        x = evt.xdata
        y = evt.ydata
        
        if self.com.cluster_calculated == 1:     
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
            
            self.slidershow.SetValue(self.selcluster)
            
            self.showClusterSpectrum()
            self.showIndvClusterImage()
            
                       
#----------------------------------------------------------------------
    def CalcClusters(self):
        nclusters = self.anlz.calculate_clusters(self.init_nclusters, self.wo_1st_pca, self.sigma_split)
        self.numclusters = nclusters
        self.ntc_clusters_found.SetValue(self.numclusters)
    
        
        
        
#----------------------------------------------------------------------     
#Show composite cluster image  
    def showClusterImage(self):       
        
        self.clusterimage = self.anlz.cluster_indices  
        
        #print self.selpca

        fig = self.ClusterImagePan.get_figure()
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        axes = fig.gca()
        
        mtplot.rcParams['font.size'] = self.fontsize
        
        
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

        fig = self.ClusterIndvImagePan.get_figure()
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        axes = fig.gca()
        
        mtplot.rcParams['font.size'] = self.fontsize        
        
        im = axes.imshow(indvclusterimage, cmap=self.clusterclrmap2, norm=self.bnorm2)
        axes.axis("off")

        self.ClusterIndvImagePan.draw()
        
        self.tc_cluster.SetValue("Cluster " + str(self.selcluster))
         
   
#----------------------------------------------------------------------     
    def showClusterSpectrum(self):
        
         
        clusterspectrum = self.anlz.clusterspectra[self.selcluster-1, ]

       
        fig = self.ClusterSpecPan.get_figure()
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        
        mtplot.rcParams['font.size'] = self.fontsize
        
        
        if self.selcluster >= self.maxclcolors:
            clcolor = self.colors[self.maxclcolors-1]
        else:
            clcolor = self.colors[self.selcluster-1]
        
        specplot = axes.plot(self.anlz.stack.ev,clusterspectrum, color = clcolor)
        
        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Optical Density')
        
        
        self.ClusterSpecPan.draw()
        
        self.tc_clustersp.SetValue("Cluster " + str(self.selcluster)+ " spectrum"  ) 
        
        
#----------------------------------------------------------------------     
    def showClusterDistanceMap(self):       
        
        mapimage = self.anlz.cluster_distances
        
        #print self.selpca

        fig = self.ClusterDistMapPan.get_figure()
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        axes = fig.gca()
        
#         divider = make_axes_locatable(axes)
#         axcb = divider.new_horizontal(size="3%", pad=0.03)  
#         fig.add_axes(axcb)  
#         axes.set_position([0.03,0.03,0.8,0.94])
               
        mtplot.rcParams['font.size'] = self.fontsize
        
        
        im = axes.imshow(mapimage, cmap=mtplot.cm.get_cmap('gray'))
        
        #cbar = axes.figure.colorbar(im, orientation='vertical',cax=axcb) 
        
        axes.axis("off")
       
        self.ClusterDistMapPan.draw()
        
#----------------------------------------------------------------------           
    def OnRemove1stpca(self, event):
        if self.remove1stpcacb.GetValue():
            self.wo_1st_pca = 1
        else: self.wo_1st_pca = 0

#----------------------------------------------------------------------           
    def OnSplitClusters(self, event):
        if self.cb_splitclusters.GetValue():
            self.sigma_split = 1
        else: self.sigma_split = 0        
        
#----------------------------------------------------------------------    
    def OnSave(self, event):     
               

        SaveWinP3().Show()
        
        
#----------------------------------------------------------------------    
    def Save(self, filename, path, spec_png = True, spec_pdf = False, spec_csv = False,
             img_png = True, img_pdf = False, 
             indimgs_png = True, indimgs_pdf = False,
             scatt_png = True, scatt_pdf = False): 
        
        self.SaveFileName = os.path.join(path,filename)
        
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas    
        mtplot.rcParams['pdf.fonttype'] = 42
   
        try: 
            
            if img_png:
                ext = 'png'
                suffix = "." + ext
            
                fig = mtplot.figure.Figure(figsize = (float(self.stk.n_rows)/10, float(self.stk.n_cols)/10))
                canvas = FigureCanvas(fig)
                fig.clf()
                fig.add_axes((0.0,0.0,1.0,1.0))
                axes = fig.gca()      
                mtplot.rcParams['font.size'] = self.fontsize        
        
                im = axes.imshow(self.clusterimage, cmap=self.clusterclrmap1, norm=self.bnorm1)
                axes.axis("off")
                
                fileName_caimg = self.SaveFileName+"_CAcimg."+ext       
                fig.savefig(fileName_caimg, dpi=ImgDpi, pad_inches = 0.0)
                
            
            if img_pdf:
                ext = 'pdf'
                suffix = "." + ext
            
                fig = mtplot.figure.Figure(figsize = (float(self.stk.n_rows)/30, float(self.stk.n_cols)/30))
                canvas = FigureCanvas(fig)
                fig.clf()
                fig.add_axes((0.0,0.0,1.0,1.0))
                axes = fig.gca()      
                mtplot.rcParams['font.size'] = self.fontsize        
        
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

                    fig = mtplot.figure.Figure(figsize =(float(self.stk.n_rows)/10, float(self.stk.n_cols)/10))
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.0,0.0,1.0,1.0))
                    axes = fig.gca()      
                    mtplot.rcParams['font.size'] = self.fontsize        
                    im = axes.imshow(indvclusterimage, cmap=self.clusterclrmap2, norm=self.bnorm2)
                    axes.axis("off")
                   
                    fileName_img = self.SaveFileName+"_CAimg_" +str(i+1)+"."+ext               
                    fig.savefig(fileName_img, dpi=ImgDpi, pad_inches = 0.0)
                
            if spec_png:
                for i in range (self.numclusters):
                   
                    clusterspectrum = self.anlz.clusterspectra[i, ]
                    fig = mtplot.figure.Figure(figsize =(PlotW, PlotH))
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.15,0.15,0.75,0.75))
                    axes = fig.gca()
                    mtplot.rcParams['font.size'] = self.fontsize
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
                    cname = 'CAspectrum_' +str(i+1)
                    self.stk.write_csv(fileName_spec, self.anlz.stack.ev, clusterspectrum, cname=cname)
                                                     
                
            ext = 'pdf'
            suffix = "." + ext
                
            if indimgs_pdf:
                for i in range (self.numclusters):
              
                    indvclusterimage = npy.zeros((self.anlz.stack.n_cols, self.anlz.stack.n_rows))+20.      
                    ind = npy.where(self.anlz.cluster_indices == i)    
                    colorcl = min(i,9)
                    indvclusterimage[ind] = colorcl

                    fig = mtplot.figure.Figure(figsize =(float(self.stk.n_rows)/30, float(self.stk.n_cols)/30))
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.0,0.0,1.0,1.0))
                    axes = fig.gca()      
                    mtplot.rcParams['font.size'] = self.fontsize        
                    im = axes.imshow(indvclusterimage, cmap=self.clusterclrmap2, norm=self.bnorm2)
                    axes.axis("off")
                   
                    fileName_img = self.SaveFileName+"_CAimg_" +str(i+1)+"."+ext               
                    fig.savefig(fileName_img, dpi=300, pad_inches = 0.0)
                
            if spec_pdf:
                for i in range (self.numclusters):
                   
                    clusterspectrum = self.anlz.clusterspectra[i, ]
                    fig = mtplot.figure.Figure(figsize =(PlotW, PlotH))
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.15,0.15,0.75,0.75))
                    axes = fig.gca()
                    mtplot.rcParams['font.size'] = self.fontsize
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
   
            wx.MessageBox('Could not save file: %s' % err, 'Error', 
                          parent=self, style=wx.OK|wx.ICON_ERROR) 
  
  
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
            wx.BeginBusyCursor()
            
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
                fig = mtplot.figure.Figure(figsize =(6.0,plotsize*nplotsrows))
                fig.subplots_adjust(wspace = 0.4, hspace = 0.4)
            else:
                fig = mtplot.figure.Figure(figsize =(3.0,2.5))
                fig.subplots_adjust(bottom = 0.2, left = 0.2)
                
            canvas = FigureCanvas(fig)
            #axes = fig.gca()
            mtplot.rcParams['font.size'] = 6   
            
                                   
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
            mtplot.rcParams['pdf.fonttype'] = 42
            fig.savefig(fileName_sct)
            
            wx.EndBusyCursor()
     
            
        except IOError, e:
            wx.EndBusyCursor()
            if e.strerror:
                err = e.strerror 
            else: 
                err = e 
   
            wx.MessageBox('Could not save file: %s' % err, 'Error', 
                          parent=self, style=wx.OK|wx.ICON_ERROR) 
            
        
#----------------------------------------------------------------------     
    def MakeColorTable(self):
        self.maxclcolors = 11
        colors_i = npy.linspace(0,self.maxclcolors,self.maxclcolors+1)

       
        self.colors=['#0000FF','#FF0000','#DFE32D','#36F200','#B366FF',
                '#FF470A','#33FFFF','#006600','#CCCC99','#993300',
                '#000000']


        
        self.clusterclrmap1=mtplot.colors.LinearSegmentedColormap.from_list('clustercm',self.colors)
     
        self.bnorm1 = mtplot.colors.BoundaryNorm(colors_i, self.clusterclrmap1.N)
        
        colors_i = npy.linspace(0,self.maxclcolors+2,self.maxclcolors+3)
        
        #use black color for clusters > maxclcolors, the other 2 colors are for background      
        colors2=['#0000FF','#FF0000','#DFE32D','#36F200','#B366FF',
                '#FF470A','#33FFFF','#006600','#CCCC99','#993300',
                '#000000','#FFFFFF','#EEEEEE']
        
        self.clusterclrmap2=mtplot.colors.LinearSegmentedColormap.from_list('clustercm2',colors2)
     
        self.bnorm2 = mtplot.colors.BoundaryNorm(colors_i, self.clusterclrmap2.N)
        

        
#----------------------------------------------------------------------     
    def OnShowScatterplots(self, evt):    
        Scatterplots(self.com, self.anlz).Show()
        
#---------------------------------------------------------------------- 
class Scatterplots(wx.Frame):

    title = "Scatter plots"
 
    def __init__(self, common, analz):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=self.title, size=(600, 550))
        
        self.SetBackgroundColour("White")
        
        ico = logos.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White") 
        
        self.colors = wx.GetApp().TopWindow.page3.colors
        
        self.com = common       
        self.fontsize = self.com.fontsize
        
        self.anlz = analz
        self.numsigpca = self.anlz.numsigpca  
        self.ncols = self.anlz.stack.n_cols
        self.nrows = self.anlz.stack.n_rows
        
        self.od_reduced = self.anlz.pcaimages[:,:,0:self.numsigpca]        
        self.od_reduced = npy.reshape(self.od_reduced, (self.ncols*self.nrows,self.numsigpca), order='F')
        
        self.clindices = self.anlz.cluster_indices
        self.clindices = npy.reshape(self.clindices, (self.ncols*self.nrows), order='F')
        self.numclusters = wx.GetApp().TopWindow.page3.numclusters
             
                     
        self.pca_y = 1
        self.pca_x = 1
        
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)  
          
        grid1 = wx.FlexGridSizer(2, 2)

        i1panel = wx.Panel(panel, -1, style = wx.SUNKEN_BORDER)
        self.ScatterPPanel = wxmpl.PlotPanel(i1panel, -1, size=(5.0, 4.0), cursor=False, crosshairs=False, location=False, zoom=False)
        
        self.slidershow_y = wx.Slider(panel, -1, self.pca_y, 1, self.numsigpca, style=wx.SL_RIGHT|wx.SL_VERTICAL|wx.SL_LABELS|wx.SL_INVERSE)
        self.slidershow_y.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnSliderScroll_y, self.slidershow_y)
        
        grid1.Add(self.slidershow_y, 0, wx.EXPAND)
        grid1.Add(i1panel, 0)
               
        
        self.slidershow_x = wx.Slider(panel, -1, self.pca_x, 1, self.numsigpca, style=wx.SL_TOP|wx.SL_LABELS|wx.SL_HORIZONTAL)
        self.slidershow_x.SetFocus()        
        self.Bind(wx.EVT_SCROLL, self.OnSliderScroll_x, self.slidershow_x)
        
        grid1.Add(wx.StaticText(panel, -1, ''))
        grid1.Add(self.slidershow_x, 0,  wx.EXPAND)
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        
        button_close = wx.Button(panel, -1, 'Close')
        button_close.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=button_close.GetId())
        hbox.Add(button_close, -1, wx.ALIGN_RIGHT|wx.RIGHT, 40)

        
        vbox.Add(grid1, 0, wx.EXPAND|wx.ALL,20)
        vbox.Add(hbox, -1, wx.ALIGN_RIGHT)
                
        panel.SetSizer(vbox)
        
        vboxtop.Add(panel,-1, wx.EXPAND|wx.ALIGN_CENTER_VERTICAL)
        
        self.SetSizer(vboxtop)
        
        self.SetPosition((220, 150))
        
        self.draw_scatterplot()
        
        
#----------------------------------------------------------------------        
    def OnSliderScroll_x(self, event):
        self.pca_x = event.GetInt()
        
        self.draw_scatterplot()

#----------------------------------------------------------------------        
    def OnSliderScroll_y(self, event):
        self.pca_y = event.GetInt()
        
        self.draw_scatterplot()        

      
#----------------------------------------------------------------------        
    def draw_scatterplot(self):
                
        x_comp = self.od_reduced[:,self.pca_x-1]
        y_comp = self.od_reduced[:,self.pca_y-1]
        
        fig = self.ScatterPPanel.get_figure()
        fig.clf()
        axes = fig.gca()
        
        mtplot.rcParams['font.size'] = self.fontsize        
        
        for i in range(self.numclusters):
            thiscluster = npy.where(self.clindices == i)
            axes.plot(x_comp[thiscluster], y_comp[thiscluster],'.',color=self.colors[i],alpha=0.5)
            
        axes.set_xlabel('Component '+str(self.pca_x))
        axes.set_ylabel('Component '+str(self.pca_y))
      
        self.ScatterPPanel.draw()
       

  
#----------------------------------------------------------------------              
    def OnClose(self, evt):
        self.Destroy()       
        
        
        
#---------------------------------------------------------------------- 
class SaveWinP3(wx.Frame):
    
    title = "Save"

    def __init__(self):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=self.title, size=(400, 330))
               
        ico = logos.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White") 

        
        self.com = wx.GetApp().TopWindow.common         
        self.fontsize = self.com.fontsize   
        
        
        path, ext = os.path.splitext(self.com.filename) 
        ext = ext[1:].lower()   
        suffix = "." + ext
        path, fn = os.path.split(self.com.filename)
        filename = fn[:-len(suffix)]
        
        self.path = self.com.path
        self.filename = filename
            
        
        
        vboxtop = wx.BoxSizer(wx.VERTICAL)
                
        panel1 = wx.Panel(self, -1)
        
        gridtop = wx.FlexGridSizer(5, 4, vgap=20, hgap=20)
    
        fontb = wx.SystemSettings_GetFont(wx.SYS_DEFAULT_GUI_FONT)
        fontb.SetWeight(wx.BOLD)
        
        
        st1 = wx.StaticText(panel1, -1, 'Save',  style=wx.ALIGN_LEFT)
        st1.SetFont(fontb)
        st2 = wx.StaticText(panel1, -1, '.pdf',  style=wx.ALIGN_LEFT)
        st2.SetFont(fontb)
        st3 = wx.StaticText(panel1, -1, '.png',  style=wx.ALIGN_LEFT)
        st3.SetFont(fontb)
        st3a = wx.StaticText(panel1, -1, '.csv',  style=wx.ALIGN_LEFT)
        st3a.SetFont(fontb)
                
        st4 = wx.StaticText(panel1, -1, '_spectra',  style=wx.ALIGN_LEFT)
        st4.SetFont(self.com.font)
        
        self.cb1 = wx.CheckBox(panel1, -1, '')
        self.cb1.SetFont(self.com.font)
        self.cb1.Set3StateValue(wx.CHK_CHECKED)
        
        self.cb2 = wx.CheckBox(panel1, -1, '')
        self.cb2.SetFont(self.com.font)
        
        self.cb2a = wx.CheckBox(panel1, -1, '')
        self.cb2a.SetFont(self.com.font)
        
        st5 = wx.StaticText(panel1, -1, '_composite_images',  style=wx.ALIGN_LEFT)
        st5.SetFont(self.com.font)
        
        self.cb3 = wx.CheckBox(panel1, -1, '')
        self.cb3.SetFont(self.com.font)
        self.cb3.Set3StateValue(wx.CHK_CHECKED)
        
        self.cb4 = wx.CheckBox(panel1, -1, '')
        self.cb4.SetFont(self.com.font)
        
        
        st6 = wx.StaticText(panel1, -1, '_individual_images',  style=wx.ALIGN_LEFT)
        st6.SetFont(self.com.font)
        
        self.cb5 = wx.CheckBox(panel1, -1, '')
        self.cb5.SetFont(self.com.font)
        self.cb5.Set3StateValue(wx.CHK_CHECKED)
        
        self.cb6 = wx.CheckBox(panel1, -1, '')
        self.cb6.SetFont(self.com.font)
        
        st7 = wx.StaticText(panel1, -1, '_scatter_plots',  style=wx.ALIGN_LEFT)
        st7.SetFont(self.com.font)
        
        self.cb7 = wx.CheckBox(panel1, -1, '')
        self.cb7.SetFont(self.com.font)
        self.cb7.Set3StateValue(wx.CHK_CHECKED)
        
        self.cb8 = wx.CheckBox(panel1, -1, '')
        self.cb8.SetFont(self.com.font)

        
        gridtop.Add(st1, 0)
        gridtop.Add(st2, 0)
        gridtop.Add(st3, 0)
        gridtop.Add(st3a, 0)
                
        gridtop.Add(st4, 0)
        gridtop.Add(self.cb1, 0)
        gridtop.Add(self.cb2, 0)              
        gridtop.Add(self.cb2a, 0)    
          
        gridtop.Add(st5, 0)
        gridtop.Add(self.cb3, 0)
        gridtop.Add(self.cb4, 0)  
        gridtop.Add(wx.StaticText(panel1, -1, ' '))
        
        gridtop.Add(st6, 0)
        gridtop.Add(self.cb5, 0)
        gridtop.Add(self.cb6, 0)  
        gridtop.Add(wx.StaticText(panel1, -1, ' '))
        
        gridtop.Add(st7, 0)
        gridtop.Add(self.cb7, 0)
        gridtop.Add(self.cb8, 0) 
        gridtop.Add(wx.StaticText(panel1, -1, ' '))
        
        panel1.SetSizer(gridtop)
        
        panel2 = wx.Panel(self, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        
        hbox0 = wx.BoxSizer(wx.HORIZONTAL)
        
        st0 = wx.StaticText(panel2, -1, 'Filename: ')
        st0.SetFont(self.com.font)
        self.tc_savefn = wx.TextCtrl(panel2, -1,  style=wx.TE_RICH, size =((100, -1)),
                                         value=self.filename)
        self.tc_savefn.SetFont(self.com.font)
        hbox0.Add(st0, 0, flag = wx.ALIGN_CENTER_VERTICAL)
        hbox0.Add(self.tc_savefn,1, wx.EXPAND|wx.LEFT, 2)         
        
        
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
                
        st1 = wx.StaticText(panel2, label='Path: ')
        st1.SetFont(self.com.font)
        self.tc_savepath = wx.TextCtrl(panel2, -1,  style=wx.TE_RICH|wx.TE_READONLY, size =((100, -1)),
                                         value=self.path)
        self.tc_savepath.SetFont(self.com.font)
        hbox1.Add(st1, 0, wx.ALIGN_CENTER_VERTICAL|wx.RIGHT, 20)
        hbox1.Add(self.tc_savepath, 1, wx.EXPAND|wx.LEFT, 2)  
        
        button_path = wx.Button(panel2, -1, 'Browse...')
        self.Bind(wx.EVT_BUTTON, self.OnBrowseDir, id=button_path.GetId())
        hbox1.Add(button_path, 0)
        
        
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        button_save = wx.Button(panel2, -1, 'Save')
        self.Bind(wx.EVT_BUTTON, self.OnSave, id=button_save.GetId())
        hbox2.Add(button_save, 0, wx.ALIGN_RIGHT)
        
        button_cancel = wx.Button(panel2, -1, 'Cancel')
        self.Bind(wx.EVT_BUTTON, self.OnCancel, id=button_cancel.GetId())
        hbox2.Add(button_cancel, 0, wx.ALIGN_RIGHT|wx.LEFT,10)
        
        vbox1.Add(hbox0, 1, wx.EXPAND)
        vbox1.Add((0,5))
        vbox1.Add(hbox1, 1, wx.EXPAND)
        vbox1.Add((0,20))
        vbox1.Add(hbox2, 0, wx.ALIGN_RIGHT)
        panel2.SetSizer(vbox1)
        
        
        vboxtop.Add(panel1, 1, wx.ALL | wx.EXPAND, 20)
        vboxtop.Add(panel2, 0, wx.ALL | wx.EXPAND, 20) 
        
        self.SetSizer(vboxtop)    
        
#----------------------------------------------------------------------        
    def OnBrowseDir(self, evt):
        
        dialog = wx.DirDialog(None, "Choose a directory",
                               style=wx.DD_DIR_MUST_EXIST,
                               defaultPath=self.path)

        if dialog.ShowModal() == wx.ID_OK:
            directory = dialog.GetPath()
            
            self.path = directory
            
            self.tc_savepath.SetValue(self.path)
            
            
                
#----------------------------------------------------------------------        
    def OnSave(self, evt):
        
        self.filename = self.tc_savefn.GetValue()
        
        sp_pdf = self.cb1.GetValue()
        sp_png = self.cb2.GetValue()
        sp_csv = self.cb2a.GetValue()
        im_pdf = self.cb3.GetValue()
        im_png = self.cb4.GetValue()
        indim_pdf = self.cb5.GetValue()
        indim_png = self.cb6.GetValue()
        scatt_pdf = self.cb7.GetValue()
        scatt_png = self.cb8.GetValue()
        
        self.Destroy() 
        wx.GetApp().TopWindow.page3.Save(self.filename, self.path,
                                         spec_png = sp_png, 
                                         spec_pdf = sp_pdf, 
                                         spec_csv = sp_csv,
                                         img_png = im_png, 
                                         img_pdf = im_pdf,
                                         indimgs_png = indim_png, 
                                         indimgs_pdf = indim_pdf,
                                         scatt_png = scatt_png,
                                         scatt_pdf = scatt_pdf)

#---------------------------------------------------------------------- 
    def OnCancel(self, evt):
        self.Destroy()  
                
    
""" ------------------------------------------------------------------------------------------------"""
class PagePCA(wx.Panel):
    def __init__(self, parent, common, data_struct, stack, anlz):
        wx.Panel.__init__(self, parent)
        
        self.SetBackgroundColour("White")
        
        self.com = common 
        self.data_struct = data_struct
        self.stk = stack       
        self.anlz = anlz
        
        self.selpca = 1       
        self.numsigpca = 2
        
        
        pw = PlotW*0.97
        ph = PlotH*0.97
        
          
        self.fontsize = self.com.fontsize
        
  
        #panel 1        
        panel1 = wx.Panel(self, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_PCAcomp = wx.TextCtrl(panel1, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_PCAcomp.SetFont(self.com.font)
        self.tc_PCAcomp.SetValue("PCA component ")
        
        hbox11 = wx.BoxSizer(wx.HORIZONTAL)
   
        i1panel = wx.Panel(panel1, -1, style = wx.SUNKEN_BORDER)
        self.PCAImagePan = wxmpl.PlotPanel(i1panel, -1, size =(ph*1.10, ph), cursor=False, crosshairs=False, location=False, zoom=False)
               
        vbox11 = wx.BoxSizer(wx.VERTICAL)               
        self.slidershow = wx.Slider(panel1, -1, self.selpca, 1, 20, style=wx.SL_LEFT|wx.SL_VERTICAL)
        self.slidershow.Disable()          
        self.slidershow.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnPCAScroll, self.slidershow)
        

        self.pcaspin = wx.SpinButton(panel1, -1, size = ((8,-1)), style=wx.SP_ARROW_KEYS)
        self.Bind(wx.EVT_SPIN_UP, self.OnPCASpinUp, self.pcaspin)
        self.Bind(wx.EVT_SPIN_DOWN, self.OnPCASpinDown, self.pcaspin)
        
        hbox11.Add(i1panel, 0)
        vbox11.Add((0,3))
        vbox11.Add(self.slidershow, 1,  wx.EXPAND) 
        vbox11.Add(self.pcaspin, 0,  wx.EXPAND)         
        hbox11.Add(vbox11, 0,  wx.EXPAND) 

        
        vbox1.Add(self.tc_PCAcomp, 1, wx.EXPAND )        
        vbox1.Add(hbox11, 0)

        panel1.SetSizer(vbox1)
        
        
                
        #panel 2
        panel2 = wx.Panel(self, -1)
        vbox2 = wx.BoxSizer(wx.VERTICAL)
                
        sizer1 = wx.StaticBoxSizer(wx.StaticBox(panel2, -1, 'PCA'), wx.VERTICAL)
        self.button_calcpca = wx.Button(panel2, -1, 'Calculate PCA', (10, 10))
        self.button_calcpca.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnCalcPCA, id=self.button_calcpca.GetId())     
        self.button_calcpca.Disable()   
        sizer1.Add(self.button_calcpca, 0, wx.EXPAND)
        self.button_savepca = wx.Button(panel2, -1, 'Save PCA Results...', (10, 10))
        self.button_savepca.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnSave, id=self.button_savepca.GetId())
        self.button_savepca.Disable()
        sizer1.Add(self.button_savepca, 0, wx.EXPAND)
        
        hbox21 = wx.BoxSizer(wx.HORIZONTAL)
        text1 = wx.StaticText(panel2, -1, 'Number of significant components',  style=wx.ALIGN_LEFT)
        text1.SetFont(self.com.font)
        self.npcaspin = wx.SpinCtrl(panel2, -1, '',  size= (60, -1), style=wx.ALIGN_LEFT)
        self.npcaspin.SetRange(1,20)
        self.Bind(wx.EVT_SPINCTRL, self.OnNPCAspin, self.npcaspin)
        hbox21.Add(text1, 0, wx.TOP, 20)
        hbox21.Add((10,0))
        hbox21.Add(self.npcaspin, 0, wx.TOP, 15)            
        sizer1.Add(hbox21, 0, wx.EXPAND)    
              
        hbox22 = wx.BoxSizer(wx.HORIZONTAL)
        text2 = wx.StaticText(panel2, -1, 'Cumulative variance', style=wx.ALIGN_LEFT)
        text2.SetFont(self.com.font)
        self.vartc = wx.StaticText(panel2, -1, '0%',  style=wx.ALIGN_LEFT)
        hbox22.Add(text2, 0)
        hbox22.Add(self.vartc, 0, wx.LEFT , 10)

        sizer1.Add(hbox22, 0)      
        
        vbox2.Add(sizer1)

        panel2.SetSizer(vbox2)

        
        #panel 3
        panel3 = wx.Panel(self, -1)
        panel3.SetBackgroundColour("White")
  
        i3panel = wx.Panel(panel3, -1, style = wx.SUNKEN_BORDER)
        self.PCASpecPan = wxmpl.PlotPanel(i3panel, -1, size =(pw, ph), cursor=False, crosshairs=False, location=False, zoom=False)
             
        vbox3 = wx.BoxSizer(wx.VERTICAL)
        self.text_pcaspec = wx.TextCtrl(panel3, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE, size=(600,-1))
        self.text_pcaspec.SetFont(self.com.font)
        self.text_pcaspec.SetValue("PCA spectrum ")        
        vbox3.Add(self.text_pcaspec, 0)
        vbox3.Add(i3panel, 0)        
        panel3.SetSizer(vbox3)
        
        
        #panel 4
        panel4 = wx.Panel(self, -1)
        panel4.SetBackgroundColour("White")        

        i4panel = wx.Panel(panel4, -1, style = wx.SUNKEN_BORDER)
        self.PCAEvalsPan = wxmpl.PlotPanel(i4panel, -1, size =(pw, ph*0.75), cursor=False, crosshairs=False, location=False, zoom=False)
        wxmpl.EVT_POINT(i4panel, self.PCAEvalsPan.GetId(), self.OnPointEvalsImage)   
        
        vbox4 = wx.BoxSizer(wx.VERTICAL)
        text4 = wx.TextCtrl(panel4, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE, size=(600,-1))
        text4.SetFont(self.com.font)
        text4.SetValue("PCA eigenvalues ")        
        vbox4.Add(text4, 0)
        vbox4.Add(i4panel, 0)
        
        panel4.SetSizer(vbox4)      
        

        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        gridtop = wx.FlexGridSizer(2, 2, vgap=10, hgap=20)
        gridtop.Add(panel2, 0, wx.LEFT|wx.TOP, 20)
        gridtop.Add(panel4, 0, wx.ALIGN_LEFT)
        
        gridtop.Add(panel1, 0)
        gridtop.Add(panel3, 0, wx.ALIGN_RIGHT)
              
        vboxtop.Add((0,10))
        vboxtop.Add(gridtop, 0, wx.LEFT, 20)
         
        
        self.SetSizer(vboxtop) 
        
        
#----------------------------------------------------------------------
    def OnCalcPCA(self, event):
       
        wx.BeginBusyCursor()
        self.calcpca = False  
        self.selpca = 1       
        self.numsigpca = 2
        self.slidershow.SetValue(self.selpca)
        
        scrollmax = npy.min([self.stk.n_ev, 20])
        self.slidershow.SetMax(scrollmax)

        try: 
            self.CalcPCA()
            self.calcpca = True
            self.loadPCAImage()
            self.loadPCASpectrum()
            self.showEvals()
            self.com.pca_calculated = 1
            wx.EndBusyCursor() 
        except:
            self.com.pca_calculated = 0
            wx.EndBusyCursor()
            wx.MessageBox("PCA not calculated.")
        
        wx.GetApp().TopWindow.refresh_widgets()

#----------------------------------------------------------------------        
    def OnNPCAspin(self, event):
        num = event.GetInt()
        self.numsigpca = num
        
        if self.com.pca_calculated == 1:      
            self.anlz.numsigpca = self.numsigpca 
        
            # cumulative variance
            var = self.anlz.variance[:self.numsigpca].sum()
            self.vartc.SetLabel(str(var.round(decimals=2)*100)+'%')
                 
       
        
#----------------------------------------------------------------------
    def CalcPCA(self):
 
        self.anlz.calculate_pca()
     
        #Scree plot criterion
        self.numsigpca = self.anlz.numsigpca
        
        self.npcaspin.SetValue(self.numsigpca)
      
        # cumulative variance
        var = self.anlz.variance[:self.numsigpca].sum()
        self.vartc.SetLabel(str(var.round(decimals=2)*100)+'%')
        

        

 
#----------------------------------------------------------------------        
    def OnPCAScroll(self, event):
        self.sel = event.GetInt()
        self.selpca = self.sel
        if self.calcpca == True:
            self.loadPCAImage()
            self.loadPCASpectrum()


#----------------------------------------------------------------------            
    def OnPCASpinUp(self, event):
        if (self.calcpca == True) and (self.selpca > 1):
            self.selpca = self.selpca - 1
            self.slidershow.SetValue(self.selpca)

            self.loadPCAImage()
            self.loadPCASpectrum()

            
#----------------------------------------------------------------------            
    def OnPCASpinDown(self, event):
        if (self.calcpca == True) and (self.selpca < self.stk.n_ev-1):
            self.selpca = self.selpca + 1
            self.slidershow.SetValue(self.selpca) 
            
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

        
        SaveWinP2().Show()


            
#----------------------------------------------------------------------    
    def Save(self, filename, path, spec_png = True, spec_pdf = False, spec_csv = False,
             img_png = True, img_pdf = False, evals_png = True, evals_pdf = False): 
        
        
        self.SaveFileName = os.path.join(path,filename)
   
        try: 

            mtplot.rcParams['pdf.fonttype'] = 42
            if evals_png:
                ext = 'png'
                suffix = "." + ext
                fileName_evals = self.SaveFileName+"_PCAevals."+ext
                            
                fig = self.PCAEvalsPan.get_figure()
                fig.savefig(fileName_evals)
                
            if evals_pdf:
                ext = 'pdf'
                suffix = "." + ext
                fileName_evals = self.SaveFileName+"_PCAevals."+ext
                            
                fig = self.PCAEvalsPan.get_figure()
                fig.savefig(fileName_evals)                

                        
            from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
            mtplot.rcParams['pdf.fonttype'] = 42
            
            ext = 'png'
            suffix = "." + ext        
                
            if img_png:
                for i in range (self.numsigpca):
              
                    self.pcaimage = self.anlz.pcaimages[:,:,i]
              
                    fig = mtplot.figure.Figure(figsize =(PlotH*1.15, PlotH))
                    canvas = FigureCanvas(fig)
                    axes = fig.gca()
                    divider = make_axes_locatable(axes)
                    ax_cb = divider.new_horizontal(size="3%", pad=0.03)  
                    fig.add_axes(ax_cb)
                    axes.set_position([0.03,0.03,0.8,0.94])
                    bound = self.anlz.pcaimagebounds[i]       
                    mtplot.rcParams['font.size'] = self.fontsize        
        
                    im = axes.imshow(self.pcaimage, cmap=mtplot.cm.get_cmap("seismic_r"), vmin = -bound, vmax = bound)
                    cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)  
                    axes.axis("off") 
                                
                    fileName_img = self.SaveFileName+"_PCA_" +str(i+1)+"."+ext
                    fig.savefig(fileName_img, bbox_inches='tight', pad_inches = 0.0)
            
            if spec_png:
                for i in range (self.numsigpca):
                
                    pcaspectrum = self.anlz.eigenvecs[:,i]
                    fig = mtplot.figure.Figure(figsize =(PlotW, PlotH))
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.15,0.15,0.75,0.75))
                    axes = fig.gca()
                    mtplot.rcParams['font.size'] = self.fontsize
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
              
                    fig = mtplot.figure.Figure(figsize =(PlotH*1.15, PlotH))
                    canvas = FigureCanvas(fig)
                    axes = fig.gca()
                    divider = make_axes_locatable(axes)
                    ax_cb = divider.new_horizontal(size="3%", pad=0.03)  
                    fig.add_axes(ax_cb)
                    axes.set_position([0.03,0.03,0.8,0.94])
                    bound = self.anlz.pcaimagebounds[i]       
                    mtplot.rcParams['font.size'] = self.fontsize        
        
                    im = axes.imshow(self.pcaimage, cmap=mtplot.cm.get_cmap("seismic_r"), vmin = -bound, vmax = bound)
                    cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)  
                    axes.axis("off") 
                                
                    fileName_img = self.SaveFileName+"_PCA_" +str(i+1)+"."+ext
                    fig.savefig(fileName_img, bbox_inches='tight', pad_inches = 0.0)
            
            if spec_pdf:
                for i in range (self.numsigpca):
                
                    self.pcaspectrum = self.anlz.eigenvecs[:,i]
                    fig = mtplot.figure.Figure(figsize =(PlotW, PlotH))
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.15,0.15,0.75,0.75))
                    axes = fig.gca()
                    mtplot.rcParams['font.size'] = self.fontsize
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
   
            wx.MessageBox('Could not save file: %s' % err, 'Error', 
                          parent=self, style=wx.OK|wx.ICON_ERROR) 
            
        
     
#----------------------------------------------------------------------      
    def showEvals(self):
        
        evalmax = npy.min([self.stk.n_ev, 40])
        self.pcaevals = self.anlz.eigenvals[0:evalmax]
        

        fig = self.PCAEvalsPan.get_figure()
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        
        mtplot.rcParams['font.size'] = self.fontsize
       
        evalsplot = axes.semilogy(npy.arange(1,evalmax+1), self.pcaevals,'b.')    
        
        axes.set_xlabel('Principal Component')
        axes.set_ylabel('Log(Eigenvalue)')
         

        self.PCAEvalsPan.draw()


#----------------------------------------------------------------------      
    def loadPCAImage(self):
        
        self.tc_PCAcomp.SetValue("PCA component " + str(self.selpca))
        self.text_pcaspec.SetValue("PCA spectrum "+ str(self.selpca))          
        
        self.pcaimage = self.anlz.pcaimages[:,:,self.selpca-1]
        
        
        fig = self.PCAImagePan.get_figure()
        fig.clf()
     
        axes = fig.gca()
    
        divider = make_axes_locatable(axes)
        ax_cb = divider.new_horizontal(size="3%", pad=0.03)  

        fig.add_axes(ax_cb)
        
        axes.set_position([0.03,0.03,0.8,0.94])
        
        bound = self.anlz.pcaimagebounds[self.selpca-1]
        
     
        im = axes.imshow(self.pcaimage, cmap=mtplot.cm.get_cmap("seismic_r"), vmin = -bound, vmax = bound)
        cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)  
    
        axes.axis("off") 
        self.PCAImagePan.draw()

        
        
#----------------------------------------------------------------------     
    def loadPCASpectrum(self):

        self.pcaspectrum = self.anlz.eigenvecs[:,self.selpca-1]
            
        
        fig = self.PCASpecPan.get_figure()
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        
        mtplot.rcParams['font.size'] = self.fontsize

        specplot = axes.plot(self.stk.ev,self.pcaspectrum)
                        
        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Optical Density')
        
        self.PCASpecPan.draw()
        
        
        
#---------------------------------------------------------------------- 
class SaveWinP2(wx.Frame):
    
    title = "Save"

    def __init__(self):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=self.title, size=(400, 330))
               
        ico = logos.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White") 

        
        self.com = wx.GetApp().TopWindow.common         
        self.fontsize = self.com.fontsize   
        
        
        path, ext = os.path.splitext(self.com.filename) 
        ext = ext[1:].lower()   
        suffix = "." + ext
        path, fn = os.path.split(self.com.filename)
        filename = fn[:-len(suffix)]
        
        self.path = self.com.path
        self.filename = filename
            
        
        
        vboxtop = wx.BoxSizer(wx.VERTICAL)
                
        panel1 = wx.Panel(self, -1)
        
        gridtop = wx.FlexGridSizer(4, 4, vgap=20, hgap=20)
    
        fontb = wx.SystemSettings_GetFont(wx.SYS_DEFAULT_GUI_FONT)
        fontb.SetWeight(wx.BOLD)
        
        
        st1 = wx.StaticText(panel1, -1, 'Save',  style=wx.ALIGN_LEFT)
        st1.SetFont(fontb)
        st2 = wx.StaticText(panel1, -1, '.pdf',  style=wx.ALIGN_LEFT)
        st2.SetFont(fontb)
        st3 = wx.StaticText(panel1, -1, '.png',  style=wx.ALIGN_LEFT)
        st3.SetFont(fontb)
        st3a = wx.StaticText(panel1, -1, '.csv',  style=wx.ALIGN_LEFT)
        st3a.SetFont(fontb)        
        
        st4 = wx.StaticText(panel1, -1, '_spectra',  style=wx.ALIGN_LEFT)
        st4.SetFont(self.com.font)
        
        self.cb1 = wx.CheckBox(panel1, -1, '')
        self.cb1.SetFont(self.com.font)
        self.cb1.Set3StateValue(wx.CHK_CHECKED)
        
        self.cb2 = wx.CheckBox(panel1, -1, '')
        self.cb2.SetFont(self.com.font)
        self.cb2a = wx.CheckBox(panel1, -1, '')
        self.cb2a.SetFont(self.com.font)
        
        st5 = wx.StaticText(panel1, -1, '_images',  style=wx.ALIGN_LEFT)
        st5.SetFont(self.com.font)
        
        self.cb3 = wx.CheckBox(panel1, -1, '')
        self.cb3.SetFont(self.com.font)
        self.cb3.Set3StateValue(wx.CHK_CHECKED)
        
        self.cb4 = wx.CheckBox(panel1, -1, '')
        self.cb4.SetFont(self.com.font)
        
        
        st6 = wx.StaticText(panel1, -1, '_eigenvals',  style=wx.ALIGN_LEFT)
        st6.SetFont(self.com.font)
        
        self.cb5 = wx.CheckBox(panel1, -1, '')
        self.cb5.SetFont(self.com.font)
        self.cb5.Set3StateValue(wx.CHK_CHECKED)
        
        self.cb6 = wx.CheckBox(panel1, -1, '')
        self.cb6.SetFont(self.com.font)

        
        gridtop.Add(st1, 0)
        gridtop.Add(st2, 0)
        gridtop.Add(st3, 0)
        gridtop.Add(st3a, 0)
                
        gridtop.Add(st4, 0)
        gridtop.Add(self.cb1, 0)
        gridtop.Add(self.cb2, 0)              
        gridtop.Add(self.cb2a, 0)  
          
        gridtop.Add(st5, 0)
        gridtop.Add(self.cb3, 0)
        gridtop.Add(self.cb4, 0)  
        gridtop.Add(wx.StaticText(panel1, -1, ' '))
        
        
        gridtop.Add(st6, 0)
        gridtop.Add(self.cb5, 0)
        gridtop.Add(self.cb6, 0)  
        gridtop.Add(wx.StaticText(panel1, -1, ' '))
        
        panel1.SetSizer(gridtop)
        
        panel2 = wx.Panel(self, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        
        hbox0 = wx.BoxSizer(wx.HORIZONTAL)
        
        st0 = wx.StaticText(panel2, -1, 'Filename: ')
        st0.SetFont(self.com.font)
        self.tc_savefn = wx.TextCtrl(panel2, -1,  style=wx.TE_RICH, size =((100, -1)),
                                         value=self.filename)
        self.tc_savefn.SetFont(self.com.font)
        hbox0.Add(st0, 0, flag = wx.ALIGN_CENTER_VERTICAL)
        hbox0.Add(self.tc_savefn,1, wx.EXPAND|wx.LEFT, 2)         
        
        
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
                
        st1 = wx.StaticText(panel2, label='Path: ')
        st1.SetFont(self.com.font)
        self.tc_savepath = wx.TextCtrl(panel2, -1,  style=wx.TE_RICH|wx.TE_READONLY, size =((100, -1)),
                                         value=self.path)
        self.tc_savepath.SetFont(self.com.font)
        hbox1.Add(st1, 0, wx.ALIGN_CENTER_VERTICAL|wx.RIGHT, 20)
        hbox1.Add(self.tc_savepath, 1, wx.EXPAND|wx.LEFT, 2)  
        
        button_path = wx.Button(panel2, -1, 'Browse...')
        self.Bind(wx.EVT_BUTTON, self.OnBrowseDir, id=button_path.GetId())
        hbox1.Add(button_path, 0)
        
        
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        button_save = wx.Button(panel2, -1, 'Save')
        self.Bind(wx.EVT_BUTTON, self.OnSave, id=button_save.GetId())
        hbox2.Add(button_save, 0, wx.ALIGN_RIGHT)
        
        button_cancel = wx.Button(panel2, -1, 'Cancel')
        self.Bind(wx.EVT_BUTTON, self.OnCancel, id=button_cancel.GetId())
        hbox2.Add(button_cancel, 0, wx.ALIGN_RIGHT|wx.LEFT,10)
        
        vbox1.Add(hbox0, 1, wx.EXPAND)
        vbox1.Add((0,5))
        vbox1.Add(hbox1, 1, wx.EXPAND)
        vbox1.Add((0,20))
        vbox1.Add(hbox2, 0, wx.ALIGN_RIGHT)
        panel2.SetSizer(vbox1)
        
        
        vboxtop.Add(panel1, 1, wx.ALL | wx.EXPAND, 20)
        vboxtop.Add(panel2, 0, wx.ALL | wx.EXPAND, 20) 
        
        self.SetSizer(vboxtop)    
        
#----------------------------------------------------------------------        
    def OnBrowseDir(self, evt):
        
        dialog = wx.DirDialog(None, "Choose a directory",
                               style=wx.DD_DIR_MUST_EXIST,
                               defaultPath=self.path)

        if dialog.ShowModal() == wx.ID_OK:
            directory = dialog.GetPath()
            
            self.path = directory
            
            self.tc_savepath.SetValue(self.path)
            
            
                
#----------------------------------------------------------------------        
    def OnSave(self, evt):
        
        self.filename = self.tc_savefn.GetValue()
        
        sp_pdf = self.cb1.GetValue()
        sp_png = self.cb2.GetValue()
        sp_csv = self.cb2a.GetValue()
        im_pdf = self.cb3.GetValue()
        im_png = self.cb4.GetValue()
        ev_pdf = self.cb5.GetValue()
        ev_png = self.cb6.GetValue()
        
        self.Destroy() 
        wx.GetApp().TopWindow.page2.Save(self.filename, self.path,
                                         spec_png = sp_png, 
                                         spec_pdf = sp_pdf, 
                                         spec_csv = sp_csv,
                                         img_png = im_png, 
                                         img_pdf = im_pdf,
                                         evals_png = ev_png, 
                                         evals_pdf = ev_pdf)

#---------------------------------------------------------------------- 
    def OnCancel(self, evt):
        self.Destroy()  
                
                
                
 

""" ------------------------------------------------------------------------------------------------"""
class PageStack(wx.Panel):
    def __init__(self, parent, common, data_struct, stack):
        wx.Panel.__init__(self, parent)
        
        self.data_struct = data_struct
        self.stk = stack
        self.com = common                  
        self.SetBackgroundColour("white")
        
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

        vbox = wx.BoxSizer(wx.VERTICAL)
        hboxT = wx.BoxSizer(wx.HORIZONTAL)
        hboxB = wx.BoxSizer(wx.HORIZONTAL)
    
        #panel 1        
        panel1 = wx.Panel(self, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        panel1.SetBackgroundColour("White")
        self.tc_imageeng = wx.TextCtrl(panel1, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_imageeng.SetFont(self.com.font)
        self.tc_imageeng.SetValue("Image at energy: ")

       
        hbox11 = wx.BoxSizer(wx.HORIZONTAL)
   
        i1panel = wx.Panel(panel1, -1, style = wx.SUNKEN_BORDER)
        self.AbsImagePanel = wxmpl.PlotPanel(i1panel, -1, size =(PlotH, PlotH), cursor=False, crosshairs=True, location=False, zoom=False)
        wxmpl.EVT_POINT(i1panel, self.AbsImagePanel.GetId(), self.OnPointAbsimage)
        
        vbox11 = wx.BoxSizer(wx.VERTICAL)                    
        self.slider_eng = wx.Slider(panel1, -1, self.iev, 0, 100, style=wx.SL_LEFT|wx.SL_VERTICAL)        
        self.slider_eng.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnScrollEng, self.slider_eng)
        
        
        self.engspin = wx.SpinButton(panel1, -1, size = ((8,-1)), style=wx.SP_ARROW_KEYS)
        self.Bind(wx.EVT_SPIN_UP, self.OnEngspinUp, self.engspin)
        self.Bind(wx.EVT_SPIN_DOWN, self.OnEngspinDown, self.engspin)
       

        vbox11.Add((0,3))
        vbox11.Add(self.slider_eng, 1,  wx.EXPAND) 
        vbox11.Add(self.engspin, 0,  wx.EXPAND)    
        hbox11.Add(i1panel, 0)
        hbox11.Add(vbox11, 0,  wx.EXPAND)              
        
        vbox1.Add(self.tc_imageeng,1, wx.LEFT | wx.TOP | wx.EXPAND, 20)        
        vbox1.Add(hbox11, 0,  wx.LEFT, 20)

        panel1.SetSizer(vbox1)
     
     
        #panel 2
        panel2 = wx.Panel(self, -1)
        vbox2 = wx.BoxSizer(wx.VERTICAL)
        panel2.SetBackgroundColour("White")
        self.tc_spec = wx.TextCtrl(panel2, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_spec.SetValue("Spectrum")
        self.tc_spec.SetFont(self.com.font)
          
        i2panel = wx.Panel(panel2, -1, style = wx.SUNKEN_BORDER)
        self.SpectrumPanel = wxmpl.PlotPanel(i2panel, -1, size=(PlotW, PlotH), cursor=False, crosshairs=False, location=False, zoom=False)
        wxmpl.EVT_POINT(i2panel, self.SpectrumPanel.GetId(), self.OnPointSpectrum)
        
        vbox2.Add(self.tc_spec, 1, wx.LEFT | wx.TOP | wx.EXPAND, 20)       
        vbox2.Add(i2panel, 0, wx.LEFT , 20)
        
        panel2.SetSizer(vbox2)
        
        
        #panel 3
        panel3 = wx.Panel(self, -1)
        sizer1 = wx.StaticBoxSizer(wx.StaticBox(panel3, -1, 'Preprocess'),orient=wx.VERTICAL)
        vbox31 = wx.BoxSizer(wx.VERTICAL)
        vbox31.Add((0,3))        
        
        self.button_align = wx.Button(panel3, -1, 'Align images...')
        self.Bind(wx.EVT_BUTTON, self.OnAlignImgs, id=self.button_align.GetId())
        self.button_align.SetFont(self.com.font)
        self.button_align.Disable()
        vbox31.Add(self.button_align, 0, wx.EXPAND) 
        
        vbox31.Add((0,2))       

        self.button_i0ffile = wx.Button(panel3, -1, 'I0 from file...')
        self.button_i0ffile.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnI0FFile, id=self.button_i0ffile.GetId())
        self.button_i0ffile.Disable()
        vbox31.Add(self.button_i0ffile, 0, wx.EXPAND)
        self.button_i0histogram = wx.Button(panel3, -1, 'I0 from histogram...')
        self.button_i0histogram.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnI0histogram, id=self.button_i0histogram.GetId())   
        self.button_i0histogram.Disable()     
        vbox31.Add(self.button_i0histogram, 0, wx.EXPAND)
        self.button_showi0 = wx.Button(panel3, -1, 'Show I0...')
        self.button_showi0.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnShowI0, id=self.button_showi0.GetId())   
        self.button_showi0.Disable()
        vbox31.Add(self.button_showi0, 0, wx.EXPAND)
        
        vbox31.Add((0,2))
        
        self.button_limitev = wx.Button(panel3, -1, 'Limit energy range...')
        self.Bind(wx.EVT_BUTTON, self.OnLimitEv, id=self.button_limitev.GetId())
        self.button_limitev.SetFont(self.com.font)
        self.button_limitev.Disable()
        vbox31.Add(self.button_limitev, 0, wx.EXPAND)
        
        self.button_subregion = wx.Button(panel3, -1, 'Clip to subregion...')
        self.Bind(wx.EVT_BUTTON, self.OnCliptoSubregion, id=self.button_subregion.GetId())
        self.button_subregion.SetFont(self.com.font)
        self.button_subregion.Disable()
        vbox31.Add(self.button_subregion, 0, wx.EXPAND)
        
        vbox31.Add((0,2))
        
            
        self.button_savestack = wx.Button(panel3, -1, 'Save preprocessed stack')
        self.button_savestack.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnSaveStack, id=self.button_savestack.GetId())
        self.button_savestack.Disable()          
        vbox31.Add(self.button_savestack, 0, wx.EXPAND)
        sizer1.Add(vbox31,1, wx.LEFT|wx.RIGHT|wx.EXPAND,2)
        panel3.SetSizer(sizer1)
        

        
        #panel 4
        panel4 = wx.Panel(self, -1)
        vbox4 = wx.BoxSizer(wx.VERTICAL)
        sb = wx.StaticBox(panel4, -1, 'Display')
        sb.SetBackgroundColour("White")
        sizer4 = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)
        sb = wx.StaticBox(panel4, -1, 'File')
        sb.SetBackgroundColour("White")
        sizer41 = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)

        hbox40 = wx.BoxSizer(wx.HORIZONTAL)    
        self.textctrl = wx.TextCtrl(panel4, -1, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        sizer41.Add(self.textctrl, 0, wx.EXPAND|wx.TOP|wx.LEFT, 5)
        hbox40.Add(sizer41, 1, wx.EXPAND)
        
        vbox41 = wx.BoxSizer(wx.VERTICAL)
        vbox41.Add((0,2))
        self.button_slideshow = wx.Button(panel4, -1, 'Play stack movie')
        self.Bind(wx.EVT_BUTTON, self.OnSlideshow, id=self.button_slideshow.GetId())
        self.button_slideshow.SetFont(self.com.font)
        self.button_slideshow.Disable()
        vbox41.Add(self.button_slideshow, 0, wx.EXPAND)
        
        self.button_save = wx.Button(panel4, -1, 'Save images...')
        self.button_save.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnSave, id=self.button_save.GetId())
        self.button_save.Disable()          
        vbox41.Add(self.button_save, 0, wx.EXPAND)
              
        hbox40.Add(vbox41, 0, wx.EXPAND|wx.LEFT|wx.RIGHT, 20)
        vbox4.Add(hbox40, 0, wx.EXPAND)
        
    
        hbox41 = wx.BoxSizer(wx.HORIZONTAL)
        sizer42 = wx.StaticBoxSizer(wx.StaticBox(panel4, -1, 'Image'),  orient=wx.VERTICAL)
        self.rb_flux = wx.RadioButton(panel4, -1, 'Flux', style=wx.RB_GROUP)
        self.rb_od = wx.RadioButton(panel4, -1, 'Optical Density')
        self.rb_flux.SetFont(self.com.font)
        self.rb_od.SetFont(self.com.font)
        self.Bind(wx.EVT_RADIOBUTTON, self.onrb_fluxod, id=self.rb_flux.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.onrb_fluxod, id=self.rb_od.GetId())
        
        self.rb_flux.Disable()
        self.rb_od.Disable()

        self.add_scale_cb = wx.CheckBox(panel4, -1, '  Scalebar')
        self.add_scale_cb.SetFont(self.com.font)
        self.Bind(wx.EVT_CHECKBOX, self.OnShowScale, self.add_scale_cb)
        self.add_scale_cb.SetValue(True)
        
        self.add_colbar_cb = wx.CheckBox(panel4, -1, '  Colorbar')
        self.add_colbar_cb.SetFont(self.com.font)
        self.Bind(wx.EVT_CHECKBOX, self.OnShowColBar, self.add_colbar_cb)

        sizer42.Add((0,10))
        sizer42.Add(self.rb_flux, 0, wx.LEFT|wx.RIGHT, 2)
        sizer42.Add((0,3))
        sizer42.Add(self.rb_od, 0, wx.LEFT|wx.RIGHT, 2)
        sizer42.Add((0,10))
        sizer42.Add(self.add_scale_cb, 0, wx.LEFT|wx.RIGHT, 2)
        sizer42.Add((0,3))
        sizer42.Add(self.add_colbar_cb, 0, wx.LEFT|wx.RIGHT, 2)
        
        hbox41.Add(sizer42, 0, wx.EXPAND)
                

        sizer43 = wx.StaticBoxSizer(wx.StaticBox(panel4, -1, 'Display settings'),  orient=wx.VERTICAL)
        sizer43.Add((0,10))
        hbox42 = wx.BoxSizer(wx.HORIZONTAL)
        hbox42.Add((10,0))

        fgs41 = wx.FlexGridSizer(3, 2, 0, 0)
        
        self.tc_min = wx.TextCtrl(panel4, -1, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE, size=(100,-1))
        self.tc_min.SetFont(self.com.font)
        self.tc_min.AppendText('Minimum: \t{0:5d}%'.format(int(100*self.brightness_min)))
        
        self.tc_max = wx.TextCtrl(panel4, -1, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE, size=(100,-1))
        self.tc_max.SetFont(self.com.font)
        self.tc_max.AppendText('Maximum:{0:5d}%'.format(int(100*self.brightness_max)))

        self.slider_brightness_min = wx.Slider(panel4, -1, self.dispbrightness_min, 0, 49, size= (50,20), style=wx.SL_HORIZONTAL)        
        self.slider_brightness_min.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnScrollBrightnessMin, self.slider_brightness_min)
        
        self.slider_brightness_max = wx.Slider(panel4, -1, self.dispbrightness_max, 50, 120, style=wx.SL_HORIZONTAL)        
        self.slider_brightness_max.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnScrollBrightnessMax, self.slider_brightness_max)        
        
        self.tc_gamma = wx.TextCtrl(panel4, -1, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE, size=(100,-1))
        self.tc_gamma.SetFont(self.com.font)
        self.tc_gamma.AppendText('Gamma:  \t{0:5.2f}'.format(self.gamma))
        
        self.slider_gamma = wx.Slider(panel4, -1, self.displaygamma, 1, 20, style=wx.SL_HORIZONTAL)        
        self.slider_gamma.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnScrollGamma, self.slider_gamma)

        
        fgs41.AddMany([(self.tc_min), (self.slider_brightness_min, 0, wx.EXPAND), (self.tc_max), 
            (self.slider_brightness_max, 0, wx.EXPAND),(self.tc_gamma), (self.slider_gamma, 0, wx.EXPAND)])
        
      
        hbox42.Add(fgs41, 0, wx.EXPAND)
        hbox42.Add((20,0))

        
        vbox43 = wx.BoxSizer(wx.VERTICAL)
        self.button_despike = wx.Button(panel4, -1, 'Despike', (100,10))
        self.button_despike.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnDespike, id=self.button_despike.GetId())   
        self.button_despike.Disable()     
        vbox43.Add(self.button_despike, 0, wx.EXPAND)
        self.button_resetdisplay = wx.Button(panel4, -1, 'Reset')
        self.button_resetdisplay.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.onResetDisplaySettings, id=self.button_resetdisplay.GetId())   
        self.button_resetdisplay.Disable()     
        vbox43.Add(self.button_resetdisplay, 0, wx.EXPAND)
        self.button_displaycolor = wx.Button(panel4, -1, '   Color Table...   ')
        self.button_displaycolor.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.onSetColorTable, id=self.button_displaycolor.GetId())   
        self.button_displaycolor.Disable()     
        vbox43.Add(self.button_displaycolor, 0, wx.EXPAND)
        
        hbox42.Add(vbox43, 0, wx.EXPAND|wx.LEFT|wx.RIGHT, 3)                   
        sizer43.Add(hbox42)
    
        hbox41.Add((2,0))
        hbox41.Add(sizer43, 1, wx.EXPAND)
        vbox4.Add((0,10))
        vbox4.Add(hbox41, 1, wx.EXPAND)

        sizer4.Add(vbox4,1, wx.EXPAND)
        
        panel4.SetSizer(sizer4)
        
        
        #panel 5
        panel5 = wx.Panel(self, -1)
        sizer5 = wx.StaticBoxSizer(wx.StaticBox(panel5, -1, 'Region of Interest'), wx.VERTICAL)
        
        vbox51 = wx.BoxSizer(wx.VERTICAL)
        #vbox51.Add((0,1))
        self.button_addROI = wx.Button(panel5, -1, 'Add ROI')
        self.button_addROI.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnAddROI, id=self.button_addROI.GetId())
        self.button_addROI.Disable()
        vbox51.Add(self.button_addROI, 0, wx.EXPAND)
        
        self.button_acceptROI = wx.Button(panel5, -1, 'Accept ROI')
        self.button_acceptROI.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnAcceptROI, id=self.button_acceptROI.GetId())   
        self.button_acceptROI.Disable()     
        vbox51.Add(self.button_acceptROI, 0, wx.EXPAND)
        
        self.button_resetROI = wx.Button(panel5, -1, 'Reset ROI')
        self.button_resetROI.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnResetROI, id=self.button_resetROI.GetId())
        self.button_resetROI.Disable()
        vbox51.Add(self.button_resetROI, 0, wx.EXPAND) 

        self.button_setROII0 = wx.Button(panel5, -1, 'Set ROI As I0')
        self.button_setROII0.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnSetROII0, id=self.button_setROII0.GetId())
        self.button_setROII0.Disable()
        vbox51.Add(self.button_setROII0, 0, wx.EXPAND)
        
        self.button_saveROIspectr = wx.Button(panel5, -1, 'Save ROI Spectrum...')
        self.button_saveROIspectr.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnSaveROISpectrum, id=self.button_saveROIspectr.GetId())   
        self.button_saveROIspectr.Disable()     
        vbox51.Add(self.button_saveROIspectr, 0, wx.EXPAND)
        
        vbox51.Add((0,1))
        
        self.button_ROIdosecalc = wx.Button(panel5, -1, 'ROI Dose Calculation...')
        self.button_ROIdosecalc.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnROI_DoseCalc, id=self.button_ROIdosecalc.GetId())   
        self.button_ROIdosecalc.Disable()     
        vbox51.Add(self.button_ROIdosecalc, 0, wx.EXPAND)       
        vbox51.Add((0,1))        
        
        self.button_spectralROI = wx.Button(panel5, -1, 'Spectral ROI...')
        self.button_spectralROI.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnSpectralROI, id=self.button_spectralROI.GetId())   
        self.button_spectralROI.Disable()     
        vbox51.Add(self.button_spectralROI, 0, wx.EXPAND)
        
        
        
        sizer5.Add(vbox51,1, wx.ALL|wx.EXPAND,2)        
        panel5.SetSizer(sizer5)

        

        hboxB.Add(panel1, 0, wx.BOTTOM | wx.TOP, 9)
        hboxB.Add(panel2, 0, wx.BOTTOM | wx.TOP, 9)
        hboxT.Add((10,0)) 
               
        hboxT.Add(panel3, 1, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 10)
        hboxT.Add(panel4, 0, wx.LEFT | wx.RIGHT |wx.TOP | wx.EXPAND,10)
        hboxT.Add(panel5, 1, wx.LEFT | wx.RIGHT |wx.TOP | wx.EXPAND,10)

        vbox.Add(hboxT, 0, wx.ALL, 5)
        
        vbox.Add(hboxB, 0,  wx.ALL, 5)
  
        self.SetSizer(vbox) 
        
        

      

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
                

                      
        fig = self.AbsImagePanel.get_figure()
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
            im = axes.imshow(image, cmap=mtplot.cm.get_cmap(self.colortable)) 
        else:
            imgmax = npy.amax(image)
            imgmin = npy.amin(image)
            if (self.gamma != 1.0) or (imgmin < 0.0):
                image = (image-imgmin)/(imgmax-imgmin)
                imgmax = 1.0
                imgmin = 0.0
                if (self.gamma != 1.0):
                    image = npy.power(image, self.gamma)
            im = axes.imshow(image, cmap=mtplot.cm.get_cmap(self.colortable), 
                             vmin=(imgmin+imgmax*self.brightness_min),vmax=imgmax*self.brightness_max)
            
        if (self.showROImask == 1) and (self.addroi == 1):
            im_red = axes.imshow(self.ROIpix_masked,cmap=mtplot.cm.get_cmap("autumn")) 
         
        axes.axis("off") 
        
        if self.show_colorbar == 1:
            cbar = axes.figure.colorbar(im, orientation='vertical',cax=axcb) 
        
        if self.show_scale_bar == 1:
            startx = int(self.stk.n_rows*0.05)
            starty = self.stk.n_cols-int(self.stk.n_cols*0.05)-self.stk.scale_bar_pixels_y
            um_string = ' $\mathrm{\mu m}$'
            microns = '$'+self.stk.scale_bar_string+' $'+um_string
            axes.text(self.stk.scale_bar_pixels_x+startx+1,starty+1, microns, horizontalalignment='left', verticalalignment='center',
                      color = 'black', fontsize=14)
            #Matplotlib has flipped scales so I'm using rows instead of cols!
            p = mtplot.patches.Rectangle((startx,starty), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                                   color = 'black', fill = True)
            axes.add_patch(p)
        
        self.AbsImagePanel.draw()
        
        self.tc_imageeng.SetValue("Image at energy: {0:5.2f} eV".format(float(self.stk.ev[self.iev])))


#----------------------------------------------------------------------          
    def loadSpectrum(self, xpos, ypos):
        
        fig = self.SpectrumPanel.get_figure()
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

        
        mtplot.rcParams['font.size'] = self.fontsize

        specplot = axes.plot(self.stk.ev,self.spectrum)
        

        
        axes.axvline(x=self.stk.ev[self.iev], color = 'g', alpha=0.5)

        
        self.SpectrumPanel.draw()
        
        self.tc_spec.SetValue("Spectrum at pixel [" +str(ypos)+", " + str(xpos)+"] or position ["+
                              str(self.stk.x_dist[xpos])+", "+ str(self.stk.y_dist[ypos])+ "]")

        
#----------------------------------------------------------------------            
    def OnScrollEng(self, event):
        self.iev = event.GetInt()
        if self.com.stack_loaded == 1:
            self.loadImage()
            self.loadSpectrum(self.ix, self.iy)

#----------------------------------------------------------------------            
    def OnEngspinUp(self, event):
        if (self.com.stack_loaded == 1) and (self.iev > 0):
            self.iev = self.iev - 1
            self.slider_eng.SetValue(self.iev)

            self.loadImage()
            self.loadSpectrum(self.ix, self.iy)
            
#----------------------------------------------------------------------            
    def OnEngspinDown(self, event):
        if (self.com.stack_loaded == 1) and (self.iev < self.stk.n_ev-1):
            self.iev = self.iev + 1
            self.slider_eng.SetValue(self.iev) 
            
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
            
            self.slider_eng.SetValue(self.iev)
            
                   
#----------------------------------------------------------------------  
    def OnPointAbsimage(self, evt):
        x = evt.xdata
        y = evt.ydata
        
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
                self.line = mtplot.lines.Line2D([x,  x], [y, y], marker = '.', color = 'red')
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
        
    def OnI0FFile(self, event):


        try: 
            wildcard = "I0 CSV files (*.csv)|*.csv| I0 files (*.xas)|*.xas|SDF I0 files (*.hdr)|*.hdr"
            dialog = wx.FileDialog(None, "Choose i0 file",
                                   wildcard=wildcard,
                                   defaultDir = defaultDir,
                                   style=wx.OPEN)
            if dialog.ShowModal() == wx.ID_OK:
                filepath_i0 = dialog.GetPath()
                self.filename = dialog.GetFilename()
                                                        
            basename, extension = os.path.splitext(self.filename)      
            
            
            if extension == '.hdr':
                wx.BeginBusyCursor()                                    

                x=self.stk.n_cols
                y=self.stk.n_rows
                z=self.iev               

                self.ix = x/2
                self.iy = y/2
                self.stk.read_sdf_i0(filepath_i0)
                self.loadSpectrum(self.ix, self.iy)
                self.loadImage()
                self.com.i0_loaded = 1
                wx.EndBusyCursor()
                
                
            elif extension == '.xas':
                wx.BeginBusyCursor()                                    

                x=self.stk.n_cols
                y=self.stk.n_rows
                z=self.iev               

                self.ix = x/2
                self.iy = y/2

                self.stk.read_stk_i0(filepath_i0, extension)

                self.loadSpectrum(self.ix, self.iy)
                self.loadImage()
                self.com.i0_loaded = 1
                wx.EndBusyCursor()

            elif extension == '.csv':
                wx.BeginBusyCursor()                             

                x=self.stk.n_cols
                y=self.stk.n_rows
                z=self.iev               

                self.ix = x/2
                self.iy = y/2

                self.stk.read_stk_i0(filepath_i0, extension)

                self.loadSpectrum(self.ix, self.iy)
                self.loadImage()
                self.com.i0_loaded = 1
                wx.EndBusyCursor()
                
            
        except:
            wx.BeginBusyCursor()  
            self.com.i0_loaded = 0        
            wx.EndBusyCursor()  
            wx.MessageBox("I0 file not loaded.")
            import sys; print sys.exc_info()
                       
        dialog.Destroy()
                          
        wx.GetApp().TopWindow.refresh_widgets()

        
        
        
#----------------------------------------------------------------------       
    def OnI0histogram(self, event):    
        wx.GetApp().TopWindow.Hide()
        ShowHistogram(self.stk).Show()
         

#----------------------------------------------------------------------       
    def I0histogramCalculated(self):    
        
        PlotFrame(self.stk.evi0hist,self.stk.i0datahist).Show()
        
         
        self.loadSpectrum(self.ix, self.iy)
        self.loadImage()
        
        self.com.i0_loaded = 1
        wx.GetApp().TopWindow.refresh_widgets()
        
#----------------------------------------------------------------------    
    def OnShowI0(self, event):     

        PlotFrame(self.stk.evi0,self.stk.i0data).Show()
 
#----------------------------------------------------------------------    
    def OnSaveStack(self, event): 
        
        wx.GetApp().TopWindow.SaveStackH5()
               
#----------------------------------------------------------------------    
    def OnSave(self, event):     
        
        
        SaveWinP1().Show()

 
   
#----------------------------------------------------------------------    
    def Save(self, filename, path, spec_png = True, spec_pdf = False, sp_csv = False, img_png = True, img_pdf = False, img_all = False): 

        self.SaveFileName = os.path.join(path,filename)
      
        
        try: 
            ext = 'png'
            suffix = "." + ext
            mtplot.rcParams['pdf.fonttype'] = 42
            
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
                wx.BeginBusyCursor()
                from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
                mtplot.rcParams['pdf.fonttype'] = 42

                for i in range (self.stk.n_ev):
                    if self.showflux:
                        #Show flux image      
                        image = self.stk.absdata[:,:,i] 
                    else:
                        #Show OD image
                        image = self.stk.od3d[:,:,i]

                    fig = mtplot.figure.Figure(figsize =(float(self.stk.n_rows)/10, float(self.stk.n_cols)/10))
                    fig.clf()
                    canvas = FigureCanvas(fig)
                    fig.add_axes((0.0,0.0,1.0,1.0))
                    axes = fig.gca()
                    im = axes.imshow(image, cmap=mtplot.cm.get_cmap(self.colortable)) 
                    axes.axis("off") 
                                
                    fileName_img = self.SaveFileName+"_imnum_" +str(i+1)+"."+ext
                    fig.savefig(fileName_img,  dpi=ImgDpi, pad_inches = 0.0)
                wx.EndBusyCursor()
                    
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
   
            wx.MessageBox('Could not save file: %s' % err, 'Error', 
                          parent=self, style=wx.OK|wx.ICON_ERROR) 


        
#----------------------------------------------------------------------    
    def OnAlignImgs(self, event):  
        
        wx.GetApp().TopWindow.Hide()
        ImageRegistration(self.com, self.stk).Show()

#----------------------------------------------------------------------    
    def OnSlideshow(self, event):  

        if (self.com.stack_loaded == 1) and (self.addroi == 0):    
            
            if (self.movie_playing == 1):
                self.movie_playing = 0
                return
            
            self.button_slideshow.SetLabel("Stop stack movie")
            old_iev =  self.iev
            self.movie_playing = 1
            
            for i in range(self.stk.n_ev):   
                wx.YieldIfNeeded() 
                if self.movie_playing == 0:
                    break
                self.iev = i                   
                self.loadImage()
                self.slider_eng.SetValue(self.iev)
                self.loadSpectrum(self.ix, self.iy)
                
                time.sleep(0.01)
                
            self.iev = old_iev 
            self.loadImage()
            self.slider_eng.SetValue(self.iev)
            self.loadSpectrum(self.ix, self.iy)
            self.button_slideshow.SetLabel("Play stack movie")
            self.movie_playing = 0
            

                        
#----------------------------------------------------------------------          
    def onrb_fluxod(self, evt):
        state = self.rb_flux.GetValue()
        
      
        if state:
            self.showflux = True
        else:        
            self.showflux = False
            
        self.ResetDisplaySettings()
        self.loadImage()

#----------------------------------------------------------------------           
    def OnShowScale(self, event):
        if self.add_scale_cb.GetValue():
            self.show_scale_bar = 1
        else: 
            self.show_scale_bar = 0
        
        if self.com.stack_loaded == 1:
            self.loadImage()
            
 #----------------------------------------------------------------------           
    def OnShowColBar(self, event):
        if self.add_colbar_cb.GetValue():
            self.show_colorbar = 1
        else: 
            self.show_colorbar = 0
        
        if self.com.stack_loaded == 1:
            self.loadImage()
             
             
#----------------------------------------------------------------------
    def onResetDisplaySettings(self, event):
        
        self.ResetDisplaySettings()
        self.loadImage()
        
#----------------------------------------------------------------------
    def OnScrollBrightnessMin(self, event):
        
        self.dispbrightness_min = event.GetInt()
        
        self.brightness_min = float(self.dispbrightness_min)/100.0
        
        self.defaultdisplay = 0.0
        
        self.tc_min.Clear()
        self.tc_min.AppendText('Minimum: \t{0:5d}%'.format(int(100*self.brightness_min)))
        
        if self.com.stack_loaded == 1:
            self.loadImage()
        
#----------------------------------------------------------------------
    def OnScrollBrightnessMax(self, event):
        
        self.dispbrightness_max = event.GetInt()
        
        self.brightness_max = float(self.dispbrightness_max)/100.0
        
        self.defaultdisplay = 0.0
        
        self.tc_max.Clear()
        self.tc_max.AppendText('Maximum:{0:5d}%'.format(int(100*self.brightness_max)))
        
        if self.com.stack_loaded == 1:
            self.loadImage()
        
#----------------------------------------------------------------------
    def OnScrollGamma(self, event):
        
        self.displaygamma = event.GetInt()
        
        self.gamma = float(self.displaygamma)/10.0  
        
        self.defaultdisplay = 0.0
        
        self.tc_gamma.Clear()
        self.tc_gamma.AppendText('Gamma:  \t{0:5.2f}'.format(self.gamma))
        
        if self.com.stack_loaded == 1:
            self.loadImage()
        
#----------------------------------------------------------------------
    def onSetColorTable(self, event):
        
        ColorTableFrame().Show()
        
#----------------------------------------------------------------------
    def ResetDisplaySettings(self):

        self.defaultdisplay = 1.0
        
        self.dispbrightness_min = 0
        self.dispbrightness_max = 100
        self.displaygamma = 10.0
        
        self.brightness_min = 0.0
        self.brightness_max = 1.0
        self.gamma = 1.0
        
        self.slider_brightness_max.SetValue(self.dispbrightness_max)
        self.slider_brightness_min.SetValue(self.dispbrightness_min) 
        self.slider_gamma.SetValue(self.displaygamma)      

        self.tc_min.Clear()
        self.tc_min.AppendText('Minimum: \t{0:5d}%'.format(int(100*self.brightness_min)))
        self.tc_max.Clear()
        self.tc_max.AppendText('Maximum:{0:5d}%'.format(int(100*self.brightness_max)))        
        self.tc_gamma.Clear()
        self.tc_gamma.AppendText('Gamma:  \t{0:5.2f}'.format(self.gamma))       
 
#----------------------------------------------------------------------        
    def OnLimitEv(self, evt):    
        LimitEv(self.com, self.stk).Show()

#----------------------------------------------------------------------        
    def OnCliptoSubregion(self, evt):    
        CliptoSubregion(self.com, self.stk).Show()
        
                
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
                
        fig = self.SpectrumPanel.get_figure()
        fig.clf()
        self.SpectrumPanel.draw()
        self.tc_spec.SetValue("Average ROI Spectrum: ")
        
        self.button_acceptROI.Disable()
        self.button_resetROI.Enable()
        self.button_ROIdosecalc.Disable() 
        wx.GetApp().TopWindow.refresh_widgets()

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

        fig = self.SpectrumPanel.get_figure()
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        
        mtplot.rcParams['font.size'] = self.fontsize

        specplot = axes.plot(self.stk.ev,self.ROIspectrum)
        
        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Optical Density')
        
        self.SpectrumPanel.draw()
        
        self.tc_spec.SetValue("Average ROI Spectrum: ")
        
#----------------------------------------------------------------------    
    def OnAcceptROI(self, evt):    
        wx.BeginBusyCursor()
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

        self.button_saveROIspectr.Enable()
        self.button_setROII0.Enable()
        self.button_ROIdosecalc.Enable() 
        wx.GetApp().TopWindow.refresh_widgets()
                
        self.loadImage()
        if (self.com.i0_loaded == 1):
            self.ShowROISpectrum()
            
        wx.EndBusyCursor()
        
    
#----------------------------------------------------------------------    
    def OnResetROI(self, evt): 
        self.addroi = 0   
        self.showROImask = 0
        self.ROIpix = None
        
        self.button_acceptROI.Disable()
        self.button_setROII0.Disable()
        self.button_resetROI.Disable()
        self.button_saveROIspectr.Disable()
        self.button_ROIdosecalc.Disable() 
        wx.GetApp().TopWindow.refresh_widgets()
        
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
        
        PlotFrame(self.stk.evi0,self.stk.i0data).Show()              
         
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
        
        self.button_acceptROI.Disable()
        self.button_setROII0.Disable()
        wx.GetApp().TopWindow.refresh_widgets()
        
#----------------------------------------------------------------------    
    def OnROI_DoseCalc(self, event):  
        
        self.CalcROISpectrum()
        
        DoseCalculation(self.com, self.stk, self.ROIspectrum).Show()
        
        
#----------------------------------------------------------------------    
    def OnSaveROISpectrum(self, event):  
               
        fileName = wx.FileSelector('Save ROI Spectrum (.csv)', default_extension='csv', 
                                   wildcard=('csv (*.csv)|*.csv'), 
                                              parent=self, flags=wx.SAVE|wx.OVERWRITE_PROMPT) 
   
        if not fileName: 
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
   
            wx.MessageBox('Could not save file: %s' % err, 'Error', 
                          parent=self, style=wx.OK|wx.ICON_ERROR) 
            
     
#----------------------------------------------------------------------        
    def OnSpectralROI(self, evt):    
        SpectralROI(self.com, self.stk).Show()
        
#----------------------------------------------------------------------        
    def OnDespike(self, evt):    

        image = self.stk.absdata[:,:,self.iev] 
       
        image = self.stk.despike(image)
        
        
        self.stk.data_struct.exchange.data = self.stk.absdata
        
        if self.com.i0_loaded:
            self.stk.calculate_optical_density()
        
        self.loadImage()
        
        
#---------------------------------------------------------------------- 
class SaveWinP1(wx.Frame):
    
    title = "Save"

    def __init__(self):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=self.title, size=(400, 330))
               
        ico = logos.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White") 

        
        self.com = wx.GetApp().TopWindow.common         
        self.fontsize = self.com.fontsize   
        
        path, ext = os.path.splitext(self.com.filename) 
        ext = ext[1:].lower()   
        suffix = "." + ext
        path, fn = os.path.split(self.com.filename)
        filename = fn[:-len(suffix)]
        
        self.path = self.com.path
        self.filename = filename
                          
        
        vboxtop = wx.BoxSizer(wx.VERTICAL)
                
        panel1 = wx.Panel(self, -1)
        
        gridtop = wx.FlexGridSizer(4, 4, vgap=20, hgap=20)
    
        fontb = wx.SystemSettings_GetFont(wx.SYS_DEFAULT_GUI_FONT)
        fontb.SetWeight(wx.BOLD)
        
        
        st1 = wx.StaticText(panel1, -1, 'Save',  style=wx.ALIGN_LEFT)
        st1.SetFont(fontb)
        st2 = wx.StaticText(panel1, -1, '.pdf',  style=wx.ALIGN_LEFT)
        st2.SetFont(fontb)
        st3 = wx.StaticText(panel1, -1, '.png',  style=wx.ALIGN_LEFT)
        st3.SetFont(fontb)
        st3a = wx.StaticText(panel1, -1, '.csv',  style=wx.ALIGN_LEFT)
        st3a.SetFont(fontb)
        
        st4 = wx.StaticText(panel1, -1, '_spectrum',  style=wx.ALIGN_LEFT)
        st4.SetFont(self.com.font)
        
        self.cb1 = wx.CheckBox(panel1, -1, '')
        self.cb1.SetFont(self.com.font)
        self.cb1.Set3StateValue(wx.CHK_CHECKED)
        
        self.cb2 = wx.CheckBox(panel1, -1, '')
        self.cb2.SetFont(self.com.font)
        
        self.cb2a = wx.CheckBox(panel1, -1, '')
        self.cb2a.SetFont(self.com.font)        
        
        st5 = wx.StaticText(panel1, -1, '_image',  style=wx.ALIGN_LEFT)
        st5.SetFont(self.com.font)
        
        self.cb3 = wx.CheckBox(panel1, -1, '')
        self.cb3.SetFont(self.com.font)
        self.cb3.Set3StateValue(wx.CHK_CHECKED)
        
        self.cb4 = wx.CheckBox(panel1, -1, '')
        self.cb4.SetFont(self.com.font)
        
        st6 = wx.StaticText(panel1, -1, 'all images',  style=wx.ALIGN_LEFT)
        st6.SetFont(self.com.font)        

        st7 = wx.StaticText(panel1, -1, ' ',  style=wx.ALIGN_LEFT)
        st7.SetFont(self.com.font) 
        
        self.cb5 = wx.CheckBox(panel1, -1, '')
        self.cb5.SetFont(self.com.font)

        
        gridtop.Add(st1, 0)
        gridtop.Add(st2, 0)
        gridtop.Add(st3, 0)
        gridtop.Add(st3a, 0)
                        
        gridtop.Add(st4, 0)
        gridtop.Add(self.cb1, 0)
        gridtop.Add(self.cb2, 0)  
        gridtop.Add(self.cb2a, 0)             
  
        gridtop.Add(st5, 0)
        gridtop.Add(self.cb3, 0)
        gridtop.Add(self.cb4, 0)  
        gridtop.Add(wx.StaticText(panel1, -1, ' '), 0)
        
        gridtop.Add(st6, 0)
        gridtop.Add(st7, 0)
        gridtop.Add(self.cb5, 0)  
        gridtop.Add(wx.StaticText(panel1, -1, ' '), 0)
        
        panel1.SetSizer(gridtop)
        
        
        panel2 = wx.Panel(self, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        
        hbox0 = wx.BoxSizer(wx.HORIZONTAL)
        
        st0 = wx.StaticText(panel2, -1, 'Filename: ')
        st0.SetFont(self.com.font)
        self.tc_savefn = wx.TextCtrl(panel2, -1,  style=wx.TE_RICH, size =((100, -1)),
                                         value=self.filename)
        self.tc_savefn.SetFont(self.com.font)
        hbox0.Add(st0, 0, flag = wx.ALIGN_CENTER_VERTICAL)
        hbox0.Add(self.tc_savefn,1, wx.EXPAND|wx.LEFT, 2)         
        
        
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
                
        st1 = wx.StaticText(panel2, label='Path: ')
        st1.SetFont(self.com.font)
        self.tc_savepath = wx.TextCtrl(panel2, -1,  style=wx.TE_RICH|wx.TE_READONLY, size =((100, -1)),
                                         value=self.path)
        self.tc_savepath.SetFont(self.com.font)
        hbox1.Add(st1, 0, wx.ALIGN_CENTER_VERTICAL|wx.RIGHT, 20)
        hbox1.Add(self.tc_savepath, 1, wx.EXPAND|wx.LEFT, 2)  
        
        button_path = wx.Button(panel2, -1, 'Browse...')
        self.Bind(wx.EVT_BUTTON, self.OnBrowseDir, id=button_path.GetId())
        hbox1.Add(button_path, 0)
        
        
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        button_save = wx.Button(panel2, -1, 'Save')
        self.Bind(wx.EVT_BUTTON, self.OnSave, id=button_save.GetId())
        hbox2.Add(button_save, 0, wx.ALIGN_RIGHT)
        
        button_cancel = wx.Button(panel2, -1, 'Cancel')
        self.Bind(wx.EVT_BUTTON, self.OnCancel, id=button_cancel.GetId())
        hbox2.Add(button_cancel, 0, wx.ALIGN_RIGHT|wx.LEFT,10)
        
        vbox1.Add(hbox0, 1, wx.EXPAND)
        vbox1.Add((0,5))
        vbox1.Add(hbox1, 1, wx.EXPAND)
        vbox1.Add((0,20))
        vbox1.Add(hbox2, 0, wx.ALIGN_RIGHT)
        panel2.SetSizer(vbox1)
        
        
        vboxtop.Add(panel1, 1, wx.ALL | wx.EXPAND, 20)
        vboxtop.Add(panel2, 0, wx.ALL | wx.EXPAND, 20) 
        
        self.SetSizer(vboxtop)    
        
#----------------------------------------------------------------------        
    def OnBrowseDir(self, evt):
        
        dialog = wx.DirDialog(None, "Choose a directory",
                               style=wx.DD_DIR_MUST_EXIST,
                               defaultPath=self.path)

        if dialog.ShowModal() == wx.ID_OK:
            directory = dialog.GetPath()
            
            self.path = directory
            
            self.tc_savepath.SetValue(self.path)
            
            
                
#----------------------------------------------------------------------        
    def OnSave(self, evt):
        
        self.filename = self.tc_savefn.GetValue()
        
        sp_pdf = self.cb1.GetValue()
        sp_png = self.cb2.GetValue()
        sp_csv = self.cb2a.GetValue()
        im_pdf = self.cb3.GetValue()
        im_png = self.cb4.GetValue()
        im_all = self.cb5.GetValue()
        
        self.Destroy() 
        wx.GetApp().TopWindow.page1.Save(self.filename, self.path,
                                         spec_png = sp_png, 
                                         spec_pdf = sp_pdf, 
                                         sp_csv = sp_csv,
                                         img_png = im_png, 
                                         img_pdf = im_pdf,
                                         img_all = im_all)

#---------------------------------------------------------------------- 
    def OnCancel(self, evt):
        self.Destroy()  
                
                     
#---------------------------------------------------------------------- 
class ShowHistogram(wx.Frame):

    title = "Histogram"

    def __init__(self, stack):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=self.title, size=(630, 700))
               
        ico = logos.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White") 
        
        self.stack = stack
        
        self.com = wx.GetApp().TopWindow.common         
        self.fontsize = self.com.fontsize
        
        self.stack.calc_histogram()
        averagefluxmax = npy.max(self.stack.histogram)
        self.histmin = 0.98*averagefluxmax
        self.histmax = averagefluxmax
        
    
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
               
        i1panel = wx.Panel(panel, -1, style = wx.SUNKEN_BORDER)
        self.HistogramPanel = wxmpl.PlotPanel(i1panel, -1, size=(6.0, 3.7), cursor=False, crosshairs=False, location=False, zoom=False)
        
        wxmpl.EVT_SELECTION(i1panel, self.HistogramPanel.GetId(), self.OnSelection)

        vbox.Add(i1panel, 0, wx.ALL, 20)
        
       
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        sizer2 = wx.StaticBoxSizer(wx.StaticBox(panel, -1, 'I0 pixels'), orient=wx.VERTICAL)
        self.textctrl = wx.TextCtrl(panel, -1, size = (565, 20), style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.textctrl.SetValue('Selection: [ {0:5.2f} kHz, {1:5.2f} kHz ]'.format(float(self.histmin), float(self.histmax)))
        self.textctrl.SetBackgroundColour("White")
        sizer2.Add(self.textctrl, 0)
        

        self.AbsImagePanel = wxmpl.PlotPanel(panel, -1, size =(1.5,1.5), cursor=False, crosshairs=False, location=False, zoom=False)
        sizer2.Add(self.AbsImagePanel, 0)
        hbox2.Add(sizer2, 0, wx.LEFT|wx.RIGHT ,20)
        vbox.Add(hbox2, 0, wx.EXPAND)      
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        
        button_ok = wx.Button(panel, -1, 'Accept')
        self.Bind(wx.EVT_BUTTON, self.OnAccept, id=button_ok.GetId())
        hbox.Add(button_ok, 1, wx.ALL,20)
        
        button_cancel = wx.Button(panel, -1, 'Cancel')
        self.Bind(wx.EVT_BUTTON, self.OnCancel, id=button_cancel.GetId())
        self.Bind(wx.EVT_CLOSE, self.OnCancel)
        hbox.Add(button_cancel, 1, wx.ALL ,20)
        
        vbox.Add(hbox, 0 )
        
        panel.SetSizer(vbox)
        
        vboxtop.Add(panel,1, wx.EXPAND )
        
        self.SetSizer(vboxtop)
        
        self.draw_histogram()
        self.draw_image()
        
        self.Centre()

        
#----------------------------------------------------------------------        
    def draw_histogram(self):
        
     
        fig = self.HistogramPanel.get_figure()
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        self.axes = fig.gca()
        
        mtplot.rcParams['font.size'] = self.fontsize
        
        histdata =  npy.reshape(self.stack.histogram, (self.stack.n_cols*self.stack.n_rows), order='F')
        
        self.n, self.bins, patches = self.axes.hist(histdata, 200, normed=1, facecolor='green', alpha=0.75)
        
        self.axes.set_xlabel('Average Flux [kHz]')
        self.axes.set_ylabel('Percentage of Pixels')

        self.axes.axvspan(self.histmin, self.histmax, facecolor='r', alpha=0.3)

        
        self.HistogramPanel.draw()
        
        
        pass
    
    
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
        
               
        fig = self.AbsImagePanel.get_figure()
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        
        
        axes = fig.gca()
        fig.patch.set_alpha(1.0) 
      
        im = axes.imshow(image, cmap=mtplot.cm.get_cmap("gray")) 

        im_red = axes.imshow(redpix_masked,cmap=mtplot.cm.get_cmap("autumn"))  
         
        axes.axis("off")  
        self.AbsImagePanel.draw()


#----------------------------------------------------------------------        
    def OnSelection(self, evt):
        
        x1, y1 = evt.x1data, evt.y1data
        x2, y2 = evt.x2data, evt.y2data
        
        self.histmin = x1
        self.histmax = x2

        self.textctrl.SetValue('Selection: [ {0:5.2f} kHz, {1:5.2f} kHz ]'.format(float(self.histmin), float(self.histmax)))
    
        self.draw_histogram()
        self.draw_image()

#----------------------------------------------------------------------        
    def OnAccept(self, evt):
        
        self.stack.i0_from_histogram(self.histmin, self.histmax)
        self.Destroy() 
        wx.GetApp().TopWindow.Show()
        wx.GetApp().TopWindow.page1.I0histogramCalculated()


                
#---------------------------------------------------------------------- 
    def OnCancel(self, evt):
        wx.GetApp().TopWindow.Show()
        self.Destroy()   


#---------------------------------------------------------------------- 
class LimitEv(wx.Frame):

    title = "Limit energy range"

    def __init__(self, common, stack):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=self.title, size=(630, 560))
               
        ico = logos.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White") 
        
        self.stack = stack
        self.com = common         
        self.fontsize = self.com.fontsize
        
        self.evlimited = 0
        self.limitevmin = 0
        self.limitevmax = self.stack.n_ev-1
        
        
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
               
        i1panel = wx.Panel(panel, -1, style = wx.SUNKEN_BORDER)
        self.SpectrumPanel = wxmpl.PlotPanel(i1panel, -1, size=(6.0, 3.7), cursor=False, crosshairs=False, location=False, zoom=False)
        
        wxmpl.EVT_SELECTION(i1panel, self.SpectrumPanel.GetId(), self.OnSelection)

        vbox.Add(i1panel, 0, wx.ALL, 20)
        
       
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        sizer2 = wx.StaticBoxSizer(wx.StaticBox(panel, -1, 'Energy'), orient=wx.VERTICAL)
        self.textctrl = wx.TextCtrl(panel, -1, size = (565, 40), style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.textctrl.SetValue(' ')
        self.textctrl.SetBackgroundColour("White")
        sizer2.Add(self.textctrl, 0)
        hbox2.Add(sizer2, 0, wx.LEFT|wx.RIGHT ,20)
        vbox.Add(hbox2, 0, wx.EXPAND)
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        
        button_ok = wx.Button(panel, -1, 'Accept')
        self.Bind(wx.EVT_BUTTON, self.OnAccept, id=button_ok.GetId())
        hbox.Add(button_ok, 1, wx.ALL,20)
        
        button_cancel = wx.Button(panel, -1, 'Cancel')
        self.Bind(wx.EVT_BUTTON, self.OnCancel, id=button_cancel.GetId())
        self.Bind(wx.EVT_CLOSE, self.OnCancel)
        hbox.Add(button_cancel, 1, wx.ALL ,20)
        
        vbox.Add(hbox, 0 )
        
        panel.SetSizer(vbox)
        
        vboxtop.Add(panel,1, wx.EXPAND )
        
        self.SetSizer(vboxtop)
        
        self.Centre()
        
        self.draw_limitev_plot()
        
      
#----------------------------------------------------------------------        
    def draw_limitev_plot(self):
        
        odtotal = self.stack.od3d.sum(axis=0)   
        odtotal = odtotal.sum(axis=0)/(self.stack.n_rows*self.stack.n_cols) 
        
        fig = self.SpectrumPanel.get_figure()
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        self.axes = fig.gca()
        
        mtplot.rcParams['font.size'] = self.fontsize

        specplot = self.axes.plot(self.stack.ev,odtotal)
        
        self.axes.set_xlabel('Photon Energy [eV]')
        self.axes.set_ylabel('Optical Density')
        
        if self.evlimited == 1:
            self.axes.axvspan(self.stack.ev[self.limitevmin], self.stack.ev[self.limitevmax], facecolor='g', alpha=0.5)

        
        self.SpectrumPanel.draw()
        
        self.textctrl.SetValue('Min energy {0:5.2f} eV\n'.format(float(self.stack.ev[self.limitevmin]))
                              + 'Max energy {0:5.2f} eV'.format(float(self.stack.ev[self.limitevmax])))

#----------------------------------------------------------------------        
    def OnSelection(self, evt):

        x1, y1 = evt.x1data, evt.y1data
        x2, y2 = evt.x2data, evt.y2data

        
        
        self.limitevmin = npy.abs(self.stack.ev-x1).argmin()
        self.limitevmax = npy.abs(self.stack.ev-x2).argmin()
        
        self.evlimited = 1
             
        self.draw_limitev_plot()
        
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
        wx.GetApp().TopWindow.page1.slider_eng.SetRange(0,self.stack.n_ev-1)
        wx.GetApp().TopWindow.page1.iev = self.stack.n_ev/2
        wx.GetApp().TopWindow.page1.slider_eng.SetValue(wx.GetApp().TopWindow.page1.iev)
        
        wx.GetApp().TopWindow.page1.loadSpectrum(wx.GetApp().TopWindow.page1.ix, wx.GetApp().TopWindow.page1.iy)
        wx.GetApp().TopWindow.page1.loadImage()
        
        self.Destroy()
        
#---------------------------------------------------------------------- 
    def OnCancel(self, evt):
        self.Destroy()
        

#---------------------------------------------------------------------- 
class CliptoSubregion(wx.Frame):

    title = "Clip to Subregion"

    def __init__(self, common, stack):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=self.title, size=(500, 600))
               
        ico = logos.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White") 
        
        self.stack = stack
        
        self.com = common         
        self.fontsize = self.com.fontsize
        
        self.new_y1 = int(self.stack.n_cols*0.10)
        self.new_y2 = int(self.stack.n_cols*0.90)
        self.new_x1 = int(self.stack.n_rows*0.10)
        self.new_x2 = int(self.stack.n_rows*0.90)
        
        self.new_ncols = self.new_y2 - self.new_y1
        self.new_nrows = self.new_x2 - self.new_x1
    
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
        
        vbox.Add((0,40))
               
       
       
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        sizer2 = wx.StaticBoxSizer(wx.StaticBox(panel, -1, 'Select new stack size'), orient=wx.VERTICAL)
        self.textctrl1 = wx.TextCtrl(panel, -1, size = (200, 20), style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.textctrl1.SetValue('Original stack size:\t{0:5d}   x{1:5d} '.format(self.stack.n_cols, self.stack.n_rows))
        self.textctrl1.SetBackgroundColour("White")
        sizer2.Add(self.textctrl1, 0)
 
        self.textctrl2 = wx.TextCtrl(panel, -1, size = (200, 20), style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.textctrl2.SetValue('New stack size:\t{0:5d}   x{1:5d} '.format(self.new_ncols, self.new_nrows))
        self.textctrl2.SetBackgroundColour("White")
        sizer2.Add(self.textctrl2, 0)       
        
        self.textctrl3 = wx.TextCtrl(panel, -1, size = (400, 20), style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.textctrl3.SetValue('Clip coordinates [[x1, x2], [y1, y2]] : [[{0:5d},{1:5d}], [{2:5d},{3:5d}]]'.format(
                                    self.new_x1, self.new_x2, self.new_y1, self.new_y2))
        self.textctrl3.SetBackgroundColour("White")
        sizer2.Add(self.textctrl3, 0)    

        self.AbsImagePanel = wxmpl.PlotPanel(panel, -1, size =(PlotH,PlotH), cursor=False, crosshairs=False, location=False, zoom=False)
        wxmpl.EVT_SELECTION(panel, self.AbsImagePanel.GetId(), self.OnSelection)
        sizer2.Add(self.AbsImagePanel, 0)
        hbox2.Add(sizer2, 0, wx.LEFT|wx.RIGHT ,20)
        vbox.Add(hbox2, 0, wx.EXPAND)              
        vbox.Add((0,10))
         
         
        hbox = wx.BoxSizer(wx.HORIZONTAL)
              
        button_ok = wx.Button(panel, -1, 'Accept')
        self.Bind(wx.EVT_BUTTON, self.OnAccept, id=button_ok.GetId())
        hbox.Add(button_ok, 1, wx.ALL,20)
        
        button_cancel = wx.Button(panel, -1, 'Cancel')
        self.Bind(wx.EVT_BUTTON, self.OnCancel, id=button_cancel.GetId())
        self.Bind(wx.EVT_CLOSE, self.OnCancel)
        hbox.Add(button_cancel, 1, wx.ALL ,20)
        
        vbox.Add(hbox, 0 )
        
        panel.SetSizer(vbox)
        
        vboxtop.Add(panel,1, wx.EXPAND | wx.LEFT, 20 )
        
        self.SetSizer(vboxtop)
        

        self.draw_image()
        
        self.Centre()

    
    
#----------------------------------------------------------------------        
    def draw_image(self):
        
   
        image = self.stack.absdata[:,:,self.stack.n_ev/2].copy() 
       
       
               
        fig = self.AbsImagePanel.get_figure()
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        
        
        axes = fig.gca()
        fig.patch.set_alpha(1.0) 
      
        im = axes.imshow(image, cmap=mtplot.cm.get_cmap("gray")) 
        
        # Draw the rectangle
        line1=mtplot.lines.Line2D([self.new_x1,self.new_x2], [self.new_y1,self.new_y1] ,color="red")
        line1.set_clip_on(False)
        axes.add_line(line1)   

        line2=mtplot.lines.Line2D([self.new_x1,self.new_x2], [self.new_y2,self.new_y2] ,color="red")
        line2.set_clip_on(False)
        axes.add_line(line2)   
        
        line3=mtplot.lines.Line2D([self.new_x1,self.new_x1], [self.new_y1,self.new_y2] ,color="red")
        line3.set_clip_on(False)
        axes.add_line(line3)   
        
        line4=mtplot.lines.Line2D([self.new_x2,self.new_x2], [self.new_y1,self.new_y2] ,color="red")
        line4.set_clip_on(False)
        axes.add_line(line4)   
         
        axes.axis("off")  
        self.AbsImagePanel.draw()


#----------------------------------------------------------------------        
    def OnSelection(self, evt):
        
        x1, y1 = evt.x1data, evt.y1data
        x2, y2 = evt.x2data, evt.y2data
        
        self.new_y1 = int(y1)
        self.new_y2 = int(y2)
        self.new_x1 = int(x1)
        self.new_x2 = int(x2)
        
        self.new_ncols = self.new_y1 - self.new_y2 + 1
        self.new_nrows = self.new_x2 - self.new_x1 + 1

        self.textctrl2.SetValue('New stack size:\t{0:5d}   x{1:5d} '.format(self.new_ncols, self.new_nrows))
        self.textctrl3.SetValue('Clip coordinates [[x1, x2], [y1, y2]]:[[{0:5d},{1:5d}], [{2:5d},{3:5d}]]'.format(
                                    self.new_x1, self.new_x2, self.new_y1, self.new_y2))
    
        self.draw_image()

#----------------------------------------------------------------------        
    def OnAccept(self, evt):
        
        #change the stack size to [x1,x2], [y1,y2] 
        #Matlab axis are inverted
        self.stack.absdata = self.stack.absdata[self.new_y2:self.new_y1+1, self.new_x1:self.new_x2+1, :]
              
        self.stack.n_cols = self.stack.absdata.shape[0]
        self.stack.n_rows = self.stack.absdata.shape[1]

        
        if self.com.i0_loaded == 1:        
            self.stack.od3d = self.stack.od3d[self.new_y2:self.new_y1+1, self.new_x1:self.new_x2+1, :]
        
            self.stack.od = self.stack.od3d.copy()
        
            self.stack.od = npy.reshape(self.stack.od, (self.stack.n_rows*self.stack.n_cols, self.stack.n_ev), order='F')
        
        #Fix the slider on Page 1! 
        wx.GetApp().TopWindow.page1.ix = self.stack.n_cols/2
        wx.GetApp().TopWindow.page1.iy = self.stack.n_rows/2
        
        self.stack.fill_h5_struct_from_stk()
        
        wx.GetApp().TopWindow.page1.loadSpectrum(wx.GetApp().TopWindow.page1.ix, wx.GetApp().TopWindow.page1.iy)
        wx.GetApp().TopWindow.page1.loadImage()
        wx.GetApp().TopWindow.page0.ShowImage()


        self.Destroy() 
        wx.GetApp().TopWindow.Show()



                
#---------------------------------------------------------------------- 
    def OnCancel(self, evt):
        wx.GetApp().TopWindow.Show()
        self.Destroy()   


        
#---------------------------------------------------------------------- 
class ImageRegistration(wx.Frame):

    title = "Image Alignment"

    def __init__(self, common, stack):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=self.title, size=(920, 750))

               
        ico = logos.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White") 
        
        self.stack = stack
        self.com = common         
        self.fontsize = self.com.fontsize
        
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
                                  
        
        #panel 1        
        panel1 = wx.Panel(self, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_imageeng = wx.TextCtrl(panel1, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_imageeng.SetFont(self.com.font)
        self.tc_imageeng.SetValue("Image at energy: ")
       
        hbox11 = wx.BoxSizer(wx.HORIZONTAL)
   
        i1panel = wx.Panel(panel1, -1, style = wx.SUNKEN_BORDER)
        self.AbsImagePanel = wxmpl.PlotPanel(i1panel, -1, size =(3.0,3.0), cursor=True, crosshairs=False, location=False, zoom=False)
        wxmpl.EVT_POINT(i1panel, self.AbsImagePanel.GetId(), self.OnPointCorrimage)

                    
        self.slider_eng = wx.Slider(panel1, -1, self.iev, 0, self.stack.n_ev-1, style=wx.SL_LEFT|wx.SL_VERTICAL)        
        self.slider_eng.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnScrollEng, self.slider_eng)

        hbox11.Add(i1panel, 0)
        hbox11.Add(self.slider_eng, 0,  wx.EXPAND)
        
        vbox1.Add(self.tc_imageeng,1, wx.EXPAND)        
        vbox1.Add(hbox11, 0)

        panel1.SetSizer(vbox1)
        
        
        
        #panel 2        
        panel2 = wx.Panel(self, -1)
        vbox2 = wx.BoxSizer(wx.VERTICAL)
        
        tc2 = wx.TextCtrl(panel2, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        tc2.SetFont(self.com.font)
        tc2.SetValue('Cross-correlation')

        i2panel = wx.Panel(panel2, -1, style = wx.SUNKEN_BORDER)
        self.CscorrPanel = wxmpl.PlotPanel(i2panel, -1, size =(2.4,2.4), cursor=False, crosshairs=True, location=False, zoom=False)
                                      
        vbox2.Add(tc2,1, wx.EXPAND)        
        vbox2.Add(i2panel, 0)

        panel2.SetSizer(vbox2)
        
        
        #panel 3
        panel3 = wx.Panel(self, -1)
        vbox3 = wx.BoxSizer(wx.VERTICAL)
        
        tc3= wx.TextCtrl(panel3, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        tc3.SetValue('Image shifts')
        tc3.SetFont(self.com.font)
          
        i4panel = wx.Panel(panel3, -1, style = wx.SUNKEN_BORDER)
        self.ShiftsPanel = wxmpl.PlotPanel(i4panel, -1, size=(4.0, 2.4), cursor=False, crosshairs=False, location=False, zoom=False)
        wxmpl.EVT_POINT(i4panel, self.ShiftsPanel.GetId(), self.OnPlotShifts)
        
        vbox3.Add(tc3, 1, wx.EXPAND)       
        vbox3.Add(i4panel, 0)
        
        panel3.SetSizer(vbox3)
        
        #panel 9
        panel9 = wx.Panel(self, -1)     
        vbox9 = wx.BoxSizer(wx.VERTICAL)   
        
        self.rb_auto = wx.RadioButton(panel9, -1, 'Automatic Alignment', style=wx.RB_GROUP)
        self.rb_man = wx.RadioButton(panel9, -1, 'Manual Alignment')
        self.rb_auto.SetFont(self.com.font)
        self.rb_man.SetFont(self.com.font)
        self.Bind(wx.EVT_RADIOBUTTON, self.Onrb_automanual, id=self.rb_auto.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.Onrb_automanual, id=self.rb_man.GetId())
        
        vbox9.Add(self.rb_auto, 1, wx.EXPAND)
        vbox9.Add((0,3))
        vbox9.Add(self.rb_man, 1, wx.EXPAND)
        panel9.SetSizer(vbox9)
        
        
        
        #panel 8
        panel8 = wx.Panel(self, -1)
        sizer8 = wx.StaticBoxSizer(wx.StaticBox(panel8, -1, 'This Image'),orient=wx.VERTICAL)
        vbox81 = wx.BoxSizer(wx.VERTICAL)
        vbox81.Add((0,3))        

        self.button_refimg = wx.Button(panel8, -1, 'Set as Reference Image')
        self.button_refimg.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.SetRefImage, id=self.button_refimg.GetId())
        vbox81.Add(self.button_refimg, 0, wx.EXPAND)
        
        self.button_remove = wx.Button(panel8, -1, 'Remove image')
        self.button_remove.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnRemoveImage, id=self.button_remove.GetId())    
        vbox81.Add(self.button_remove, 0, wx.EXPAND)
        
        
        
        sizer8.Add(vbox81,1, wx.LEFT|wx.RIGHT|wx.EXPAND,2)
        panel8.SetSizer(sizer8)
        
        
        
        #panel 4
        panel4 = wx.Panel(self, -1)
        sizer1 = wx.StaticBoxSizer(wx.StaticBox(panel4, -1, 'Automatic Alignment'),orient=wx.VERTICAL)
        vbox41 = wx.BoxSizer(wx.VERTICAL)
        vbox41.Add((0,3))        

        
        self.button_register = wx.Button(panel4, -1, 'Calculate image shifts')
        self.button_register.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnCalcRegistration, id=self.button_register.GetId())   
        self.button_register.Disable()     
        vbox41.Add(self.button_register, 0, wx.EXPAND)
        
        vbox41.Add((0,8))
        
        self.button_subregion = wx.Button(panel4, -1, 'Select subregion on reference')
        self.button_subregion.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnSelectSubregion, id=self.button_subregion.GetId())   
        self.button_subregion.Disable()     
        vbox41.Add(self.button_subregion, 0, wx.EXPAND)
        
        self.button_delsubregion = wx.Button(panel4, -1, 'Remove subregion selection')
        self.button_delsubregion.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnDeleteSubregion, id=self.button_delsubregion.GetId())   
        self.button_delsubregion.Disable()     
        vbox41.Add(self.button_delsubregion, 0, wx.EXPAND)
        
        
        vbox41.Add((0,10))
        
        hbox42 = wx.BoxSizer(wx.HORIZONTAL)
        text1 = wx.StaticText(panel4, label=' Max shift [pixels]: ')
        self.tc_maxshift = wx.lib.intctrl.IntCtrl(panel4, size=( 70, -1 ), 
                                                  style = wx.TE_CENTRE, 
                                                  value = 0, 
                                                  limited = True )
        self.tc_maxshift.SetMin(0)
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnSetMaxShift, self.tc_maxshift)
        
        
        hbox42.Add((5,0))
        hbox42.Add(text1, 0)
        hbox42.Add(self.tc_maxshift,0)
        vbox41.Add(hbox42,0, wx.EXPAND)
        
        
        self.tc_shift = wx.TextCtrl(panel4, -1, style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_shift.SetFont(self.com.font)
        vbox41.Add(self.tc_shift, 1, wx.EXPAND|wx.TOP|wx.LEFT, 5)
        
        
        self.tc_shift.AppendText('X shift: {0:5.2f} pixels\n'.format(self.yshifts[self.iev]))
        self.tc_shift.AppendText('Y shift: {0:5.2f} pixels'.format(self.xshifts[self.iev]))
        
        
        vbox41.Add((0,3))
        
        hbox41 = wx.BoxSizer(wx.HORIZONTAL)
        hbox41.Add((5,0))
        self.showcscor_cb = wx.CheckBox(panel4, -1, 'Show Cross-correlation')
        self.showcscor_cb.SetFont(self.com.font)
        self.Bind(wx.EVT_CHECKBOX, self.OnShowCCorr, self.showcscor_cb)
        hbox41.Add(self.showcscor_cb, 0)
        hbox41.Add((5,0))
        vbox41.Add(hbox41,0, wx.EXPAND)
        vbox41.Add((0,3))

        
        
        sizer1.Add(vbox41,1, wx.LEFT|wx.RIGHT|wx.EXPAND,2)
        panel4.SetSizer(sizer1)
        
        
        #panel 5        
        panel5 = wx.Panel(self, -1)
        vbox5 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_refimg = wx.TextCtrl(panel5, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_refimg.SetFont(self.com.font)
        self.tc_refimg.SetValue('Reference image')
        
        i3panel = wx.Panel(panel5, -1, style = wx.SUNKEN_BORDER)
        self.RefImagePanel = wxmpl.PlotPanel(i3panel, -1, size =(3.0,3.0), cursor=True, crosshairs=False, location=False, zoom=False)
        wxmpl.EVT_POINT(i3panel, self.RefImagePanel.GetId(), self.OnPointRefimage)
        
        wxmpl.EVT_SELECTION(i3panel, self.RefImagePanel.GetId(), self.OnSelection)
        
        vbox5.Add(self.tc_refimg,1, wx.EXPAND)        
        vbox5.Add(i3panel, 0)

        panel5.SetSizer(vbox5)
        
        
        #panel 6
        panel6 = wx.Panel(self, -1)
        sizer6 = wx.StaticBoxSizer(wx.StaticBox(panel6, -1, 'Manual Alignment'),orient=wx.VERTICAL)
        vbox61 = wx.BoxSizer(wx.VERTICAL)
        vbox61.Add((0,3))        
                
        self.button_manalign = wx.Button(panel6, -1, 'Pick a point on reference image')
        self.button_manalign.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnPickRefPoint, id=self.button_manalign.GetId())
        self.button_manalign.Disable()
        vbox61.Add(self.button_manalign, 0, wx.EXPAND)
        
        self.button_pick2ndpoint = wx.Button(panel6, -1, 'This image: click on same point')
        self.button_pick2ndpoint.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnPickCorrPoint, id=self.button_pick2ndpoint.GetId())
        self.button_pick2ndpoint.Disable()
        vbox61.Add(self.button_pick2ndpoint, 0, wx.EXPAND)
        
        self.button_applyman = wx.Button(panel6, -1, 'Apply manual shifts')
        self.button_applyman.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnApplyManShifts, id=self.button_applyman.GetId())
        self.button_applyman.Disable()
        vbox61.Add(self.button_applyman, 0, wx.EXPAND)
        
        
        self.textctrl_ms = wx.TextCtrl(panel6, -1, style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.textctrl_ms.SetFont(self.com.font)
        vbox61.Add(self.textctrl_ms, 1, wx.EXPAND|wx.TOP, 10)
        
        vbox61.Add((0,3))
        
        self.textctrl_ms.AppendText('X manual shift: \n')
        self.textctrl_ms.AppendText('Y manual shift: ')
        
        
        
        sizer6.Add(vbox61,1, wx.LEFT|wx.RIGHT|wx.EXPAND,2)
        panel6.SetSizer(sizer6)
        
        #panel 7
        panel7 = wx.Panel(self, -1)
        vbox7 = wx.BoxSizer(wx.VERTICAL)
        
        self.button_saveimg = wx.Button(panel7, -1, 'Save image shifts plot')
        self.button_saveimg.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnSaveShiftsPlot, id=self.button_saveimg.GetId())   
        self.button_saveimg.Disable()
        vbox7.Add(self.button_saveimg, 0, wx.EXPAND)
        
        self.button_crop = wx.Button(panel7, -1, 'Crop aligned images')
        self.button_crop.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnCropShifts, id=self.button_crop.GetId())   
        self.button_crop.Disable()
        vbox7.Add(self.button_crop, 0, wx.EXPAND)

        self.button_saveshifts = wx.Button(panel7, -1, 'Save image shifts')
        self.button_saveshifts.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnSaveShifts, id=self.button_saveshifts.GetId())
        self.button_saveshifts.Disable()
        vbox7.Add(self.button_saveshifts, 0, wx.EXPAND)
        
        self.button_loadshifts = wx.Button(panel7, -1, 'Load image shifts')
        self.button_loadshifts.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnLoadShifts, id=self.button_loadshifts.GetId())
        vbox7.Add(self.button_loadshifts, 0, wx.EXPAND)
                
        self.button_accept = wx.Button(panel7, -1, 'Accept changes')
        self.button_accept.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnAccept, id=self.button_accept.GetId())
        self.button_accept.Disable()
        vbox7.Add(self.button_accept, 0, wx.EXPAND)
        
        self.button_close = wx.Button(panel7, -1, 'Dismiss changes')
        self.button_close.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=self.button_close.GetId())
        self.Bind(wx.EVT_CLOSE, self.OnClose)


        vbox7.Add(self.button_close, 0, wx.EXPAND)
        
        panel7.SetSizer(vbox7)
        
        
               
        
        hboxtop = wx.BoxSizer(wx.HORIZONTAL)
        
        vboxL = wx.BoxSizer(wx.VERTICAL)
        vboxR = wx.BoxSizer(wx.VERTICAL)
        
        vboxL.Add((0,10))
        vboxL.Add(panel8, 0, wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM, 15)
        vboxL.Add(panel9, 0, wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM, 15)
        vboxL.Add(panel4, 1, wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM, 15)
        vboxL.Add(panel6, 0, wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM, 15)
        vboxL.Add(panel7, 1, wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM, 15)
        
        hboxRT = wx.BoxSizer(wx.HORIZONTAL)
        hboxRB = wx.BoxSizer(wx.HORIZONTAL)
        
        hboxRT.Add(panel1, 0, wx.RIGHT, 25)
        hboxRT.Add(panel5)
        
        hboxRB.Add(panel3, 0, wx.RIGHT, 10)
        hboxRB.Add(panel2)
        
        vboxR.Add((0,20))
        vboxR.Add(hboxRT)
        vboxR.Add((0,30))
        vboxR.Add(hboxRB)
        
        hboxtop.Add(vboxL)
        hboxtop.Add((20,0))
        hboxtop.Add(vboxR)
        
        self.SetSizer(hboxtop)
        
        self.Centre()
        
        self.ShowImage()

        
        
#----------------------------------------------------------------------        
    def ShowImage(self):
               
        image = self.aligned_stack[:,:,self.iev]
            
            
        fig = self.AbsImagePanel.get_figure()
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        axes = fig.gca()  
             

        im = axes.imshow(image, cmap=mtplot.cm.get_cmap('gray')) 
       
        if self.man_align == 2:           
            lx=mtplot.lines.Line2D([self.man_yref-self.man_ys[self.iev],
                                    self.man_yref-self.man_ys[self.iev]], 
                                    [0,self.stack.n_cols],color='red')
            ly=mtplot.lines.Line2D([0,self.stack.n_rows], 
                                   [self.man_xref-self.man_xs[self.iev],
                                    self.man_xref-self.man_xs[self.iev]] ,color='red')
            axes.add_line(lx)
            axes.add_line(ly)
            
        axes.axis("off") 
        self.AbsImagePanel.draw()
        
        
        self.tc_imageeng.SetValue('Image at energy: {0:5.2f} eV'.format(float(self.stack.ev[self.iev])))
        
        self.tc_shift.Clear()
        self.tc_shift.AppendText('X shift: {0:5.2f} pixels\n'.format(self.yshifts[self.iev]))
        self.tc_shift.AppendText('Y shift: {0:5.2f} pixels'.format(self.xshifts[self.iev]))

        
        if (self.man_align == 2):                  
            self.textctrl_ms.Clear()
            self.textctrl_ms.AppendText('X manual shift:  {0:5.2f}  pixels\n'.format(self.man_ys[self.iev]))
            self.textctrl_ms.AppendText('Y manual shift:  {0:5.2f}  pixels'.format(self.man_xs[self.iev]))
            

        
        
#----------------------------------------------------------------------            
    def OnScrollEng(self, event):
        self.iev = event.GetInt()
        
        self.ShowImage()
            

#----------------------------------------------------------------------        
    def SetRefImage(self, event):
        
        self.ref_image_index = self.iev
               
        self.ref_image = self.aligned_stack[:,:,self.iev].copy()
        
        self.ShowRefImage()
        self.have_ref_image = 1
        
        self.UpdateWidgets()
        

#----------------------------------------------------------------------        
    def ShowRefImage(self):

        fig = self.RefImagePanel.get_figure()
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        axes = fig.gca()
        
    
        im = axes.imshow(self.ref_image, cmap=mtplot.cm.get_cmap('gray')) 
        
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
        
            patch = patches.PathPatch(path, facecolor='red')
            patch.set_alpha(0.3)
            axes.add_patch(patch)
         
        if self.man_align > 0:           
            lx=mtplot.lines.Line2D([self.man_yref,self.man_yref], [0,self.stack.n_cols],color='green')
            ly=mtplot.lines.Line2D([0,self.stack.n_rows], [self.man_xref,self.man_xref] ,color='green')
            axes.add_line(lx)
            axes.add_line(ly)
    
        axes.axis("off") 

        self.RefImagePanel.draw()
        
        self.tc_refimg.SetValue('Reference image at energy: {0:5.2f} eV'.format(float(self.stack.ev[self.ref_image_index])))

#----------------------------------------------------------------------          
    def Onrb_automanual(self, evt):
        state = self.rb_auto.GetValue()
        
      
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
            fig = self.ShiftsPanel.get_figure()
            fig.clf()
            self.ShiftsPanel.draw()
            
            fig = self.CscorrPanel.get_figure()
            fig.clf()
            self.CscorrPanel.draw()
            
            self.ShowImage()
            self.ShowRefImage()
            
        self.UpdateWidgets()
        
#----------------------------------------------------------------------        
    def ShowCrossCorrelation(self, ccorr, xshift, yshift):

        fig = self.CscorrPanel.get_figure()
    
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        axes = fig.gca()
    
        im = axes.imshow(ccorr, cmap=mtplot.cm.get_cmap('gray')) 
        
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
                   
        lx=mtplot.lines.Line2D([ycenter,ycenter], [xl,xr],color='green')
        ly=mtplot.lines.Line2D([yl,yr], [xcenter,xcenter] ,color='green')
        axes.add_line(lx)
        axes.add_line(ly)
         
    
        axes.axis("off") 

        self.CscorrPanel.draw()        
        
#----------------------------------------------------------------------           
    def OnShowCCorr(self, event):
        if self.showcscor_cb.GetValue():
            self.showccorr = 1
        else: 
            self.showccorr = 0


#----------------------------------------------------------------------           
    def OnSetMaxShift(self, event):
        ctl = event.GetEventObject()
        self.maxshift = ctl.GetValue()      

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
        
        
        wx.GetApp().TopWindow.page1.slider_eng.SetRange(0,self.stack.n_ev-1)
        wx.GetApp().TopWindow.page1.iev = self.stack.n_ev/2
        wx.GetApp().TopWindow.page1.slider_eng.SetValue(wx.GetApp().TopWindow.page1.iev)
        
        wx.GetApp().TopWindow.page1.loadSpectrum(wx.GetApp().TopWindow.page1.ix, wx.GetApp().TopWindow.page1.iy)
        wx.GetApp().TopWindow.page1.loadImage() 

#----------------------------------------------------------------------            
    def OnCalcRegistration(self, event):
        
        wx.BeginBusyCursor()
        
        #Subregion selection on a reference image
        if self.subregion == 0:
            self.sr_x1 = 0
            self.sr_x2 = 0
            self.sr_y1 = 0
            self.sr_y2 = 0            
            referenceimage = self.ref_image            
        else:    
            referenceimage = self.ref_image[self.sr_y1:self.sr_y2, self.sr_x1:self.sr_x2]  
          

        for i in range(self.stack.n_ev):

            if self.subregion == 0:
                img2 = self.aligned_stack[:,:,i]  
            else:
                img2 = self.aligned_stack[self.sr_y1:self.sr_y2, self.sr_x1:self.sr_x2, i]  
               
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
        self.slider_eng.SetValue(self.iev)
        
        self.UpdateWidgets()
        
            
        wx.EndBusyCursor()
        
#----------------------------------------------------------------------            
    def OnCropShifts(self, event):
        
        wx.BeginBusyCursor()
        
        self.aligned_stack = self.stack.crop_registed_images(self.aligned_stack, 
                                                             self.xshifts,
                                                             self.yshifts)
        
        
        self.iev = 0
        self.ShowImage()
        self.slider_eng.SetValue(self.iev)
        
        wx.EndBusyCursor()
        
        
#----------------------------------------------------------------------            
    def OnPickRefPoint(self, event):    
        self.man_align = 1   

    
#----------------------------------------------------------------------  
    def OnPointRefimage(self, evt):
        
        if self.man_align == 0:
            pass   
                
        x = evt.xdata
        y = evt.ydata
        
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
    def OnPickCorrPoint(self, event):    

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
                
            
            self.textctrl_ms.Clear()
            self.textctrl_ms.AppendText('X manual shift:  {0:5.2f}  pixels\n'.format(self.man_ys[self.iev]))
            self.textctrl_ms.AppendText('Y manual shift:  {0:5.2f}  pixels'.format(self.man_xs[self.iev]))
            
            self.iev = self.iev + 1
            if self.iev > (self.stack.n_ev-1):
                self.iev = 0
            
            self.slider_eng.SetValue(self.iev)
                    
           
            self.ShowImage()
            self.UpdateWidgets()
            
#----------------------------------------------------------------------            
    def OnApplyManShifts(self, event):
        
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
        
        self.textctrl_ms.Clear()
        self.textctrl_ms.AppendText('X manual shift: \n')
        self.textctrl_ms.AppendText('Y manual shift: ')
        
        self.UpdateWidgets()
        

        self.ShowImage()
        
        
      
                    
#----------------------------------------------------------------------            
    def OnAccept(self, event):
        
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
        
        
        wx.GetApp().TopWindow.page1.slider_eng.SetRange(0,self.stack.n_ev-1)
        wx.GetApp().TopWindow.page1.iev = self.stack.n_ev/2
        wx.GetApp().TopWindow.page1.slider_eng.SetValue(wx.GetApp().TopWindow.page1.iev)
        
        wx.GetApp().TopWindow.page1.ix = int(self.stack.n_cols/2)
        wx.GetApp().TopWindow.page1.iy = int(self.stack.n_rows/2)
        
        wx.GetApp().TopWindow.page1.loadSpectrum(wx.GetApp().TopWindow.page1.ix, wx.GetApp().TopWindow.page1.iy)
        wx.GetApp().TopWindow.page1.loadImage()
        
        wx.GetApp().TopWindow.Show()
        self.Destroy()

        
#----------------------------------------------------------------------              
    def OnClose(self, evt):

        wx.GetApp().TopWindow.Show()
        self.Destroy()

        
        
        
#----------------------------------------------------------------------          
    def PlotShifts(self):
        
        fig = self.ShiftsPanel.get_figure()
        fig.clf()
        
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        
        
        mtplot.rcParams['font.size'] = self.fontsize
        
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
            
            self.slider_eng.SetValue(self.iev)
            
#----------------------------------------------------------------------  
    def OnSaveShiftsPlot(self, evt):
        wildcard = 'Portable Network Graphics (*.png)|*.png|Adobe PDF Files (*.pdf)|*.pdf|All files (*.*)|*.*'               
        self.SaveFileName = wx.FileSelector('Save Image Shifts Plot', default_extension='png', 
                                   wildcard=wildcard, parent=self, flags=wx.SAVE|wx.OVERWRITE_PROMPT) 
   
        if not self.SaveFileName: 
            return 
        
        path, ext = os.path.splitext(self.SaveFileName) 
        ext = ext[1:].lower() 
        
        try: 

            mtplot.rcParams['pdf.fonttype'] = 42
            if ext == 'png':
                            
                fig = self.ShiftsPanel.get_figure()
                fig.savefig(self.SaveFileName)
                
            if ext =='pdf':

                            
                fig = self.ShiftsPanel.get_figure()
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
        
        self.button_delsubregion.Enable()
        
#----------------------------------------------------------------------            
    def OnDeleteSubregion(self, event):
        
        self.subregion = 0
        self.sr_x1 = 0
        self.sr_x2 = 0
        self.sr_y1 = 0
        self.sr_y2 = 0
        
        self.ShowRefImage()
        
                    
#----------------------------------------------------------------------        
    def OnSelection(self, evt):
        
        if (self.man_align > 0) and (self.subregion == 0):
            pass
        
        x1, y1 = evt.x1data, evt.y1data
        x2, y2 = evt.x2data, evt.y2data
        
        self.sr_x1 = int(x1)
        self.sr_x2 = int(x2)
        self.sr_y1 = int(y2)
        self.sr_y2 = int(y1)
        
        self.ShowRefImage()

#----------------------------------------------------------------------  
    def OnSaveShifts(self, evt):
        
        wildcard = "CSV files (*.csv)|*.csv"
        dialog = wx.FileDialog(None, "Please select an alignment file (.csv)",
                               wildcard=wildcard,
                               style=wx.SAVE)
        if dialog.ShowModal() == wx.ID_OK:
            filepath = dialog.GetPath()
            
            file = open(filepath, 'w')
            print>>file, '*********************  Alignment file  ********************'
            print>>file, '***  for ', self.com.filename
            print>>file, '***  ev, xshift, yshift'           
            for ie in range(self.stack.n_ev):
                print>>file, '%.6f, %.6f, %.6f' %(self.stack.ev[ie], self.xshifts[ie], self.yshifts[ie])
        
            file.close()            
        
#----------------------------------------------------------------------  
    def OnLoadShifts(self, evt):
        wildcard = "CSV files (*.csv)|*.csv"
        dialog = wx.FileDialog(None, "Please select an alignment file (.csv)",
                               wildcard=wildcard,
                               style=wx.OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            filepath = dialog.GetPath()


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
            self.slider_eng.SetValue(self.iev)
            
            self.UpdateWidgets()                    
           
            
#----------------------------------------------------------------------  
    def UpdateWidgets(self):
        
        if self.auto:
            self.button_manalign.Disable()
            if self.have_ref_image == 1:
                self.button_register.Enable()
                self.button_subregion.Enable()
            else:
                self.button_register.Disable()
                self.button_subregion.Disable()
                self.button_delsubregion.Disable()
            
            if self.regist_calculated == 1:
                self.button_crop.Enable()
                self.button_accept.Enable()
                self.button_saveshifts.Enable()
                self.button_saveimg.Enable()
            else:
                self.button_crop.Disable()
                self.button_accept.Disable() 
                self.button_saveshifts.Disable()
                self.button_saveimg.Disable()
                           
        else:
            self.button_register.Disable()
            self.button_delsubregion.Disable()
            self.button_subregion.Disable()
            
            if self.have_ref_image == 1:
                self.button_manalign.Enable()
            else:
                self.button_manalign.Disable()
            
            if self.regist_calculated == 1:
                self.button_crop.Enable()
                self.button_accept.Enable()
                self.button_saveshifts.Enable()
            else:
                self.button_crop.Disable()
                self.button_accept.Disable() 
                self.button_saveshifts.Disable()
                           
            if self.man_align == 0:
                self.button_pick2ndpoint.Disable()
                self.button_applyman.Disable()
            elif self.man_align == 1:
                self.button_pick2ndpoint.Enable() 
            elif self.man_align == 2:
                self.button_applyman.Enable()       
               
        
        
        
#---------------------------------------------------------------------- 
class SpectralROI(wx.Frame):

    title = "Spectral Regions of Interest"

    def __init__(self, common, stack):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=self.title, size=(630, 720))
               
        ico = logos.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White") 
        
        self.stack = stack
        self.com = common
               
        self.fontsize = self.com.fontsize
        
        self.imin = 0
        self.imax = 0
        self.i0min = 0
        self.i0max = 0
        self.iselected = 0
        self.i0selected = 0        
        
                
        self.odtotal = self.stack.od3d.sum(axis=0)   
        self.odtotal = self.odtotal.sum(axis=0)/(self.stack.n_rows*self.stack.n_cols) 
        
        self.image_i0 = npy.zeros((self.stack.n_cols, self.stack.n_rows))
        self.image_i = npy.zeros((self.stack.n_cols, self.stack.n_rows))
        self.odthickmap = npy.zeros((self.stack.n_cols, self.stack.n_rows))
        
    
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
        text = wx.StaticText(panel, 0, 'First select I0 region below the edge, then select I region above the edge:')
        
        i1panel = wx.Panel(panel, -1, style = wx.SUNKEN_BORDER)
        self.SpectrumPanel = wxmpl.PlotPanel(i1panel, -1, size=(6.0, 3.7), cursor=False, crosshairs=False, location=False, zoom=False)
        
        wxmpl.EVT_SELECTION(i1panel, self.SpectrumPanel.GetId(), self.OnSelection)

        vbox.Add(text, 0, wx.TOP|wx.LEFT, 20)
        vbox.Add(i1panel, 0, wx.BOTTOM|wx.LEFT, 20)
        
       
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        
        
        
        sizer2 = wx.StaticBoxSizer(wx.StaticBox(panel, -1, 'Selected Spectral Regions'), orient=wx.VERTICAL)
        text = wx.StaticText(panel, 0, 'I Selection (red): ')
        self.textctrl1 = wx.TextCtrl(panel, -1, size = (328, 20), style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.textctrl1.SetValue('[  ]' )
        self.textctrl1.SetBackgroundColour("White")
        sizer2.Add((0,10))
        sizer2.Add(text, 1, wx.EXPAND|wx.LEFT, 10)
        sizer2.Add(self.textctrl1, 1, wx.EXPAND|wx.LEFT, 10)
        
        text = wx.StaticText(panel, 0, 'I0 Selection (green): ')
        self.textctrl2 = wx.TextCtrl(panel, -1, size = (328, 20), style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.textctrl2.SetValue('[  ]' )
        self.textctrl2.SetBackgroundColour("White")
        sizer2.Add((0,10))
        sizer2.Add(text, 1, wx.EXPAND|wx.LEFT, 10)
        sizer2.Add(self.textctrl2, 1, wx.EXPAND|wx.LEFT, 10)
        
        text = wx.StaticText(panel, 0, 'Optical density map')
        vbox.Add(text, 0, wx.LEFT|wx.RIGHT,20)
        vbox.Add((0,2))

        i2panel = wx.Panel(panel, -1, style = wx.SUNKEN_BORDER)
        self.ODMImagePanel = wxmpl.PlotPanel(i2panel, -1, size =(2.0,2.0), cursor=False, crosshairs=False, location=False, zoom=False)
        hbox2.Add(i2panel, 0, wx.LEFT|wx.RIGHT,20)
        hbox2.Add(sizer2, 1, wx.EXPAND|wx.LEFT|wx.RIGHT,20)
        
        vbox.Add(hbox2, 0, wx.EXPAND)      
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        
        button_save = wx.Button(panel, -1, 'Save')
        self.Bind(wx.EVT_BUTTON, self.OnSave, id=button_save.GetId())
        hbox.Add(button_save, 1, wx.ALL,20)
        
        button_cancel = wx.Button(panel, -1, 'Dismiss')
        self.Bind(wx.EVT_BUTTON, self.OnCancel, id=button_cancel.GetId())
        self.Bind(wx.EVT_CLOSE, self.OnCancel)
        hbox.Add(button_cancel, 1, wx.ALL ,20)
        
        vbox.Add(hbox, 0 )
        
        panel.SetSizer(vbox)
        
        vboxtop.Add(panel,1, wx.EXPAND )
        
        self.SetSizer(vboxtop)
        
        self.draw_spectrum()
        
        
#----------------------------------------------------------------------        
    def draw_spectrum(self):

        
        fig = self.SpectrumPanel.get_figure()
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        self.axes = fig.gca()
        
        mtplot.rcParams['font.size'] = self.fontsize

        specplot = self.axes.plot(self.stack.ev,self.odtotal)
        
        if self.i0selected == 1:
            self.axes.axvspan(self.stack.ev[self.i0min], self.stack.ev[self.i0max], facecolor='g', alpha=0.5)
            
        if self.iselected == 1:
            self.axes.axvspan(self.stack.ev[self.imin], self.stack.ev[self.imax], facecolor='r', alpha=0.5)

        
        self.axes.set_xlabel('Photon Energy [eV]')
        self.axes.set_ylabel('Optical Density')
        

        self.SpectrumPanel.draw()
        
        self.SetPosition((220, 150))
    
#----------------------------------------------------------------------        
    def draw_image(self):

               
        fig = self.ODMImagePanel.get_figure()
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        
        axes = fig.gca()
        divider = make_axes_locatable(axes)
        axcb = divider.new_horizontal(size="3%", pad=0.03)  

        fig.add_axes(axcb)
        
        axes.set_position([0.03,0.03,0.8,0.94])
        
        
        im = axes.imshow(self.odthickmap, cmap=mtplot.cm.get_cmap("gray")) 

        cbar = axes.figure.colorbar(im, orientation='vertical',cax=axcb) 

        #Show Scale Bar
        startx = int(self.stack.n_rows*0.05)
        starty = self.stack.n_cols-int(self.stack.n_cols*0.05)-self.stack.scale_bar_pixels_y
        um_string = ' $\mathrm{\mu m}$'
        microns = '$'+self.stack.scale_bar_string+' $'+um_string
        axes.text(self.stack.scale_bar_pixels_x+startx+1,starty+1, microns, horizontalalignment='left', verticalalignment='center',
                  color = 'black', fontsize=14)
        #Matplotlib has flipped scales so I'm using rows instead of cols!
        p = mtplot.patches.Rectangle((startx,starty), self.stack.scale_bar_pixels_x, self.stack.scale_bar_pixels_y,
                               color = 'black', fill = True)
        axes.add_patch(p)

         
        axes.axis("off")  
        self.ODMImagePanel.draw()


#----------------------------------------------------------------------        
    def OnSelection(self, evt):
        
        x1, y1 = evt.x1data, evt.y1data
        x2, y2 = evt.x2data, evt.y2data

        
        if (self.i0selected == 1) and (self.iselected ==1):
            self.i0selected = 0
            self.iselected = 0
        
        if self.i0selected == 0:       
            self.i0min = npy.abs(self.stack.ev - x1).argmin()
            self.i0max = npy.abs(self.stack.ev - x2).argmin()
            
            self.image_i0 = npy.sum(self.stack.absdata[:, :, self.i0min:self.i0max+1], axis=2)/(self.i0max+1-self.i0min)

            self.textctrl1.SetValue('Selection: [ '+str(self.stack.ev[self.i0min]) + ' eV, '+ str(self.stack.ev[self.i0max])+' eV ]' )
            self.i0selected = 1
            
        elif self.iselected == 0:
            self.imin = npy.abs(self.stack.ev - x1).argmin()
            self.imax = npy.abs(self.stack.ev - x2).argmin()
                       
            self.image_i = npy.sum(self.stack.absdata[:, :, self.imin:self.imax+1], axis=2)/(self.imax+1-self.imin)
            
            self.textctrl2.SetValue('Selection: [ '+str(self.stack.ev[self.imin]) + ' eV, '+ str(self.stack.ev[self.imax])+' eV ]' )
            self.iselected = 1        
            
        if (self.i0selected == 1) and (self.iselected ==1):              
            nonzeroind = self.image_i0.nonzero()
            self.odthickmap = npy.zeros((self.stack.n_cols, self.stack.n_rows))
            self.odthickmap[nonzeroind] = - npy.log(self.image_i[nonzeroind]/self.image_i0[nonzeroind])
            self.draw_image()
            
        
        self.draw_spectrum()

#----------------------------------------------------------------------        
    def OnSave(self, evt):
        #Save images
                       
        fileName = wx.FileSelector('Save OD Map', default_extension='png', 
                                   wildcard=('Portable Network Graphics (*.png)|*.png|' 
                                             + 'Adobe PDF Files (*.pdf)|*.pdf|All files (*.*)|*.*'), 
                                              parent=self, flags=wx.SAVE|wx.OVERWRITE_PROMPT) 
   
        if not fileName: 
            return 

        path, ext = os.path.splitext(fileName) 
        ext = ext[1:].lower() 
        
       
        if ext != 'png' and ext != 'pdf': 
            error_message = ( 
                  'Only the PNG and PDF image formats are supported.\n' 
                 'A file extension of `png\' or `pdf\' must be used.') 
            wx.MessageBox(error_message, 'Error - Could not save file.', 
                  parent=self, style=wx.OK|wx.ICON_ERROR) 
            return 
   
        try: 

            mtplot.rcParams['pdf.fonttype'] = 42
            
            fig = self.ODMImagePanel.get_figure()
            fig.savefig(fileName)

            
        except IOError, e:
            if e.strerror:
                err = e.strerror 
            else: 
                err = e 
   
            wx.MessageBox('Could not save file: %s' % err, 'Error', 
                          parent=self, style=wx.OK|wx.ICON_ERROR) 
            

                
#---------------------------------------------------------------------- 
    def OnCancel(self, evt):
        wx.GetApp().TopWindow.Show()
        self.Destroy()   

        
#---------------------------------------------------------------------- 
class DoseCalculation(wx.Frame):

    title = "Dose Calculation"

    def __init__(self, common, stack, ROIspectrum):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=self.title, size=(415, 270))
               
        ico = logos.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White") 
        
        self.stack = stack
        self.com = common
        self.ROIspectrum = ROIspectrum
               
        self.fontsize = self.com.fontsize
        
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        panel1 = wx.Panel(self, -1)
        
        gridtop = wx.FlexGridSizer(3, 2, vgap=20, hgap=20)
        #gridtop = wx.FlexGridSizer(4, 2, vgap=20, hgap=20)
        
        fontb = wx.SystemSettings_GetFont(wx.SYS_DEFAULT_GUI_FONT)
        fontb.SetWeight(wx.BOLD)
        
        
        st1 = wx.StaticText(panel1, -1, 'Detector efficiency [%]:',  style=wx.ALIGN_LEFT)
        st1.SetFont(fontb)
        st2 = wx.StaticText(panel1, -1, 'I region composition:',  style=wx.ALIGN_LEFT)
        st2.SetFont(fontb)
#        st3 = wx.StaticText(panel1, -1, 'Xray absorption length:',  style=wx.ALIGN_LEFT)
#        st3.SetFont(fontb)
        st4 = wx.StaticText(panel1, -1, 'Dose [Gray]:',  style=wx.ALIGN_LEFT)
        st4.SetFont(fontb)

        
        self.tc_1 = wx.TextCtrl(panel1, -1, size=((200,-1)), style=wx.TE_RICH|wx.VSCROLL, 
                                         value='30')
        
        self.tc_2 = wx.TextCtrl(panel1, -1, size=((200,-1)), style=wx.TE_RICH|wx.VSCROLL, 
                                         value='')
        
#        self.tc_3 = wx.TextCtrl(panel1, -1, size=((200,-1)), style=wx.TE_RICH|wx.VSCROLL|wx.TE_READONLY, 
#                                         value=' ')   
        
        self.tc_4 = wx.TextCtrl(panel1, -1, size=((200,-1)), style=wx.TE_RICH|wx.VSCROLL|wx.TE_READONLY, 
                                         value=' ')
        
      
        gridtop.Add(st1, 0)
        gridtop.Add( self.tc_1, 0)
        gridtop.Add(st2, 0)
        gridtop.Add( self.tc_2, 0)
#        gridtop.Add(st3, 0)
#        gridtop.Add( self.tc_3, 0)
        gridtop.Add(st4, 0)        
        gridtop.Add( self.tc_4, 0)             
          
        panel1.SetSizer(gridtop)

        button_calcdose = wx.Button(self, -1, 'Calculate Dose')
        self.Bind(wx.EVT_BUTTON, self.OnCalcDose, id=button_calcdose.GetId())
                
        button_cancel = wx.Button(self, -1, 'Dismiss')
        self.Bind(wx.EVT_BUTTON, self.OnCancel, id=button_cancel.GetId())


        vboxtop.Add(panel1, 1, wx.LEFT| wx.RIGHT|wx.TOP|wx.EXPAND, 20)
        vboxtop.Add(button_calcdose, 0, wx.LEFT| wx.RIGHT | wx.EXPAND, 20) 
        vboxtop.Add(button_cancel, 0, wx.LEFT| wx.RIGHT | wx.BOTTOM | wx.EXPAND, 20) 
        
        self.SetSizer(vboxtop)    
        
              
#---------------------------------------------------------------------- 
    def CalcDose(self):
                      
        
        try:
            detector_eff = 0.01*float(self.tc_1.GetValue())
        except:
            print 'Please enter numeric number for detector efficiency.'
            return
            
        
        i_composition = str(self.tc_2.GetValue())
        
        dose = 0.
        
        Chenke = henke.henke()
        
        #Check if composition array is recognizable
        try:
            z_array, atwt = Chenke.compound(i_composition,1.0)
        except:
            print 'Composition string error: Please re-enter composition string.'
            return
        
        
        dose = Chenke.dose_calc(self.stack, i_composition, self.ROIspectrum, self.stack.i0data, detector_eff)
        
        self.tc_4.SetValue(str(dose))
        
        return
        
        
#---------------------------------------------------------------------- 
    def OnCalcDose(self, evt):

        wx.BeginBusyCursor()
        self.CalcDose()
        wx.EndBusyCursor()
        
#---------------------------------------------------------------------- 
    def OnCancel(self, evt):

        self.Destroy()  
                
            
#---------------------------------------------------------------------- 
class PlotFrame(wx.Frame):

    def __init__(self, datax, datay, title = "I0 data"):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title = title, size=(630, 500))
        
        self.datax = datax
        self.datay = datay
        self.title = title
        
        ico = logos.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White")
        
        
        self.com = wx.GetApp().TopWindow.common         
        self.fontsize = self.com.fontsize
        
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
               
        i1panel = wx.Panel(panel, -1, style = wx.SUNKEN_BORDER)
        self.PlotPanel = wxmpl.PlotPanel(i1panel, -1, size=(6.0, 3.7), cursor=False, crosshairs=False, location=False, zoom=False)

        vbox.Add(i1panel, 0, wx.ALL, 20)
        

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        button_save = wx.Button(panel, -1, 'Save Spectrum')
        self.Bind(wx.EVT_BUTTON, self.OnSave, id=button_save.GetId())
        hbox.Add(button_save, 0, wx.ALL | wx.ALIGN_RIGHT ,20)        
                
        button_close = wx.Button(panel, -1, 'Close')
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=button_close.GetId())
        hbox.Add(button_close, 0, wx.ALL | wx.ALIGN_RIGHT ,20)
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        
        vbox.Add(hbox)
        

        panel.SetSizer(vbox)
        
        vboxtop.Add(panel,1, wx.EXPAND )
        
        self.SetSizer(vboxtop)
        
        #self.Centre()
        self.SetPosition((220, 150))

        
        self.draw_plot(datax,datay)
        
      
#----------------------------------------------------------------------        
    def draw_plot(self, datax, datay):
        
        
        fig = self.PlotPanel.get_figure()
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        self.axes = fig.gca()
        
        mtplot.rcParams['font.size'] = self.fontsize
        
        plot = self.axes.plot(datax,datay)
        
        self.axes.set_xlabel('Photon Energy [eV]')
        self.axes.set_ylabel('I0 Flux')
        

        
        self.PlotPanel.draw()
        
#----------------------------------------------------------------------
    def OnSave(self, event):


        try: 
            wildcard = "Csv files (*.csv)|*.csv"
            dialog = wx.FileDialog(None, "Save as .csv", wildcard=wildcard,
                                    style=wx.SAVE|wx.OVERWRITE_PROMPT)

            if dialog.ShowModal() == wx.ID_OK:
                filepath = dialog.GetPath()
                            
            wx.BeginBusyCursor()                  
            self.Save(filepath)    
            wx.EndBusyCursor()      

        except:

            wx.EndBusyCursor()
            wx.MessageBox("Could not save .csv file.")
                   
        dialog.Destroy()

        
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
            print>>f, '%.6f, %.6f' %(self.datax[ie], self.datay[ie])
        
        f.close()
    
        
#----------------------------------------------------------------------              
    def OnClose(self, evt):
        
        self.Destroy()

        
        
        
#---------------------------------------------------------------------- 
class ColorTableFrame(wx.Frame):

    def __init__(self):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title = "Pick Color Table", size=(200, 430))
        
        ico = logos.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White")
        
        
        self.com = wx.GetApp().TopWindow.common         
        self.fontsize = self.com.fontsize
        
        self.colors= ["gray","jet","autumn","bone", "cool","copper", "flag","hot","hsv","pink",
                      "prism","spring","summer","winter", "spectral"]
        
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
               
        sizer1 = wx.StaticBoxSizer(wx.StaticBox(panel, -1, 'Color Tables'),  orient=wx.VERTICAL)
        self.rb_grey = wx.RadioButton(panel, -1, self.colors[0], style=wx.RB_GROUP)
        self.rb_jet = wx.RadioButton(panel, -1, self.colors[1])
        self.rb_autumn = wx.RadioButton(panel, -1, self.colors[2])
        self.rb_bone = wx.RadioButton(panel, -1, self.colors[3])
        self.rb_cool = wx.RadioButton(panel, -1, self.colors[4])
        self.rb_copper = wx.RadioButton(panel, -1, self.colors[5])
        self.rb_flag = wx.RadioButton(panel, -1, self.colors[6])
        self.rb_hot = wx.RadioButton(panel, -1, self.colors[7])
        self.rb_hsv = wx.RadioButton(panel, -1, self.colors[8])
        self.rb_pink = wx.RadioButton(panel, -1, self.colors[9])
        self.rb_prism = wx.RadioButton(panel, -1, self.colors[10])
        self.rb_spring = wx.RadioButton(panel, -1, self.colors[11])
        self.rb_summer = wx.RadioButton(panel, -1, self.colors[12])
        self.rb_winter = wx.RadioButton(panel, -1, self.colors[13])
        self.rb_spectral = wx.RadioButton(panel, -1, self.colors[14])
        
        
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
        

        self.Bind(wx.EVT_RADIOBUTTON, self.OnColorTable, id=self.rb_grey.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.OnColorTable, id=self.rb_jet.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.OnColorTable, id=self.rb_autumn.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.OnColorTable, id=self.rb_bone.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.OnColorTable, id=self.rb_cool.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.OnColorTable, id=self.rb_copper.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.OnColorTable, id=self.rb_flag.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.OnColorTable, id=self.rb_hot.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.OnColorTable, id=self.rb_hsv.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.OnColorTable, id=self.rb_pink.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.OnColorTable, id=self.rb_prism.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.OnColorTable, id=self.rb_spring.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.OnColorTable, id=self.rb_summer.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.OnColorTable, id=self.rb_winter.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.OnColorTable, id=self.rb_spectral.GetId())
        
        sizer1.Add(self.rb_grey, 1, wx.EXPAND)
        sizer1.Add(self.rb_jet, 1, wx.EXPAND)
        sizer1.Add(self.rb_autumn, 1, wx.EXPAND)
        sizer1.Add(self.rb_bone, 1, wx.EXPAND)
        sizer1.Add(self.rb_cool, 1, wx.EXPAND)
        sizer1.Add(self.rb_copper, 1, wx.EXPAND)
        sizer1.Add(self.rb_flag, 1, wx.EXPAND)
        sizer1.Add(self.rb_hot, 1, wx.EXPAND)
        sizer1.Add(self.rb_hsv, 1, wx.EXPAND)
        sizer1.Add(self.rb_pink, 1, wx.EXPAND)
        sizer1.Add(self.rb_prism, 1, wx.EXPAND)
        sizer1.Add(self.rb_spring, 1, wx.EXPAND)
        sizer1.Add(self.rb_summer, 1, wx.EXPAND)
        sizer1.Add(self.rb_winter, 1, wx.EXPAND)
        sizer1.Add(self.rb_spectral, 1, wx.EXPAND)

        vbox.Add(sizer1, 1, wx.EXPAND | wx.ALL, 20)
        

        
        button_close = wx.Button(panel, -1, 'Close')
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=button_close.GetId())
        vbox.Add(button_close, 0, wx.RIGHT | wx.BOTTOM | wx.ALIGN_RIGHT ,20)
        

        panel.SetSizer(vbox)
        
        vboxtop.Add(panel,1, wx.EXPAND )
        
        self.SetSizer(vboxtop)
        
        self.ct_dict[ wx.GetApp().TopWindow.page1.colortable].SetValue(True)
        
        self.SetPosition((250, 150))
        
#----------------------------------------------------------------------          
    def OnColorTable(self, event):
        
        radioSelected = event.GetEventObject()

        wx.GetApp().TopWindow.page1.colortable = radioSelected.GetLabel()
        
        wx.GetApp().TopWindow.page1.loadImage()
        
#----------------------------------------------------------------------              
    def OnClose(self, evt):
        self.Destroy()
        
        
""" ------------------------------------------------------------------------------------------------"""
class PageLoadData(wx.Panel):
    def __init__(self, parent, common, data_struct, stack):
        wx.Panel.__init__(self, parent)
        
        self.data_struct = data_struct
        self.stk = stack
        self.com = common                  
        self.SetBackgroundColour("white")
        
        self.filename = " "
       
        self.fontsize = self.com.fontsize
        
        self.iev = 0

   
        
        #panel 1
        panel1 = wx.Panel(self, -1)
        sizer1 = wx.StaticBoxSizer(wx.StaticBox(panel1, -1, 'Load Data Stack'),orient=wx.VERTICAL)
        vbox1 = wx.BoxSizer(wx.VERTICAL)   

        self.button_hdf5 = wx.Button(panel1, -1, 'Load HDF5 Stack (*.hdf5)', size=((200,-1)))
        self.button_hdf5.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnLoadHDF5, id=self.button_hdf5.GetId())
        vbox1.Add(self.button_hdf5, 0, wx.EXPAND)
        
        self.button_sdf = wx.Button(panel1, -1, 'Load SDF Stack (*.hdr)')
        self.button_sdf.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnLoadSDF, id=self.button_sdf.GetId())   
        vbox1.Add(self.button_sdf, 0, wx.EXPAND)
        
        self.button_stk = wx.Button(panel1, -1, 'Load STK Stack (*.stk)')
        self.button_stk.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnLoadSTK, id=self.button_stk.GetId())   
        vbox1.Add(self.button_stk, 0, wx.EXPAND)
        
        self.button_xrm = wx.Button(panel1, -1, 'Load XRM Image (*.xrm)')
        self.Bind(wx.EVT_BUTTON, self.OnLoadXRM, id=self.button_xrm.GetId())
        self.button_xrm.SetFont(self.com.font)
        vbox1.Add(self.button_xrm, 0, wx.EXPAND)
        
        self.button_txrm = wx.Button(panel1, -1, 'Load TXRM Stack (*.txrm)')
        self.Bind(wx.EVT_BUTTON, self.OnLoadTXRM, id=self.button_txrm.GetId())
        self.button_txrm.SetFont(self.com.font)
        vbox1.Add(self.button_txrm, 0, wx.EXPAND)    
        
        self.button_tif = wx.Button(panel1, -1, 'Load TIF Stack (*.tif)')
        self.Bind(wx.EVT_BUTTON, self.OnLoadTIF, id=self.button_tif.GetId())
        self.button_tif.SetFont(self.com.font)
        vbox1.Add(self.button_tif, 0, wx.EXPAND)      

        sizer1.Add(vbox1,1, wx.ALL|wx.EXPAND,10)
        panel1.SetSizer(sizer1)
        
        
        #panel 2
        panel2 = wx.Panel(self, -1)
        sizer2 = wx.StaticBoxSizer(wx.StaticBox(panel2, -1, 'Build a stack from a set of files'),orient=wx.VERTICAL)
        vbox2 = wx.BoxSizer(wx.VERTICAL)


        self.button_sm = wx.Button(panel2, -1, 'Build a stack from a set of NetCDF (*.sm) files')
        self.button_sm.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnBuildStack, id=self.button_sm.GetId())
        vbox2.Add(self.button_sm, 0, wx.EXPAND)
        
        self.button_xrm_list = wx.Button(panel2, -1, 'Build a stack from a set of XRM (*.xrm) files')
        self.button_xrm_list.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnBuildStack, id=self.button_xrm_list.GetId())
        vbox2.Add(self.button_xrm_list, 0, wx.EXPAND)
        
        sizer2.Add(vbox2,1, wx.ALL|wx.EXPAND,10)
        panel2.SetSizer(sizer2)
        

        #panel 3
        panel3 = wx.Panel(self, -1)
        sizer3 = wx.StaticBoxSizer(wx.StaticBox(panel3, -1, 'File'),orient=wx.VERTICAL)
        vbox3 = wx.BoxSizer(wx.VERTICAL)
        vbox3.Add((0,3))      
  
        self.tc_file = wx.TextCtrl(panel3, -1, size=((700,30)), style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        vbox3.Add(self.tc_file, 0, wx.EXPAND)

        
        sizer3.Add(vbox3,1, wx.ALL|wx.EXPAND,10)
        panel3.SetSizer(sizer3)

        
        #panel 4
        panel4 = wx.Panel(self, -1)
        sizer4 = wx.StaticBoxSizer(wx.StaticBox(panel4, -1, 'Path'),orient=wx.VERTICAL)
        vbox4 = wx.BoxSizer(wx.VERTICAL)
        vbox4.Add((0,3))      
  
        self.tc_path = wx.TextCtrl(panel4, -1, size=((700,30)), style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        vbox4.Add(self.tc_path, 0, wx.EXPAND)

        
        sizer4.Add(vbox4,1, wx.ALL|wx.EXPAND,10)
        panel4.SetSizer(sizer4)
        
        #panel 5     
        panel5 = wx.Panel(self, -1)
        vbox5 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_imageeng = wx.TextCtrl(panel5, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE, size=(300,-1))
        self.tc_imageeng.SetFont(self.com.font)
        self.tc_imageeng.SetValue("Image at energy: ")
        
        hbox51 = wx.BoxSizer(wx.HORIZONTAL)
   
        i1panel = wx.Panel(panel5, -1, style = wx.SUNKEN_BORDER)
        self.AbsImagePanel = wxmpl.PlotPanel(i1panel, -1, size =(PlotH*.8, PlotH*.8), cursor=False, crosshairs=False, location=False, zoom=False)
        
        vbox51 = wx.BoxSizer(wx.VERTICAL)
        self.slider_eng = wx.Slider(panel5, 0, self.iev, 0, 100, style=wx.SL_LEFT|wx.SL_VERTICAL)        
        self.slider_eng.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnScrollEng, self.slider_eng)

        self.engspin = wx.SpinButton(panel5, -1, size = ((8,-1)), style=wx.SP_ARROW_KEYS)
        self.Bind(wx.EVT_SPIN_UP, self.OnEngspinUp, self.engspin)
        self.Bind(wx.EVT_SPIN_DOWN, self.OnEngspinDown, self.engspin)
        

        hbox51.Add(i1panel, 0)
        vbox51.Add((0,3))
        vbox51.Add(self.slider_eng, 1,  wx.EXPAND) 
        vbox51.Add(self.engspin, 0,  wx.EXPAND)         
        hbox51.Add(vbox51, 0,  wx.EXPAND) 
        
        vbox5.Add(self.tc_imageeng, 0)        
        vbox5.Add(hbox51, 0)

        panel5.SetSizer(vbox5)
        
        
        
        hboxT = wx.BoxSizer(wx.HORIZONTAL)
        vboxT1 = wx.BoxSizer(wx.VERTICAL)
        vboxT2 = wx.BoxSizer(wx.VERTICAL)

        hboxT1 = wx.BoxSizer(wx.HORIZONTAL)
        hboxT2 = wx.BoxSizer(wx.HORIZONTAL)
        hboxT3 = wx.BoxSizer(wx.HORIZONTAL)
        hboxT4 = wx.BoxSizer(wx.HORIZONTAL)                  
          
        hboxT1.Add((50,0))
        hboxT1.Add(panel1, 0, wx.BOTTOM | wx.TOP, 10)
              
        hboxT4.Add((50,0)) 
        hboxT4.Add(panel2, 0, wx.BOTTOM | wx.TOP, 10) 


        vboxT1.Add(hboxT1, 0, wx.ALL, 5)
        vboxT1.Add(hboxT4, 0, wx.ALL, 5)
        
        hboxT.Add(vboxT1, 0,  wx.ALL, 5)
        hboxT.Add((100,0)) 
        hboxT.Add(panel5, 0, wx.ALL, 5) 
        
        hboxT2.Add((50,0))  
        hboxT2.Add(panel3, 1, wx.TOP | wx.EXPAND, 10)  
        hboxT3.Add((50,0))
        hboxT3.Add(panel4, 1, wx.TOP | wx.EXPAND, 10)
                
        vboxT2.Add((0,30))
        vboxT2.Add(hboxT, 0, wx.ALL, 5)
        vboxT2.Add(hboxT2, 0, wx.ALL, 5)
        vboxT2.Add(hboxT3, 0,  wx.ALL, 5)
        

        self.SetSizer(vboxT2) 
        
        
#----------------------------------------------------------------------          
    def OnLoadHDF5(self, event):

        wildcard =  "HDF5 files (*.hdf5)|*.hdf5|SDF files (*.hdr)|*.hdr|STK files (*.stk)|*.stk|TXRM (*.txrm)|*.txrm|XRM (*.xrm)|*.xrm|TIF (*.tif)|*.tif" 
        wx.GetApp().TopWindow.LoadStack(wildcard)
        
#----------------------------------------------------------------------          
    def OnLoadSDF(self, event):

        wildcard =  "SDF files (*.hdr)|*.hdr|HDF5 files (*.hdf5)|*.hdf5|STK files (*.stk)|*.stk|TXRM (*.txrm)|*.txrm|XRM (*.xrm)|*.xrm|TIF (*.tif)|*.tif" 
        wx.GetApp().TopWindow.LoadStack(wildcard)
                
#----------------------------------------------------------------------          
    def OnLoadSTK(self, event):

        wildcard =  "STK files (*.stk)|*.stk|HDF5 files (*.hdf5)|*.hdf5|SDF files (*.hdr)|*.hdr|TXRM (*.txrm)|*.txrm|XRM (*.xrm)|*.xrm|TIF (*.tif)|*.tif" 
        wx.GetApp().TopWindow.LoadStack(wildcard)
        
#----------------------------------------------------------------------          
    def OnLoadTXRM(self, event):

        wildcard =  "TXRM (*.txrm)|*.txrm|HDF5 files (*.hdf5)|*.hdf5|SDF files (*.hdr)|*.hdr|STK files (*.stk)|*.stk|XRM (*.xrm)|*.xrm|TIF (*.tif)|*.tif" 
        wx.GetApp().TopWindow.LoadStack(wildcard)

#----------------------------------------------------------------------          
    def OnLoadXRM(self, event):

        wildcard =  "XRM (*.xrm)|*.xrm|HDF5 files (*.hdf5)|*.hdf5|SDF files (*.hdr)|*.hdr|STK files (*.stk)|*.stk|TXRM (*.txrm)|*.txrm|TIF (*.tif)|*.tif" 
        wx.GetApp().TopWindow.LoadStack(wildcard)

#----------------------------------------------------------------------          
    def OnLoadTIF(self, event):

        wildcard =  "TIF (*.tif)|*.tif|XRM (*.xrm)|*.xrm|HDF5 files (*.hdf5)|*.hdf5|SDF files (*.hdr)|*.hdr|STK files (*.stk)|*.stk|TXRM (*.txrm)|*.txrm" 
        wx.GetApp().TopWindow.LoadStack(wildcard)
                
#----------------------------------------------------------------------          
    def OnBuildStack(self, event):

        wx.GetApp().TopWindow.BuildStack()
        
#----------------------------------------------------------------------            
    def OnScrollEng(self, event):
        self.iev = event.GetInt()

        if self.com.stack_loaded == 1:
            self.ShowImage()
            
#----------------------------------------------------------------------            
    def OnEngspinUp(self, event):
        if (self.com.stack_loaded == 1) and (self.iev > 0):
            self.iev = self.iev - 1
            self.slider_eng.SetValue(self.iev)

            self.ShowImage()

            
#----------------------------------------------------------------------            
    def OnEngspinDown(self, event):
        if (self.com.stack_loaded == 1) and (self.iev < self.stk.n_ev-1):
            self.iev = self.iev + 1
            self.slider_eng.SetValue(self.iev)

            self.ShowImage()
            
#----------------------------------------------------------------------        
    def ShowImage(self):
        
        image = self.stk.absdata[:,:,int(self.iev)].copy() 

        fig = self.AbsImagePanel.get_figure()
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        axes = fig.gca()
        fig.patch.set_alpha(1.0)
        
        im = axes.imshow(image, cmap=mtplot.cm.get_cmap("gray")) 
        
        if wx.GetApp().TopWindow.page1.show_scale_bar == 1:
            #Show Scale Bar
            startx = int(self.stk.n_rows*0.05)
            starty = self.stk.n_cols-int(self.stk.n_cols*0.05)-self.stk.scale_bar_pixels_y
            um_string = ' $\mathrm{\mu m}$'
            microns = '$'+self.stk.scale_bar_string+' $'+um_string
            axes.text(self.stk.scale_bar_pixels_x+startx+1,starty+1, microns, horizontalalignment='left', verticalalignment='center',
                      color = 'black', fontsize=14)
            #Matplotlib has flipped scales so I'm using rows instead of cols!
            p = mtplot.patches.Rectangle((startx,starty), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                                   color = 'black', fill = True)
            axes.add_patch(p)
            
       
        axes.axis("off")      
        self.AbsImagePanel.draw()
        
        self.tc_imageeng.SetValue('Image at energy: {0:5.2f} eV'.format(float(self.stk.ev[self.iev])))
        

        
#----------------------------------------------------------------------        
    def ShowInfo(self, filename, filepath):
        
        self.ShowImage()
        
        self.tc_file.SetValue(filename)
        self.tc_path.SetValue(filepath)
        
        
            
""" ------------------------------------------------------------------------------------------------"""
class MainFrame(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, size=(Winsizex, Winsizey))
        
        
      
        ico = logos.getlogo_2l_32Icon()
        self.SetIcon(ico)

        
        self.initToolbar()
        
        self.data_struct = data_struct.h5()
        self.stk = data_stack.data(self.data_struct)
        self.anlz = analyze.analyze(self.stk)
#        self.nnma = nnma.nnma(self.stk, self.data_struct, self.anlz)
        self.common = common()
        
        self.SetFont(self.common.font)
              

        # Here we create a panel and a notebook on the panel
        p = wx.Panel(self)
        nb = wx.Notebook(p, style=wx.BORDER_STATIC)

        # create the page windows as children of the notebook
        self.page0 = PageLoadData(nb, self.common, self.data_struct, self.stk)
        self.page1 = PageStack(nb, self.common, self.data_struct, self.stk)
        self.page2 = PagePCA(nb, self.common, self.data_struct, self.stk, self.anlz)
        self.page3 = PageCluster(nb, self.common, self.data_struct, self.stk, self.anlz)
        self.page4 = PageSpectral(nb, self.common, self.data_struct, self.stk, self.anlz)
        self.page7 = None


        # add the pages to the notebook with the label to show on the tab
        nb.AddPage(self.page0, "Load Data")
        nb.AddPage(self.page1, "Preprocess Data")
        nb.AddPage(self.page2, "PCA")
        nb.AddPage(self.page3, "Cluster Analysis")
        nb.AddPage(self.page4, "Spectral Analysis")
        # Only add NNMA pages if option "--nnma" is given in command line
        options, extraParams = getopt.getopt(sys.argv[1:], '', ['wx', 'batch', 'nnma', 'ica', 'keyeng'])
        for opt, arg in options:
            if opt in '--nnma':
                if verbose: print "Running with NNMA."
                self.page5Notebook = wx.Notebook(nb)    # page 5 is a wx.Notebook, with different pages for different NNMA analysis results
                #self.page5 = PageNNMA(nb, self.common, self.data_struct, self.stk, self.nnma)
                self.page5a = PageNNMAOptDensity(self.page5Notebook, self.common, self.data_struct, self.stk, self.nnma)
                self.page5b = PageNNMASpectra(self.page5Notebook, self.common, self.data_struct, self.stk, self.nnma)
                self.page5c = PageNNMAThickness(self.page5Notebook, self.common, self.data_struct, self.stk, self.nnma)
                nb.AddPage(self.page5Notebook, "NNMA Analysis")
                self.page5Notebook.AddPage(self.page5a, 'NNMA optical density')
                self.page5Notebook.AddPage(self.page5b, 'NNMA spectra')
                self.page5Notebook.AddPage(self.page5c, 'NNMA thickness maps')

            if opt in '--ica':
                if verbose: print "Running with ICA."
                self.page6 = PageICA(nb, self.common, self.data_struct, self.stk, self.anlz)
                nb.AddPage(self.page6, "ICA")
                
            
#             if opt in '--keyeng':
#                 if verbose: print "Running with KeyEng."
#                 self.page7 = PageKeyEng(nb, self.common, self.data_struct, self.stk, self.anlz)
#                 nb.AddPage(self.page7, "Key Energies")

          

        self.page7 = PageKeyEng(nb, self.common, self.data_struct, self.stk, self.anlz)
        nb.AddPage(self.page7, "Key Energies")
                            
                
                
        # finally, put the notebook in a sizer for the panel to manage
        # the layout
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(nb, 1, wx.EXPAND)
        p.SetSizer(sizer)
  
        self.Centre()
        self.Show(True) 
        
        
        
#----------------------------------------------------------------------        
    def initToolbar(self):

        self.toolbar = self.CreateToolBar(style= wx.TB_HORIZONTAL | wx.NO_BORDER | wx.TB_FLAT | wx.TB_TEXT)
        self.toolbar.SetToolBitmapSize((16,16))
         
        open_ico = wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN, wx.ART_TOOLBAR, (16,16))
        openTool = self.toolbar.AddSimpleTool(wx.ID_OPEN, open_ico, "Open stack .hdf5 or .stk", "Open stack .hdf5 or .stk")
        self.Bind(wx.EVT_MENU, self.onBrowse, openTool)
        
        open_sl_ico = wx.Image(resource_path(os.path.join('images', 'open-sl.png')), wx.BITMAP_TYPE_ANY).ConvertToBitmap()
        #open_sl_ico = logos.getslBitmap()
        openslTool = self.toolbar.AddSimpleTool(101, open_sl_ico, "Open directory with file sequence", "Open direcory with file sequence")   
        self.Bind(wx.EVT_MENU, self.onOpenSL, openslTool)
        
        save_ico = wx.ArtProvider.GetBitmap(wx.ART_FILE_SAVE, wx.ART_TOOLBAR, (16,16))
        saveTool = self.toolbar.AddSimpleTool(wx.ID_SAVE, save_ico, "Save results to .hdf5", "Save results to .hdf5")
        self.toolbar.EnableTool(wx.ID_SAVE, False)       
        self.Bind(wx.EVT_MENU, self.onSaveResultsToH5, saveTool)     
        
        saveas_ico = wx.ArtProvider.GetBitmap(wx.ART_FILE_SAVE_AS, wx.ART_TOOLBAR, (16,16))
        saveasTool = self.toolbar.AddSimpleTool(wx.ID_SAVEAS, saveas_ico, "Save as .hdf5", "Save as .hdf5")
        self.Bind(wx.EVT_MENU, self.onSaveAsH5, saveasTool)
        self.toolbar.EnableTool(wx.ID_SAVEAS, False)
        
        help_ico = wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, wx.ART_TOOLBAR, (16,16))
        helpTool = self.toolbar.AddSimpleTool(wx.ID_HELP, help_ico, "About", "About")
        self.Bind(wx.EVT_MENU, self.onAbout, helpTool)
    
               
 
        self.toolbar.Realize()


#----------------------------------------------------------------------
    def onBrowse(self, event):
        
        self.LoadStack()
 
        
#----------------------------------------------------------------------
    def LoadStack(self, wildcard = ''):
        """
        Browse for a stack file:
        """

        try:
            if wildcard == '':
                wildcard =  "HDF5 files (*.hdf5)|*.hdf5|SDF files (*.hdr)|*.hdr|STK files (*.stk)|*.stk|TXRM (*.txrm)|*.txrm|XRM (*.xrm)|*.xrm|TIF (*.tif)|*.tif" 
            dialog = wx.FileDialog(None, "Choose a file", style=wx.OPEN)
            
            dialog.SetWildcard(wildcard)
            if dialog.ShowModal() == wx.ID_OK:
                filepath = dialog.GetPath()
                self.page1.filename = dialog.GetFilename()
                directory = dialog.GetDirectory()
                defaultDir = directory
            else:
                return
            
            wx.BeginBusyCursor() 
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
            self.page0.slider_eng.SetRange(0,self.stk.n_ev-1)
            self.page0.iev = self.iev
            self.page0.slider_eng.SetValue(self.iev)
                     
            self.page1.slider_eng.SetRange(0,self.stk.n_ev-1)
            self.page1.iev = self.iev
            self.page1.slider_eng.SetValue(self.iev)
            
        
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
            self.page1.textctrl.SetValue(self.page1.filename)
            
            self.page0.ShowInfo(self.page1.filename, directory)
            

            wx.EndBusyCursor()
                
        except:

            self.common.stack_loaded = 0 
            self.common.i0_loaded = 0
            self.new_stack_refresh()
                               
            wx.EndBusyCursor()
            
            wx.MessageBox("Image stack not loaded.")
            import sys; print sys.exc_info()
                   
        dialog.Destroy()
        self.refresh_widgets()
        
#----------------------------------------------------------------------
    def onOpenSL(self, event):

        self.BuildStack()
                
#----------------------------------------------------------------------
    def BuildStack(self):
        """
        Browse for .sm files
        """

        try:
            dialog = wx.DirDialog(None, "Choose a directory",
                                   style=wx.DD_DIR_MUST_EXIST|wx.DD_CHANGE_DIR)
            #dialog.SetWildcard(wildcard)
            if dialog.ShowModal() == wx.ID_OK:
                directory = dialog.GetPath()
                
            self.common.path = directory

            StackListFrame(directory, self.common, self.stk, self.data_struct).Show()
            
        except:
            print 'Error could not build stack list.'
            self.common.stack_loaded = 0 
            self.common.i0_loaded = 0
            self.new_stack_refresh()
            self.refresh_widgets()
                               
            wx.MessageBox(".sm files not loaded.")
            import sys; print sys.exc_info()
        

#----------------------------------------------------------------------
    def onSaveAsH5(self, event):
        self.SaveStackH5()
        
        
#----------------------------------------------------------------------
    def SaveStackH5(self):

        """
        Browse for .hdf5 file
        """

        try: 
            wildcard = "HDF5 files (*.hdf5)|*.hdf5"
            dialog = wx.FileDialog(None, "Save as .hdf5", wildcard=wildcard,
                                    style=wx.SAVE|wx.OVERWRITE_PROMPT)

            if dialog.ShowModal() == wx.ID_OK:
                filepath = dialog.GetPath()
                self.page1.filename = dialog.GetFilename()
                dir = dialog.GetDirectory()
                
                self.common.path = dir
                self.common.filename = self.page1.filename

            wx.BeginBusyCursor()                  
            self.stk.write_h5(filepath, self.data_struct)    
            wx.EndBusyCursor()      

        except:

            wx.EndBusyCursor()
            wx.MessageBox("Could not save HDF5 file.")
                   
        dialog.Destroy()
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
            wildcard = "HDF5 files (*.hdf5)|*.hdf5"
            dialog = wx.FileDialog(None, "Save as .hdf5", wildcard=wildcard,
                                    style=wx.SAVE|wx.OVERWRITE_PROMPT)

            if dialog.ShowModal() == wx.ID_OK:
                filepath = dialog.GetPath()
                self.page1.filename = dialog.GetFilename()
                dir = dialog.GetDirectory()
                
                self.common.path = dir
                self.common.filename = self.page1.filename

            wx.BeginBusyCursor()                  
            self.stk.write_results_h5(filepath, self.data_struct, self.anlz)    
            wx.EndBusyCursor()      

        except:

            wx.EndBusyCursor()
            wx.MessageBox("Could not save HDF5 file.")
                   
        dialog.Destroy()
        self.refresh_widgets()
        
        return
          
#----------------------------------------------------------------------
    def onAbout(self, event):
        AboutFrame().Show()
        return
        
#----------------------------------------------------------------------        
    def refresh_widgets(self):
        

        if self.common.stack_loaded == 0:
            self.page1.button_i0ffile.Disable()
            self.page1.button_i0histogram.Disable() 
            self.page1.button_save.Disable() 
            self.page1.button_savestack.Disable()
            self.page1.button_align.Disable()
            self.page1.button_slideshow.Disable()
            self.page1.button_addROI.Disable()
            self.page1.button_spectralROI.Disable()
            self.page1.button_resetdisplay.Disable() 
            self.page1.button_despike.Disable()   
            self.page1.button_displaycolor.Disable()
            self.toolbar.EnableTool(wx.ID_SAVE, False)
            self.toolbar.EnableTool(wx.ID_SAVEAS, False)
        else:
            self.page1.button_i0ffile.Enable()
            self.page1.button_i0histogram.Enable() 
            self.page1.button_save.Enable() 
            self.page1.button_savestack.Enable()  
            self.page1.button_align.Enable()  
            self.page1.button_slideshow.Enable()
            self.page1.button_addROI.Enable()  
            self.page1.button_spectralROI.Enable()
            self.page1.button_resetdisplay.Enable() 
            self.page1.button_despike.Enable() 
            self.page1.button_displaycolor.Enable()
            self.toolbar.EnableTool(wx.ID_SAVE, True) 
            self.toolbar.EnableTool(wx.ID_SAVEAS, True)  
            
            
        if self.common.i0_loaded == 0:
            self.page1.button_limitev.Disable()
            self.page1.button_subregion.Disable()
            self.page1.button_showi0.Disable() 
            self.page1.rb_flux.Disable()
            self.page1.rb_od.Disable()
            self.page2.button_calcpca.Disable()
            self.page4.button_loadtspec.Disable()
            self.page4.button_addflat.Disable()
        else:
            self.page1.button_limitev.Enable()
            self.page1.button_subregion.Enable()
            self.page1.button_showi0.Enable()
            self.page1.rb_flux.Enable()
            self.page1.rb_od.Enable()   
            self.page2.button_calcpca.Enable() 
            self.page4.button_loadtspec.Enable()
            self.page4.button_addflat.Enable()   
            
            
            
        if self.common.pca_calculated == 0:      
            self.page2.button_savepca.Disable()
            self.page2.slidershow.Disable() 
            self.page3.button_calcca.Disable()
            self.page4.rb_fit.Disable()
            if self.page7: self.page7.button_calckeng.Disable()
        else:
            self.page2.button_savepca.Enable()
            self.page2.slidershow.Enable()
            self.page3.button_calcca.Enable()  
            self.page4.rb_fit.Enable()  
            if self.page7: self.page7.button_calckeng.Enable()       
            
        if self.common.cluster_calculated == 0:   
            self.page3.button_scatterplots.Disable()
            self.page3.button_savecluster.Disable()
            self.page3.slidershow.Disable()
            self.page4.button_addclspec.Disable()
        else:
            self.page3.button_scatterplots.Enable()
            self.page3.button_savecluster.Enable()  
            self.page3.slidershow.Enable()
            self.page4.button_addclspec.Enable()
            
        if self.common.spec_anl_calculated == 0:
            self.page4.button_removespec.Disable()
            self.page4.button_movespdown.Disable()
            self.page4.button_movespup.Disable()
            self.page4.button_save.Disable()
            self.page4.button_showrgb.Disable() 
        else:
            self.page4.button_removespec.Enable()
            self.page4.button_movespdown.Enable()
            self.page4.button_movespup.Enable()
            self.page4.button_save.Enable() 
            self.page4.button_showrgb.Enable()          
                  
            
        self.page1.ResetDisplaySettings()
            
            
            
#----------------------------------------------------------------------        
    def new_stack_refresh(self):
        
        
        self.common.i0_loaded = 0
        self.common.pca_calculated = 0
        self.common.cluster_calculated = 0
        
        self.refresh_widgets()
        
        #page 1
        self.page1.rb_flux.SetValue(True)
        self.page1.rb_od.SetValue(False)
        self.page1.showflux = True
        
        fig = self.page1.SpectrumPanel.get_figure()
        fig.clf()
        self.page1.SpectrumPanel.draw()
        self.page1.tc_spec.SetValue("Spectrum at point: ")       
        
        fig = self.page1.AbsImagePanel.get_figure()
        fig.clf()
        self.page1.AbsImagePanel.draw()        
        self.page1.tc_imageeng.SetValue("Image at energy: ")
        
        self.page1.textctrl.SetValue(' ')
        
        self.page1.ResetDisplaySettings()
        
        
        #page 2
        fig = self.page2.PCAEvalsPan.get_figure()
        fig.clf()
        self.page2.PCAEvalsPan.draw()
        
        fig = self.page2.PCAImagePan.get_figure()
        fig.clf()
        self.page2.PCAImagePan.draw()
        
        fig = self.page2.PCASpecPan.get_figure()
        fig.clf()
        self.page2.PCASpecPan.draw()
        
        self.page2.vartc.SetLabel('0%')
        self.page2.npcaspin.SetValue(1)
        self.page2.tc_PCAcomp.SetValue("PCA component ")
        self.page2.text_pcaspec.SetValue("PCA spectrum ")
        
        self.page2.selpca = 1       
        self.page2.numsigpca = 2
        self.page2.slidershow.SetValue(self.page2.selpca)
        
        #page 3
        fig = self.page3.ClusterImagePan.get_figure()
        fig.clf()
        self.page3.ClusterImagePan.draw()
        
        fig = self.page3.ClusterIndvImagePan.get_figure()
        fig.clf()
        self.page3.ClusterIndvImagePan.draw()
        
        fig = self.page3.ClusterSpecPan.get_figure()
        fig.clf()
        self.page3.ClusterSpecPan.draw()
        
        fig = self.page3.ClusterDistMapPan.get_figure()
        fig.clf()
        self.page3.ClusterDistMapPan.draw()       
       
        self.page3.selcluster = 1
        self.page3.slidershow.SetValue(self.page3.selcluster)
        self.page3.numclusters = 5
        self.page3.nclusterspin.SetValue(self.page3.numclusters)
        self.page3.tc_cluster.SetValue("Cluster ")
        self.page3.tc_clustersp.SetValue("Cluster spectrum")
        self.page3.wo_1st_pca = 0
        self.page3.remove1stpcacb.SetValue(False)

        
        #page 4
        self.page4.ClearWidgets()
        
        #page 7
        if self.page7:
            self.page7.button_calckeng.Disable()
            fig = self.page7.KESpecPan.get_figure()
            fig.clf()
            self.page7.KESpecPan.draw()         
            fig = self.page7.AbsImagePanel.get_figure()
            fig.clf()
            self.page7.AbsImagePanel.draw()   
            self.page7.lc_1.DeleteAllItems()   
            self.page7.keyenergies = []
            self.page7.keyengs_calculated = 0     
        
#---------------------------------------------------------------------- 
class StackListFrame(wx.Frame):

    def __init__(self, filepath, com, stack, data_struct):
        
        self.data_struct = data_struct
        self.stk = stack
        self.common = com
        
           
            
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title = "Stack File List", size=(410, 510))
            
        ico = logos.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.filepath = filepath
        
        self.have1st = 0
        self.havelast = 0
        self.file1st = ' '
        self.filelast = ' '
        
        self.filetype = ''
        
        
        self.com = wx.GetApp().TopWindow.common         
        self.fontsize = self.com.fontsize
        
        self.SetBackgroundColour("White")
        
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        
        panel1 = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
        
        
        self.textt = wx.TextCtrl(panel1,-1, 'Select first stack file', 
                                 style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.textt.SetFont(self.com.font)
        
        
        vbox.Add(self.textt, 0, wx.EXPAND)

        
        self.filelist = wx.wx.ListCtrl(panel1,-1,style=wx.LC_REPORT)
        self.filelist.InsertColumn(0,"File list", width=150)
        self.filelist.InsertColumn(1,"X", width=50)
        self.filelist.InsertColumn(2,"Y", width=50)
        self.filelist.InsertColumn(3,"eV", width=70)
        
        self.filelist.SetFont(self.com.font)
        self.filelist.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnFileList)
        vbox.Add(self.filelist, 1, wx.EXPAND)
        
        
        self.tc_first = wx.TextCtrl(panel1,-1, 'First stack file: ', size = (300,20),
                                    style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_last = wx.TextCtrl(panel1,-1, 'Last stack file: ', size = (300,20),
                                   style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_first.SetFont(self.com.font)
        self.tc_last.SetFont(self.com.font)
        
        vbox.Add(self.tc_first,0, wx.TOP, 10)
        vbox.Add(self.tc_last, 0, wx.TOP, 10)
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        
        self.button_accept = wx.Button(panel1, -1, 'Accept')
        self.button_accept.Disable()
        self.button_accept.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnAccept, id=self.button_accept.GetId())
        hbox.Add(self.button_accept, 0, wx.ALL|wx.ALIGN_RIGHT, 10)
        
        button_cancel = wx.Button(panel1, -1, 'Cancel')
        self.Bind(wx.EVT_BUTTON, self.OnCancel, id=button_cancel.GetId())
        button_cancel.SetFont(self.com.font)
        hbox.Add(button_cancel, 0, wx.ALL|wx.ALIGN_RIGHT,10)
        
        vbox.Add(hbox, 0, wx.ALIGN_RIGHT )
      
        panel1.SetSizer(vbox)
        
        vboxtop.Add(panel1, 1, wx.ALL|wx.EXPAND,20)
        
        self.SetSizer(vboxtop)
        self.Centre()
        
        self.ShowFileList()
        
        
#----------------------------------------------------------------------           
    def ShowFileList(self):    
        
                
        self.sm_files = [x for x in os.listdir(self.filepath) if x.endswith('.sm')]
        
        if self.sm_files:
            
            self.filetype = 'sm'
            
            try: 
                from netCDF4 import Dataset
                import file_sm_netcdf
                self.sm = file_sm_netcdf.sm(data_struct)
            
            except:
                wx.MessageBox("Could not import netCDF4 library.")
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
            
            import file_xrm
            self.xrm = file_xrm.xrm()

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
    def OnFileList(self, evt):
        
        if (self.have1st == 1) and (self.havelast==1):
            self.have1st = 0
            self.havelast = 0
            self.button_accept.Disable()
            self.tc_first.SetValue('First stack file: ')
            self.tc_last.SetValue('Last stack file: ')
            
        
        fn =  evt.GetText()
        if self.have1st == 0:
            self.tc_first.SetValue('First stack file: ' + fn) 
            self.textt.SetValue('Select last stack file')
            self.file1st = fn
            self.have1st = 1
        elif self.havelast == 0:
            self.tc_last.SetValue('Last stack file: ' + fn) 
            self.textt.SetValue('Select first stack file')         
            self.filelast = fn
            self.button_accept.Enable()
            self.havelast = 1

        
        
#----------------------------------------------------------------------        
    def OnAccept(self, evt):
        
        wx.BeginBusyCursor()
        
        wx.GetApp().TopWindow.new_stack_refresh()
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
                   

        wx.GetApp().TopWindow.page1.iev = int(self.stk.n_ev/3)
   
        wx.GetApp().TopWindow.ix = int(self.stk.n_cols/2)
        wx.GetApp().TopWindow.iy = int(self.stk.n_rows/2)
                        
        self.common.stack_loaded = 1
        
        wx.GetApp().TopWindow.refresh_widgets()
        wx.GetApp().TopWindow.page1.ResetDisplaySettings()
        wx.GetApp().TopWindow.page1.filename = filelist[0]
        wx.GetApp().TopWindow.page1.textctrl.SetValue(filelist[0])
 
        wx.GetApp().TopWindow.page0.slider_eng.SetRange(0,self.stk.n_ev-1)
        wx.GetApp().TopWindow.page0.iev = self.stk.n_ev/2
        wx.GetApp().TopWindow.page0.slider_eng.SetValue(wx.GetApp().TopWindow.page1.iev)       
        
        wx.GetApp().TopWindow.page1.slider_eng.SetRange(0,self.stk.n_ev-1)
        wx.GetApp().TopWindow.page1.iev = self.stk.n_ev/2
        wx.GetApp().TopWindow.page1.slider_eng.SetValue(wx.GetApp().TopWindow.page1.iev)
        
        wx.GetApp().TopWindow.page1.loadSpectrum(wx.GetApp().TopWindow.page1.ix, wx.GetApp().TopWindow.page1.iy)
        wx.GetApp().TopWindow.page1.loadImage()
        
        wx.GetApp().TopWindow.page0.ShowInfo(filelist[0], self.filepath)
        
        wx.EndBusyCursor()
        self.Destroy()
        
#---------------------------------------------------------------------- 
    def OnCancel(self, evt):
        self.Destroy()
               
        
        
#---------------------------------------------------------------------- 
class AboutFrame(wx.Frame):

    def __init__(self):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title = "About Mantis", size=(420, 700))
        
        #ico = logos.getlogo_2l_32Icon()
        #self.SetIcon(ico)
        
        
        self.com = wx.GetApp().TopWindow.common         
        self.fontsize = self.com.fontsize
        
        self.SetBackgroundColour("White")
        
        
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        
        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
        
        img = wx.Image(resource_path(os.path.join('images', 'Mantis_logo_about.png')), wx.BITMAP_TYPE_ANY).ConvertToBitmap()
        #img = logos.getMantis_logo_aboutImage()
        self.imageCtrl = wx.StaticBitmap(panel, wx.ID_ANY, img)
 
        vbox.Add(self.imageCtrl, 0, wx.ALL, 2)

        font1 = wx.Font(12, wx.SWISS, wx.NORMAL, wx.NORMAL)
        text1 = wx.StaticText(panel, 0, "www.2ndlookconsulting.com")
        text1.SetFont(font1)   
        text1.SetForegroundColour((53,159,217)) 
        
        font2 = wx.Font(12, wx.SWISS, wx.NORMAL, wx.NORMAL)
        text2 = wx.StaticText(panel, 0, '''Mantis 1.15''')  
        text2.SetFont(font2)        

        font3 = wx.Font(8, wx.SWISS, wx.NORMAL, wx.NORMAL)
        text3 = wx.StaticText(panel, 0, '''
Developed by Mirna Lerotic, based on earlier programs by Mirna 
Lerotic and Chris Jacobsen. Initial development supported by 
Argonne National Laboratory LDRD 2010-193-R1 9113. ''')  
        text3.SetFont(font3)          

        
        font4 = wx.Font(8, wx.SWISS, wx.NORMAL, wx.NORMAL)
        text4 = wx.StaticText(panel, 0, '''        
Mantis is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published 
by the Free Software Foundation, either version 3 of the License, 
or any later version.

Mantis is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty 
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU General Public License for more details 
http://www.gnu.org/licenses/.''')  
        text4.SetFont(font4)          


        vbox .Add((0,30))
        vbox.Add(text1,0, wx.LEFT, 50)  
        vbox.Add((0,10))    
        vbox.Add(text2,0, wx.LEFT, 50)   
        vbox.Add(text3,0, wx.LEFT, 50)   
        vbox.Add(text4,0, wx.LEFT, 50)  
        vbox.Add((0,10)) 

        
        
        button_close = wx.Button(panel, 0, 'Close')
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=button_close.GetId())
        vbox.Add(button_close, 0, wx.LEFT | wx.RIGHT | wx.ALIGN_RIGHT ,20)
        

        panel.SetSizer(vbox)
        
        vboxtop.Add(panel,0, wx.EXPAND )
        
        self.SetSizer(vboxtop)
        
        self.SetPosition((250, 150))
        
#----------------------------------------------------------------------              
    def OnClose(self, evt):
        self.Destroy()
        
    
""" ------------------------------------------------------------------------------------------------"""    
def show_splash():
    # create, show and return the splash screen
    bitmap = logos.getMantis_logo_splashBitmap()
    
    splash = wx.SplashScreen(bitmap, wx.SPLASH_CENTRE_ON_SCREEN|wx.SPLASH_NO_TIMEOUT, 3000, None, -1)
    splash.Show()
    return splash


""" ------------------------------------------------------------------------------------------------"""
def main():
    app = wx.App()
    #splash = show_splash()
   
    #time.sleep(1)
    frame = MainFrame(None, -1, 'Mantis')
    frame.Show()

    #splash.Destroy()
    app.MainLoop()

if __name__ == '__main__':
    main()

