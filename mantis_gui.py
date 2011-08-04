
from __future__ import division
import wx
import matplotlib as mtplot 
mtplot.interactive( True )
mtplot.use( 'WXAgg', warn=False )

import wxmpl
import numpy as npy
import os.path
from mpl_toolkits.axes_grid import make_axes_locatable
import time

import data_struct
import data_stack
import analyze
import logos



Winsizex = 1000
Winsizey = 740



#----------------------------------------------------------------------
class common:
    def __init__(self):
        self.fontsize = 8
        
        self.stack_loaded = 0
        self.i0_loaded = 0
        self.pca_calculated = 0
        self.cluster_calculated = 0
        self.spec_anl_calculated = 0

        self.font = wx.SystemSettings_GetFont(wx.SYS_DEFAULT_GUI_FONT)
        

""" ------------------------------------------------------------------------------------------------"""
class PageSpectral(wx.Panel):
    def __init__(self, parent, common, data_struct, stack, anlz):
        wx.Panel.__init__(self, parent)
        
        self.data_struct = data_struct
        self.stk = stack
        self.com = common
        self.anlz = anlz
        
        self.i_tspec = 1
        self.showraw = False
        self.show_scale_bar = 0
        
        self.SetBackgroundColour("White")
           
        self.fontsize = self.com.fontsize        
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        hboxT = wx.BoxSizer(wx.HORIZONTAL)
        hboxB = wx.BoxSizer(wx.HORIZONTAL)
    
        #panel 1        
        panel1 = wx.Panel(self, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_spmap = wx.TextCtrl(panel1, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_spmap.SetFont(self.com.font)
        self.tc_spmap.SetValue("Spectrum composition map")

        
        self.MapPanel = wxmpl.PlotPanel(panel1, -1, size =(3.5,3.5), cursor=False, crosshairs=False, location=False, zoom=False)                            
  
        vbox1.Add((0,10))
        vbox1.Add(self.tc_spmap,1, wx.LEFT | wx.EXPAND, 20)        
        vbox1.Add(self.MapPanel, 0,  wx.LEFT, 20)

        panel1.SetSizer(vbox1)
     
     
        #panel 2
        panel2 = wx.Panel(self, -1)
        vbox2 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_tspec = wx.TextCtrl(panel2, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_tspec.SetFont(self.com.font)
        self.tc_tspec.SetValue("Target Spectrum: ")
        hbox11 = wx.BoxSizer(wx.HORIZONTAL)         
        
        self.TSpectrumPanel = wxmpl.PlotPanel(panel2, -1, size=(5.65, 3.5), cursor=False, crosshairs=False, location=False, zoom=False)

        self.slider_tspec = wx.Slider(panel2, -1, 1, 1, 5, style=wx.SL_LEFT )        
        self.slider_tspec.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnTSScroll, self.slider_tspec)

        hbox11.Add(self.TSpectrumPanel, 0)
        hbox11.Add(self.slider_tspec, 0,  wx.EXPAND)
          
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

        self.button_save = wx.Button(panel3, -1, 'Save', (10,10))
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
        vbox4.Add(sizer41, 0, wx.EXPAND)
        
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
        self.rb_fit.SetValue(True)

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
         
        self.tc_speclist = wx.TextCtrl(panel5, -1, style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_speclist.SetFont(self.com.font)
        hbox51.Add(self.tc_speclist, 1, wx.EXPAND)
        sizer5.Add(hbox51,1, wx.ALL|wx.EXPAND,2)        
        panel5.SetSizer(sizer5)

        
        hboxB.Add(panel2, 0, wx.BOTTOM | wx.TOP, 9)
        hboxB.Add(panel1, 0, wx.BOTTOM | wx.TOP, 9)
        hboxT.Add((10,0)) 
               
        hboxT.Add(panel3, 1, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 10)
        hboxT.Add(panel4, 3, wx.LEFT | wx.RIGHT |wx.TOP | wx.EXPAND,10)
        hboxT.Add(panel5, 1, wx.LEFT | wx.RIGHT |wx.TOP | wx.EXPAND,10)

        vbox.Add(hboxT, 0, wx.ALL, 5)
        
        vbox.Add(hboxB, 0, wx.LEFT | wx.RIGHT, 5)
  
        self.SetSizer(vbox) 
        
        
        
#----------------------------------------------------------------------
    def OnTSpecFromFile(self, event):
        

        try: 
            wildcard = "Spectrum files (*.xas)|*.xas"
            dialog = wx.FileDialog(None, "Choose Spectrum file",
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
            
        except:
            wx.EndBusyCursor()  
            wx.MessageBox("Spectrum file not loaded.")
                                   
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

        try: 
            wx.BeginBusyCursor() 
            self.anlz.add_cluster_target_spectra()
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
            wx.MessageBox("CLuster spectra not loaded.")
 
                                                        
        wx.GetApp().TopWindow.refresh_widgets()
        
        
#----------------------------------------------------------------------
    def OnSave(self, event):
        
        wildcard = 'Portable Network Graphics (*.png)|*.png|Adobe PDF Files (*.pdf)|*.pdf|All files (*.*)|*.*'
      
        fileName = wx.FileSelector('Save', default_extension='png', wildcard = wildcard,
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
            
            suffix = "." + ext
                
            
            mtplot.rcParams['pdf.fonttype'] = 42
            
            
            from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas   
            
            colors=['#FF0000','#000000','#FFFFFF']
            spanclrmap=mtplot.colors.LinearSegmentedColormap.from_list('spancm',colors)
                       
            for i in range (self.anlz.n_target_spectra):
              
                #Save composition maps
                if self.showraw == True:
                    tsmapimage = self.anlz.target_svd_maps[:,:,i]
                else:
                    tsmapimage = self.anlz.target_pcafit_maps[:,:,i] 
  
                fig = mtplot.figure.Figure(figsize =(3.5,3.5))
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
                    p = mtplot.patches.Rectangle((5,self.stk.n_cols-10), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                                                 color = 'white', fill = True)
                    axes.add_patch(p)     
     
                im = axes.imshow(tsmapimage, cmap=spanclrmap, vmin = -bound, vmax = bound)
                cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)  
    
                axes.axis("off") 
                
                   
                fileName_img = fileName[:-len(suffix)]+"_TSmap_" +str(i+1)+"."+ext               
                fig.savefig(fileName_img)
                
                #Save spectra images
                tspectrum = self.anlz.target_spectra[i, :]
                tspectrumfit = self.anlz.target_pcafit_spectra[i, :]
            
        
                fig = mtplot.figure.Figure(figsize =(5.65, 3.5))
                canvas = FigureCanvas(fig)
                fig.clf()
                fig.add_axes((0.15,0.15,0.75,0.75))
                axes = fig.gca()
        
                mtplot.rcParams['font.size'] = self.fontsize

                specplot = axes.plot(self.stk.ev,tspectrum, color='black')
                specplot = axes.plot(self.stk.ev,tspectrumfit, color='green')
                        
                axes.set_xlabel('Photon Energy [eV]')
                axes.set_ylabel('Optical Density')

                fileName_spec = fileName[:-len(suffix)]+"_Tspectrum_" +str(i+1)+"."+ext
                fig.savefig(fileName_spec)                
                
                
            
        except IOError, e:
            if e.strerror:
                err = e.strerror 
            else: 
                err = e 
   
            wx.MessageBox('Could not save file: %s' % err, 'Error', 
                          parent=self, style=wx.OK|wx.ICON_ERROR) 
    
#----------------------------------------------------------------------        
    def OnTSScroll(self, event):
        sel = event.GetInt()
        self.i_tspec = sel
        if self.com.spec_anl_calculated == 1:
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
        self.tc_speclist.Clear() 
        self.tc_spfitlist.Clear()
        
        self.com.spec_anl_calculated = 0
        self.i_tspec = 1
        self.showraw = False
        self.rb_fit.SetValue(True)
        
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
        
        self.tc_speclist.Clear()     
        
        for i in range(self.anlz.n_target_spectra):
            self.tc_speclist.AppendText(self.anlz.tspec_names[i])
            
            if i < self.anlz.n_target_spectra:
                self.tc_speclist.AppendText('\n')
        self.tc_speclist.SetInsertionPoint(0)
        
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
            um_string = '$\mu m$'
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
        
        tspectrumfit = self.anlz.target_pcafit_spectra[self.i_tspec-1, :]
            
        
        fig = self.TSpectrumPanel.get_figure()
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        
        mtplot.rcParams['font.size'] = self.fontsize

        specplot = axes.plot(self.stk.ev,tspectrum, color='black')
        specplot = axes.plot(self.stk.ev,tspectrumfit, color='green')
                        
        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Optical Density')
        
        self.TSpectrumPanel.draw()
        
        self.tc_tspec.SetValue("Target Spectrum: " + 
                               self.anlz.tspec_names[self.i_tspec-1])
        
        self.textctrl_sp.Clear()
        
        self.textctrl_sp.AppendText('Common Name: '+ 
                                    self.anlz.tspec_names[self.i_tspec-1]+'\n')
        self.textctrl_sp.AppendText('RMS Error: '+ str('{0:7.5f}').format(self.anlz.target_rms[self.i_tspec-1]))
        
        self.ShowFitWeights()

        
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
        self.numclusters = 5
        self.wo_1st_pca = 0
        
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
        self.nclusterspin = wx.SpinCtrl(panel1, -1, '',  size= (60, -1), style=wx.ALIGN_LEFT)
        self.nclusterspin.SetRange(2,20)
        self.nclusterspin.SetValue(self.numclusters)
        self.Bind(wx.EVT_SPINCTRL, self.OnNClusterspin, self.nclusterspin)
        hbox11.Add((20,0))
        hbox11.Add(text1, 0, wx.TOP, 20)
        hbox11.Add((10,0))
        hbox11.Add(self.nclusterspin, 0, wx.TOP, 15)  
        hbox11.Add((20,0))       
            
        
        hbox12 = wx.BoxSizer(wx.HORIZONTAL)
        hbox12.Add((20,0))
        self.remove1stpcacb = wx.CheckBox(panel1, -1, 'Reduce thickness effects')
        self.remove1stpcacb.SetFont(self.com.font)
        self.Bind(wx.EVT_CHECKBOX, self.OnRemove1stpca, self.remove1stpcacb)
        hbox12.Add(self.remove1stpcacb, 0, wx.EXPAND|wx.TOP, 15)
        hbox12.Add((20,0))
        
        
        sizer1.Add((0,10)) 
        sizer1.Add(self.button_calcca, 0, wx.EXPAND|wx.LEFT|wx.RIGHT, 20)
        sizer1.Add(hbox11, 0, wx.EXPAND)
        sizer1.Add(hbox12, 0, wx.EXPAND)
        
        sizer1.Add((0,10))        
        sizer1.Add(wx.StaticLine(panel1), 0, wx.ALL|wx.EXPAND, 5)        
        sizer1.Add((0,10)) 
      
        sizer1.Add(self.button_scatterplots, 0, wx.EXPAND|wx.LEFT|wx.RIGHT, 20)
        sizer1.Add(self.button_savecluster, 0, wx.EXPAND|wx.LEFT|wx.RIGHT, 20)       
        sizer1.Add((0,10))
         
        panel1.SetSizer(sizer1)
        
        
                
        #panel 2        
        panel2 = wx.Panel(self, -1)
        vbox2 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_clustercomp = wx.TextCtrl(panel2, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_clustercomp.SetValue("Composite cluster image")        
        self.tc_clustercomp.SetFont(self.com.font)
        self.ClusterImagePan = wxmpl.PlotPanel(panel2, -1, size =(3.40,3.40), cursor=False, crosshairs=True, location=False, zoom=False)                              
        wxmpl.EVT_POINT(self, self.ClusterImagePan.GetId(), self.OnPointClusterImage)   
        vbox2.Add(self.tc_clustercomp, 0, wx.EXPAND) 
        vbox2.Add(self.ClusterImagePan, 0)   

        panel2.SetSizer(vbox2)
        
        
        #panel 3 
        panel3 = wx.Panel(self, -1)
        vbox3 = wx.BoxSizer(wx.VERTICAL)
        
        fgs = wx.FlexGridSizer(2, 3, 0, 0)
        
        self.tc_cluster = wx.TextCtrl(panel3, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_cluster.SetValue("Cluster ")
        self.tc_cluster.SetFont(self.com.font)
        hbox31 = wx.BoxSizer(wx.HORIZONTAL)
        self.ClusterIndvImagePan = wxmpl.PlotPanel(panel3, -1, size =(2.63,2.63), cursor=False, crosshairs=False, location=False, zoom=False)
    
        self.slidershow = wx.Slider(panel3, -1, self.selcluster, 1, 20, style=wx.SL_LEFT)   
        self.slidershow.Disable()    
        self.slidershow.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnClusterScroll, self.slidershow)
        
        text3 = wx.StaticText(panel3, -1, 'Cluster Distance Map',  style=wx.ALIGN_LEFT)
        text3.SetFont(self.com.font)
        self.ClusterDistMapPan = wxmpl.PlotPanel(panel3, -1, size =(2.63,2.63), cursor=False, crosshairs=False, location=False, zoom=False)
        
          
        fgs.AddMany([(self.tc_cluster), (wx.StaticText(panel3, -1, ' ')), (text3, 0, wx.LEFT, 15), 
                     (self.ClusterIndvImagePan), (self.slidershow, 0,  wx.EXPAND), (self.ClusterDistMapPan, 0, wx.LEFT, 15)])

        vbox3.Add(fgs)

        panel3.SetSizer(vbox3)
        
        
        #panel 4 
        panel4 = wx.Panel(self, -1)
        vbox4 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_clustersp = wx.TextCtrl(panel4, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_clustersp.SetValue("Cluster spectrum")
        self.tc_clustersp.SetFont(self.com.font)        

        self.ClusterSpecPan = wxmpl.PlotPanel(panel4, -1, size =(5.7,3.40), cursor=False, crosshairs=False, location=False, zoom=False)
        
        vbox4.Add(self.tc_clustersp, 0, wx.EXPAND)        
        vbox4.Add(self.ClusterSpecPan, 0)

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
        
        a=1
        if a==1:
        
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
        self.numclusters = num
        
               
 
#----------------------------------------------------------------------        
    def OnClusterScroll(self, event):
        sel = event.GetInt()
        self.selcluster = sel
        if self.com.cluster_calculated == 1:
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
        self.anlz.calculate_clusters_kmeans(self.numclusters, self.wo_1st_pca)
        
        
        
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
        
        mtplot.rcParams['font.size'] = self.fontsize
        
        
        im = axes.imshow(mapimage, cmap=mtplot.cm.get_cmap('gray'))
        axes.axis("off")
       
        self.ClusterDistMapPan.draw()
        
#----------------------------------------------------------------------           
    def OnRemove1stpca(self, event):
        if self.remove1stpcacb.GetValue():
            self.wo_1st_pca = 1
        else: self.wo_1st_pca = 0
        
        #recalculate CA!
        
        
#----------------------------------------------------------------------    
    def OnSave(self, event):     
               
        fileName = wx.FileSelector('Save', default_extension='png', 
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
            
            suffix = "." + ext
            
            fileName_evals = fileName[:-len(suffix)]+"_CAcimg."+ext          
            
            mtplot.rcParams['pdf.fonttype'] = 42
            
            fig = self.ClusterImagePan.get_figure()
            fig.savefig(fileName_evals)
            
            
            from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas          
            for i in range (self.numclusters):
              
                indvclusterimage = npy.zeros((self.anlz.stack.n_cols, self.anlz.stack.n_rows))+20.      
                ind = npy.where(self.anlz.cluster_indices == i)    
                colorcl = min(i,9)
                indvclusterimage[ind] = colorcl

                fig = mtplot.figure.Figure(figsize =(3.5,3.5))
                canvas = FigureCanvas(fig)
                fig.add_axes((0.02,0.02,0.96,0.96))
                axes = fig.gca()      
                mtplot.rcParams['font.size'] = self.fontsize        
                im = axes.imshow(indvclusterimage, cmap=self.clusterclrmap2, norm=self.bnorm2)
                axes.axis("off")
                   
                fileName_img = fileName[:-len(suffix)]+"_CAimg_" +str(i+1)+"."+ext               
                fig.savefig(fileName_img)
                
                
                clusterspectrum = self.anlz.clusterspectra[i, ]
                fig = mtplot.figure.Figure(figsize =(5.9,3.5))
                canvas = FigureCanvas(fig)
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

                fileName_spec = fileName[:-len(suffix)]+"_CAspectrum_" +str(i+1)+"."+ext
                fig.savefig(fileName_spec)                
                
                
            
        except IOError, e:
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

       
        self.colors=['#0000FF','#FF0000','#FFFF00','#33FF33','#B366FF',
                '#FF470A','#33FFFF','#006600','#CCCC99','#993300',
                '#000000']
        
        self.clusterclrmap1=mtplot.colors.LinearSegmentedColormap.from_list('clustercm',self.colors)
     
        self.bnorm1 = mtplot.colors.BoundaryNorm(colors_i, self.clusterclrmap1.N)
        
        colors_i = npy.linspace(0,self.maxclcolors+2,self.maxclcolors+3)
        
        #use black color for clusters > maxclcolors, the other 2 colors are for background      
        colors2=['#0000FF','#FF0000','#FFFF00','#33FF33','#B366FF',
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

        self.ScatterPPanel = wxmpl.PlotPanel(panel, -1, size=(5.0, 4.0), cursor=False, crosshairs=False, location=False, zoom=False)
        
        self.slidershow_y = wx.Slider(panel, -1, self.pca_y, 1, self.numsigpca, style=wx.SL_RIGHT|wx.SL_LABELS|wx.SL_INVERSE)
        self.slidershow_y.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnSliderScroll_y, self.slidershow_y)
        
        grid1.Add(self.slidershow_y, 0, wx.EXPAND)
        grid1.Add(self.ScatterPPanel, 0)
               
        
        self.slidershow_x = wx.Slider(panel, -1, self.pca_x, 1, self.numsigpca, style=wx.SL_TOP|wx.SL_LABELS)
        self.slidershow_x.SetFocus()        
        self.Bind(wx.EVT_SCROLL, self.OnSliderScroll_x, self.slidershow_x)
        
        grid1.Add(wx.StaticText(panel, -1, ''))
        grid1.Add(self.slidershow_x, 0,  wx.EXPAND)
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.button_savescatt = wx.Button(panel, -1, 'Save')
        self.button_savescatt.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnSaveScatt, id=self.button_savescatt.GetId())
        hbox.Add(self.button_savescatt, -1, wx.ALIGN_RIGHT|wx.RIGHT, 20)
        
        button_close = wx.Button(panel, -1, 'Close')
        button_close.SetFont(self.com.font)
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=button_close.GetId())
        hbox.Add(button_close, -1, wx.ALIGN_RIGHT|wx.RIGHT, 40)

        
        vbox.Add(grid1, 0, wx.EXPAND|wx.ALL,20)
        vbox.Add(hbox, -1, wx.ALIGN_RIGHT)
                
        panel.SetSizer(vbox)
        
        vboxtop.Add(panel,-1, wx.EXPAND|wx.ALIGN_CENTER_VERTICAL)
        
        self.SetSizer(vboxtop)
        
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
    def OnSaveScatt(self, event): 
          
               
        fileName = wx.FileSelector('Save Plots', default_extension='png', 
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
            wx.BeginBusyCursor()
            
            suffix = "." + ext
            
            nplots = 0
            for ip in range(self.numsigpca):
                for jp in range(self.numsigpca):
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
            for ip in range(self.numsigpca):
                for jp in range(self.numsigpca):
                    if jp >= (ip+1):
                        

                        x_comp = self.od_reduced[:,ip]
                        y_comp = self.od_reduced[:,jp]
                        if nplots > 1 :
                            axes = fig.add_subplot(nplotsrows,2, pplot)
                        else:
                            axes = fig.add_subplot(1,1,1)
                            
                        pplot = pplot+1
                        
                        for i in range(self.numclusters):
                            thiscluster = npy.where(self.clindices == i)
                            axes.plot(x_comp[thiscluster], y_comp[thiscluster],'.',color=self.colors[i],alpha=0.5)
                        axes.set_xlabel('Component '+str(ip+1))
                        axes.set_ylabel('Component '+str(jp+1))
                            
    
            fileName_sct = fileName[:-len(suffix)]+"_CAscatterplot_" +str(i+1)+"."+ext
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
    def OnClose(self, evt):
        self.Close(True)       
    
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
        
          
        self.fontsize = self.com.fontsize
        
  
        #panel 1        
        panel1 = wx.Panel(self, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_PCAcomp = wx.TextCtrl(panel1, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_PCAcomp.SetFont(self.com.font)
        self.tc_PCAcomp.SetValue("PCA component ")
        
        hbox11 = wx.BoxSizer(wx.HORIZONTAL)
   
        self.PCAImagePan = wxmpl.PlotPanel(panel1, -1, size =(3.9,3.4), cursor=False, crosshairs=False, location=False, zoom=False)
                              
        self.slidershow = wx.Slider(panel1, -1, self.selpca, 1, 20, style=wx.SL_LEFT)
        self.slidershow.Disable()          
        self.slidershow.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnPCAScroll, self.slidershow)

        hbox11.Add(self.PCAImagePan, 0)
        hbox11.Add(self.slidershow, 0,  wx.EXPAND)
        
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
        
        self.PCASpecPan = wxmpl.PlotPanel(panel3, -1, size =(5.4,3.4), cursor=False, crosshairs=False, location=False, zoom=False)
             
        vbox3 = wx.BoxSizer(wx.VERTICAL)
        self.text_pcaspec = wx.TextCtrl(panel3, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.text_pcaspec.SetFont(self.com.font)
        self.text_pcaspec.SetValue("PCA spectrum ")        
        vbox3.Add(self.text_pcaspec, 0)
        vbox3.Add(self.PCASpecPan, 0)        
        panel3.SetSizer(vbox3)
        
        
        #panel 4
        panel4 = wx.Panel(self, -1)
             
        self.PCAEvalsPan = wxmpl.PlotPanel(panel4, -1, size =(5.4,2.6), cursor=False, crosshairs=False, location=False, zoom=False)
        
        vbox4 = wx.BoxSizer(wx.VERTICAL)
        text4 = wx.TextCtrl(panel4, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        text4.SetFont(self.com.font)
        text4.SetValue("PCA eigenvalues ")        
        vbox4.Add(text4, 0)
        vbox4.Add(self.PCAEvalsPan, 0)
        
        panel4.SetSizer(vbox4)      
        

        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        gridtop = wx.FlexGridSizer(2, 2, vgap=10, hgap=20)
        gridtop.Add(panel2, 0, wx.LEFT|wx.TOP, 20)
        gridtop.Add(panel4, 0, wx.ALIGN_RIGHT)
        
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
     
        #Keiser criterion
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
            self.showEvals()
            
#----------------------------------------------------------------------    
    def OnSave(self, event):     
        
        wildcard = 'Portable Network Graphics (*.png)|*.png|Adobe PDF Files (*.pdf)|*.pdf|All files (*.*)|*.*'               
        fileName = wx.FileSelector('Save Plot', default_extension='png', 
                                   wildcard=wildcard, parent=self, flags=wx.SAVE|wx.OVERWRITE_PROMPT) 
   
        if not fileName: 
            return 

        path, ext = os.path.splitext(fileName) 
        ext = ext[1:].lower() 
        
       
        if ext != 'png' and ext != 'pdf': 
            error_message = ( 
                  'Only the PNG and PDF image formats are supported.\n' 
                 'A file extension of `png\' or `pdf\' must be used.') 
            wx.MessageBox(error_message, 'Error - Could not save file.'+error_message, 
                  parent=self, style=wx.OK|wx.ICON_ERROR) 
            return 
   
        try: 
            
            suffix = "." + ext
            fileName_evals = fileName[:-len(suffix)]+"_PCAevals."+ext
            mtplot.rcParams['pdf.fonttype'] = 42
            
            fig = self.PCAEvalsPan.get_figure()
            fig.savefig(fileName_evals)
            
            from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
            
            
            for i in range (self.numsigpca):
              
                self.pcaimage = self.anlz.pcaimages[:,:,i]
              
                fig = mtplot.figure.Figure(figsize =(4.0,3.5))
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
                                
                fileName_img = fileName[:-len(suffix)]+"_PCA_" +str(i+1)+"."+ext
                fig.savefig(fileName_img)
            
                
                self.pcaspectrum = self.anlz.eigenvecs[:,i]
                fig = mtplot.figure.Figure(figsize =(5.65,3.5))
                canvas = FigureCanvas(fig)
                fig.add_axes((0.15,0.15,0.75,0.75))
                axes = fig.gca()
                mtplot.rcParams['font.size'] = self.fontsize
                specplot = axes.plot(self.stk.ev,self.pcaspectrum)    
                axes.set_xlabel('Photon Energy [eV]')
                axes.set_ylabel('Optical Density')
                
                fileName_spec = fileName[:-len(suffix)]+"_PCAspectrum_" +str(i+1)+"."+ext
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
        
        self.pcaevals = self.anlz.eigenvals[0:40]
        

        fig = self.PCAEvalsPan.get_figure()
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        
        mtplot.rcParams['font.size'] = self.fontsize
       
        evalsplot = axes.semilogy(npy.arange(1,41), self.pcaevals,'b.')    
        
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
        self.sel = 50
        self.showflux = True
        self.fontsize = self.com.fontsize
        
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
        
        self.show_scale_bar = 0

        vbox = wx.BoxSizer(wx.VERTICAL)
        hboxT = wx.BoxSizer(wx.HORIZONTAL)
        hboxB = wx.BoxSizer(wx.HORIZONTAL)
    
        #panel 1        
        panel1 = wx.Panel(self, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_imageeng = wx.TextCtrl(panel1, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_imageeng.SetFont(self.com.font)
        self.tc_imageeng.SetValue("Image at energy: ")

       
        hbox11 = wx.BoxSizer(wx.HORIZONTAL)
   
        self.AbsImagePanel = wxmpl.PlotPanel(panel1, -1, size =(3.5,3.5), cursor=False, crosshairs=True, location=False, zoom=False)
        wxmpl.EVT_POINT(panel1, self.AbsImagePanel.GetId(), self.OnPointAbsimage)
                              
        self.slider_eng = wx.Slider(panel1, -1, self.sel, 0, 100, style=wx.SL_LEFT )        
        self.slider_eng.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnScrollEng, self.slider_eng)

        hbox11.Add(self.AbsImagePanel, 0)
        hbox11.Add(self.slider_eng, 0,  wx.EXPAND)
        
        vbox1.Add(self.tc_imageeng,1, wx.LEFT | wx.TOP | wx.EXPAND, 20)        
        vbox1.Add(hbox11, 0,  wx.LEFT, 20)

        panel1.SetSizer(vbox1)
     
     
        #panel 2
        panel2 = wx.Panel(self, -1)
        vbox2 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_spec = wx.TextCtrl(panel2, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_spec.SetValue("Spectrum at point: ")
        self.tc_spec.SetFont(self.com.font)
          
        
        self.SpectrumPanel = wxmpl.PlotPanel(panel2, -1, size=(5.75, 3.5), cursor=False, crosshairs=False, location=False, zoom=False)
        wxmpl.EVT_POINT(panel2, self.SpectrumPanel.GetId(), self.OnPointSpectrum)
        
        vbox2.Add(self.tc_spec, 1, wx.LEFT | wx.TOP | wx.EXPAND, 20)       
        vbox2.Add(self.SpectrumPanel, 0, wx.LEFT , 20)
        
        panel2.SetSizer(vbox2)
        
        
        #panel 3
        panel3 = wx.Panel(self, -1)
        sizer1 = wx.StaticBoxSizer(wx.StaticBox(panel3, -1, 'Preprocess'),orient=wx.VERTICAL)
        vbox31 = wx.BoxSizer(wx.VERTICAL)
        vbox31.Add((0,10)) 
        
        self.button_limitev = wx.Button(panel3, -1, 'Limit energy range...')
        self.Bind(wx.EVT_BUTTON, self.OnLimitEv, id=self.button_limitev.GetId())
        self.button_limitev.SetFont(self.com.font)
        self.button_limitev.Disable()
        vbox31.Add(self.button_limitev, 0, wx.EXPAND)
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
        self.button_save = wx.Button(panel3, -1, 'Save')
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
        sb.SetBackgroundColour("White")
        sizer4 = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)
        sb = wx.StaticBox(panel4, -1, 'File')
        sb.SetBackgroundColour("White")
        sizer41 = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)

        self.textctrl = wx.TextCtrl(panel4, -1, style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        sizer41.Add(self.textctrl, 1, wx.EXPAND|wx.TOP|wx.LEFT, 5)
        vbox4.Add(sizer41, 0, wx.EXPAND)
        
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

        self.add_scale_cb = wx.CheckBox(panel4, -1, '  Scale')
        self.add_scale_cb.SetFont(self.com.font)
        self.Bind(wx.EVT_CHECKBOX, self.OnShowScale, self.add_scale_cb)

        sizer42.Add((0,3))
        sizer42.Add(self.rb_flux)
        sizer42.Add((0,5))
        sizer42.Add(self.rb_od)
        sizer42.Add((0,10))
        sizer42.Add(self.add_scale_cb)
        
        hbox41.Add(sizer42, 0, wx.EXPAND)
                

        sizer43 = wx.StaticBoxSizer(wx.StaticBox(panel4, -1, 'Display settings'),  orient=wx.VERTICAL)
        hbox42 = wx.BoxSizer(wx.HORIZONTAL)
        hbox42.Add((15,0))

        fgs41 = wx.FlexGridSizer(3, 2, 2, 5)
        min = wx.StaticText(panel4, label="Minimum")
        min.SetFont(self.com.font)
        max = wx.StaticText(panel4, label="Maximum")
        max.SetFont(self.com.font)
        self.slider_brightness_min = wx.Slider(panel4, -1, self.dispbrightness_min, 0, 49, style=wx.SL_HORIZONTAL)        
        self.slider_brightness_min.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnScrollBrightnessMin, self.slider_brightness_min)
        
        self.slider_brightness_max = wx.Slider(panel4, -1, self.dispbrightness_max, 50, 100, style=wx.SL_HORIZONTAL)        
        self.slider_brightness_max.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnScrollBrightnessMax, self.slider_brightness_max)        
        
        gamma = wx.StaticText(panel4, label="Gamma")
        gamma.SetFont(self.com.font)
        self.slider_gamma = wx.Slider(panel4, -1, self.displaygamma, 1, 20, style=wx.SL_HORIZONTAL)        
        self.slider_gamma.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnScrollGamma, self.slider_gamma)

        
        fgs41.AddMany([(min), (self.slider_brightness_min, 0, wx.EXPAND), (max), 
            (self.slider_brightness_max, 0, wx.EXPAND),(gamma), (self.slider_gamma, 0, wx.EXPAND)])
      
        hbox42.Add(fgs41, 0, wx.EXPAND)
        hbox42.Add((20,0))

        
        vbox43 = wx.BoxSizer(wx.VERTICAL)
        self.button_despike = wx.Button(panel4, -1, 'Despike')
        self.button_despike.SetFont(self.com.font)
        #self.Bind(wx.EVT_BUTTON, self.onDespike, id=self.button_despike.GetId())   
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
        
        hbox42.Add(vbox43, 1, wx.EXPAND)                   
        sizer43.Add(hbox42)
    
        hbox41.Add((10,0))
        hbox41.Add(sizer43, 1, wx.EXPAND)
        vbox4.Add(hbox41, 1, wx.EXPAND)

        sizer4.Add(vbox4,1, wx.EXPAND)
        
        panel4.SetSizer(sizer4)
        
        
        #panel 5
        panel5 = wx.Panel(self, -1)
        sizer5 = wx.StaticBoxSizer(wx.StaticBox(panel5, -1, 'Region of Interest'), wx.VERTICAL)
        
        vbox51 = wx.BoxSizer(wx.VERTICAL)
        vbox51.Add((0,2))
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
        
        vbox51.Add((0,8))        
        
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
        hboxT.Add(panel4, 3, wx.LEFT | wx.RIGHT |wx.TOP | wx.EXPAND,10)
        hboxT.Add(panel5, 1, wx.LEFT | wx.RIGHT |wx.TOP | wx.EXPAND,10)

        vbox.Add(hboxT, 0, wx.ALL, 5)
        
        vbox.Add(hboxB, 0,  wx.ALL, 5)
  
        self.SetSizer(vbox) 
        
        

      

#----------------------------------------------------------------------        
    def loadImage(self):
               
        self.image = 0     
        
        if self.defaultdisplay == 1.0:
            #use a pointer to the data not a copy
            if self.showflux:
                #Show flux image      
                self.image = self.stk.absdata[:,:,self.iev]#.copy() 
            else:
                #Show OD image
                self.image = self.stk.od3d[:,:,self.iev]#.copy()
        else:   
            #Adjustment to the data display setting has been made so make a copy
            if self.showflux:
                self.image = self.stk.absdata[:,:,self.iev].copy() 
            else:
                self.image = self.stk.od3d[:,:,self.iev].copy() 
                


                      
        fig = self.AbsImagePanel.get_figure()
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        axes = fig.gca()
        fig.patch.set_alpha(1.0)
        
        if (self.line != None) and (self.addroi == 1):
            axes.add_line(self.line)
            
         

        if self.defaultdisplay == 1.0:
            im = axes.imshow(self.image, cmap=mtplot.cm.get_cmap(self.colortable)) 
        else:
            imgmax = npy.amax(self.image)
            imgmin = npy.amin(self.image)
            if (self.gamma != 1.0) or (imgmin < 0.0):
                self.image = (self.image-imgmin)/(imgmax-imgmin)
                imgmax = 1.0
                imgmin = 0.0
                if (self.gamma != 1.0):
                    self.image = npy.power(self.image, self.gamma)
            im = axes.imshow(self.image, cmap=mtplot.cm.get_cmap(self.colortable), 
                             vmin=(imgmin+imgmax*self.brightness_min),vmax=imgmax*self.brightness_max)
            
        if (self.showROImask == 1) and (self.addroi == 1):
            im_red = axes.imshow(self.ROIpix_masked,cmap=mtplot.cm.get_cmap("autumn")) 
         
        axes.axis("off") 
        

        
        if self.show_scale_bar == 1:
            um_string = '$\mu m$'
            microns = '$'+self.stk.scale_bar_string+' $'+um_string
            axes.text(self.stk.scale_bar_pixels_x+10,self.stk.n_cols-9, microns, horizontalalignment='left', verticalalignment='center',
                      color = 'black', fontsize=14)
            #Matplotlib has flipped scales so I'm using rows instead of cols!
            p = mtplot.patches.Rectangle((5,self.stk.n_cols-10), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                                   color = 'black', fill = True)
            axes.add_patch(p)
        
        self.AbsImagePanel.draw()
        
        self.tc_imageeng.SetValue("Image at energy: " +str(self.stk.ev[self.iev])+" eV")

#----------------------------------------------------------------------          
    def loadSpectrum(self, xpos, ypos):

        self.spectrum = self.stk.od3d[xpos,ypos, :]

        fig = self.SpectrumPanel.get_figure()
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()
        
        mtplot.rcParams['font.size'] = self.fontsize

        specplot = axes.plot(self.stk.ev,self.spectrum)
        
        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Optical Density')
        
        axes.axvline(x=self.stk.ev[self.iev], color = 'g', alpha=0.5)

        
        self.SpectrumPanel.draw()
        
        self.tc_spec.SetValue("Spectrum at point: [" +str(ypos)+", " + str(xpos)+"] ")

        
#----------------------------------------------------------------------            
    def OnScrollEng(self, event):
        self.sel = event.GetInt()
        self.iev = self.sel
        if self.com.stack_loaded == 1:
            self.loadImage()
            
                   
#----------------------------------------------------------------------  
    def OnPointSpectrum(self, evt):
        x = evt.xdata
        y = evt.ydata
        
        if (self.com.i0_loaded == 1) and (self.addroi == 0):      
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
            
                   
#----------------------------------------------------------------------  
    def OnPointAbsimage(self, evt):
        x = evt.xdata
        y = evt.ydata
        
        if (self.com.i0_loaded == 1) and (self.addroi == 0):      
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
            wildcard = "I0 files (*.xas)|*.xas"
            dialog = wx.FileDialog(None, "Choose i0 file",
                                   wildcard=wildcard,
                                   style=wx.OPEN)
            if dialog.ShowModal() == wx.ID_OK:
                filepath_i0 = dialog.GetPath()
                self.filename = dialog.GetFilename()
                                                        
            wx.BeginBusyCursor()                                    
                        
            x=self.stk.n_cols
            y=self.stk.n_rows
            z=self.iev               

            self.ix = x/2
            self.iy = y/2

            self.stk.read_stk_i0(filepath_i0)
            
            self.loadSpectrum(self.ix, self.iy)
            self.loadImage()
            self.com.i0_loaded = 1
            wx.EndBusyCursor()
            
        except:
            wx.BeginBusyCursor()  
            self.com.i0_loaded = 0        
            wx.EndBusyCursor()  
            wx.MessageBox("I0 file not loaded.")
            
                       
        dialog.Destroy()
                          
        wx.GetApp().TopWindow.refresh_widgets()

        
        
#----------------------------------------------------------------------       
    def OnI0histogram(self, event):    
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
    def OnSave(self, event):     
               
        fileName = wx.FileSelector('Save Plot', default_extension='png', 
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
            
            suffix = "." + ext
            fileName_spec = fileName[:-len(suffix)]+"_spectrum."+ext

            mtplot.rcParams['pdf.fonttype'] = 42
            
            fig = self.SpectrumPanel.get_figure()
            fig.savefig(fileName_spec)

            
            fileName_img = fileName[:-len(suffix)]+"_" +str(self.stk.ev[self.iev])+"eV."+ext
            fig = self.AbsImagePanel.get_figure()
            fig.savefig(fileName_img)

            
        except IOError, e:
            if e.strerror:
                err = e.strerror 
            else: 
                err = e 
   
            wx.MessageBox('Could not save file: %s' % err, 'Error', 
                          parent=self, style=wx.OK|wx.ICON_ERROR) 


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
    def onResetDisplaySettings(self, event):
        
        self.ResetDisplaySettings()
        self.loadImage()
        
#----------------------------------------------------------------------
    def OnScrollBrightnessMin(self, event):
        
        self.dispbrightness_min = event.GetInt()
        
        self.brightness_min = float(self.dispbrightness_min)/100.0
        
        self.defaultdisplay = 0.0
        
        if self.com.stack_loaded == 1:
            self.loadImage()
        
#----------------------------------------------------------------------
    def OnScrollBrightnessMax(self, event):
        
        self.dispbrightness_max = event.GetInt()
        
        self.brightness_max = float(self.dispbrightness_max)/100.0
        
        self.defaultdisplay = 0.0
        
        if self.com.stack_loaded == 1:
            self.loadImage()
        
#----------------------------------------------------------------------
    def OnScrollGamma(self, event):
        
        self.displaygamma = event.GetInt()
        
        self.gamma = float(self.displaygamma)/10.0  
        
        self.defaultdisplay = 0.0
        
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
 
#----------------------------------------------------------------------        
    def OnLimitEv(self, evt):    
        LimitEv(self.com, self.stk).Show()
        
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
    def OnSaveROISpectrum(self, event):  
               
        fileName = wx.FileSelector('Save ROI Spectrum (.xas)', default_extension='xas', 
                                   wildcard=('XAS (*.xas)|*.xas'), 
                                              parent=self, flags=wx.SAVE|wx.OVERWRITE_PROMPT) 
   
        if not fileName: 
            return 

        path, ext = os.path.splitext(fileName) 
        ext = ext[1:].lower() 
   
        try:
            if (self.com.i0_loaded == 1):
                self.stk.write_xas(fileName, self.stk.ev, self.ROIspectrum)
            else:
                self.CalcROI_I0Spectrum()
                self.stk.write_xas(fileName, self.stk.ev, self.ROIspectrum)
                     
                
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
               
        self.HistogramPanel = wxmpl.PlotPanel(panel, -1, size=(6.0, 3.7), cursor=False, crosshairs=False, location=False, zoom=False)
        
        wxmpl.EVT_SELECTION(panel, self.HistogramPanel.GetId(), self.OnSelection)

        vbox.Add(self.HistogramPanel, 0, wx.ALL, 20)
        
       
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        sizer2 = wx.StaticBoxSizer(wx.StaticBox(panel, -1, 'I0 pixels'), orient=wx.VERTICAL)
        self.textctrl = wx.TextCtrl(panel, -1, size = (565, 20), style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.textctrl.SetValue('Selection: [ '+str(self.histmin) + ' kHz, '+ str(self.histmax)+' kHz ]' )
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
        hbox.Add(button_cancel, 1, wx.ALL ,20)
        
        vbox.Add(hbox, 0 )
        
        panel.SetSizer(vbox)
        
        vboxtop.Add(panel,1, wx.EXPAND )
        
        self.SetSizer(vboxtop)
        
        self.draw_histogram()
        self.draw_image()
        
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
        
   
        self.image = self.stack.absdata[:,:,self.stack.n_ev/2].copy() 
       
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
      
        im = axes.imshow(self.image, cmap=mtplot.cm.get_cmap("gray")) 

        im_red = axes.imshow(redpix_masked,cmap=mtplot.cm.get_cmap("autumn"))  
         
        axes.axis("off")  
        self.AbsImagePanel.draw()


#----------------------------------------------------------------------        
    def OnSelection(self, evt):
        
        x1, y1 = evt.x1data, evt.y1data
        x2, y2 = evt.x2data, evt.y2data
        
        self.histmin = x1
        self.histmax = x2

        self.textctrl.SetValue('Selection: [ '+str(self.histmin) + ' kHz, '+ str(self.histmax)+' kHz ]' )

        self.draw_histogram()
        self.draw_image()

#----------------------------------------------------------------------        
    def OnAccept(self, evt):
        
        self.stack.i0_from_histogram(self.histmin, self.histmax)
        self.Close(True) 
        wx.GetApp().TopWindow.page1.I0histogramCalculated()

                
#---------------------------------------------------------------------- 
    def OnCancel(self, evt):
        self.Close(True)   

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
        self.limitevmax = self.stack.original_n_ev-1
        
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
               
        self.SpectrumPanel = wxmpl.PlotPanel(panel, -1, size=(6.0, 3.7), cursor=False, crosshairs=False, location=False, zoom=False)
        
        wxmpl.EVT_SELECTION(panel, self.SpectrumPanel.GetId(), self.OnSelection)

        vbox.Add(self.SpectrumPanel, 0, wx.ALL, 20)
        
       
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
        hbox.Add(button_cancel, 1, wx.ALL ,20)
        
        vbox.Add(hbox, 0 )
        
        panel.SetSizer(vbox)
        
        vboxtop.Add(panel,1, wx.EXPAND )
        
        self.SetSizer(vboxtop)
        
        self.draw_limitev_plot()
        
      
#----------------------------------------------------------------------        
    def draw_limitev_plot(self):
        
        odtotal = self.stack.original_od3d.sum(axis=0)   
        odtotal = odtotal.sum(axis=0)/(self.stack.original_n_rows*self.stack.original_n_cols) 
        
        fig = self.SpectrumPanel.get_figure()
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        self.axes = fig.gca()
        
        mtplot.rcParams['font.size'] = self.fontsize

        specplot = self.axes.plot(self.stack.original_ev,odtotal)
        
        self.axes.set_xlabel('Photon Energy [eV]')
        self.axes.set_ylabel('Optical Density')
        
        if self.evlimited == 1:
            self.axes.axvspan(self.stack.original_ev[self.limitevmin], self.stack.original_ev[self.limitevmax], facecolor='g', alpha=0.5)

        
        self.SpectrumPanel.draw()
        
        self.textctrl.SetValue("Min energy " + str(self.stack.original_ev[self.limitevmin]) + " eV\n"
                              + "Max energy " + str(self.stack.original_ev[self.limitevmax]) + " eV")
        
        
#----------------------------------------------------------------------        
    def OnSelection(self, evt):

        x1, y1 = evt.x1data, evt.y1data
        x2, y2 = evt.x2data, evt.y2data

        
        
        self.limitevmin = npy.abs(self.stack.original_ev-x1).argmin()
        self.limitevmax = npy.abs(self.stack.original_ev-x2).argmin()
        
        self.evlimited = 1
             
        self.draw_limitev_plot()
        
#----------------------------------------------------------------------        
    def OnAccept(self, evt):
        #change the energy range to limitevmin-limitev-max  
        #print self.stack.n_ev, self.stack.ev.shape
        self.stack.n_ev = self.limitevmax+1-self.limitevmin
        self.stack.ev = self.stack.original_ev[self.limitevmin:self.limitevmax+1]
        
        #print self.stack.n_ev, self.stack.ev.shape
        
        
        self.stack.absdata = self.stack.original_absdata[:,:,self.limitevmin:self.limitevmax+1]
        
        self.stack.od3d = self.stack.original_od3d[:,:,self.limitevmin:self.limitevmax+1]
        
        self.stack.od = self.stack.od3d.copy()
        
        self.stack.od = npy.reshape(self.stack.od, (self.stack.n_rows*self.stack.n_cols, self.stack.n_ev), order='F')
        
        #Fix the slider on Page 1! 
        wx.GetApp().TopWindow.page1.slider_eng.SetRange(0,self.stack.n_ev-1)
        wx.GetApp().TopWindow.page1.iev = self.stack.n_ev/2
        wx.GetApp().TopWindow.page1.slider_eng.SetValue(wx.GetApp().TopWindow.page1.iev)
        
        wx.GetApp().TopWindow.page1.loadSpectrum(wx.GetApp().TopWindow.page1.ix, wx.GetApp().TopWindow.page1.iy)
        wx.GetApp().TopWindow.page1.loadImage()
        
        self.Close(True)
        
#---------------------------------------------------------------------- 
    def OnCancel(self, evt):
        self.Close(True)
        
        
               
#---------------------------------------------------------------------- 
class SpectralROI(wx.Frame):

    title = "Spectral Regions of Interest"

    def __init__(self, common, stack):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=self.title, size=(630, 700))
               
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
        
                
        self.odtotal = self.stack.original_od3d.sum(axis=0)   
        self.odtotal = self.odtotal.sum(axis=0)/(self.stack.n_rows*self.stack.n_cols) 
        
        self.image_i0 = npy.zeros((self.stack.n_cols, self.stack.n_rows))
        self.image_i = npy.zeros((self.stack.n_cols, self.stack.n_rows))
        self.odthickmap = npy.zeros((self.stack.n_cols, self.stack.n_rows))
        
    
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
               
        self.SpectrumPanel = wxmpl.PlotPanel(panel, -1, size=(6.0, 3.7), cursor=False, crosshairs=False, location=False, zoom=False)
        
        wxmpl.EVT_SELECTION(panel, self.SpectrumPanel.GetId(), self.OnSelection)

        vbox.Add(self.SpectrumPanel, 0, wx.ALL, 20)
        
       
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

        self.ODMImagePanel = wxmpl.PlotPanel(panel, -1, size =(2.0,2.0), cursor=False, crosshairs=False, location=False, zoom=False)
        hbox2.Add(self.ODMImagePanel, 0, wx.LEFT|wx.RIGHT,20)
        hbox2.Add(sizer2, 1, wx.EXPAND|wx.LEFT|wx.RIGHT,20)
        
        vbox.Add(hbox2, 0, wx.EXPAND)      
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        
        button_ok = wx.Button(panel, -1, 'Save')
        self.Bind(wx.EVT_BUTTON, self.OnAccept, id=button_ok.GetId())
        hbox.Add(button_ok, 1, wx.ALL,20)
        
        button_cancel = wx.Button(panel, -1, 'Dismiss')
        self.Bind(wx.EVT_BUTTON, self.OnCancel, id=button_cancel.GetId())
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
    
#----------------------------------------------------------------------        
    def draw_image(self):

               
        fig = self.ODMImagePanel.get_figure()
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        
        
        axes = fig.gca()
        
        im = axes.imshow(self.odthickmap, cmap=mtplot.cm.get_cmap("gray")) 

         
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
    def OnAccept(self, evt):
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
        self.Close(True)   

        
    
#---------------------------------------------------------------------- 
class PlotFrame(wx.Frame):

    def __init__(self, datax, datay):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title = "I0 data", size=(630, 500))
        
        ico = logos.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White")
        
        
        self.com = wx.GetApp().TopWindow.common         
        self.fontsize = self.com.fontsize
        
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
               
        self.PlotPanel = wxmpl.PlotPanel(panel, -1, size=(6.0, 3.7), cursor=False, crosshairs=False, location=False, zoom=False)

        vbox.Add(self.PlotPanel, 0, wx.ALL, 20)
        

        
        button_close = wx.Button(panel, -1, 'Close')
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=button_close.GetId())
        vbox.Add(button_close, 0, wx.ALL | wx.ALIGN_RIGHT ,20)
        

        panel.SetSizer(vbox)
        
        vboxtop.Add(panel,1, wx.EXPAND )
        
        self.SetSizer(vboxtop)
        
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
    def OnClose(self, evt):
        self.Close(True)
        
        
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
        
#----------------------------------------------------------------------          
    def OnColorTable(self, event):
        
        radioSelected = event.GetEventObject()

        wx.GetApp().TopWindow.page1.colortable = radioSelected.GetLabel()
        
        wx.GetApp().TopWindow.page1.loadImage()
        
#----------------------------------------------------------------------              
    def OnClose(self, evt):
        self.Hide()
        
        
            
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
        self.common = common()
        
        self.SetFont(self.common.font)
              

        # Here we create a panel and a notebook on the panel
        p = wx.Panel(self)
        nb = wx.Notebook(p, style=wx.BORDER_STATIC)

        # create the page windows as children of the notebook
        self.page1 = PageStack(nb, self.common, self.data_struct, self.stk)
        self.page2 = PagePCA(nb, self.common, self.data_struct, self.stk, self.anlz)
        self.page3 = PageCluster(nb, self.common, self.data_struct, self.stk, self.anlz)
        self.page4 = PageSpectral(nb, self.common, self.data_struct, self.stk, self.anlz)

        # add the pages to the notebook with the label to show on the tab
        nb.AddPage(self.page1, "Image Stack")
        nb.AddPage(self.page2, "PCA")
        nb.AddPage(self.page3, "Cluster Analysis")
        nb.AddPage(self.page4, "Spectral Analysis")
        

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
        openTool = self.toolbar.AddSimpleTool(wx.ID_OPEN, open_ico, "Open", "Open stack")
        self.Bind(wx.EVT_MENU, self.onBrowse, openTool)
        
        save_ico = wx.ArtProvider.GetBitmap(wx.ART_FILE_SAVE, wx.ART_TOOLBAR, (16,16))
        saveTool = self.toolbar.AddSimpleTool(wx.ID_SAVE, save_ico, "Save", "Save analysis")
        self.toolbar.EnableTool(wx.ID_SAVE, False)       
        #self.Bind(wx.EVT_MENU, self.onBrowse, saveTool)     
        
        saveas_ico = wx.ArtProvider.GetBitmap(wx.ART_FILE_SAVE_AS, wx.ART_TOOLBAR, (16,16))
        saveasTool = self.toolbar.AddSimpleTool(wx.ID_SAVEAS, saveas_ico, "Convert .stk to .hdf5", "Convert .stk to .hdf5")
        self.Bind(wx.EVT_MENU, self.onSaveAsH5, saveasTool)
        self.toolbar.EnableTool(wx.ID_SAVEAS, False)
        
        help_ico = wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, wx.ART_TOOLBAR, (16,16))
        helpTool = self.toolbar.AddSimpleTool(wx.ID_HELP, help_ico, "About", "About")
        self.Bind(wx.EVT_MENU, self.onAbout, helpTool)
        
               
 
        self.toolbar.Realize()

#----------------------------------------------------------------------
    def onBrowse(self, event):
        """
        Browse for .hdf5 or .stk file
        """
        a=1
        if a==1:
        #try: 
            wildcard =  "HDF5 files (*.hdf5)|*.hdf5|STK files (*.stk)|*.stk" 
            dialog = wx.FileDialog(None, "Choose a file",
                                   style=wx.OPEN)
            dialog.SetWildcard(wildcard)
            if dialog.ShowModal() == wx.ID_OK:
                filepath = dialog.GetPath()
                self.page1.filename = dialog.GetFilename()
                                  
            basename, extension = os.path.splitext(self.page1.filename)      
            
            if extension == '.stk':
                wx.BeginBusyCursor()     
            
                if self.common.stack_loaded == 1:
                    self.new_stack_refresh()  
                    self.data_struct = data_struct.h5()
                    self.stk = data_stack.data(self.data_struct)
                    self.anlz = analyze.analyze(self.stk)                 
                self.stk.read_stk(filepath)        
                self.page1.slider_eng.SetRange(0,self.stk.n_ev-1)
                self.iev = int(self.stk.n_ev/3)
                self.page1.iev = self.iev
                self.page1.slider_eng.SetValue(self.iev)
            
                x=self.stk.n_cols
                y=self.stk.n_rows
                z=self.iev               
                self.page1.imgrgb = npy.zeros(x*y*3,dtype = "uint8")        
                self.page1.maxval = npy.amax(self.stk.imagestack)
            
                self.ix = int(x/2)
                self.iy = int(y/2)
                
                self.page1.ix = self.ix
                self.page1.iy = self.iy
                        
                self.common.stack_loaded = 1
                
                self.page1.ResetDisplaySettings()
                self.page1.loadImage()
                self.page1.textctrl.SetValue(self.page1.filename)
                

                wx.EndBusyCursor()
                
            if extension == '.hdf5':
                wx.BeginBusyCursor()     
                
                if self.common.stack_loaded == 1:
                    self.new_stack_refresh()  
                    self.data_struct = data_struct.h5()
                    self.stk = data_stack.data(self.data_struct)
                    self.anlz = analyze.analyze(self.stk)                 
            
                self.stk.read_h5(filepath)
                         
                self.page1.slider_eng.SetRange(0,self.stk.n_ev-1)
                self.iev = self.stk.n_ev/2
                self.page1.iev = self.iev
                self.page1.slider_eng.SetValue(self.iev)
            
                x=self.stk.n_cols
                y=self.stk.n_rows
                z=self.iev               
                self.page1.imgrgb = npy.zeros(x*y*3,dtype = "uint8")        
                self.page1.maxval = npy.amax(self.stk.absdata)
            
                self.ix = x/2
                self.iy = y/2
                        
                self.common.stack_loaded = 1
                self.common.i0_loaded = 1
                
                self.page1.ix = self.ix
                self.page1.iy = self.iy
                
                self.page1.ResetDisplaySettings()
                self.page1.loadImage()
                self.page1.loadSpectrum(self.ix, self.iy)
                self.page1.textctrl.SetValue(self.page1.filename)
                

                wx.EndBusyCursor()

#        except:
#
#            self.common.stack_loaded = 0 
#            self.common.i0_loaded = 0
#                               
#            wx.EndBusyCursor()
#            wx.MessageBox("Image stack not loaded.")
                   
        dialog.Destroy()
        self.refresh_widgets()
        
#----------------------------------------------------------------------
    def onSaveAsH5(self, event):

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
    def onAbout(self, event):
        AboutFrame().Show()
        return
        
#----------------------------------------------------------------------        
    def refresh_widgets(self):
        
        #page 1
        if self.common.stack_loaded == 0:
            self.page1.button_i0ffile.Disable()
            self.page1.button_i0histogram.Disable() 
            self.page1.button_save.Disable() 
            self.page1.button_addROI.Disable()
            self.page1.button_spectralROI.Disable()
            self.page1.button_resetdisplay.Disable() 
            self.page1.button_despike.Disable()   
            self.page1.button_displaycolor.Disable()
        else:
            self.page1.button_i0ffile.Enable()
            self.page1.button_i0histogram.Enable() 
            self.page1.button_save.Enable()     
            self.page1.button_addROI.Enable()  
            self.page1.button_spectralROI.Enable()
            self.page1.button_resetdisplay.Enable() 
            #self.page1.button_despike.Enable() 
            self.page1.button_displaycolor.Enable()
            
        if self.common.i0_loaded == 0:
            self.page1.button_limitev.Disable()
            self.page1.button_showi0.Disable() 
            self.page1.rb_flux.Disable()
            self.page1.rb_od.Disable()
            self.page2.button_calcpca.Disable()
            self.toolbar.EnableTool(wx.ID_SAVEAS, False)
        else:
            self.page1.button_limitev.Enable()
            self.page1.button_showi0.Enable()
            self.page1.rb_flux.Enable()
            self.page1.rb_od.Enable()   
            self.page2.button_calcpca.Enable() 
            self.toolbar.EnableTool(wx.ID_SAVEAS, True)  
            
            
        if self.common.pca_calculated == 0:      
            self.page2.button_savepca.Disable()
            self.page2.slidershow.Disable() 
            self.page3.button_calcca.Disable()
            self.page4.button_loadtspec.Disable()
            self.page4.button_addflat.Disable()
        else:
            self.page2.button_savepca.Enable()
            self.page2.slidershow.Enable()
            self.page3.button_calcca.Enable()
            self.page4.button_loadtspec.Enable()
            self.page4.button_addflat.Enable()            
            
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
        else:
            self.page4.button_removespec.Enable()
            self.page4.button_movespdown.Enable()
            self.page4.button_movespup.Enable()
            self.page4.button_save.Enable()          
                  
            
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
            
#---------------------------------------------------------------------- 
class AboutFrame(wx.Frame):

    def __init__(self):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title = "About Mantis", size=(410, 510))
        
        ico = logos.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        
        self.com = wx.GetApp().TopWindow.common         
        self.fontsize = self.com.fontsize
        
        self.SetBackgroundColour("White")
        
        
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        
        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
        
        img = logos.getMantis_logo_aboutImage()
        self.imageCtrl = wx.StaticBitmap(panel, wx.ID_ANY, wx.BitmapFromImage(img))
 
        vbox.Add(self.imageCtrl, 0, wx.ALL, 2)

        
        font3 = wx.Font(12, wx.SWISS, wx.NORMAL, wx.NORMAL)
        text3 = wx.StaticText(panel, 0, '''Mantis 1.0  
Developed by Mirna Lerotic''')
        text3.SetFont(font3)
             
        font4 = wx.Font(12, wx.SWISS, wx.NORMAL, wx.NORMAL)
        text4 = wx.StaticText(panel, 0, "www.2ndlook.co")
        text4.SetFont(font4)   
        text4.SetForegroundColour((53,159,217)) 

        vbox .Add((0,30))
  
        vbox.Add(text3,0, wx.LEFT, 50)     
        vbox.Add((0,10)) 
        vbox.Add(text4,0, wx.LEFT, 50) 
        
        
        button_close = wx.Button(panel, 0, 'Close')
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=button_close.GetId())
        vbox.Add(button_close, 0, wx.ALL | wx.ALIGN_RIGHT ,20)
        

        panel.SetSizer(vbox)
        
        vboxtop.Add(panel,0, wx.EXPAND )
        
        self.SetSizer(vboxtop)
        
#----------------------------------------------------------------------              
    def OnClose(self, evt):
        self.Close(True)
        
    
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
    splash = show_splash()
   
    time.sleep(1)
    frame = MainFrame(None, -1, 'Mantis')
    frame.Show()

    splash.Destroy()
    app.MainLoop()

if __name__ == '__main__':
    main()

