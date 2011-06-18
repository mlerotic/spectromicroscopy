from __future__ import division
import wx
import matplotlib as mtplot 
mtplot.interactive( True )
mtplot.use( 'WXAgg',warn=False )

import wxmpl as mpl 
import numpy as npy
import os.path
from mpl_toolkits.axes_grid import make_axes_locatable


import x1a_stk
import analyze
import logo_2l_32
import logo_2ndlook 

Winsizex = 1000
Winsizey = 750





#----------------------------------------------------------------------
class common:
    def __init__(self):
        self.fontsize = 8
        
        self.stack_loaded = 0
        self.i0_loaded = 0
        self.pca_calculated = 0
        self.cluster_calculated = 0
        

""" ------------------------------------------------------------------------------------------------"""
class PageCluster(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        
        self.SetBackgroundColour("White")
        
        self.stk = parent.GetParent().GetParent().stk
        self.anlz = parent.GetParent().GetParent().anlz
        
        self.selcluster = 1
        self.numclusters = 5
        
        self.wo_1st_pca = 0
        
        self.MakeColorTable()
        
        self.com = wx.GetApp().TopWindow.common      
        self.fontsize = self.com.fontsize
        
        
        #panel 1
        panel1 = wx.Panel(self, -1)
        sizer1 = wx.StaticBoxSizer(wx.StaticBox(panel1, -1, 'Cluster analysis'), wx.VERTICAL)
        self.button_calcca = wx.Button(panel1, -1, 'Calculate Clusters', (10, 10))
        self.Bind(wx.EVT_BUTTON, self.OnCalcClusters, id=self.button_calcca.GetId())   
        self.button_calcca.Disable()     
        self.button_scatterplots = wx.Button(panel1, -1, 'Show scatter plots', (10, 10))
        self.Bind(wx.EVT_BUTTON, self.OnShowScatterplots, id=self.button_scatterplots.GetId())
        self.button_scatterplots.Disable()
        self.button_savecluster = wx.Button(panel1, -1, 'Save CA Results', (10, 10))
        self.Bind(wx.EVT_BUTTON, self.OnSave, id=self.button_savecluster.GetId())
        self.button_savecluster.Disable()
        
        
        hbox11 = wx.BoxSizer(wx.HORIZONTAL)
        text1 = wx.StaticText(panel1, -1, 'Number of clusters',  style=wx.ALIGN_LEFT)
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
        self.ClusterImagePan = mpl.PlotPanel(panel2, -1, size =(3.45,3.45), cursor=False, crosshairs=True, location=False, zoom=False)                              
        mpl.EVT_POINT(self, self.ClusterImagePan.GetId(), self.OnPointClusterImage)   
        vbox2.Add(self.tc_clustercomp, 0, wx.EXPAND) 
        vbox2.Add(self.ClusterImagePan, 0)   

        panel2.SetSizer(vbox2)
        
        
        #panel 3 
        panel3 = wx.Panel(self, -1)
        vbox3 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_cluster = wx.TextCtrl(panel3, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_cluster.SetValue("Cluster ")
        hbox31 = wx.BoxSizer(wx.HORIZONTAL)
        self.ClusterIndvImagePan = mpl.PlotPanel(panel3, -1, size =(2.65,2.65), cursor=False, crosshairs=False, location=False, zoom=False)
    
        self.slidershow = wx.Slider(panel3, -1, self.selcluster, 1, 20, style=wx.SL_LEFT)   
        self.slidershow.Disable()    
        self.slidershow.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnClusterScroll, self.slidershow)
          
        hbox31.Add(self.ClusterIndvImagePan, 0)
        hbox31.Add(self.slidershow, 0,  wx.EXPAND)
        
        vbox3.Add(self.tc_cluster, 0)        
        vbox3.Add(hbox31, 0)

        panel3.SetSizer(vbox3)
        
        
        #panel 4 
        panel4 = wx.Panel(self, -1)
        vbox4 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_clustersp = wx.TextCtrl(panel4, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_clustersp.SetValue("Cluster spectrum")
        
        self.ClusterSpecPan = mpl.PlotPanel(panel4, -1, size =(5.7,3.45), cursor=False, crosshairs=False, location=False, zoom=False)
        
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
        try: 
            self.CalcClusters()
            self.calcclusters = True
            
            self.selcluster = 1
            self.slidershow.SetValue(self.selcluster)
            self.slidershow.SetMax(self.numclusters)
            
            self.showClusterImage()
            self.showClusterSpectrum()
            self.showIndvClusterImage()     
            self.com.cluster_calculated = 1       
            wx.EndBusyCursor() 
        except:
            self.com.cluster_calculated = 0
            wx.EndBusyCursor()      
            
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
        colorcl = min(self.selcluster-1,9)
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

        specplot = axes.plot(self.anlz.stack.ev,clusterspectrum)
        
        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Optical Density')
        
        
        self.ClusterSpecPan.draw()
        
        self.tc_clustersp.SetValue("Cluster " + str(self.selcluster)+ " spectrum"  ) 
        
#----------------------------------------------------------------------           
    def OnRemove1stpca(self, event):
        if self.remove1stpcacb.GetValue():
            self.wo_1st_pca = 1
        else: self.wo_1st_pca = 0
        
        #recalculate CA!
        
        
#----------------------------------------------------------------------    
    def OnSave(self, event):     
               
        fileName = wx.FileSelector('Save Plot', default_extension='png', 
                                   wildcard=('Portable Network Graphics (*.png)|*.png|' 
                                             + 'Encapsulated Postscript (*.eps)|*.eps|All files (*.*)|*.*'), 
                                              parent=self, flags=wx.SAVE|wx.OVERWRITE_PROMPT) 
   
        if not fileName: 
            return 

        path, ext = os.path.splitext(fileName) 
        ext = ext[1:].lower() 
        
       
        if ext != 'png' and ext != 'eps': 
            error_message = ( 
                  'Only the PNG and EPS image formats are supported.\n' 
                 'A file extension of `png\' or `eps\' must be used.') 
            wx.MessageBox(error_message, 'Error - plotit', 
                  parent=self, style=wx.OK|wx.ICON_ERROR) 
            return 
   
        try: 
            
            suffix = "." + ext
            
            fileName_evals = fileName[:-len(suffix)]+"_CAcimg."+ext
            self.ClusterImagePan.print_figure(fileName_evals)
            
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
                canvas.print_figure(fileName_img)
                
                
                clusterspectrum = self.anlz.clusterspectra[i, ]
                fig = mtplot.figure.Figure(figsize =(5.9,3.5))
                canvas = FigureCanvas(fig)
                fig.add_axes((0.15,0.15,0.75,0.75))
                axes = fig.gca()
                mtplot.rcParams['font.size'] = self.fontsize
                specplot = axes.plot(self.anlz.stack.ev,clusterspectrum)
        
                axes.set_xlabel('Photon Energy [eV]')
                axes.set_ylabel('Optical Density')

                fileName_spec = fileName[:-len(suffix)]+"_CAspectrum_" +str(i+1)+"."+ext
                canvas.print_figure(fileName_spec)
                
                
                
            
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
        
        colors_i = npy.linspace(0,self.maxclcolors,self.maxclcolors+2)
        
        #use black color for clusters > maxclcolors        
        colors2=['#0000FF','#FF0000','#FFFF00','#33FF33','#B366FF',
                '#FF470A','#33FFFF','#006600','#CCCC99','#993300',
                '#000000','#FFFFFF','#EEEEEE']
        
        self.clusterclrmap2=mtplot.colors.LinearSegmentedColormap.from_list('clustercm2',colors2)
     
        self.bnorm2 = mtplot.colors.BoundaryNorm(colors_i, self.clusterclrmap2.N)
        
#----------------------------------------------------------------------     
    def OnShowScatterplots(self, evt):    
        Scatterplots().Show()
        
#---------------------------------------------------------------------- 
class Scatterplots(wx.Frame):

    title = "Scatter plots"
 
    def __init__(self):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=self.title, size=(600, 550))
        
        self.SetBackgroundColour("White")
        
        ico = logo_2l_32.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White") 
        
        self.colors = wx.GetApp().TopWindow.page3.colors
        
        self.com = wx.GetApp().TopWindow.common       
        self.fontsize = self.com.fontsize
        
        self.anlz = wx.GetApp().TopWindow.page3.anlz
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

        self.ScatterPPanel = mpl.PlotPanel(panel, -1, size=(5.0, 4.0), cursor=False, crosshairs=False, location=False, zoom=False)
        
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
        self.Bind(wx.EVT_BUTTON, self.OnSaveScatt, id=self.button_savescatt.GetId())
        hbox.Add(self.button_savescatt, -1, wx.ALIGN_RIGHT|wx.RIGHT, 20)
        
        button_close = wx.Button(panel, -1, 'Close')
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
          
               
        fileName = wx.FileSelector('Save Plot', default_extension='png', 
                                   wildcard=('Portable Network Graphics (*.png)|*.png|' 
                                             + 'Encapsulated Postscript (*.eps)|*.eps|All files (*.*)|*.*'), 
                                              parent=self, flags=wx.SAVE|wx.OVERWRITE_PROMPT) 
        
   
        if not fileName: 
            return 

        path, ext = os.path.splitext(fileName) 
        ext = ext[1:].lower() 
        
       
        if ext != 'png' and ext != 'eps': 
            error_message = ( 
                  'Only the PNG and EPS image formats are supported.\n' 
                 'A file extension of `png\' or `eps\' must be used.') 
            wx.MessageBox(error_message, 'Error - plotit', 
                  parent=self, style=wx.OK|wx.ICON_ERROR) 
            return 
   
        try: 
            
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
            canvas.print_figure(fileName_sct)
                

            
        except IOError, e:
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
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        
        self.SetBackgroundColour("White")
        
        self.stk = parent.GetParent().GetParent().stk       
        self.anlz = parent.GetParent().GetParent().anlz
        
        self.selpca = 1       
        self.numsigpca = 2
        
        self.com = wx.GetApp().TopWindow.common   
        self.fontsize = self.com.fontsize
        
  
        #panel 1        
        panel1 = wx.Panel(self, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_PCAcomp = wx.TextCtrl(panel1, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_PCAcomp.SetValue("PCA component ")
        
        hbox11 = wx.BoxSizer(wx.HORIZONTAL)
   
        self.PCAImagePan = mpl.PlotPanel(panel1, -1, size =(3.9,3.4), cursor=False, crosshairs=False, location=False, zoom=False)
                              
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
        self.Bind(wx.EVT_BUTTON, self.OnCalcPCA, id=self.button_calcpca.GetId())     
        self.button_calcpca.Disable()   
        sizer1.Add(self.button_calcpca, 0, wx.EXPAND)
        self.button_savepca = wx.Button(panel2, -1, 'Save PCA Results', (10, 10))
        self.Bind(wx.EVT_BUTTON, self.OnSave, id=self.button_savepca.GetId())
        self.button_savepca.Disable()
        sizer1.Add(self.button_savepca, 0, wx.EXPAND)
        
        hbox21 = wx.BoxSizer(wx.HORIZONTAL)
        text1 = wx.StaticText(panel2, -1, 'Number of significant components',  style=wx.ALIGN_LEFT)
        self.npcaspin = wx.SpinCtrl(panel2, -1, '',  size= (60, -1), style=wx.ALIGN_LEFT)
        self.npcaspin.SetRange(1,20)
        self.Bind(wx.EVT_SPINCTRL, self.OnNPCAspin, self.npcaspin)
        hbox21.Add(text1, 0, wx.TOP, 20)
        hbox21.Add((10,0))
        hbox21.Add(self.npcaspin, 0, wx.TOP, 15)            
        sizer1.Add(hbox21, 0, wx.EXPAND)    
              
        hbox22 = wx.BoxSizer(wx.HORIZONTAL)
        text2 = wx.StaticText(panel2, -1, 'Cumulative variance', style=wx.ALIGN_LEFT)
        self.vartc = wx.StaticText(panel2, -1, '0%',  style=wx.ALIGN_LEFT)
        hbox22.Add(text2, 0)
        hbox22.Add(self.vartc, 0, wx.LEFT , 10)

        sizer1.Add(hbox22, 0)      
        
        vbox2.Add(sizer1)

        panel2.SetSizer(vbox2)

        
        #panel 3
        panel3 = wx.Panel(self, -1)
        
        self.PCASpecPan = mpl.PlotPanel(panel3, -1, size =(5.4,3.4), cursor=False, crosshairs=False, location=False, zoom=False)
             
        vbox3 = wx.BoxSizer(wx.VERTICAL)
        self.text_pcaspec = wx.TextCtrl(panel3, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.text_pcaspec.SetValue("PCA spectrum ")        
        vbox3.Add(self.text_pcaspec, 0)
        vbox3.Add(self.PCASpecPan, 0)        
        panel3.SetSizer(vbox3)
        
        
        #panel 4
        panel4 = wx.Panel(self, -1)
             
        self.PCAEvalsPan = mpl.PlotPanel(panel4, -1, size =(5.4,2.6), cursor=False, crosshairs=False, location=False, zoom=False)
        
        vbox4 = wx.BoxSizer(wx.VERTICAL)
        text4 = wx.TextCtrl(panel4, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
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
 
        
        self.anlz.setdata(self.stk)
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
               
        fileName = wx.FileSelector('Save Plot', default_extension='png', 
                                   wildcard=('Portable Network Graphics (*.png)|*.png|' 
                                             + 'Encapsulated Postscript (*.eps)|*.eps|All files (*.*)|*.*'), 
                                              parent=self, flags=wx.SAVE|wx.OVERWRITE_PROMPT) 
   
        if not fileName: 
            return 

        path, ext = os.path.splitext(fileName) 
        ext = ext[1:].lower() 
        
       
        if ext != 'png' and ext != 'eps': 
            error_message = ( 
                  'Only the PNG and EPS image formats are supported.\n' 
                 'A file extension of `png\' or `eps\' must be used.') 
            wx.MessageBox(error_message, 'Error - plotit', 
                  parent=self, style=wx.OK|wx.ICON_ERROR) 
            return 
   
        try: 
            
            suffix = "." + ext
            fileName_evals = fileName[:-len(suffix)]+"_PCAevals."+ext
            
            self.PCAEvalsPan.print_figure(fileName_evals)
            
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
                canvas.print_figure(fileName_img)
                
                
                
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
                canvas.print_figure(fileName_spec)
                
                
                
            
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
        #self.PCAImagePan.Refresh()        
        
""" ------------------------------------------------------------------------------------------------"""
class PageStack(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        
        self.SetBackgroundColour("White")
        
        self.stk = parent.GetParent().GetParent().stk
        self.filename = " "
       
        self.ix = 0
        self.iy = 0
        self.iev = 50  
        self.sel = 50
        self.showflux = True
        
        self.com = wx.GetApp().TopWindow.common      
        self.fontsize = self.com.fontsize
 

        vbox = wx.BoxSizer(wx.VERTICAL)
        hboxT = wx.BoxSizer(wx.HORIZONTAL)
        hboxB = wx.BoxSizer(wx.HORIZONTAL)
    
        #panel 1        
        panel1 = wx.Panel(self, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_imageeng = wx.TextCtrl(panel1, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_imageeng.SetValue("Image at energy: ")
        hbox11 = wx.BoxSizer(wx.HORIZONTAL)
   
        self.AbsImagePanel = mpl.PlotPanel(panel1, -1, size =(3.5,3.5), cursor=False, crosshairs=True, location=False, zoom=False)
        mpl.EVT_POINT(panel1, self.AbsImagePanel.GetId(), self.on_point_absimage)
                              
        self.slider = wx.Slider(panel1, -1, self.sel, 0, 100, style=wx.SL_LEFT )        
        self.slider.SetFocus()
        self.Bind(wx.EVT_SCROLL, self.OnScroll)

        hbox11.Add(self.AbsImagePanel, 0)
        hbox11.Add(self.slider, 0,  wx.EXPAND)
        
        vbox1.Add(self.tc_imageeng,1, wx.LEFT | wx.TOP | wx.EXPAND, 20)        
        vbox1.Add(hbox11, 0,  wx.LEFT, 20)

        panel1.SetSizer(vbox1)
     
        #panel 2

        panel2 = wx.Panel(self, -1)
        vbox2 = wx.BoxSizer(wx.VERTICAL)
        
        self.tc_spec = wx.TextCtrl(panel2, 0, style=wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.tc_spec.SetValue("Spectrum at point: ")
          
        
        self.SpectrumPanel = mpl.PlotPanel(panel2, -1, size=(5.75, 3.5), cursor=False, crosshairs=False, location=False, zoom=False)

        vbox2.Add(self.tc_spec, 1, wx.LEFT | wx.TOP | wx.EXPAND, 20)       
        vbox2.Add(self.SpectrumPanel, 0, wx.LEFT , 20)
        
        panel2.SetSizer(vbox2)
        
        
        #panel 3
        panel3 = wx.Panel(self, -1)
        sizer1 = wx.StaticBoxSizer(wx.StaticBox(panel3, -1, 'Preprocess'), wx.VERTICAL)
        self.button_limitev = wx.Button(panel3, -1, 'Limit energy range')
        self.Bind(wx.EVT_BUTTON, self.OnLimitEv, id=self.button_limitev.GetId())
        self.button_limitev.Disable()
        sizer1.Add(self.button_limitev, 0, )
        self.button_i0ffile = wx.Button(panel3, -1, 'I0 from file')
        self.Bind(wx.EVT_BUTTON, self.OnI0FFile, id=self.button_i0ffile.GetId())
        self.button_i0ffile.Disable()
        sizer1.Add(self.button_i0ffile, 0, wx.EXPAND)
        self.button_i0histogram = wx.Button(panel3, -1, 'I0 from histogram')
        self.Bind(wx.EVT_BUTTON, self.OnI0histogram, id=self.button_i0histogram.GetId())   
        self.button_i0histogram.Disable()     
        sizer1.Add(self.button_i0histogram, 0, wx.EXPAND)
        self.button_showi0 = wx.Button(panel3, -1, 'Show I0')
        self.Bind(wx.EVT_BUTTON, self.OnShowI0, id=self.button_showi0.GetId())   
        self.button_showi0.Disable()
        sizer1.Add(self.button_showi0, 0, wx.EXPAND)
        self.button_save = wx.Button(panel3, -1, 'Save')
        self.Bind(wx.EVT_BUTTON, self.OnSave, id=self.button_save.GetId())
        self.button_save.Disable()          
        sizer1.Add(self.button_save, 0, wx.EXPAND)
        panel3.SetSizer(sizer1)
        

        
        #panel 4
        panel4 = wx.Panel(self, -1)
        sizer4 = wx.StaticBoxSizer(wx.StaticBox(panel4, -1, 'Display', size =(500,-1)), orient=wx.VERTICAL)
        vbox3 = wx.BoxSizer(wx.VERTICAL)
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        sizer11 = wx.StaticBoxSizer(wx.StaticBox(panel4, -1, 'File', size =(500,-1)), orient=wx.VERTICAL)
        self.textctrl = wx.TextCtrl(panel4, -1, style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        sizer11.Add(self.textctrl, 1, wx.EXPAND)
        hbox1.Add(sizer11, 1, wx.EXPAND)
        
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        sizer21 = wx.StaticBoxSizer(wx.StaticBox(panel4, -1, 'Image', size =(500,-1)),  orient=wx.VERTICAL)
        self.rb_flux = wx.RadioButton(panel4, -1, 'Flux', style=wx.RB_GROUP)
        self.rb_od = wx.RadioButton(panel4, -1, 'Optical Density')
        self.Bind(wx.EVT_RADIOBUTTON, self.onrb_fluxod, id=self.rb_flux.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.onrb_fluxod, id=self.rb_od.GetId())
        
        self.rb_flux.Disable()
        self.rb_od.Disable()

        sizer21.Add(self.rb_flux)
        sizer21.Add(self.rb_od)
        hbox2.Add(sizer21, 1, wx.EXPAND)

        
        vbox3.Add(hbox1, 1, wx.LEFT | wx.RIGHT | wx.EXPAND, 10)
        vbox3.Add(hbox2, 1, wx.LEFT | wx.RIGHT | wx.EXPAND, 10)
        sizer4.Add(vbox3,1, wx.EXPAND)

        
        panel4.SetSizer(sizer4)
        

        hboxB.Add(panel1, 0, wx.BOTTOM | wx.TOP, 9)
        hboxB.Add(panel2, 0, wx.BOTTOM | wx.TOP, 9)
        hboxT.Add((10,0))        
        hboxT.Add(panel3, 1, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 10)
        hboxT.Add(panel4, 3.5, wx.LEFT | wx.RIGHT |wx.TOP | wx.EXPAND,10)
        
        vbox.Add(hboxT, 0, wx.ALL, 5)
        vbox.Add(hboxB, 0,  wx.ALL, 5)
  
        self.SetSizer(vbox) 
        

      

#----------------------------------------------------------------------        
    def loadImage(self):
        
        self.image = 0         
        if self.showflux:  
            #Show flux image      
            self.image = self.stk.absdata[:,:,self.iev]#.copy() 
        else:
            #Show OD image
            self.image = self.stk.od3d[:,:,self.iev]#.copy() 
              
               
        fig = self.AbsImagePanel.get_figure()
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))
        axes = fig.gca()
      
        im = axes.imshow(self.image, cmap=mtplot.cm.get_cmap("gray"))   
        axes.axis("off")  
        self.AbsImagePanel.draw()
        
        #self.textctrl.SetValue(self.filename+"\n"+str(self.stk.ev[self.iev])+" eV")
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
        
        self.SpectrumPanel.draw()
        
        self.tc_spec.SetValue("Spectrum at point: [" +str(xpos)+", " + str(ypos)+"] ")

        
#----------------------------------------------------------------------            
    def OnScroll(self, event):
        self.sel = event.GetInt()
        self.iev = self.sel
        if self.com.stack_loaded == 1:
            self.loadImage()
                   
#----------------------------------------------------------------------  
    def on_point_absimage(self, evt):
        x = evt.xdata
        y = evt.ydata
        
        if self.com.i0_loaded == 1:      
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
        
#----------------------------------------------------------------------
        
    def OnI0FFile(self, event):

        try: 
            wildcard = "STK files (*.xas)|*.xas"
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
            self.com.i0_loaded = 0        
            wx.EndBusyCursor()  
            wx.MessageBox("I0 file not loaded.")
            
                       
        dialog.Destroy()
                          
        wx.GetApp().TopWindow.refresh_widgets()

        
        
#----------------------------------------------------------------------       
    def OnI0histogram(self, event):    
        ShowHistogram().Show()
         

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
                                             + 'Encapsulated Postscript (*.eps)|*.eps|All files (*.*)|*.*'), 
                                              parent=self, flags=wx.SAVE|wx.OVERWRITE_PROMPT) 
   
        if not fileName: 
            return 

        path, ext = os.path.splitext(fileName) 
        ext = ext[1:].lower() 
        
       
        if ext != 'png' and ext != 'eps': 
            error_message = ( 
                  'Only the PNG and EPS image formats are supported.\n' 
                 'A file extension of `png\' or `eps\' must be used.') 
            wx.MessageBox(error_message, 'Error - plotit', 
                  parent=self, style=wx.OK|wx.ICON_ERROR) 
            return 
   
        try: 
            
            suffix = "." + ext
            fileName_spec = fileName[:-len(suffix)]+"_spectrum."+ext
          
            self.SpectrumPanel.print_figure(fileName_spec)
            
            fileName_img = fileName[:-len(suffix)]+"_" +str(self.stk.ev[self.iev])+"eV."+ext
            
            self.AbsImagePanel.print_figure(fileName_img)
            
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
            
        self.loadImage()

        
 
#----------------------------------------------------------------------    
        
    def OnLimitEv(self, evt):    
        LimitEv().Show()
        
#---------------------------------------------------------------------- 
class ShowHistogram(wx.Frame):

    title = "Histogram"

    def __init__(self):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=self.title, size=(630, 700))
               
        ico = logo_2l_32.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White") 
        
        self.stack = wx.GetApp().TopWindow.page1.stk
        
        self.com = wx.GetApp().TopWindow.common         
        self.fontsize = self.com.fontsize
        
        self.stack.calc_histogram()
        averagefluxmax = npy.max(self.stack.histogram)
        self.histmin = 0.98*averagefluxmax
        self.histmax = averagefluxmax
        
    
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
               
        self.HistogramPanel = mpl.PlotPanel(panel, -1, size=(6.0, 3.7), cursor=False, crosshairs=False, location=False, zoom=False)
        
        mpl.EVT_SELECTION(panel, self.HistogramPanel.GetId(), self.OnSelection)

        vbox.Add(self.HistogramPanel, 0, wx.ALL, 20)
        
       
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        sizer2 = wx.StaticBoxSizer(wx.StaticBox(panel, -1, 'I0 pixels'), orient=wx.VERTICAL)
        self.textctrl = wx.TextCtrl(panel, -1, size = (565, 20), style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH|wx.BORDER_NONE)
        self.textctrl.SetValue('Selection: [ '+str(self.histmin) + ' kHz, '+ str(self.histmax)+' kHz ]' )
        self.textctrl.SetBackgroundColour("White")
        sizer2.Add(self.textctrl, 0)
        

        self.AbsImagePanel = mpl.PlotPanel(panel, -1, size =(1.5,1.5), cursor=False, crosshairs=False, location=False, zoom=False)
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

    def __init__(self):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=self.title, size=(630, 560))
               
        ico = logo_2l_32.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White") 
        
        self.stack = wx.GetApp().TopWindow.page1.stk
        
        self.com = wx.GetApp().TopWindow.common         
        self.fontsize = self.com.fontsize
        
        self.evlimited = 0
        self.limitevmin = 0
        self.limitevmax = self.stack.original_n_ev-1
        
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
               
        self.SpectrumPanel = mpl.PlotPanel(panel, -1, size=(6.0, 3.7), cursor=False, crosshairs=False, location=False, zoom=False)
        
        mpl.EVT_SELECTION(panel, self.SpectrumPanel.GetId(), self.OnSelection)

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
        wx.GetApp().TopWindow.page1.slider.SetRange(0,self.stack.n_ev-1)
        wx.GetApp().TopWindow.page1.iev = self.stack.n_ev/2
        wx.GetApp().TopWindow.page1.slider.SetValue(wx.GetApp().TopWindow.page1.iev)
        
        wx.GetApp().TopWindow.page1.loadSpectrum(wx.GetApp().TopWindow.page1.ix, wx.GetApp().TopWindow.page1.iy)
        wx.GetApp().TopWindow.page1.loadImage()
        
        self.Close(True)
        
#---------------------------------------------------------------------- 
    def OnCancel(self, evt):
        self.Close(True)
        
    
#---------------------------------------------------------------------- 
class PlotFrame(wx.Frame):

    def __init__(self, datax, datay):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title = "I0 data", size=(630, 500))
        
        ico = logo_2l_32.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        self.SetBackgroundColour("White")
        
        
        self.com = wx.GetApp().TopWindow.common         
        self.fontsize = self.com.fontsize
        
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
               
        self.PlotPanel = mpl.PlotPanel(panel, -1, size=(6.0, 3.7), cursor=False, crosshairs=False, location=False, zoom=False)

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
        
            
""" ------------------------------------------------------------------------------------------------"""
class MainFrame(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, size=(Winsizex, Winsizey))
        
      
        ico = logo_2l_32.getlogo_2l_32Icon()
        self.SetIcon(ico)

        
        self.initToolbar()
        
        self.stk = x1a_stk.x1astk()
        self.anlz = analyze.analyze()
        
        self.common = common()
                      

        # Here we create a panel and a notebook on the panel
        p = wx.Panel(self)
        nb = wx.Notebook(p, style=wx.BORDER_STATIC)

        # create the page windows as children of the notebook
        self.page1 = PageStack(nb)
        self.page2 = PagePCA(nb)
        self.page3 = PageCluster(nb)

        # add the pages to the notebook with the label to show on the tab
        nb.AddPage(self.page1, "Image Stack")
        nb.AddPage(self.page2, "PCA")
        nb.AddPage(self.page3, "Cluster Analysis")

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
        #self.Bind(wx.EVT_MENU, self.onBrowse, openTool)     
        
        
        help_ico = wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, wx.ART_TOOLBAR, (16,16))
        helpTool = self.toolbar.AddSimpleTool(wx.ID_HELP, help_ico, "About", "About")
        self.Bind(wx.EVT_MENU, self.onAbout, helpTool)
        
               
 
        self.toolbar.Realize()

#----------------------------------------------------------------------
    def onBrowse(self, event):
        """
        Browse for .stk file
        """
        try: 
            wildcard = "STK files (*.stk)|*.stk"
            dialog = wx.FileDialog(None, "Choose a file",
                                    wildcard=wildcard,
                                    style=wx.OPEN)
            if dialog.ShowModal() == wx.ID_OK:
                            filepath = dialog.GetPath()
                            self.page1.filename = dialog.GetFilename()
                                                        
            wx.BeginBusyCursor()     
            
            if self.common.stack_loaded == 1:
                self.new_stack_refresh()  
                self.stk = x1a_stk.x1astk() 
                self.anlz = analyze.analyze()                   
            
            self.stk.read_stk(filepath)
                         
            self.page1.slider.SetRange(0,self.stk.n_ev-1)
            self.iev = self.stk.n_ev/2
            self.page1.iev = self.iev
            self.page1.slider.SetValue(self.iev)
            
            x=self.stk.n_cols
            y=self.stk.n_rows
            z=self.iev               
            self.page1.imgrgb = npy.zeros(x*y*3,dtype = "uint8")        
            self.page1.maxval = npy.amax(self.stk.imagestack)
            
            self.ix = x/2
            self.iy = y/2
                        
            self.common.stack_loaded = 1
                
            self.page1.loadImage()
            self.page1.textctrl.SetValue(self.page1.filename)

            wx.EndBusyCursor()
            
        except:
            self.common.stack_loaded = 0 
            self.common.i0_loaded = 0
                               
            wx.EndBusyCursor()
            wx.MessageBox("Image stack not loaded.")
                   
        dialog.Destroy()
        self.refresh_widgets()
       
#----------------------------------------------------------------------
    def onAbout(self, event):
        AboutFrame().Show()
        
#----------------------------------------------------------------------        
    def refresh_widgets(self):
        
        #page 1
        if self.common.stack_loaded == 0:
            self.page1.button_i0ffile.Disable()
            self.page1.button_i0histogram.Disable() 
            self.page1.button_save.Disable() 
        else:
            self.page1.button_i0ffile.Enable()
            self.page1.button_i0histogram.Enable() 
            self.page1.button_save.Enable()             
            
        if self.common.i0_loaded == 0:
            self.page1.button_limitev.Disable()
            self.page1.button_showi0.Disable() 
            self.page1.rb_flux.Disable()
            self.page1.rb_od.Disable()
            self.page2.button_calcpca.Disable()
        else:
            self.page1.button_limitev.Enable()
            self.page1.button_showi0.Enable()
            self.page1.rb_flux.Enable()
            self.page1.rb_od.Enable()   
            self.page2.button_calcpca.Enable()   
            
            
        if self.common.pca_calculated == 0:      
            self.page2.button_savepca.Disable()
            self.page2.slidershow.Disable() 
            self.page3.button_calcca.Disable()
        else:
            self.page2.button_savepca.Enable()
            self.page2.slidershow.Enable()
            self.page3.button_calcca.Enable()
            
            
        if self.common.cluster_calculated == 0:   
            self.page3.button_scatterplots.Disable()
            self.page3.button_savecluster.Disable()
            self.page3.slidershow.Disable()
        else:
            self.page3.button_scatterplots.Enable()
            self.page3.button_savecluster.Enable()  
            self.page3.slidershow.Enable()          
            
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
        
        self.page3.selcluster = 1
        self.page3.slidershow.SetValue(self.page3.selcluster)
        self.page3.numclusters = 5
        self.page3.nclusterspin.SetValue(self.page3.numclusters)
        self.page3.tc_cluster.SetValue("Cluster ")
        self.page3.tc_clustersp.SetValue("Cluster spectrum")
        self.page3.wo_1st_pca = 0
        self.page3.remove1stpcacb.SetValue(False)
            
#---------------------------------------------------------------------- 
class AboutFrame(wx.Frame):

    def __init__(self):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title = "About", size=(350, 400))
        
        ico = logo_2l_32.getlogo_2l_32Icon()
        self.SetIcon(ico)
        
        
        self.com = wx.GetApp().TopWindow.common         
        self.fontsize = self.com.fontsize
        
        self.SetBackgroundColour("White")
        
        
        vboxtop = wx.BoxSizer(wx.VERTICAL)
        
        
        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
        
        img = logo_2ndlook.getlogo_2ndlookImage()              
        self.imageCtrl = wx.StaticBitmap(panel, wx.ID_ANY, wx.BitmapFromImage(img))
 
        vbox.Add(self.imageCtrl, 0, wx.ALL, 20)

        
        font1 = wx.Font(12, wx.SWISS, wx.NORMAL, wx.NORMAL)
        text1 = wx.StaticText(panel, 0, "Mantis")
        text1.SetFont(font1)
        
        font2 = wx.Font(10, wx.SWISS, wx.NORMAL, wx.LIGHT)
        text2 = wx.StaticText(panel, 0, "version 1.0")
        text2.SetFont(font2)
        
        font3 = wx.Font(12, wx.SWISS, wx.NORMAL, wx.NORMAL)
        text3 = wx.StaticText(panel, 0, "Developed by Mirna Lerotic")
        text3.SetFont(font3)
             
        font4 = wx.Font(12, wx.SWISS, wx.NORMAL, wx.NORMAL)
        text4 = wx.StaticText(panel, 0, "www.2ndlook.co")
        text4.SetFont(font4)   
        text4.SetForegroundColour((53,159,217)) 

        vbox .Add((0,10))
        vbox.Add(text1,0, wx.LEFT, 40)
        vbox.Add(text2,0, wx.LEFT, 40)  
        vbox.Add((0,25))      
        vbox.Add(text3,0, wx.LEFT, 40)     
        vbox.Add((0,10)) 
        vbox.Add(text4,0, wx.LEFT, 120) 
        
        
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

if __name__ == '__main__':
    app = wx.App()
    MainFrame(None, -1, 'Mantis')
    app.MainLoop()
