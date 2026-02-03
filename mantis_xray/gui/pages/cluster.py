
import os
import numpy as np
from PIL import Image
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtCore import Qt
import matplotlib
import matplotlib.cm
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
# from mpl_toolkits.axes_grid1 import make_axes_locatable

from ...core.constants import PlotH, PlotW, ImgDpi
from ..widgets import SpecFig, ImgFig
from ..dialogs.scatter import Scatterplots
from ..dialogs.save import SaveWinP3

class PageCluster(QtWidgets.QWidget):
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
        sizer1 = QtWidgets.QGroupBox('Cluster analysis')
        vbox1 = QtWidgets.QVBoxLayout()


        self.button_calcca = QtWidgets.QPushButton('Calculate Clusters')
        self.button_calcca.clicked.connect( self.OnCalcClusters)
        self.button_calcca.setEnabled(False)
        vbox1.addWidget(self.button_calcca)
        self.button_scatterplots = QtWidgets.QPushButton('Show scatter plots...')
        self.button_scatterplots.clicked.connect( self.OnShowScatterplots)
        self.button_scatterplots.setEnabled(False)
        self.button_savecluster = QtWidgets.QPushButton('Save CA Results...')
        self.button_savecluster.clicked.connect( self.OnSave)
        self.button_savecluster.setEnabled(False)

        vbox1.addStretch(1)


        hbox11 = QtWidgets.QHBoxLayout()
        text1 = QtWidgets.QLabel(self)
        text1.setText('Number of clusters')
        hbox11.addWidget(text1)
        self.nclusterspin = QtWidgets.QSpinBox()
        self.nclusterspin.setRange(2,20)
        self.nclusterspin.setValue(self.init_nclusters)
        self.nclusterspin.valueChanged[int].connect(self.OnNClusterspin)
        hbox11.addWidget(text1)
        hbox11.addWidget(self.nclusterspin)

        vbox1.addLayout(hbox11)

        hbox12 = QtWidgets.QHBoxLayout()
        text1a = QtWidgets.QLabel(self)
        text1a.setText("Number of clusters found")
        hbox12.addWidget(text1a)

        self.ntc_clusters_found = QtWidgets.QLabel(self)
        self.ntc_clusters_found.setText(str(self.numclusters))
        hbox12.addWidget(self.ntc_clusters_found)

        vbox1.addLayout(hbox12)


        hbox13 = QtWidgets.QHBoxLayout()
        self.remove1stpcacb = QtWidgets.QCheckBox('Reduce thickness effects', self)
        self.remove1stpcacb.stateChanged.connect(self.OnRemove1stpca)
        hbox13.addWidget(self.remove1stpcacb)

        vbox1.addLayout(hbox13)

        hbox14 = QtWidgets.QHBoxLayout()
        self.cb_splitclusters = QtWidgets.QCheckBox('Divide clusters with large Sigma', self)
        self.cb_splitclusters.stateChanged.connect(self.OnSplitClusters)
        hbox14.addWidget(self.cb_splitclusters)

        vbox1.addLayout(hbox14)

        hbox14a = QtWidgets.QHBoxLayout()
        tc1 = QtWidgets.QLabel(self)
        tc1.setText("PC scaling factor")
        hbox14a.addWidget(tc1)
        self.ntc_pcscaling = QtWidgets.QLineEdit(self)
        self.ntc_pcscaling.setFixedWidth(65)
        self.ntc_pcscaling.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.ntc_pcscaling.setAlignment(QtCore.Qt.AlignRight)

        self.ntc_pcscaling.setText(str(self.pcscalingfactor))
        hbox14a.addWidget(self.ntc_pcscaling)
        hbox14a.addStretch(1)
        vbox1.addLayout(hbox14a)

        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.HLine)
        line.setFrameShadow(QtWidgets.QFrame.Sunken)


        vbox1.addStretch(1)
        vbox1.addWidget(line)
        vbox1.addStretch(1)
        vbox1.addWidget(self.button_scatterplots)
        vbox1.addWidget(self.button_savecluster)


        sizer1.setLayout(vbox1)



        #panel 2
        vbox2 = QtWidgets.QVBoxLayout()

        tc_clustercomp = QtWidgets.QLabel(self)
        tc_clustercomp.setText("Composite cluster image")
        vbox2.addWidget(tc_clustercomp)

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.clusterimgfig = Figure((PlotH, PlotH))
        self.ClusterImagePan = FigureCanvas(self.clusterimgfig)
        self.ClusterImagePan.mpl_connect('button_press_event', self.OnPointClusterImage)
        self.ClusterImagePan.setParent(self)
        fbox.addWidget(self.ClusterImagePan)
        frame.setLayout(fbox)
        vbox2.addWidget(frame)

        #panel 3
        vbox3 = QtWidgets.QVBoxLayout()
        fgs = QtWidgets.QGridLayout()

        self.tc_cluster = QtWidgets.QLabel(self)
        self.tc_cluster.setText("Cluster ")
        fgs.addWidget(self.tc_cluster, 0, 0, QtCore .Qt. AlignLeft)

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        #self.clusterindvimgfig = Figure((PlotH*0.73, PlotH*0.73))
        self.clusterindvimgfig = Figure((PlotH, PlotH))
        self.ClusterIndvImagePan = FigureCanvas(self.clusterindvimgfig)
        self.ClusterIndvImagePan.setParent(self)
        fbox.addWidget(self.ClusterIndvImagePan)
        frame.setLayout(fbox)
        fgs.addWidget(frame, 1, 0, QtCore .Qt. AlignLeft)

        self.slidershow = QtWidgets.QScrollBar(QtCore.Qt.Vertical)
        self.slidershow.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slidershow.setEnabled(False)
        self.slidershow.valueChanged[int].connect(self.OnClusterScroll)
        self.slidershow.setRange(1, 20)
        fgs.addWidget(self.slidershow, 1, 1, QtCore .Qt. AlignLeft)


        text3 = QtWidgets.QLabel(self)
        text3.setText('Cluster Error Map')
        fgs.addWidget(text3, 0, 2, QtCore .Qt. AlignLeft)
        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()
        self.clusterdistmapfig = Figure((PlotH, PlotH))
        self.ClusterDistMapPan = FigureCanvas(self.clusterdistmapfig)
        self.ClusterDistMapPan.setParent(self)
        fbox.addWidget(self.ClusterDistMapPan)
        frame.setLayout(fbox)
        fgs.addWidget(frame, 1, 2, QtCore .Qt. AlignLeft)

        vbox3.addLayout(fgs)




        #panel 4
        vbox4 = QtWidgets.QVBoxLayout()

        self.tc_clustersp = QtWidgets.QLabel(self)
        self.tc_clustersp.setText("Cluster spectrum")
        vbox4.addWidget(self.tc_clustersp)

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.clusterspecfig = Figure((PlotW, PlotH))
        self.ClusterSpecPan = FigureCanvas(self.clusterspecfig)
        self.ClusterSpecPan.setParent(self)

        fbox.addWidget(self.ClusterSpecPan)
        frame.setLayout(fbox)
        vbox4.addWidget(frame)


        #panel 5
        sizer5 = QtWidgets.QGroupBox('Display')
        vbox5 = QtWidgets.QVBoxLayout()

        hbox51 = QtWidgets.QHBoxLayout()
        self.showallspectracb = QtWidgets.QCheckBox('Show all spectra', self)
        self.showallspectracb.stateChanged.connect(self.OnShowallspectra)
        hbox51.addWidget(self.showallspectracb)

        vbox5.addLayout(hbox51)

        sizer5.setLayout(vbox5)



        vboxtop = QtWidgets.QVBoxLayout()
        hboxtopL = QtWidgets.QHBoxLayout()
        vboxtopL = QtWidgets.QVBoxLayout()
        vboxtopL.addWidget(sizer1)
        vboxtopL.addWidget(sizer5)
        hboxtopL.addLayout(vboxtopL)
        hboxtopL.addStretch(1)

        gridsizertop = QtWidgets.QGridLayout()
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

        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
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
            QtWidgets.QApplication.restoreOverrideCursor()

        except Exception as e:
            self.com.cluster_calculated = 0
            QtWidgets.QApplication.restoreOverrideCursor()
            print(e)

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

        if (x == None) or (y == None):
            return

            try:
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


                self.selcluster = self.anlz.cluster_indices[self.ix,self.iy] + 1

                self.slidershow.setValue(self.selcluster)

                self.showClusterSpectrum()
                self.showIndvClusterImage()
            except Exception:
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

        indvclusterimage = self.anlz.cluster_indices.copy()
        indvclusterimage[indvclusterimage!=self.selcluster-1] = 20.

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



        im = axes.imshow(np.rot90(mapimage), cmap=matplotlib.colormaps["gray"])

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
#----------------------------------------------------------------------
    def OnSave(self, event):

        savewin = SaveWinP3(self.window())
        savewin.show()


#----------------------------------------------------------------------
    def OnShowScatterplots(self, event):
        scattplwin = Scatterplots(self.window(), self.com, self.anlz)
        scattplwin.show()
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
                #ToDo: Recently throws an error. Possible conflict in module PIL. "Cannot handle this data type: (1, 1), <i8"
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

        except IOError as e:
            if e.strerror:
                err = e.strerror
            else:
                err = e
            print(e)
            QtWidgets.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err)



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
            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))

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
                            axes = fig.add_subplot(int(nplotsrows),int(2), int(pplot))
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

            QtWidgets.QApplication.restoreOverrideCursor()


        except IOError as e:
            QtWidgets.QApplication.restoreOverrideCursor()
            if e.strerror:
                err = e.strerror
            else:
                err = e

            QtWidgets.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err)
