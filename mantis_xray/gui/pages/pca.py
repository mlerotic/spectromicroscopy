
import os
import numpy as np
from PIL import Image
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtCore import Qt
import matplotlib
import matplotlib.cm
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable

from ..dialogs.save import SaveWinP2
from ...core.constants import PlotH, PlotW

class PagePCA(QtWidgets.QWidget):
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
        vbox1 = QtWidgets.QVBoxLayout()

        self.tc_PCAcomp = QtWidgets.QLabel(self)
        self.tc_PCAcomp.setText("PCA component ")
        vbox1.addWidget(self.tc_PCAcomp)

        hbox11 = QtWidgets.QHBoxLayout()

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()
        self.pcaimgfig = Figure((PlotH*1.10, PlotH))
        self.PCAImagePan = FigureCanvas(self.pcaimgfig)
        self.PCAImagePan.setParent(self)

        fbox.addWidget(self.PCAImagePan)
        frame.setLayout(fbox)
        hbox11.addWidget(frame)

        self.slidershow = QtWidgets.QScrollBar(QtCore.Qt.Vertical)
        self.slidershow.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slidershow.setRange(1,20)
        self.slidershow.setEnabled(False)
        self.slidershow.valueChanged[int].connect(self.OnPCAScroll)


        hbox11.addWidget(self.slidershow)
        vbox1.addLayout(hbox11)


        self.slider_theta = QtWidgets.QScrollBar(QtCore.Qt.Horizontal)
        self.slider_theta.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_theta.valueChanged[int].connect(self.OnScrollTheta)
        self.slider_theta.setRange(0, 100)
        self.slider_theta.setMinimumWidth(250)
        self.slider_theta.setVisible(False)
        self.tc_imagetheta = QtWidgets.QLabel(self)
        self.tc_imagetheta.setText("4D Data Angle: ")
        self.tc_imagetheta.setVisible(False)
        hbox51 = QtWidgets.QHBoxLayout()
        hbox51.addWidget(self.tc_imagetheta)
        hbox51.addStretch(1)
        hbox51.addWidget(self.slider_theta)
        hbox51.addStretch(1)
        vbox1.addLayout(hbox51)



        #panel 2
        vbox2 = QtWidgets.QVBoxLayout()
        vbox2.setContentsMargins(20,20,20,20)
        sizer2 = QtWidgets.QGroupBox('PCA')
        vbox21 = QtWidgets.QVBoxLayout()

        self.button_calcpca = QtWidgets.QPushButton('Calculate PCA')
        self.button_calcpca.clicked.connect( self.OnCalcPCA)
        self.button_calcpca.setEnabled(False)
        vbox21.addWidget(self.button_calcpca)
        self.button_savepca = QtWidgets.QPushButton('Save PCA Results...')
        self.button_savepca.clicked.connect( self.OnSave)
        self.button_savepca.setEnabled(False)
        vbox21.addWidget(self.button_savepca)


        hbox21 = QtWidgets.QHBoxLayout()
        text1 = QtWidgets.QLabel(self)
        text1.setText('Number of significant components')


        self.npcaspin = QtWidgets.QSpinBox()
        self.npcaspin.setRange(1,20)
        self.npcaspin.valueChanged[int].connect(self.OnNPCAspin)

        hbox21.addWidget(text1)
        hbox21.addWidget(self.npcaspin)
        vbox21.addLayout(hbox21)

        hbox22 = QtWidgets.QHBoxLayout()
        text2 = QtWidgets.QLabel(self)
        text2.setText( 'Cumulative variance')
        self.vartc = QtWidgets.QLabel(self)
        self.vartc.setText('0%')
        hbox22.addWidget(text2)
        hbox22.addWidget(self.vartc)

        vbox21.addLayout(hbox22)

        self.button_movepcup = QtWidgets.QPushButton('Move PC up')
        self.button_movepcup.clicked.connect( self.OnMovePCUP)
        self.button_movepcup.setEnabled(False)
        vbox21.addWidget(self.button_movepcup)


#         line = QtWidgets.QFrame()
#         line.setFrameShape(QtWidgets.QFrame.HLine)
#         line.setFrameShadow(QtWidgets.QFrame.Sunken)
#         vbox21.addWidget(line)

        #Todo: Repair 4d pca
        self.button_calcpca4D = QtWidgets.QPushButton('Calculate PCA for all angles')
        self.button_calcpca4D.clicked.connect( self.OnCalcPCA4D)
        self.button_calcpca4D.setEnabled(False)
        self.button_calcpca4D.setVisible(False)
        vbox21.addWidget(self.button_calcpca4D)


        sizer2.setLayout(vbox21)
        vbox2.addStretch(1)
        vbox2.addWidget(sizer2)
        vbox2.addStretch(3)


        #panel 3
        vbox3 = QtWidgets.QVBoxLayout()



        self.text_pcaspec = QtWidgets.QLabel(self)
        self.text_pcaspec.setText("PCA spectrum ")
        vbox3.addWidget(self.text_pcaspec)

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.pcaspecfig = Figure((PlotW, PlotH))
        self.PCASpecPan = FigureCanvas(self.pcaspecfig)
        self.PCASpecPan.setParent(self)

        fbox.addWidget(self.PCASpecPan)
        frame.setLayout(fbox)
        vbox3.addWidget(frame)


        #panel 4
        vbox4 = QtWidgets.QVBoxLayout()

        text4 = QtWidgets.QLabel(self)
        text4.setText("PCA eigenvalues ")
        vbox4.addWidget(text4)

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.pcaevalsfig = Figure((PlotW, PlotH*0.75))
        self.PCAEvalsPan = FigureCanvas(self.pcaevalsfig)
        self.PCAEvalsPan.setParent(self)
        self.PCAEvalsPan.mpl_connect('button_press_event', self.OnPointEvalsImage)

        fbox.addWidget(self.PCAEvalsPan)
        frame.setLayout(fbox)
        vbox4.addWidget(frame)



        vboxtop = QtWidgets.QVBoxLayout()
        gridsizertop = QtWidgets.QGridLayout()

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

        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
        self.calcpca = False
        self.selpca = 1
        self.numsigpca = 2
        self.slidershow.setValue(self.selpca)

        scrollmax = np.min([self.stk.n_ev, 20])
        self.slidershow.setMaximum(scrollmax)

        #try:
        self.CalcPCA()
        self.calcpca = True
        self.loadPCAImage()
        self.loadPCASpectrum()
        self.showEvals()
        self.com.pca_calculated = 1
        QtWidgets.QApplication.restoreOverrideCursor()
        #except:
        #    pass
        #    self.com.pca_calculated = 0
        #    QtWidgets.QApplication.restoreOverrideCursor()
        #    QtWidgets.QMessageBox.warning(self, 'Error', 'PCA not calculated.')

        self.window().refresh_widgets()


#----------------------------------------------------------------------
    def OnCalcPCA4D(self, event):

        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
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

            QtWidgets.QApplication.restoreOverrideCursor()
        except:
            self.com.pca_calculated = 0
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, 'Error', 'PCA not calculated.')




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

        savewin = SaveWinP2(self)
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

                    im = axes.imshow(np.rot90(self.pcaimage), cmap=matplotlib.colormaps["seismic_r"], vmin = -bound, vmax = bound)
                    cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)
                    axes.axis("off")

                    fileName_img = self.SaveFileName+"_PCA_" +str(i+1)+"."+ext
                    fig.savefig(fileName_img, pad_inches = 0.0)

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

                    im = axes.imshow(np.rot90(self.pcaimage), cmap=matplotlib.colormaps["seismic_r"], vmin = -bound, vmax = bound)
                    cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)
                    axes.axis("off")

                    fileName_img = self.SaveFileName+"_PCA_" +str(i+1)+"."+ext
                    fig.savefig(fileName_img, pad_inches = 0.0)

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

                    im = axes.imshow(np.rot90(self.pcaimage), cmap=matplotlib.colormaps["seismic_r"], vmin = -bound, vmax = bound)
                    cbar = axes.figure.colorbar(im, orientation='vertical',cax=ax_cb)
                    axes.axis("off")

                    fileName_img = self.SaveFileName+"_PCA_" +str(i+1)+"."+ext
                    fig.savefig(fileName_img, pad_inches = 0.0)

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

        except IOError as e:
            if e.strerror:
                err = e.strerror
            else:
                err = e
            print(err)
            QtWidgets.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err)


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


        im = axes.imshow(np.rot90(self.pcaimage), cmap=matplotlib.colormaps["seismic_r"], vmin = -bound, vmax = bound)
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
