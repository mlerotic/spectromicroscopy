
import numpy as np
from PyQt5 import QtWidgets, QtGui, QtCore
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

class Scatterplots(QtWidgets.QDialog):

    def __init__(self, parent,  common, analz):
        QtWidgets.QWidget.__init__(self, parent)

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



        vbox = QtWidgets.QVBoxLayout()

        grid1 = QtWidgets.QGridLayout()

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()
        self.scattplfig = Figure((5.0, 4.8))
        self.ScatterPPanel = FigureCanvas(self.scattplfig)
        self.ScatterPPanel.setParent(self)
        fbox.addWidget(self.ScatterPPanel)
        frame.setLayout(fbox)

        self.slidershow_y = QtWidgets.QSlider(QtCore.Qt.Vertical)
        self.slidershow_y.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slidershow_y.setRange(1, self.numsigpca)
        self.slidershow_y.setValue(self.pca_y)
        self.slidershow_y.valueChanged[int].connect(self.OnSliderScroll_y)

        grid1.addWidget(self.slidershow_y, 0, 0)
        grid1.addWidget(frame, 0, 1)


        self.slidershow_x = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.slidershow_x.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slidershow_x.setRange(1, self.numsigpca)
        self.slidershow_x.setValue(self.pca_x)
        self.slidershow_x.valueChanged[int].connect(self.OnSliderScroll_x)

        #grid1.addWidget(wx.StaticText(panel, -1, ''))
        grid1.addWidget(self.slidershow_x, 1,  1)

        hbox = QtWidgets.QVBoxLayout()

        button_close = QtWidgets.QPushButton('Close')
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
