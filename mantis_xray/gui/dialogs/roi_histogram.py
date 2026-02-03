
from PyQt5 import QtWidgets, QtCore, QtGui
import numpy as np
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class ROIHistogram(QtWidgets.QDialog):

    def __init__(self, parent, image, n_cols, n_rows):
        QtWidgets.QWidget.__init__(self, parent)

        self.parent = parent


        self.resize(600, 600)
        self.setWindowTitle('Histogram')

        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)

        self.image = image

        self.n_cols = n_cols
        self.n_rows = n_rows


        imgmax = np.max(self.image)
        self.histmin = 0.3*imgmax
        self.histmax = 0.6*imgmax


        vbox = QtWidgets.QVBoxLayout()

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.histfig = Figure((6.0, 4.2))
        self.HistogramPanel = FigureCanvas(self.histfig)
        self.HistogramPanel.setParent(self)
        self.HistogramPanel.mpl_connect('button_press_event', self.OnSelection1)
        self.HistogramPanel.mpl_connect('button_release_event', self.OnSelection2)


        fbox.addWidget(self.HistogramPanel)
        frame.setLayout(fbox)
        vbox.addWidget(frame)

        hbox = QtWidgets.QHBoxLayout()
        vbox1 = QtWidgets.QVBoxLayout()
        sizer1 = QtWidgets.QGroupBox('ROI selection')


        tc1 = QtWidgets.QLabel(self)
        tc1.setText('Select an ROI region by clicking and dragging a mouse on the histogram.')
        vbox1.addWidget(tc1)

        self.textctrl1 = QtWidgets.QLabel(self)
        self.textctrl1.setText('Min =  {0:04.3f}'.format(float(self.histmin)))
        vbox1.addWidget(self.textctrl1)

        self.textctrl2 = QtWidgets.QLabel(self)
        self.textctrl2.setText('Max =  {0:04.3f}'.format( float(self.histmax)))
        vbox1.addWidget(self.textctrl2)

        sizer1.setLayout(vbox1)
        hbox.addWidget(sizer1)

        self.absimgfig = Figure((4,4))
        self.AbsImagePanel = FigureCanvas(self.absimgfig)
        self.AbsImagePanel.setParent(self)

        hbox.addWidget(self.AbsImagePanel,0,QtCore .Qt. AlignLeft)

        vbox.addLayout(hbox)

        hbox2 = QtWidgets.QHBoxLayout()
        button_ok = QtWidgets.QPushButton('Accept')
        button_ok.clicked.connect(self.OnAccept)
        hbox2.addWidget(button_ok)

        button_cancel = QtWidgets.QPushButton('Cancel')
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

        histogram = self.image.flatten()

        self.n, self.bins, patches = self.axes.hist(histogram, 200, normed=1, facecolor='green', alpha=0.75)

        self.axes.set_xlabel('Optical Density')
        self.axes.set_ylabel('Percentage of Pixels')

        self.patch = self.axes.axvspan(self.histmin, self.histmax, facecolor='r', alpha=0.3)


        self.HistogramPanel.draw()




#----------------------------------------------------------------------
    def draw_image(self):


        image = self.image.copy()

        fluxmin = self.histmin
        fluxmax = self.histmax

        hist_indices = np.where((fluxmin<image)&(image<fluxmax))

        redpix = np.zeros((self.n_cols,self.n_rows))
        redpix[hist_indices] = 255

        redpix = np.ma.array(redpix)

        redpix_masked =  np.ma.masked_values(redpix, 0)


        fig = self.absimgfig
        fig.clf()
        fig.add_axes(((0.0,0.0,1.0,1.0)))


        axes = fig.gca()
        fig.patch.set_alpha(1.0)

        im = axes.imshow(np.rot90(image), cmap=matplotlib.colormaps["gray"])

        im_red = axes.imshow(np.rot90(redpix_masked),cmap=matplotlib.colormaps["autumn"])

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

        self.textctrl1.setText('Min =  {0:04.3f}'.format(float(self.histmin)))
        self.textctrl2.setText('Max =  {0:04.3f}'.format( float(self.histmax)))

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

        self.parent.MakeHistogramROI(self.histmin, self.histmax)

        self.close()
