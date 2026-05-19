
import numpy as np
from PyQt5 import QtWidgets, QtGui, QtCore
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

class ShowMapHistogram(QtWidgets.QDialog):

    def __init__(self, parent, common, analz):
        QtWidgets.QWidget.__init__(self, parent)

        self.parent = parent

        self.limit = 1

        self.histmax = None


        self.resize(600, 500)
        self.setWindowTitle('Histogram')

        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)

        self.com = common
        self.anlz = analz

        vbox = QtWidgets.QVBoxLayout()

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.histfig = Figure((6.0, 4.2))
        self.HistogramPanel = FigureCanvas(self.histfig)
        self.HistogramPanel.setParent(self)
        self.HistogramPanel.mpl_connect('button_press_event', self.OnClick)


        fbox.addWidget(self.HistogramPanel)
        frame.setLayout(fbox)
        vbox.addWidget(frame)



        vbox1 = QtWidgets.QVBoxLayout()
        sizer1 = QtWidgets.QGroupBox('Histogram Cutoff')

        st = QtWidgets.QLabel(self)
        st.setText('Select a cutoff values on the histogram. All the values outside the defined limits will be set to cutoff limit value.')

        vbox1.addWidget(st)

        hbox1 = QtWidgets.QHBoxLayout()
        self.rb_min = QtWidgets.QRadioButton( 'Lower Limit', self)
        self.rb_max = QtWidgets.QRadioButton('Upper Limit',self)
        self.rb_min.setChecked(True)
        self.rb_min.toggled.connect(self.OnRb_limit)


        hbox1.addWidget(self.rb_min)
        hbox1.addWidget(self.rb_max)
        hbox1.addStretch (1)
        vbox1.addLayout(hbox1)

        self.tl_cutmin = QtWidgets.QLabel(self)
        self.tl_cutmin.setText('Lower Cutoff Value: ')
        vbox1.addWidget(self.tl_cutmin)

        self.tl_cutmax = QtWidgets.QLabel(self)
        self.tl_cutmax.setText('Upper Cutoff Value: ')
        vbox1.addWidget(self.tl_cutmax)

        sizer1.setLayout(vbox1)
        vbox.addWidget(sizer1)

        hbox2 = QtWidgets.QHBoxLayout()
        self.button_ok = QtWidgets.QPushButton('Accept')
        self.button_ok.clicked.connect(self.OnAccept)
        self.button_ok.setEnabled(False)
        hbox2.addWidget(self.button_ok)

        button_cancel = QtWidgets.QPushButton('Cancel')
        button_cancel.clicked.connect(self.close)
        hbox2.addWidget(button_cancel)

        vbox.addLayout(hbox2)

        self.setLayout(vbox)

        if len(self.anlz.original_svd_maps4D) == 0:
            self.draw_histogram()
        else:
            self.draw_histogram4D()




#----------------------------------------------------------------------
    def draw_histogram(self):


        fig = self.histfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        self.axes = fig.gca()

        #target_svd_maps;target_pcafit_maps
        if self.parent.tab_spec.showraw == True:
            self.histogram = self.anlz.original_svd_maps
        else:
            self.histogram = self.anlz.original_fit_maps


        histdata = np.reshape(self.histogram, (self.anlz.stack.n_cols*self.anlz.stack.n_rows*self.anlz.n_target_spectra), order='F')

        self.n, self.bins, patches = self.axes.hist(histdata, 200, normed=1, facecolor='green', alpha=0.75)

        self.axes.set_xlabel('Thickness per Pixel in Spectral Maps')
        self.axes.set_ylabel('Percentage of Pixels')

        self.HistogramPanel.draw()


#----------------------------------------------------------------------
    def draw_histogram4D(self):


        fig = self.histfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        self.axes = fig.gca()

        #target_svd_maps;target_pcafit_maps
        if self.parent.tab_spec.showraw == True:
            self.histogram = self.anlz.original_svd_maps4D
        else:
            self.histogram = self.anlz.original_fit_maps4D


        histdata = np.reshape(self.histogram, (self.anlz.stack.n_cols*self.anlz.stack.n_rows*self.anlz.n_target_spectra*self.anlz.stack.n_theta), order='F')

        self.n, self.bins, patches = self.axes.hist(histdata, 200, normed=1, facecolor='green', alpha=0.75)

        self.axes.set_xlabel('Thickness per Pixel in Spectral Maps')
        self.axes.set_ylabel('Percentage of Pixels')

        self.HistogramPanel.draw()

#----------------------------------------------------------------------
    def OnClick(self, evt):

        x1 = evt.xdata

        if x1 == None:
            return

        if self.limit == 1:
            self.tl_cutmin.setText('Lower Cutoff Value: '+str(x1))

            self.histmin = x1

        else:
            self.tl_cutmax.setText('Upper Cutoff Value: '+str(x1))

            self.histmax = x1

        self.button_ok.setEnabled(True)


#----------------------------------------------------------------------
    def OnRb_limit(self, enabled):

        state = enabled

        if state:
            self.limit = 1
        else:
            self.limit = 2


#----------------------------------------------------------------------
    def OnAccept(self, evt):

        if self.parent.tab_spec.showraw == True:
            self.anlz.svd_map_threshold(self.histmin, self.histmax, svd=True)
        else:
            self.anlz.svd_map_threshold(self.histmin, self.histmax, pca=True)
        self.parent.tab_spec.loadTargetMap()
        self.close()
