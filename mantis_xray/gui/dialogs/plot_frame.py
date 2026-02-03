from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class PlotFrame(QtWidgets.QDialog):
    def __init__(self, parent, datax, datay, title = "I0 data"):
        QtWidgets.QWidget.__init__(self, parent)

        self.parent = parent

        self.title = title

        self.datax = datax
        self.datay = datay

        self.resize(630, 500)
        self.setWindowTitle(title)

        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)



        vbox = QtWidgets.QVBoxLayout()

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.plotfig = Figure((6.0, 4.2))
        self.PlotPanel = FigureCanvas(self.plotfig)
        self.PlotPanel.setParent(self)

        fbox.addWidget(self.PlotPanel)
        frame.setLayout(fbox)
        vbox.addWidget(frame)


        hbox = QtWidgets.QHBoxLayout()

        button_save = QtWidgets.QPushButton('Save Spectrum')
        button_save.clicked.connect(self.OnSave)
        hbox.addWidget(button_save)

        button_close = QtWidgets.QPushButton('Close')
        button_close.clicked.connect(self.close)
        hbox.addWidget(button_close)

        vbox.addLayout(hbox)

        self.setLayout(vbox)


        self.draw_plot(datax,datay)


#----------------------------------------------------------------------
    def draw_plot(self, datax, datay):


        fig = self.plotfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        self.axes = fig.gca()


        plot = self.axes.plot(datax,datay)

        self.axes.set_xlabel('Photon Energy [eV]')
        self.axes.set_ylabel('I0 Flux')



        self.PlotPanel.draw()

#----------------------------------------------------------------------
    def OnSave(self, event):


        try:
            wildcard = "CSV files (*.csv)"

            filepath, _filter = QtWidgets.QFileDialog.getSaveFileName(self, 'Save plot as .csv (.csv)', '', wildcard)

            filepath = str(filepath)
            if filepath == '':
                return


            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
            self.Save(filepath)
            QtWidgets.QApplication.restoreOverrideCursor()

        except:

            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, 'Error', "Could not save .csv file.")

        self.close()


        return

#----------------------------------------------------------------------
    def Save(self, filename):

        f = open(filename, 'w')
        print('*********************  X-ray Absorption Data  ********************', file=f)
        print('*', file=f)
        print('* Formula: ', file=f)
        print('* Common name: ', self.title, file=f)
        print('* Edge: ', file=f)
        print('* Acquisition mode: ', file=f)
        print('* Source and purity: ', file=f)
        print('* Comments: ', file=f)
        print('* Delta eV: ', file=f)
        print('* Min eV: ', file=f)
        print('* Max eV: ', file=f)
        print('* Y axis: ', file=f)
        print('* Contact person: ', file=f)
        print('* Write date: ', file=f)
        print('* Journal: ', file=f)
        print('* Authors: ', file=f)
        print('* Title: ', file=f)
        print('* Volume: ', file=f)
        print('* Issue number: ', file=f)
        print('* Year: ', file=f)
        print('* Pages: ', file=f)
        print('* Booktitle: ', file=f)
        print('* Editors: ', file=f)
        print('* Publisher: ', file=f)
        print('* Address: ', file=f)
        print('*--------------------------------------------------------------', file=f)
        dim = self.datax.shape
        n=dim[0]
        for ie in range(n):
            print('{0:06.2f}, {1:06f}'.format(self.datax[ie], self.datay[ie]), file=f)

        f.close()
