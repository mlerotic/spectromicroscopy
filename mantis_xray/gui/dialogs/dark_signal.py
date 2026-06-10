from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt

class DarkSignal(QtWidgets.QDialog):

    def __init__(self, parent, common, stack):
        QtWidgets.QWidget.__init__(self, parent)

        self.parent = parent

        self.stack = stack
        self.com = common

        self.resize(300, 170)
        self.setWindowTitle('Dark Signal Correction')

        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)


        vboxtop = QtWidgets.QVBoxLayout()


        gridtop = QtWidgets.QGridLayout()


        #fontb = wx.SystemSettings_GetFont(wx.SYS_DEFAULT_GUI_FONT)
        #fontb.SetWeight(wx.BOLD)


        st1 = QtWidgets.QLabel(self)
        st1.setText('Dark Signal Value:')


        self.ntc_ds = QtWidgets.QLineEdit(self)
        self.ntc_ds.setFixedWidth(150)
        self.ntc_ds.setValidator(QtGui.QDoubleValidator(-99999, 99999, 2, self))
        self.ntc_ds.setAlignment(QtCore.Qt.AlignRight)

        self.ntc_ds.setText(str(0.0))


        gridtop.addWidget(st1, 0,0)
        gridtop.addWidget(self.ntc_ds, 0,1)


        button_dscalc = QtWidgets.QPushButton('Subtract Dark Signal')
        button_dscalc.clicked.connect(self.OnDSCalc)

        button_cancel = QtWidgets.QPushButton('Dismiss')
        button_cancel.clicked.connect(self.close)


        vboxtop.addLayout(gridtop)
        vboxtop.addWidget(button_dscalc)
        vboxtop.addWidget(button_cancel)

        self.setLayout(vboxtop)


#----------------------------------------------------------------------
    def OnDSCalc(self):

        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))

        darksig = 0.0

        try:
            value = self.ntc_ds.text()
            darksig = float(value)
        except:
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, 'Error', 'Please enter numeric number for dark signal.')
            return


        self.stack.absdata = self.stack.absdata - darksig

        if self.com.i0_loaded == 1:
            self.stack.calculate_optical_density()

        self.stack.fill_h5_struct_from_stk()
        if self.com.i0_loaded == 1:
            self.stack.fill_h5_struct_normalization()



        self.parent.tab_prep.absimgfig.loadNewImageWithROI()
        self.parent.tab_load.absimgfig.loadNewImage()
        self.parent.tab_prep.specfig.ClearandReload()
        self.parent.window().refresh_widgets()



        QtWidgets.QApplication.restoreOverrideCursor()

        self.close()

        return
