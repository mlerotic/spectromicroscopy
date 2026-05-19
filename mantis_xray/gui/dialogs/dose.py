import numpy as np
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtCore import Qt
from ...analysis import henke

class DoseCalculation(QtWidgets.QDialog):

    def __init__(self, parent,  stack, ROIspectrum, roi_spectra=None):
        QtWidgets.QWidget.__init__(self, parent)

        self.parent = parent

        self.resize(300, 170)
        self.setWindowTitle('Dose Calculation')

        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)

        self.stack = stack
        self.ROIspectrum = ROIspectrum
        self.roi_spectra = roi_spectra or []



        vboxtop = QtWidgets.QVBoxLayout()


        gridtop = QtWidgets.QGridLayout()


        #fontb = wx.SystemSettings_GetFont(wx.SYS_DEFAULT_GUI_FONT)
        #fontb.SetWeight(wx.BOLD)


        st1 = QtWidgets.QLabel(self)
        st1.setText('Detector efficiency [%]:')
        #st1.SetFont(fontb)
        st2 = QtWidgets.QLabel(self)
        st2.setText('I region composition:')
        #st2.SetFont(fontb)
        st_roi = QtWidgets.QLabel(self)
        st_roi.setText('ROI selection:')
#        st3 = QtWidgets.QLabel(self,'Xray absorption length:')
#        st3.SetFont(fontb)
        st4 = QtWidgets.QLabel(self)
        st4.setText('Dose [Gray]:')
        #st4.SetFont(fontb)


        self.tc_1 = QtWidgets.QLineEdit(self)
        self.tc_1.setText('30')

        self.tc_2 = QtWidgets.QLineEdit(self)

        self.cb_roi = QtWidgets.QComboBox(self)
        if self.roi_spectra:
            self.cb_roi.addItems([label for label, _ in self.roi_spectra])
        else:
            self.cb_roi.addItem('Current ROI')

#        self.tc_3 = wx.TextCtrl(panel1, -1, size=((200,-1)), style=wx.TE_RICH|wx.VSCROLL|wx.TE_READONLY,
#                                         value=' ')

        self.tc_4 = QtWidgets.QLabel(self)


        gridtop.addWidget(st1, 0,0)
        gridtop.addWidget( self.tc_1, 0,1)
        gridtop.addWidget(st2, 1,0)
        gridtop.addWidget( self.tc_2, 1,1)
        gridtop.addWidget(st_roi, 2,0)
        gridtop.addWidget(self.cb_roi, 2,1)
#        gridtop.addWidget(st3, 0)
#        gridtop.addWidget( self.tc_3, 0)
        gridtop.addWidget(st4, 3,0)
        gridtop.addWidget( self.tc_4, 3,1)



        button_calcdose = QtWidgets.QPushButton('Calculate Dose')
        button_calcdose.clicked.connect(self.OnCalcDose)

        button_cancel = QtWidgets.QPushButton('Dismiss')
        button_cancel.clicked.connect(self.close)


        vboxtop.addLayout(gridtop)
        vboxtop.addWidget(button_calcdose)
        vboxtop.addWidget(button_cancel)

        self.setLayout(vboxtop)


#----------------------------------------------------------------------
    def CalcDose(self):


        try:
            detector_eff = 0.01*float(self.tc_1.text())
        except:
            QtWidgets.QMessageBox.warning(self, 'Error', 'Please enter numeric number for detector efficiency.')
            print('Please enter numeric number for detector efficiency.')
            return


        i_composition = str(self.tc_2.text())

        selected_spectrum = np.asarray(self.ROIspectrum, dtype=float)
        if self.roi_spectra:
            idx = self.cb_roi.currentIndex()
            if idx < 0 or idx >= len(self.roi_spectra):
                QtWidgets.QMessageBox.warning(self, 'Error', 'Please select a valid ROI.')
                return
            selected_spectrum = np.asarray(self.roi_spectra[idx][1], dtype=float)

        dose = 0.

        Chenke = henke.henke()

        #Check if composition array is recognizable
        try:
            z_array, atwt = Chenke.compound(i_composition,1.0)
        except:
            QtWidgets.QMessageBox.warning(self, 'Error', "Please enter new compound.")
            return

        try:
            dose = Chenke.dose_calc(self.stack, i_composition, selected_spectrum, self.stack.i0data, detector_eff)
        except Exception as e:
            QtWidgets.QMessageBox.warning(self, 'Error', "Could not calculate dose: {}".format(e))
            return

        self.tc_4.setText(str(dose))

        return


#----------------------------------------------------------------------
    def OnCalcDose(self, evt):

        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
        self.CalcDose()
        QtWidgets.QApplication.restoreOverrideCursor()


#----------------------------------------------------------------------
