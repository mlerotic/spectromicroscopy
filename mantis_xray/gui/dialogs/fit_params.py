
from PyQt5 import QtWidgets, QtGui

class FitParams(QtWidgets.QDialog):

    def __init__(self, parent, title, fitparams, readonly, initialization):
        QtWidgets.QWidget.__init__(self, parent)

        self.parent = parent

        self.base = fitparams[0]
        self.stepfitparams = fitparams[1]
        self.gauss_fp_a = fitparams[2]
        self.gauss_fp_m = fitparams[3]
        self.gauss_fp_s = fitparams[4]

        self.init = initialization


        self.resize(600, 450)
        self.setWindowTitle(title)

        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)


        vboxtop = QtWidgets.QVBoxLayout()

        #panel
        sizer = QtWidgets.QGroupBox('Fit Parameters')
        vbox = QtWidgets.QVBoxLayout()


        fgs1 = QtWidgets.QGridLayout()



        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Amplitude')
        fgs1.addWidget(textctrl, 0, 1)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Position')
        fgs1.addWidget(textctrl, 0, 2)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('FWHM')
        fgs1.addWidget(textctrl, 0, 3)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Base')
        fgs1.addWidget(textctrl, 0, 4)


        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Step 1')
        fgs1.addWidget(textctrl, 1, 0)

        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Step 2')
        fgs1.addWidget(textctrl, 2, 0)

        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Amplitude')
        fgs1.addWidget(textctrl, 3, 1)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Mu')
        fgs1.addWidget(textctrl, 3, 2)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Sigma')
        fgs1.addWidget(textctrl, 3, 3)


        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 1')
        fgs1.addWidget(textctrl, 4, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 2')
        fgs1.addWidget(textctrl, 5, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 3')
        fgs1.addWidget(textctrl, 6, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 4')
        fgs1.addWidget(textctrl, 7, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 5')
        fgs1.addWidget(textctrl, 8, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 6')
        fgs1.addWidget(textctrl, 9, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 7')
        fgs1.addWidget(textctrl, 10, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 8')
        fgs1.addWidget(textctrl, 11, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 9')
        fgs1.addWidget(textctrl, 12, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 10')
        fgs1.addWidget(textctrl, 13, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 11')
        fgs1.addWidget(textctrl, 14, 0)
        textctrl =  QtWidgets.QLabel(self)
        textctrl.setText('Gauss 12')
        fgs1.addWidget(textctrl, 15, 0)


        self.le_sa1 = QtWidgets.QLineEdit(self)
        self.le_sa1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sa1.setReadOnly(readonly)
        self.le_sa1.setText(str(self.stepfitparams[0]))
        fgs1.addWidget(self.le_sa1, 1, 1)

        self.le_sp1 = QtWidgets.QLineEdit(self)
        self.le_sp1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sp1.setReadOnly(readonly)
        self.le_sp1.setText(str(self.stepfitparams[1]))
        fgs1.addWidget(self.le_sp1, 1, 2)

        self.le_sf1 = QtWidgets.QLineEdit(self)
        self.le_sf1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sf1.setReadOnly(readonly)
        self.le_sf1.setText(str(self.stepfitparams[2]))
        fgs1.addWidget(self.le_sf1, 1, 3)

        self.le_base = QtWidgets.QLineEdit(self)
        self.le_base.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_base.setReadOnly(readonly)
        self.le_base.setText(str(self.base))
        fgs1.addWidget(self.le_base, 1, 4)


        self.le_sa2 = QtWidgets.QLineEdit(self)
        self.le_sa2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sa2.setReadOnly(readonly)
        self.le_sa2.setText(str(self.stepfitparams[3]))
        fgs1.addWidget(self.le_sa2, 2, 1)

        self.le_sp2 = QtWidgets.QLineEdit(self)
        self.le_sp2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sp2.setReadOnly(readonly)
        self.le_sp2.setText(str(self.stepfitparams[4]))
        fgs1.addWidget(self.le_sp2, 2, 2)

        self.le_sf2 = QtWidgets.QLineEdit(self)
        self.le_sf2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        self.le_sf2.setReadOnly(readonly)
        self.le_sf2.setText(str(self.stepfitparams[5]))
        fgs1.addWidget(self.le_sf2, 2, 3)



        self.le_amplitudes = []
        self.le_mu = []
        self.le_sigma = []

        self.le_ga1 = QtWidgets.QLineEdit(self)
        self.le_ga1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga1, 4, 1)
        self.le_amplitudes.append(self.le_ga1)

        self.le_gm1 = QtWidgets.QLineEdit(self)
        self.le_gm1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm1, 4, 2)
        self.le_mu.append(self.le_gm1)

        self.le_gs1 = QtWidgets.QLineEdit(self)
        self.le_gs1.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs1, 4, 3)
        self.le_sigma.append(self.le_gs1)

        self.le_ga2 = QtWidgets.QLineEdit(self)
        self.le_ga2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga2, 5, 1)
        self.le_amplitudes.append(self.le_ga2)

        self.le_gm2 = QtWidgets.QLineEdit(self)
        self.le_gm2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm2, 5, 2)
        self.le_mu.append(self.le_gm2)

        self.le_gs2 = QtWidgets.QLineEdit(self)
        self.le_gs2.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs2, 5, 3)
        self.le_sigma.append(self.le_gs2)

        self.le_ga3 = QtWidgets.QLineEdit(self)
        self.le_ga3.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga3, 6, 1)
        self.le_amplitudes.append(self.le_ga3)

        self.le_gm3 = QtWidgets.QLineEdit(self)
        self.le_gm3.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm3, 6, 2)
        self.le_mu.append(self.le_gm3)

        self.le_gs3 = QtWidgets.QLineEdit(self)
        self.le_gs3.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs3, 6, 3)
        self.le_sigma.append(self.le_gs3)

        self.le_ga4 = QtWidgets.QLineEdit(self)
        self.le_ga4.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga4, 7, 1)
        self.le_amplitudes.append(self.le_ga4)

        self.le_gm4 = QtWidgets.QLineEdit(self)
        self.le_gm4.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm4, 7, 2)
        self.le_mu.append(self.le_gm4)

        self.le_gs4 = QtWidgets.QLineEdit(self)
        self.le_gs4.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs4, 7, 3)
        self.le_sigma.append(self.le_gs4)

        self.le_ga5 = QtWidgets.QLineEdit(self)
        self.le_ga5.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga5, 8, 1)
        self.le_amplitudes.append(self.le_ga5)

        self.le_gm5 = QtWidgets.QLineEdit(self)
        self.le_gm5.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm5, 8, 2)
        self.le_mu.append(self.le_gm5)

        self.le_gs5 = QtWidgets.QLineEdit(self)
        self.le_gs5.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs5, 8, 3)
        self.le_sigma.append(self.le_gs5)

        self.le_ga6 = QtWidgets.QLineEdit(self)
        self.le_ga6.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga6, 9, 1)
        self.le_amplitudes.append(self.le_ga6)

        self.le_gm6 = QtWidgets.QLineEdit(self)
        self.le_gm6.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm6, 9, 2)
        self.le_mu.append(self.le_gm6)

        self.le_gs6 = QtWidgets.QLineEdit(self)
        self.le_gs6.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs6, 9, 3)
        self.le_sigma.append(self.le_gs6)

        self.le_ga7 = QtWidgets.QLineEdit(self)
        self.le_ga7.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga7, 10, 1)
        self.le_amplitudes.append(self.le_ga7)

        self.le_gm7 = QtWidgets.QLineEdit(self)
        self.le_gm7.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm7, 10, 2)
        self.le_mu.append(self.le_gm7)

        self.le_gs7 = QtWidgets.QLineEdit(self)
        self.le_gs7.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs7, 10, 3)
        self.le_sigma.append(self.le_gs7)

        self.le_ga8 = QtWidgets.QLineEdit(self)
        self.le_ga8.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga8, 11, 1)
        self.le_amplitudes.append(self.le_ga8)

        self.le_gm8 = QtWidgets.QLineEdit(self)
        self.le_gm8.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm8, 11, 2)
        self.le_mu.append(self.le_gm8)

        self.le_gs8 = QtWidgets.QLineEdit(self)
        self.le_gs8.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs8, 11, 3)
        self.le_sigma.append(self.le_gs8)

        self.le_ga9 = QtWidgets.QLineEdit(self)
        self.le_ga9.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga9, 12, 1)
        self.le_amplitudes.append(self.le_ga9)

        self.le_gm9 = QtWidgets.QLineEdit(self)
        self.le_gm9.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm9, 12, 2)
        self.le_mu.append(self.le_gm9)

        self.le_gs9 = QtWidgets.QLineEdit(self)
        self.le_gs9.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs9, 12, 3)
        self.le_sigma.append(self.le_gs9)

        self.le_ga10 = QtWidgets.QLineEdit(self)
        self.le_ga10.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga10, 13, 1)
        self.le_amplitudes.append(self.le_ga10)

        self.le_gm10 = QtWidgets.QLineEdit(self)
        self.le_gm10.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm10, 13, 2)
        self.le_mu.append(self.le_gm10)

        self.le_gs10 = QtWidgets.QLineEdit(self)
        self.le_gs10.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs10, 13, 3)
        self.le_sigma.append(self.le_gs10)

        self.le_ga11 = QtWidgets.QLineEdit(self)
        self.le_ga11.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga11, 14, 1)
        self.le_amplitudes.append(self.le_ga11)

        self.le_gm11 = QtWidgets.QLineEdit(self)
        self.le_gm11.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm11, 14, 2)
        self.le_mu.append(self.le_gm11)

        self.le_gs11 = QtWidgets.QLineEdit(self)
        self.le_gs11.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs11, 14, 3)
        self.le_sigma.append(self.le_gs11)

        self.le_ga12 = QtWidgets.QLineEdit(self)
        self.le_ga12.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_ga12, 15, 1)
        self.le_amplitudes.append(self.le_ga12)

        self.le_gm12 = QtWidgets.QLineEdit(self)
        self.le_gm12.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gm12, 15, 2)
        self.le_mu.append(self.le_gm12)

        self.le_gs12 = QtWidgets.QLineEdit(self)
        self.le_gs12.setValidator(QtGui.QDoubleValidator(0, 99999, 2, self))
        fgs1.addWidget(self.le_gs12, 15, 3)
        self.le_sigma.append(self.le_gs12)

        for i in range(12):
            self.le_amplitudes[i].setReadOnly(readonly)
            self.le_mu[i].setReadOnly(readonly)
            self.le_sigma[i].setReadOnly(readonly)

            self.le_amplitudes[i].setText(str(self.gauss_fp_a[i]))
            self.le_mu[i].setText(str(self.gauss_fp_m[i]))
            self.le_sigma[i].setText(str(self.gauss_fp_s[i]))

        vbox.addLayout(fgs1)
        sizer.setLayout(vbox)

        vboxtop.addWidget(sizer)

        if initialization:
            bname = 'Confirm'
        else:
            bname = 'Close'
        self.button_close = QtWidgets.QPushButton(bname)
        self.button_close.clicked.connect( self.OnClose)
        vboxtop.addWidget(self.button_close)


        self.setLayout(vboxtop)


#----------------------------------------------------------------------
    def OnClose(self, event):


        if self.init:

            self.base = float(self.le_base.text())


            self.stepfitparams[0] = float(self.le_sa1.text())
            self.stepfitparams[1] = float(self.le_sp1.text())
            self.stepfitparams[2] = float(self.le_sf1.text())
            self.stepfitparams[3] = float(self.le_sa2.text())
            self.stepfitparams[4] = float(self.le_sp2.text())
            self.stepfitparams[5] = float(self.le_sf2.text())

            for i in range(12):

                self.gauss_fp_a[i] = float(self.le_amplitudes[i].text())
                self.gauss_fp_m[i] = float(self.le_mu[i].text())
                self.gauss_fp_s[i] = float(self.le_sigma[i].text())


            self.parent.page6.SetInitFitParams(self.base, self.stepfitparams, self.gauss_fp_a, self.gauss_fp_m, self.gauss_fp_s)

        self.close()
