
from PyQt5 import QtWidgets, QtGui
import os
from ...helpers import resource_path
from ... import __version__ as version

class AboutFrame(QtWidgets.QDialog):

    def __init__(self, parent = None, title='About'):
        QtWidgets.QWidget.__init__(self, parent)

        self.resize(360, 620)
        self.setFixedSize(360, 620)
        self.setWindowTitle('About Mantis')

        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)

        vbox = QtWidgets.QVBoxLayout()

        self.image = QtGui.QImage(resource_path(os.path.join('images','Mantis_logo_about.png')))

        self.imageLabel = QtWidgets.QLabel()
        self.imageLabel.setBackgroundRole(QtGui.QPalette.Base)

        self.imageLabel.setPixmap(QtGui.QPixmap.fromImage(self.image))
        vbox.addWidget(self.imageLabel)



        text1 = QtWidgets.QLabel(self)
        text1.setText("www.2ndlookconsulting.com")
        text1.setStyleSheet('color: rgb(53,159,217);font-size: 14pt; font-family: SansSerif;')
        #text1.setFont(QtGui.QFont('SansSerif', 14))

        #font2 = wx.Font(12, wx.SWISS, wx.NORMAL, wx.NORMAL)
        text2 = QtWidgets.QLabel(self)
        text2.setText('Mantis '+version)
        text2.setStyleSheet('color: rgb(0,0,0);font-size: 12pt')
        #text2.SetFont(font2)

        #font3 = wx.Font(8, wx.SWISS, wx.NORMAL, wx.NORMAL)
        text3 = QtWidgets.QLabel(self)
        text3.setText( '''
Developed by Mirna Lerotic, based on earlier programs by Mirna
Lerotic and Chris Jacobsen. Initial development supported by
Argonne National Laboratory LDRD 2010-193-R1 9113. ''')
        text3.setStyleSheet('color: rgb(0,0,0)')
        #text3.SetFont(font3)


        #font4 = wx.Font(8, wx.SWISS, wx.NORMAL, wx.NORMAL)
        text4 = QtWidgets.QLabel(self)
        text4.setText( '''
Mantis is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License,
or any later version.

Mantis is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details
http://www.gnu.org/licenses/.''')
        text4.setStyleSheet('color: rgb(0,0,0)')
        #text4.SetFont(font4)

        vbox.addStretch(1)
        hbox = QtWidgets.QHBoxLayout()
        hbox.addStretch(1)
        vbox2 = QtWidgets.QVBoxLayout()
        vbox2.addWidget(text1)
        vbox2.addStretch(1)
        vbox2.addWidget(text2)
        vbox2.addWidget(text3)
        vbox2.addWidget(text4)
        hbox.addLayout(vbox2)
        hbox.addStretch(1)
        vbox.addLayout(hbox)
        vbox.addStretch(1)

        vbox.addStretch(1)

        esc_shortcut = QtWidgets.QShortcut(QtGui.QKeySequence("Esc"), self)
        esc_shortcut.activated.connect(self.reject)



        self.setLayout(vbox)
