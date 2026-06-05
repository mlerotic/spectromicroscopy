
from PyQt5 import QtWidgets, QtGui
import os, webbrowser
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
        text1.linkActivated.connect(self.openUrl)
        text1.setText("<A href='https://www.2ndlookconsulting.com'>www.2ndlookconsulting.com</A>")
        text1.setStyleSheet('color: rgb(53,159,217);font-size: 14pt; font-family: SansSerif;')
        text1.setTextFormat(1)

        #font2 = wx.Font(12, wx.SWISS, wx.NORMAL, wx.NORMAL)
        text2 = QtWidgets.QLabel(self)
        text2.setText('Mantis '+version)
        text2.setStyleSheet('color: rgb(0,0,0);font-size: 12pt')
        #text2.SetFont(font2)

        #font3 = wx.Font(8, wx.SWISS, wx.NORMAL, wx.NORMAL)
        text3 = QtWidgets.QLabel(self)
        text3.linkActivated.connect(self.openUrl)
        text3.setText( '''
Current development and maintenance by Benjamin Watts and <BR>
Jan-David Förster, with contributions from <A href='https://github.com/mlerotic/spectromicroscopy/graphs/contributors?all=1'>others</A>.<BR>
Originally developed by Mirna Lerotic, based on earlier programs by<BR>
Mirna Lerotic and Chris Jacobsen. Initial development supported by<BR>
Argonne National Laboratory LDRD 2010-193-R1 9113. ''')
        text3.setStyleSheet('color: rgb(0,0,0)')
        text3.setTextFormat(1)


        #font4 = wx.Font(8, wx.SWISS, wx.NORMAL, wx.NORMAL)
        text4 = QtWidgets.QLabel(self)
        text4.linkActivated.connect(self.openUrl)
        text4.setText( '''
Mantis is free software: you can redistribute it and/or modify<BR>
it under the terms of the GNU General Public License as published<BR>
by the Free Software Foundation, either version 3 of the License,<BR>
or any later version.<BR>
<BR>
Mantis is distributed in the hope that it will be useful,<BR>
but WITHOUT ANY WARRANTY; without even the implied warranty<BR>
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.<BR>
See the <A href='https://www.gnu.org/licenses/'>GNU General Public License</A> for more details.<BR>''')
        text4.setStyleSheet('color: rgb(0,0,0)')
        text4.setTextFormat(1)

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

    def openUrl(self,url):
        print(url)
        webbrowser.open(url, new=1, autoraise=True)
