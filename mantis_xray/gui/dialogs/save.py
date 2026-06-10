import os
import numpy as np
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtCore import Qt
from PyQt5 import uic
from PIL import Image
from skimage import exposure
from skimage.util import img_as_ubyte
import matplotlib.figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

class SaveWin(QtWidgets.QDialog):

    def __init__(self, parent, com, stk, i0_mode=False):
        QtWidgets.QWidget.__init__(self, parent)
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dialogsave.ui'), self)
        self.parent = parent
        self.setWindowTitle('Save')
        self.com = com
        self.stk = stk

        self.cb11.setChecked(True)
        self.cb21.setChecked(True)

        if i0_mode:
            # Batch export not applicable when saving I0 spectrum
            self.label_6.setVisible(False)
            self.cb32.setVisible(False)
            self.cb32.setChecked(False)
            self.cb35.setVisible(False)
            self.cb35.setChecked(False)

        if not hasattr(parent, 'specfig'):
            self.label_spectrum.setVisible(False)
            self.label_csv.setVisible(False)
            self.cb11.setVisible(False)
            self.cb11.setChecked(False)
            self.cb12.setVisible(False)
            self.cb13.setVisible(False)
            self.cb14.setVisible(False)

        path, ext = os.path.splitext(self.com.filename) # currently empty?
        ext = ext[1:].lower()
        suffix = "." + ext
        path, fn = os.path.split(self.com.filename)
        self.filename = fn[:-len(suffix)]
        self.tc_savefn.setText(self.filename)

        self.path = self.com.path
        self.tc_savepath.setText(self.path)
        self.button_path.clicked.connect(self.OnBrowseDir)
        self.buttonBox.accepted.connect(self.OnSave)
        self.buttonBox.rejected.connect(self.close)


#----------------------------------------------------------------------
    def OnBrowseDir(self, evt):

        directory = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose a directory", self.path, QtWidgets.QFileDialog.ShowDirsOnly|QtWidgets.QFileDialog.ReadOnly)



        if directory == '':
            return

        directory = str(directory)
        self.com.path = directory

        self.path = directory

        self.tc_savepath.setText(self.path)



#----------------------------------------------------------------------
    def OnSave(self):

        self.filename = str(self.tc_savefn.text())

        sp_pdf = self.cb11.isChecked()
        sp_png = self.cb12.isChecked()
        sp_svg = self.cb13.isChecked()
        sp_csv = self.cb14.isChecked()
        im_pdf = self.cb21.isChecked()
        im_png = self.cb22.isChecked()
        im_svg = self.cb23.isChecked()
        im_tif = self.cb25.isChecked()
        im_all = self.cb32.isChecked()
        im_all_tif = self.cb35.isChecked()

        self.close()
        self.Save(self.filename, self.path,
                                         spec_png = sp_png,
                                         spec_pdf = sp_pdf,
                                         spec_svg = sp_svg,
                                         sp_csv = sp_csv,
                                         img_png = im_png,
                                         img_pdf = im_pdf,
                                         img_svg = im_svg,
                                         img_tif = im_tif,
                                         img_all = im_all,
                                         img_all_tif = im_all_tif)


#----------------------------------------------------------------------
    def Save(self, filename, path, spec_png = True, spec_pdf = False, spec_svg = False, sp_csv = False,
             img_png = True, img_pdf = False, img_svg = False, img_tif = False, img_all = False, img_all_tif = False):

        self.SaveFileName = os.path.join(path,filename)

        try:
            ext = 'png'
            suffix = "." + ext

            if spec_png:
                fileName_spec = self.SaveFileName+"_spectrum."+ext
                fig = self.parent.specfig
                fig.SaveFig(fileName_spec)

            if img_png:

                fileName_img = self.SaveFileName+"_" +str(self.stk.ev[self.parent.iev])+"eV."+ext
                fig = self.parent.absimgfig
                fig.SaveFig(fileName_img)

            #Save all images in the stack
            if img_all:
                QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))

                for i in range (self.stk.n_ev):
                    if self.parent.showflux:
                        #Show flux image
                        image = self.stk.absdata[:,:,i]
                    else:
                        #Show OD image
                        image = self.stk.od3d[:,:,i]
                    image = img_as_ubyte(exposure.rescale_intensity(image))
                    img = Image.fromarray(np.rot90(image))
                    fileName_img = self.SaveFileName + "_" + str(self.stk.ev[i])+"eV_"+ str(i + 1) + "." + ext
                    img.save(fileName_img)
                QtWidgets.QApplication.restoreOverrideCursor()

            #Save all images in the stack
            if img_all_tif:
                QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))

                for i in range (self.stk.n_ev):
                    if self.parent.showflux:
                        #Show flux image
                        image = self.stk.absdata[:,:,i]
                    else:
                        #Show OD image
                        image = self.stk.od3d[:,:,i]

                    fileName_img = self.SaveFileName+"_" +str(self.stk.ev[i])+"eV_" +str(i+1)+".tif"
                    img = Image.fromarray(np.rot90(image))
                    img.save(fileName_img)

                QtWidgets.QApplication.restoreOverrideCursor()

            ext = 'pdf'
            suffix = "." + ext

            if spec_pdf:
                fileName_spec = self.SaveFileName+"_spectrum."+ext
                fig = self.parent.specfig
                fig.SaveFig(fileName_spec)

            if img_pdf:
                fileName_img = self.SaveFileName+"_" +str(self.stk.ev[self.parent.iev])+"eV."+ext
                fig = self.parent.absimgfig
                fig.SaveFig(fileName_img)

            if sp_csv:
                evdata = self.parent.specfig.pi.items[1].xData.tolist()
                fileName_spec = self.SaveFileName + "_spectrum.csv"
                if self.parent.button_showi0.isChecked():
                    name = "I0 data"
                    data = (self.parent.specfig.pi.items[-1].yData.tolist())
                else:
                    name = "ROI spectrum"
                    plot_num = len(self.parent.specfig.pi.items)
                    data = []
                    for i in range(plot_num - 1):
                        data.append(self.parent.specfig.pi.items[i + 1].yData.tolist())
                self.stk.write_csv(fileName_spec, evdata, data, cname=name)

            ext = 'svg'
            suffix = "." + ext

            if spec_svg:
                fileName_spec = self.SaveFileName+"_spectrum."+ext
                fig = self.parent.specfig
                fig.SaveFig(fileName_spec)

            if img_svg:
                fileName_img = self.SaveFileName+"_" +str(self.stk.ev[self.parent.iev])+"eV."+ext
                fig = self.parent.absimgfig
                fig.SaveFig(fileName_img)

            if img_tif:
                fileName_img = self.SaveFileName+"_" +str(self.stk.ev[self.parent.iev])+"eV.tif"
                if self.parent.showflux:
                    image = self.stk.absdata[:,:,self.parent.iev]
                else:
                    image = self.stk.od3d[:,:,self.parent.iev]
                img1 = Image.fromarray(np.rot90(image))
                img1.save(fileName_img)


        except IOError as e:
            if e.strerror:
                err = e.strerror
            else:
                err = e

            QtWidgets.QMessageBox.warning(self,'Error','Could not save file: %s' % err)


class SaveWinP2(QtWidgets.QDialog):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dialogsavepca.ui'), self)
        self.parent = parent
        self.setWindowTitle('Save')
        self.com = self.parent.common

        self.cb11.setChecked(True)
        self.cb21.setChecked(True)

        path, ext = os.path.splitext(self.com.filename) # currently empty?
        ext = ext[1:].lower()
        suffix = "." + ext
        path, fn = os.path.split(self.com.filename)
        self.filename = fn[:-len(suffix)]
        self.tc_savefn.setText(self.filename)

        self.path = self.com.path
        self.tc_savepath.setText(self.path)
        self.button_path.clicked.connect(self.OnBrowseDir)
        self.buttonBox.accepted.connect(self.OnSave)
        self.buttonBox.rejected.connect(self.close)


#----------------------------------------------------------------------
    def OnBrowseDir(self, evt):

        directory = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose a directory", self.path, QtWidgets.QFileDialog.ShowDirsOnly|QtWidgets.QFileDialog.ReadOnly)



        if directory == '':
            return

        directory = str(directory)
        self.com.path = directory

        self.path = directory

        self.tc_savepath.setText(self.path)



#----------------------------------------------------------------------
    def OnSave(self, evt=None):

        self.filename = str(self.tc_savefn.text())

        sp_pdf = self.cb11.isChecked()
        sp_png = self.cb12.isChecked()
        sp_svg = self.cb13.isChecked()
        sp_csv = self.cb14.isChecked()
        im_pdf = self.cb21.isChecked()
        im_png = self.cb22.isChecked()
        im_svg = self.cb23.isChecked()
        im_tif = self.cb25.isChecked()
        ev_pdf = self.cb31.isChecked()
        ev_png = self.cb32.isChecked()
        ev_svg = self.cb33.isChecked()

        self.close()
        self.parent.Save(self.filename, self.path,
                               spec_png = sp_png,
                               spec_pdf = sp_pdf,
                               spec_svg = sp_svg,
                               spec_csv = sp_csv,
                               img_png = im_png,
                               img_pdf = im_pdf,
                               img_svg = im_svg,
                               img_tif = im_tif,
                               evals_png = ev_png,
                               evals_pdf = ev_pdf,
                               evals_svg = ev_svg)


class SaveWinP3(QtWidgets.QDialog):

    def __init__(self, parent, page=None):
        QtWidgets.QWidget.__init__(self, parent)

        self.parent = parent
        self._ca_page = page  # optional explicit page; if None, uses self.parent.tab_clus


        self.resize(400, 300)
        self.setWindowTitle('Save')

        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)


        self.com = self.parent.common

        path, ext = os.path.splitext(self.com.filename)
        ext = ext[1:].lower()
        suffix = "." + ext
        path, fn = os.path.split(self.com.filename)
        filename = fn[:-len(suffix)]

        self.path = self.com.path
        self.filename = filename


        vboxtop = QtWidgets.QVBoxLayout()
        vboxtop.setContentsMargins(20,20,20,20)

        gridtop = QtWidgets.QGridLayout()
        gridtop.setVerticalSpacing(20)

        fontb = QtGui.QFont()
        fontb.setBold(True)

        st1 = QtWidgets.QLabel(self)
        st1.setText('Save')
        st1.setFont(fontb)
        st2 = QtWidgets.QLabel(self)
        st2.setText('.pdf')
        st2.setFont(fontb)
        st3 = QtWidgets.QLabel(self)
        st3.setText('.png')
        st3.setFont(fontb)
        st4 = QtWidgets.QLabel(self)
        st4.setText('.svg')
        st4.setFont(fontb)
        st5 = QtWidgets.QLabel(self)
        st5.setText('.csv')
        st5.setFont(fontb)
        st_tif = QtWidgets.QLabel(self)
        st_tif.setText('.tif (data)')
        st_tif.setFont(fontb)

        st6 = QtWidgets.QLabel(self)
        st6.setText('_spectrum')

        self.cb11 = QtWidgets.QCheckBox('', self)
        self.cb12 = QtWidgets.QCheckBox('', self)
        self.cb12.setChecked(True)
        self.cb13 = QtWidgets.QCheckBox('', self)
        self.cb14 = QtWidgets.QCheckBox('', self)

        st7 = QtWidgets.QLabel(self)
        st7.setText('_composite_images')

        self.cb21 = QtWidgets.QCheckBox('', self)
        self.cb22 = QtWidgets.QCheckBox('', self)
        self.cb22.setChecked(True)
        self.cb23 = QtWidgets.QCheckBox('', self)
        self.cb24 = QtWidgets.QCheckBox('', self)

        st8 = QtWidgets.QLabel(self)
        st8.setText('_individual_images')

        self.cb31 = QtWidgets.QCheckBox('', self)
        self.cb32 = QtWidgets.QCheckBox('', self)
        self.cb32.setChecked(True)
        self.cb33 = QtWidgets.QCheckBox('', self)
        self.cb34 = QtWidgets.QCheckBox('', self)

        st9 = QtWidgets.QLabel(self)
        st9.setText('_scatter_plots')

        self.cb41 = QtWidgets.QCheckBox('', self)
        self.cb42 = QtWidgets.QCheckBox('', self)
        self.cb42.setChecked(True)
        self.cb43 = QtWidgets.QCheckBox('', self)


        gridtop.addWidget(st1, 0, 0)
        gridtop.addWidget(st2, 0, 1)
        gridtop.addWidget(st3, 0, 2)
        gridtop.addWidget(st4, 0, 3)
        gridtop.addWidget(st5, 0, 4)
        gridtop.addWidget(st_tif, 0, 5)

        gridtop.addWidget(st6, 1, 0)
        gridtop.addWidget(self.cb11, 1, 1)
        gridtop.addWidget(self.cb12, 1, 2)
        gridtop.addWidget(self.cb13, 1, 3)
        gridtop.addWidget(self.cb14, 1, 4)

        gridtop.addWidget(st7, 2, 0)
        gridtop.addWidget(self.cb21, 2, 1)
        gridtop.addWidget(self.cb22, 2, 2)
        gridtop.addWidget(self.cb23, 2, 3)
        gridtop.addWidget(self.cb24, 2, 5)

        gridtop.addWidget(st8, 3, 0)
        gridtop.addWidget(self.cb31, 3, 1)
        gridtop.addWidget(self.cb32, 3, 2)
        gridtop.addWidget(self.cb33, 3, 3)
        gridtop.addWidget(self.cb34, 3, 5)

        gridtop.addWidget(st9, 4, 0)
        gridtop.addWidget(self.cb41, 4, 1)
        gridtop.addWidget(self.cb42, 4, 2)
        gridtop.addWidget(self.cb43, 4, 3)


        vboxtop.addStretch(1)
        vboxtop.addLayout(gridtop)
        vboxtop.addStretch(2)


        hbox0 = QtWidgets.QHBoxLayout()

        stf = QtWidgets.QLabel(self)
        stf.setText('Filename:\t')
        self.tc_savefn = QtWidgets.QLineEdit(self)
        self.tc_savefn.setText(self.filename)

        hbox0.addWidget(stf)
        hbox0.addWidget(self.tc_savefn)

        hbox1 = QtWidgets.QHBoxLayout()

        stp = QtWidgets.QLabel(self)
        stp.setText('Path:  \t')
        self.tc_savepath = QtWidgets.QLineEdit(self)
        self.tc_savepath.setReadOnly(True)
        self.tc_savepath.setText(self.path)
        self.tc_savepath.setMinimumWidth(100)
        hbox1.addWidget(stp)
        hbox1.addWidget(self.tc_savepath)

        button_path = QtWidgets.QPushButton('Browse...')
        button_path.clicked.connect(self.OnBrowseDir)
        hbox1.addWidget(button_path)


        hbox2 = QtWidgets.QHBoxLayout()
        button_save = QtWidgets.QPushButton('Save')
        button_save.clicked.connect(self.OnSave)
        hbox2.addWidget(button_save)

        button_cancel = QtWidgets.QPushButton('Cancel')
        button_cancel.clicked.connect(self.close)
        hbox2.addWidget(button_cancel)

        vboxtop.addLayout(hbox0)
        vboxtop.addLayout(hbox1)
        vboxtop.addStretch(1)
        vboxtop.addLayout(hbox2)



        self.setLayout(vboxtop)

#----------------------------------------------------------------------
    def OnBrowseDir(self, evt):

        directory = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose a directory", self.path, QtWidgets.QFileDialog.ShowDirsOnly|QtWidgets.QFileDialog.ReadOnly)



        if directory == '':
            return

        directory = str(directory)
        self.com.path = directory

        self.path = directory

        self.tc_savepath.setText(self.path)



#----------------------------------------------------------------------
    def OnSave(self, evt):

        self.filename = str(self.tc_savefn.text())

        sp_pdf = self.cb11.isChecked()
        sp_png = self.cb12.isChecked()
        sp_svg = self.cb13.isChecked()
        sp_csv = self.cb14.isChecked()
        im_pdf = self.cb21.isChecked()
        im_png = self.cb22.isChecked()
        im_svg = self.cb23.isChecked()
        im_tif = self.cb24.isChecked()
        indim_pdf = self.cb31.isChecked()
        indim_png = self.cb32.isChecked()
        indim_svg = self.cb33.isChecked()
        indim_tif = self.cb34.isChecked()
        scatt_pdf = self.cb41.isChecked()
        scatt_png = self.cb42.isChecked()
        scatt_svg = self.cb43.isChecked()

        self.close()
        if self._ca_page is not None:
            self._ca_page.Save_CA(self.filename, self.path,
                                             spec_png = sp_png,
                                             spec_pdf = sp_pdf,
                                             spec_svg = sp_svg,
                                             spec_csv = sp_csv,
                                             img_png = im_png,
                                             img_pdf = im_pdf,
                                             img_svg = im_svg,
                                             img_tif = im_tif,
                                             indimgs_png = indim_png,
                                             indimgs_pdf = indim_pdf,
                                             indimgs_svg = indim_svg,
                                             indimgs_tif = indim_tif,
                                             scatt_png = scatt_png,
                                             scatt_pdf = scatt_pdf,
                                             scatt_svg = scatt_svg)
        else:
            self.parent.tab_clus.Save_CA(self.filename, self.path,
                                             spec_png = sp_png,
                                             spec_pdf = sp_pdf,
                                             spec_svg = sp_svg,
                                             spec_csv = sp_csv,
                                             img_png = im_png,
                                             img_pdf = im_pdf,
                                             img_svg = im_svg,
                                             img_tif = im_tif,
                                             indimgs_png = indim_png,
                                             indimgs_pdf = indim_pdf,
                                             indimgs_svg = indim_svg,
                                             indimgs_tif = indim_tif,
                                             scatt_png = scatt_png,
                                             scatt_pdf = scatt_pdf,
                                             scatt_svg = scatt_svg)


class SaveWinP4(QtWidgets.QDialog):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)

        self.parent = parent


        self.resize(400, 300)
        self.setWindowTitle('Save')

        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)


        self.com = self.parent.common

        path, ext = os.path.splitext(self.com.filename)
        ext = ext[1:].lower()
        suffix = "." + ext
        path, fn = os.path.split(self.com.filename)
        filename = fn[:-len(suffix)]

        self.path = self.com.path
        self.filename = filename


        vboxtop = QtWidgets.QVBoxLayout()
        vboxtop.setContentsMargins(20,20,20,20)

        gridtop = QtWidgets.QGridLayout()
        gridtop.setVerticalSpacing(20)

        fontb = QtGui.QFont()
        fontb.setBold(True)

        st1 = QtWidgets.QLabel(self)
        st1.setText('Save')
        st1.setFont(fontb)
        st2 = QtWidgets.QLabel(self)
        st2.setText('.pdf')
        st2.setFont(fontb)
        st3 = QtWidgets.QLabel(self)
        st3.setText('.png')
        st3.setFont(fontb)
        st4 = QtWidgets.QLabel(self)
        st4.setText('.svg')
        st4.setFont(fontb)
        st5 = QtWidgets.QLabel(self)
        st5.setText('.csv')
        st5.setFont(fontb)
        st8 = QtWidgets.QLabel(self)
        st8.setText('.tif (data)')
        st8.setFont(fontb)

        st6 = QtWidgets.QLabel(self)
        st6.setText('_spectrum')

        self.cb11 = QtWidgets.QCheckBox('', self)
        self.cb11.setChecked(True)
        self.cb12 = QtWidgets.QCheckBox('', self)
        self.cb13 = QtWidgets.QCheckBox('', self)
        self.cb14 = QtWidgets.QCheckBox('', self)

        st7 = QtWidgets.QLabel(self)
        st7.setText('_images')

        self.cb21 = QtWidgets.QCheckBox('', self)
        self.cb21.setChecked(True)
        self.cb22 = QtWidgets.QCheckBox('', self)
        self.cb23 = QtWidgets.QCheckBox('', self)
        self.cb24 = QtWidgets.QCheckBox('', self)



        gridtop.addWidget(st1, 0, 0)
        gridtop.addWidget(st2, 0, 1)
        gridtop.addWidget(st3, 0, 2)
        gridtop.addWidget(st4, 0, 3)
        gridtop.addWidget(st5, 0, 4)
        gridtop.addWidget(st8, 0, 5)

        gridtop.addWidget(st6, 1, 0)
        gridtop.addWidget(self.cb11, 1, 1)
        gridtop.addWidget(self.cb12, 1, 2)
        gridtop.addWidget(self.cb13, 1, 3)
        gridtop.addWidget(self.cb14, 1, 4)

        gridtop.addWidget(st7, 2, 0)
        gridtop.addWidget(self.cb21, 2, 1)
        gridtop.addWidget(self.cb22, 2, 2)
        gridtop.addWidget(self.cb23, 2, 3)
        gridtop.addWidget(self.cb24, 2, 5)



        vboxtop.addStretch(1)
        vboxtop.addLayout(gridtop)
        vboxtop.addStretch(2)


        hbox0 = QtWidgets.QHBoxLayout()

        stf = QtWidgets.QLabel(self)
        stf.setText('Filename:\t')
        self.tc_savefn = QtWidgets.QLineEdit(self)
        self.tc_savefn.setText(self.filename)

        hbox0.addWidget(stf)
        hbox0.addWidget(self.tc_savefn)

        hbox1 = QtWidgets.QHBoxLayout()

        stp = QtWidgets.QLabel(self)
        stp.setText('Path:  \t')
        self.tc_savepath = QtWidgets.QLineEdit(self)
        self.tc_savepath.setReadOnly(True)
        self.tc_savepath.setText(self.path)
        self.tc_savepath.setMinimumWidth(100)
        hbox1.addWidget(stp)
        hbox1.addWidget(self.tc_savepath)

        button_path = QtWidgets.QPushButton('Browse...')
        button_path.clicked.connect(self.OnBrowseDir)
        hbox1.addWidget(button_path)


        hbox2 = QtWidgets.QHBoxLayout()
        button_save = QtWidgets.QPushButton('Save')
        button_save.clicked.connect(self.OnSave)
        hbox2.addWidget(button_save)

        button_cancel = QtWidgets.QPushButton('Cancel')
        button_cancel.clicked.connect(self.close)
        hbox2.addWidget(button_cancel)

        vboxtop.addLayout(hbox0)
        vboxtop.addLayout(hbox1)
        vboxtop.addStretch(1)
        vboxtop.addLayout(hbox2)



        self.setLayout(vboxtop)

#----------------------------------------------------------------------
    def OnBrowseDir(self, evt):

        directory = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose a directory", self.path, QtWidgets.QFileDialog.ShowDirsOnly|QtWidgets.QFileDialog.ReadOnly)



        if directory == '':
            return

        directory = str(directory)
        self.com.path = directory

        self.path = directory

        self.tc_savepath.setText(self.path)



#----------------------------------------------------------------------
    def OnSave(self, evt):

        self.filename = str(self.tc_savefn.text())

        sp_pdf = self.cb11.isChecked()
        sp_png = self.cb12.isChecked()
        sp_svg = self.cb13.isChecked()
        sp_csv = self.cb14.isChecked()
        im_pdf = self.cb21.isChecked()
        im_png = self.cb22.isChecked()
        im_svg = self.cb23.isChecked()
        im_tif = self.cb24.isChecked()


        self.close()
        self.parent.tab_spec.Save(self.filename, self.path,
                                         spec_png = sp_png,
                                         spec_pdf = sp_pdf,
                                         spec_svg = sp_svg,
                                         spec_csv = sp_csv,
                                         img_png = im_png,
                                         img_pdf = im_pdf,
                                         img_svg = im_svg,
                                         img_tif = im_tif)


class SaveWinP5(QtWidgets.QDialog):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)

        self.parent = parent


        self.resize(400, 300)
        self.setWindowTitle('Save')

        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)


        self.com = self.parent.common

        _path, ext = os.path.splitext(self.com.filename)
        ext = ext[1:].lower()
        suffix = "." + ext
        _path, fn = os.path.split(self.com.filename)
        filename = fn[:-len(suffix)]

        self.path = self.com.path
        self.filename = filename


        vboxtop = QtWidgets.QVBoxLayout()
        vboxtop.setContentsMargins(20,20,20,20)

        gridtop = QtWidgets.QGridLayout()
        gridtop.setVerticalSpacing(20)

        fontb = QtGui.QFont()
        fontb.setBold(True)

        st1 = QtWidgets.QLabel(self)
        st1.setText('Save')
        st1.setFont(fontb)
        st2 = QtWidgets.QLabel(self)
        st2.setText('.pdf')
        st2.setFont(fontb)
        st3 = QtWidgets.QLabel(self)
        st3.setText('.png')
        st3.setFont(fontb)
        st4 = QtWidgets.QLabel(self)
        st4.setText('.svg')
        st4.setFont(fontb)
        st5 = QtWidgets.QLabel(self)
        st5.setText('.csv')
        st5.setFont(fontb)
        st9 = QtWidgets.QLabel(self)
        st9.setText('.tif (data)')
        st9.setFont(fontb)

        st6 = QtWidgets.QLabel(self)
        st6.setText('_spectrum')

        self.cb11 = QtWidgets.QCheckBox('', self)
        self.cb11.setChecked(True)
        self.cb12 = QtWidgets.QCheckBox('', self)
        self.cb13 = QtWidgets.QCheckBox('', self)
        self.cb14 = QtWidgets.QCheckBox('', self)

        st7 = QtWidgets.QLabel(self)
        st7.setText('_map')

        self.cb21 = QtWidgets.QCheckBox('', self)
        self.cb21.setChecked(True)
        self.cb22 = QtWidgets.QCheckBox('', self)
        self.cb23 = QtWidgets.QCheckBox('', self)
        self.cb24 = QtWidgets.QCheckBox('', self)

        st8 = QtWidgets.QLabel(self)
        st8.setText('_costfunction')
        self.cb31 = QtWidgets.QCheckBox('', self)
        self.cb32 = QtWidgets.QCheckBox('', self)
        self.cb33 = QtWidgets.QCheckBox('', self)


        gridtop.addWidget(st1, 0, 0)
        gridtop.addWidget(st2, 0, 1)
        gridtop.addWidget(st3, 0, 2)
        gridtop.addWidget(st4, 0, 3)
        gridtop.addWidget(st5, 0, 4)
        gridtop.addWidget(st9, 0, 5)

        gridtop.addWidget(st6, 1, 0)
        gridtop.addWidget(self.cb11, 1, 1)
        gridtop.addWidget(self.cb12, 1, 2)
        gridtop.addWidget(self.cb13, 1, 3)
        gridtop.addWidget(self.cb14, 1, 4)

        gridtop.addWidget(st7, 2, 0)
        gridtop.addWidget(self.cb21, 2, 1)
        gridtop.addWidget(self.cb22, 2, 2)
        gridtop.addWidget(self.cb23, 2, 3)
        gridtop.addWidget(self.cb24, 2, 5)

        gridtop.addWidget(st8, 3, 0)
        gridtop.addWidget(self.cb31, 3, 1)
        gridtop.addWidget(self.cb32, 3, 2)
        gridtop.addWidget(self.cb33, 3, 3)

        vboxtop.addStretch(1)
        vboxtop.addLayout(gridtop)
        vboxtop.addStretch(2)


        hbox0 = QtWidgets.QHBoxLayout()

        stf = QtWidgets.QLabel(self)
        stf.setText('Filename:\t')
        self.tc_savefn = QtWidgets.QLineEdit(self)
        self.tc_savefn.setText(self.filename)

        hbox0.addWidget(stf)
        hbox0.addWidget(self.tc_savefn)

        hbox1 = QtWidgets.QHBoxLayout()

        stp = QtWidgets.QLabel(self)
        stp.setText('Path:  \t')
        self.tc_savepath = QtWidgets.QLineEdit(self)
        self.tc_savepath.setReadOnly(True)
        self.tc_savepath.setText(self.path)
        self.tc_savepath.setMinimumWidth(100)
        hbox1.addWidget(stp)
        hbox1.addWidget(self.tc_savepath)

        button_path = QtWidgets.QPushButton('Browse...')
        button_path.clicked.connect(self.OnBrowseDir)
        hbox1.addWidget(button_path)


        hbox2 = QtWidgets.QHBoxLayout()
        button_save = QtWidgets.QPushButton('Save')
        button_save.clicked.connect(self.OnSave)
        hbox2.addWidget(button_save)

        button_cancel = QtWidgets.QPushButton('Cancel')
        button_cancel.clicked.connect(self.close)
        hbox2.addWidget(button_cancel)

        vboxtop.addLayout(hbox0)
        vboxtop.addLayout(hbox1)
        vboxtop.addStretch(1)
        vboxtop.addLayout(hbox2)



        self.setLayout(vboxtop)

#----------------------------------------------------------------------
    def OnBrowseDir(self, evt):

        directory = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose a directory", self.path, QtWidgets.QFileDialog.ShowDirsOnly|QtWidgets.QFileDialog.ReadOnly)



        if directory == '':
            return

        directory = str(directory)
        self.com.path = directory

        self.path = directory

        self.tc_savepath.setText(self.path)



#----------------------------------------------------------------------
    def OnSave(self, evt):

        self.filename = str(self.tc_savefn.text())

        sp_pdf = self.cb11.isChecked()
        sp_png = self.cb12.isChecked()
        sp_svg = self.cb13.isChecked()
        sp_csv = self.cb14.isChecked()
        im_pdf = self.cb21.isChecked()
        im_png = self.cb22.isChecked()
        im_svg = self.cb23.isChecked()
        im_tif = self.cb24.isChecked()
        cf_pdf = self.cb31.isChecked()
        cf_png = self.cb32.isChecked()
        cf_svg = self.cb33.isChecked()

        self.close()
        self.parent.tab_nnma.Save(self.filename, self.path,
                               spec_png = sp_png,
                               spec_pdf = sp_pdf,
                               spec_svg = sp_svg,
                               spec_csv = sp_csv,
                               map_png = im_png,
                               map_pdf = im_pdf,
                               map_svg = im_svg,
                               map_tif = im_tif,
                               costf_png = cf_png,
                               costf_pdf = cf_pdf,
                               costf_svg = cf_svg)
