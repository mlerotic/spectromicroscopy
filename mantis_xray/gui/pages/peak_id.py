
import os
import numpy as np
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtCore import Qt
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from ...core.constants import PlotH, PlotW

class PagePeakID(QtWidgets.QWidget):
    def __init__(self, common, data_struct, stack, anlz):
        super(PagePeakID, self).__init__()

        self.initUI(common, data_struct, stack, anlz)

#----------------------------------------------------------------------
    def initUI(self, common, data_struct, stack, anlz):


        self.com = common
        self.data_struct = data_struct
        self.stk = stack
        self.anlz = anlz



        pw = PlotW*0.8
        ph = PlotH*0.8


        self.peak_engs = [285.0, 286.5, 288.5, 289.5, 290.5]
        self.peak_names = ['aromatic', 'phenolic', 'carboxyl', 'alkyl', 'carbonyl']

        self.i_spectrum = 1
        self.i_eng = 0

        self.n_spectra = 0

        self.spectrum_loaded = 0


        #panel 1
        vbox1 = QtWidgets.QVBoxLayout()

        self.tc_1 = QtWidgets.QLabel(self)
        self.tc_1.setText("X-ray Spectrum")
        hbox11 = QtWidgets.QHBoxLayout()

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.kespecfig = Figure((pw, ph))
        self.KESpecPan = FigureCanvas(self.kespecfig)
        self.KESpecPan.mpl_connect('button_press_event', self.OnPointSpectrum)
        fbox.addWidget(self.KESpecPan)
        frame.setLayout(fbox)

        self.slider_spec = QtWidgets.QScrollBar(QtCore.Qt.Vertical)
        self.slider_spec.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_spec.setEnabled(False)
        self.slider_spec.valueChanged[int].connect(self.OnSScroll)
        self.slider_spec.setRange(1, 5)

        hbox11.addWidget(frame)
        hbox11.addWidget(self.slider_spec)


        vbox1.addStretch(1)
        vbox1.addWidget(self.tc_1)
        vbox1.addLayout(hbox11)
        vbox1.addStretch(1)




        #panel 2
        sizer2 = QtWidgets.QGroupBox('Load spectra')
        vbox2 = QtWidgets.QVBoxLayout()

        self.button_addclspec = QtWidgets.QPushButton('Add cluster spectra')
        self.button_addclspec.clicked.connect( self.OnAddClusterSpectra)
        self.button_addclspec.setEnabled(False)
        vbox2.addWidget(self.button_addclspec)


        self.button_loadtspec = QtWidgets.QPushButton('Load spectra')
        self.button_loadtspec.clicked.connect(self.OnSpecFromFile)
        vbox2.addWidget(self.button_loadtspec)


        self.button_save = QtWidgets.QPushButton('Save images...')
        self.button_save.clicked.connect( self.OnSave)
        self.button_save.setEnabled(False)
        vbox2.addWidget(self.button_save)

        sizer2.setLayout(vbox2)


        #panel 3
        vbox3 = QtWidgets.QVBoxLayout()

        t1 = QtWidgets.QLabel(self)
        t1.setText("Peak ID Energies")

        self.lc_1 = QtWidgets.QListWidget()
        self.lc_1.itemClicked.connect(self.OnEngListClick)
        self.lc_1.setMinimumSize(200, 400)

        vbox3.addWidget(t1)
        vbox3.addWidget(self.lc_1)


        #panel 4
        vbox4 = QtWidgets.QVBoxLayout()

        self.tc_imageeng = QtWidgets.QLabel(self)
        self.tc_imageeng.setText("Image at peak energy: ")
        vbox4.addWidget(self.tc_imageeng)

        hbox41 = QtWidgets.QHBoxLayout()

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.absimgfig = Figure((ph, ph))
        self.AbsImagePanel = FigureCanvas(self.absimgfig)

        fbox.addWidget(self.AbsImagePanel)
        frame.setLayout(fbox)
        hbox41.addWidget(frame)

        self.slider_eng = QtWidgets.QScrollBar(QtCore.Qt.Vertical)
        self.slider_eng.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider_eng.valueChanged[int].connect(self.OnScrollEng)
        self.slider_eng.setRange(0, len(self.peak_names)-1)
        hbox41.addWidget(self.slider_eng)
        hbox41.addStretch(1)
        vbox4.addLayout(hbox41)





        vboxtop1 = QtWidgets.QVBoxLayout()

        vboxtop1.addStretch(1)
        vboxtop1.addWidget(sizer2)
        vboxtop1.addStretch(1)
        vboxtop1.addLayout(vbox3)

        vboxtop2 = QtWidgets.QVBoxLayout()
        vboxtop2.addStretch(1)
        vboxtop2.addLayout(vbox1)
        vboxtop2.addStretch(1)
        vboxtop2.addLayout(vbox4)


        hboxtop = QtWidgets.QHBoxLayout()
        hboxtop.addStretch(1)
        hboxtop.addLayout(vboxtop1)
        hboxtop.addStretch(1)
        hboxtop.addLayout(vboxtop2)
        hboxtop.addStretch(1)


        hboxtop.setContentsMargins(20,20,20,20)
        self.setLayout(hboxtop)



        self.ShowPeakEngs()


#----------------------------------------------------------------------
    def OnSpecFromFile(self, event):
        #ToDo: Clean-up needed. One dialog for all pages is sufficient.

        try:
        #if True:
            wildcard = "ASCII files (*.csv *.txt *.xas);;CSV spectrum (*.csv);;TXT spectrum (*.txt);;XAS spectrum (*.xas)"
            directory = self.com.path
            filepath, _filter = QtWidgets.QFileDialog.getOpenFileName(self, 'Choose Spectrum file', directory, wildcard)

            filepath = str(filepath)
            if filepath == '':
                return

            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))

            self.anlz.load_xraypeakfit_spectra(filename=filepath)
            if not self.anlz.xrayfitsp_to_add:
                raise IndexError
            if self.window().tab_pfit != None:
                self.com.xpf_loaded = 1
                for i in range(self.anlz.xrayfitsp_to_add):
                    self.window().tab_pfit.spectrumfitted.append(0)
                    self.window().tab_pfit.firstinit.append(1)
                    self.window().tab_pfit.initdone.append(0)
                    self.window().tab_pfit.nsteps.append(1)
                    self.window().tab_pfit.npeaks.append(4)
                    self.window().tab_pfit.fits.append(0)
                    self.window().tab_pfit.fits_sep.append(0)

                self.window().tab_pfit.i_spec = self.anlz.n_xrayfitsp
                self.window().tab_pfit.slider_spec.setMaximum(self.anlz.n_xrayfitsp)
                self.window().tab_pfit.slider_spec.setValue(self.window().tab_pfit.i_spec)
                self.window().tab_pfit.nstepsspin.setValue(self.window().tab_pfit.nsteps[self.window().tab_pfit.i_spec-1])
                self.window().tab_pfit.ngaussspin.setValue(self.window().tab_pfit.npeaks[self.window().tab_pfit.i_spec-1])
                self.window().tab_pfit.loadSpectrum()
                self.window().tab_pfit.updatewidgets()

            self.spectrum_loaded = 1
            self.updatewidgets()


            QtWidgets.QApplication.restoreOverrideCursor()
        except IndexError:
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, 'Error', 'Spectrum file not loaded. Data invalid or zero.')

        except TypeError:
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, 'Error', 'Spectrum file not loaded. Unknown format.')

        except Exception as error:
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, 'Error', "Unexpected Error: "+ str(error) +" "+ str(type(error)))

        self.window().refresh_widgets()

        self.slider_spec.setMaximum(self.n_spectra)
        self.slider_spec.setEnabled(True)

        self.ShowODPlot()


#----------------------------------------------------------------------
    def OnAddClusterSpectra(self, event):

        try:

            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))


            for i in range(self.anlz.nclusters):
                self.anlz.load_xraypeakfit_clusterspectrum(i)

                if self.window().tab_pfit != None:
                    self.window().tab_pfit.spectrumfitted.append(0)
                    self.window().tab_pfit.firstinit.append(1)
                    self.window().tab_pfit.initdone.append(0)
                    self.window().tab_pfit.nsteps.append(1)
                    self.window().tab_pfit.npeaks.append(4)
                    self.window().tab_pfit.fits.append(0)
                    self.window().tab_pfit.fits_sep.append(0)


            self.spectrum_loaded = 1


            self.com.xpf_loaded = 1
            if self.window().tab_pfit != None:
                self.window().tab_pfit.i_spec = self.anlz.n_xrayfitsp
                self.window().tab_pfit.slider_spec.setMaximum(self.anlz.n_xrayfitsp)
                self.window().tab_pfit.slider_spec.setValue(self.window().tab_pfit.i_spec)
                self.window().tab_pfit.nstepsspin.setValue(self.window().tab_pfit.nsteps[self.window().tab_pfit.i_spec-1])
                self.window().tab_pfit.ngaussspin.setValue(self.window().tab_pfit.npeaks[self.window().tab_pfit.i_spec-1])
                self.window().tab_pfit.loadSpectrum()
                self.window().tab_pfit.updatewidgets()


            QtWidgets.QApplication.restoreOverrideCursor()

        except:
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, 'Error', 'Cluster spectra not loaded.')




        self.window().refresh_widgets()

        self.updatewidgets()

        self.ShowODPlot()
        self.ShowImage()

#----------------------------------------------------------------------
    def OnSScroll(self, value):


        sel = value
        self.i_spectrum = sel


        self.ShowODPlot()
        self.ShowImage()

#----------------------------------------------------------------------
    def OnScrollEng(self, value):
        self.i_eng = value


        self.ShowImage()
        self.ShowODPlot()


#----------------------------------------------------------------------
    def OnPointSpectrum(self, evt):
        x = evt.xdata

        if (self.com.i0_loaded != 0) or (self.spectrum_loaded == 1):
            if x < self.stk.ev[0]:
                sel_ev = 0
            elif x > self.stk.ev[self.stk.n_ev-1]:
                sel_ev = self.stk.n_ev-1
            else:
                indx = np.abs(self.stk.ev - x).argmin()
                sel_ev = indx


            self.i_eng=(np.abs(self.peak_engs-self.stk.ev[sel_ev])).argmin()

            self.ShowImage()
            self.ShowODPlot()

            self.slider_eng.setValue(self.i_eng)

#----------------------------------------------------------------------
    def ShowODPlot(self):

        if (self.spectrum_loaded == 0) and (self.com.i0_loaded == 0):
            return

        if (self.com.i0_loaded == 0):
            spectrum = self.anlz.xrayfitspectra[self.i_spectrum-1, :]
            self.tc_1.setText("X-ray Spectrum: " +
                           self.anlz.xfspec_names[self.i_spectrum-1])
        else:
            if (self.i_spectrum == 1):
                spectrum = self.stk.od3d.sum(axis=0)
                spectrum = spectrum.sum(axis=0)/(self.stk.n_rows*self.stk.n_cols)

                self.tc_1.setText("X-ray Spectrum: Average OD")
            else:
                spectrum = self.anlz.xrayfitspectra[self.i_spectrum-2, :]
                self.tc_1.setText("X-ray Spectrum: " +
                               self.anlz.xfspec_names[self.i_spectrum-2])


        fig = self.kespecfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        axes = fig.gca()


        specplot = axes.plot(self.stk.ev,spectrum)

        for i in range(len(self.peak_engs)):
            if (self.peak_engs[i]>self.stk.ev[0]) and (self.peak_engs[i]<self.stk.ev[-1]):
                axes.axvline(x=self.peak_engs[i], color = 'g', alpha=0.5)

        if (self.peak_engs[self.i_eng]>self.stk.ev[0]) and (self.peak_engs[self.i_eng]<self.stk.ev[-1]):
            axes.axvline(x=self.peak_engs[self.i_eng], color = 'g', alpha=1.0)

        axes.set_xlabel('Photon Energy [eV]')
        axes.set_ylabel('Optical Density')

        self.KESpecPan.draw()

#----------------------------------------------------------------------
    def OnEngListClick(self):
        item = self.lc_1.currentRow()

        self.i_eng = item

        #self.ShowPeakEngs()
        self.ShowImage()
        self.ShowODPlot()

#----------------------------------------------------------------------
    def ShowPeakEngs(self):

        self.lc_1.clear()

        for i in range(len(self.peak_engs)):
            self.lc_1.addItem('{0:8.2f}'.format(self.peak_engs[i])+'\t'+self.peak_names[i])


#----------------------------------------------------------------------
    def ShowImage(self):

        if self.com.i0_loaded == 0:
            return

        iev=(np.abs(self.stk.ev-self.peak_engs[self.i_eng])).argmin()
        image = self.stk.absdata[:,:,int(iev)].copy()

        fig = self.absimgfig
        fig.clf()
        fig.add_axes(((0.0,0.0,1.0,1.0)))
        axes = fig.gca()
        fig.patch.set_alpha(1.0)

        im = axes.imshow(np.rot90(image), cmap=matplotlib.colormaps["gray"])


        #Show Scale Bar
        startx = int(self.stk.n_rows*0.05)
        starty = int(self.stk.n_cols*0.05)+self.stk.scale_bar_pixels_y
        um_string = ' $\mathrm{\mu m}$'
        microns = '$'+self.stk.scale_bar_string+' $'+um_string
        axes.text(self.stk.scale_bar_pixels_x+startx+1,starty+1, microns, horizontalalignment='left', verticalalignment='center',
                  color = 'black', fontsize=14)
        #Matplotlib has flipped scales so I'm using rows instead of cols!
        p = matplotlib.patches.Rectangle((startx,starty), self.stk.scale_bar_pixels_x, self.stk.scale_bar_pixels_y,
                               color = 'black', fill = True)
        axes.add_patch(p)


        axes.axis("off")
        self.AbsImagePanel.draw()

        self.tc_imageeng.setText('Image at peak energy: {0:5.2f} eV'.format(float(self.stk.ev[iev])))

        self.lc_1.setCurrentRow(self.i_eng)

#----------------------------------------------------------------------
    def updatewidgets(self):

        if self.spectrum_loaded:
            self.n_spectra = self.anlz.n_xrayfitsp
        else:
            self.n_spectra = 0

        if self.com.i0_loaded == 1:
            self.n_spectra += 1

        self.i_spectrum = self.n_spectra

        self.slider_spec.setMaximum(self.n_spectra)
        self.slider_spec.setEnabled(True)
        self.slider_spec.setValue(self.i_spectrum)
        self.ShowODPlot()
        self.ShowImage()


#----------------------------------------------------------------------
    def OnSave(self, event):


        #Save images
        wildcard = "Portable Network Graphics (*.png);;Adobe PDF Files (*.pdf);;"

        fileName, _filter = QtWidgets.QFileDialog.getSaveFileName(self, 'Save OD Plot with Key Energies', '', wildcard)

        fileName = str(fileName)
        if fileName == '':
            return

        path, ext = os.path.splitext(fileName)
        ext = ext[1:].lower()


        if ext != 'png' and ext != 'pdf':
            error_message = (
                  'Only the PNG and PDF image formats are supported.\n'
                 'A file extension of `png\' or `pdf\' must be used.')

            QtWidgets.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % error_message)
            return



        try:

            matplotlib.rcParams['pdf.fonttype'] = 42

            fig = self.kespecfig
            fig.savefig(fileName)


        except IOError as e:
            if e.strerror:
                err = e.strerror
            else:
                err = e


            QtWidgets.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err)



        return

