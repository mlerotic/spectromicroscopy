import os
import numpy as np
from PyQt5 import QtWidgets, QtCore, QtGui, uic
from PyQt5.QtCore import pyqtSignal
import pyqtgraph as pg
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable


class SpectralImageMap(QtWidgets.QDialog):
    def __init__(self, parent, common, stack):
        QtWidgets.QWidget.__init__(self, parent)
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'showspectralroi.ui'), self)
        self.parent = parent
        self.stack = stack
        self.com = common
        self.iev = 0
        self.itheta = 0

        self.setWindowTitle('Spectral Image Map')
        self.parent.pbSelfromSpec.setEnabled(False)
        self.button_close.clicked.connect(self.closeEvent)

        if self.com.stack_loaded == 1:
            self.SetupPlot()

    def OnSelectionChanged(self):
        self.region_i0.blockSignals(True)
        if self.idx_selected:
            self.region_i0.setRegion([self.stack.ev[min(self.idx_selected)], self.stack.ev[max(self.idx_selected)]])
            self.region_i0.blockSignals(False)
        return

    def GenerateSpectrum(self, evselection):
        if self.com.i0_loaded == 1:
            if self.com.stack_4d == 1:
                total = self.stack.od4d[:, :, :, int(self.itheta)].copy()
            else:
                total = self.stack.od3d[:, :, :].copy()

        else:
            if self.com.stack_4d == 1:
                # t = [self.stack.theta[i] for i in self.thetaidx_selected]
                # self.label_theta_range.setText(
                #     "Theta range: [ " + str(min(t, default=0)) + "° .. " + str(
                #         max(t, default=0)) + "° ], # values: " + str(
                #         len(t)))
                total = self.stack.stack4D[:, :, :, int(self.itheta)].copy()
            else:
                total = self.stack.absdata[:, :, :].copy()
        total = total.sum(axis=(0,1)) #/ (int(self.box.size().x()) * int(self.box.size().y()))
        x = self.stack.ev
        y = total
        return (x, y)

    def SetupPlot(self):
        x, y = self.GenerateSpectrum(list(range(self.stack.n_ev)))
        self.spectrum_plotwidget.setBackground("w")

        self.region_i0 = pg.LinearRegionItem(brush=QtGui.QColor('#88beaed4'),hoverBrush=QtGui.QColor('#ccbeaed4'), bounds=[np.min(x), np.max(x)])
        self.region_i0.setZValue(10)
        self.region_i = pg.LinearRegionItem(brush=QtGui.QColor('#887fc97f'),hoverBrush=QtGui.QColor('#cc7fc97f'), bounds=[np.min(x), np.max(x)])
        self.region_i.setZValue(11)

        plot = self.spectrum_plotwidget
        plot.setBackground("w")
        plot.addItem(self.region_i0, ignoreBounds=False)
        plot.addItem(self.region_i, ignoreBounds=False)
        plot.setMouseEnabled(x=True, y=True)
        plot.showGrid(y=True)

        plot.showAxis("top", show=True)
        plot.showAxis("right", show=True)
        by = plot.getAxis("right")
        bx = plot.getAxis("top")
        by.setStyle(showValues=False, tickLength=0)
        bx.setStyle(showValues=False, tickLength=0)
        ay = plot.getAxis("left")
        ax = plot.getAxis("bottom")

        ax.setLabel(text="Photon energy [eV]")
        if self.com.i0_loaded:
            ay.setLabel(text="Sum of optical densities")
        else:
            ay.setLabel(text="Photon flux [cps]")

        self.plotitem = plot.plot(x, y, pen=pg.mkPen(color="b", width=2),symbolBrush=(0,0,255),symbolPen=('k'))
        tristate_list = [item for row in self.stack.shifts for item in row][1::3]
        self.ColorizeScatterDots(tristate_list)
        self.plotitem.setZValue(1)
        self.region_i0.setRegion((min(x), min(x)))
        self.region_i.setRegion((max(x), max(x)))
        self.label_i0 = pg.InfLineLabel(self.region_i0.lines[0], text="I0 / pre-edge", movable=False,angle=90, position=0.05, anchors=(0,0),color=(0, 0, 0))
        self.label_i =  pg.InfLineLabel(self.region_i.lines[1], text="I / on-edge", movable=False,angle=90, position=0.75, anchors=(0,0),color=(0, 0, 0))
        self.region_i0.sigRegionChanged.connect(lambda region: self.UpdateZvalue(region))
        self.region_i.sigRegionChanged.connect(lambda region: self.UpdateZvalue(region))
        self.region_i0.sigRegionChangeFinished.connect(lambda region: self.UpdateSelection(region))
        self.region_i.sigRegionChangeFinished.connect(lambda region: self.UpdateSelection(region))
    # ----------------------------------------------------------------------
    # Colorize scatter dots according to the pre/on-edge color
    def ColorizeScatterDots(self,list):
        brushes = [QtGui.QColor('#ffbeaed4') if val == -1 else QtGui.QColor('blue') if val == 0 else QtGui.QColor('#ff7fc97f') for val in list]
        self.plotitem.scatter.setBrush(brushes,update=True)
    # ----------------------------------------------------------------------
    # Bring active region to front
    def UpdateZvalue(self,region):
        otherregion = {self.region_i: self.region_i0, self.region_i0: self.region_i}
        region.setZValue(11)
        otherregion[region].setZValue(10)
    # ----------------------------------------------------------------------
    def getDataClosestToRegion(self, region, plotitem, snapregion=True):
        otherregion = {self.region_i : self.region_i0, self.region_i0: self.region_i}
        otherregion = otherregion[region]
        minidx, maxidx = region.getRegion()
        data = plotitem.getData()[0]
        index = lambda x: np.argmin(np.abs(data - x))
        minidx = index(minidx)+1 if minidx > data[index(minidx)] else index(minidx)
        maxidx = index(maxidx)-1 if maxidx < data[index(maxidx)] else index(maxidx)
        mindata = (data[minidx] + data[max((minidx-1,0))])/2
        maxdata = (data[maxidx] + data[min((maxidx+1,np.argmax(data)))])/2
        #push regions. regions may not overlap!
        if region.zValue() > otherregion.zValue() : #choose the active/just dragged region
            othermin, othermax = otherregion.getRegion()
            if othermin <= mindata <= othermax:
                otherregion.setBounds((othermin, mindata))
                #otherregion.setRegion([othermin, mindata])
            elif othermin <= maxdata <= othermax:
                otherregion.setBounds((maxdata, othermax))
        if snapregion:
            region.blockSignals(True)
            otherregion.setBounds((min(data), max(data)))
            region.setRegion([mindata, maxdata])  # snap region to data points
            region.blockSignals(False)
            #self.parent.OnSelectionChanged()
        return minidx, maxidx, mindata, maxdata

    def UpdateSelection(self,region):
        otherregion = {self.region_i : self.region_i0, self.region_i0: self.region_i}
        otherregion = otherregion[region]
        if region == self.region_i0:
            mini, maxi, mindi, maxdi =      self.getDataClosestToRegion(otherregion,self.plotitem)
            mini0, maxi0, mindi0, maxdi0 =  self.getDataClosestToRegion(region,self.plotitem)
        else:
            mini, maxi, mindi, maxdi =      self.getDataClosestToRegion(region,self.plotitem)
            mini0, maxi0, mindi0, maxdi0 =  self.getDataClosestToRegion(otherregion,self.plotitem)
        #self.region_i0.blockSignals(True)
        #self.region_i.blockSignals(True)
        #self.region_i0.setBounds((min(self.plotitem.getData()[0]),mindi))
        #self.region_i.setBounds((maxdi0,max(self.plotitem.getData()[0])))
        #self.region_i0.blockSignals(False)
        #self.region_i.blockSignals(False)

        qlist = self.parent.MapSelectWidget1
        for row in range(qlist.count()):
            if row in range(mini0, maxi0+1):
                self.stack.shifts[row][1]= -1
                qlist.item(row).setBackground(QtGui.QColor('#beaed4'))
            elif row in range(mini, maxi+1):
                self.stack.shifts[row][1]= 1
                qlist.item(row).setBackground(QtGui.QColor('#7fc97f'))
            else:
                self.stack.shifts[row][1]= 0
                qlist.item(row).setBackground(QtGui.QColor(0, 0, 0, 0))
        tristate_list = [item for row in self.stack.shifts for item in row][1::3]
        self.ColorizeScatterDots(tristate_list)
        self.parent.OnSelectionChanged()
    # ----------------------------------------------------------------------
    def closeEvent(self, event):
        self.close()
        self.parent.pbSelfromSpec.setEnabled(True)
# ----------------------------------------------------------------------
class SpectralROI(QtWidgets.QDialog):

    def __init__(self, parent,  common, stack):
        QtWidgets.QWidget.__init__(self, parent)

        self.parent = parent

        self.stack = stack
        self.com = common

        self.resize(630, 700)
        self.setWindowTitle('Spectral Regions of Interest')

        pal = QtGui.QPalette()
        self.setAutoFillBackground(True)
        pal.setColor(QtGui.QPalette.Window,QtGui.QColor('white'))
        self.setPalette(pal)

        self.stack = stack
        self.com = common


        self.imin = 0
        self.imax = 0
        self.i0min = 0
        self.i0max = 0
        self.iselected = 0
        self.i0selected = 0


        self.odtotal = self.stack.od3d.sum(axis=0)
        self.odtotal = self.odtotal.sum(axis=0)/(self.stack.n_rows*self.stack.n_cols)

        self.image_i0 = np.zeros((self.stack.n_cols, self.stack.n_rows))
        self.image_i = np.zeros((self.stack.n_cols, self.stack.n_rows))
        self.odthickmap = np.zeros((self.stack.n_cols, self.stack.n_rows))



        vbox = QtWidgets.QVBoxLayout()
        text = QtWidgets.QLabel(self)
        text.setText('First select I0 region below the edge, then select I region above the edge:')
        vbox.addWidget(text)

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.specfig = Figure((6.0, 4.2))
        self.SpectrumPanel = FigureCanvas(self.specfig)
        self.SpectrumPanel.setParent(self)
        self.SpectrumPanel.mpl_connect('button_press_event', self.OnSelection1)
        self.SpectrumPanel.mpl_connect('button_release_event', self.OnSelection2)

        fbox.addWidget(self.SpectrumPanel)
        frame.setLayout(fbox)
        vbox.addWidget(frame)


        hbox2 = QtWidgets.QHBoxLayout()


        sizer2 =QtWidgets.QGroupBox('Selected Spectral Regions')
        sizer2.setMinimumWidth(350)
        vbox2 = QtWidgets.QVBoxLayout()
        text = QtWidgets.QLabel(self)
        text.setText('I Selection (red): ')
        self.textctrl1 = QtWidgets.QLabel(self)
        self.textctrl1.setText('[  ]' )
        vbox2.addWidget(text)
        vbox2.addWidget(self.textctrl1)

        text = QtWidgets.QLabel(self)
        text.setText('I0 Selection (green): ')
        self.textctrl2 = QtWidgets.QLabel(self)
        self.textctrl2.setText('[  ]' )
        vbox2.addWidget(text)
        vbox2.addWidget(self.textctrl2)
        sizer2.setLayout(vbox2)

        text = QtWidgets.QLabel(self)
        text.setText('Optical density map')
        vbox.addWidget(text)

        frame = QtWidgets.QFrame()
        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Sunken)
        fbox = QtWidgets.QHBoxLayout()

        self.odmfig = Figure((2.4,2.4))
        self.ODMImagePanel = FigureCanvas(self.odmfig)
        self.ODMImagePanel.setParent(self)

        fbox.addWidget(self.ODMImagePanel, stretch=0)
        frame.setLayout(fbox)
        hbox2.addWidget(frame, stretch=0)
        hbox2.addStretch(1)
        hbox2.addWidget(sizer2)

        vbox.addLayout(hbox2)

        hbox = QtWidgets.QHBoxLayout()

        button_save = QtWidgets.QPushButton('Save')
        button_save.clicked.connect(self.OnSave)
        hbox.addWidget(button_save)

        button_cancel = QtWidgets.QPushButton('Dismiss')
        button_cancel.clicked.connect(self.close)
        hbox.addWidget(button_cancel)

        vbox.addLayout(hbox)


        self.setLayout(vbox)

        self.draw_spectrum()


#----------------------------------------------------------------------
    def draw_spectrum(self):


        fig = self.specfig
        fig.clf()
        fig.add_axes((0.15,0.15,0.75,0.75))
        self.axes = fig.gca()


        specplot = self.axes.plot(self.stack.ev,self.odtotal)

        if self.i0selected == 1:
            self.axes.axvspan(self.stack.ev[self.i0min], self.stack.ev[self.i0max], facecolor='g', alpha=0.5)

        if self.iselected == 1:
            self.axes.axvspan(self.stack.ev[self.imin], self.stack.ev[self.imax], facecolor='r', alpha=0.5)


        self.axes.set_xlabel('Photon Energy [eV]')
        self.axes.set_ylabel('Optical Density')


        self.SpectrumPanel.draw()


#----------------------------------------------------------------------
    def draw_image(self):


        fig = self.odmfig
        fig.clf()
        fig.add_axes((0.02,0.02,0.96,0.96))

        axes = fig.gca()
        divider = make_axes_locatable(axes)
        axcb = divider.new_horizontal(size="3%", pad=0.03)

        fig.add_axes(axcb)

        axes.set_position([0.03,0.03,0.8,0.94])


        im = axes.imshow(np.rot90(self.odthickmap), cmap=matplotlib.colormaps["gray"])

        cbar = axes.figure.colorbar(im, orientation='vertical',cax=axcb)

        #Show Scale Bar
        startx = int(self.stack.n_rows*0.05)
        starty = self.stack.n_cols-int(self.stack.n_cols*0.05)-self.stack.scale_bar_pixels_y
        um_string = ' $\mathrm{\mu m}$'
        microns = '$'+self.stack.scale_bar_string+' $'+um_string
        axes.text(self.stack.scale_bar_pixels_x+startx+1,starty+1, microns, horizontalalignment='left', verticalalignment='center',
                  color = 'white', fontsize=14)
        #Matplotlib has flipped scales so I'm using rows instead of cols!
        p = matplotlib.patches.Rectangle((startx,starty), self.stack.scale_bar_pixels_x, self.stack.scale_bar_pixels_y,
                               color = 'white', fill = True)
        axes.add_patch(p)


        axes.axis("off")
        self.ODMImagePanel.draw()


#----------------------------------------------------------------------
    def OnSelection1(self, evt):

        x1 = evt.xdata

        self.button_pressed = True
        self.patch = None
        self.conn = self.SpectrumPanel.mpl_connect('motion_notify_event', self.OnSelectionMotion)

        if x1 == None:
            return

        self.x1 = x1

#----------------------------------------------------------------------
    def OnSelection2(self, evt):

        x2 = evt.xdata


        self.button_pressed = False
        self.SpectrumPanel.mpl_disconnect(self.conn)

        if x2 == None:
            return

        x1 = self.x1


        if (self.i0selected == 1) and (self.iselected ==1):
            self.i0selected = 0
            self.iselected = 0

        if self.i0selected == 0:
            self.i0min = np.abs(self.stack.ev - x1).argmin()
            self.i0max = np.abs(self.stack.ev - x2).argmin()

            self.image_i0 = np.sum(self.stack.absdata[:, :, self.i0min:self.i0max+1], axis=2)/(self.i0max+1-self.i0min)

            self.textctrl1.setText('Selection: [ '+str(self.stack.ev[self.i0min]) + ' eV, '+ str(self.stack.ev[self.i0max])+' eV ]' )
            self.i0selected = 1

        elif self.iselected == 0:
            self.imin = np.abs(self.stack.ev - x1).argmin()
            self.imax = np.abs(self.stack.ev - x2).argmin()

            self.image_i = np.sum(self.stack.absdata[:, :, self.imin:self.imax+1], axis=2)/(self.imax+1-self.imin)

            self.textctrl2.setText('Selection: [ '+str(self.stack.ev[self.imin]) + ' eV, '+ str(self.stack.ev[self.imax])+' eV ]' )
            self.iselected = 1

        if (self.i0selected == 1) and (self.iselected ==1):
            nonzeroind = self.image_i0.nonzero()
            self.odthickmap = np.zeros((self.stack.n_cols, self.stack.n_rows))
            self.odthickmap[nonzeroind] = - np.log(self.image_i[nonzeroind]/self.image_i0[nonzeroind])
            self.draw_image()


        self.draw_spectrum()


#----------------------------------------------------------------------
    def OnSelectionMotion(self, event):

        x2 = event.xdata

        if x2 == None:
            return

        x1 = self.x1


        fig = self.specfig

        axes = fig.gca()

        if self.patch != None:
            self.patch.remove()
        self.patch = self.axes.axvspan(x1, x2, facecolor='w', alpha=0.5)

        self.SpectrumPanel.draw()


#----------------------------------------------------------------------
    def OnSave(self, evt):
        #Save images

        wildcard = "Portable Network Graphics (*.png);;Adobe PDF Files (*.pdf);;TIFF File (*.tif);;"

        fileName, _filter = QtWidgets.QFileDialog.getSaveFileName(self, 'Save OD Map', '', wildcard)

        fileName = str(fileName)
        if fileName == '':
            return

        path, ext = os.path.splitext(fileName)
        ext = ext[1:].lower()



        if ext != 'png' and ext != 'pdf' and ext != 'tif':
            error_message = (
                  'Only the PNG and PDF image formats are supported.\n'
                 'A file extension of `png\' or `pdf\' must be used.')
            QtWidgets.QMessageBox.warning(self, 'Error', 'Error - Could not save file.')
            return


        if ext == 'tif':
            from PIL import Image
            img1 = Image.fromarray(self.odthickmap)
            img1.save(fileName)
        else:

            try:
                matplotlib.rcParams['pdf.fonttype'] = 42

                fig = self.odmfig
                fig.savefig(fileName)


            except IOError as e:
                if e.strerror:
                    err = e.strerror
                else:
                    err = e

                QtWidgets.QMessageBox.warning(self, 'Error', 'Could not save file: %s' % err)
#-----------------------------------------------------------------------
