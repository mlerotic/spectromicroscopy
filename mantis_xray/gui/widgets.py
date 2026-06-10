
import pyqtgraph as pg
from PyQt5 import QtCore, QtGui, QtWidgets
import numpy as np
import os
import re
from scipy.interpolate import interp1d
from scipy import ndimage
from ..helpers import PDFExporter

MANTIS_HEX = [
    '#6699DD', '#EE7733', '#ABCC44', '#99DDFF', '#FFAABB',
    '#BAAA00', '#AB2622', '#44BB99', '#AA4499', '#EEDD89'
]


def build_mantis_lut(num_colors):
    n = max(int(num_colors), 1)
    base = np.array([QtGui.QColor(c).getRgb()[:3] for c in MANTIS_HEX], dtype=np.uint8)
    idx = np.arange(n) % len(base)
    return base[idx]

# Note: These classes heavily depend on 'parent' which seems to be a page or something with specific attributes.
# In a proper refactor, we should define interfaces or pass specific data/signals instead of 'parent'.
# For now, we move them as-is to separate file, but they still are tightly coupled.

class SpecFig():
    def __init__(self,parent, plotwidget):
        self.parent = parent
        self.plot = plotwidget
        self.pi = self.plot.getPlotItem()
        self.plot.setBackground("w")
        self.pi.layout.setSpacing(12)
        self.pi.layout.setContentsMargins(10,10,40,10)
        self.pi.showGrid(y=True)
        self.pi.showAxis("top", show=True)
        self.pi.showAxis("right", show=True)
        by = self.pi.getAxis("right")
        bx = self.pi.getAxis("top")
        by.setStyle(showValues=False,tickLength=0)
        bx.setStyle(showValues=False,tickLength=0)
        self.ay = self.pi.getAxis("left")
        self.ay.setLabel(text="counts")

        self.ax = self.pi.getAxis("bottom")
        self.ax.setLabel(text="Photon energy [eV]")
        self.pi.setTitle("")
        self.pi.addLegend()
        self.LineIndicator = pg.InfiniteLine(angle=90, movable=True, markers=None,
                                              pen=pg.mkPen(color=QtGui.QColor(0, 0, 0, 128), width=1.5, style=QtCore.Qt.DashLine))
        self.LineIndicatorLabel = pg.InfLineLabel(self.LineIndicator, " ")
        try:
            self.parent.button_lockspectrum.clicked.connect(self.OnLockSpectrum)
            self.parent.button_clearspecfig.clicked.connect(self.ClearandReload)
            self.parent.button_clearlastroi.clicked.connect(self.ClearLast)
            self.parent.button_mergeroi.clicked.connect(self.mergeROI)
            self.parent.button_subtractroi.clicked.connect(self.subtractROI)
        except AttributeError: #e.g., PagePCACluster does not have roi manipulation buttons.
            pass

        # canvas.scene().sigMouseClicked.connect(lambda _: self.setFocus())
        self.plot.scene().sigMouseClicked.connect(lambda event: self.setActive())

    def setActive(self):
        self.parent.window().active_widget = self  # Returns the active widget for CTRL +C functionality

    def ClearLast(self):
        #self.plot.blockSignals(True)
        i = 0
        if self.parent.ROIShapeBox.currentText() == "Histogram":
            i = 2
        roiitems = [image for image in self.parent.absimgfig.imageplot.items if isinstance(image, pg.ImageItem)]
        try:
            if roiitems[-1] != self.parent.absimgfig.ROImask and len(self.parent.absimgfig.imageplot.items) > 2:
                self.pi.removeItem(self.pi.items[-1 - i])
                self.parent.absimgfig.imageplot.removeItem(roiitems[-1])
        except IndexError:  # if previously removed, ignore
            pass
        #self.plot.blockSignals(False)

    def ClearandReload(self, fig=None):
        self.roicolor = (0, 0, 255, 255)
        self.plot.blockSignals(True)
        self.pi.setMouseEnabled(x=True, y=True)
        # Qt clicked(bool) may pass a bool positional arg; ignore it here.
        if fig is None or isinstance(fig, bool) or not hasattr(fig, "imageplot"):
            fig = self.parent.absimgfig
        iterator = len(self.pi.items)-1
        # The expression "for item in self.pi.items: " does not work! Instead we count the items and iterate through them
        for i in range(iterator,-1,-1):
                self.pi.removeItem(self.pi.items[i])  # remove spectra
                try:    # remove locked ROIs from imageplot
                    roiitem = fig.imageplot.items[i+1]# +1 because transparent selection "ROImask" exists.
                    if isinstance(roiitem, pg.ImageItem) and len(fig.imageplot.items) > 2:
                        fig.imageplot.removeItem(roiitem)
                except IndexError:  # if previously removed, ignore
                    pass
        try:
            self.LineIndicator.sigPositionChangeFinished.disconnect()
            self.LineIndicator.sigPositionChanged.disconnect()
            self.plot.sigRangeChanged.disconnect()
            self.plot.scene().sigMouseClicked.disconnect()
            self.parent.absimgfig.roi.sigRegionChanged.disconnect()
        except:
            pass
        self.plot.blockSignals(False)
        self.setCaption()
        self.loadNewSpectrum()
        vb = self.pi.items[1].getViewBox()
        vb.enableAutoRange()
        if hasattr(self.parent, "ROIShapeBox"):
            if self.parent.ROIShapeBox.currentText() != "Lasso":
                self.parent.absimgfig.OnROIVisibility(self.parent.ROIvisibleCheckBox.checkState())
        # Clear stale lock fallback data after explicit clear-all.
        self.last_locked_roi_spectrum = None


    def ClearforHistogram(self):
        self.plot.blockSignals(True)
        self.pi.setMouseEnabled(x=False, y=False)
        iterator = len(self.pi.items)-1
        # The expression "for item in self.pi.items: " does not work! Instead we count the items and iterate through them
        for i in range(iterator,-1,-1):
                self.pi.removeItem(self.pi.items[i])  # remove spectra
                try:    # remove locked ROIs from imageplot
                    roiitem = self.parent.absimgfig.imageplot.items[i+1]# +1 because transparent selection "ROImask" exists.
                    if isinstance(roiitem, pg.ImageItem) and len(self.parent.absimgfig.imageplot.items) > 3:
                        self.parent.absimgfig.imageplot.removeItem(roiitem)
                except IndexError:  # if previously removed, ignore
                    pass
        try:
            self.LineIndicator.sigPositionChangeFinished.disconnect()
            self.LineIndicator.sigPositionChanged.disconnect()
            self.plot.sigRangeChanged.disconnect()
            self.plot.scene().sigMouseClicked.disconnect()
            self.parent.absimgfig.roi.sigRegionChanged.disconnect()
        except:
            pass
        self.plot.blockSignals(False)

    def setPlotItemVisibility(self, show):
        iterator = len(self.pi.items) - 1
        if show:
            for i in range(iterator, -1, -1):
                self.pi.items[i].show()
        else:
            for i in range(iterator, -1, -1):
                self.pi.items[i].hide()


    def toggleI0Spectrum(self):
        if self.parent.button_showi0.isChecked():
            self.formatAxesLabels(type="showi0")
            self.parent.ROIvisibleCheckBox.setEnabled(False)
            self.parent.ROIShapeBox.setEnabled(False)
            self.parent.button_lockspectrum.setEnabled(False)
            self.parent.button_clearspecfig.setEnabled(False)
            self.parent.button_clearlastroi.setEnabled(False)
            self.parent.button_mergeroi.setEnabled(False)
            self.parent.button_subtractroi.setEnabled(False)

            self.setPlotItemVisibility(False)
            x, y = (self.parent.stk.evi0, self.parent.stk.i0data)
            curve = pg.PlotCurveItem(x,y, pen=({'color': "#ff7700", 'width': 2}),
                                     skipFiniteCheck=True, name="I0")
            self.pi.addItem(curve)
            #Show I0 region and hide roi selection
            indices = np.where(self.parent.stk.i0_mask)
            self.drawROImask(None,indices= indices,color=(255,119,0,255))
            self.parent.absimgfig.OnROIVisibility(QtCore.Qt.Unchecked)
            self.parent.absimgfig.ROImask.show()

        else:
            self.setPlotItemVisibility(True)
            if len(self.pi.items) > 2:
                self.pi.removeItem( self.pi.items[-1])
            self.formatAxesLabels()
            self.parent.ROIvisibleCheckBox.setEnabled(True)
            self.parent.ROIShapeBox.setEnabled(True)
            self.parent.button_lockspectrum.setEnabled(True)
            self.parent.button_clearspecfig.setEnabled(True)
            self.parent.button_clearlastroi.setEnabled(True)
            self.parent.button_mergeroi.setEnabled(True)
            self.parent.button_subtractroi.setEnabled(True)

            # Restore ROI
            self.parent.absimgfig.OnROIVisibility(self.parent.ROIvisibleCheckBox.checkState())
            self.updatePlotData()


    def OnI0Histogram(self):
        self.parent.OnShowMean()
        self.parent.button_meanflux.setChecked(True)
        self.parent.ROIvisibleCheckBox.setEnabled(False)
        self.parent.button_lockspectrum.setEnabled(False)
        self.parent.button_clearspecfig.setEnabled(False)
        self.parent.button_clearlastroi.setEnabled(False)
        self.parent.button_mergeroi.setEnabled(False)
        self.parent.button_subtractroi.setEnabled(False)
        self.roicolor = (255,119,0,255)
        self.parent.absimgfig.OnROIShapeChanged("Histogram")
        #self.parent.ROIShapeBox.setCurrentText("Histogram")
        self.parent.ROIShapeBox.setStyleSheet("color: #ff7700;");
        self.parent.label_roitype.setText("I0 type")
        self.parent.label_roitype.setStyleSheet("color: #ff7700;");
        self.parent.button_i0.disconnect()
        self.parent.button_i0.setText("Accept I0")
        self.parent.button_i0.setStyleSheet("color: #ff7700;");
        self.parent.button_i0.clicked.connect( self.OnI0Accept)

    # ----------------------------------------------------------------------
    def OnI0Accept(self, evt):
        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
        self.parent.OnShowMean()
        self.parent.button_meanflux.setChecked(False)
        self.parent.ROIvisibleCheckBox.setEnabled(True)
        self.parent.button_lockspectrum.setEnabled(True)
        self.parent.button_clearspecfig.setEnabled(True)
        self.parent.button_clearlastroi.setEnabled(True)
        self.parent.button_mergeroi.setEnabled(True)
        self.parent.button_subtractroi.setEnabled(True)
        self.parent.stk.i0_mask = np.sum(self.parent.absimgfig.ROIrgba, axis=2)[:, :] > 0
        self.I0Update()

    # ----------------------------------------------------------------------
    def OnI0Reset(self):

        self.parent.stk.reset_i0()
        self.parent.com.i0_loaded = 0
        self.roicolor = (0,0,255,255)
        self.parent.showflux = True
        # self.rb_flux.setChecked(True)

        self.parent.absimgfig.loadNewImageWithROI()
        self.ClearandReload()
        self.parent.window().refresh_widgets()

    # ----------------------------------------------------------------------
    def I0Update(self):
        bool = np.where(self.parent.stk.i0_mask == True)
        if bool[0].size and bool[1].size:
            result = self.parent.stk.i0_from_histogram(bool)
            if result is False: # Aborted due to duplicates
                QtWidgets.QApplication.restoreOverrideCursor()
                return

            self.parent.I0histogramCalculated()
            self.parent.button_i0.disconnect()
            self.parent.button_i0.setText("Reset I0")
            self.parent.button_i0.setStyleSheet("");  # pass an empty string to return to default style
            self.parent.ROIShapeBox.setStyleSheet("");
            self.parent.button_i0.clicked.connect(self.OnI0Reset)
            self.parent.label_roitype.setText("ROI type")
            self.parent.label_roitype.setStyleSheet("");
            QtWidgets.QApplication.restoreOverrideCursor()
        else:
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self.parent, 'Error', 'I0 region is empty!')

    # ----------------------------------------------------------------------
    def I0Cancel(self):
        self.parent.button_i0.disconnect()
        self.roicolor = (0,0,255,255)
        self.parent.showflux = True
        self.parent.button_i0.setStyleSheet("");  # pass an empty string to return to default style
        self.parent.ROIShapeBox.setStyleSheet("");
        self.parent.label_roitype.setText("ROI type")
        self.parent.label_roitype.setStyleSheet("");
        self.parent.button_i0.setText("Select I0")
        self.parent.button_i0.clicked.connect(self.OnI0Histogram)
        self.parent.button_i0ffile.setEnabled(True)
        self.parent.button_prenorm.setEnabled(True)
        self.parent.button_refimgs.setEnabled(True)

    # ----------------------------------------------------------------------
    def GetNextROINumberandColor(self):
        #MANTiS unique Light qualitative color scheme
        lut = ['#6699DD','#EE7733','#ABCC44','#99DDFF','#FFAABB','#BAAA00','#AB2622','#44BB99','#AA4499','#EEDD89']
        hues = len(lut)
        index = len(self.pi.items)
        if self.parent.ROIShapeBox.currentText() == "Histogram":
            index = index - 2

        name = "ROI " + str(index - 1)
        color = QtGui.QColor(lut[(index - 2) % hues])
        return name, color

    def OnLockSpectrum(self):
        name, color = self.GetNextROINumberandColor()
        color.setAlpha(200)
        curve = pg.PlotCurveItem(pen=({'color': color, 'width': 2}), skipFiniteCheck=True,
                                 name=name)
        curve.hide()
        x = self.pi.items[1].xData
        y = self.pi.items[1].yData
        self.pi.addItem(curve)
        if self.parent.ROIShapeBox.currentText() == "Histogram":
            self.pi.items.remove(curve)
            self.pi.items.insert(-2, curve)
        else:
            curve.show()
        self.parent.absimgfig.addLockedROI(color)
        curve.setData(x,y)
        self.last_locked_roi_spectrum = (np.asarray(x).copy(), np.asarray(y).copy())

    def get_roi_spectra_for_dose(self):
        spectra = []

        def add_if_valid(label, y_data):
            if y_data is None:
                return
            arr = np.asarray(y_data, dtype=float)
            if arr.ndim == 1 and arr.size == self.parent.stk.n_ev:
                spectra.append((label, arr.copy()))

        # Current interactive ROI spectrum (main curve).
        current_y = None
        if len(self.pi.items) > 1 and isinstance(self.pi.items[1], pg.PlotCurveItem):
            current_y = getattr(self.pi.items[1], "yData", None)
        if current_y is None:
            # Recompute from the active ROI if plot curve data are temporarily unavailable.
            data = self.prefilterData()
            mask = self.createROImask()
            _x, current_y = self.getSpecfromROI(data, mask)
        add_if_valid("Current ROI", current_y)

        # Locked ROI spectra
        for item in self.pi.items:
            if not isinstance(item, pg.PlotCurveItem):
                continue
            name = None
            if hasattr(item, "opts") and isinstance(item.opts, dict):
                name = item.opts.get("name")
            if isinstance(name, str) and name.startswith("ROI "):
                add_if_valid(name, getattr(item, "yData", None))

        unique = {}
        for label, spectrum in spectra:
            unique[label] = spectrum

        ordered = []
        if "Current ROI" in unique:
            ordered.append(("Current ROI", unique.pop("Current ROI")))

        def roi_sort_key(label):
            match = re.match(r"ROI\s+(\d+)", label)
            return int(match.group(1)) if match else 10**9

        for label in sorted(unique.keys(), key=roi_sort_key):
            ordered.append((label, unique[label]))

        return ordered

    def removeLast2ROI(self,i,roiitems):
        self.pi.removeItem(self.pi.items[-1 - i])
        self.pi.removeItem(self.pi.items[-1 - i])
        self.parent.absimgfig.imageplot.removeItem(roiitems[-1])
        self.parent.absimgfig.imageplot.removeItem(roiitems[-2])

    def mergeROI(self):
        i = 0
        if self.parent.ROIShapeBox.currentText() == "Histogram":
            i = 2
        self.plot.blockSignals(True)
        roiitems = [image for image in self.parent.absimgfig.imageplot.items if isinstance(image, pg.ImageItem)]
        try:
            if isinstance(roiitems[-2], pg.ImageItem) and len(roiitems) > 3:
                b1 = np.sum(roiitems[-1].image, axis = 2)[:,:] > 0
                b2 = np.sum(roiitems[-2].image, axis = 2)[:,:] > 0
                boolmask = ~np.logical_or(b1, b2)
                indices = np.where(boolmask == False)
                self.removeLast2ROI(i,roiitems)
                roi = np.zeros([*boolmask.shape, 4], dtype=np.uint8)
                name, color = self.GetNextROINumberandColor()
                color.setAlpha(200)
                roi[indices] = color.getRgb()
                lockedroi = pg.ImageItem(image=roi, border="k", opacity=0.5)
                self.parent.absimgfig.imageplot.addItem(lockedroi, ignoreBounds=True)

                data = self.prefilterData()
                x, y = self.getSpecfromROI(data,boolmask)
                curve = pg.PlotCurveItem(pen=({'color': color, 'width': 2}),
                                         skipFiniteCheck=True, name=name)
                curve.hide()
                self.pi.addItem(curve)
                if self.parent.ROIShapeBox.currentText() == "Histogram":
                    self.pi.items.remove(curve)
                    self.pi.items.insert(-2, curve)
                else:
                    curve.show()
                curve.setData(x, y)
        except IndexError:
            pass
        self.plot.blockSignals(False)

    def subtractROI(self):
        i = 0
        if self.parent.ROIShapeBox.currentText() == "Histogram":
            i = 2
        self.plot.blockSignals(True)
        roiitems = [image for image in self.parent.absimgfig.imageplot.items if isinstance(image, pg.ImageItem)]
        try:
            if roiitems[-2] != self.parent.absimgfig.ROImask and len(self.parent.absimgfig.imageplot.items) > 2:
                b1 = np.sum(roiitems[-1].image, axis = 2)[:,:] > 0
                b2 = np.sum(roiitems[-2].image, axis = 2)[:,:] > 0
                #boolmask = ~np.logical_and(b1, b2) # intersect
                boolmask = ~np.logical_and(~b1, b2)  # subtract
                if not np.any(~boolmask):
                    self.ClearLast()
                    self.ClearLast() #!
                    return
                indices = np.where(boolmask == False)
                self.removeLast2ROI(i,roiitems)
                roi = np.zeros([*boolmask.shape, 4], dtype=np.uint8)
                name, color = self.GetNextROINumberandColor()
                color.setAlpha(200)
                roi[indices] = color.getRgb()
                lockedroi = pg.ImageItem(image=roi, border="k", opacity=0.5)
                self.parent.absimgfig.imageplot.addItem(lockedroi, ignoreBounds=True)

                data = self.prefilterData()
                x, y = self.getSpecfromROI(data,boolmask)
                curve = pg.PlotCurveItem(pen=({'color': color, 'width': 2}),
                                         skipFiniteCheck=True, name=name)
                curve.hide()
                self.pi.addItem(curve)
                if self.parent.ROIShapeBox.currentText() == "Histogram":
                    self.pi.items.remove(curve)
                    self.pi.items.insert(-2, curve)
                else:
                    curve.show()
                curve.setData(x, y)
            elif roiitems[-2] == self.parent.absimgfig.ROImask:
                self.ClearLast()
        except IndexError:
            pass
        self.plot.blockSignals(False)

    def loadNewSpectrum(self):
        if self.pi.items:
            try:
                self.ClearandReload()
            except IndexError:
                pass
        curve = pg.PlotCurveItem(pen=({'color': pg.intColor(6), 'width': 2}), skipFiniteCheck=True)

        self.pi.addItem(self.LineIndicator, ignoreBounds=True)
        self.pi.addItem(curve)
        try:
            self.LineIndicator.sigPositionChangeFinished.disconnect()
            self.LineIndicator.sigPositionChanged.disconnect()
            self.plot.sigRangeChanged.disconnect()
            self.plot.scene().sigMouseClicked.disconnect()
            self.parent.absimgfig.roi.sigRegionChanged.disconnect()
        except:
            pass
        if hasattr(self.parent, "ROIvisibleCheckBox"):
            self.parent.ROIvisibleCheckBox.setEnabled(True)
        if hasattr(self.parent, "ROIShapeBox"):
            self.parent.absimgfig.OnROIShapeChanged(self.parent.ROIShapeBox.currentText())
            if self.parent.ROIShapeBox.currentText() == "Lasso":
                self.pi.items[1].hide()
                self.parent.absimgfig.OnROIVisibility(self.parent.ROIvisibleCheckBox.checkState())
        if hasattr(self.parent, "selpca"):
            self.updatePlotDataOnPCA()
        self.LineIndicator.addMarker("o")
        self.dot = self.LineIndicator.markers[0][0]
        self.LineIndicator.setZValue(10)
        self.ypos = self.getIntersectionY()

        self.LineIndicator.sigPositionChanged.connect(self.OnUpdateLineIndicator)
        self.LineIndicator.sigPositionChangeFinished.connect(self.SnapIndicatorToEV)
        self.plot.sigRangeChanged.connect(self.OnUpdateLineIndicator)
        self.plot.scene().sigMouseClicked.connect(self.OnMouseClick)
        self.OnUpdateLineIndicator()

    def updatePlotData(self):
        data = self.prefilterData()
        mask = self.createROImask()
        self.drawROImask(mask, color = self.roicolor)
        x, y = self.getSpecfromROI(data,self.parent.absimgfig.boolmask)
        if len(self.pi.items) < 2 or not isinstance(self.pi.items[1], pg.PlotCurveItem):
            return
        self.plot.blockSignals(True)
        self.pi.items[1].setData(x, y)
        if self.pi.items[1].xData is not None and self.pi.items[1].yData is not None and len(self.pi.items[1].xData) > self.parent.iev:
            self.LineIndicator.setPos(QtCore.QPointF(self.pi.items[1].xData[self.parent.iev], self.pi.items[1].yData[self.parent.iev]))
        self.plot.blockSignals(False)

    def updatePlotDataOnPCA(self):
        if hasattr(self, "roicolor"):
            color = self.roicolor
            self.plot.blockSignals(True)
            #self.formatAxesLabels(type=type)
            try:
                self.pi.removeItem(self.pi.items[1])
            except:
                pass
            #self.region.sigRegionChanged.disconnect()
            self.setPlotItemVisibility(True)
            #data = self.prefilterData()
            #mask = self.createROImask()
            #self.drawROImask(mask, color=color)
            #x, y = self.getSpecfromROI(data,self.parent.absimgfig.boolmask)
            #print(self.parent.anlz.eigenvecs[:,0])
            curve = pg.PlotCurveItem(self.parent.stk.ev , self.parent.anlz.eigenvecs[:,self.parent.selpca], pen=({'color': color, 'width': 2}), skipFiniteCheck=True, name='Component {}'.format(self.parent.selpca + 1))
            self.pi.setMouseEnabled(x=True, y=True)
            self.pi.addItem(curve)
            self.pi.items.remove(curve)
            self.pi.items.insert(1, curve)
            self.plot.blockSignals(False)
            if hasattr(self, "dot"):
                self.OnUpdateLineIndicator()

    def updatePlotDataOnROIShapeChange(self,type):
        color = self.roicolor
        self.plot.blockSignals(True)
        self.formatAxesLabels(type=type)
        self.pi.removeItem(self.pi.items[1])
        #self.region.sigRegionChanged.disconnect()
        try:
            self.pi.removeItem(self.region)
            self.pi.removeItem(self.histogram)
        except:
            pass

        if type == "Histogram":
            self.setPlotItemVisibility(False)
            self.region = pg.LinearRegionItem(brush=color, bounds=[np.min(self.parent.stk.hist_data_x), np.max(self.parent.stk.hist_data_x)])
            #self.region.setZValue(10)
            self.region.setOpacity(0.3)
            self.pi.addItem(self.region, ignoreBounds=False)
            curve = pg.PlotCurveItem(pen=({'color': color, 'width': 2}), skipFiniteCheck=True)
            self.histogram = pg.PlotCurveItem(self.parent.stk.hist_data_x, self.parent.stk.hist_data_y, pen=({'color': color, 'width': 1}), skipFiniteCheck=True, name="histogram",
                                     stepMode=True, fillLevel=0, brush=color)
            self.pi.setMouseEnabled(x=False, y=False)
            self.pi.addItem(self.histogram)
            curve.hide()
        else:
            self.setPlotItemVisibility(True)
            data = self.prefilterData()
            mask = self.createROImask()
            self.drawROImask(mask, color=color)
            x, y = self.getSpecfromROI(data,self.parent.absimgfig.boolmask)
            curve = pg.PlotCurveItem(x, y, pen=({'color': color, 'width': 2}), skipFiniteCheck=True)
            self.pi.setMouseEnabled(x=True, y=True)
        self.pi.addItem(curve)
        self.pi.items.remove(curve)
        self.pi.items.insert(1, curve)
        self.plot.blockSignals(False)
        if hasattr(self, "dot"):
            self.OnUpdateLineIndicator()

        def update(region):
            self.region.setZValue(10)
            minX, maxX = region
            i0_indices = np.where((minX <= self.parent.stk.histogram) & (self.parent.stk.histogram <= maxX))
            data = self.prefilterData()
            self.drawROImask(None,indices = i0_indices, color=color)
            x, y = self.getSpecfromROI(data,self.parent.absimgfig.boolmask)
            self.pi.items[1].setData(x, y)
        if type == "Histogram":
            self.region.sigRegionChanged.connect(lambda: update(self.region.getRegion()))
            self.region.setRegion((self.parent.stk.histmin, self.parent.stk.histmax))

    def setCaption(self):
        fn = ""
        if self.parent.com.fntocaption and hasattr(self.parent, "tc_file"):
            fn = self.parent.tc_file.text()
        self.pi.setTitle("<center>{}</center>".format(fn))

    def formatAxesLabels(self,type=None):
        if type == "Histogram":
            #self.pi.getViewBox().setLimits()#yMin=0, yMax=np.max(y))
            self.ax.setLabel(text="Average Flux")
            self.ay.setLabel(text="log<sub>10</sub> (Number of pixels)")
            return
        elif type == "showi0":
            self.ay.setLabel(text="Flux in selected I0 area [counts]")
            return
        self.ax.setLabel(text="Photon energy [eV]")
        if self.parent.com.i0_loaded:
            self.ay.setLabel(text="Optical density per px inside ROI")
        else:
            self.ay.setLabel(text="Flux per px inside ROI [counts]")

    def getIntersectionY(self):
        x_newgrid = self.pi.items[1].xData
        y_newgrid = self.pi.items[1].yData

        if len(self.pi.items[1].xData)> 2:
            func = interp1d(x_newgrid, self.pi.items[1].yData)
            x_newgrid = np.linspace(min(x_newgrid), max(x_newgrid), num=min(5000,(40*len(x_newgrid))), endpoint=True)
            y_newgrid = func(x_newgrid)
        diff = np.abs(x_newgrid - self.LineIndicator.value())
        idx = np.argmin(diff)

        vb = self.pi.items[1].getViewBox()
        ypos = 1 + (y_newgrid[idx] - vb.viewRect().bottom()) / (
                    vb.viewRect().bottom() - vb.viewRect().top())
        return ypos

    def SnapIndicatorToEV(self):
        diff = np.abs(self.pi.items[1].xData - self.LineIndicator.value())
        idx = np.argmin(diff)
        self.parent.slider_eng.blockSignals(True)
        self.parent.OnScrollEng(idx)
        self.parent.slider_eng.blockSignals(False)

    def OnUpdateLineIndicator(self):
        if hasattr(self.parent, "absimgfig") and self.parent.absimgfig.currentroishape == "Histogram":
            return
        self.ypos = self.getIntersectionY()
        self.LineIndicator.markers = [(self.dot, self.ypos, 10)]
        self.LineIndicator.update()
        self.LineIndicatorLabel.setPosition(self.ypos)
        self.LineIndicatorLabel.setFormat(" "+str(round(self.LineIndicator.value(),2)) + " eV ")

    def OnMouseClick(self,e):
        if e.double():
            vb = self.pi.items[1].getViewBox()
            pos = vb.mapSceneToView(e.scenePos()).x()
            self.LineIndicator.blockSignals(True)
            self.LineIndicator.setPos(pos)
            self.LineIndicator.blockSignals(False)
            self.LineIndicator.sigPositionChangeFinished.emit(self)

    def GetRegion(self,box):
        left = int(box.pos().x())
        right = left + int(box.size().x())
        bottom = int(box.pos().y())
        top = bottom + int(box.size().y())
        left = max(0,left)
        bottom = max(0,bottom)
        top = min(self.parent.stk.n_rows,max(0,top))
        right = min(self.parent.stk.n_cols,max(0,right))
        return (left,right,top,bottom)

    def createROImask(self):
        if isinstance(self.parent.absimgfig.roi, type(None)):
            return None
        cols = self.parent.stk.n_cols
        rows = self.parent.stk.n_rows
        angle = self.parent.absimgfig.roi.angle()
        offsetx = 0
        offsety = 0
        boolmask = np.full((cols , rows), True)
        if self.parent.absimgfig.roi.boundingRect().width() > 1000:
            print("ROI wider than 1000 px is not supported. Refer to _getArrayRegionForArbitraryShape().")
            return None
        mask, coords = self.parent.absimgfig.roi.getArrayRegion(boolmask,self.parent.absimgfig.imageitem, axes=(0,1), returnMappedCoords=True)
        if mask.size == 0:
            return None
        # calculate offsets to avoid mismatch of roi and selected pixels
        if angle:
            offsetx = np.clip(np.sin(0.785398) * np.cos(np.radians(angle) + 0.785398) - 0.5, -1, 0)
            offsety = np.clip(np.sin(0.785398) * np.cos(np.radians(angle) - 0.785398) - 0.5, -1, 0)
        x = np.rint(coords[0] + offsetx).astype(int).flatten()
        y = np.rint(coords[1] + offsety).astype(int).flatten()
        # Delete indices outside image region
        if isinstance(self.parent.absimgfig.roi, pg.RectROI):
            idcs = np.where((x >= cols) | (x < 0) | (y >= rows) | (y < 0))
        else:
            mask = mask.astype(bool).flatten()
            idcs = np.where((x > cols) | (x < 0) | (y > rows) | (y < 0) | ~mask)
        x = np.delete(x, idcs)
        y = np.delete(y, idcs)
        boolmask[x, y] = False
        # Handle edge cases and gaps. Fill holes in data due to rounding errors by applying and reverting a binary dilation.
        boolmask = np.pad(boolmask, ((1, 1), (1, 1)), 'constant', constant_values=True)
        boolmask = ndimage.binary_dilation(~ndimage.binary_dilation(~boolmask))[1:-1, 1:-1]
        return boolmask

    def getSpecfromROI(self, data, boolmask):
        cols = self.parent.stk.n_cols
        rows = self.parent.stk.n_rows
        ev = self.parent.stk.n_ev
        if boolmask is None:
            return [self.parent.stk.ev[i] for i in list(range(ev))], np.zeros(ev)
        valid_pixel_count = (cols * rows) - np.count_nonzero(boolmask)
        if not valid_pixel_count:
            return [self.parent.stk.ev[i] for i in list(range(ev))], np.zeros(ev)
        mask = np.broadcast_to(np.expand_dims(boolmask, axis=2), (cols, rows, ev))
        spectrum = np.ma.array(data, mask=mask).sum(axis=(0,1)) / max(valid_pixel_count,1)
        x = [self.parent.stk.ev[i] for i in list(range(ev))]
        y = [spectrum[i] for i in list(range(ev))]
        return x, y

    def drawROImask(self, boolmask, indices=None, color = (0,0,255,255)):
        roirgba = self.parent.absimgfig.ROIrgba
        roirgba[:, :] = [0, 0, 0, 0]
        roimask = self.parent.absimgfig.ROImask
        if not indices:
            cols = self.parent.stk.n_cols
            rows = self.parent.stk.n_rows
            valid_pixel_count = (cols * rows) - np.count_nonzero(boolmask)
            if not valid_pixel_count or boolmask is None:
                self.parent.absimgfig.boolmask[:,:] = True
            else:
                self.parent.absimgfig.boolmask = boolmask
            indices = ~self.parent.absimgfig.boolmask
        else:
            self.parent.absimgfig.boolmask[:, :] = True
            self.parent.absimgfig.boolmask[indices] = False
        roirgba[indices] = color
        roimask.setImage(roirgba)

    def prefilterData(self):
        if self.parent.com.i0_loaded:
            #self.cb_od_per_px.setVisible(True)
            #if self.cb_od_per_px.isChecked():
            if self.parent.com.stack_4d:
                data = self.parent.stk.od4d[:, :, :, int(self.parent.itheta)].copy()
            else:
                data = self.parent.stk.od3d
        else:
            if self.parent.com.stack_4d == 1:
                #t = [self.parent.stk.theta[i] for i in self.parent.itheta]
                #self.label_theta_range.setText(
                #    "Theta range: [ " + str(min(t, default=0)) + "° .. " + str(
                #        max(t, default=0)) + "° ], # values: " + str(
                #        len(t)))
                data = self.parent.stk.stack4D[:, :, :, int(self.parent.itheta)]
            else:
                data = self.parent.stk.absdata
        # self.label_spatial_range.setText("Stack size: [ "+str(int(self.box.size().x()))+" x "+str(int(self.box.size().y()))+" ] px²")
        # self.label_ev_range.setText(
        #     "Energy range: [ " + str(min(x, default=0)) + " .. " + str(max(x, default=0)) + " ] eV, # values: "+ str(len(x)))
        return data

    def OnCopy(self):
        # self.exp = pg.exporters.ImageExporter(self.pi) # just plot
        exp = pg.exporters.ImageExporter(self.plot.scene()) # plot and axes, i.e., complete viewbox
        exp.export(copy=True)
        return

    def SaveFig(self,fileName):
        fileName = str(fileName)
        if fileName == '':
            return
        path, ext = os.path.splitext(fileName)
        ext = ext[1:].lower()

        if ext == 'svg':
            exp = pg.exporters.SVGExporter(self.plot.scene())
        elif ext== 'pdf':
            exp = PDFExporter(self.plot.scene())
        else:
            exp = pg.exporters.ImageExporter(self.plot.scene())
        if ext in ['tif','png','jpg','svg','pdf']:
            exp.export(fileName)

class ImgFig():
    def __init__(self,parent,canvas):
        self.scale = 0.000001
        self.parent = parent
        canvas.setBackground("w") # canvas is a pg.GraphicsLayoutWidget
        self.imageplot = canvas.addPlot()
        self.imageplot.setMouseEnabled(x=False, y=False)
        self.imageitem = pg.ImageItem(border="k")
        self.scalebar = pg.ScaleBar(size=1, suffix="m")
        self.imageplot.setAspectLocked(lock=True, ratio=1)
        self.imageplot.showAxis("top", show=True)
        self.imageplot.showAxis("bottom", show=True)
        self.imageplot.showAxis("left", show=True)
        self.imageplot.showAxis("right", show=True)
        self.ay1 = self.imageplot.getAxis("left")
        by1 = self.imageplot.getAxis("right")
        self.ax1 = self.imageplot.getAxis("bottom")
        bx1 = self.imageplot.getAxis("top")
        self.ay1.setLabel(text="y",units="px")
        self.ay1.enableAutoSIPrefix(enable=True)
        self.ax1.setLabel(text="x",units="px")
        self.ax1.enableAutoSIPrefix(enable=True)
        self.ay1.setStyle(tickLength=8)
        self.ax1.setStyle(tickLength=8)
        by1.setStyle(showValues=False,tickLength=0)
        bx1.setStyle(showValues=False,tickLength=0)
        self.imageplot.setTitle("no data")
        self.map = "gray"
        cm = pg.colormap.get(self.map, source="matplotlib")
        self.bar = pg.ColorBarItem(values=(0, 1), colorMap=cm, rounding=0.0001)  # init color bar
        self.mousepressed = False
        #canvas.scene().sigMouseClicked.connect(lambda _: self.setFocus())
        canvas.scene().sigMouseClicked.connect(lambda event: self.setActive())
    
    # ... (Rest of ImgFig methods would go here, need to extract carefully)
    # Since I cannot see all ImgFig methods in the previous view, I will just put placeholders and methods I saw.
    # Actually I should read more of the file to get full ImgFig.

    def setActive(self):
        self.parent.window().active_widget = self  # Returns the active widget for CTRL +C functionality

    def loadPCAImage(self):
        self.clear()
        self.parent.selpca = 0
        self.loadPCAData()

    def loadNewImage(self):
        self.clear()
        self.parent.iev = 0
        self.loadData()

    def loadNewImageWithROI(self):
        self.parent.stk.calc_histogram()
        self.clear()
        self.parent.iev = 0
        self.loadData()
        self.currentroishape = "Lasso"
        self.addROI((0,0),(self.imageitem.boundingRect().width(), self.imageitem.boundingRect().height()), self.currentroishape)
        self.parent.ROIShapeBox.blockSignals(True)
        self.parent.ROIShapeBox.setCurrentText(self.currentroishape)
        self.parent.ROIShapeBox.blockSignals(False)

    def clear(self):
        self.imageplot.removeItem(self.scalebar)
        self.imageplot.removeItem(self.imageitem)
        self.imageplot.clear()

    def OnROIVisibility(self, state):
        if state == QtCore.Qt.Checked:
            if not isinstance(self.roi, type(None)):
                self.roi.show()
                self.parent.specfig.pi.items[1].show()  # curve
                self.parent.specfig.pi.items[0].show()  # line & dot
            try:
                self.parent.specfig.histogram.show()
                self.parent.specfig.region.show()
            except:
                pass
            self.ROImask.show()
            self.parent.ROIShapeBox.setEnabled(True)
            self.parent.button_lockspectrum.setEnabled(True)
            self.parent.button_clearspecfig.setEnabled(True)
            self.parent.button_clearlastroi.setEnabled(True)
            self.parent.button_mergeroi.setEnabled(True)
            self.parent.button_subtractroi.setEnabled(True)
        else:
            if not isinstance(self.roi, type(None)):
                self.roi.hide()
                self.parent.specfig.pi.items[1].hide()
                self.parent.specfig.pi.items[0].hide()
            try:
                self.parent.specfig.histogram.hide()
                self.parent.specfig.region.hide()
            except:
                pass
            self.ROImask.hide()
            vb = self.parent.specfig.pi.getViewBox()
            vb.updateAutoRange()
            self.parent.ROIShapeBox.setEnabled(False)
            self.parent.button_lockspectrum.setEnabled(False)

    def OnROIShapeChanged(self, shape):
        self.parent.ROIShapeBox.blockSignals(True)
        self.parent.ROIShapeBox.setCurrentText(shape)
        self.parent.ROIShapeBox.blockSignals(False)
        #self.parent.specfig.plot.blockSignals(True)
        if isinstance(self.roi, pg.PolyLineROI):
            pos = np.rint(self.roi.pos())
            state = self.roi.getState()
            points = state['points']
            try:
                minx= min(points, key=lambda item: item[0])[0]
                miny = min(points, key=lambda item: item[1])[1]
                x = (max(points, key=lambda item: item[0])[0] -
                    minx)
                y = (max(points, key=lambda item: item[1])[1] -
                     miny)
                pos = QtCore.QPointF(np.rint(pos[0]+minx),np.rint(pos[1]+miny))
            except ValueError: # of no points found expand to image dimension.
                x,y = (self.imageitem.boundingRect().width(), self.imageitem.boundingRect().height())
                pos = QtCore.QPointF(0,0)
                pass
            size = np.rint((x,y))
        elif isinstance(self.roi, type(None)):
            pos = QtCore.QPointF(0, 0)
            size = np.rint((self.imageitem.boundingRect().width(), self.imageitem.boundingRect().height()))
        else:
            pos = np.rint(self.roi.pos())
            size = np.rint(self.roi.size())
        try:
            self.roi.disconnect()
        except:
            pass
        self.imageplot.removeItem(self.ROImask)
        self.imageplot.removeItem(self.roi)
        self.currentroishape = shape
        self.addROI(pos, size, shape)
        self.parent.specfig.updatePlotDataOnROIShapeChange(shape)

    def onMousePress(self,e):
        if not self.parent.button_showi0.isChecked() and self.parent.ROIvisibleCheckBox.isChecked():
            self.mousepressed = True
    def onMouseRelease(self,e):
        if self.mousepressed:
            self.mousepressed = False
            if len(self.roi.handles) > 2:
                self.proxy.disconnect()
                self.roi.addSegment(self.roi.handles[-1]['item'], self.roi.handles[0]['item'])
                self.imageitem.mousePressEvent = self.imageplot.mousePressEvent
                self.imageitem.mouseReleaseEvent = self.imageplot.mouseReleaseEvent
    def onMouseMoved(self,e):
        pos = self.vb.mapSceneToView(e[0])
        roipos = pos-self.vb.mapFromViewToItem(self.roi,pos)
        if self.mousepressed and self.vb.itemBoundingRect(self.imageitem).contains(pos):
            pos = (np.rint(pos.x()-roipos.x()), np.rint(pos.y()-roipos.y()))
            try:
                if self.roi.handles[-1]['pos'] != QtCore.QPointF(*pos): # if same point as before, do not add handle
                    self.roi.addFreeHandle(pos)
                    self.roi.addSegment(self.roi.handles[-2]['item'], self.roi.handles[-1]['item'])
            except IndexError: # first handle
                self.parent.specfig.pi.items[1].show()
                self.parent.specfig.pi.items[0].show()
                self.roi.addFreeHandle(pos)

    def addROI(self,pos, size, shape):
        self.imageitem.mousePressEvent = self.imageplot.mousePressEvent
        self.imageitem.mouseReleaseEvent = self.imageplot.mouseReleaseEvent
        kwargs= {'pen': (5, 8), 'handlePen' : QtGui.QPen(QtGui.QColor(255, 0, 128, 255)), 'resizable' : True, 'removable' : False, 'movable' : True, 'scaleSnap' : True, 'translateSnap' : True}
        selection = {"Rectangle": pg.RectROI(pos,size,**kwargs), "Circle": pg.CircleROI(pos,size,**kwargs),
                     "Ellipse": pg.EllipseROI(pos,size,**kwargs), "Polygon": pg.PolyLineROI(positions= [(0,0),(0,size[1]),(size[0],size[1]),(size[0],0)], pos=pos, closed=True,**kwargs),
                     "Lasso": pg.PolyLineROI(positions= [], pos=pos, closed=True,**kwargs), "Histogram":None}
        self.roi = selection[shape]
        if not shape == "Histogram":
            self.imageplot.addItem(self.roi, ignoreBounds=True)
        if shape == "Lasso":
            try:
                self.parent.specfig.pi.items[1].hide()
                self.parent.specfig.pi.items[0].hide()
            except IndexError: # if first roi, spectra are not existing at this point
              pass
            self.roi.handlePen = QtGui.QPen(QtGui.QColor(0, 0, 0, 0)) # Make handles invisible
            self.imageitem.mousePressEvent = self.onMousePress
            self.imageitem.mouseReleaseEvent = self.onMouseRelease
            self.proxy = pg.SignalProxy(self.vb.scene().sigMouseMoved, rateLimit=15, slot=self.onMouseMoved)
        else:
            try:
                self.parent.specfig.pi.items[1].show()
                self.parent.specfig.pi.items[0].show()
            except IndexError: # if first roi, spectra are not existing at this point
              pass
        self.ROImask = pg.ImageItem(border="k", opacity=0.5)
        self.imageplot.addItem(self.ROImask)
        self.ROIrgba = np.zeros([*self.imageitem.image.shape, 4], dtype=np.uint8)
        self.boolmask = np.full((self.parent.stk.n_cols, self.parent.stk.n_rows), True)
        # items are appended to the items list. relocate the two freshly added items, otherwise roi removal fails.
        self.imageplot.items.remove(self.ROImask)
        self.imageplot.items.insert(1, self.ROImask)
        if not shape == "Histogram":
            self.imageplot.items.remove(self.roi)
            self.imageplot.items.insert(1, self.roi)
            self.roi.setZValue(10)  # make sure ROI is drawn above image
            self.roi.sigRegionChanged.connect(self.parent.specfig.updatePlotData)

    def addLockedROI(self,color):
        roi = np.zeros([*self.imageitem.image.shape, 4], dtype=np.uint8)
        indices = np.where(self.boolmask == False)
        roi[indices] = color.getRgb()
        lockedroi = pg.ImageItem(image=roi, border="k", opacity=0.5)

        self.imageplot.addItem(lockedroi, ignoreBounds=True)
        lockedroi.setZValue(11)  # make sure ROI is drawn above image
        # Locked ROIs are already valid selection state; enable ROI actions now.
        for button_name in ("button_saveROIspectr", "button_setROII0", "button_ROIdosecalc"):
            button = getattr(self.parent, button_name, None)
            if button is not None:
                button.setEnabled(True)
        if self.parent.ROIShapeBox.currentText() == "Lasso": # remove lasso ROI after function call
            self.OnROIShapeChanged("Lasso")

    def loadPCAData(self): # Called when fresh data are loaded.
        try:
            self.vb.sigRangeChanged.disconnect()
        except:
            pass
        self.imageplot.addItem(self.imageitem)
        self.vb = self.imageitem.getViewBox()
        rightlabel = self.bar.getAxis("right")
        #if self.parent.com.i0_loaded:
        #    rightlabel.setLabel(text="OD", units="")
        #else:
        #    rightlabel.setLabel(text="counts", units="")
        #self.parent.slider_eng.blockSignals(True)
        self.parent.slider_cl.setRange(0, min(19, self.parent.stk.n_ev - 1))
        self.bar.setImageItem(self.imageitem, insert_in=self.imageplot)
        self.OnColormapChange(map=self.parent.CMMapBox.currentText(),num_colors=self.parent.StepSpin.value(),fliplut=self.parent.ColorFlipCheckBox.isChecked())
        self.parent.OnScrollCl(self.parent.selpca) # Plot image & set Scrollbar
        try:
            self.OnMetricScale(self.parent.MetricCheckBox.isChecked(), True, False)
        except AttributeError:
            self.OnMetricScale(False, True, self.parent.SquarePxCheckBox.isChecked())
        self.OnShowScale()
        self.vb.sigRangeChanged.connect(lambda: self.OnUpdateScale(self.parent.ScalebarCheckBox.isChecked()))


    def loadData(self): # Called when fresh data are loaded.
        try:
            self.vb.sigRangeChanged.disconnect()
        except:
            pass
        self.imageplot.addItem(self.imageitem)
        self.vb = self.imageitem.getViewBox()
        rightlabel = self.bar.getAxis("right")
        if self.parent.com.i0_loaded:
            rightlabel.setLabel(text="OD", units="")
        else:
            rightlabel.setLabel(text="counts", units="")
        #self.parent.slider_eng.blockSignals(True)
        self.parent.slider_eng.setRange(0, self.parent.stk.n_ev - 1)

        self.bar.setImageItem(self.imageitem, insert_in=self.imageplot)
        self.OnColormapChange(map=self.parent.CMMapBox.currentText(),num_colors=self.parent.StepSpin.value(),fliplut=self.parent.ColorFlipCheckBox.isChecked())
        self.parent.OnScrollEng(self.parent.iev) # Plot image & set Scrollbar
        try:
            self.OnMetricScale(self.parent.MetricCheckBox.isChecked(), True, False)
        except AttributeError:
            self.OnMetricScale(False, True, self.parent.SquarePxCheckBox.isChecked())
        self.OnShowScale()
        self.vb.sigRangeChanged.connect(lambda: self.OnUpdateScale(self.parent.ScalebarCheckBox.isChecked()))

    def draw(self,image,setlabel=True,setpcalabel=None,setlut=False,levels=False):
        if setlut:
            self.OnColormapChange(map=self.parent.CMMapBox.currentText(),num_colors=self.parent.StepSpin.value(),fliplut=self.parent.ColorFlipCheckBox.isChecked())
        self.imageitem.setImage(image)
        if setlabel:
            self.setCaption(setpcalabel)
        if levels:
            min,max = levels
        else:
            min = np.nanmin(image)  # ignoring nans
            max = np.nanmax(image)
        if not np.isnan(min) and not np.isnan(max):
            self.bar.setLevels(low=min, high=max)

    def setCaption(self,setpcalabel=None,custom_title=False):
        fn = ""
        if self.parent.com.fntocaption and hasattr(self.parent, "tc_file"):
            fn = self.parent.tc_file.text() + " - "
        mean = getattr(self.parent, "mean_visible", None)
        if mean:
            if self.parent.com.i0_loaded == 0:
                str = "Mean flux"
            elif self.parent.com.i0_loaded == 1:
                str = "Mean OD"
            self.imageplot.setTitle("<center>{}{}</center>".format(fn, str))
            return
        elif setpcalabel is None:
            if self.parent.com.stack_4d == 1:
                self.imageplot.setTitle("<center>{0}Image at {1:5.2f} eV and {2:5.1f}°</center>".format(fn,float(self.parent.stk.ev[self.parent.iev]),
                                                                                              float(self.parent.stk.theta[
                                                                                                        self.parent.itheta])))
            else:
                self.imageplot.setTitle("<center>{0}Image at {1:5.2f} eV</center>".format(fn,float(self.parent.stk.ev[self.parent.iev])))

        elif setpcalabel is True:
            # In cluster display mode the left panel always shows the cluster
            # composite (or error map), not a PCA component – use the correct title.
            if getattr(self.parent, 'cluster_display_mode', False):
                if getattr(self.parent, 'showerrormap', 0):
                    self.imageplot.setTitle("<center>{}Cluster Error Map</center>".format(fn))
                else:
                    if self.parent.com.stack_4d == 1 and hasattr(self.parent.stk, 'theta') and len(self.parent.stk.theta) > self.parent.itheta:
                        self.imageplot.setTitle("<center>{}Cluster composite at {:5.1f}\u00b0</center>".format(
                            fn, float(self.parent.stk.theta[self.parent.itheta])))
                    else:
                        self.imageplot.setTitle("<center>{}Cluster composite</center>".format(fn))
            else:
                self.imageplot.setTitle("<center>{}Component {:02d}</center>".format(fn, int(self.parent.selpca) + 1))
        #if custom_title:
        #    self.imageplot.setTitle("<center>{}{}</center>".format(fn, custom_title))


    def OnMetricScale(self, setmetric= True, zeroorigin= True, square= False):
        if self.parent.com.stack_loaded == 1:
            if setmetric==True:
                self.parent.SquarePxCheckBox.setVisible(False)
                self.parent.ZeroOriginCheckBox.setVisible(True)
                self.imageplot.setAspectLocked(lock=True, ratio=1)
                #self.p2.setAspectLocked(lock=True, ratio=1)
                if not zeroorigin:
                    x_start = self.parent.stk.x_start*self.scale
                    y_start = self.parent.stk.y_start*self.scale
                else:
                    x_start = 0
                    y_start = 0
                self.ay1.setLabel(text="y", units="m")
                self.ax1.setLabel(text="x", units="m")
                #self.ay2.setLabel(text="y", units="m")
                #self.ax2.setLabel(text="x", units="m")

                self.imageitem.setRect(QtCore.QRectF(x_start, y_start, self.scale*self.parent.stk.n_cols*self.parent.stk.x_pxsize, self.scale*self.parent.stk.n_rows*self.parent.stk.y_pxsize))
                #if hasattr(self, "OD"):
                #    self.m_item.setRect(QtCore.QRectF(x_start, y_start, self.scale*np.shape(self.OD)[0]*self.stk.x_pxsize, self.scale*np.shape(self.OD)[1]*self.stk.y_pxsize))
                #    self.setCrosshair()
            else:
                try:
                    self.parent.ZeroOriginCheckBox.setVisible(False)
                except AttributeError:
                    pass
                self.parent.SquarePxCheckBox.setVisible(True)
                if square == True:
                    aspect = 1
                    self.parent.ScalebarCheckBox.setVisible(False)
                else:
                    self.parent.ScalebarCheckBox.setVisible(True)
                    aspect = self.parent.stk.x_pxsize/self.parent.stk.y_pxsize
                    #print(aspect)
                self.imageplot.setAspectLocked(lock=True, ratio=aspect)
                self.ay1.setLabel(text="y", units="px")
                self.ax1.setLabel(text="x", units="px")
                self.imageitem.setRect(QtCore.QRectF(0, 0, self.parent.stk.n_cols, self.parent.stk.n_rows))

    def OnCatChanged(self):
        self.parent.CMMapBox.blockSignals(True)
        self.parent.CMMapBox.clear()
        self.parent.CMMapBox.blockSignals(False)
        self.parent.CMMapBox.addItems(self.parent.cmaps[self.parent.CMCatBox.currentIndex()][1])

    def OnColormapChange(self, map="gray", num_colors=256, fliplut=False):
        self.map = map
        if self.map == "Mantis":
            lut = build_mantis_lut(num_colors)
            if fliplut:
                lut = np.ascontiguousarray(lut[::-1])
        else:
            luttup = (0,1)
            if fliplut:
                luttup = (1,0)
            cm = pg.colormap.get(self.map, source="matplotlib")
            lut = cm.getLookupTable(*luttup, num_colors)
        if self.parent.com.stack_loaded == 1:
            try:
                lut = np.ascontiguousarray(lut)
                self.imageitem.setLookupTable(lut)
                lut = np.expand_dims(lut, axis=1)
                qimg = pg.functions.ndarray_to_qimage(lut, QtGui.QImage.Format.Format_RGB888)
                self.bar.bar.setPixmap(QtGui.QPixmap.fromImage(qimg).scaled(1, 256))
            except AttributeError:
                self.bar.bar.setLookupTable(lut)
                self.imageitem.setLookupTable(lut)

    def OnShowScale(self):
        suffix = "m"
        self.scalebar.text.setText(pg.siFormat(float(self.parent.stk.scale_bar_string)*self.scale, suffix=suffix))
        self.scalebar.setParentItem(self.imageplot.getViewBox())
        self.scalebar.anchor((1, 1), (1, 1), offset=(-20, -20))
        self.scalebar.hide()
        self.OnUpdateScale(self.parent.ScalebarCheckBox.isChecked())

    def OnUpdateScale(self, set):
        if hasattr(self.parent.stk, "scale_bar_string"):
            if not set or self.parent.SquarePxCheckBox.isChecked():
                self.scalebar.hide()
                return
            if not hasattr(self.parent, "MetricCheckBox"):
                self.scalebar.size = float(self.parent.stk.scale_bar_string) / self.parent.stk.x_pxsize
            elif not self.parent.MetricCheckBox.isChecked():
                self.scalebar.size = float(self.parent.stk.scale_bar_string) / self.parent.stk.x_pxsize
            elif self.parent.MetricCheckBox.isChecked():
                self.scalebar.size = float(self.parent.stk.scale_bar_string)*self.scale
            self.scalebar.updateBar()
            self.scalebar.show()

    def OnCopy(self):
        # self.exp = pg.exporters.ImageExporter(self.imageitem) # just image
        exp = pg.exporters.ImageExporter(self.imageplot)  # image and axes, i.e., complete viewbox
        exp.export(copy=True)
        return

    def SaveFig(self,fileName):
        fileName = str(fileName)
        if fileName == '':
            return
        path, ext = os.path.splitext(fileName)
        ext = ext[1:].lower()

        if ext == 'svg':

            exp = pg.exporters.SVGExporter(self.imageplot)
            # The SVG output is clean.
            # Display errors (line widths, etc.) likely result from external software not properly handling SVG.

        elif ext== 'pdf':
            exp = PDFExporter(self.imageplot)
        else:
            exp = pg.exporters.ImageExporter(self.imageplot)
        if ext in ['tif','png','jpg','svg','pdf']:
            exp.export(fileName)
