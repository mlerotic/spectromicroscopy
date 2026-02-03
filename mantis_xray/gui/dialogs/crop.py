from PyQt5 import QtWidgets, QtCore, QtGui, uic
from PyQt5.QtCore import pyqtSignal
import pyqtgraph as pg
import numpy as np
import os

class MultiCrop(QtWidgets.QDialog, QtWidgets.QGraphicsScene):
    evlistchanged = pyqtSignal([object])
    thetalistchanged = pyqtSignal([object])
    def __init__(self, parent, common, stack):
        QtWidgets.QWidget.__init__(self, parent)
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'showmulticrop.ui'), self)
        self.parent = parent
        self.stack = stack
        self.com = common
        self.iev = 0
        self.itheta = 0
        self.itheta = 0
        self.sort_mode = 'index'
        self.sort_asc = True
        #self.stack.absdata_shifted_cropped = self.stack.absdata_shifted.copy()

        self.poolthread = QtCore.QThread()
        self.aligned = False
        self.button_ok.setEnabled(True)

        self.setWindowTitle('Stack Cropping')
        self.pglayout = pg.GraphicsLayout(border=None)
        self.canvas.setBackground("w") # canvas is a pg.GraphicsView widget
        self.canvas.setCentralWidget(self.pglayout)
        self.vb = self.pglayout.addViewBox()
        self.vb.setAspectLocked()
        self.i_item = pg.ImageItem(border="k",parent= self)

        self.vb.setMouseEnabled(x=False, y=False)
        self.vb.addItem(self.i_item, ignoreBounds=False)

        self.button_ok.clicked.connect(self.OnAccept)
        self.button_cancel.clicked.connect(self.OnCancel)

        if self.com.stack_loaded == 1:
            self.cb_od_per_px.setVisible(False)
            self.cb_od_per_px.setChecked(True)
            self.cb_od_per_px.stateChanged.connect(self.RedrawPlots)
            self.label_theta_range.setVisible(False)
            self.slider_theta.setVisible(False)
            self.cb_remove_theta.setVisible(False)
            self.groupBox_theta.setVisible(False)
            if self.com.stack_4d == 1:
                self.label_theta_range.setVisible(True)
                self.slider_theta.setVisible(True)
                self.cb_remove_theta.setVisible(True)
                self.groupBox_theta.setVisible(True)
                self.slider_theta.setRange(0, self.stack.n_theta - 1)
                self.slider_theta.valueChanged[int].connect(self.OnScrollTheta)
                self.SetupListTheta()
            # self.maskedvals = [True] * int(self.stack.n_ev)
            # self.spinBoxError.setEnabled(False)
            self.slider_eng.sliderPressed.connect(self.ShowImage)
            self.slider_eng.sliderReleased.connect(self.ShowImage)
            self.slider_eng.valueChanged[int].connect(self.OnScrollEng)
            self.slider_eng.setRange(0, self.stack.n_ev - 1)
            self.pb_selectall.clicked.connect(self.OnSelectAll)
            self.pb_clearall.clicked.connect(self.OnClearAll)
            self.evlistchanged.connect(lambda row: self.qListChangeHandler(row, "energy"))
            self.thetalistchanged.connect(lambda row: self.qListChangeHandler(row, "theta"))
            #self.ev_widget.itemClicked.connect(lambda item: self.OnItemClicked(item))
            self.ev_widget.mousePressEvent = self.mouseEventOnEVList
            self.ev_widget.mouseMoveEvent = self.mouseEventOnEVList
            self.theta_widget.mousePressEvent = self.mouseEventOnThetaList
            self.theta_widget.mouseMoveEvent = self.mouseEventOnThetaList
            #self.ev_widget.itemSelectionChanged.connect(lambda item: self.OnItemClicked(item))
            self.SetupSortCheckbox()
            self.SetupListEV()
            self.OnScrollEng(0)
            self.SetupROI()
            self.SetupPlot()

    def mouseEventOnEVList(self, e):
        if e.type() == QtCore.QEvent.MouseMove or e.type() == QtCore.QEvent.MouseButtonPress:
            qlist = self.ev_widget
            pos = qlist.mapFromGlobal(QtGui.QCursor.pos())
            row = qlist.indexAt(pos).row()
            item = qlist.itemAt(pos)
            #print(row,self.latest_row)
            if row >= 0:
                if e.type() != QtCore.QEvent.MouseMove or row != self.latest_row:
                    if e.buttons() == QtCore.Qt.LeftButton:
                        qlist.setCurrentRow(row)
                        self.evlistchanged.emit(item)
                        self.latest_row = row
        return
    def mouseEventOnThetaList(self, e):
        if e.type() == QtCore.QEvent.MouseMove or e.type() == QtCore.QEvent.MouseButtonPress:
            qlist = self.theta_widget
            pos = qlist.mapFromGlobal(QtGui.QCursor.pos())
            row = qlist.indexAt(pos).row()
            item = qlist.itemAt(pos)
            #print(row,self.latest_row)
            if row >= 0:
                if e.type() != QtCore.QEvent.MouseMove or row != self.latest_row:
                    if e.buttons() == QtCore.Qt.LeftButton:
                        qlist.setCurrentRow(row)
                        self.thetalistchanged.emit(item)
                        self.latest_row = row
        return
    def OnSelectionChanged(self):
        self.RedrawNewPlot()
        #self.UpdateIndices()
        self.region.blockSignals(True)
        if self.idx_selected:
            self.region.setRegion([self.stack.ev[min(self.idx_selected)], self.stack.ev[max(self.idx_selected)]])
            self.region.blockSignals(False)
            self.spectrum_plotwidget.setXRange(*self.region.getRegion())
            self.spectrum_plotwidget.setYRange(np.min(self.pi_new.yData),np.max(self.pi_new.yData))
        return
    def UpdateIndices(self):
        self.idx_selected = sorted([i.data(QtCore.Qt.UserRole) for i in self.ev_selected])
        if self.com.stack_4d:
            self.thetaidx_selected = sorted([i.data(QtCore.Qt.UserRole) for i in self.theta_selected])
    def RedrawPlots(self):
        x,y = self.GenerateSpectrum(list(range(self.stack.n_ev)))
        self.pi.setData(x,y)
        self.OnSelectionChanged()
    def RedrawNewPlot(self):
        self.UpdateIndices()
        x,y = self.GenerateSpectrum(self.idx_selected)
        self.pi_new.setData(x,y)
        if self.idx_selected:
            self.region.show()
    def qListChangeHandler(self,row, dimension):
        if dimension == "theta":
            selection = self.theta_selected
            widget = self.theta_widget
        elif dimension == "energy":
            selection = self.ev_selected
            widget = self.ev_widget

        if row in selection:
            selection.remove(row)
            row.setBackground(QtGui.QColor(0, 0, 0, 0))
        else:
            selection.append(row)
            row.setBackground(QtGui.QColor('#beaed4'))
        if dimension == "theta":
            self.OnScrollTheta(widget.row(row))
        elif dimension == "energy":
            self.OnScrollEng(widget.row(row))
        self.OnSelectionChanged()

    def SetupPlot(self):
        plot = self.spectrum_plotwidget
        plot.setBackground("w")
        plot.setMouseEnabled(x=False, y=False)
        plot.showGrid(y=True)

        plot.showAxis("top", show=True)
        plot.showAxis("right", show=True)
        by = plot.getAxis("right")
        bx = plot.getAxis("top")
        by.setStyle(showValues=False,tickLength=0)
        bx.setStyle(showValues=False,tickLength=0)
        self.ay = plot.getAxis("left")
        ax = plot.getAxis("bottom")
        ax.setLabel(text="Photon energy [eV]")
        x,y = self.GenerateSpectrum(list(range(self.stack.n_ev)))

        self.region = pg.LinearRegionItem(brush=[255,0,0,45],bounds=[np.min(x),np.max(x)])
        plot.addItem(self.region, ignoreBounds=False)
        self.region.setZValue(10)

        self.spectrum_plotwidget.setBackground("w")

        self.pi = plot.plot(x, y, pen=pg.mkPen(color=0.8, width=2))
        self.pi_new = plot.plot(x, y, pen=pg.mkPen(color="b", width=2))
        self.refmarker = pg.InfiniteLine(angle=90, movable=False,
                                         pen=pg.mkPen(color="b", width=2, style=QtCore.Qt.DashLine))
        plot.addItem(self.refmarker, ignoreBounds=True)
        self.region.setRegion((min(x),max(x)))
        self.region.sigRegionChangeFinished.connect(lambda region: self.UpdateEVRegion(region))
    # ----------------------------------------------------------------------
    def getDataClosestToRegion(self,region,plotitem,snapregion=False):
        minidx, maxidx = region.getRegion()
        data = plotitem.getData()[0]
        index = lambda x: np.argmin(np.abs(data - x))
        minidx = index(minidx)
        maxidx = index(maxidx)
        if minidx == maxidx:
            minidx = 0
            maxidx = np.argmax(data)
        mindata = data[minidx]
        maxdata = data[maxidx]
        if snapregion:
            region.blockSignals(True)
            region.setRegion([mindata, maxdata])  # snap region to data points
            region.blockSignals(False)
        return minidx, maxidx, mindata, maxdata
    
    def UpdateEVRegion(self, region):
        min_idx, max_idx,*_  = self.getDataClosestToRegion(region,self.pi, True)
        y_vals = self.pi.yData[min_idx:max_idx + 1]
        x_vals = self.pi.xData[min_idx:max_idx + 1]
        self.spectrum_plotwidget.setRange(yRange=[np.min(y_vals), np.max(y_vals)], xRange=[np.min(x_vals), np.max(x_vals)], disableAutoRange=True, padding=0.05)
        if min_idx not in self.idx_selected:
            for idx in range(int(min_idx),np.min(self.idx_selected)):
                self.ev_selected.append(self.ev_widget.item(idx))
                self.ev_widget.item(idx).setBackground(QtGui.QColor('#beaed4'))
        else:
            for idx in range(np.min(self.idx_selected),int(min_idx)):
                try:
                    self.ev_selected.remove(self.ev_widget.item(idx))
                except ValueError:
                    pass
                self.ev_widget.item(idx).setBackground(QtGui.QColor(0, 0, 0, 0))
        if max_idx not in self.idx_selected:
            for idx in range(np.max(self.idx_selected)+1,int(max_idx)+1):
                try:
                    self.ev_selected.append(self.ev_widget.item(idx))
                except ValueError:
                    pass
                self.ev_widget.item(idx).setBackground(QtGui.QColor('#beaed4'))
        else:
            for idx in range(int(max_idx)+1,np.max(self.idx_selected)+1):
                try:
                    self.ev_selected.remove(self.ev_widget.item(idx))
                except ValueError:
                    pass
                self.ev_widget.item(idx).setBackground(QtGui.QColor(0, 0, 0, 0))
        self.RedrawNewPlot()
    def OnSelectAll(self):
        self.ev_widget.clear()
        self.SetupListEV()
        if self.com.stack_4d:
            self.theta_widget.clear()
            self.SetupListTheta()
        self.region.setRegion((min(self.stack.ev), max(self.stack.ev)))
        self.region.show()
        self.RedrawNewPlot()
    def OnClearAll(self):
        #self.ev_widget.clear()
        #self.SetupListEV()
        self.region.hide()
        for idx in self.idx_selected:
            self.ev_widget.item(idx).setBackground(QtGui.QColor(0, 0, 0, 0))
        self.ev_selected = []
        self.idx_selected = []
        if self.com.stack_4d:
            for idx in self.thetaidx_selected:
                self.theta_widget.item(idx).setBackground(QtGui.QColor(0, 0, 0, 0))
            self.theta_selected = []
            self.thetaidx_selected = []
        self.RedrawNewPlot()
    def SetupListTheta(self):
        self.theta_selected = []
        self.thetaidx_selected = []
        for i,e in enumerate(self.stack.theta): # Fill QList with energies
            #self.stk.shifts.append([1,0,(0.0,0.0)]) #checked [0,1]; pre, post, undefined state for map [-1,1,0],(xshift [float],yshift [float])
            item = QtWidgets.QListWidgetItem(str(int(i)).zfill(4)+"     at     " + format(e, '.1f') + "°")
            item.setData(QtCore.Qt.UserRole, i)
            self.theta_widget.addItem(item)
            self.theta_selected.append(item)
            self.thetaidx_selected.append(i)
            item.setBackground(QtGui.QColor('#beaed4'))
            item.setForeground(QtGui.QColor(0, 0, 0, 128))
    def SetupListEV(self):
        self.ev_widget.setSortingEnabled(False)
        # Preserve selection state if possible
        preserved = set()
        select_all = True
        if hasattr(self, 'idx_selected') and self.idx_selected:
             preserved = set(self.idx_selected)
             select_all = False
             
        self.ev_widget.clear()
        self.ev_selected = []
        self.idx_selected = []
        
        data = list(enumerate(self.stack.ev))
        
        if hasattr(self, 'cb_sort_asc') and self.cb_sort_asc.isChecked():
            data.sort(key=lambda x: x[1])
               
        for i,e in data:
            item = QtWidgets.QListWidgetItem(str(int(i)).zfill(4)+"     at     " + format(e, '.2f') + " eV")
            item.setData(QtCore.Qt.UserRole, i)
            self.ev_widget.addItem(item)
            
            if select_all or i in preserved:
                self.ev_selected.append(item)
                # self.idx_selected.append(i) # Will be updated via UpdateIndices or manually here?
                # Original code updated both manually. But UpdateIndices is reliable.
                item.setBackground(QtGui.QColor('#beaed4'))
                item.setForeground(QtGui.QColor(0, 0, 0, 128))
            else:
                 item.setBackground(QtGui.QColor(0, 0, 0, 0))
        
        self.UpdateIndices()

    def SetupSortCheckbox(self):
        layout = self.ev_widget.parentWidget().layout()
        if layout:
             container = QtWidgets.QWidget()
             hl = QtWidgets.QHBoxLayout(container)
             hl.setContentsMargins(0,0,0,0)
             
             self.cb_sort_asc = QtWidgets.QCheckBox("Sort energies")
             self.cb_sort_asc.setChecked(True)
             self.cb_sort_asc.stateChanged.connect(self.OnSortChange)
             
             hl.addWidget(self.cb_sort_asc)
             hl.addStretch()
             
             idx = layout.indexOf(self.ev_widget)
             if idx >= 0:

                 if isinstance(layout, QtWidgets.QBoxLayout):
                    layout.insertWidget(idx, container)
                 elif isinstance(layout, QtWidgets.QGridLayout):
                     try:
                        index = layout.indexOf(self.ev_widget)
                        row, col, rowSpan, colSpan = layout.getItemPosition(index)
                        layout.removeWidget(self.ev_widget)
                        layout.addWidget(container, row, col, 1, colSpan)
                        layout.addWidget(self.ev_widget, row + 1, col, rowSpan, colSpan)
                     except:
                        pass
    
    def OnSortChange(self):
        self.SetupListEV()

    def OnScrollTheta(self, value):
        self.slider_theta.setValue(value)
        self.ResetAllItems(self.theta_widget)
        
        # Map visual row (value) to data index
        item = self.theta_widget.item(value)
        if item:
            self.itheta = item.data(QtCore.Qt.UserRole)
        else:
            self.itheta = value

        #self.stack.absdata = self.stack.stack4D[:,:,:,self.itheta].copy()
        self.ShowImage()
        self.theta_widget.setCurrentRow(value) # Set visual row
        if self.com.stack_loaded == 1:
            if item:
                 item.setForeground(QtGui.QColor(0, 0, 0, 255))
        self.RedrawPlots()
    def OnScrollEng(self, value):
        # value is visual row index from slider or list
        self.slider_eng.setValue(value)
        self.ResetAllItems(self.ev_widget)
        
        # Map visual row to data index
        item = self.ev_widget.item(value)
        if item:
            self.iev = item.data(QtCore.Qt.UserRole)
        else:
            self.iev = value
            
        self.ShowImage()
        try:
            self.refmarker.setValue(self.stack.ev[self.iev])
        except:
            pass
        self.ev_widget.setCurrentRow(value) # Set visual selection
        if self.com.stack_loaded == 1:
            if item:
                item.setForeground(QtGui.QColor(0, 0, 0, 255))
    def ShowImage(self):
        if self.com.stack_4d == 1:
            self.i_item.setImage(self.stack.stack4D[:, :, int(self.iev),int(self.itheta)])
            self.groupBox.setTitle(str('Stack Browser | Image at {0:5.2f} eV and {1:5.1f}°').format(float(self.stack.ev[self.iev]),float(self.stack.theta[self.itheta]), ))
        else:
            if self.com.i0_loaded:
                self.i_item.setImage(self.stack.od3d[:, :, int(self.iev)])
            else:
                self.i_item.setImage(self.stack.absdata[:, :, int(self.iev)])
            self.groupBox.setTitle(str('Stack Browser | Image at {0:5.2f} eV').format(float(self.stack.ev[self.iev])))

    ## Setup a ROI for an alignment rectangle. By default the whole image area is used.
    def SetupROI(self):
        #self.stack.absdata = self.stack.absdata.copy()
        self.box = pg.RectROI(self.i_item.boundingRect().topLeft(), self.i_item.boundingRect().bottomRight(),
                              pen=(5, 8), handlePen=QtGui.QPen(QtGui.QColor(255, 0, 128, 255)), centered=False,
                              sideScalers=False, removable=False, scaleSnap=True, translateSnap=True,
                              maxBounds=self.i_item.boundingRect())
        self.vb.addItem(self.box, ignoreBounds=False)
        self.box.sigRegionChangeFinished.connect(self.RedrawPlots)
        self.box.sigRegionChangeStarted.connect(self.OnBoxChanging)
        self.button_rstroi.clicked.connect(self.OnResetROI)
    ## The ROI is limited to the visible image area. OnMouseMoveOutside handles the behavior when
    def OnResetROI(self):
        self.box.setPos(0, 0, update=False, finish=False)
        self.box.setSize(self.i_item.boundingRect().bottomRight() - self.box.pos(), update=True, snap=True, finish=True)
        self.box.show()
    def OnBoxChanging(self):
        self.boxsize = self.box.size()
        self.proxy = pg.SignalProxy(self.vb.scene().sigMouseMoved, rateLimit=30, slot=self.OnMouseMoveOutside)
    def GetRegion(self):
        try:
            self.proxy.disconnect()
        except AttributeError:
            pass
        left = int(self.box.pos().x())
        right = left + int(self.box.size().x())
        bottom = int(self.box.pos().y())
        top = bottom + int(self.box.size().y())
        return (left,right,top,bottom)

    def GenerateSpectrum(self, evselection):
        left,right,top,bottom = self.GetRegion()
        if self.com.i0_loaded:
            self.cb_od_per_px.setVisible(True)
            if self.cb_od_per_px.isChecked():
                self.ay.setLabel(text="Optical density per px")
            else:
                self.ay.setLabel(text="Sum of optical densities in ROI")
            if self.com.stack_4d:
                total = self.stack.od4d[left:right, bottom:top, :, int(self.itheta)].copy()
            else:
                total = self.stack.od3d[left:right, bottom:top, :].copy()
        else:
            self.ay.setLabel(text="Photon flux per px [cps]")
            if self.com.stack_4d == 1:
                t = [self.stack.theta[i] for i in self.thetaidx_selected]
                self.label_theta_range.setText(
                    "Theta range: [ " + str(min(t, default=0)) + "° .. " + str(
                        max(t, default=0)) + "° ], # values: " + str(
                        len(t)))
                total = self.stack.stack4D[left:right, bottom:top, :, int(self.itheta)].copy()
            else:
                total = self.stack.absdata[left:right, bottom:top, :].copy()
        if self.cb_od_per_px.isChecked():
            total = total.sum(axis=(0,1)) / (int(self.box.size().x()) * int(self.box.size().y()))
        else:
            total = total.sum(axis=(0,1))
        x = [self.stack.ev[i] for i in evselection]
        y = [total[i] for i in evselection]
        self.label_spatial_range.setText("Stack size: [ "+str(int(self.box.size().x()))+" x "+str(int(self.box.size().y()))+" ] px²")
        self.label_ev_range.setText(
            "Energy range: [ " + str(min(x, default=0)) + " .. " + str(max(x, default=0)) + " ] eV, # values: "+ str(len(x)))
        return (x, y)

    def ResetAllItems(self,widget):
        for i in range(widget.count()):
            widget.item(i).setForeground(QtGui.QColor(0, 0, 0, 128))

    def OnMouseMoveOutside(self, ev):
        mousepos = self.vb.mapSceneToView(ev[0])
        if not self.vb.itemBoundingRect(self.i_item).contains(mousepos):
            maxrect = self.i_item.boundingRect().bottomRight()
            # if bounds exceeded
            out_x = max((mousepos-maxrect).x(),0)
            out_y = max((mousepos-maxrect).y(), 0)
            if self.box.size() != self.boxsize: # prevents taking action when the box is just dragged and not resized
                if out_x and out_y:
                    self.box.setSize(self.i_item.boundingRect().bottomRight()-self.box.pos(), update=True, snap=True, finish=False)
                elif out_x:
                    self.box.setSize([self.i_item.boundingRect().right()-self.box.pos().x(),mousepos.y()-self.box.pos().y()], update=True, snap=True, finish=False)
                elif out_y:
                    self.box.setSize([mousepos.x()-self.box.pos().x(),self.i_item.boundingRect().bottom()-self.box.pos().y()], update=True, snap=True, finish=False)
    # ----------------------------------------------------------------------
    def OnCancel(self, evt):
        self.close()
    # ----------------------------------------------------------------------
    def OnAccept(self, evt):
        if self.cb_croptoroi.isChecked():
            left, right, top, bottom = self.GetRegion()
        else:
            left, right, top, bottom = (None,None,None,None)

        if self.cb_remove_evs.isChecked():
            selection = []
            for row in range(self.ev_widget.count()):
                item = self.ev_widget.item(row)
                if item in self.ev_selected:
                    selection.append(item.data(QtCore.Qt.UserRole))
            
            if len(selection) == 0:
                QtWidgets.QMessageBox.warning(self, 'Error', 'Please select at least one energy value!')
                return
            self.stack.n_ev = np.array(len(selection))
            self.stack.ev = self.stack.ev[selection]
            self.stack.data_dwell = self.stack.data_dwell[selection]
        else:
            selection = list(range(self.stack.n_ev))


        self.stack.absdata = self.stack.absdata[left:right,bottom:top,selection]
        self.stack.n_cols = self.stack.absdata.shape[0]
        self.stack.n_rows = self.stack.absdata.shape[1]
        self.parent.page1.ix = int(self.stack.n_cols/2)
        self.parent.page1.iy = int(self.stack.n_rows/2)
        self.stack.scale_bar()

        if self.com.stack_4d:
            if self.cb_remove_theta.isChecked():
                thetas = []
                for row in range(self.theta_widget.count()):
                    item = self.theta_widget.item(row)
                    if item in self.theta_selected:
                        thetas.append(item.data(QtCore.Qt.UserRole))

                if len(thetas) == 0:
                    QtWidgets.QMessageBox.warning(self, 'Error', 'Please select at least one theta value!')
                    return
                self.stack.n_theta = len(thetas)
                self.stack.theta = self.stack.theta[thetas]

            else:
                thetas = list(range(self.stack.n_theta))
            self.stack.stack4D = self.stack.stack4D[left:right, bottom:top, selection, :]
            self.stack.stack4D = self.stack.stack4D[:,:, :, thetas]
        if self.com.i0_loaded:
            self.stack.i0_mask = self.stack.i0_mask[left:right,bottom:top]
            if self.com.stack_4d:
                self.stack.od4d = self.stack.od4d[left:right,bottom:top,selection, thetas]
            else:
                self.stack.od3d =  self.stack.od3d[left:right,bottom:top,selection]
                self.stack.od = self.stack.od3d.copy()
                self.stack.od = np.reshape(self.stack.od, (self.stack.n_rows * self.stack.n_cols, self.stack.n_ev),
                                       order='F')
            self.parent.page1.specfig.I0Update()

        self.stack.fill_h5_struct_from_stk()
        if self.com.i0_loaded == 1:
            self.stack.fill_h5_struct_normalization()

        # Fix the slider on Page 1!
        if self.com.stack_4d:
            self.parent.page1.slider_theta.setRange(0, self.stack.n_theta - 1)
            self.parent.page1.itheta = 0
            self.parent.page1.slider_theta.blockSignals(True)
            self.parent.page1.slider_theta.setValue(int(self.parent.page1.itheta))
            self.parent.page1.slider_theta.blockSignals(False)

            self.parent.page0.slider_theta.setRange(0, self.stack.n_theta - 1)
            self.parent.page0.itheta = 0
            self.parent.page0.slider_theta.blockSignals(True)
            self.parent.page0.slider_theta.setValue(int(self.parent.page1.itheta))
            self.parent.page0.slider_theta.blockSignals(False)

        #self.parent.page1.slider_eng.setRange(0, self.stack.n_ev - 1)
        #self.parent.page1.iev = 0
        #self.parent.page1.slider_eng.setValue(int(self.parent.page1.iev))

        #self.parent.page0.slider_eng.setRange(0, self.stack.n_ev - 1)
        #self.parent.page0.iev = 0
        #self.parent.page0.slider_eng.setValue(int(self.parent.page1.iev))

#        self.parent.page1.loadSpectrum(self.parent.page1.ix, self.parent.page1.iy)
        #self.parent.page1.Clear()
        self.parent.page1.absimgfig.loadNewImageWithROI()
        self.parent.page0.absimgfig.loadNewImage()
        self.parent.page1.specfig.ClearandReload()
        self.parent.window().refresh_widgets()
        #if showmaptab:
        #    self.parent.page9.Clear()
        #    self.parent.page9.loadData()

        self.close()
