import os
import numpy as np
from PyQt5 import QtWidgets, QtCore, QtGui, uic
from PyQt5.QtCore import pyqtSignal
import matplotlib.pyplot as plt
from scipy import ndimage
import pyqtgraph as pg

from ..widgets import ImgFig
from .spectral_dialogs import SpectralImageMap

class ShowODMap(QtWidgets.QWidget):
    qlistchanged = pyqtSignal([tuple])
    def __init__(self, parent, common, data_struct, stack):
        super(ShowODMap, self).__init__()
        dir_path = os.path.dirname(os.path.realpath(__file__))
        uic.loadUi(os.path.join(dir_path,'showodmap.ui'), self)
        self.show()
        self.cmaps = [('Perceptually Uniform Sequential', [
            'viridis', 'plasma', 'inferno', 'magma']),
                      ('Sequential', [
                          'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                          'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                          'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']),
                                            ('Sequential (2)', [
                          'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
                          'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
                          'hot', 'afmhot', 'gist_heat', 'copper']),
                      ('Diverging', [
                          'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
                          'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']),
                      ('Qualitative', [
                          'Pastel1', 'Pastel2', 'Paired', 'Accent',
                          'Dark2', 'Set1', 'Set2', 'Set3',
                          'tab10', 'tab20', 'tab20b', 'tab20c']),
                      ('Miscellaneous', [
                          'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
                          'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv',
                          'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar'])]
        self.initUI(parent, common, data_struct, stack)
        #    self.Clear()
        self.odimgfig.loadNewImage()
        self.loadData()

#-----------------------------------------------------------------------
    def initUI(self, parent, common, data_struct, stack):
        self.scale = 0.000001
        self.data_struct = data_struct
        self.stk = stack
        self.com = common
        self.parent = parent
        self.iev = 0
        self.latest_row = -1
        self.xoffset = 0
        self.yoffset = 0

        self.odimgfig = ImgFig(self, self.canvas)
        self.slider_eng.valueChanged[int].connect(self.OnScrollEng)
        self.unitlabel = self.odimgfig.bar.getAxis("right")

        self.pbExpData.clicked.connect(self.OnSaveData)
        self.pbExpImg.clicked.connect(self.OnSaveImage)
        self.pbCopy.clicked.connect(self.OnCopy)

        self.MetricCheckBox.toggled.connect(lambda: self.odimgfig.OnMetricScale(self.MetricCheckBox.isChecked(), self.ZeroOriginCheckBox.isChecked(),self.SquarePxCheckBox.isChecked()))
        self.ZeroOriginCheckBox.toggled.connect(lambda: self.odimgfig.OnMetricScale(self.MetricCheckBox.isChecked(), self.ZeroOriginCheckBox.isChecked(),self.SquarePxCheckBox.isChecked()))
        self.SquarePxCheckBox.toggled.connect(lambda: self.odimgfig.OnMetricScale(self.MetricCheckBox.isChecked(), self.ZeroOriginCheckBox.isChecked(),self.SquarePxCheckBox.isChecked()))
        self.SquarePxCheckBox.setVisible(False)
        self.ScalebarCheckBox.toggled.connect(lambda: self.odimgfig.OnUpdateScale(self.ScalebarCheckBox.isChecked()))

        #self.CropCheckBox.toggled.connect(lambda: self.OnCropCB(self.CropCheckBox.isChecked()))
        #self.cropflag = True
        #self.ShiftLabel.setText("x = %0.1f \ny = %0.1f" % (0, 0))
        #self.ODHighSpinBox.valueChanged.connect(lambda: self.setODlimits(self.ODLowSpinBox.value(),self.ODHighSpinBox.value()))
        #self.ODLowSpinBox.valueChanged.connect(lambda: self.setODlimits(self.ODLowSpinBox.value(),self.ODHighSpinBox.value()))
        self.pbRSTOD.clicked.connect(lambda: self.ShowMap(self.prelst, self.postlst))
        # self.pbClrShifts.clicked.connect(self.OnClrShifts)
        self.pbClrSel.clicked.connect(self.OnClrSelection)

        self.CMCatBox.addItems([self.cmaps[0][0],self.cmaps[1][0],self.cmaps[2][0],self.cmaps[3][0],self.cmaps[4][0],self.cmaps[5][0]])
        self.CMMapBox.addItems(self.cmaps[2][1])
        self.CMCatBox.setCurrentIndex(2)
        self.CMMapBox.setCurrentIndex(14)
        self.CMCatBox.currentIndexChanged.connect(self.odimgfig.OnCatChanged)
        self.CMMapBox.currentIndexChanged.connect(lambda: self.odimgfig.OnColormapChange(map=self.CMMapBox.currentText(),num_colors=self.StepSpin.value(),fliplut=self.ColorFlipCheckBox.isChecked()))
        self.StepSpin.valueChanged.connect(lambda: self.odimgfig.OnColormapChange(map=self.CMMapBox.currentText(),num_colors=self.StepSpin.value(),fliplut=self.ColorFlipCheckBox.isChecked()))
        self.ColorFlipCheckBox.toggled.connect(lambda: self.odimgfig.OnColormapChange(map=self.CMMapBox.currentText(), num_colors=self.StepSpin.value(),fliplut=self.ColorFlipCheckBox.isChecked()))
        self.filterSpinBox.valueChanged.connect(lambda: self.ShowMap(self.prelst, self.postlst))
        self.filterSpinBox.setEnabled(False)
        self.rb_filterlee.setEnabled(False)
        self.rb_filteruniform.setEnabled(False)
        self.rb_filterlee.toggled.connect(lambda: self.ShowMap(self.prelst, self.postlst))
        self.MapSelectWidget1.mousePressEvent = self.mouseEventOnQList
        self.MapSelectWidget1.mouseMoveEvent = self.mouseEventOnQList

        self.parent.tab_prep.button_spectralROI.setEnabled(False)

    def OnScrollEng(self, value):
        self.iev = value
        if self.com.stack_loaded == 1:
            self.slider_eng.setValue(value)
            if self.com.i0_loaded == 0:
                # Show flux image
                image = self.stk.absdata[:, :, self.iev]  # .copy()
            else:
                # Show OD image
                image = self.stk.od3d[:, :, self.iev]  # .copy()
            self.odimgfig.draw(image)

    def mouseEventOnQList(self, e):
        if e.type() == QtCore.QEvent.MouseMove or e.type() == QtCore.QEvent.MouseButtonPress:
            qlist = self.MapSelectWidget1
            pos = qlist.mapFromGlobal(QtGui.QCursor.pos())
            row = qlist.indexAt(pos).row()
            #print(row,self.latest_row)
            if self.stk.shifts:
                params = [self.stk.shifts[row][0],self.stk.shifts[row][1],self.stk.shifts[row][2]]
                if row >= 0:
                        if e.type() != QtCore.QEvent.MouseMove or row != self.latest_row:
                            if e.buttons() == QtCore.Qt.RightButton:
                                qlist.setCurrentRow(row)
                                params[1] = 1
                            elif e.buttons() == QtCore.Qt.LeftButton:
                                qlist.setCurrentRow(row)
                                params[1] = -1
                            self.qlistchanged.emit((row,params))
                            self.latest_row = row
        return
    def qListChangeHandler(self,paramtup):
        qlist = self.MapSelectWidget1
        row, params = paramtup
        if self.stk.shifts[row][1] == params[1]:
            #print("deselect it!")
            qlist.item(row).setBackground(QtGui.QColor(0, 0, 0, 0))
            self.stk.shifts[row][1] = 0
        elif self.stk.shifts[row][1] != params[1] and params[1] == 1: # select right
            #print("select right!")
            qlist.item(row).setBackground(QtGui.QColor('#7fc97f'))
            self.stk.shifts[row][1] = 1
        elif self.stk.shifts[row][1] != params[1] and params[1] == -1: # select left
            #print("select left!")
            qlist.item(row).setBackground(QtGui.QColor('#beaed4'))
            self.stk.shifts[row][1] = -1
        self.OnSelectionChanged()
        #print(row, params)
    def OnSelectionChanged(self):
        self.prelst = [index for index, value in enumerate([x[1] for x in  self.stk.shifts]) if value == -1]
        self.postlst = [index for index, value in enumerate([x[1] for x in  self.stk.shifts]) if value == 1]
        #print(self.prelst,self.postlst)

        if len(self.prelst) == 0 or len(self.postlst) == 0:
            if len(self.prelst + self.postlst) == 0:
                self.pbClrSel.setEnabled(False)
            else:
                self.pbClrSel.setEnabled(True)

            #self.odimgfig.loadData()
            self.odimgfig.OnColormapChange()
            self.OnScrollEng(self.MapSelectWidget1.currentRow())
            #            self.cm.clear()
            self.odimgfig.imageplot.titleLabel.setText("<center>Select at least one pre- and post-edge image!</center>",size='10pt')
            if self.com.i0_loaded == 0:
                self.unitlabel.setLabel(text="counts", units="")
            #self.ODHighSpinBox.setEnabled(False)
            #self.ODLowSpinBox.setEnabled(False)
            self.pbRSTOD.setEnabled(False)
            self.filterSpinBox.setEnabled(False)
            self.rb_filterlee.setEnabled(False)
            self.rb_filteruniform.setEnabled(False)
            self.pbExpData.setEnabled(False)
            self.pbExpImg.setEnabled(False)

        else:
            self.pbClrSel.setEnabled(True)
            #self.InfWarning = False
            self.unitlabel.setLabel(text="OD", units="")
            self.ShowMap(self.prelst,self.postlst)
            # self.OnScrollEng(self.MapSelectWidget1.currentRow())
            #self.OnMetricScale(self.MetricCheckBox.isChecked(), self.ZeroOriginCheckBox.isChecked(),
            #                        self.SquarePxCheckBox.isChecked())
        return
    # ----------------------------------------------------------------------
    def OnSaveData(self,event):
        #Save Data
        wildcard = "TIFF Float32 File (*.tif);; TXT Image (*.txt);;"

        fileName, _filter = QtWidgets.QFileDialog.getSaveFileName(self, 'Save OD Map', '', wildcard)

        fileName = str(fileName)
        if fileName == '':
            return
        path, ext = os.path.splitext(fileName)
        ext = ext[1:].lower()
        if ext == '':
            ext = _filter.split()[0].lower()[0:3]
            fileName = fileName+'.'+ext
        print(ext)
        if ext != 'tif' and ext != 'txt':
            error_message = (
                  'Only the TIF and TXT data formats are supported.\n'
                 'A file extension of `tif\' or `txt\' must be used.')
            QtWidgets.QMessageBox.warning(self, 'Error', 'Error - Could not save file.')
            return

        if ext == 'tif':
            from PIL import Image
            img1 = Image.fromarray(np.rot90(self.OD))
            img1.save(fileName)
        if ext == 'txt':
            np.savetxt(fileName, np.rot90(self.OD), delimiter='\t', newline='\n',fmt='%.5f')
    def OnCopy(self):
        self.exp = pg.exporters.ImageExporter(self.odimgfig.imageplot)
        self.exp.export(copy=True)
        return
    def OnSaveImage(self, event):
        wildcard = "TIFF (*.tif);;PNG (*.png);;JPG (*.jpg);;SVG (*.svg);;"
        fileName, _filter = QtWidgets.QFileDialog.getSaveFileName(self, 'Save OD Map', '', wildcard)

        fileName = str(fileName)
        if fileName == '':
            return
        path, ext = os.path.splitext(fileName)
        ext = ext[1:].lower()
        if ext == '':
            ext = _filter.split()[0].lower()[0:3]
            fileName = fileName+'.'+ext

        if ext == 'svg':
            exp = pg.exporters.SVGExporter(self.odimgfig.imageplot)
        else:
            exp = pg.exporters.ImageExporter(self.odimgfig.imageplot)
        if ext in ['tif','png','jpg','svg']:
            exp.export(fileName)

    def OnCatChanged(self):
        self.CMMapBox.blockSignals(True)
        self.CMMapBox.clear()
        self.CMMapBox.blockSignals(False)
        self.CMMapBox.addItems(self.cmaps[self.CMCatBox.currentIndex()][1])

    def setODlimits(self,low, high):
        #self.ODHighSpinBox.setMaximum(self.ODmax)
        #self.ODHighSpinBox.setMinimum(low)
        #self.ODLowSpinBox.setMaximum(high)
        #self.ODLowSpinBox.setMinimum(self.ODmin)
        OD = np.clip(self.OD,low,high)
        #self.OnColormap(map=self.CMMapBox.currentText(), colors=self.StepSpin.value())
        #self.setODbar(low, high)
        self.odimgfig.draw(OD,False,True)
        #self.m_item.setImage(OD)
    # def setCrosshair(self):
    #     if hasattr(self, 'vLine'):
    #         self.p2.removeItem(self.vLine)
    #         self.p2.removeItem(self.vLine)
    #         self.canvas.removeItem(self.ODlabel)
    #
    #     self.vLine = pg.InfiniteLine(angle=90, movable=False, pen="0000")
    #     self.hLine = pg.InfiniteLine(angle=0, movable=False,pen="0000")
    #     self.ODlabel = pg.LabelItem(justify="left")
    #     self.canvas.addItem(self.ODlabel)
    #     self.p2.addItem(self.vLine, ignoreBounds=True)
    #     self.p2.addItem(self.hLine, ignoreBounds=True)
    #     self.proxy = pg.SignalProxy(self.p2.scene().sigMouseMoved, rateLimit=30, slot=self.OnMouseHover)

    # def OnMouseHover(self,evt):
    #     pos = evt[0]
    #     #print(pos)
    #     if hasattr(self, 'OD'):
    #         if self.p2.getViewBox().itemBoundingRect(self.m_item).contains(self.p2.getViewBox().mapSceneToView(pos)):
    #             self.vLine.setPen("r")
    #             self.hLine.setPen("r")
    #             self.vLine.setZValue(1000)
    #             self.hLine.setZValue(1000)
    #             mousePoint = self.p2.getViewBox().mapSceneToView(pos)
    #             if self.MetricCheckBox.isChecked():
    #                 x_min = (self.p2.getViewBox().itemBoundingRect(self.m_item).x() / (self.scale * self.stk.x_pxsize)) # minimum x value in px
    #                 y_min = (self.p2.getViewBox().itemBoundingRect(self.m_item).y() / (self.scale * self.stk.y_pxsize)) # minimum y value in px
    #                 x_off = x_min - int(x_min)
    #                 y_off = y_min - int(y_min)
    #                 x_pos = mousePoint.x() / (self.scale * self.stk.x_pxsize) +1 #x mousepos in number of rows
    #                 y_pos = mousePoint.y() / (self.scale * self.stk.y_pxsize) +1 #y mousepos in number of rows
    #
    #                 if mousePoint.x() < 0:
    #                     xi = int(x_pos - x_off-1)
    #                 else:
    #                     xi = int(x_pos - x_off)
    #                 if mousePoint.y() < 0:
    #                     yi = int(y_pos - y_off -1)
    #                 else:
    #                     yi = int(y_pos - y_off)
    #
    #                 x = (xi - 0.5 + x_off) * self.scale * self.stk.x_pxsize
    #                 y = (yi - 0.5 + y_off) * self.scale * self.stk.y_pxsize
    #                 pt = QtCore.QPointF(x,y)
    #                 self.ODlabel.setText("<span style='font-size: 10pt; color: red'> x = %0.3f  µm, <span style='color: red'> y = %0.3f µm</span>, <span style='color: red'> OD = %0.2f</span>" % (x/self.scale - 0.5 * self.stk.x_pxsize, y/self.scale - 0.5 * self.stk.y_pxsize, self.OD[xi-1-int(x_min),yi-1-int(y_min)]))
    #             else:
    #                 x = int(mousePoint.x()) + 0.5
    #                 y = int(mousePoint.y()) + 0.5
    #                 pt = QtCore.QPointF(x,y)
    #                 self.ODlabel.setText("<span style='font-size: 10pt; color: red'> x = %0.0f, <span style='color: red'> y = %0.0f</span>, <span style='color: red'> OD = %0.2f</span>" % (x, y,self.OD[int(x),int(y)]))
    #             #print(x,y)
    #             self.vLine.setPos(x)
    #             self.hLine.setPos(y)
    #         else:
    #             self.ODlabel.setText("")
    #             self.vLine.setZValue(-1000)
    #             self.hLine.setZValue(-1000)
    #             self.vLine.setPen("0000")
    #             self.hLine.setPen("0000")

    # def setODbar(self,min=None,max=None):
    #     self.cm.setRange(xRange=[0,1], yRange=[min,max], update=False, disableAutoRange=True,padding=0)
    #     self.cmimg.setRect(QtCore.QRectF(0,min,1,max-min))

    # ToDo: Restore fine-alignment by arrow keys
    # def keyPressEvent(self, e):
    #     modifiers = QtWidgets.QApplication.keyboardModifiers()
    #     if modifiers == QtCore.Qt.ShiftModifier:
    #         noshift = False
    #     elif modifiers == QtCore.Qt.KeypadModifier:
    #         noshift = True
    #     elif modifiers == (QtCore.Qt.KeypadModifier |
    #                        QtCore.Qt.ShiftModifier):
    #         noshift = False
    #     else:
    #         noshift = True
    #     if e.key() == 67 and (e.modifiers() & QtCore.Qt.ControlModifier):
    #         self.OnCopy()
    #     if e.key() == Qt.Key_Up or (e.key() == QtCore.Qt.Key_8):
    #         if noshift:
    #             self.pbUU.click()
    #         else:
    #             self.pbU.click()
    #     elif e.key() == Qt.Key_Down or (e.key() == QtCore.Qt.Key_2):
    #         if noshift:
    #             self.pbDD.click()
    #         else:
    #             self.pbD.click()
    #     elif e.key() == Qt.Key_Left or (e.key() == QtCore.Qt.Key_4):
    #         if noshift:
    #             self.pbLL.click()
    #         else:
    #             self.pbL.click()
    #     elif e.key() == Qt.Key_Right or (e.key() == QtCore.Qt.Key_6):
    #         if noshift:
    #             self.pbRR.click()
    #         else:
    #             self.pbR.click()
    #     elif e.key() == QtCore.Qt.Key_Home:
    #         if noshift:
    #             self.pbLLUU.click()
    #         else:
    #             self.pbLU.click()
    #     elif e.key() == QtCore.Qt.Key_PageUp:
    #         if noshift:
    #             self.pbRRUU.click()
    #         else:
    #             self.pbRU.click()
    #     elif e.key() == QtCore.Qt.Key_End:
    #         if noshift:
    #             self.pbLLDD.click()
    #         else:
    #             self.pbLD.click()
    #     elif e.key() == QtCore.Qt.Key_PageDown:
    #         if noshift:
    #             self.pbRRDD.click()
    #         else:
    #             self.pbRD.click()
    #     elif e.key() == QtCore.Qt.Key_Clear:
    #             self.pbRST.click()

#     def Clear(self):
#         try:
#             self.slider_eng.valueChanged.disconnect()
#             self.qlistchanged.disconnect()
#             self.pbSelfromSpec.disconnect()
#         except:
#             pass
#         self.prelst = []
#         self.postlst =[]
#         self.stk.shifts = []
#         self.stk.absdata_shifted= []
# #        self.p1.clear()
#         self.p2.clear()
#         self.cm.clear()
#         self.MapSelectWidget1.clear()

    def loadData(self): # Called when fresh data are loaded.
        self.stk.shifts = []
        #self.slider_eng.setRange(0, self.stk.n_ev - 1)
        #self.ODHighSpinBox.setEnabled(False)
        #self.ODLowSpinBox.setEnabled(False)
        #self.pbRSTOD.setEnabled(False)
        self.filterSpinBox.setEnabled(False)
        self.pbExpData.setEnabled(False)
        self.pbExpImg.setEnabled(False)
        self.pbClrSel.setEnabled(False)
        self.stk.absdata_shifted = self.stk.absdata.copy()
        #self.p1.addItem(self.i_item)
        for i,e in enumerate(self.stk.ev): # Fill QList with energies
            self.stk.shifts.append([1,0,(0.0,0.0)]) #checked [0,1]; pre, post, undefined state for map [-1,1,0],(xshift [float],yshift [float])
            item = QtWidgets.QListWidgetItem(str(int(i)).zfill(3)+"     at     " + format(e, '.2f') + " eV     "+"+0.0"+"    +0.0")
            self.MapSelectWidget1.addItem(item)
        #self.slider_eng.valueChanged[int].connect(self.OnScrollEng)
        self.qlistchanged.connect(self.qListChangeHandler)
        self.odimgfig.OnColormapChange()

        #self.OnScrollEng(0) # Plot first image & set Scrollbar
        #self.OnMetricScale(self.MetricCheckBox.isChecked(), True, False)
        self.pbSelfromSpec.setEnabled(True)
        self.pbSelfromSpec.clicked.connect(self.OnSelfromSpec)
    # def UpdateEntry(self,row):
    #     self.MapSelectWidget1.item(row).setText(str(int(row)).zfill(3)+"     at     " + format(self.stk.ev[row], '.2f') + " eV     "+format(self.stk.shifts[row][2][0], '+.1f')+"    "+format(self.stk.shifts[row][2][1], '+.1f'))
    #     #self.MapSelectWidget1.addItem(self.MapSelectWidget1.item(row))
    #     self.i_item.setImage(self.Shift(row))
    # def ResetAllItems(self,widget):
    #     for i in range(widget.count()):
    #         widget.item(i).setForeground(QtGui.QColor(0, 0, 0, 128))
    # def OnScrollEng(self, value):
    #     self.slider_eng.blockSignals(True)
    #     self.slider_eng.setValue(value)
    #     self.slider_eng.blockSignals(False)
    #     self.MapSelectWidget1.setCurrentRow(value)
    #     self.ResetAllItems(self.MapSelectWidget1)
    #     self.iev = value
    #     self.p1.titleLabel.item.setTextWidth(self.p1.width() * 0.7)
    #     self.p2.titleLabel.item.setTextWidth(self.p2.width() * 0.7)
    #     #self.canvas.resizeEvent(None)
    #     if self.com.stack_loaded == 1:
    #         self.p1.titleLabel.setText("<center>Image at {0:5.2f} eV</center>".format(float(self.stk.ev[self.iev])),size='10pt')
    #         self.i_item.setImage(self.Shift(int(self.iev)))
    #         self.MapSelectWidget1.item(value).setForeground(QtGui.QColor(0, 0, 0, 255))

    # def OnCropCB(self, value=True):
    #     if self.com.stack_loaded == 1:
    #         if value == True:
    #             self.cropflag = True
    #         else:
    #             self.cropflag = False
    #         self.ShowMap(self.prelst, self.postlst)

    def OnClrSelection(self):
        for row in self.prelst + self.postlst:
            self.MapSelectWidget1.item(row).setBackground(QtGui.QColor(0, 0, 0, 0))
            self.stk.shifts[row][1] = 0
        self.OnSelectionChanged()

    def OnSelfromSpec(self):
        #print(len(self.stk.shifts)) #ToDo: Allow multiple map windows. Task: prevent that stk.shifts is appended each time the window is opened.
        spectralimgmap = SpectralImageMap(self, self.com, self.stk)
        spectralimgmap.show()

    # def OnClrShifts(self):
    #     for row in [index for index, value in enumerate([x[2] for x in  self.stk.shifts]) if value != (0.0,0.0)]:
    #         self.stk.shifts[row].pop(2)  # remove tuple
    #         self.stk.shifts[row].insert(2, (0.0, 0.0))
    #         self.stk.absdata_shifted[:, :, row] = self.stk.absdata[:, :, row]
    #         self.MapSelectWidget1.item(row).setText(
    #             str(int(row)).zfill(3) + "     at     " + format(self.stk.ev[row], '.2f') + " eV     " + format(
    #                 self.stk.shifts[row][2][0], '+.1f') + "    " + format(self.stk.shifts[row][2][1], '+.1f'))
    #     if hasattr(self, 'prelst') and hasattr(self, 'postlst'):
    #         if len(self.prelst) != 0 and len(self.postlst) != 0:
    #             self.ShowMap(self.prelst,self.postlst)
    #     if self.stk.shifts:
    #         self.UpdateEntry(self.MapSelectWidget1.currentRow())
    #
    # def setShifts(self,shift_x, shift_y):
    #     #print("setshifts called")
    #     # if hasattr(self, "OD"):
    #     if self.stk.shifts:
    #         row = self.MapSelectWidget1.currentRow()
    #         xoffset, yoffset = self.stk.shifts[row][2] # current offset stored as tuple in table stk.shifts
    #         #print(self.stk.shifts[self.MapSelectWidget1.currentRow()])
    #         #yoffset = self.stk.shifts[self.MapSelectWidget1.currentRow()][2][1]
    #         if shift_x == 0 and shift_y == 0: # if reset button pressed
    #             if xoffset != 0 or yoffset != 0:
    #                 self.stk.shifts[row].pop(2) # remove tuple
    #                 self.stk.shifts[row].insert(2,(0.0,0.0))
    #             else:
    #                 print("Reset has no effect")
    #                 return
    #         else:
    #             self.stk.shifts[row].pop(2)
    #             xoffset = round(xoffset + shift_x,1)
    #             yoffset = round(yoffset + shift_y,1)
    #             self.stk.shifts[row].insert(2,(xoffset, yoffset))
    #         #current_img = self.stk.absdata[:, :, self.MapSelectWidget1.currentRow()]
    #         #print(type(self.stk.absdata), self.stk.shifts[self.MapSelectWidget1.currentRow()])
    #         #self.Shift(row)
    #         self.UpdateEntry(row)
    #         if hasattr(self, 'prelst') and hasattr(self, 'postlst'):
    #             if len(self.prelst) != 0 and len(self.postlst) != 0:
    #                 self.ShowMap(self.prelst,self.postlst)
    # def Shift(self,row):
    #     #current_img = self.stk.absdata_shifted[:, :, row]
    #     original_img = self.stk.absdata[:, :, row]
    #     xoffset, yoffset = self.stk.shifts[row][2] # x and y offsets of current image
    #     if xoffset == 0 and yoffset == 0:
    #         self.stk.absdata_shifted[:, :, row] = original_img # replace with original if no shift is applied
    #         #self.ShiftLabel.setText("x = %0.1f \ny = %0.1f" % (self.xoffset, self.yoffset))
    #         #return original_img
    #     else:
    #         shifted = ndimage.fourier_shift(np.fft.fft2(original_img), [float(xoffset), float(yoffset)])
    #         shifted = np.fft.ifft2(shifted)
    #         shifted_real = shifted.real
    #         self.stk.absdata_shifted[:, :, row] = shifted_real
    #     return self.stk.absdata_shifted[:,:, row]

            #self.ShiftLabel.setText("x = %0.1f \ny = %0.1f" % (xoffset, yoffset))
    def CalcODMap(self,im_idx1,im_idx2):
        if len(im_idx1) == 1 and len(im_idx2) == 1:
            im1 = self.stk.absdata_shifted[:, :, int(im_idx1[0])]
            im2 = self.stk.absdata_shifted[:, :, int(im_idx2[0])]
        else:
            im1 = np.mean([self.stk.absdata_shifted[:, :, i] for i in im_idx1], axis=0)
            im2 = np.mean([self.stk.absdata_shifted[:, :, i] for i in im_idx2], axis=0)
        OD = np.log(im1/im2)
        # if self.cropflag:
        #     shiftlst = [self.stk.shifts[i][2] for i in im_idx1 + im_idx2]
        #     #print(shiftlst)
        #     # print(max(shiftlst, key=itemgetter(0))[0]) # possible alternative to lambda function
        #     #Calculate crops
        #     l = int(np.floor(min(shiftlst, key=lambda item: item[0])[0]))
        #     cl = l if l < 0 else None
        #     r = int(np.ceil(max(shiftlst, key=lambda item: item[0])[0]))
        #     t = int(np.ceil(max(shiftlst, key=lambda item: item[1])[1]))
        #     b = int(np.floor(min(shiftlst, key=lambda item: item[1])[1]))
        #     cb = b if b < 0 else None
        #     # Crop map
        #     OD = OD[r:cl,t:cb]
        #     #print(r,cl,t,cb)
        inf_idx = np.where(np.isinf(OD))
        nan_idx = np.where(np.isnan(OD))
        if np.any(inf_idx) or np.any(nan_idx):
            #if not self.InfWarning:
            #    self.InfWarning = True
            #    QtWidgets.QMessageBox.warning(self, 'Warning!', "The OD map contained infinite or nan values. Please note that they have been zeroed.")
            OD[inf_idx] = 0 #infinite values get replaced by zero. This is an ugly work around, but with nans the image is not displayed correctly.
            OD[nan_idx] = 0
        return OD

    # def calcBinSize(self,i,N):
    #     return int(round(256*(i+1)/N) - round(256*i/N))

    # def OnColormap(self,map="afmhot", colors=256):
    #     if hasattr(self, "OD"):
    #         colormap = cm.get_cmap(map, colors)
    #         colormap = colormap(np.arange(colors))
    #         cm_lst = [[colormap[idx][0], colormap[idx][1], colormap[idx][2], colormap[idx][3]] for idx in range(np.shape(colormap)[0])] #convert to r,g,b,a list
    #         cm_lst = [item for sub in [[cm_lst[i]]*self.calcBinSize(i,colors) for i in range(colors)] for item in sub] #fills 256 bins as equal as possible with n colors
    #         cm_array = np.array([np.asarray(cm_lst)]) #vertical colorbar
    #         cm_lst.extend((cm_lst[-1],cm_lst[-1],cm_lst[-1]))
    #         lut = np.asarray(cm_lst)
    #         lut = (lut * 255).view(np.ndarray) #lut for OD map
    #         self.cmimg.setImage(cm_array)
    #         self.m_item.setLookupTable(lut)
    #         if hasattr(self, "OD"):
    #             self.setODbar(self.ODmin, self.ODmax)

    def ShowMap(self, preidx, postidx):
        # self.p2.clear()
        # self.p2.addItem(self.m_item)
        # self.cm.clear()
        # self.cm.addItem(self.cmimg)

        #self.setCrosshair()
        try:
            ## Optional sorting switched off for convenience. Allows maps to range to negative OD
            #selection = preidx + postidx
            # selection.sort()
            if len(preidx + postidx) == 2:
                self.odimgfig.imageplot.titleLabel.setText("<center>Binary map from energies " + str(round(self.stk.ev[preidx[0]], 2)) + " and " + str(
                round(self.stk.ev[postidx[0]], 2)) + " eV</center>",
                                           size='10pt')
            elif len(preidx + postidx) <= 6:
                self.odimgfig.imageplot.titleLabel.setText("<center>Map from energies " + str([round(self.stk.ev[e], 2) for e in preidx]).strip('[]') + " and " +
                                 str([round(self.stk.ev[e], 2) for e in postidx]).strip('[]') + " eV</center>",
                                           size='10pt')
            else:
                self.odimgfig.imageplot.titleLabel.setText("<center>Map from "+ str(len(preidx + postidx)) +" energies: " + str(round(self.stk.ev[preidx[0]], 2)) + ' ... ' + str(round(self.stk.ev[preidx[-1]], 2)) + " and "
                                 + str(round(self.stk.ev[postidx[0]], 2)) + ' ... ' + str(round(self.stk.ev[postidx[-1]], 2)) + " eV</center>",
                                           size='10pt')
            self.OD = self.CalcODMap(preidx, postidx)
            if self.filterSpinBox.value() > 1:
                if self.rb_filteruniform.isChecked():
                    self.OD = ndimage.filters.uniform_filter(self.OD, size=self.filterSpinBox.value(), mode='nearest')
                elif self.rb_filterlee.isChecked():
                    self.OD = self.leeFilter(self.OD, size=self.filterSpinBox.value())
            self.ODmin = np.min(self.OD)
            self.ODmax = np.max(self.OD)
            #self.ODHighSpinBox.setEnabled(True)
            #self.ODLowSpinBox.setEnabled(True)
            #self.ODHighSpinBox.blockSignals(True)
            #self.ODLowSpinBox.blockSignals(True)
            #self.ODHighSpinBox.setMaximum(self.ODmax)
            #self.ODHighSpinBox.setMinimum(self.ODmin)
            #self.ODLowSpinBox.setMaximum(self.ODmax)
            #self.ODLowSpinBox.setMinimum(self.ODmin)
            #self.ODLowSpinBox.setValue(self.ODmin)
            #self.ODHighSpinBox.setValue(self.ODmax)
            self.setODlimits(self.ODmin, self.ODmax)
            #self.ODHighSpinBox.blockSignals(False)
            #self.ODLowSpinBox.blockSignals(False)
            self.pbRSTOD.setEnabled(True)
            self.filterSpinBox.setEnabled(True)
            self.rb_filterlee.setEnabled(True)
            self.rb_filteruniform.setEnabled(True)
            self.pbExpData.setEnabled(True)
            self.pbExpImg.setEnabled(True)
        #     #self.pglayout.layout.setColumnMaximumWidth(1, self.p1.width())
        except IndexError:
            pass
        #     self.p2.clear()
        #     self.cm.clear()
            #self.p2.setTitle("Please select a second image!")
            #print("Select a second image!")

    def leeFilter(self, img, size):
        mean = ndimage.filters.uniform_filter(img, (size, size),mode='nearest')
        sqrmean = ndimage.filters.uniform_filter(img ** 2, (size, size),mode='nearest')
        var = sqrmean - mean ** 2
        totvar = np.var(img)

        weights = var / (var + totvar)
        img = mean + weights * (img - mean)
        return img

    def closeEvent(self, event):
        #self.Clear()
        self.close()
        self.parent.tab_prep.button_spectralROI.setEnabled(True)
#-----------------------------------------------------------------------
