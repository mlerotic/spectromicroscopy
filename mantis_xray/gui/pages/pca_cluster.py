import os
import csv
import numpy as np
import pyqtgraph as pg
from PIL import Image
from PyQt5 import QtWidgets, uic, QtCore, QtGui
from ..widgets import SpecFig, ImgFig
from ..dialogs.save import SaveWinP2, SaveWinP3
from ..dialogs.scatter import Scatterplots
from ...helpers import PDFExporter

class PagePCACluster(QtWidgets.QWidget):
    def __init__(self, common, data_struct, stack, anlz):
        super(PagePCACluster, self).__init__()
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pagepca.ui'), self)
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

        self.initUI(common, data_struct, stack, anlz)

#----------------------------------------------------------------------
    def initUI(self, common, data_struct, stack, anlz):

        self.spectrum_plotwidget.setBackground("w")
        self.pca_eigenvals_plotwidget.setBackground("w")

        self.data_struct = data_struct
        self.stk = stack
        self.com = common
        self.common = common  # alias for SaveWinP2 compatibility
        self.anlz = anlz
        self.calcpca = False
        self.selpca = 0

        self.filename = " "

        self.ix = 0
        self.iy = 0
        self.iev = 0
        self.itheta = 0
        self.showflux = True
        self.show_colorbar = 0

        self.dispbrightness_min = 0
        self.dispbrightness_max = 100
        self.displaygamma = 10.0
        self.defaultdisplay = 1.0

        self.brightness_min = 0.0
        self.brightness_max = 1.0
        self.gamma = 1.0

        self.colortable = "gray"

        self.addroi = 0
        self.showROImask = 0
        self.line = None
        self.ROIpix = None

        self.show_scale_bar = 1
        self.white_scale_bar = 0

        self.movie_playing = 0
        self.mean_visible = 0

        # Cluster state
        self.cluster_display_mode = False
        self.selcluster = 0
        self.numclusters = 0
        self.init_nclusters = 5
        self.wo_1st_pca = 0
        self.sigma_split = 0
        self.showallspectra = 0
        self.showerrormap = 0
        self._cluster_display_opts = None
        self.pcscalingfactor = 0.0

        self.climgfig = ImgFig(self, self.canvas)
        self.specfig = SpecFig(self, self.spectrum_plotwidget)
        self.setupEigenvaluePlot()

        self.slider_cl.valueChanged[int].connect(self.OnScrollCl)
        self.slider_theta.valueChanged[int].connect(self.OnScrollTheta)
        self.slider_theta.setVisible(False)
        self.slider_theta.setEnabled(False)

        self.pb_copy_img.clicked.connect(self.climgfig.OnCopy)
        self.pb_copy_img.setEnabled(False)
        self.pb_copy_specimg.clicked.connect(self.OnCopyPCAPlot)
        self.pb_copy_specimg.setEnabled(False)

        self.SquarePxCheckBox.toggled.connect(lambda: self.climgfig.OnMetricScale(False, True, self.SquarePxCheckBox.isChecked()))
        self.SquarePxCheckBox.setVisible(True)
        self.ScalebarCheckBox.toggled.connect(lambda: self.climgfig.OnUpdateScale(self.ScalebarCheckBox.isChecked()))

        self.CMCatBox.addItems([self.cmaps[0][0], self.cmaps[1][0], self.cmaps[2][0], self.cmaps[3][0], self.cmaps[4][0], self.cmaps[5][0]])
        self.CMMapBox.addItems(self.cmaps[3][1])
        self.CMCatBox.currentIndexChanged.connect(self.climgfig.OnCatChanged)
        self.CMMapBox.currentIndexChanged.connect(lambda: self.climgfig.OnColormapChange(map=self.CMMapBox.currentText(), num_colors=self.StepSpin.value(), fliplut=self.ColorFlipCheckBox.isChecked()))
        self.CMCatBox.setCurrentIndex(3)
        self.CMMapBox.setCurrentIndex(5)
        self.StepSpin.valueChanged.connect(lambda: self.climgfig.OnColormapChange(map=self.CMMapBox.currentText(), num_colors=self.StepSpin.value(), fliplut=self.ColorFlipCheckBox.isChecked()))
        self.ColorFlipCheckBox.toggled.connect(lambda: self.climgfig.OnColormapChange(map=self.CMMapBox.currentText(), num_colors=self.StepSpin.value(), fliplut=self.ColorFlipCheckBox.isChecked()))
        self.CMMapBox.currentIndexChanged.connect(self.OnClusterColormapChanged)
        self.ColorFlipCheckBox.toggled.connect(self.OnClusterColormapChanged)
        self.StepSpin.valueChanged.connect(self.OnClusterColormapChanged)

        self.button_calcpca.clicked.connect(self.OnCalcPCA)
        self.button_calcpca.setEnabled(False)

        self.button_savepca.clicked.connect(self.OnSave)
        self.button_savepca.setEnabled(False)

        self.button_movepcup.clicked.connect(self.OnMovePCUP)
        self.button_movepcup.setEnabled(False)

        self.npcaspin.valueChanged[int].connect(self.OnNPCAspin)

        self.button_toggle_pca_plot.toggled.connect(self.OnTogglePCAPlot)
        self.button_toggle_pca_plot.setEnabled(False)

        # Cluster UI wiring
        self.button_calcca.clicked.connect(self.OnCalcClusters)
        self.button_calcca.setEnabled(False)
        self.nclusterspin.setValue(self.init_nclusters)
        self.nclusterspin.valueChanged[int].connect(self.OnNClusterspin)

        # Dynamically add the missing cluster option used by refresh/reset logic.
        self.showallspectracb = QtWidgets.QCheckBox('Show all spectra', self)
        self.showallspectracb.setChecked(False)
        self.showallspectracb.setEnabled(False)
        self.showallspectracb.stateChanged.connect(self.OnShowallspectra)
        if hasattr(self, 'gridLayout_5') and self.gridLayout_5 is not None:
            # Reflow rows to avoid overlap with existing controls from the .ui file.
            grid = self.gridLayout_5
            grid.removeWidget(self.remove1stpcacb)
            grid.removeWidget(self.cb_splitclusters)
            grid.removeWidget(self.ntc_pcscaling)
            grid.removeWidget(self.ntc_pcscaling_2)
            grid.addWidget(self.showallspectracb, 2, 0, 1, 2)
            grid.addWidget(self.remove1stpcacb, 3, 0, 1, 2)
            grid.addWidget(self.cb_splitclusters, 4, 0, 1, 2)
            grid.addWidget(self.ntc_pcscaling, 5, 0)

            # Replace text input with QDoubleSpinBox for PC scaling factor
            self.ntc_pcscaling_2 = QtWidgets.QDoubleSpinBox()
            self.ntc_pcscaling_2.setRange(0.0, 1.0)
            self.ntc_pcscaling_2.setSingleStep(0.1)
            self.ntc_pcscaling_2.setValue(0.0)
            self.ntc_pcscaling_2.setDecimals(1)
            grid.addWidget(self.ntc_pcscaling_2, 5, 1)

        # Cluster block must stay inactive until PCA results exist.
        self.nclusterspin.setEnabled(False)
        self.remove1stpcacb.setEnabled(False)
        self.cb_splitclusters.setEnabled(False)
        self.ntc_pcscaling_2.setEnabled(False)

        # Set PC scaling factor initial value to 0.0
        self.ntc_pcscaling_2.setValue(0.0)

        self.npcaspin.setMaximumWidth(50)
        self.nclusterspin.setMaximumWidth(50)
        self.ntc_pcscaling_2.setMaximumWidth(50)
        self.remove1stpcacb.stateChanged.connect(self.OnRemove1stpca)
        self.cb_splitclusters.stateChanged.connect(self.OnSplitClusters)
        self.button_scatterplots.clicked.connect(self.OnShowScatterplots)
        self.button_scatterplots.setEnabled(False)
        self.button_savecluster.clicked.connect(self.OnSaveCluster)
        self.button_savecluster.setEnabled(False)
        self.button_errormap.setCheckable(True)
        self.button_errormap.setChecked(False)
        self.button_errormap.toggled.connect(self.OnShowErrorMap)
        self.button_errormap.setEnabled(False)
        self._setErrorMapButtonText(False)
        self.canvas.scene().sigMouseClicked.connect(self.OnPointClusterImage)

        self.UpdatePCAPlotView()

#----------------------------------------------------------------------
    def setupEigenvaluePlot(self):
        self.evals_plot_item = self.pca_eigenvals_plotwidget.getPlotItem()
        self.evals_plot_item.layout.setSpacing(12)
        self.evals_plot_item.layout.setContentsMargins(10, 10, 40, 10)
        self.evals_plot_item.showGrid(y=True)
        self.evals_plot_item.showAxis("top", show=True)
        self.evals_plot_item.showAxis("right", show=True)
        right_axis = self.evals_plot_item.getAxis("right")
        top_axis = self.evals_plot_item.getAxis("top")
        right_axis.setStyle(showValues=False, tickLength=0)
        top_axis.setStyle(showValues=False, tickLength=0)
        self.evals_plot_item.getAxis("bottom").setLabel(text="Principal component")
        self.evals_plot_item.getAxis("left").setLabel(text="log\u2081\u2080(Eigenvalue)")
        # Keep log mode OFF – ScatterPlotItem does not honour pyqtgraph's
        # automatic log-axis transformation.  We pre-compute log10 ourselves.
        self.evals_plot_item.setLogMode(x=False, y=False)
        self.evals_plot_item.addLegend()

#----------------------------------------------------------------------
    def _caption_prefix(self):
        if self.com.fntocaption and hasattr(self, "tc_file"):
            fn = self.tc_file.text().strip()
            if fn:
                return fn + " - "
        return ""

    def _cluster_composite_title(self):
        prefix = self._caption_prefix()
        if self.com.stack_4d == 1 and hasattr(self.stk, 'theta') and len(self.stk.theta) > self.itheta:
            return "{}Cluster composite at {:5.1f}\u00b0".format(prefix, float(self.stk.theta[self.itheta]))
        return "{}Cluster composite".format(prefix)

#----------------------------------------------------------------------
    def UpdatePCAPlotView(self):
        # Cluster mode has a fixed right-hand view: cluster spectra.
        if self.cluster_display_mode:
            self.pca_plot_stack.setCurrentWidget(self.pca_eigenvalues_page)
            self.button_toggle_pca_plot.setEnabled(False)
            self.button_movepcup.setEnabled(False)
            self.button_toggle_pca_plot.setText("Show eigenvalues")
            self.groupBox_6.setTitle("")
            return

        showing_right = self.button_toggle_pca_plot.isChecked()
        self.pca_plot_stack.setCurrentWidget(
            self.pca_eigenvalues_page if showing_right else self.pca_spectrum_page
        )
        self.button_toggle_pca_plot.setEnabled(self.com.pca_calculated == 1)
        self.button_toggle_pca_plot.setText(
            "Show PCA spectrum" if showing_right else "Show eigenvalues"
        )
        self.groupBox_6.setTitle("")

        # Move PC Up is only valid while looking at PCA eigenvalues.
        self.button_movepcup.setEnabled(self.com.pca_calculated == 1 and showing_right)

#----------------------------------------------------------------------
    def showEvals(self):
        self.evals_plot_item.clear()
        # Log mode stays OFF – ScatterPlotItem ignores pyqtgraph's log-axis
        # transform.  We pass pre-computed log10 values ourselves so that
        # the dots land at the correct positions.
        self.evals_plot_item.setLogMode(x=False, y=False)
        self.evals_plot_item.getAxis("bottom").setLabel(text="Principal component")
        self.evals_plot_item.getAxis("left").setLabel(text="log\u2081\u2080(Eigenvalue)")
        self.evals_plot_item.setTitle("")

        if self.anlz.pca_calculated != 1:
            return

        evalmax = min(self.stk.n_ev, 40, len(self.anlz.eigenvals))
        if evalmax <= 0:
            return

        x = np.arange(1, evalmax + 1)
        y_raw = np.asarray(self.anlz.eigenvals[:evalmax], dtype=float)
        # Use absolute values; clamp away from zero to keep log10 defined
        y_raw = np.clip(np.abs(y_raw), 1e-10, None)
        y = np.log10(y_raw)   # display values (log-space)

        selected_idx = int(np.clip(self.selpca, 0, evalmax - 1))

        scatter = pg.ScatterPlotItem(
            x=x,
            y=y,
            data=np.arange(evalmax),
            pen=None,
            brush=pg.mkBrush(31, 119, 180, 180),
            size=10,
        )
        scatter.sigClicked.connect(self.OnPointEvalsImage)
        self.evals_plot_item.addItem(scatter)

        selected = pg.ScatterPlotItem(
            x=[x[selected_idx]],
            y=[y[selected_idx]],
            pen=pg.mkPen('#d62728', width=1.5),
            brush=pg.mkBrush('#d62728'),
            size=14,
        )
        self.evals_plot_item.addItem(selected)

        y_min = float(np.min(y))
        y_max = float(np.max(y))
        y_span = max(y_max - y_min, 0.5)  # at least half a decade visible
        self.evals_plot_item.setXRange(0.5, evalmax + 0.5, padding=0)
        self.evals_plot_item.setYRange(y_min - 0.05 * y_span, y_max + 0.15 * y_span, padding=0)

#----------------------------------------------------------------------
    def OnMovePCUP(self):
        # selpca is 0-based in this tab
        if self.selpca == 0:
            return

        self.anlz.move_pc_up(self.selpca)

        self.selpca -= 1
        # setValue triggers OnScrollCl which updates image, spectrum,
        # eigenvalue plot and the view title in one go.
        self.slider_cl.setValue(self.selpca)


#----------------------------------------------------------------------
    def OnNPCAspin(self, value):
        self.numsigpca = value

        if self.com.pca_calculated == 1:
            self.anlz.numsigpca = self.numsigpca

        self._update_variance_label()

    def _get_active_variance(self):
        variance = getattr(self.anlz, 'variance', None)
        if variance is not None and np.size(variance) > 0:
            return np.asarray(variance)

        variance_4d = getattr(self.anlz, 'variance4D', None)
        if isinstance(variance_4d, (list, tuple)) and len(variance_4d) > 0:
            theta_idx = int(getattr(self, 'itheta', 0))
            theta_idx = int(np.clip(theta_idx, 0, len(variance_4d) - 1))
            return np.asarray(variance_4d[theta_idx])

        return np.asarray([])

    def _update_variance_label(self):
        variance = self._get_active_variance()
        if variance.size == 0:
            self.vartc.setText('N/A')
            return
        n_comp = int(np.clip(self.numsigpca, 0, variance.size))
        cumulative = variance[:n_comp].sum()
        self.vartc.setText(str(np.round(cumulative, 2) * 100) + '%')

#----------------------------------------------------------------------
    def OnTogglePCAPlot(self, checked):
        if self.cluster_display_mode:
            return

        self.UpdatePCAPlotView()
        if checked:
            if self.anlz.pca_calculated == 1:
                self.showEvals()
        else:
            if self.anlz.pca_calculated == 1:
                self.specfig.updatePlotDataOnPCA()

#----------------------------------------------------------------------
    def OnPointEvalsImage(self, item, points, event=None):
        if self.com.pca_calculated != 1 or not points:
            return

        selected_component = points[0].data()
        if selected_component is None:
            return

        self.slider_cl.setValue(int(selected_component))

#----------------------------------------------------------------------
    def OnCopyPCAPlot(self):
        if self.button_toggle_pca_plot.isChecked():
            exp = pg.exporters.ImageExporter(self.pca_eigenvals_plotwidget.scene())
            exp.export(copy=True)
            return
        self.specfig.OnCopy()

#----------------------------------------------------------------------
    def OnCalcPCA(self, event):

        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
        # Re-running PCA invalidates previously calculated cluster outputs.
        self.com.cluster_calculated = 0
        self.button_scatterplots.setEnabled(False)
        self.button_savecluster.setEnabled(False)
        self.button_errormap.setEnabled(False)
        self.showallspectracb.setEnabled(False)
        self.calcpca = False
        self.selpca = 1
        self.numsigpca = 2
        # Leaving cluster mode when PCA is recalculated
        self.cluster_display_mode = False
        #self.slidershow.setValue(self.selpca)

        #scrollmax = np.min([self.stk.n_ev, 20])
        #self.slidershow.setMaximum(scrollmax)

        #try:
        if self.com.stack_4d == 1:
            self.CalcPCA4D()
        else:
            self.CalcPCA()
        self.calcpca = True
        #self.loadPCAImage()
        #self.loadPCASpectrum()
        #self.showEvals()
        self.com.pca_calculated = 1
        QtWidgets.QApplication.restoreOverrideCursor()
        # Reset to diverging/RdBu – the appropriate colormap for signed PCA images.
        # Done here (not in UpdatePCAPlotView) so manual changes are not overridden
        # on every scroll/toggle, only when a fresh PCA is calculated.
        self.CMCatBox.setCurrentIndex(3)   # Diverging
        self.CMMapBox.setCurrentIndex(5)   # RdBu
        self.ColorFlipCheckBox.setChecked(False)
        self.StepSpin.setEnabled(True)     # unlock – was locked during cluster mode
        self.StepSpin.setValue(256)        # continuous gradient for PCA
        # Re-enable and reconfigure slider for PCA component scrolling
        scrollmax = max(self.numsigpca - 1, 0)
        self.slider_cl.blockSignals(True)
        self.slider_cl.setRange(0, scrollmax)
        self.slider_cl.setValue(0)
        self.slider_cl.blockSignals(False)
        self.slider_cl.setEnabled(True)
        self.showerrormap = 0
        self.button_errormap.blockSignals(True)
        self.button_errormap.setChecked(False)
        self.button_errormap.blockSignals(False)
        self._setErrorMapButtonText(False)
        self.ntc_pcscaling_2.setValue(0.0)
        self._cluster_display_opts = None
        self._setClusterColorbarHandlesEnabled(True)
        self._resetClusterColorbarTicks()
        self.specfig.ClearandReload(self.climgfig)
        if self.com.stack_4d == 1:
            self.slider_theta.setVisible(True)
            self.slider_theta.setEnabled(True)
            self.slider_theta.blockSignals(True)
            self.slider_theta.setRange(0, self.stk.n_theta - 1)
            self.slider_theta.setValue(self.itheta)
            self.slider_theta.blockSignals(False)
        else:
            self.slider_theta.setVisible(False)
            self.slider_theta.setEnabled(False)
        self.climgfig.loadPCAImage()
        self.showEvals()
        self.UpdatePCAPlotView()
        #except:
        #    pass
        #    self.com.pca_calculated = 0
        #    QtWidgets.QApplication.restoreOverrideCursor()
        #   QtWidgets.QMessageBox.warning(self, 'Error', 'PCA not calculated.')

        self.window().refresh_widgets()

#----------------------------------------------------------------------
    def CalcPCA(self):

        self.anlz.calculate_pca()

        #Scree plot criterion
        self.numsigpca = self.anlz.numsigpca

        self.npcaspin.setValue(self.numsigpca)

        self._update_variance_label()

    def _apply_theta_pca_state(self, theta_idx):
        pcaimages_4d = getattr(self.anlz, 'pcaimages4D', None)
        if not isinstance(pcaimages_4d, (list, tuple)) or len(pcaimages_4d) == 0:
            return False
        theta_idx = int(np.clip(theta_idx, 0, len(pcaimages_4d) - 1))
        self.anlz.pcaimages = self.anlz.pcaimages4D[theta_idx]
        self.anlz.eigenvals = self.anlz.eigenvals4D[theta_idx]
        self.anlz.eigenvecs = self.anlz.eigenvecs4D[theta_idx]
        self.anlz.variance = self.anlz.variance4D[theta_idx]
        self.anlz.pcaimagebounds = self.anlz.pcaimagebounds4D[theta_idx]
        self.itheta = theta_idx
        return True

#----------------------------------------------------------------------
    def CalcPCA4D(self):

        if self.com.stack_4d == 1:
            self.anlz.calculate_pca_4D()
        else:
            return

        #Scree plot criterion
        self.numsigpca = self.anlz.numsigpca
        self.npcaspin.setValue(self.numsigpca)
        self._apply_theta_pca_state(self.itheta)
        self.com.pca4D_calculated = 1
        self._update_variance_label()

        if self.com.spec_anl4D_calculated == 1:
            self.anlz.calculate_targetmaps_4D()

#-----------------------------------------------------------------------
    def OnScrollTheta(self, value):
        if self.com.stack_4d != 1:
            return

        self.itheta = int(value)
        self.slider_theta.blockSignals(True)
        self.slider_theta.setValue(self.itheta)
        self.slider_theta.blockSignals(False)

        if self.com.pca4D_calculated == 1 and self._apply_theta_pca_state(self.itheta):
            current_component = int(np.clip(self.selpca, 0, max(self.numsigpca - 1, 0)))
            if self.com.cluster_calculated == 1:
                self._reset_cluster_display_to_pca(current_component)
            if self.com.pca_calculated == 1:
                self._update_variance_label()
                self.OnScrollCl(current_component)
                self.window().refresh_widgets()

        # Keep theta slider behavior independent across tabs.

#-----------------------------------------------------------------------
    def OnScrollCl(self, value):
        if self.cluster_display_mode:
            self.selcluster = value
            # Left panel always remains the composite; scrollbar only affects right panel (spectrum)
            self.showClusterSpectrum()
            self.UpdatePCAPlotView()
            return

        self.selpca = value
        if self.anlz.pca_calculated == 1:
            self.slider_cl.setValue(value)
            image = self.anlz.pcaimages[:, :, self.selpca]
            self.climgfig.draw(image, levels=(-self.anlz.pcaimagebounds[self.selpca], self.anlz.pcaimagebounds[self.selpca]), setpcalabel=True)
            if self.com.stack_4d == 1 and hasattr(self.stk, 'theta') and len(self.stk.theta) > self.itheta:
                title = "<center>{}Component {:02d} at {:5.1f}\u00b0</center>".format(
                    self._caption_prefix(), int(self.selpca) + 1, float(self.stk.theta[self.itheta]))
                self.climgfig.imageplot.setTitle(title)
            self.specfig.updatePlotDataOnPCA()
            self.showEvals()
            self.UpdatePCAPlotView()

#----------------------------------------------------------------------
    def OnCalcClusters(self, event):
        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
        try:
            self.pcscalingfactor = self.ntc_pcscaling_2.value()
        except Exception:
            self.pcscalingfactor = 0.0

        try:
            self.CalcClusters()
            self.com.cluster_calculated = 1
            self.selcluster = 0
            self.cluster_display_mode = True

            # Set slider range for clusters
            self.slider_cl.blockSignals(True)
            self.slider_cl.setRange(0, self.numclusters - 1)
            self.slider_cl.setValue(0)
            self.slider_cl.blockSignals(False)

            # Start each CA run in composite-spectrum mode by default.
            self.showallspectracb.blockSignals(True)
            self.showallspectracb.setChecked(True)
            self.showallspectracb.blockSignals(False)
            self.showallspectra = 1
            self.slider_cl.setEnabled(False)

            # Switch colormap to Qualitative/Set1 with one step per cluster so
            # the image, colorbar and spectra all share the same palette.
            self.CMCatBox.setCurrentIndex(4)           # Qualitative
            self.StepSpin.setValue(self.numclusters)   # one discrete colour per cluster
            self.CMMapBox.setCurrentIndex(5)           # Set1
            self.ColorFlipCheckBox.setChecked(True)    # legacy CA style: flipped colours enabled
            self.StepSpin.setEnabled(False)            # locked to numclusters

            self.showClusterImage()
            self.showClusterSpectrum()
            self.UpdatePCAPlotView()
            self._setClusterColorbarHandlesEnabled(False)

            QtWidgets.QApplication.restoreOverrideCursor()
        except Exception as e:
            self.com.cluster_calculated = 0
            QtWidgets.QApplication.restoreOverrideCursor()
            print(e)
            QtWidgets.QMessageBox.warning(self, 'Error', 'Cluster analysis failed: {}'.format(e))

        self.window().refresh_widgets()

#----------------------------------------------------------------------
    def CalcClusters(self):
        nclusters = self.anlz.calculate_clusters(
        self.init_nclusters,
        self.wo_1st_pca,
        self.sigma_split,
        pcscalingfactor=self.pcscalingfactor
        )
        self.numclusters = nclusters
        self.ntc_clusters_found.setText(str(self.numclusters))

#----------------------------------------------------------------------
    def showClusterImage(self):
        """Display the cluster-index image (always composite) using the current Set1 colormap.
        cluster_indices values (0..numclusters-1) map directly to LUT entries.
        """
        if self.numclusters <= 0:
            return

        cluster_indices = self.anlz.cluster_indices.astype(float)
        self.climgfig.draw(
            cluster_indices,
            setlabel=False,
            # Integer cluster ids (0..N-1) become centered bins in the colorbar.
            levels=(-0.5, self.numclusters - 0.5),
            # Always rebuild the cluster colormap LUT so that any temporary
            # per-cluster binary LUT set during Save_CA() export is correctly reset.
            setlut=True,
        )
        self._setColorbarAxisLabel("Cluster", "")
        self._setClusterColorbarTicks()

        title = self._cluster_composite_title()
        self.climgfig.imageplot.setTitle("<center>{}</center>".format(title))

#----------------------------------------------------------------------
    def showClusterDistanceMap(self):
        """Display the cluster error map in the left panel."""
        if self.numclusters <= 0:
            return

        mapimage = np.asarray(self.anlz.cluster_distances, dtype=float)
        if mapimage.size == 0:
            return

        finite = mapimage[np.isfinite(mapimage)]
        if finite.size:
            low = float(np.nanmin(finite))
            high = float(np.nanmax(finite))
        else:
            low, high = 0.0, 1.0

        if high <= low:
            high = low + 1.0

        self.climgfig.draw(
            mapimage,
            setlabel=False,
            levels=(low, high),
        )
        self._resetClusterColorbarTicks()
        self._setColorbarAxisLabel("Cluster error", "a.u.")
        self.climgfig.OnColormapChange(
            map=self.CMMapBox.currentText(),
            num_colors=self.StepSpin.value(),
            fliplut=self.ColorFlipCheckBox.isChecked(),
        )
        prefix = self._caption_prefix()
        self.climgfig.imageplot.setTitle("<center>{}Cluster Error Map</center>".format(prefix))

#----------------------------------------------------------------------
    def OnShowErrorMap(self, checked):
        if self.com.cluster_calculated != 1:
            self.button_errormap.blockSignals(True)
            self.button_errormap.setChecked(False)
            self.button_errormap.blockSignals(False)
            self.showerrormap = 0
            self._setErrorMapButtonText(False)
            return

        self.showerrormap = 1 if checked else 0
        if checked:
            self._captureClusterDisplayOptions()
            self.CMCatBox.setCurrentIndex(2)  # Sequential (2)
            self.CMMapBox.setCurrentText('gray')
            self.ColorFlipCheckBox.setChecked(False)
            self.StepSpin.setEnabled(True)
            self.StepSpin.setValue(256)
            self.showClusterDistanceMap()
        else:
            self._restoreClusterDisplayOptions()
            self.showClusterImage()
            self._setErrorMapButtonText(bool(checked))

#----------------------------------------------------------------------
    def _setErrorMapButtonText(self, active):
        self.button_errormap.setText("Show cluster image" if active else "Show error map")

#----------------------------------------------------------------------
    def _captureClusterDisplayOptions(self):
        if self._cluster_display_opts is None:
            self._cluster_display_opts = {
                'cat_index': self.CMCatBox.currentIndex(),
                'map_name': self.CMMapBox.currentText(),
                'flip': self.ColorFlipCheckBox.isChecked(),
                'steps': self.StepSpin.value(),
                'steps_enabled': self.StepSpin.isEnabled(),
            }

#----------------------------------------------------------------------
    def _restoreClusterDisplayOptions(self):
        if not self._cluster_display_opts:
            return
        opts = self._cluster_display_opts
        self.CMCatBox.setCurrentIndex(opts['cat_index'])
        self.CMMapBox.setCurrentText(opts['map_name'])
        self.ColorFlipCheckBox.setChecked(opts['flip'])
        self.StepSpin.setValue(opts['steps'])
        self.StepSpin.setEnabled(opts['steps_enabled'])
        self._cluster_display_opts = None

#----------------------------------------------------------------------
    def _setClusterColorbarHandlesEnabled(self, enabled):
        # ColorBarItem creates a LinearRegionItem only in interactive mode.
        region = getattr(self.climgfig.bar, "region", None)
        if region is None:
            return
        region.setVisible(enabled)
        region.setMovable(enabled)

#----------------------------------------------------------------------
    def _setClusterColorbarTicks(self):
        axis = self.climgfig.bar.getAxis("right")
        ticks = [(i, str(i + 1)) for i in range(self.numclusters)]
        axis.setTicks([ticks])

#----------------------------------------------------------------------
    def _setColorbarAxisLabel(self, text, units=""):
        axis = self.climgfig.bar.getAxis("right")
        axis.setLabel(text=text, units=units)

#----------------------------------------------------------------------
    def _resetClusterColorbarTicks(self):
        axis = self.climgfig.bar.getAxis("right")
        axis.setTicks(None)

#----------------------------------------------------------------------
    def _reset_cluster_display_to_pca(self, current_component):
        # Leaving cluster mode on theta change: restore PCA plot/lut/colorbar defaults.
        self.com.cluster_calculated = 0
        self.cluster_display_mode = False
        self.selcluster = 0
        self.showallspectra = 0
        self.showerrormap = 0

        self.button_errormap.blockSignals(True)
        self.button_errormap.setChecked(False)
        self.button_errormap.blockSignals(False)
        self._setErrorMapButtonText(False)

        self.CMCatBox.setCurrentIndex(3)   # Diverging
        self.CMMapBox.setCurrentIndex(5)   # RdBu
        self.ColorFlipCheckBox.setChecked(False)
        self.StepSpin.setEnabled(True)
        self.StepSpin.setValue(256)
        self._setClusterColorbarHandlesEnabled(True)
        self._resetClusterColorbarTicks()

        self.slider_cl.blockSignals(True)
        self.slider_cl.setRange(0, max(self.numsigpca - 1, 0))
        self.slider_cl.setValue(current_component)
        self.slider_cl.blockSignals(False)
        self.slider_cl.setEnabled(True)

#----------------------------------------------------------------------
    def _refreshClusterSpectrumIfVisible(self):
        # In cluster mode, the right panel is fixed to cluster spectra.
        if self.cluster_display_mode and self.com.cluster_calculated == 1:
            self.showClusterSpectrum()
            return

        # Fallback for non-cluster mode behavior.
        if self.button_toggle_pca_plot.isChecked():
            self.showClusterSpectrum()

#----------------------------------------------------------------------
    def OnClusterColormapChanged(self, *_):
        self._refreshClusterSpectrumIfVisible()

#----------------------------------------------------------------------
    def OnPointClusterImage(self, event):
        if self.cluster_display_mode is False or self.com.cluster_calculated != 1:
            return
        if not self.climgfig.imageplot.sceneBoundingRect().contains(event.scenePos()):
            return

        vb = self.climgfig.imageplot.getViewBox()
        pos = vb.mapSceneToView(event.scenePos())
        x = int(np.floor(pos.x()))
        y = int(np.floor(pos.y()))
        if (x < 0) or (x >= self.stk.n_cols) or (y < 0) or (y >= self.stk.n_rows):
            return

        sel = int(self.anlz.cluster_indices[x, y])
        if (sel < 0) or (sel >= self.numclusters):
            return

        # Clicking a cluster should display its spectrum immediately.
        self.button_toggle_pca_plot.blockSignals(True)
        self.button_toggle_pca_plot.setChecked(True)
        self.button_toggle_pca_plot.blockSignals(False)
        self.UpdatePCAPlotView()
        if self.slider_cl.value() == sel:
            self.selcluster = sel
            if self.showerrormap:
                self.showClusterDistanceMap()
            else:
                self.showClusterImage()
            self._refreshClusterSpectrumIfVisible()
            return
        self.slider_cl.setValue(sel)

#----------------------------------------------------------------------
    def _getClusterColor(self, cluster_index):
        # Keep spectrum colors consistent with the active image LUT settings.
        try:
            map_name = self.CMMapBox.currentText()
            num_lut = max(2, int(self.StepSpin.value()))
            lut_range = (1.0, 0.0) if self.ColorFlipCheckBox.isChecked() else (0.0, 1.0)
            cm = pg.colormap.get(map_name, source='matplotlib')
            lut = cm.getLookupTable(*lut_range, num_lut)
            lut_index = int(np.clip(np.floor(((cluster_index + 0.5) / max(self.numclusters, 1)) * num_lut), 0, num_lut - 1))
            return QtGui.QColor(int(lut[lut_index, 0]), int(lut[lut_index, 1]), int(lut[lut_index, 2]))
        except Exception:
            return QtGui.QColor('blue')

#----------------------------------------------------------------------
    def showClusterSpectrum(self):
        """Draw the selected cluster spectrum using the same Set1 colour as the image."""
        self.evals_plot_item.clear()
        self.evals_plot_item.setLogMode(x=False, y=False)
        self.evals_plot_item.getAxis("bottom").setLabel(text="Photon Energy [eV]")

        if self.com.cluster_calculated != 1:
            return

        i = int(np.clip(self.selcluster, 0, self.numclusters - 1))

        prefix = self._caption_prefix()

        if self.showallspectra == 0:
            self.evals_plot_item.getAxis("left").setLabel(text="Optical Density")
            color = self._getClusterColor(i)
            pen = pg.mkPen(color=color, width=2)
            spec = self.anlz.clusterspectra[i]
            self.evals_plot_item.addItem(
                pg.PlotCurveItem(self.stk.ev, spec, pen=pen, name='Cluster {}'.format(i + 1))
            )
            title = "<center>{}Cluster {} spectrum</center>".format(prefix, i + 1)
        else:
            self.evals_plot_item.getAxis("left").setLabel(text="Normalized Optical Density")
            # Draw non-selected first, selected cluster last for visibility.
            draw_order = [j for j in range(self.numclusters) if j != i] + [i]
            normalized_specs = []
            for j in draw_order:
                spec_j = np.asarray(self.anlz.clusterspectra[j], dtype=float)
                # Vertically shift each spectrum so its maximum lands at +1.
                # This is a pure translation (spec_j - nanmax + 1), so the
                # curve shape and all amplitude differences are preserved exactly.
                # Works correctly for positive, negative, and mixed-sign spectra
                # without any scaling, sign flip, or shape distortion.
                if spec_j.size:
                    s_max = np.nanmax(spec_j)
                    if np.isfinite(s_max):
                        spec_j = spec_j - s_max + 1.0
                color = self._getClusterColor(j)
                self.evals_plot_item.addItem(
                    pg.PlotCurveItem(self.stk.ev, spec_j, pen=pg.mkPen(color=color, width=2),
                                     name='Cluster {}'.format(j + 1))
                )
                normalized_specs.append(spec_j)
            spec = np.asarray(normalized_specs, dtype=float)
            title = "<center>{}Normalized cluster spectra</center>".format(prefix)

        # Explicitly reset view range; otherwise the plot can keep the
        # principal-component x-range from the scree/eigenvalue view.
        try:
            ev = np.asarray(self.stk.ev, dtype=float)
            if ev.size > 1:
                x_min = float(np.nanmin(ev))
                x_max = float(np.nanmax(ev))
                if np.isfinite(x_min) and np.isfinite(x_max) and x_max > x_min:
                    self.evals_plot_item.setXRange(x_min, x_max, padding=0.02)
            y = np.asarray(spec, dtype=float)
            if y.size > 0:
                y_min = float(np.nanmin(y))
                y_max = float(np.nanmax(y))
                if np.isfinite(y_min) and np.isfinite(y_max):
                    y_span = max(y_max - y_min, 1e-9)
                    self.evals_plot_item.setYRange(y_min - 0.05 * y_span, y_max + 0.1 * y_span, padding=0)
        except Exception:
            # Keep plotting robust even if range computation fails.
            pass
        self.evals_plot_item.setTitle(title)

#----------------------------------------------------------------------
    def OnNClusterspin(self, value):
        self.init_nclusters = value

#----------------------------------------------------------------------
    def OnRemove1stpca(self, state):
        self.wo_1st_pca = 1 if state == QtCore.Qt.Checked else 0

#----------------------------------------------------------------------
    def OnSplitClusters(self, state):
        self.sigma_split = 1 if state == QtCore.Qt.Checked else 0

#----------------------------------------------------------------------
    def OnShowallspectra(self, state):
        self.showallspectra = 1 if state == QtCore.Qt.Checked else 0
        # Composite mode shows all clusters; selected-cluster scrolling is not meaningful.
        self.slider_cl.setEnabled(self.showallspectra == 0)
        # Left panel always shows composite, so no need to update it here.
        # Only update the right panel (spectrum view).
        self._refreshClusterSpectrumIfVisible()
        self.UpdatePCAPlotView()

#----------------------------------------------------------------------
    def OnShowScatterplots(self, event=None):
        scattplwin = Scatterplots(self.window(), self.com, self.anlz, page=self)
        scattplwin.show()

#----------------------------------------------------------------------
    def OnSaveCluster(self, event=None):
        savewin = SaveWinP3(self.window(), page=self)
        savewin.show()

#----------------------------------------------------------------------
    def OnSave(self, event):

        savewin = SaveWinP2(self)
        savewin.show()

#----------------------------------------------------------------------
    def Save(self, filename, path, spec_png=True, spec_pdf=False, spec_svg=False, spec_csv=False,
             img_png=True, img_pdf=False, img_svg=False, img_tif=False,
             evals_png=True, evals_pdf=False, evals_svg=False):
        """Save PCA results (images, spectra, eigenvalue plot) using live pyqtgraph widgets."""
        self.SaveFileName = os.path.join(path, filename)
        win = self.window()
        try:
            # Suppress all screen repaints for the duration of the export so
            # the UI does not flicker while we switch pages / iterate components.
            win.setUpdatesEnabled(False)
            saved_page = self.pca_plot_stack.currentWidget()

            # --- Eigenvalue plot — switch page so widget has full size ---
            def _save_evals(ext):
                fn = self.SaveFileName + '_PCAevals.' + ext
                if ext == 'svg':
                    exp = pg.exporters.SVGExporter(self.evals_plot_item)
                elif ext == 'pdf':
                    exp = PDFExporter(self.evals_plot_item)
                else:
                    exp = pg.exporters.ImageExporter(self.evals_plot_item)
                    exp.parameters()['width'] = int(self.pca_eigenvals_plotwidget.width())
                exp.export(fn)

            if evals_png or evals_pdf or evals_svg:
                self.pca_plot_stack.setCurrentWidget(self.pca_eigenvalues_page)
                QtWidgets.QApplication.processEvents()
                if evals_png:
                    _save_evals('png')
                if evals_pdf:
                    _save_evals('pdf')
                if evals_svg:
                    _save_evals('svg')

            # --- PCA images and spectra — iterate through components ---
            need_img  = img_png or img_pdf or img_svg or img_tif
            need_spec = spec_png or spec_pdf or spec_svg

            if need_img or need_spec or spec_csv:
                # Ensure the spectrum page is visible so specfig has proper pixel dimensions
                if need_spec:
                    self.pca_plot_stack.setCurrentWidget(self.pca_spectrum_page)
                    QtWidgets.QApplication.processEvents()

                saved_selpca = self.selpca
                for i in range(self.numsigpca):
                    self.selpca = i
                    image = self.anlz.pcaimages[:, :, i]
                    self.climgfig.draw(image,
                                       levels=(-self.anlz.pcaimagebounds[i],
                                               self.anlz.pcaimagebounds[i]),
                                       setpcalabel=True)
                    self.specfig.updatePlotDataOnPCA()
                    QtWidgets.QApplication.processEvents()

                    base_img  = self.SaveFileName + '_PCA_' + str(i + 1)
                    base_spec = self.SaveFileName + '_PCAspectrum_' + str(i + 1)

                    if img_png:
                        self.climgfig.SaveFig(base_img + '.png')
                    if img_pdf:
                        self.climgfig.SaveFig(base_img + '.pdf')
                    if img_svg:
                        self.climgfig.SaveFig(base_img + '.svg')
                    if img_tif:
                        img1 = Image.fromarray(np.rot90(self.anlz.pcaimages[:, :, i]))
                        img1.save(base_img + '.tif')

                    if spec_png:
                        self.specfig.SaveFig(base_spec + '.png')
                    if spec_pdf:
                        self.specfig.SaveFig(base_spec + '.pdf')
                    if spec_svg:
                        self.specfig.SaveFig(base_spec + '.svg')

                # Write all PCA spectra into a single multi-column CSV after the loop
                if spec_csv:
                    fn = self.SaveFileName + '_PCAspectra.csv'
                    with open(fn, 'w', newline='') as f:
                        writer = csv.writer(f, quoting=csv.QUOTE_MINIMAL)
                        for line in [
                            ['*********************  X-ray Absorption Data  ********************'],
                            ['*'],
                            ['* Formula: '],
                            ['* Common name: PCA Spectra'],
                            ['* Edge: '], ['* Acquisition mode: '],
                            ['* Source and purity: '], ['* Comments: '],
                            ['* Delta eV: '], ['* Min eV: '], ['* Max eV: '],
                            ['* Y axis: '], ['* Contact person: '],
                            ['* Write date: '], ['* Journal: '],
                            ['* Authors: '], ['* Title: '],
                            ['* Volume: '], ['* Issue number: '],
                            ['* Year: '], ['* Pages: '],
                            ['* Booktitle: '], ['* Editors: '],
                            ['* Publisher: '], ['* Address: '],
                            ['*--------------------------------------------------------------'],
                        ]:
                            writer.writerow(line)
                        col_names = ['photon energy'] + [
                            'Component {}'.format(i + 1) for i in range(self.numsigpca)]
                        writer.writerow(col_names)
                        for ie in range(len(self.stk.ev)):
                            row = ['{:06.6f}'.format(self.stk.ev[ie])]
                            for i in range(self.numsigpca):
                                row.append('{:09.6f}'.format(self.anlz.eigenvecs[ie, i]))
                            writer.writerow(row)

                # restore original component
                self.selpca = saved_selpca
                image = self.anlz.pcaimages[:, :, saved_selpca]
                self.climgfig.draw(image,
                                   levels=(-self.anlz.pcaimagebounds[saved_selpca],
                                           self.anlz.pcaimagebounds[saved_selpca]),
                                   setpcalabel=True)
                self.specfig.updatePlotDataOnPCA()

            # restore original stacked page
            self.pca_plot_stack.setCurrentWidget(saved_page)

        except IOError as e:
            QtWidgets.QMessageBox.warning(self, 'Error',
                                          'Could not save file: %s' % (e.strerror or e))
        finally:
            win.setUpdatesEnabled(True)
            win.update()

#----------------------------------------------------------------------
    def Save_CA(self, filename, path, spec_png=True, spec_pdf=False, spec_svg=False, spec_csv=False,
                img_png=True, img_pdf=False, img_svg=False, img_tif=False,
                indimgs_png=True, indimgs_pdf=False, indimgs_svg=False, indimgs_tif=False,
                scatt_png=True, scatt_pdf=False, scatt_svg=False):
        """Save Cluster Analysis results (composite/individual images, spectra, scatter plots)."""
        self.SaveFileName = os.path.join(path, filename)
        win = self.window()

        def _pg_export(item, widget, fn):
            ext = os.path.splitext(fn)[1][1:].lower()
            if ext == 'svg':
                exp = pg.exporters.SVGExporter(item)
            elif ext == 'pdf':
                exp = PDFExporter(item)
            else:
                exp = pg.exporters.ImageExporter(item)
                if widget is not None:
                    exp.parameters()['width'] = int(widget.width())
            exp.export(fn)

        try:
            win.setUpdatesEnabled(False)
            saved_page = self.pca_plot_stack.currentWidget()

            # --- Composite cluster image — live climgfig widget ---
            if img_png or img_pdf or img_svg:
                self.showClusterImage()
                QtWidgets.QApplication.processEvents()
                if img_png:
                    self.climgfig.SaveFig(self.SaveFileName + '_CAcimg.png')
                if img_pdf:
                    self.climgfig.SaveFig(self.SaveFileName + '_CAcimg.pdf')
                if img_svg:
                    self.climgfig.SaveFig(self.SaveFileName + '_CAcimg.svg')
            if img_tif:
                img1 = Image.fromarray(np.rot90(self.anlz.cluster_indices.astype(np.int32)))
                img1.save(self.SaveFileName + '_CAcimg.tif')

            # --- Individual cluster images — climgfig with per-cluster 2-colour LUT ---
            if indimgs_png or indimgs_pdf or indimgs_svg or indimgs_tif:
                for i in range(self.numclusters):
                    base = self.SaveFileName + '_CAimg_' + str(i + 1)
                    if indimgs_tif:
                        # Raw mask as float tif (no rendering)
                        indv = (self.anlz.cluster_indices == i).astype(np.float32)
                        Image.fromarray(np.rot90(indv)).save(base + '.tif')
                    if indimgs_png or indimgs_pdf or indimgs_svg:
                        qc = self._getClusterColor(i)
                        lut = np.array([[255, 255, 255],
                                        [qc.red(), qc.green(), qc.blue()]], dtype=np.uint8)
                        mask = (self.anlz.cluster_indices == i).astype(float)
                        self.climgfig.imageitem.setLookupTable(lut)
                        self.climgfig.draw(mask, setlut=False, setlabel=False, levels=(0, 1))
                        # Per-cluster export: hide colorbar, set cluster-number title
                        self.climgfig.bar.setVisible(False)
                        prefix = self._caption_prefix()
                        self.climgfig.imageplot.setTitle(
                            "<center>{}Cluster {:d}</center>".format(prefix, i + 1))
                        QtWidgets.QApplication.processEvents()
                        if indimgs_png:
                            self.climgfig.SaveFig(base + '.png')
                        if indimgs_pdf:
                            self.climgfig.SaveFig(base + '.pdf')
                        if indimgs_svg:
                            self.climgfig.SaveFig(base + '.svg')
                # Restore composite cluster image + LUT + colorbar visibility
                self.climgfig.bar.setVisible(True)
                self.showClusterImage()
                QtWidgets.QApplication.processEvents()

            # --- Cluster spectra — iterate through clusters via live evals_plot_item ---
            need_spec = spec_png or spec_pdf or spec_svg
            if need_spec or spec_csv:
                saved_selcluster = self.selcluster
                saved_showallspectra = self.showallspectra

                if need_spec:
                    self.pca_plot_stack.setCurrentWidget(self.pca_eigenvalues_page)
                    QtWidgets.QApplication.processEvents()
                    self.showallspectra = 0
                    for i in range(self.numclusters):
                        self.selcluster = i
                        self.showClusterSpectrum()
                        QtWidgets.QApplication.processEvents()
                        base_spec = self.SaveFileName + '_CAspectrum_' + str(i + 1)
                        if spec_png:
                            _pg_export(self.evals_plot_item, self.pca_eigenvals_plotwidget,
                                       base_spec + '.png')
                        if spec_pdf:
                            _pg_export(self.evals_plot_item, self.pca_eigenvals_plotwidget,
                                       base_spec + '.pdf')
                        if spec_svg:
                            _pg_export(self.evals_plot_item, self.pca_eigenvals_plotwidget,
                                       base_spec + '.svg')

                if spec_csv:
                    # Write all cluster spectra into a single CSV with one column per cluster.
                    fn = self.SaveFileName + '_CAspectra.csv'
                    with open(fn, 'w', newline='') as f:
                        writer = csv.writer(f, quoting=csv.QUOTE_MINIMAL)
                        # Standard MANTiS header block
                        for line in [
                            ['*********************  X-ray Absorption Data  ********************'],
                            ['*'],
                            ['* Formula: '],
                            ['* Common name: Cluster Analysis Spectra'],
                            ['* Edge: '], ['* Acquisition mode: '],
                            ['* Source and purity: '], ['* Comments: '],
                            ['* Delta eV: '], ['* Min eV: '], ['* Max eV: '],
                            ['* Y axis: '], ['* Contact person: '],
                            ['* Write date: '], ['* Journal: '],
                            ['* Authors: '], ['* Title: '],
                            ['* Volume: '], ['* Issue number: '],
                            ['* Year: '], ['* Pages: '],
                            ['* Booktitle: '], ['* Editors: '],
                            ['* Publisher: '], ['* Address: '],
                            ['*--------------------------------------------------------------'],
                        ]:
                            writer.writerow(line)
                        # Column headers
                        col_names = ['photon energy'] + [
                            'Cluster {}'.format(i + 1) for i in range(self.numclusters)]
                        writer.writerow(col_names)
                        # Data rows — one row per energy point
                        for ie in range(len(self.stk.ev)):
                            row = ['{:06.6f}'.format(self.stk.ev[ie])]
                            for i in range(self.numclusters):
                                row.append('{:09.6f}'.format(self.anlz.clusterspectra[i, ie]))
                            writer.writerow(row)

                # Restore cluster/spectrum state
                self.selcluster = saved_selcluster
                self.showallspectra = saved_showallspectra
                self.showallspectracb.blockSignals(True)
                self.showallspectracb.setChecked(saved_showallspectra == 1)
                self.showallspectracb.blockSignals(False)
                self.showClusterSpectrum()

            # --- Scatter plots — offscreen pg.GraphicsLayoutWidget ---
            if scatt_png or scatt_pdf or scatt_svg:
                od_reduced = self.anlz.pcaimages[:, :, 0:self.anlz.numsigpca]
                od_reduced = np.reshape(od_reduced,
                                        (self.stk.n_cols * self.stk.n_rows, self.anlz.numsigpca),
                                        order='F')
                clindices = np.reshape(self.anlz.cluster_indices,
                                       (self.stk.n_cols * self.stk.n_rows), order='F')

                ip_jp_pairs = [(ip, jp)
                               for ip in range(self.anlz.numsigpca)
                               for jp in range(self.anlz.numsigpca) if jp >= ip + 1]
                nplots = len(ip_jp_pairs)
                if nplots == 0:
                    pass
                else:
                    ncols_g = min(2, nplots)
                    nrows_g = max(int(np.ceil(nplots / ncols_g)), 1)

                    def _save_scatt(ext):
                        gl = pg.GraphicsLayoutWidget()
                        gl.setBackground('w')
                        gl.resize(600 if nplots > 1 else 300,
                                  max(250 * nrows_g, 250))
                        for k, (ip, jp) in enumerate(ip_jp_pairs):
                            row_g = k // ncols_g
                            col_g = k % ncols_g
                            plot = gl.addPlot(row=row_g, col=col_g)
                            plot.setLabel('bottom', 'Component ' + str(ip + 1))
                            plot.setLabel('left', 'Component ' + str(jp + 1))
                            for ci in range(self.numclusters):
                                thiscluster = np.where(clindices == ci)
                                x = od_reduced[thiscluster, ip].flatten()
                                y = od_reduced[thiscluster, jp].flatten()
                                qc = self._getClusterColor(ci)
                                scatter = pg.ScatterPlotItem(
                                    x=x, y=y, size=3, pen=pg.mkPen(None),
                                    brush=pg.mkBrush(qc.red(), qc.green(), qc.blue(), 140))
                                plot.addItem(scatter)
                        gl.show()
                        QtWidgets.QApplication.processEvents()
                        gl.hide()
                        fn = self.SaveFileName + '_CAscatterplot.' + ext
                        _pg_export(gl.scene(), None, fn)
                        if ext != 'svg':
                            # For ImageExporter, force a reasonable output width
                            pass
                        gl.close()

                    if scatt_png:
                        _save_scatt('png')
                    if scatt_pdf:
                        _save_scatt('pdf')
                    if scatt_svg:
                        _save_scatt('svg')

            # Restore original stacked page
            self.pca_plot_stack.setCurrentWidget(saved_page)

        except IOError as e:
            QtWidgets.QMessageBox.warning(self, 'Error',
                                          'Could not save file: %s' % (e.strerror or e))
        finally:
            win.setUpdatesEnabled(True)
            win.update()

#----------------------------------------------------------------------
    def reset_after_new_stack(self):
        """Clear PCA_new state so a newly loaded dataset starts from a clean view."""
        self.calcpca = False
        self.selpca = 0
        self.numsigpca = 2
        self.itheta = 0

        # Reset cluster state
        self.cluster_display_mode = False
        self.selcluster = 0
        self.numclusters = 0

        self.vartc.setText('0%')
        self.npcaspin.setValue(1)
        self.slider_theta.setVisible(False)
        self.slider_theta.setEnabled(False)
        self.ntc_clusters_found.setText('0')
        self.ntc_pcscaling_2.setValue(0.0)
        self.StepSpin.setEnabled(True)
        self._setClusterColorbarHandlesEnabled(True)
        self._resetClusterColorbarTicks()

        self.evals_plot_item.clear()
        self.spectrum_plotwidget.getPlotItem().clear()

        self.button_toggle_pca_plot.blockSignals(True)
        self.button_toggle_pca_plot.setChecked(False)
        self.button_toggle_pca_plot.setEnabled(False)
        self.button_toggle_pca_plot.blockSignals(False)
        self.button_movepcup.setEnabled(False)
        self.button_calcca.setEnabled(False)
        self.button_scatterplots.setEnabled(False)
        self.button_savecluster.setEnabled(False)
        self.showallspectracb.blockSignals(True)
        self.showallspectracb.setChecked(False)
        self.showallspectracb.blockSignals(False)
        self.showallspectracb.setEnabled(False)
        self.showallspectra = 0
        self.showerrormap = 0
        self.button_errormap.setEnabled(False)
        self.button_errormap.blockSignals(True)
        self.button_errormap.setChecked(False)
        self.button_errormap.blockSignals(False)
        self._setErrorMapButtonText(False)
        self._cluster_display_opts = None

        self.UpdatePCAPlotView()
