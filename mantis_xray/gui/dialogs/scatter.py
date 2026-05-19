import numpy as np
import pyqtgraph as pg
from PyQt5 import QtWidgets, QtCore, QtGui


class Scatterplots(QtWidgets.QDialog):
    """Scatter-plot dialog showing PCA scores coloured by cluster assignment.

    Parameters
    ----------
    parent : QWidget
        Main window (passed through to QDialog).
    common : common
        Shared application state.
    analz : analysis object
        Provides ``pcaimages``, ``cluster_indices``, ``numsigpca``.
    page : PagePCACluster, optional
        Reference to the calling page; used to obtain per-cluster colours via
        ``page._getClusterColor(i)``.
    """

    def __init__(self, parent, common, analz, page=None):
        super().__init__(parent)

        self.com = common
        self.anlz = analz
        self.page = page

        self.numsigpca = max(self.anlz.numsigpca, 1)
        self.ncols = self.anlz.stack.n_cols
        self.nrows = self.anlz.stack.n_rows

        # Flatten PCA score images: shape (npx, numsigpca)
        self.od_reduced = np.reshape(
           self.anlz.pcaimages[:, :, :self.numsigpca],
           (self.ncols * self.nrows, self.numsigpca),
           order='F',
        )
        self.clindices = np.reshape(
           self.anlz.cluster_indices, (self.ncols * self.nrows), order='F'
        )
        self.numclusters = (
           self.page.numclusters if page is not None else int(self.clindices.max()) + 1
        )

        # 0-based component indices
        self.pca_x = 0
        self.pca_y = 1 if self.numsigpca > 1 else 0

        self._build_ui()
        self._update_plot()

    # ------------------------------------------------------------------
    def _build_ui(self):
        self.setWindowTitle('Scatter plots')
        self.resize(580, 540)

        # ── plot widget ──────────────────────────────────────────────
        self.plot_widget = pg.PlotWidget(background='w')
        self.plot_widget.showGrid(x=True, y=True, alpha=0.3)
        self.plot_item = self.plot_widget.getPlotItem()
        self.plot_item.showAxis('top',   show=True)
        self.plot_item.showAxis('right', show=True)
        self.plot_item.getAxis('top').setStyle(showValues=False, tickLength=0)
        self.plot_item.getAxis('right').setStyle(showValues=False, tickLength=0)

        # ── x slider (horizontal, below plot) ────────────────────────
        self.label_x = QtWidgets.QLabel('X: PC 1')
        self.label_x.setAlignment(QtCore.Qt.AlignCenter)

        self.slider_x = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.slider_x.setRange(1, self.numsigpca)
        self.slider_x.setValue(self.pca_x + 1)
        self.slider_x.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.slider_x.setTickInterval(1)
        self.slider_x.valueChanged[int].connect(self._on_slider_x)

        # ── y slider (vertical, left of plot) ────────────────────────
        self._yl_proxy = _RotatedLabel('Y: PC 2')

        self.slider_y = QtWidgets.QSlider(QtCore.Qt.Vertical)
        self.slider_y.setRange(1, self.numsigpca)
        self.slider_y.setValue(self.pca_y + 1)
        self.slider_y.setTickPosition(QtWidgets.QSlider.TicksLeft)
        self.slider_y.setTickInterval(1)
        self.slider_y.valueChanged[int].connect(self._on_slider_y)

        # ── assemble layout ──────────────────────────────────────────
        # Center area: [y-slider | plot]
        center = QtWidgets.QHBoxLayout()
        center.setSpacing(4)

        left_col = QtWidgets.QVBoxLayout()
        left_col.addWidget(self._yl_proxy, alignment=QtCore.Qt.AlignHCenter)
        left_col.addWidget(self.slider_y)
        center.addLayout(left_col)
        center.addWidget(self.plot_widget, stretch=1)

        # Bottom row: [x-slider | label]
        bottom = QtWidgets.QHBoxLayout()
        bottom.addWidget(self.slider_x, stretch=1)
        bottom.addWidget(self.label_x)

        main = QtWidgets.QVBoxLayout(self)
        main.addLayout(center, stretch=1)
        main.addLayout(bottom)

        esc_shortcut = QtWidgets.QShortcut(QtGui.QKeySequence("Esc"), self)
        esc_shortcut.activated.connect(self.reject)

    # ------------------------------------------------------------------
    def _on_slider_x(self, value):
        self.pca_x = value - 1
        self.label_x.setText('X: PC {}'.format(value))
        self._update_plot()

    def _on_slider_y(self, value):
        self.pca_y = value - 1
        self._yl_proxy.setText('Y: PC {}'.format(value))
        self._update_plot()

    # ------------------------------------------------------------------
    def _get_color(self, cluster_index):
        """Return a (R, G, B) tuple for the given cluster index."""
        if self.page is not None:
           try:
               qc = self.page._getClusterColor(cluster_index)
               return (qc.red(), qc.green(), qc.blue())
           except Exception:
               pass
        # Fallback: tab10-like palette
        _fallback = [
           (31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40),
           (148, 103, 189), (140, 86, 75), (227, 119, 194), (127, 127, 127),
           (188, 189, 34), (23, 190, 207),
        ]
        return _fallback[cluster_index % len(_fallback)]

    # ------------------------------------------------------------------
    def _update_plot(self):
        self.plot_item.clear()

        x_all = self.od_reduced[:, self.pca_x]
        y_all = self.od_reduced[:, self.pca_y]

        for i in range(self.numclusters):
           mask = self.clindices == i
           if not np.any(mask):
               continue
           r, g, b = self._get_color(i)
           brush = pg.mkBrush(r, g, b, 140)
           scatter = pg.ScatterPlotItem(
               x=x_all[mask], y=y_all[mask],
               pen=pg.mkPen(None), brush=brush, size=5,
           )
           self.plot_item.addItem(scatter)

        self.plot_item.setLabel('bottom', 'PC {}'.format(self.pca_x + 1))
        self.plot_item.setLabel('left',   'PC {}'.format(self.pca_y + 1))


# ---------------------------------------------------------------------------
class _RotatedLabel(QtWidgets.QWidget):
    """Compact widget that renders a text label rotated 90° counter-clockwise."""
    def __init__(self, text='', parent=None):
        super().__init__(parent)
        self._text = text
        self.setFixedWidth(20)

    def setText(self, text):
        self._text = text
        self.update()

    def sizeHint(self):
        fm = QtGui.QFontMetrics(self.font())
        return QtCore.QSize(20, fm.horizontalAdvance(self._text) + 8)

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        painter.translate(self.width() / 2, self.height())
        painter.rotate(-90)
        fm = QtGui.QFontMetrics(painter.font())
        painter.drawText(-fm.horizontalAdvance(self._text) // 2, fm.ascent() // 2, self._text)
        painter.end()
