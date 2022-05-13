def resource_path(relative_path):
    import os, sys

    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(__file__))
    return os.path.join(base_path, relative_path)


def check_for_updates(current_version):
    import urllib.request, re
    from importlib.metadata import version as version_check
    from packaging.version import parse as parse_version
    # Scrape the version string from the PyPI:mantis-xray RSS feed
    try:
        with urllib.request.urlopen("https://pypi.org/rss/project/mantis_xray/releases.xml") as init_file:
            pypi_rss = init_file.read()
        #print(pypi_rss)
        pypi_list = re.findall(r"(?:\<title\>\s*)([\d\.]+)(?:\s*\</title\>)", pypi_rss.decode())[0]
        print("Latest package on PyPI is version {0}".format(pypi_list))
    except:
        pass
    
    # Scrape version string from the code in the github repository
    try:
        with urllib.request.urlopen("https://raw.githubusercontent.com/mlerotic/spectromicroscopy/master/mantis_xray/__init__.py") as init_file:
            github_init = init_file.read()
        #print(github_init)
        github_latest = re.search(r"(?:__version__*\s=*\s)['|\"]+([\d\.]+)", github_init.decode()).group(1)
        print("Current default (master) code is version {0}".format(github_latest))
    except:
        pass
    try:
        with urllib.request.urlopen(
                "https://raw.githubusercontent.com/mlerotic/spectromicroscopy/development/mantis_xray/__init__.py") as init_file:
            github_init = init_file.read()
        # print(github_init)
        github_latest = re.search(r"(?:__version__*\s=*\s)['|\"]+([\d\.]+)", github_init.decode()).group(1)
        print("Current development code is version {0}".format(github_latest))
    except:
        pass
    # PyQt5 & pyqtgraph version check
    if parse_version(version_check('pyqt5')) >= parse_version('5.15.6'):
        print("PyQt version in use is {0}".format(version_check("pyqt5")))
    else:
        print("PyQt version in use is {0}. Please consider updating to > 5.15.6 for full functionality.".format(version_check("pyqt5")))
    if parse_version(version_check('pyqtgraph')) >= parse_version('0.12.2'):
        print("PyQtGraph version in use is {0}".format(version_check("pyqtgraph")))
    else:
        print("PyQtGraph version in use is {0}. Please consider updating to > 0.12.2 for full functionality.".format(version_check("pyqtgraph")))



# PDF Exporter adopted from Orange https://orangedatamining.com/
# https://github.com/biolab/orange-widget-base/blob/master/orangewidget/utils/PDFExporter.py
# Potential integration into pyqtgraph in the future is under discussion: https://github.com/pyqtgraph/pyqtgraph/issues/1455#issuecomment-734299674
from pyqtgraph.exporters.Exporter import Exporter

from PyQt5 import QtCore
from PyQt5.QtWidgets import QGraphicsItem, QApplication
from PyQt5.QtGui import QPainter, QPdfWriter, QPageSize
from PyQt5.QtCore import QMarginsF, Qt, QSizeF, QRectF


class PDFExporter(Exporter):
    """A pdf exporter for pyqtgraph graphs. Based on pyqtgraph's
     ImageExporter.
     There is a bug in Qt<5.12 that makes Qt wrongly use a cosmetic pen
     (QTBUG-68537). Workaround: do not use completely opaque colors.
     There is also a bug in Qt<5.12 with bold fonts that then remain bold.
     To see it, save the OWNomogram output."""

    def __init__(self, item):
        Exporter.__init__(self, item)
        if isinstance(item, QGraphicsItem):
            scene = item.scene()
        else:
            scene = item
        bgbrush = scene.views()[0].backgroundBrush()
        bg = bgbrush.color()
        if bgbrush.style() == Qt.NoBrush:
            bg.setAlpha(0)
        self.background = bg

    def export(self, filename=None):
        pw = QPdfWriter(filename)
        dpi = int(QApplication.primaryScreen().logicalDotsPerInch())
        pw.setResolution(dpi)
        pw.setPageMargins(QMarginsF(0, 0, 0, 0))
        pw.setPageSize( # Tested with pyqt5-15.5.6
            QPageSize(QSizeF(self.getTargetRect().size()) / dpi * 25.4,
                      QPageSize.Millimeter))
        painter = QPainter(pw)
        try:
            self.setExportMode(True, {'antialias': True,
                                      'background': self.background,
                                      'painter': painter})
            painter.setRenderHint(QPainter.Antialiasing, True)
            if QtCore.QT_VERSION >= 0x050D00:
                painter.setRenderHint(QPainter.LosslessImageRendering, True)
            self.getScene().render(painter,
                                   QRectF(self.getTargetRect()),
                                   QRectF(self.getSourceRect()))
        finally:
            self.setExportMode(False)
        painter.end()
    
