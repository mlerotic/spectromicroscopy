def resource_path(relative_path):
    import os, sys

    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(__file__))
    return os.path.join(base_path, relative_path)


def check_for_updates(current_version):
    import urllib.request, re
    from urllib.error import URLError
    timeout=.06
    message = 'Socket timed out. Check your internet connection.\nVersion checks skipped.'
    attempt = 0
    max_attempt = 3
    while attempt < max_attempt:
        # Scrape the version string from the PyPI:mantis-xray RSS feed
        try:
            with urllib.request.urlopen("https://pypi.org/rss/project/mantis_xray/releases.xml", timeout=timeout) as init_file:
                pypi_rss = init_file.read()
            #print(pypi_rss)
            pypi_list = re.findall(r"(?:\<title\>\s*)([\d\.]+)(?:\s*\</title\>)", pypi_rss.decode())[0]
            print("Latest package on PyPI is version {0}".format(pypi_list))
            attempt = 0
            break
        except URLError:
            attempt+=1
            print("Connection attempt " + str(attempt)+"/"+str(max_attempt)+" unsuccessful.")
            if attempt == 3:
                print(message)
                return
        except:
            pass
    
    # Scrape version string from the code in the github repository
    while attempt < max_attempt:
        try:
            with urllib.request.urlopen("https://raw.githubusercontent.com/mlerotic/spectromicroscopy/master/mantis_xray/__init__.py", timeout=timeout) as init_file:
                github_init = init_file.read()
            #print(github_init)
            github_latest = re.search(r"(?:__version__*\s=*\s)['|\"]+([\d\.]+)", github_init.decode()).group(1)
            print("Current default (master) code is version {0}".format(github_latest))
            attempt = 0
            break
        except URLError:
            attempt+=1
            print("Connection attempt " + str(attempt)+"/"+str(max_attempt)+" unsuccessful.")
            if attempt == 3:
                print(message)
                return
        except:
            pass
    # while attempt < max_attempt:
    #     try:
    #         with urllib.request.urlopen(
    #                 "https://raw.githubusercontent.com/mlerotic/spectromicroscopy/development/mantis_xray/__init__.py", timeout=timeout) as init_file:
    #             github_init = init_file.read()
    #         # print(github_init)
    #         github_latest = re.search(r"(?:__version__*\s=*\s)['|\"]+([\d\.]+)", github_init.decode()).group(1)
    #         print("Current development code is version {0}".format(github_latest))
    #         attempt = 0
    #         break
    #     except URLError:
    #         attempt+=1
    #         print("Connection attempt " + str(attempt)+"/"+str(max_attempt)+" unsuccessful.")
    #         if attempt == 3:
    #             print(message)
    #             return
    #     except:
    #         pass

def print_dependency_versions():
    from importlib.metadata import version as version_check, PackageNotFoundError
    from packaging.version import parse as parse_version
    import re
    print("Dependency versions:")
    for P in ['PyQt5>=5.15.9','numpy', 'scipy>=1.11.4', 'matplotlib>=3.6.0', 'h5py', 'Pillow', 'lxml', 'pyqtgraph>=0.13.7', "scikit-image>=0.19.1", "xdrlib3"]: #copy list from ../setup.py
        p = re.split('[><=]',P)
        # print(p[0], len(p))
        try:
            v_c = version_check(p[0].lower())
            print("\t{0} == {1}".format(p[0],v_c),end=' ')
            if len(p)>1:
                op_str = P[len(p[0]):-len(p[-1])]
                if   op_str == '<' :        check = parse_version(v_c) <  parse_version(p[-1])
                elif op_str == '>' :        check = parse_version(v_c) >  parse_version(p[-1])
                elif op_str in ['=','=='] : check = parse_version(v_c) == parse_version(p[-1])
                elif op_str == '!=':        check = parse_version(v_c) != parse_version(p[-1])
                elif op_str in ['>=','=>']: check = parse_version(v_c) >= parse_version(p[-1])
                elif op_str in ['<=','=>']: check = parse_version(v_c) <= parse_version(p[-1])
                else : check = True
                if not check:
                    print(" (upgrade to {0} for full functionality)".format(P[len(p[0]):]))
                else:
                    print("")
            else:
                print('')
        except PackageNotFoundError:
            pass

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
    
