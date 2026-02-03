# coding: utf8
#
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
#
#   Copyright (C) 2013 Mirna Lerotic, 2nd Look
#   http://2ndlookconsulting.com
#   License: GNU GPL v3
#
#   Mantis is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   any later version.
#
#   Mantis is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details <http://www.gnu.org/licenses/>.


from __future__ import absolute_import

import sys
import os
from PyQt5 import QtWidgets

from .gui import MainFrame
from .helpers import resource_path

def main():
    app = QtWidgets.QApplication(sys.argv)
    
    # Load stylesheet if available
    qsspath = os.path.join(os.path.dirname(__file__), 'stylesheet_global.qss')
    if os.path.exists(qsspath):
        with open(qsspath, "r") as stylesheet:
            app.setStyleSheet(stylesheet.read())
            
    ex = MainFrame()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
