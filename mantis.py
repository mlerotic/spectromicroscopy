# 
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
# 
#   Copyright (C) 2011 Mirna Lerotic, 2nd Look
#   http://2ndlook.co/products.html
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

import sys
import os
import getopt


import mantis_qt
import mantis_wx

""" ------------------------------------------------------------------------------------------------"""
                        
def main():
    
    verbose = True
    
    run_qt = 1
    run_wx = 0
    run_cl = 0
    arguments = sys.argv[1:]
    options, extraParams = getopt.getopt(arguments, '', ['wx', 'cl', 'nnma', 'ica', 'keyeng'])
    for opt, arg in options:
        if opt in '--wx':
            run_qt = 0
            run_wx = 1
        if opt in '--cl':
            run_qt = 0
            run_cl = 1            
            
    if run_qt == 1:
        m_qt = mantis_qt.main()
    elif run_wx == 1:
        m_wx = mantis_wx.main()
    elif run_cl == 1:
        wdir = 'D:/Mantis/x1adata/bacterium'
        wdir.strip('"') 
        if "'" in wdir: wdir = wdir[1:-1]
        
        wdir = os.path.normpath(wdir)
        
        if verbose: print 'working directory =', wdir
            
        if not os.path.exists(wdir):
            print 'Error - Directory ', wdir, ' does not exist. Please specify working directory.'
            return
    
    sys.exit()


if __name__ == '__main__':
    main()
    
