from __future__ import print_function
#
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
#
#   Copyright (C) 2011 Mirna Lerotic, 2nd Look
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

import numpy as np
import os


title = 'text spectrum'
extension = ['*.csv','*.txt','*.xas']
read_types = ['spectrum']
write_types = ['spectrum']


def identify(filename):
    try:
        with open(file_name, 'tr') as check_file:  # try open file in text mode
            check_file.read()
            return True
    except:
        return False

def GetFileStructure(FileName):
    return None


#----------------------------------------------------------------------
def read(filename, self, selection=None, *args, **kwargs):

    f = open(str(filename),'r')
    
    elist = []
    ilist = []    

    for line in f:
        if line.startswith(("*","%","#")):
            pass
        else:
            e, i = [float (x) for x in line.split(',')] 
            elist.append(e)
            ilist.append(i)
            
    self.evi0 = np.array(elist)
    self.i0data = np.array(ilist) 
            
    f.close()
    
    self.i0_dwell = None
       
    return
    
#----------------------------------------------------------------------

def write(filename, data, energies, title=None ):
    if title is None:
        title = 'Spectrum'
    
    with open(filename, 'w') as f:
        print('********************* '+title+'  ********************', file=f)
        print('*', file=f)
        print('* ev, intensity', file=f)
        for ie in range(len(energies)):
            print('{0:06.2f}, {1:06f}'.format(energies[ie], data[ie]), file=f)



