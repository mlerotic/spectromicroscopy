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


from __future__ import division

import numpy as np
import scipy.interpolate
import os



title = 'STK'
extension = ['*.stk']
read_types = ['spectrum','image','stack']
write_types = []


def identify(filename):
        
    try:
        #Binary file - the only way to check if the file is ok to check the extension
        bn, extension = os.path.splitext(filename) 
        if extension.lower() == '.stk':
            return True
        else:
            return False
    except:
        return False


def read_old( filepath, data_stk):
    data_stk.read_h5(filepath)
    return


def GetFileStructure(FileName):
    return None


#-----------------------------------------------------------------------
def read(filename, self, selection=None, *args, **kwargs):
    f = open(str(filename),'rb')
    data = np.fromfile(f, np.int32, 3)
    data.byteswap(True)
    
    self.n_cols = data[0]
    self.n_rows = data[1]
    self.n_ev = data[2]
    
    #print 'self.n_cols, self.n_rows, self.n_ev', self.n_cols, self.n_rows, self.n_ev
    
    self.x_dist = np.fromfile(f, np.float32, self.n_cols)
    self.x_dist.byteswap(True)   
    
    
    self.y_dist = np.fromfile(f, np.float32, self.n_rows)
    self.y_dist.byteswap(True)     
    

    self.ev = np.fromfile(f, np.float32, self.n_ev)
    self.ev.byteswap(True)   
    
    msec = np.fromfile(f, np.float32, self.n_ev)
    msec.byteswap(True) 
     
    self.data_dwell = msec
    
    imagestack = np.fromfile(f, np.float32, self.n_cols*self.n_rows*self.n_ev)
    imagestack.byteswap(True)     
           
    self.absdata = np.empty((self.n_cols, self.n_rows, self.n_ev))
            
    self.absdata = np.reshape(imagestack, (self.n_cols, self.n_rows, self.n_ev), order='F')       

    f.close()
    
    self.fill_h5_struct_from_stk()
    
  
    return

#-----------------------------------------------------------------------
# Read x-ray absorption spectrum
def read_ascii(self, filename):

    spectrum_common_names = [' ']
    elist = []
    ilist = []
    names = []
    # Check the first character of the line and skip if not a number
    allowedchars = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '-', '.']
    with open(str(filename), 'r') as f:
        for line in f:
            if line.startswith('*'):
                if 'Common name' in line:
                    spectrum_common_names = [line.split(':')[-1].strip()]
            elif line[0] not in allowedchars:
                names = line.split(',')
                if names[0] == "photon energy":
                    names = names[1:]
                    #print(names)
                continue
            else:
                e, i = [x for x in line.split(',',1)]
                elist.append(e)
                ilist.append([float(a) for a in i.split(',')])
    if not elist:
        return
    spectrum_evdata = np.array([float(ev) for ev in elist])
    spectrum_data = np.transpose(np.array(ilist))
    if spectrum_evdata[-1] < spectrum_evdata[0]:
        spectrum_evdata = spectrum_evdata[::-1]
        spectrum_data = np.flip(spectrum_data,1)

    if spectrum_common_names[0] == ' ':
        spectrum_common_names = [os.path.splitext(os.path.basename(str(filename)))[0]]
    elif spectrum_common_names[0] == 'ROI spectrum' and len(names) > 1:
        spectrum_common_names = names
        #print(spectrum_common_name)

    return spectrum_evdata, spectrum_data, spectrum_common_names

#-----------------------------------------------------------------------
class x1astk:
    def __init__(self):
        pass
    
#-----------------------------------------------------------------------
    def read_stk(self, filename):
        f = open(str(filename),'rb')
        data = np.fromfile(f, np.int32, 3)
        data.byteswap(True)
        
        self.n_cols = data[0]
        self.n_rows = data[1]
        self.n_ev = data[2]
        
        #print 'self.n_cols, self.n_rows, self.n_ev', self.n_cols, self.n_rows, self.n_ev
        
        self.x_dist = np.fromfile(f, np.float32, self.n_cols)
        self.x_dist.byteswap(True)   
        
        
        self.y_dist = np.fromfile(f, np.float32, self.n_rows)
        self.y_dist.byteswap(True)     
        

        self.ev = np.fromfile(f, np.float32, self.n_ev)
        self.ev.byteswap(True)   
        
        msec = np.fromfile(f, np.float32, self.n_ev)
        msec.byteswap(True) 
         
        self.data_dwell = msec
        
        imagestack = np.fromfile(f, np.float32, self.n_cols*self.n_rows*self.n_ev)
        imagestack.byteswap(True)     
               
        self.absdata = np.empty((self.n_cols, self.n_rows, self.n_ev))
                
        self.absdata = np.reshape(imagestack, (self.n_cols, self.n_rows, self.n_ev), order='F')       

        f.close()
        
#         self.original_n_cols = self.n_cols.copy()
#         self.original_n_rows = self.n_rows.copy()
#         self.original_n_ev = self.n_ev.copy()
#         self.original_ev = self.ev.copy()
#         self.original_absdata = self.absdata.copy()
      
        return
