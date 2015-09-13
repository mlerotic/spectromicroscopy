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

#----------------------------------------------------------------------
class x1astk:
    def __init__(self):
        pass
    
#----------------------------------------------------------------------
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

#----------------------------------------------------------------------   
    def read_stk_i0_xas(self, filename):

        f = open(str(filename),'r')
        
        elist = []
        ilist = []    
    
        for line in f:
            if line.startswith("*"):
                pass
            else:
                e, i = [float (x) for x in line.split()] 
                elist.append(e)
                ilist.append(i)
                
        self.evi0 = np.array(elist)
        self.i0data = np.array(ilist) 
                
        f.close()
        
        self.i0_dwell = None
           
    
        return
    

#----------------------------------------------------------------------   
    def read_stk_i0_csv(self, filename):

        f = open(str(filename),'r')
        
        elist = []
        ilist = []    
    
        for line in f:
            if line.startswith("*"):
                pass
            else:
                e, i = [float (x) for x in line.split(',')] 
                elist.append(e)
                ilist.append(i)
                
        self.evi0 = np.array(elist)
        self.i0data = np.array(ilist) 
                
        f.close()
        
        self.i0_dwell = None
           
           
#----------------------------------------------------------------------
    def read_bim(self, filename):
        
        f = open(str(filename),'rb')
        data = np.fromfile(f, np.uint32, 6)

        nmotpos = data[0]
        ndatatype = data[1]
        ndate = data[2]
        naxisnames = data[3]

        self.n_cols = data[4]
        self.n_rows = data[5]
        self.n_ev = 1
         

        angles = np.fromfile(f, np.float64, 1)
        pixelsize = np.fromfile(f, np.float32, 1)
        
        self.x_dist = np.arange(np.float(self.n_cols))*pixelsize
        self.y_dist = np.arange(np.float(self.n_rows))*pixelsize
        
        data = np.fromfile(f, np.uint32, 2)
        
        hbin = data[0]
        vbin = data[1]
        
        energy = np.fromfile(f, np.float64, 1)
        motpos = np.fromfile(f, np.float32, nmotpos)  
        axisnames = np.fromfile(f, np.uint8, naxisnames)
        exposuretime = np.fromfile(f, np.float32, 1)
        nimages = np.fromfile(f, np.uint32, 1)
        
        data = np.fromfile(f, np.uint8, ndatatype)
        data = np.fromfile(f, np.uint8, ndate)
        
        npix = self.n_cols*self.n_rows
        
        imagestack = np.fromfile(f, np.float32, npix)
        
        self.absdata = np.reshape(imagestack, (self.n_cols, self.n_rows, self.n_ev), order='C') 
        
        
        fn = os.path.basename(str(filename))        
        fnlist = fn.split('_')
        ind = fnlist.index('eV')
        
        self.ev = [float(fnlist[ind-1])] 
          
        self.data_dwell = np.zeros((self.n_ev))+exposuretime
         

        f.close()
        
    
        return            
        
        
        
        
#----------------------------------------------------------------------
    def read_bim_info(self, filename):
        
        f = open(str(filename),'rb')
        data = np.fromfile(f, np.uint32, 6)

        n_cols = data[4]
        n_rows = data[5]
       
        f.close()
   
        fn = os.path.basename(str(filename))        
        fnlist = fn.split('_')
        ind = fnlist.index('eV')
        
        ev = [float(fnlist[ind-1])] 

 
        return n_cols, n_rows, ev     
        
        
       
            
    