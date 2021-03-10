# 
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
# 
#   Copyright (C) 2015 Mirna Lerotic, 2nd Look
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
import os

          
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
def read_bim_info(filename):
    
    f = open(str(filename),'rb')
    data = np.fromfile(f, np.uint32, 6)

    n_cols = data[4]
    n_rows = data[5]
   
    f.close()

    fn = os.path.basename(str(filename))        
    fnlist = fn.split('_')
    ind = fnlist.index('eV')
    
    ev = float(fnlist[ind-1])


    return n_cols, n_rows, ev     
       
        
#----------------------------------------------------------------------    
def read_bim_list(self, filelist, filepath, ds): 
    
    #Fill the common stack data 
    file1 = os.path.join(filepath, filelist[0])   
    
    f = open(str(file1),'rb') 
    data = np.fromfile(f, np.uint32, 6)

    nmotpos = data[0]
    ndatatype = data[1]
    ndate = data[2]
    naxisnames = data[3]
                
    ncols = data[4]
    nrows = data[5]
    npix = ncols*nrows
    nev = len(filelist)
    absdata = np.zeros((ncols, nrows, nev))
    ev = np.zeros((nev))
    
    
    angles = np.fromfile(f, np.float64, 1)
    pixelsize = np.fromfile(f, np.float32, 1)
    
    x_dist = np.arange(np.float(ncols))*pixelsize
    y_dist = np.arange(np.float(nrows))*pixelsize
    
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
      
    dwell_msec = exposuretime
     
    f.close()
    
   
    
    #Read the image data        
    for i in range(len(filelist)):
        fn = filelist[i]
        filename = os.path.join(filepath, fn)   

        f = open(str(filename),'rb')
        data = np.fromfile(f, np.uint32, 6)

        angles = np.fromfile(f, np.float64, 1)
        pixelsize = np.fromfile(f, np.float32, 1)
        
        data = np.fromfile(f, np.uint32, 2)
        
        energy = np.fromfile(f, np.float64, 1)
        motpos = np.fromfile(f, np.float32, nmotpos)  
        axisnames = np.fromfile(f, np.uint8, naxisnames)
        exposuretime = np.fromfile(f, np.float32, 1)
        nimages = np.fromfile(f, np.uint32, 1)
        
        data = np.fromfile(f, np.uint8, ndatatype)
        data = np.fromfile(f, np.uint8, ndate)
        
        imagestack = np.fromfile(f, np.float32, npix)
        
        absdata[:,:,i] = np.reshape(imagestack, (ncols, nrows), order='C') 
        
        fn = os.path.basename(str(filename))        
        fnlist = fn.split('_')
        ind = fnlist.index('eV')
        
        
        ev[i] = float(fnlist[ind-1])

        f.close()            
        
        
    if ev[-1]<ev[0]:
        ev = ev[::-1]
        absdata = absdata[:,:, ::-1]
        
    
    #Fill the data structure with data: 
    ds.implements = 'information:exchange:spectromicroscopy'
    ds.version = '1.0'
    ds.information.comment = 'Converted from .bim file list',
    
    import datetime 
    now = datetime.datetime.now()
    ds.information.file_creation_datetime = now.strftime("%Y-%m-%dT%H:%M")

    ds.information.experimenter.name = ''
    ds.information.sample.name = ''
    
    ds.exchange.data = absdata
    ds.exchange.data_signal = 1
    ds.exchange.data_axes='x:y'
    
    ds.exchange.energy=ev
    ds.exchange.energy_units = 'ev'
    
    
    #Since we do not have a scanning microscope we fill the x_dist and y_dist from pixel_size
    x_dist = np.arange(np.float(ncols))*pixelsize
    y_dist = np.arange(np.float(nrows))*pixelsize

    ds.exchange.x = x_dist
    ds.exchange.x_units = 'um'
    ds.exchange.y = y_dist
    ds.exchange.y_units = 'um'     
    
    
    
    self.data_dwell = np.ones((nev))*exposuretime
    ds.spectromicroscopy.data_dwell = self.data_dwell
    
    return

