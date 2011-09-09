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


from __future__ import division
import os.path
import numpy as np
import scipy as scy
import netCDF4 

import data_struct


#----------------------------------------------------------------------
class sm:
    def __init__(self, ds):
        self.ds = ds

    
#----------------------------------------------------------------------    
    def read_sm_header(self, filename):
        
        try:
            f = netCDF4.Dataset(filename, 'r', format='NETCDF4')
        
        
#        print f.file_format
#        print f.groups
#        print f.variables
        #print f.dimensions
#        
#        for dimname, dimobj in f.dimensions.items():
#            print dimname, len(dimobj), dimobj.isunlimited()
            
            ncols = len(f.dimensions['n_fast_pixels'])
            nrows = len(f.dimensions['n_slow_pixels'])
            ev = getattr(f,'start_ev')
            x_center_um = getattr(f,'x_center_um')
            y_center_um = getattr(f,'y_center_um')
            
#        
#        for name in f.ncattrs():
#            print 'Global attr', name, '=', getattr(f,name)
            
        #print netCDF4.chartostring(getattr(f,'detector_name'))
        
        #n_cols = getattr(f, )
  
            f.close()
        
            return 1, ncols, nrows, ev
        except:
            return -1, 0,0,0
        
#----------------------------------------------------------------------    
    def read_sm_list(self, filelist, filepath, ds): 
        
        #Fill the common stack data 
        file1 = os.path.join(filepath, filelist[0])   
        fcdf = netCDF4.Dataset(file1, 'r', format='NETCDF4')   

            
        ncols = len(fcdf.dimensions['n_fast_pixels'])
        nrows = len(fcdf.dimensions['n_slow_pixels'])
        npix = ncols*nrows
        nev = len(filelist)
        absdata = np.zeros((ncols, nrows, nev))
        ev = np.zeros((nev))
        
        fdevicen = netCDF4.chartostring(getattr(fcdf,'fast_device_name'))
        sdevicen = netCDF4.chartostring(getattr(fcdf,'slow_device_name'))
        
        x_center_um = -1
        y_center_um = -1
        if fdevicen == 'XPZT':
            x_center_um = getattr(fcdf,'x_center_um')
        if sdevicen == 'YPZT':
            y_center_um = getattr(fcdf,'y_center_um')   
            
        if (x_center_um == -1) or (y_center_um == -1):
            print 'Error in stack build: Not piezo scan'
            return
        
        x_pixel_um = getattr(fcdf,'fast_pixel_um')
        y_pixel_um = getattr(fcdf,'slow_pixel_um')
        
        x_direction = getattr(fcdf,'fast_direction')
        y_direction = getattr(fcdf,'slow_direction')
                
        x_dist = (np.arange(npix)-np.fix(npix/2.0)) * x_direction * x_pixel_um + x_center_um
        y_dist = (np.arange(npix)-np.fix(npix/2.0)) * y_direction * y_pixel_um + y_center_um         
        
        
        self.ds.implements = 'base, exchange, spectromicroscopy'
        self.ds.version = '0.9'
        self.ds.file_creation_datetime = netCDF4.chartostring(getattr(fcdf,'systime'))
        self.ds.comment = netCDF4.chartostring(getattr(fcdf,'comments1'))
        
        self.ds.add_experimenter(name = netCDF4.chartostring(getattr(fcdf,'scientist')))

        self.ds.add_sample(name = netCDF4.chartostring(getattr(fcdf,'sample')))
        
        self.ds.add_exchange(collection_datetime = netCDF4.chartostring(getattr(fcdf,'systime')))
       
        
        fcdf.close()
        
        
        
        #Read the image data        
        for i in range(len(filelist)):
            fn = filelist[i]
            filename = os.path.join(filepath, fn)   
            fcdf = netCDF4.Dataset(filename, 'r', format='NETCDF4')   
            
            image = fcdf.variables['image_data'][:]
            absdata[:,:,i] = image[:,:,0]
            eng = getattr(fcdf,'start_ev')
            ev[i] = eng
                
            fcdf.close()
            
        if ev[-1]<ev[0]:
            ev = ev[::-1]
            absdata = absdata[:,:, ::-1]
            
        
        self.ds.exchange.add_detector(data=absdata,
                                      signal = 1, 
                                      axes='x:y', 
                                      energy=ev, 
                                      energy_units = 'ev')
        
        self.ds.exchange.detector[0].add_dimscale(key = 'x', units = 'um', data = x_dist)
        self.ds.exchange.detector[0].add_dimscale(key = 'y', units = 'um', data = y_dist)
        
        return