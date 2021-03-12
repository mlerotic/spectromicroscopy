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
import os.path
import numpy as np
import scipy as scy
import netCDF4 

from mantis_xray import data_struct


#----------------------------------------------------------------------    
def read_sm_header(filename):
    
    try:
        print (filename)
        f = netCDF4.Dataset(filename, 'r', format='NETCDF4')

        
        ncols = len(f.dimensions['n_fast_pixels'])
        nrows = len(f.dimensions['n_slow_pixels'])
        ev = getattr(f,'start_ev')
        x_center_um = getattr(f,'x_center_um')
        y_center_um = getattr(f,'y_center_um')


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
    
    #fdevicen = netCDF4.chartostring(getattr(fcdf,'fast_device_name'))
    #sdevicen = netCDF4.chartostring(getattr(fcdf,'slow_device_name'))
    
    fdevicenl = getattr(fcdf,'fast_device_name')
    sdevicenl = getattr(fcdf,'slow_device_name')
    
    fdevicen = chr(fdevicenl[0]) + chr(fdevicenl[1]) + chr(fdevicenl[2]) + chr(fdevicenl[3])
    sdevicen = chr(sdevicenl[0]) + chr(sdevicenl[1]) + chr(sdevicenl[2]) + chr(sdevicenl[3])
    
    x_center_um = -1
    y_center_um = -1
    if fdevicen == 'XPZT':
        x_center_um = getattr(fcdf,'x_center_um')
    if sdevicen == 'YPZT':
        y_center_um = getattr(fcdf,'y_center_um')   
        
    if (x_center_um == -1) or (y_center_um == -1):
        print ('Error in stack build: Not piezo scan')
        return
    
    x_pixel_um = getattr(fcdf,'fast_pixel_um')
    y_pixel_um = getattr(fcdf,'slow_pixel_um')
    
    x_direction = getattr(fcdf,'fast_direction')
    y_direction = getattr(fcdf,'slow_direction')
    
    dwell_msec = getattr(fcdf,'dwell_msec')
            
    x_dist = (np.arange(npix)-np.fix(npix/2.0)) * x_direction * x_pixel_um + x_center_um
    y_dist = (np.arange(npix)-np.fix(npix/2.0)) * y_direction * y_pixel_um + y_center_um         
    
    
    ds.implements = 'information:exchange:spectromicroscopy'
    ds.version = '1.0'
    ds.information.file_creation_datetime = ' '#netCDF4.chartostring(getattr(fcdf,'systime'))
    ds.information.comment = ' '#netCDF4.chartostring(getattr(fcdf,'comments1'))
    
    ds.information.experimenter.name = ' '#netCDF4.chartostring(getattr(fcdf,'scientist'))

    ds.information.sample.name = ' '#netCDF4.chartostring(getattr(fcdf,'sample'))
    
    
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
        
    
    ds.exchange.data = absdata
    ds.exchange.data_signal = 1
    ds.exchange.data_axes='x:y'
    
    ds.exchange.energy=ev
    ds.exchange.energy_units = 'ev'
    
    
    ds.exchange.x = x_dist
    ds.exchange.x_units = 'um'
    ds.exchange.y = y_dist
    ds.exchange.y_units = 'um'        
    ds.spectromicroscopy.data_dwell = dwell_msec
    
    return

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
        
        #fdevicen = netCDF4.chartostring(getattr(fcdf,'fast_device_name'))
        #sdevicen = netCDF4.chartostring(getattr(fcdf,'slow_device_name'))
        
        fdevicenl = getattr(fcdf,'fast_device_name')
        sdevicenl = getattr(fcdf,'slow_device_name')
        
        fdevicen = chr(fdevicenl[0]) + chr(fdevicenl[1]) + chr(fdevicenl[2]) + chr(fdevicenl[3])
        sdevicen = chr(sdevicenl[0]) + chr(sdevicenl[1]) + chr(sdevicenl[2]) + chr(sdevicenl[3])
        
        x_center_um = -1
        y_center_um = -1
        if fdevicen == 'XPZT':
            x_center_um = getattr(fcdf,'x_center_um')
        if sdevicen == 'YPZT':
            y_center_um = getattr(fcdf,'y_center_um')   
            
        if (x_center_um == -1) or (y_center_um == -1):
            print ('Error in stack build: Not piezo scan')
            return
        
        x_pixel_um = getattr(fcdf,'fast_pixel_um')
        y_pixel_um = getattr(fcdf,'slow_pixel_um')
        
        x_direction = getattr(fcdf,'fast_direction')
        y_direction = getattr(fcdf,'slow_direction')
        
        dwell_msec = getattr(fcdf,'dwell_msec')
                
        x_dist = (np.arange(npix)-np.fix(npix/2.0)) * x_direction * x_pixel_um + x_center_um
        y_dist = (np.arange(npix)-np.fix(npix/2.0)) * y_direction * y_pixel_um + y_center_um         
        
        
        self.ds.implements = 'information:exchange:spectromicroscopy'
        self.ds.version = '1.0'
        self.ds.information.file_creation_datetime = ' '#netCDF4.chartostring(getattr(fcdf,'systime'))
        self.ds.information.comment = ' '#netCDF4.chartostring(getattr(fcdf,'comments1'))
        
        self.ds.information.experimenter.name = ' '#netCDF4.chartostring(getattr(fcdf,'scientist'))

        self.ds.information.sample.name = ' '#netCDF4.chartostring(getattr(fcdf,'sample'))
        
        
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
            
        
        self.ds.exchange.data = absdata
        self.ds.exchange.data_signal = 1
        self.ds.exchange.data_axes='x:y'
        
        self.ds.exchange.energy=ev
        self.ds.exchange.energy_units = 'ev'
        
        
        self.ds.exchange.x = x_dist
        self.ds.exchange.x_units = 'um'
        self.ds.exchange.y = y_dist
        self.ds.exchange.y_units = 'um'        
        self.ds.spectromicroscopy.data_dwell = dwell_msec
        
        return
