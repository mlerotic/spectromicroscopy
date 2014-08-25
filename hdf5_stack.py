# 
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
# 
#   Copyright (C) 2011 Mirna Lerotic - 2nd Look
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
import scipy.ndimage
import h5py 
import datetime 


#----------------------------------------------------------------------
class h5data:
    def __init__(self):
        
        pass
        
#----------------------------------------------------------------------
    def check_h5_format(self, filename):
        
        nxs_format = False
        # Open HDF5 file
        f = h5py.File(filename, 'r')   
        
        #Check is exchange is in the file
        #Exchange HDF5 Group
        if 'entry1' in f:      
            nxs_format = True
            
        f.close()
            
        return nxs_format

#---------------------------------------------------------------------- 
    def read_h5(self, filename, data_struct):    
        
        # Open HDF5 file
        f = h5py.File(filename, 'r') 
        
        if 'entry1' in f:
            entry1Grp = f['entry1']     
               
            if 'Counter1' in  entry1Grp:
                counter1Grp = entry1Grp['Counter1']
                
                dsdata = counter1Grp['data']
                self.absdata = dsdata[...]
                self.absdata = np.transpose(self.absdata)
                
                dseng = counter1Grp['photon_energy']                           
                self.ev = dseng[...]

                dsx = counter1Grp['sample_x']                           
                self.x_dist = dsx[...]       
                
                dsy = counter1Grp['sample_y']                           
                self.y_dist = dsy[...]          
                
                dsdd = counter1Grp['count_time']                           
                self.data_dwell = dsdd[...]    
                  
        f.close()
        
        dims = np.int32(self.absdata.shape)
               
        self.n_cols = dims[0]
        self.n_rows = dims[1]
        self.n_ev = dims[2]
        
      
                
        return
    
