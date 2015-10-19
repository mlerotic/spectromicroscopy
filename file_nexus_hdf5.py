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
        regions = 0
        # Open HDF5 file
        f = h5py.File(filename, 'r')   
        
        #Check is exchange is in the file
        #Exchange HDF5 Group
        if 'entry1' in f:      
            nxs_format = True
            
        n_regions = 0
        for i in range(1,21):
            if 'entry{0}'.format(i) in f:
                n_regions += 1
            
        f.close()
            
        return nxs_format, n_regions

#---------------------------------------------------------------------- 
    def read_h5(self, filename, data_struct, loadregion):    
        
        # Open HDF5 file
        f = h5py.File(filename, 'r') 
        
        n_regions = 0
        for i in range(1,21):
            if 'entry{0}'.format(i) in f:
                n_regions += 1
                
               
        if (n_regions == 1) or (loadregion == 1):        
            if 'entry1' in f:
                entry1Grp = f['entry1']     
                   
                if 'Counter1' in entry1Grp:
                    counter1Grp = entry1Grp['Counter1']
                elif 'Counter0' in  entry1Grp:
                    counter1Grp = entry1Grp['Counter0']
                else:
                    print 'ERROR reading NEXUS HDF5 file.'
                    return
                    
                dsdata = counter1Grp['data']
                self.absdata = dsdata[...]
                self.absdata = np.transpose(self.absdata, axes=(2,1,0))

                
                if 'E' in counter1Grp:
                    dseng = counter1Grp['E']
                elif 'photon_energy' in counter1Grp:
                    dseng = counter1Grp['photon_energy']                           
                self.ev = dseng[...]
                
                #The image was rotated to show correct orientation so we have
                #to swap x and y distance    
                if 'sample_y' in counter1Grp:
                    dsx = counter1Grp['sample_y']        
                elif 'Y' in counter1Grp:
                    dsx = counter1Grp['Y']                                         
                self.y_dist = dsx[...]       
                
                if 'sample_x' in counter1Grp:
                    dsy = counter1Grp['sample_x']      
                elif 'X' in counter1Grp:
                    dsy = counter1Grp['X']                       
                self.x_dist = dsy[...]          
                
                dsdd = counter1Grp['count_time']                           
                self.data_dwell = dsdd[...]    
                      
            f.close()
            
            dims = np.int32(self.absdata.shape)
                   
            self.n_cols = dims[0]
            self.n_rows = dims[1]
            self.n_ev = dims[2]

        elif loadregion > 1:        
            if 'entry{0}'.format(loadregion) in f:
                entry1Grp = f['entry{0}'.format(loadregion)]     
                   
                if 'Counter1' in  entry1Grp:
                    counter1Grp = entry1Grp['Counter1']
                elif 'Counter0' in  entry1Grp:
                    counter1Grp = entry1Grp['Counter0']
                else:
                    print 'ERROR reading NEXUS HDF5 file.'
                    return
                    
                dsdata = counter1Grp['data']
                self.absdata = dsdata[...]
                self.absdata = np.transpose(self.absdata, axes=(2,1,0))
                
                dseng = counter1Grp['photon_energy']                           
                self.ev = dseng[...]
                
                #The image was rotated to show correct orientation so we have
                #to swap x and y distance    
                dsx = counter1Grp['sample_y']                           
                self.y_dist = dsx[...]       
                
                dsy = counter1Grp['sample_x']                           
                self.x_dist = dsy[...]          
                
                dsdd = counter1Grp['count_time']                           
                self.data_dwell = dsdd[...]    
                      
            f.close()
            
            dims = np.int32(self.absdata.shape)
                   
            self.n_cols = dims[0]
            self.n_rows = dims[1]
            self.n_ev = dims[2]
                        
            
        elif loadregion == 0: #Combined stack 
            
            #Multi region stack
            adata = []
            idimc = []
            idimr = []
            nc = 0
            nr = 0
            xd = []
            yd = []

            
            for i in range(1,n_regions+1):
                if 'entry{0}'.format(i) in f:      
                    entry1Grp = f['entry{0}'.format(i)]     
                       
                    if 'Counter1' in  entry1Grp:
                        counter1Grp = entry1Grp['Counter1']
                    elif 'Counter0' in  entry1Grp:
                        counter1Grp = entry1Grp['Counter0']
                    else:
                        print 'ERROR reading NEXUS HDF5 file.'
                        return
                        
                    dsdata = counter1Grp['data']
                    idata = dsdata[...]
                    idata = np.transpose(idata, axes=(2,1,0))
                    adata.append(idata)
                    
                    dims = np.int32(idata.shape)
               
                    inc = dims[0]
                    inr = dims[1]
                    self.n_ev = dims[2]
                    
                    if inc > nc:
                        nc = inc
                    nr += inr
                    idimc.append(inc)
                    idimr.append(inr)
                    
                    dseng = counter1Grp['photon_energy']                           
                    self.ev = dseng[...]
    
                    #The image was rotated to show correct orientation so we have
                    #to swap x and y distance  
                    dsy = counter1Grp['sample_y']                           
                    yd.append(dsy[...])      
                    
                    dsx = counter1Grp['sample_x']                           
                    xd.append(dsx[...])          
                    
                    dsdd = counter1Grp['count_time']                           
                    self.data_dwell = dsdd[...]   
                        
                       
            self.n_cols = nc
            self.n_rows = nr
            self.absdata = np.zeros((nc,nr,self.n_ev))
            self.x_dist = np.zeros(nc)
            self.y_dist = np.zeros(nr)
            thisr = 0
            for i in range(n_regions):
#                 self.absdata[0:idimc[i],thisr:idimr[i]+thisr, :] = adata[i]
#                 thisr+=idimr[i]
                self.absdata[0:idimc[i],thisr:idimr[i]+thisr, :] = adata[i]
                thisr+=idimr[i]
                
            thisr = 0
            for i in range(n_regions):            
                if nc == idimc[i]:
                    self.x_dist = xd[i]
                self.y_dist[thisr:idimr[i]+thisr]=yd[i]
                thisr+=idimr[i]
                
            
                          
        return
    
