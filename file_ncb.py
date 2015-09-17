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
import os

#----------------------------------------------------------------------
class Cncb:
    def __init__(self):
        pass
    
   
           
#----------------------------------------------------------------------
    def read_ncb(self, filename):
        
        print 'Reading stack:', filename
        
        basename, extension = os.path.splitext(filename) 
        dat_fn = basename + '.dat'
        

        f = open(str(dat_fn),'r')
        lines = f.readlines()
        temp = lines[0].split()
        
        self.n_cols = int(temp[0])
        self.n_rows = int(temp[1])
        scale = float(temp[2])
    
        #print self.n_cols, self.n_rows, scale
        
        temp = lines[1].split()
        x_start = float(temp[0])
        x_stop = float(temp[1])
        
        temp = lines[2].split()
        y_start = float(temp[0])
        y_stop = float(temp[1])

        self.n_ev = int(lines[3])

        #print self.n_ev
        
    
        self.ev = np.zeros((self.n_ev))
        self.data_dwell = np.zeros((self.n_ev))
        filename_list = []
        for i in range(self.n_ev):
            self.ev[i] = float(lines[4+i])
            
        #print self.ev

        for i in range(self.n_ev):
            self.ev[i] = float(lines[4+i])

        for i in range(self.n_ev):
            temp = lines[4+self.n_ev+i].split()
            
            filename_list.append(temp[0])
            self.data_dwell[i] = float(temp[2])

    
        f.close()

        
        
        
        if scale < 0 :
            dataformat = np.float32
        else:
            dataformat = np.int16
            
        f = open(str(filename),'rb')
        big_array = np.fromfile(f, dataformat, self.n_cols*self.n_rows*self.n_ev)
        f.close()    
        
        
        if scale < 0 and scale != 1 :
            image_stack = big_array.astype(np.float)/scale
            print "data rescaled by ", 1./scale
        else:
            image_stack = big_array.astype(np.float)
            

        if (x_start > x_stop) :
            image_stack = image_stack[::-1,:,:]
            t = x_start 
            x_start = x_stop 
            x_stop = t
            print 'x data reversed'

        if (y_start > y_stop) :
            image_stack = image_stack[:,::-1,:]
            t = y_start 
            y_start = y_stop 
            y_stop = t
            print 'y data reversed'

            
        xstep = (x_stop-x_start)/(self.n_cols-1)
        self.y_dist = np.arange(x_start,x_stop+xstep, xstep)
        
        ystep = (y_stop-y_start)/(self.n_rows-1)
        self.x_dist = np.arange(y_start,y_stop+ystep, ystep)              

   
        self.absdata = np.empty((self.n_cols, self.n_rows, self.n_ev))
                
        self.absdata = np.reshape(image_stack, (self.n_cols, self.n_rows, self.n_ev), order='F')    
        
        #Transpose the array to show correct orientation
        self.absdata = np.transpose(self.absdata, axes=(1,0,2))
        t = self.n_cols
        self.n_cols = self.n_rows
        self.n_rows = t
        

           
        
        print self.n_cols, self.n_rows
        print self.absdata.shape
  
        return
