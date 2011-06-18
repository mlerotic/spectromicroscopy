'''
Created on Oct 25, 2010

@author: Mirna Lerotic
'''

import numpy as np
import scipy.interpolate

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
        
        #print self.n_cols, self.n_rows, self.n_ev
        
        x_dist = np.fromfile(f, np.float32, self.n_cols)
        x_dist.byteswap(True)   
        
        y_dist = np.fromfile(f, np.float32, self.n_rows)
        y_dist.byteswap(True)     

        self.ev = np.fromfile(f, np.float32, self.n_ev)
        self.ev.byteswap(True)   
        
        msec = np.fromfile(f, np.float32, self.n_ev)
        msec.byteswap(True)   
        
        self.imagestack = np.fromfile(f, np.float32, self.n_cols*self.n_rows*self.n_ev)
        self.imagestack.byteswap(True)     
               
        self.absdata = np.empty((self.n_cols, self.n_rows, self.n_ev))
                
        self.absdata = np.reshape(self.imagestack, (self.n_cols, self.n_rows, self.n_ev), order='F')       

        f.close()
        
        self.original_n_cols = self.n_cols.copy()
        self.original_n_rows = self.n_rows.copy()
        self.original_n_ev = self.n_ev.copy()
        self.original_ev = self.ev.copy()
        self.original_absdata = self.absdata.copy()
      
        return

#----------------------------------------------------------------------   
    def read_stk_i0(self, filename):
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
    
        self.calculate_optical_density()    
        
    
        return
    
#----------------------------------------------------------------------   
    def calc_histogram(self):
        #calculate average flux for each pixel
        self.averageflux = np.mean(self.absdata,axis=2)
        
        self.histogram = self.averageflux
        
        
#----------------------------------------------------------------------   
    def i0_from_histogram(self, fluxmin, fluxmax):

        self.evi0hist = self.ev.copy()

        self.i0datahist = np.zeros(self.n_ev)

        i0_indices = np.where((fluxmin<self.averageflux)&(self.averageflux<fluxmax))
        if i0_indices:
            invnumel = 1./self.averageflux[i0_indices].shape[0]
            for ie in range(self.n_ev):  
                thiseng_abs = self.absdata[:,:,ie]
                self.i0datahist[ie] = np.sum(thiseng_abs[i0_indices])*invnumel


        self.evi0 = self.ev.copy()
        self.i0data = self.i0datahist 

        self.calculate_optical_density()    
    
    
        return    
    
#----------------------------------------------------------------------   
# Normalize the data: calculate optical density matrix D 
    def calculate_optical_density(self):
        n_pixels = self.n_cols*self.n_rows
        self.od = np.empty((self.n_cols, self.n_rows, self.n_ev))
        
        #little hack to deal with rounding errors
        self.evi0[self.evi0.size-1] += 0.001
        
        fi0int = scipy.interpolate.interp1d(self.evi0,self.i0data, kind='cubic', bounds_error=False, fill_value=0.0)      
        i0 = fi0int(self.ev)
        
        #zero out all negative values in the image stack
        negative_indices = np.where(self.imagestack < 0)
        if negative_indices:
            self.imagestack[negative_indices] = 0.01
                

        for i in range(self.n_ev):
            self.od[:,:,i] = - np.log(self.absdata[:,:,i]/i0[i])
        
        #clean up the result
        nan_indices = np.where(np.isfinite(self.od) == False)
        if nan_indices:
            self.od[nan_indices] = 0
            
        self.od3d = self.od.copy()
            
        #Optical density matrix is rearranged into n_pixelsxn_ev
        self.od = np.reshape(self.od, (n_pixels, self.n_ev), order='F')
        
        self.original_od3d = self.od3d.copy()

           
        
        return
            
        
        
        
        
        
        
        
       
            
    