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
        
        self.x_dist = np.fromfile(f, np.float32, self.n_cols)
        self.x_dist.byteswap(True)   
        
        self.y_dist = np.fromfile(f, np.float32, self.n_rows)
        self.y_dist.byteswap(True)     

        self.ev = np.fromfile(f, np.float32, self.n_ev)
        self.ev.byteswap(True)   
        
        self.msec = np.fromfile(f, np.float32, self.n_ev)
        self.msec.byteswap(True)   
        
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
        print 'tusam'
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
           
    
        return
    

            
        
        
        
        
        
        
        
       
            
    