'''
Created on Jun 20, 2011

@author: Mirna
'''

import numpy as np
import scipy.interpolate
import h5py 
import datetime 

import x1a_stk
import aps_hdf5

#----------------------------------------------------------------------
class data(x1a_stk.x1astk,aps_hdf5.h5):
    def __init__(self):
        x1a_stk.x1astk.__init__(self)
        aps_hdf5.h5.__init__(self)

        pass
#----------------------------------------------------------------------   
    def new_data(self):
        self.n_cols = 0
        self.n_rows = 0
        self.n_ev = 0
               
        self.x_dist = 0       
        self.y_dist = 0

        self.ev = 0 
        self.msec = 0
        self.imagestack = 0             
        self.absdata = 0

        self.original_n_cols = 0
        self.original_n_rows = 0
        self.original_n_ev = 0
        self.original_ev = 0
        self.original_absdata = 0
        
        self.i0data =0
        self.evi0 = 0
        
        self.od = 0
        self.od3d = 0
        self.original_absdata = 0
        self.original_ev = 0
        self.original_n_cols = 0
        self.original_n_ev = 0
        self.original_n_rows = 0
        self.original_od3d = 0
           
#----------------------------------------------------------------------   
    def read_stk_i0(self, filename):
        x1a_stk.x1astk.read_stk_i0(self,filename)
        self.calculate_optical_density()

#---------------------------------------------------------------------- 
    def read_stk(self, filename):    
        self.new_data()  
        x1a_stk.x1astk.read_stk(self, filename)
  
#---------------------------------------------------------------------- 
    def read_h5(self, filename):    
        self.new_data()  
        aps_hdf5.h5.read_h5(self, filename)
        
        self.calculate_optical_density()
        
#----------------------------------------------------------------------
    def convert_stk_to_h5(self, filename):
        f = h5py.File(filename, 'w')  
        sm = f.create_group('spectromicroscopy')
        
        now = datetime.datetime.now()
        #print "Current date and time using strftime:"
        #print now.strftime("%Y-%m-%d %H:%M")
               
        sm.attrs['file_creation_date'] = now.strftime("%Y-%m-%d")
        sm.attrs['version'] = 0
        sm.attrs['comment'] = 'Converted from a x1a .stk file'
        
     
        dataset = sm.create_dataset('data_stack', data=self.absdata)
        ds = sm.create_group('dataset')
        ds.attrs['n_columns'] = self.n_cols
        ds.attrs['n_rows'] =  self.n_rows
        ds.attrs['n_energies'] = self.n_ev
        ds.attrs['n_detector_elements'] = 1
        ds.attrs['ev'] = self.ev
        ds.attrs['x_dist'] = self.x_dist
        ds.attrs['y_dist'] = self.y_dist
        ds.attrs['msec'] = self.msec
        ds.attrs['type'] = '3dstack'
        

        exp = sm.create_group('experimenter')
        exp.attrs['name'] = ' '  
        exp.attrs['role'] = ' '            
        exp.attrs['affiliation'] = ' '    
        exp.attrs['address'] = ' '         
        exp.attrs['phone'] = ' '           
        exp.attrs['email'] = ' '           
        exp.attrs['facility_user_id'] = ' '         
        
        sample = sm.create_group('sample')
        sample.attrs['name'] = ' '        
        sample.attrs['description'] = ' '   
        sample.attrs['preparation_date'] = ' '           
        sample.attrs['chemical_formula'] = ' '            

                   
        i0gp = sm.create_group('i0/')
        i0gp.attrs['i0data'] = self.i0data           
        i0gp.attrs['evi0'] = self.evi0
        
        f.close()
        
        
        return
    
    
#----------------------------------------------------------------------   
    def calc_histogram(self):
        #calculate average flux for each pixel
        self.averageflux = np.mean(self.absdata,axis=2)
        self.histogram = self.averageflux
        
        return
        
        
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
    def set_i0(self, i0data, evdata):

        self.evi0 = evdata
        self.i0data = i0data 

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
    
#----------------------------------------------------------------------   
    def write_xas(self, filename, evdata, i0data):
        print filename
    
        return  