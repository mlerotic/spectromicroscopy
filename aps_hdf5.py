'''
Created on Jun 12, 2011

/ [root HDF5 group]
    spectromicroscopy/[HDF5 group]
            file_creation_date [string - ISO 8601 format]
            version [string]
            comment [string]
            experimenter/ [HDF5 group]
                    name [string]
                    role [string]
                    affiliation [string]
                    address [string]
                    phone [string]
                    email [string]
                    facility_user_id [string]
            sample/ [HDF5 group]
                    name [string]
                    description [string]
                    preparation_date [string - ISO 8601 format]
                    chemical_formula [string - abbreviated CIF format]
            data_stack [float] - n-dim dataset 
            dataset/ [HDF5 group]
                    type [string]  - spectrum, 3D stack, multielement ...
                    n_columns [integer]
                    n_rows [integer]
                    n_energies [integer]
                    n_detector_elements [integer]
                    x_dist [float]
                    y_dist [float]
                    ev [float]
                    msec [float]
            i0/ [HDF5 group]
                    i0data [float]
                    evi0[float]



@author: Mirna Lerotic
'''
import numpy as np
import scipy as scy
import h5py 


#----------------------------------------------------------------------
class h5:
    def __init__(self):
        pass
    
#----------------------------------------------------------------------
    def read_h5(self, filename):
        
        
        f = h5py.File(filename, 'r') 

        sm = f['spectromicroscopy']       
        filedate = sm.attrs['file_creation_date']        
        version = sm.attrs['version']        
        comment = sm.attrs['comment'] 
          
        dataset = sm['data_stack']
        self.absdata = dataset[...]
        
        ds = sm['dataset']
        self.n_cols = ds.attrs['n_columns'] 
        self.n_rows =  ds.attrs['n_rows'] 
        self.n_ev = ds.attrs['n_energies']  
        self.n_det = ds.attrs['n_detector_elements'] 
        self.ev = ds.attrs['ev'] 
        self.x_dist = ds.attrs['x_dist'] 
        self.y_dist = ds.attrs['y_dist']
        self.msec = ds.attrs['msec'] 
        type = ds.attrs['type'] 
        
        exp = sm['experimenter']
        exp_name = exp.attrs['name']   
        exp_role = exp.attrs['role']            
        exp_affil = exp.attrs['affiliation']     
        esp_addrs = exp.attrs['address']         
        exp_phone = exp.attrs['phone']           
        exp_email = exp.attrs['email']           
        exp_facuid = exp.attrs['facility_user_id']         
        
        sample = sm['sample']
        samp_name = sample.attrs['name']        
        samp_desc = sample.attrs['description']  
        samp_prepdate = sample.attrs['preparation_date']            
        samp_chemfor = sample.attrs['chemical_formula']            
                   
        i0gp = sm['i0/']
        self.i0data = i0gp.attrs['i0data']            
        self.evi0 = i0gp.attrs['evi0']
        
        verbose = 0
        
        if verbose == 1: 
            print  'filename ', filename
            print list(sm)
            print 'File creation date ', filedate
            print 'Version ', version
            print 'Comment ', comment
            print 'Data array shape: ', self.absdata.shape
            print 'n columns ', self.n_cols
            print 'n_rows ', self.n_rows
            print 'n_ev ', self.n_ev
            print 'n_detector ', self.n_det
            print 'ev array ', self.ev
            #print 'x dist ', self.x_dist
            #print 'y_dist ', self.y_dist
            print 'type ', type 
            print 'i0 data ', self.i0data
            print 'evi0 ', self.evi0
        
                
        f.close()
        
        self.imagestack = np.empty((self.n_cols*self.n_rows*self.n_ev))
        self.imagestack = np.reshape(self.absdata, (self.n_cols*self.n_rows*self.n_ev), order='F')
        
        return
    


