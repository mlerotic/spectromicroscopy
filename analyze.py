'''
Created on Nov 22, 2010

@author: Mirna Lerotic
'''
from __future__ import division

import numpy as np
import scipy.interpolate
from scipy.cluster.vq import  kmeans2

import warnings
warnings.simplefilter('ignore', DeprecationWarning)



#----------------------------------------------------------------------
class analyze:
    def __init__(self):
        pass

#----------------------------------------------------------------------   
    def setdata(self, stkdata):
        self.stack = stkdata
        
            
#----------------------------------------------------------------------   
# Calculate pca 
    def calculate_pca(self):
        #covariance matrix
        n_pix = self.stack.n_cols*self.stack.n_rows

        od = self.stack.od
        
        #normalize od spectra - not used in pca_gui.pro
        norms = np.apply_along_axis(np.linalg.norm, 1, od)
        odn = np.zeros((n_pix, self.stack.n_ev))
        for i in range(n_pix):
            odn[i,:] = od[i,:]/np.linalg.norm(od[i,:])
            
       
        covmatrix = np.dot(od.T,od)
  
        self.pcaimages = np.zeros((self.stack.n_cols, self.stack.n_rows, self.stack.n_ev))
        self.pcaimagebounds = np.zeros((self.stack.n_ev))
        

        try:
            self.eigenvals, self.eigenvecs = np.linalg.eigh(covmatrix)

            #sort the eigenvals and eigenvecs       
            perm = np.argsort(-np.abs(self.eigenvals))
            self.eigenvals = self.eigenvals[perm]
            self.eigenvecs = self.eigenvecs[:,perm]
            

            self.pcaimages = np.dot(od,self.eigenvecs)            
            
            #calculate eigenimages
            self.pcaimages = np.reshape(self.pcaimages, (self.stack.n_cols, self.stack.n_rows, self.stack.n_ev), order='F')
                
            #Find bounds for displaying color-tables
            for i in range(self.stack.n_ev):
                min_val = np.min(self.pcaimages[:,:,i])
                max_val = np.max(self.pcaimages[:,:,i])
                self.pcaimagebounds[i] = np.max((np.abs(min_val), np.abs(max_val)))

            
            #calculate variance captured by the pca components
            self.variance = self.eigenvals.copy()
            
            totalvar = self.variance.sum()
            
            self.variance = self.variance/totalvar
            
            #Keiser criterion - select number of significant components - evals > 1
            sigpca = np.extract(self.eigenvals>1,self.eigenvals)
            self.numsigpca = sigpca.size
 
                        
        except:
            print "pca not converging"
            

        return    
    
#----------------------------------------------------------------------   
# Find clusters using kmeans
    def calculate_clusters_kmeans(self, nclusters, remove1stpca = 0):
        #Reduced data matrix od_reduced(n_pixels,n_significant_components)
        #od_reduced = np.zeros((self.stack.n_cols, self.stack.n_rows, self.numsigpca))
        

        npixels = self.stack.n_cols * self.stack.n_rows
        
        inverse_n_pixels = 1./float(npixels)
        
        if remove1stpca == 0 : 
            od_reduced = self.pcaimages[:,:,0:self.numsigpca]
            od_reduced = np.reshape(od_reduced, (npixels,self.numsigpca), order='F')
        else: 
            od_reduced = self.pcaimages[:,:,1:self.numsigpca]
            od_reduced = np.reshape(od_reduced, (npixels,self.numsigpca-1), order='F')
        
        

        #kmeans(obs,k_or_guess,iter=20,thresh=1e-5)
        
        indx = np.zeros(npixels)

        res, indx = kmeans2(od_reduced, nclusters, iter=100, minit = 'points' )
        
        indx = np.reshape(indx, (self.stack.n_cols, self.stack.n_rows), order='F')
        
        self.clustersizes = np.zeros((nclusters,), dtype=np.int)

        
        for i in range(nclusters):
            clind = np.where(indx == i)
            self.clustersizes[i] = indx[clind].shape[0]
                    
        #sort the data with the cluster with the most members first  
        count_indices = np.argsort(self.clustersizes)
        count_indices = count_indices[::-1]
                
        self.cluster_indices = np.zeros((self.stack.n_cols, self.stack.n_rows), dtype=np.int)
              
        self.clusterspectra = np.zeros((nclusters, self.stack.n_ev))
               
        for i in range(nclusters):
            clind = np.where(indx == count_indices[i])
            self.cluster_indices[clind] = i
            self.clustersizes[i] = self.cluster_indices[clind].shape[0]

            for ie in range(self.stack.n_ev):  
                thiseng_od = self.stack.od3d[:,:,ie]
                self.clusterspectra[i,ie] = np.sum(thiseng_od[clind])/self.clustersizes[i]
     

      
        
#----------------------------------------------------------------------   
# Find clusters using EM clustering
    def calculate_clusters_em(self, nclusters):
        #Reduced data matrix od_reduced(n_pixels,n_significant_components)
        #od_reduced = np.zeros((self.stack.n_cols, self.stack.n_rows, self.numsigpca))
        

        npixels = self.stack.n_cols * self.stack.n_rows
        
        inverse_n_pixels = 1./float(npixels)
        
        od_reduced = self.pcaimages[:,:,0:self.numsigpca]
        
        od_reduced = np.reshape(od_reduced, (npixels,self.numsigpca), order='F')
        

        #kmeans(obs,k_or_guess,iter=20,thresh=1e-5)
        
        self.indx = np.zeros(npixels)

        res, self.indx = kmeans2(od_reduced,5)
        
        self.indx = np.reshape(self.indx, (self.stack.n_cols, self.stack.n_rows), order='F')        

        
            

