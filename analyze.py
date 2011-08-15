'''
Created on Nov 22, 2010

@author: Mirna Lerotic
'''
from __future__ import division

import numpy as np
import scipy.interpolate
import scipy.spatial
import scipy.ndimage
from scipy.cluster.vq import  kmeans2

import warnings
warnings.simplefilter('ignore', DeprecationWarning)



#----------------------------------------------------------------------
class analyze:
    def __init__(self, stkdata):

        self.stack = stkdata
       
        self.target_spectra = 0
        self.tspectrum_loaded = 0
        self.n_target_spectra = 0
        self.tspec_names = []
        
#----------------------------------------------------------------------   
# Calculate pca 
    def delete_data(self):
        
        self.target_spectra = 0
        self.tspectrum_loaded = 0
        self.n_target_spectra = 0
        self.tspec_names = []
        
        self.pcaimages = 0
        self.pcaimagebounds = 0
        self.eigenvals = 0
        self.eigenvecs = 0
        
        self.cluster_distances = 0
        self.clustersizes = 0
        self.cluster_indices = 0
        self.clusterspectra = 0
                    
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
            
            if self.numsigpca > 5:
                self.numsigpca = 5
 
                        
        except:
            print "pca not converging"
            

        return    
    
#----------------------------------------------------------------------   
# Find clusters using kmeans
    def calculate_clusters_kmeans(self, nclusters, remove1stpca = 0):
        #Reduced data matrix od_reduced(n_pixels,n_significant_components)
        #od_reduced = np.zeros((self.stack.n_cols, self.stack.n_rows, self.numsigpca))
        
        self.nclusters = nclusters
                
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

        clustercentroids, indx = kmeans2(od_reduced, nclusters, iter=100, minit = 'points' )
        
        
        #calculate cluster distances
        self.cluster_distances = np.zeros((self.stack.n_cols*self.stack.n_rows))
        for i in range(npixels):
            clind = indx[i]
            self.cluster_distances[i] = scipy.spatial.distance.euclidean(od_reduced[i,:],clustercentroids[clind,:])
            
        self.cluster_distances = np.reshape(self.cluster_distances, (self.stack.n_cols, self.stack.n_rows), order='F')
                    
     
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
        
        
#-----------------------------------------------------------------------------
# Spectral analysis
# This routine reads in a mu spectrum in units of inverse microns.
# The spectrum is interpolated onto the energy range of the stack,
# and loaded into the matrix target_spectra(pca_gui_par.n_targets,n_ev).
# If there is a PCA calculation done, we find the fits to the
# target spectra from the components.
    def read_target_spectrum(self, filename = '', flat = False):
        # Load spectrum from a file
        spectrum_evdata = 0
        spectrum_data = 0
        spectrum_common_name = ' '
        if flat == False:
            spectrum_evdata, spectrum_data, spectrum_common_name = self.stack.read_xas(filename)
                                    
            # Map this spectrum onto our energy range - interpolate to ev
            ftspec = scipy.interpolate.interp1d(spectrum_evdata, spectrum_data, kind='cubic', bounds_error=False)      
            target_spectrum = np.reshape(ftspec(self.stack.ev), (1,self.stack.n_ev))
            
            #fix the edges if needed 
            if self.stack.ev[0]<spectrum_evdata[0]:
                indx = np.where(self.stack.ev<spectrum_evdata[0])
                target_spectrum[0,indx] = spectrum_data[0]
            if self.stack.ev[-1]>spectrum_evdata[-1]:
                indx = np.where(self.stack.ev>spectrum_evdata[-1])
                target_spectrum[0,indx] = spectrum_data[-1]
                
        else:
            target_spectrum = np.ones((1,self.stack.n_ev))
            spectrum_common_name = 'Flat'

     
        if self.tspectrum_loaded == 0:
            self.target_spectra = target_spectrum
            self.tspectrum_loaded = 1
            self.n_target_spectra += 1
        else:
            self.target_spectra = np.vstack((self.target_spectra,target_spectrum))
            self.n_target_spectra += 1
        self.tspec_names.append(spectrum_common_name)
        
        self.fit_target_spectra()
        self.calc_svd_maps()
        
        
            
#-----------------------------------------------------------------------------          
    def add_cluster_target_spectra(self):
        # Load spectrum from a file or cluster spectra 

        for i in range(self.nclusters):
   
            target_spectrum = self.clusterspectra[i,:]
       
            if self.tspectrum_loaded == 0:
                self.target_spectra = target_spectrum
                self.tspectrum_loaded = 1
                self.n_target_spectra += 1
            else:
                self.target_spectra = np.vstack((self.target_spectra,target_spectrum))
                self.n_target_spectra += 1
            self.tspec_names.append('Cluster '+str(i+1))
        
        self.fit_target_spectra()
        self.calc_svd_maps()  
        
#-----------------------------------------------------------------------------          
    def remove_spectrum(self, i_spec):   
        
        if self.n_target_spectra > 1:
            self.target_spectra = np.delete(self.target_spectra, i_spec, axis=0)
            
            del self.tspec_names[i_spec]
            self.n_target_spectra -= 1
        
            self.fit_target_spectra()
            self.calc_svd_maps()  
        else:
            self.target_spectra = 0
            self.tspectrum_loaded = 0
            self.n_target_spectra = 0
            
#-----------------------------------------------------------------------------          
    def move_spectrum(self, old_position, new_position):   
        
        temp_target_spectra = self.target_spectra.copy()
        temp_target_spectra[old_position,:] = self.target_spectra[new_position,:]
        temp_target_spectra[new_position,:] = self.target_spectra[old_position,:]
        self.target_spectra = temp_target_spectra
        
        temp_tspec_name = self.tspec_names[new_position]
        self.tspec_names[new_position] = self.tspec_names[old_position]
        self.tspec_names[old_position] = temp_tspec_name
                
        self.fit_target_spectra()
        self.calc_svd_maps()  

        
#-----------------------------------------------------------------------------
# This routine calculates:
#   - the transformation matrix T as target_spectrumfit_coeffs, and
#     the fits to the target spectra as target_fittedspectra
#   - the inverse of T
#   - the eigenvector matrix C(S_abstract,N) by transposing the matrix
#       CT(N,S_abstract)=evecs(N,S_abstract)
#   - the target maps t(P,S_targets) as targetfit_maps
    def fit_target_spectra(self):
#       We want to find T(S_physical,S_abstract) which is
#       CT(N,S_abstract)##mu(S_physical,N).  But mu(S_physical,N) is
#       known here as target_spectra(S_physical,N), and
#       CT(N,S_abstract) is just a limited version of evecs(N,S).
#       We will call T(S_physical,S_abstract) by the name
#       target_spectrumfit_coeffs(S_physical,S_abstract).
        CT = self.eigenvecs[:,0:self.numsigpca]
        self.target_pcafit_coeffs = np.dot(self.target_spectra, CT )
      
#       Now we get the target spectra as approximated from our
#       components with
#       mu(S_physical,N)=C(S_abstract,N)##T(S_physical,S_abstract).
        self.target_pcafit_spectra = np.dot(self.target_pcafit_coeffs, CT.T)

#       To get the maps, we need to find R(P,Sbar_abstract)
#       from ct(N,Sbar_abstract)##d(P,N), and we also
#       need to invert the transformation matrix.  Start by
#       finding the singular value decomposition of the
#       transformation matrix.
        U, s, V = np.linalg.svd(self.target_pcafit_coeffs, full_matrices=False)
                        
#       This gives T^{-1}(Sbar_abstract,S_physical)
        t_inverse = np.dot(np.dot(V.T, np.linalg.inv(np.diag(s))), U.T)
     
#       This is R(P,Sbar_abstract)=CT(N,Sbar_abstract)##D(P,N)
        r_matrix = np.dot(self.stack.od, CT)

#       and this gives us the maps as
#       t(P,S_physical) = T^{-1}(Sbar_abstract,S_physical)##R(P,Sbar_abstract)
#       but in fact it is t(P,n_targets)!
        self.target_pcafit_maps = np.dot(r_matrix, t_inverse)
        
        self.target_pcafit_maps = np.reshape(self.target_pcafit_maps, 
                                             (self.stack.n_cols, self.stack.n_rows, self.n_target_spectra), order='F')  
        
        #Find fit errors
        self.target_rms = (self.target_spectra-self.target_pcafit_spectra)**2
        self.target_rms = np.sqrt(np.sum(self.target_rms, axis=1)/self.stack.n_ev)
        

    
#-----------------------------------------------------------------------------
# This routine calculates the SVD composition maps
#   1. The optical density is calculated from the stack, and the
#      data matrix D of dimensions (pixels,energies) is formed
#   2. mu_inverse is calculated using SVD
#   3. svd_maps is calculated from multiplying mu_inverse ## d
    def calc_svd_maps(self, usefittedspectra = False):
        
        if usefittedspectra:
            U, s, V = np.linalg.svd(self.target_spectra, full_matrices=False)
        else:
            U, s, V = np.linalg.svd(self.target_pcafit_spectra, full_matrices=False)
      
        mu_inverse = t_inverse = np.dot(np.dot(V.T, np.linalg.inv(np.diag(s))), U.T)
        self.target_svd_maps = np.dot(self.stack.od, mu_inverse)
        
        self.target_svd_maps = np.reshape(self.target_svd_maps, 
                                             (self.stack.n_cols, self.stack.n_rows, self.n_target_spectra), order='F') 
        
