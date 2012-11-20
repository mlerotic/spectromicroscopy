# 
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
# 
#   Copyright (C) 2011 Mirna Lerotic, 2nd Look
#   http://2ndlook.co/products.html
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
import scipy.spatial
import scipy.ndimage
from scipy.cluster.vq import  kmeans2

import scipy as sp
mmult = sp.dot

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
        

        #try:
        if True:
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
                min_val = np.amin(self.pcaimages[:,:,i])
                max_val = np.amax(self.pcaimages[:,:,i])
                self.pcaimagebounds[i] = np.amax((np.abs(min_val), np.abs(max_val)))

            
            #calculate variance captured by the pca components
            self.variance = self.eigenvals.copy()
            
            totalvar = self.variance.sum()
            
            self.variance = self.variance/totalvar
            
            #Scree plot - find an elbow in the curve - between 1 and 20 components
            maxpoints = 25
            #Find a line between first (x1, y1) and last point (x2, y2) and calculate distances:
            y2 = np.log(self.eigenvals[maxpoints])
            x2 = maxpoints
            y1 = np.log(self.eigenvals[0])
            x1 = 0
            
            #Calculate distances between all the points and the line x1 and x2 are points on the line and x0 are eigenvals
            distance = np.zeros((25))
            for i in range(25):
                y0 = np.log(self.eigenvals[i])
                x0=i
                distance[i] = np.abs((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1))/np.math.sqrt((x2-x1)**2+(y2-y1)**2)  
            
            #Point with the largest distance is the "elbow"
            sigpca = np.argmax(distance)

            self.numsigpca = sigpca + 1
            

                        
#        except:
#            print "pca not converging"
            

        return    
    
#----------------------------------------------------------------------   
# Find clusters 
    def calculate_clusters(self, nclusters, remove1stpca = 0):
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
       

       
        indx = np.zeros(npixels)

        clustercentroids, indx = kmeans2(od_reduced, nclusters, iter=200, minit = 'points' )
        
              
       
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
     
        #Calculate SSE Sum of Squared errors
        indx = np.reshape(self.cluster_indices, (npixels), order='F')
        self.sse = np.zeros((npixels))
        for i in range(npixels):
            clind = indx[i]
            self.sse[i] = np.sum(np.square(self.stack.od[i,:]-self.clusterspectra[clind,:]))         
           
        self.sse = np.reshape(self.sse, (self.stack.n_cols, self.stack.n_rows), order='F')
        


        #Check the validity of cluster analysis and if needed add another cluster
        new_cluster_indices = self.cluster_indices.copy()
        new_nclusters = nclusters
        recalc_clusters = False
        
        for i in range(nclusters):
            clind = np.where(self.cluster_indices == i)
            cl_sse_mean = np.mean(self.sse[clind])
            cl_see_std = np.std(self.sse[clind])
            #print i, cl_sse_mean, cl_see_std 
            
            sigma9 = cl_sse_mean+9*cl_see_std
            maxsse = np.max(self.sse[clind])
            if (maxsse > sigma9): 
                #print 'have new cluster', max, sigma6
                recalc_clusters = True
                sse_helper = np.zeros((self.stack.n_cols, self.stack.n_rows), dtype=np.int)
                sse_helper[clind] = self.sse[clind]
                newcluster_ind = np.where(sse_helper > sigma9)
                new_cluster_indices[newcluster_ind] = new_nclusters
                new_nclusters += 1
                        
        
        if recalc_clusters == True:
            nclusters = new_nclusters
            self.cluster_indices = new_cluster_indices
            self.clusterspectra = np.zeros((nclusters, self.stack.n_ev))
            self.clustersizes = np.zeros((nclusters,), dtype=np.int)
            for i in range(nclusters):
                clind = np.where(self.cluster_indices == i)
                self.clustersizes[i] = self.cluster_indices[clind].shape[0]

                for ie in range(self.stack.n_ev):  
                    thiseng_od = self.stack.od3d[:,:,ie]
                    self.clusterspectra[i,ie] = np.sum(thiseng_od[clind])/self.clustersizes[i]
         
            #Calculate SSE Sum of Squared errors
            indx = np.reshape(self.cluster_indices, (npixels), order='F')
            self.sse = np.zeros((npixels))
            for i in range(npixels):
                clind = indx[i]
                self.sse[i] = np.sum(np.square(self.stack.od[i,:]-self.clusterspectra[clind,:]))         
               
            self.sse = np.reshape(self.sse, (self.stack.n_cols, self.stack.n_rows), order='F')
        
        
        self.cluster_distances = self.sse
        
        return int(nclusters)
        
    
        
#----------------------------------------------------------------------   
# Find clusters using K-means clustering
    def calculate_clusters_kmeans(self, N_clusters, data):

        N_Variables = data.shape[1] # same as n_significant_components
        N_Samples = data.shape[0] # same as n_pixels
        N_Iterations = 10
        
 
        #Initial 'learning rate'.
        LearningRates = 0.3-0.2*np.arange(N_Iterations)/(N_Iterations-1)
        
        # Normal random cluster weights.
        cluster_weights = (-1.+ 2*np.random.random((N_Variables, N_clusters)))

        
        cluster_distances = np.zeros((N_Samples))
        cluster_histogram = np.zeros((N_clusters), dtype=np.int32)
        Tempcluster_indices = np.zeros((N_Samples), dtype=np.int32)        
        cluster_indices = np.zeros((N_Samples), dtype=np.int32)                          
                        
                        
        WorkCol = np.ones((1,N_clusters))


        Metric = np.zeros((N_clusters))
        # Start by picking a percentage of the pixels at random,
        # and using them to start finding our cluster centers
        n_random_pixels = int(0.20*float(N_Samples))
        random_sample_indices = np.arange((N_Samples))
        np.random.shuffle(random_sample_indices)

        this_learning_rate = LearningRates[0]

        for i_sample in range(n_random_pixels):
            Sample = random_sample_indices[i_sample]
            Vector =  np.outer(data[Sample, :], WorkCol) - cluster_weights
            
            #Calculate distances
            for i_cluster in range(N_clusters):
                Metric[i_cluster] = np.sqrt(np.dot(Vector[:, i_cluster],Vector[:, i_cluster]).conj())
                
            MinIndex = np.argmin(Metric)
            cluster_weights[:,MinIndex] = this_learning_rate* Vector[:,MinIndex] + cluster_weights[:,MinIndex]

            #cluster_indices[random_sample_indices[i_sample]] = MinIndex


    
        # Random ordering of sample indices.
        random_ordered = np.arange((N_Samples))
        np.random.shuffle(random_ordered)


        New_RMSDistanceIterations = np.zeros((N_Iterations))
        for i_iteration in range(N_Iterations):
      
            this_max_distance = float(0)
            this_learning_rate = LearningRates[i_iteration]
            for Sample in range (N_Samples):
                RSample = random_ordered[Sample]
                # In our case the data array is
                # d_reduced(n_significant_components,n_pixels), and we have
                # WorkCol(1,n_clusters).  Calculate
                # Vector(n_significant_components,n_clusters) by multiplying
                # all the components for this pixel by n_clusters values of
                # 1 to pick them off, and then subtracting from the result
                # the current guess of the weights (the cluster centers).
                temp = np.outer(data[RSample, :], WorkCol)
                Vector =  temp - cluster_weights

                #Calculate distances
                for i_cluster in range(N_clusters):
                    Metric[i_cluster] = np.sqrt(np.dot(Vector[:, i_cluster],Vector[:, i_cluster]))
                
                MinIndex = np.argmin(Metric)
                MinMetric = Metric[MinIndex]
                cluster_weights[:,MinIndex] = this_learning_rate* Vector[:,MinIndex] + cluster_weights[:,MinIndex]

                cluster_distances[RSample] = MinMetric
                
                if (i_iteration == (N_Iterations-1)) :
                    Tempcluster_indices[RSample] = MinIndex
                    cluster_histogram[MinIndex] = cluster_histogram[MinIndex]+1

            # Since we're talking about distances from the cluster
            # center, which is in some ways an average of pixel positions,
            # we use (N_Samples-1) in the denominator.
            #New_RMSDistanceIterations[i_iteration] = np.sqrt(np.sum(cluster_distances**2)/float(N_Samples-1))

        # Next we sort the data with the cluster with the most members
        # first
        count_indices = np.flipud(np.argsort(cluster_histogram))
        cluster_histogram = cluster_histogram[count_indices]
        cluster_indices = np.zeros((N_Samples), dtype=np.int32)     
        ClustersFound = 0L
        for i_cluster in range(N_clusters):
            i_temp_cluster = count_indices[i_cluster]
            these_pixels = np.where(Tempcluster_indices == i_temp_cluster)
            if Tempcluster_indices[these_pixels].any():
                cluster_indices[these_pixels] = i_cluster
                ClustersFound = ClustersFound + 1
    
#        cluster_histogram = cluster_histogram[0:(ClustersFound-1)]
#        cluster_weights = cluster_weights[:, 0:(ClustersFound-1)]   
#
#        # Recalculate the cluster centers to be equal to the average of the
#        # pixel weights. 
#        for i_cen in range(N_clusters):
#            cluster_members = np.where(cluster_indices == i_cen)
#            n_mem = cluster_indices[cluster_members].size
##            if n_mem>0:
##                WorkRow2=np.ones((n_mem))
##                cluster_weights[:,i_cen]=np.outer(data[cluster_members,:],WorkRow2)/n_mem

      
        return cluster_weights.T, cluster_indices
       
        
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
            spectrum_evdata, spectrum_data, spectrum_common_name = self.stack.read_csv(filename)
                                    
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
        



#-----------------------------------------------------------------------------
# Find key energies by finding peaks and valleys in significant pca spectra
    def calc_key_engs(self, threshold):
 
        key_engs = []
               
        for i in range(self.numsigpca):
            pcaspectrum = self.eigenvecs[:,i]
            pmax,pmin = self.find_peaks(pcaspectrum, threshold, x = self.stack.ev)
            
            for i in range(len(pmin)):
                key_engs.append(pmin[i][0])
            for i in range(len(pmax)):
                key_engs.append(pmax[i][0])
                
        key_engs = np.array(key_engs)
        
        #Sort the energies and remove double entries
        key_engs = np.unique(key_engs)
                
        return key_engs


#-----------------------------------------------------------------------------
#Peakfinder
    def find_peaks(self, v, delta, x = None):
        """
        Converted from MATLAB script at http://billauer.co.il/peakdet.html by endolith
        https://gist.github.com/250860
        Returns two arrays
        function [maxtab, mintab]=peakdet(v, delta, x)
        %PEAKDET Detect peaks in a vector
        % [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
        % maxima and minima ("peaks") in the vector V.
        % MAXTAB and MINTAB consists of two columns. Column 1
        % contains indices in V, and column 2 the found values.
        %
        % With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
        % in MAXTAB and MINTAB are replaced with the corresponding
        % X-values.
        %
        % A point is considered a maximum peak if it has the maximal
        % value, and was preceded (to the left) by a value lower by
        % DELTA.
        % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
        % This function is released to the public domain; Any use is allowed.
        """
        maxtab = []
        mintab = []
           
        if x is None:
            x = np.arange(len(v))
        
        v = np.asarray(v)
        
        if len(v) != len(x):
            print 'Input vectors v and x must have same length'
            return -1
        
        if not np.isscalar(delta):
            print 'Input argument delta must be a scalar'
            return -1
        
        if delta <= 0:
            print 'Input argument delta must be positive'
        
        mn, mx = np.Inf, -np.Inf
        mnpos, mxpos = np.NaN, np.NaN
        
        lookformax = True
        
        for i in np.arange(len(v)):
            this = v[i]
            if this > mx:
                mx = this
                mxpos = x[i]
            if this < mn:
                mn = this
                mnpos = x[i]
            
            if lookformax:
                if this < mx-delta:
                    maxtab.append((mxpos, mx))
                    mn = this
                    mnpos = x[i]
                    lookformax = False
            else:
                if this > mn+delta:
                    mintab.append((mnpos, mn))
                    mx = this
                    mxpos = x[i]
                    lookformax = True
    
        return np.array(maxtab), np.array(mintab)




#----------------------------------------------------------------------   
# Calculate Fast Independent Component Analysis; FASTICA uses Hyvarinen's 
# fixed-point algorithm 
# A. Hyvarinen. Fast and Robust Fixed-Point Algorithms for Independent 
#   Component Analysis. IEEE Transactions on Neural Networks 10(3):626-634, 1999.

    def calculate_fastica(self, mixedsig, numOfIC):

        mixedsig = mixedsig.transpose()
        
        print 'msig', mixedsig.shape

        # Remove the mean and check the data
        mixedmean = np.mean(mixedsig, axis=1)
        mixedsig = mixedsig - mixedmean[:,np.newaxis]

        Dim = mixedsig.shape[0]
        NumOfSampl = mixedsig.shape[1]
        
        print 'Dim, NumOfSampl',Dim,NumOfSampl


        # Default values for optional parameters
        verbose  = True

        # Default values for 'pcamat' parameters
        firstEig          = 1
        lastEig           = Dim
        interactivePCA    = 'off'

        # Default values for 'fpica' parameters
        approach          = 'defl'
        g                 = 'pow3'
        finetune          = 'off'
        a1                = 1
        a2                = 1
        myy               = 1
        stabilization     = 'off'
        epsilon           = 0.0001
        maxNumIterations  = 1000
        maxFinetune       = 5
        initState         = 'rand'
        guess             = 0
        sampleSize        = 1
        displayMode       = 'off'
        displayInterval   = 1


        # Parameters for fastICA 
        b_verbose = True
        

        # print information about data
        if b_verbose:
            print 'Number of signals:', Dim
            print 'Number of samples: ', NumOfSampl


        # Check if the data has been entered the wrong way,
        # but warn only... it may be on purpose

        if (Dim > NumOfSampl):
            if b_verbose:
                print 'Warning: '
                print 'The signal matrix may be oriented in the wrong way.'
                print 'In that case transpose the matrix.'

    
        # Calculating PCA
        # We already have the PCA data

        if b_verbose:
            print 'Values for PCA calculations supplied.\n'
            print 'PCA calculations not needed.\n' 
    
        # PCA was already calculated:
        D = np.identity(self.numsigpca)*self.eigenvals[0:self.numsigpca] 
        E = self.eigenvecs[:,0:self.numsigpca]


        # Calculate the whitening
           
        # Calculate the whitening and dewhitening matrices (these handle
        # dimensionality simultaneously).
        whiteningMatrix = mmult(np.linalg.inv (np.sqrt(D)), E.transpose())
        dewhiteningMatrix = mmult(E, np.sqrt(D))
        
        print 'wd=', whiteningMatrix.shape, dewhiteningMatrix.shape
        
        # Project to the eigenvectors of the covariance matrix.
        # Whiten the samples and reduce dimension simultaneously.
        if b_verbose:
            print 'Whitening...'
            whitesig =  np.dot(whiteningMatrix,mixedsig)
            print 'whitesig', whitesig.shape
            
        # Just some security...
        if np.sum(np.imag(whitesig)) != 0:
            print 'Whitened vectors have imaginary values.'
            

        # Calculating the ICA
  
        # Check some parameters
        # The dimension of the data may have been reduced during PCA calculations.
        # The original dimension is calculated from the data by default, and the
        # number of IC is by default set to equal that dimension.
  
        Dim = whitesig.shape[0]
  
        # The number of IC's must be less or equal to the dimension of data
        if numOfIC > Dim:
            numOfIC = Dim
            # Show warning only if verbose = 'on' and user supplied a value for 'numOfIC'
            if b_verbose:
                print 'Warning: estimating only ,', numOfIC, '  independent components' 
                print '(Cannot estimate more independent components than dimension of data)'

  
        # Calculate the ICA with fixed point algorithm.
        A, W = self.calc_fpica (whitesig,  whiteningMatrix, dewhiteningMatrix, approach, 
          numOfIC, g, finetune, a1, a2, myy, stabilization, epsilon, 
          maxNumIterations, maxFinetune, initState, guess, sampleSize, 
          displayMode, displayInterval, verbose)
        
        print 'A,W', A.shape, W.shape
        # Check for valid return
        if W.any():
            # Add the mean back in.
            if b_verbose:
                print 'Adding the mean back to the data.'
    
            icasig = mmult(W, mixedsig) + mmult(mmult(W, mixedmean), np.ones((1, NumOfSampl)))

   
        else:
            icasig = []
            
        print icasig.shape
 
        return icasig

#----------------------------------------------------------------------   
    def getSamples(self, max, percentage):
        Samples = (np.random.random((max,)) < percentage).nonzero()
        return Samples
        
#----------------------------------------------------------------------   
# Fixed point ICA
# This function is adapted from Hyvarinen's fixed point algorithm Matlab version
# [A, W] = fpica(whitesig, whiteningMatrix, dewhiteningMatrix, approach,
#        numOfIC, g, finetune, a1, a2, mu, stabilization, epsilon, 
#        maxNumIterations, maxFinetune, initState, guess, sampleSize,
#        displayMode, displayInterval, verbose);
# 
# Perform independent component analysis using Hyvarinen's fixed point
# algorithm. Outputs an estimate of the mixing matrix A and its inverse W.
#
# whitesig                              :the whitened data as row vectors
# whiteningMatrix                       :whitening matrix
# dewhiteningMatrix                     :dewhitening matrix
# approach      [ 'symm' | 'defl' ]     :the approach used (deflation or symmetric)
# numOfIC       [ 0 - Dim of whitesig ] :number of independent components estimated
# g             [ 'pow3' | 'tanh' |     :the nonlinearity used
#                 'gaus' | 'skew' ]     
# finetune      [same as g + 'off']     :the nonlinearity used in finetuning.
# a1                                    :parameter for tuning 'tanh'
# a2                                    :parameter for tuning 'gaus'
# mu                                    :step size in stabilized algorithm
# stabilization [ 'on' | 'off' ]        :if mu < 1 then automatically on
# epsilon                               :stopping criterion
# maxNumIterations                      :maximum number of iterations 
# maxFinetune                           :maximum number of iterations for finetuning
# initState     [ 'rand' | 'guess' ]    :initial guess or random initial state. See below
# guess                                 :initial guess for A. Ignored if initState = 'rand'
# sampleSize    [ 0 - 1 ]               :percentage of the samples used in one iteration
# displayMode   [ 'signals' | 'basis' | :plot running estimate
#                 'filters' | 'off' ]
# displayInterval                       :number of iterations we take between plots
# verbose       [ 'on' | 'off' ]        :report progress in text format

    def calc_fpica(self, X, whiteningMatrix, dewhiteningMatrix, approach, 
            numOfIC, g, finetune, a1, a2, myy, stabilization, 
            epsilon, maxNumIterations, maxFinetune, initState, 
            guess, sampleSize, displayMode, displayInterval, 
            b_verbose):
        

        vectorSize = X.shape[0]
        numSamples = X.shape[1]
        
        # Checking the value for approach
        if approach == 'symm':
            approachMode = 1
        elif approach == 'defl':
            approachMode = 2    
        else:
            print 'Illegal value for parameter approach:', approach
            return
        if b_verbose:
            print 'Used approach:', approach     
            
        #Checking the value for numOfIC
        if vectorSize < numOfIC:
            print 'Must have numOfIC <= Dimension!'     
            return
        
        # Checking the sampleSize
        if sampleSize > 1:
            sampleSize = 1
            if b_verbose:
                print 'Warning: Setting sampleSize to 1.\n'
   
        elif sampleSize < 1:
            if (sampleSize * numSamples) < 1000:
                sampleSize = np.min(1000/numSamples, 1)
                if b_verbose:
                    print 'Warning: Setting ampleSize to ',sampleSize,' samples=', np.floor(sampleSize * numSamples)

        print 'sample size', sampleSize
        if  b_verbose and (sampleSize < 1):
            print 'Using about  ',sampleSize*100,' of the samples in random order in every step.'
       
        
        # Checking the value for nonlinearity.

        if g == 'pow3':
            gOrig = 10
        elif g =='tanh':
            gOrig = 20
        elif g == 'gauss':
            gOrig = 30
        elif g == 'skew':
            gOrig = 40
        else:
            print 'Illegal value for parameter g: ', g

        if sampleSize != 1:
            gOrig = gOrig + 2

        if myy != 1:
            gOrig = gOrig + 1


        if b_verbose:
            print 'Used nonlinearity: ', g


        finetuningEnabled = 1
        if finetune == 'pow3':
            gFine = 10 + 1
        elif finetune == 'tanh':
            gFine = 20 + 1
        elif finetune == 'gauss':
            gFine = 30 + 1
        elif finetune == 'skew':
            gFine = 40 + 1
        elif finetune == 'off':
            if myy != 1:
                gFine = gOrig
            else :
                gFine = gOrig + 1
  
            finetuningEnabled = 0
        else:
            print 'Illegal value for parameter finetune :', finetune
            return

        if b_verbose and finetuningEnabled:
            print 'Finetuning enabled, nonlinearity: ', finetune


        if stabilization == 'on':
            stabilizationEnabled = 1
        elif stabilization == 'off':
            if myy != 1:
                stabilizationEnabled = 1
            else:
                stabilizationEnabled = 0

        else:
            print 'Illegal value for parameter stabilization: ', stabilization 


        if b_verbose and stabilizationEnabled:
            print 'Using stabilized algorithm.'
       
        # Some other parameters
        myyOrig = myy
        # When we start fine-tuning we'll set myy = myyK * myy
        myyK = 0.01
        # How many times do we try for convergence until we give up.
        failureLimit = 5


        usedNlinearity = gOrig
        stroke = 0
        notFine = 1
        long = 0
        
        
        # Checking the value for initial state.

        if initState == 'rand':
            initialStateMode = 0;
        elif initState == 'guess':
            if guess.shape[0] != whiteningMatrix.shape[1]:
                initialStateMode = 0
                if b_verbose:
                    print 'Warning: size of initial guess is incorrect. Using random initial guess.'
    
            else:
                initialStateMode = 1
                if guess.shape[0] < numOfIC:
                    if b_verbose:
                        print 'Warning: initial guess only for first ',guess.shape[0],' components. Using random initial guess for others.' 
    
                    guess[:, guess.shape[0] + 1:numOfIC] = np.random.uniform(-0.5,0.5,(vectorSize,numOfIC-guess.shape[0]))
                elif guess.shape[0]>numOfIC:
                    guess=guess[:,1:numOfIC]
                    print 'Warning: Initial guess too large. The excess column are dropped.'
  
                if b_verbose:
                    print 'Using initial guess.'

        else:
            print 'Illegal value for parameter initState:', initState
            return        
        
        
        # Checking the value for display mode.

        if (displayMode =='off') or (displayMode == 'none'):
            usedDisplay = 0
        elif (displayMode =='on') or (displayMode == 'signals'):
            usedDisplay = 1
            if (b_verbose and (numSamples > 10000)):
                print 'Warning: Data vectors are very long. Plotting may take long time.'
 
            if (b_verbose and (numOfIC > 25)):
                print 'Warning: There are too many signals to plot. Plot may not look good.'
 
        elif (displayMode =='basis'):
            usedDisplay = 2
            if (b_verbose and (numOfIC > 25)):
                print 'Warning: There are too many signals to plot. Plot may not look good.'
 
        elif (displayMode =='filters'):
            usedDisplay = 3
            if (b_verbose and (vectorSize > 25)):
                print 'Warning: There are too many signals to plot. Plot may not look good.'

        else:
            print 'Illegal value for parameter displayMode:', displayMode
            return


        # The displayInterval can't be less than 1...
        if displayInterval < 1:
            displayInterval = 1

        # Start ICA calculation
        if b_verbose:
            print 'Starting ICA calculation...'
            
            
        # SYMMETRIC APPROACH
        if approachMode == 1:
            print 'Symmetric approach under construction'
            return
        
        
        # DEFLATION APPROACH
        elif approachMode == 2:      
            
            
            B = np.zeros((numOfIC, numOfIC))
  
            # The search for a basis vector is repeated numOfIC times.
            round = 0  
            numFailures = 0 
                      
            while (round < numOfIC):   
                myy = myyOrig
                usedNlinearity = gOrig
                stroke = 0
                notFine = 1
                long = 0
                endFinetuning = 0
                
                
                # Show the progress...
                if b_verbose:
                    print 'IC :', round
        
                # Take a random initial vector of length 1 and orthogonalize it
                # with respect to the other vectors.
                if initialStateMode == 0:
                    w = np.random.standard_normal((vectorSize,))
                    
                elif initialStateMode == 1:
                    w=mmult(whiteningMatrix,guess[:,round])
                    
       
                w = w - mmult(mmult(B, B.T), w) 
                norm =  sp.sqrt((w*w).sum())
                w = w / norm
                
        
                wOld = np.zeros(w.shape)
                wOld2 = np.zeros(w.shape)
                    
                    
                # This is the actual fixed-point iteration loop.
                #    for i = 1 : maxNumIterations + 1
                i = 1
                gabba = 1
                while (i <= (maxNumIterations + gabba)):
                    if (usedDisplay > 0):
                        print 'display'
          
                    #Project the vector into the space orthogonal to the space
                    # spanned by the earlier found basis vectors. Note that we can do
                    # the projection with matrix B, since the zero entries do not
                    # contribute to the projection.
                    w =  w - mmult(mmult(B, B.T), w)
                    norm =  sp.sqrt((w*w).sum())
                    w = w / norm
    
                    if notFine:
                        if i == (maxNumIterations + 1):   
                            if b_verbose:
                                print 'Component number',round,'  did not converge in ',maxNumIterations, 'iterations.'
          
                            round = round - 1
                            numFailures = numFailures + 1
                            if numFailures > failureLimit:
                                if b_verbose:
                                    print 'Too many failures to converge ', numFailures,' Giving up.'
            
                                if round == 0:
                                    A=[]
                                    W=[]
                                    return
                            break                    
                    else:
                        if i >= endFinetuning:
                            #So the algorithm will stop on the next test...
                            wOld = w.copy()
                            
                    # Show the progress...
                    if b_verbose:
                        print '.'
                        
                    # Test for termination condition. Note that the algorithm has
                    # converged if the direction of w and wOld is the same, this
                    # is why we test the two cases.
                    normm = sp.sqrt(((w - wOld)*(w - wOld)).sum())
                    normp = sp.sqrt(((w + wOld)*(w + wOld)).sum())
                    conv = min((normm), normp)
                    if  (conv < epsilon):
                        if finetuningEnabled and notFine:
                            if b_verbose:
                                print 'Initial convergence, fine-tuning: '
                            notFine = 0
                            gabba = maxFinetune
                            wOld = np.zeros(w.shape)
                            wOld2 = np.zeros(w.shape)
                            usedNlinearity = gFine
                            myy = myyK * myyOrig
          
                            endFinetuning = maxFinetune + i
                       
                        else:
                            numFailures = 0
                            # Save the vector
                            B[:, round] = w.copy()

          
                            # Calculate the de-whitened vector.
                            A = np.dot(dewhiteningMatrix, w)
                            # Calculate ICA filter.
                            W = np.dot(w.transpose(), whiteningMatrix)
                            
                            # Show the progress...
                            if b_verbose:
                                print 'computed ( ',i,' steps ) '
                                
                            break
                    elif stabilizationEnabled:
                                
                        if (not stroke) and (np.linalg.norm(w - wOld2) < epsilon or np.linalg.norm(w + wOld2) < epsilon):
                            stroke = myy
                            if b_verbose:
                                print 'Stroke!'
                            myy = .5*myy
                            if np.mod(usedNlinearity,2) == 0:
                                usedNlinearity = usedNlinearity + 1
        
                        elif stroke:
                            myy = stroke
                            stroke = 0
                            if (myy == 1) and (np.mod(usedNlinearity,2) != 0):
                                usedNlinearity = usedNlinearity - 1
         
                        elif (notFine) and (not long) and (i > maxNumIterations / 2):
                            if b_verbose:
                                print 'Taking long (reducing step size) '
                                long = 1
                                myy = .5*myy
                                if np.mod(usedNlinearity,2) == 0:
                                    usedNlinearity = usedNlinearity + 1
    
          
                    wOld2 = wOld
                    wOld = w      
                                
                    # pow3
                    if usedNlinearity == 10:
                        u = mmult(X.T, w)
                        w = mmult(X, u*u*u)/numSamples - 3.*w          

                    elif usedNlinearity == 11:
                        u = mmult(X.T, w)
                        EXGpow3 = mmult(X, u*u*u)/numSamples
                        Beta = mmult(w.T, EXGpow3)
                        w = w - myy * (EXGpow3 - Beta*w)/(3-Beta)         
                                       
                    elif usedNlinearity == 12:
                        Xsub = self._get_rsamples(X)
                        u = mmult(Xsub.T, w)
                        w = mmult(Xsub, u*u*u)/Xsub.shape[1] - 3.*w                        

                    elif usedNlinearity == 13:
                        Xsub=X[:,self.getSamples(numSamples, sampleSize)]                        
                        u = mmult(Xsub.T, w)
                        EXGpow3 = mmult(Xsub, u*u*u)/Xsub.shape[1]
                        Beta = mmult(w.T, EXGpow3)
                        w = w - myy * (EXGpow3 - Beta*w)/(3-Beta)
                        
                    # tanh
                    elif usedNlinearity == 20:
                        u = mmult(X.T, w)
                        tang = sp.tanh(a1 * u)
                        temp = mmult((1. - tang*tang).sum(axis=0), w)
                        w = (mmult(X, tang) - a1*temp)/numSamples
                        
                    elif usedNlinearity == 21:
                        u = mmult(X.T, w)
                        tang = sp.tanh(a1 * u)
                        Beta = mmult(u.T, tang)
                        temp = (1. - tang*tang).sum(axis=0)
                        w = w-myy*((mmult(X, tang)-Beta*w)/(a1*temp-Beta))
                                                
                    elif usedNlinearity == 22:
                        Xsub=X[:,self.getSamples(numSamples, sampleSize)]
                        u = mmult(Xsub.T, w)
                        tang = sp.tanh(a1 * u)
                        temp = mmult((1. - tang*tang).sum(axis=0), w)
                        w = (mmult(Xsub, tang) - a1*temp)/Xsub.shape[1]
                                                
                    elif usedNlinearity == 23:
                        Xsub=X[:,self.getSamples(numSamples, sampleSize)]
                        u = mmult(Xsub.T, w)
                        tang = sp.tanh(a1 * u)
                        Beta = mmult(u.T, tang)
                        w = w - myy * ((mmult(Xsub, tang)-Beta*w) /
                                      (a1*(1. - tang*tang).sum(axis=0) -
                                       Beta))                        

                    # gauss
                    elif usedNlinearity == 30:
                        # This has been split for performance reasons.
                        u = mmult(X.T, w)
                        u2 = u*u
                        ex = sp.exp(-a2*u2*0.5)
                        gauss =  u*ex
                        dgauss = (1. - a2 *u2)*ex
                        w = (mmult(X, gauss)-mmult(dgauss.sum(axis=0), w))/numSamples                        
                        

                    elif usedNlinearity == 31:
                        u = mmult(X.T, w)
                        u2 = u*u
                        ex = sp.exp(-a2*u2*0.5)
                        gauss =  u*ex
                        dgauss = (1. - a2 *u2)*ex
                        Beta = mmult(u.T, gauss)
                        w = w - myy*((mmult(X, gauss)-Beta*w)/
                                    (dgauss.sum(axis=0)-Beta))                        
                        
                    elif usedNlinearity == 32:
                        Xsub=X[:,self.getSamples(numSamples, sampleSize)]
                        u = mmult(Xsub.T, w)
                        u2 = u*u
                        ex = sp.exp(-a2*u2*0.5)
                        gauss =  u*ex
                        dgauss = (1. - a2 *u2)*ex
                        w = (mmult(Xsub, gauss)-
                             mmult(dgauss.sum(axis=0), w))/Xsub.shape[1]                        
                        
                    elif usedNlinearity == 33:
                        Xsub=X[:,self.getSamples(numSamples, sampleSize)]
                        u = mmult(Xsub.T, w)
                        u2 = u*u
                        ex = sp.exp(-a2*u2*0.5)
                        gauss =  u*ex
                        dgauss = (1. - a2 *u2)*ex
                        Beta = mmult(u.T, gauss)
                        w = w - myy*((mmult(Xsub, gauss)-Beta*w)/
                                    (dgauss.sum(axis=0)-Beta))                        
                        
                    # skew
                    elif usedNlinearity == 40:
                        u = mmult(X.T, w)
                        w = mmult(X, u*u)/numSamples                       
                        
                    elif usedNlinearity == 41:
                        u = mmult(X.T, w)
                        EXGskew = mmult(X, u*u) / numSamples
                        Beta = mmult(w.T, EXGskew)
                        w = w - myy * (EXGskew - mmult(Beta, w))/(-Beta)                        
                        
                    elif usedNlinearity == 42:                 
                        Xsub=X[:,self.getSamples(numSamples, sampleSize)]
                        u = mmult(Xsub.T, w)
                        w = mmult(Xsub, u*u)/Xsub.shape[1]

                    elif usedNlinearity == 43:
                        Xsub=X[:,self.getSamples(numSamples, sampleSize)]
                        u = mmult(Xsub.T, w)
                        EXGskew = mmult(Xsub, u*u) / Xsub.shape[1]
                        Beta = mmult(w.T, EXGskew)
                        w = w - myy * (EXGskew - Beta*w)/(-Beta)
                        
                    else:
                        print 'Code for desired nonlinearity not found!'
                        return
                                 
                    # Normalize the new w.
                    norm = sp.sqrt((w*w).sum())
                    w = w / norm
                    i = i + 1  
                round = round + 1
            
            if b_verbose:
                print 'Done.'
            

        # In the end let's check the data for some security
        if A.imag.any():
            if b_verbose:
                print 'Warning: removing the imaginary part from the result.'
            A = A.real
            W = W.imag

        return A, W
                
                      
