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
        
        

#----------------------------------------------------------------------   
# Calculate Fast Independent Component Analysis; FASTICA uses Hyvarinen's 
# fixed-point algorithm 
# A. Hyvarinen. Fast and Robust Fixed-Point Algorithms for Independent 
#   Component Analysis. IEEE Transactions on Neural Networks 10(3):626-634, 1999.

    def calculate_fastica(self, mixedsig, numOfIC):

        mixedsig = mixedsig.transpose()

        # Remove the mean and check the data
        mixedmean = np.mean(mixedsig, axis=0)
        mixedsig = mixedsig - mixedmean*np.ones((1,mixedsig.shape[1]))

        Dim = mixedsig.shape[0]
        NumOfSampl = mixedsig.shape[1]


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
        whiteningMatrix = np.dot(np.linalg.inv (np.sqrt(D)), E.transpose())
        dewhiteningMatrix = np.dot(E, np.sqrt(D))
        
        # Project to the eigenvectors of the covariance matrix.
        # Whiten the samples and reduce dimension simultaneously.
        if b_verbose:
            print 'Whitening...'
            whitesig =  np.dot(whiteningMatrix,mixedsig)
            
        # Just some security...
        if np.sum(np.imag(whitesig)) != 0:
            print 'Whitened vectors have imaginary values.'
            

        # Calculating the ICA
  
        # Check some parameters
        # The dimension of the data may have been reduced during PCA calculations.
        # The original dimension is calculated from the data by default, and the
        # number of IC is by default set to equal that dimension.
  
        Dim = whitesig.shape[1]
  
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
            B = np.zeros((vectorSize))
  
            # The search for a basis vector is repeated numOfIC times.
            round = 1
  
            numFailures = 0 
                      
        while (round <= numOfIC):   
            myy = myyOrig
            usedNlinearity = gOrig
            stroke = 0
            notFine = 1
            long = 0
            endFinetuning = 0
            
            
            # Show the progress...
            if b_verbose:
                print 'IC :', round
    
            # Take a random initial vector of lenght 1 and orthogonalize it
            # with respect to the other vectors.
            if initialStateMode == 0:
                w = np.random.standard_normal((vectorSize,))
                
            elif initialStateMode == 1:
                w=np.dot(whiteningMatrix,guess[:,round])
   
            w = w - np.dot(B , np.dot( B.transpose(), w))     
            w = w / np.linalg.norm(w)
    
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
                w = w - np.dot(B , np.dot( B.transpose(), w))     
                w = w / np.linalg.norm(w)

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
                if (np.linalg.norm(w - wOld) < epsilon) or (np.linalg.norm(w + wOld) < epsilon):
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
                        B[:, round] = w
      
                        # Calculate the de-whitened vector.
                        A[:,round] = np.dot(dewhiteningMatrix, w)
                        # Calculate ICA filter.
                        W[round,:] = np.dot(w.transpose() * whiteningMatrix)
                        
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

      
                wOld2 = wOld.copy()
                wOld = w.copy()         
                            
                            
