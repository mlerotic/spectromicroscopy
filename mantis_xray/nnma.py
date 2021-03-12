import numpy as np
import scipy as sp
from scipy.integrate import trapz
import sys, time

from . import analyze

class nnma():

    def __init__(self, stkdata):

        self.stack = stkdata
    
        # Set default values for various NNMA parameters
        self.kNNMA = 5            # use no. of significant PCA components, otherwise default = 5              
        self.maxIters = 10         # default no. of iterations
        self.deltaErrorThreshold = 1e-3        # default threshold for change in reconstruction error
        self.initMatrices = 'Random'  # default initialization for matrices
    
        self.lambdaSparse = 0.        # sparseness regularization param
        self.lambdaClusterSim = 0.    # cluster spectra similarity regularization param
        self.lambdaSmooth = 0.       # smoothness regularization param


    # Fill initial matrices depending on user-specified initialization method
    def setParameters(self, kNNMA = 5, 
                            maxIters = 10,
                            deltaErrorThreshold = 1e-3,
                            initMatrices = 'Random',
                            lambdaSparse = 0.,
                            lambdaClusterSim = 0.,
                            lambdaSmooth = 0.):
        
        # Set the values for various NNMA parameters
        self.kNNMA = kNNMA                  
        self.maxIters = maxIters         
        self.deltaErrorThreshold = deltaErrorThreshold        
        self.initMatrices = initMatrices  
    
        self.lambdaSparse = lambdaSparse       
        self.lambdaClusterSim = lambdaClusterSim    
        self.lambdaSmooth = lambdaSmooth      


    def setClusterSpectra(self, clusterspectra):
        
        self.clusterspectra = clusterspectra.T.copy() 
        
        
    def setStandardsSpectra(self, standspectra):
        
        self.standspectra = standspectra

#---------------------------------------------------------------------------------------
# Define some helper functions for NNMA
#---------------------------------------------------------------------------------------

    # Fill initial matrices depending on user-specified initialization method
    def fillInitMatrices(self):
        print ('Init method:', self.initMatrices)
        if self.initMatrices == 'Random':
            muInit = np.random.rand(self.nEnergies, self.kNNMA)
            tInit = np.random.rand(self.kNNMA, self.nPixels)
        elif self.initMatrices == "Cluster":
            muInit = self.clusterspectra
            tInit = np.random.rand(self.kNNMA, self.nPixels)  # use SVD on mu to find t instead?
        elif self.initMatrices == "Standards":
            muInit = self.standspectra
            tInit = np.random.rand(self.kNNMA, self.nPixels)  # use SVD on mu to find t instead?
        return muInit, tInit

    # Update t
    def tUpdate(self, mu, t):
        tUpdateFactor = np.dot(mu.T, self.OD) / ( np.dot(mu.T, np.dot(mu, t)) + self.lambdaSparse + 1e-9 )
        tUpdated = t * tUpdateFactor
        return tUpdated

    # Update mu
    def muUpdate(self, mu, tUpdated):
        dCostSmooth_dMu = self.calcDCostSmooth_dMu(mu)
        muDiff = mu - self.muCluster
        muUpdateFactor = np.dot(self.OD, tUpdated.T) / ( np.dot(mu, np.dot(tUpdated, tUpdated.T))
                         + self.lambdaSmooth*dCostSmooth_dMu + self.lambdaClusterSim*2*muDiff + 1e-9 )
        muUpdated = mu * muUpdateFactor
        # Normalize each column of mu
        return muUpdated
  
    # Calculate sparseness contribution (JSparse) to cost function
    def calcCostSparse(self, t):
        costSparse = np.sum(np.sum(np.abs(t)))
        return costSparse

    # Calculate cluster spectra similarity contribution (JClusterSim) to cost function
    def calcCostClusterSim(self, mu):
        muDiff = mu - self.muCluster
        costClusterSim = (np.linalg.norm(muDiff))**2 
        return costClusterSim

    # Calculate smoothness contribution (JSmooth) to cost function
    def calcCostSmooth(self, mu):
        return 0.
  
    # Calculate dJSmooth/dMu needed in mu update algorithm
    def calcDCostSmooth_dMu(self, mu):
        return 0.
 
    # Calculate integral of each column in mu for normalization
    def calcMuColNorm(self, mu, muRefNorm):
        for k in range(self.kNNMA):
            muNorm = trapz(mu[:, k], x=self.energies)
            mu[:, k] =  (mu[:, k] / muNorm) * muRefNorm[k]
        return mu
 
    # Calculate current total cost function, fill cost function array
    def calcCostFn(self, mu, t, count):
        D = np.dot(mu, t)
        costDataMatch = 0.5 * (np.linalg.norm(self.OD - D))**2
        costSparse = self.calcCostSparse(t)
        costClusterSim = self.calcCostClusterSim(mu)
        costSmooth = self.calcCostSmooth(mu)
        costTotal = ( costDataMatch + self.lambdaSparse*costSparse 
		  + self.lambdaClusterSim*costClusterSim + self.lambdaSmooth*costSmooth )
        if count > 0:
            deltaError = self.costFnArray[count, 1] - self.costFnArray[count-1, 1]
        elif count == 0:
            deltaError = -1e-13
        self.costFnArray[count, :] = np.array([costTotal, deltaError, 
	          costSparse, costClusterSim, costSmooth])
        return costTotal, deltaError


#---------------------------------------------------------------------------------------
# Calculate NNMA
#---------------------------------------------------------------------------------------
    def calcNNMA(self, initmatrices = 'Random'):
      
        print ('calculating nnma')
        
        self.initMatrices = initmatrices

        self.OD = self.stack.od.copy()
        self.energies = self.stack.ev
        self.nEnergies = self.stack.n_ev
        self.nCols = self.stack.n_cols
        self.nRows = self.stack.n_rows
        self.nPixels = self.nCols * self.nRows
    
        # Transpose optical density matrix since NNMA needs dim NxP
        self.OD = self.OD.T
        # Zero out negative values in OD matrix
        negInd = np.where(self.OD < 0.)
        if negInd: self.OD[negInd] = 0.

        self.muRecon = np.zeros((self.nEnergies, self.kNNMA))
        self.tRecon = np.zeros((self.kNNMA, self.nPixels))
        self.DRecon = np.zeros((self.nEnergies, self.nPixels))
        self.costFnArray = np.zeros((self.maxIters+1, 5))
        self.costFnArray[0, 0] = 1e99
        self.costTotal = 0.    # stores current value of total cost function
        self.deltaError = 0    # stores current value of cost function change

        # If doing cluster spectra similarity regularization, import cluster spectra
        if initmatrices == 'Cluster':
            self.muCluster = self.clusterspectra
        else:
            self.muCluster = 0.

        self.timeTaken = 0.    # stores time taken to complete NNMA analysis

        # Initialize matrices
        muInit, tInit = self.fillInitMatrices()
        self.muinit = muInit.copy()
        muCurrent = muInit
        tCurrent = tInit
        count = 0
        costCurrent, deltaErrorCurrent = self.calcCostFn(muCurrent, tCurrent, count)

        # Start NNMA
        startTime = time.time()
  
        while ((count < self.maxIters) and (self.deltaError < self.deltaErrorThreshold)):
            # Store values from previous iterations before updating
            self.tRecon = tCurrent
            self.muRecon = muCurrent
            self.costTotal = costCurrent
            self.deltaError = deltaErrorCurrent

            # Now do NNMA update
            tUpdated = self.tUpdate(muCurrent, tCurrent)
            muUpdated = self.muUpdate(muCurrent, tUpdated)
            # Zero out any negative values in t and mu
            negIndT = np.where(tUpdated < 0.)
            if negIndT: tUpdated[negIndT] = 0.
            negIndMu = np.where(muUpdated < 0.)
            if negIndMu: muUpdated[negIndMu] = 0.
            tCurrent = tUpdated
            muCurrent = muUpdated

            # Calculate cost function
            costCurrent, deltaErrorCurrent = self.calcCostFn(muCurrent, tCurrent, count)

            count = count + 1
            print ('Iteration number {0}/{1}'.format(count,self.maxIters))

        endTime = time.time()
        self.timeTaken = endTime - startTime
        if count < self.maxIters:
            self.maxIters = count
            
        self.tRecon = np.reshape(self.tRecon, (self.kNNMA, self.nCols, self.nRows), order='F')
        
        
        return 1
