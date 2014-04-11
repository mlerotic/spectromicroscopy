import numpy as np
import scipy as sp
import sys, time

import analyze

class nnma():

  def __init__(self, stkdata, kNNMA=5):

    self.stack = stkdata
    
    # Set default values for various NNMA parameters
    self.kNNMA = kNNMA            # use no. of significant PCA components, otherwise default = 5              
    self.maxIters = 10         # default no. of iterations
    self.deltaErrorThreshold = 1e-3        # default threshold for change in reconstruction error
    self.initMatrices = 'Random'  # default initialization for matrices
    
    self.lambdaSparse = 0.        # sparseness regularization param
    self.lambdaClusterSim = 0.    # cluster spectra similarity regularization param
    self.lambdaSmooth = 0.       # smoothness regularization param


#---------------------------------------------------------------------------------------
# Define some helper functions for NNMA
#---------------------------------------------------------------------------------------

  # Fill initial matrices depending on user-specified initialization method
  def fillInitMatrices(self):
    if self.initMatrices == 'Random':
      muInit = np.random.rand(self.nEnergies, self.kNNMA)
      tInit = np.random.rand(self.kNNMA, self.nPixels)
    elif self.initMatrices == "Cluster":
      muInit = analyze.clusterspectra.T 
      tInit = np.random.rand(self.kNNMA, self.nPixels)  # use SVD on mu to find t instead?
    elif self.initMatrices == "Standards":
      # Read from file???
      pass
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
      muNorm = scipy.integrate.trapz(mu[:, k], x=self.energies)
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
  def calcNNMA(self):
      
    print 'calculating nnma'

    self.OD = self.stack.od.copy()
    self.energies = self.stack.ev
    self.nEnergies = self.stack.n_ev
    self.nCols = self.stack.n_cols
    self.nRows = self.stack.n_rows
    self.nPixels = self.nCols * self.nRows
    
    # Need to do some data reshaping
    self.OD = self.OD.T
    for n in range(self.nEnergies):
      self.OD[n, :] = ( (self.OD[n, :].reshape(self.nCols, self.nRows)).T ).flatten()

    self.muRecon = np.zeros((self.nEnergies, self.kNNMA))
    self.tRecon = np.zeros((self.kNNMA, self.nPixels))
    self.DRecon = np.zeros((self.nEnergies, self.nPixels))
    self.costFnArray = np.zeros((self.maxIters+1, 5))
    self.costFnArray[0, 0] = 1e99
    self.costTotal = 0.    # stores current value of total cost function
    self.deltaError = 0    # stores current value of cost function change

    # If doing cluster spectra similarity regularization, import cluster spectra
    if self.lambdaClusterSim != 0.:
      self.muCluster = analyze.clusterspectra.T
    else:
      self.muCluster = 0.

    self.timeTaken = 0.    # stores time taken to complete NNMA analysis

    # Initialize matrices
    muInit, tInit = self.fillInitMatrices()
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
      tCurrent = tUpdated
      muCurrent = muUpdated

      # Calculate cost function
      costCurrent, deltaErrorCurrent = self.calcCostFn(muCurrent, tCurrent, count)

      count = count + 1
      print count

    endTime = time.time()
    self.timeTaken = endTime - startTime
    if count < self.maxIters:
      self.maxIters = count
