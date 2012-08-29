'''
Created on Oct 14, 2011
@author: Rachel Mak
'''

import numpy as np
import scipy as sp
import sys, time
from pylab import *
import csv

#----------------------------------------------------------------------
# All NNMA algorithms share the same following template (i.e., initialization, iteration)
#----------------------------------------------------------------------
class NNMATemplate():

  def initFactors(self, D, kComponents, muInit=None, tInit=None):
    print("Inside initFactors()")
    nEnergies, nPixels = D.shape
    print("nEnergies = ", nEnergies)
    print("nPixels = ", nPixels)
    # if no initial matrices are given, generate ones with random numbers
    if muInit is None:
      muInit = np.random.rand(nEnergies, kComponents)
      #mu = np.ones((nEnergies, kComponents), np.float)
    if tInit is None:
      tInit = np.random.rand(kComponents, nPixels)
      #t = np.ones((kComponents, nPixels), np.float)
    print("initial muInit = ", muInit)
    print("initial tInit = ", tInit)
    print("self.stepsizeT = ", self.stepsizeT)
    return muInit, tInit

  # Method __call__ is where the actual NNMA iteration takes place
  def __call__(self, D, energies, kNNMA, maxIters, mu0=None, verbose=True):
    print("In NNMATemplate.__call__()")
    mu, t = self.initFactors(D, kNNMA, muInit=mu0)
    #mu, t = self.initFactors(D, kNNMA, muInit=muCluster)		# give option to seed initial spectra with cluster spectra
    normD = np.linalg.norm(D)
    print("norm(D) = ", normD)
    eps = 1e-3	# tolerance level for iterations
    maxCount = maxIters	# max number of iterations
    print("maxCount = ", maxCount)
    count = 0		# counter for number of iterations
    objOld = 1e99	# diff between previous and current iterations (initially set to some large number)

    nEnergies, nPixels = D.shape
    muAllIterations = np.zeros([kNNMA, nEnergies, maxCount]) 	# array to store mu at each iteration
    costFunction = np.zeros([maxCount, 2])	# array to store [|D-D'|, delta|D-D'|]

    while True:

      print("count = ", count)
      for k in range(kNNMA):
	muAllIterations[k, :, count] = mu[:, k]		# store current reconstruction of mu

      # Save mu at each iteration to compile into movie
      self.saveMuImage(energies, mu, count)

      muUpdated, tUpdated, distUpdated = self.update(mu, t, D)	# <----- where NNMA updates take place
      if np.any(np.isnan(muUpdated)) or np.any(np.isinf(muUpdated)) or \
	 np.any(np.isnan(tUpdated)) or np.any(np.isinf(tUpdated)):
	if verbose: print "RESTART"
	mu, t = self.initFactors(D, kNNMA, muInit=m0)
	count = 0
      #count += 1
      print("dist = ", distUpdated)
      obj = distUpdated / normD 	# normalize cost function by dividing by its norm
      print("obj = ", obj)
      deltaObj = obj - objOld
      print("deltaObj = ", deltaObj)
      if deltaObj > 0:		# if cost function increased, try decreasing stepsize
        print("deltObj > 0")
        print("self.stepsizeT = ", self.stepsizeT)
        self.stepsizeT = self.stepsizeT / 2
        if self.stepsizeT < 1e-12:
          print("Algorithm has converged; stepsize = ", self.stepsizeT)
          break
      else:
        if deltaObj > -1e-6:	# if cost function decreased too slowly, increase stepsize
	  print("Increasing self.stepsizeT")
	  self.stepsizeT = self.stepsizeT * 1.5
	  print("self.stepsizeT = ", self.stepsizeT)
        count += 1
        mu = muUpdated
        t = tUpdated
        costFunction[count-1, 0] = obj
        costFunction[count-1, 1] = deltaObj 
        objOld = obj

      if (count >= maxCount) and (-eps < deltaObj <= eps):
	print("count = ", count)
	print("(eps, deltaObj) = ", eps, deltaObj)
	break
      elif (count >= maxCount):
	print("Finished specified number of iterations, but deltaObj < eps not satisfied.")
	break
      elif (self.stepsizeT < 1e-12):
        print("Stepsize is too small, breaking...")
        break

      #objFunctionWriter = csv.writer(open(dir + 'objFunction.csv', 'wb'), delimiter=',')
      #objFunctionWriter.writerow([obj, deltaObj])


    #print("mu = ", mu)
    #print("t = ", t)

    # Plot convergence of mu spectra
    figure(4)
    for k in range(kNNMA):
      subplot(kNNMA+1, 1, k+1)
      imshow(muAllIterations[k, :, :], cmap=cm.gist_rainbow)
    subplot(kNNMA+1, 1, kNNMA+1) 	# on last subplot, show cost function
    plot(arange(0, maxCount, 1), costFunction[:, 0], 'b')
    xlabel('Iteration number')
    ylabel('Cost function |D - D\'|')
    np.savetxt("costFunction.txt", costFunction)

    # Calculate and plot power spectra for mu (to quantify how smooth the mu's are)
    figure(5)
    muPower = np.zeros([nEnergies, kNNMA])
    for k in range(kNNMA):
      muPower[:,k] = np.power(np.abs(sp.fft(mu[:,k])), 2)	# need normalization?  Divide by nEnergies?
      subplot(kNNMA, 1, k)
      plot(muPower[:,k], 'r')
    
    DRecon = np.dot(mu, t)

    # Save results of reconstructed mu and t into a csv file
    #dirName = "/home/rmak/scratch/"
    #muWriter = csv.writer(open(dirName + 'mu.csv', 'wb'), delimiter=',')
    #for i in range(nEnergies):
    #  muWriter.writerow([energies[i], mu[i, :]])
    #tWriter = csv.writer(open(dir + 't.csv', 'wb'), delimiter=',')
    #  tWriter.writerow(t[row,:])

    np.savetxt('mu.txt', mu) 

    print("End of NNMATemplate.__call__()")
    return mu, t, DRecon

#----------------------------------------------------------------------
class factorizedNNMA(NNMATemplate):

  def __init__(self, muUpdate, tUpdate, dist):
    print("In factorizedNNMA.__init__()")
    self.muUpdate = muUpdate
    self.tUpdate = tUpdate
    self.dist = dist
    self.stepsizeT = 0.1	# t stepsize for gradient descent

  def update(self, mu, t, D):
    print("Inside update() in factorizedNNMA")
    muUpdated = self.muUpdate(mu, t, D)
    tUpdated = self.tUpdate(muUpdated, t, D, stepsizeT=self.stepsizeT)
    dist = self.dist(D, muUpdated, tUpdated)
    return muUpdated, tUpdated, dist

  # save mu reconstruction at each iteration so we can compile it into a movie
  def saveMuImage(self, energies, muRecon, count):
    print("In saveMuImage()")
    kComp = muRecon.shape[1]

    for k in range(kComp):
      plot(energies, muRecon[:, k], 'b')
      xlabel('Energy (eV)')
      ylabel(r'$\mu$' + str(k + 1))
      fileName = "mu" + str(k + 1) + "_" + ("%04d.png" % count)
      savefig("/home/rmak/scratch/movieImages/" + fileName)
      clf()


#----------------------------------------------------------------------
class nnma(factorizedNNMA):

  def __init__(self, stkdata, data_struct, pcaAnalz):

    self.stack = stkdata
    self.data_struct = data_struct
    self.PCAAnalz = pcaAnalz
    self.kNNMA = 5	# set some default guess value for number of chemical components
    self.maxIters = 100	# number of NNMA iterations
    self.algoNNMA = "Basic"	# set default NNMA algorithm
    self.initMatrices = 'Random'

    self.kPCA = 3	# number of principal components from PCA <-- set this as GUI button?
    self.sparsenessT = 0.8      # sparseness of t matrix, in [0, 1], 0 being most sparse

    self.NNMA = factorizedNNMA(self.muMultUpdate, self.tMultUpdate, self.frobDist)	# Basic NNMA multiplicative update with Frobenius distance measure
    self.NNMASparse = factorizedNNMA(self.muMultUpdate, self.tSparseUpdate, self.frobDist)
    #self.NNMASparse.stepsizeT = 0.1
    #print("self.NNMASparse.stepsizeT = ", self.NNMASparse.stepsizeT)

#----------------------------------------------------------------------
# Some functions for use in NNMA algorithms   

# Note: In Python, * and / are element-wise operations.
#	For matrix multiplication, use numpy.dot 
#----------------------------------------------------------------------
# Basic multiplicative update for mu
#----------------------------------------------------------------------
  def muMultUpdate(self, mu, t, D, **param):
    print("In nnma.muMultUpdate()")
    #print("mu before update = ", mu)
    updateFactor = np.dot(D, t.T) / ( np.dot(mu, np.dot(t, t.T)) + 1e-9 )
    muUpdated = mu * updateFactor
    #print("muUpdated = ", muUpdated)
    #print("t = ", t)
    return muUpdated 

#----------------------------------------------------------------------
# Basic multiplicative update for t
#----------------------------------------------------------------------
  def tMultUpdate(self, mu, t, D, **param):
    #updateFactor = np.dot(mu.T, D) / ( np.dot(mu.T, np.dot(mu, t)) + 1e-9 )
    updateFactor = np.dot(D.T, mu).T / ( np.dot(np.dot(mu.T, mu), t) + 1e-9 )
    tUpdated = t * updateFactor
    return tUpdated

#----------------------------------------------------------------------
# Frobenius distance metric
#----------------------------------------------------------------------
  def frobDist(self, D, mu, t):
    dist = np.linalg.norm(D - np.dot(mu, t))
    return dist

#----------------------------------------------------------------------
# Projection operator to enforce sparseness, from Hoyer (2004)
# x is the vector we wish to project onto non-negative space with specified constraints
# L1 and L2 norms are set to achieve desired sparseness (and non-negativity)
#----------------------------------------------------------------------
  def sparsenessProjector(self, x, L1, L2):
    print("In sparsenessProjector()")
    N = len(x)
    print("N = ", N)
    print("x.shape = ", x.shape)
    x = x.reshape(N, 1)
    # Start by projecting point to sum constraint hyperplane
    s = x + (L1 - np.sum(x)) / N
    print("s = ", s)
    zeroCoeff = []	# initially, no elements are assumed to be zero
    iters = 0
    while True:
      print("Inside sparsenessProjector while True loop")
      if len(zeroCoeff) == N:
        print("Reached maximal iterations to find s")
        break
      midPoint = np.ones((N, 1), float) * L1 / (N - len(zeroCoeff))
      print("midPoint = ", midPoint)
      midPoint[zeroCoeff] = 0
      print("midPoint = ", midPoint)
      w = s - midPoint
      a = np.sum(np.power(w, 2)) + 0j
      b = 2 * np.dot(w.T, midPoint)
      c = np.sum(np.power(midPoint, 2)) - L2
      print("a = ", a)
      print("b = ", b)
      print("c = ", c)
      alpha = (-b + (np.sqrt(np.power(b, 2) - 4 * a * c)).real) / (2 * a)	# solve quadratic eqn to find alpha such that s satisfies L2 norm constraint
      s = alpha * w + midPoint		# <--- is it alpha*w + s (as in Hoyer's Matlab code), or alpha*w + midPoint (as in Hoyer's paper)???
      print("s = ", s)

      print("alpha = ", alpha)
      if (alpha > 0) == False or np.isnan(alpha):	# alpha should be > 0
        print("breaking...; iters = ", iters)
        break

      elif np.all(s >= 0):		# if all elements of v are +ve, solution is found
        iters = iters + 1
        print("All s >= 0; iters = ", iters)
        break

      else:	# if some elements are -ve, then set them to zero and iterate again to find soln
        print("In else statement")
        iters = iters + 1
        zeroCoeff = np.where(s <= 0)[0]
        s[zeroCoeff] = 0
        sumTemp = np.sum(s)
        s = s + (L1 - sumTemp) / (N - len(zeroCoeff))
        s[zeroCoeff] = 0

    if (np.abs(s.imag)).max() > 1e-10:
      print("Error: imaginary values in v!")

    print("About to return s")
    return s


#----------------------------------------------------------------------
# Updates with sparseness constraints
# sparsenessT - sparseness of t in [0, 1]
# L1 and L2 norms set to achieve desired sparseness set by self.sparsenessT 
#----------------------------------------------------------------------
  def tSparseUpdate(self, mu, t, D, stepsizeT, L1=1., L2=1.):

    print("In tSparseUpdate()")
    print("stepsizeT = ", stepsizeT)

    # First take a step in direction of negative gradient
    tUpdated = t - stepsizeT * np.dot(mu.T, (np.dot(mu, t) - D))
    print("t = ", t)
    print("mu = ", mu)
    print("tUpdated = ", tUpdated)

    # Now project onto constraint space that satisfies desired sparseness
    L2 = 1.0	# rows of t normalized to L2
    nPixels = D.shape[1]
    L1 = (np.sqrt(nPixels) - (np.sqrt(nPixels) - 1) * self.sparsenessT) * L2
    #L1 = (np.sqrt(nPixels) - (np.sqrt(nPixels) - 1) * self.sparsenessT)
    print("L1 = ", L1)
    for i in range(self.kNNMA):
      print("i = ", i)
      tUpdated[i, :] = (self.sparsenessProjector(tUpdated[i, :].T, L1, L2)).T

    print("tUpdated.shape = ", tUpdated.shape)
    print("t = ", t)
    print("norms are: ")
    for i in range(self.kNNMA):
      print(np.linalg.norm(t[i,:]))

    return tUpdated

#----------------------------------------------------------------------
# Run the NNMA analysis
#----------------------------------------------------------------------
  def calcNNMA(self, algoName, kComponents=5, verbose=1):

    print("Doing nnma.calcNNMA:")
    print("kComponents = ", self.kNNMA)
    print("maxIters = ", self.maxIters)
    print("sparsenessT = ", self.sparsenessT)

    self.nEnergies = self.stack.n_ev
    self.nCols = self.stack.n_cols
    self.nRows = self.stack.n_rows
    self.nPixels = self.nCols * self.nRows
    self.energies = self.stack.ev.copy()

## -----------------------------------------------------------------------------------------------
#
#    # Here we're going to test using OD with a subset of original data
#    #self.kPCA = self.PCAAnalz.numsigpca 	# use num significant components calculated from PCA, instead of inputting manually as above
#    #print("self.PCAAnalz.numsigpca = ", self.PCAAnalz.numsigpca)
#    self.PCAEigenvecsSubset = self.PCAAnalz.eigenvecs[:, 0:self.kPCA]
#    print("self.PCAEigenvecsSubset.shape = ", self.PCAEigenvecsSubset.shape)
#    self.PCAImagesSubset = self.PCAAnalz.pcaimages2D[0:self.kPCA, :]
#    print("self.PCAImagesSubset.shape = ", self.PCAImagesSubset.shape)
#    self.ODSubset = np.dot(self.PCAEigenvecsSubset, self.PCAImagesSubset)
#    print("self.ODSubset.shape = ", self.ODSubset.shape)
#    D = self.ODSubset
#    self.OD = D.reshape(self.nEnergies, self.nRows, self.nCols)
#
## ----------------------------------------------------------------------------------------------

    self.OD = self.stack.absdata.copy()	# absdata in HDF5 file actually  holds transmission data (not absorption!)
    self.i0 = self.stack.i0data
    self.ODTemp = np.zeros((self.nEnergies, self.nRows, self.nCols))
    print("self.ODTemp.shape = ", self.ODTemp.shape)
    print("self.OD.shape = ", self.OD.shape)
    for i in range(self.nEnergies):		# transpose so that OD is in C-order: energy is leftmost (slowest-changing) index
      self.ODTemp[i, :, :] = self.OD[:, :, i].T
      #self.ODTemp[i, :, :] = self.OD[:, :, i]
      self.ODTemp[i, :, :] = -np.log(self.ODTemp[i, :, :] / self.i0[i])	# do this if absdata is actually transmission
    neg = np.where(self.ODTemp < 0.)	# find where OD is negative
    self.ODTemp[neg] = 0.		# and set to zero
    print("neg = ", neg)
    self.OD = self.ODTemp
    print("self.OD.shape = ", self.OD.shape)
    self.ODRecon = np.zeros((self.nEnergies, self.nCols * self.nRows), dtype=float)	# initialize matrix which holds reconstructed OD
    self.ODError = np.zeros((self.nEnergies, self.nCols * self.nRows), dtype=float)	# initialize matrix which holds difference between experimental and reconstructed OD
    self.mu = np.zeros((self.nEnergies, kComponents), dtype=float)
    self.t = np.zeros((kComponents, self.nPixels), dtype=float)

    # Rearrange 3D OD matrix into 2D
    D = np.zeros((self.nEnergies, self.nPixels), float)
    for i in range(self.nEnergies):
      D[i, :] = self.OD[i, :, :].flatten()

    if (self.kNNMA < 1) or (self.kNNMA > self.nEnergies) or (self.kNNMA > self.nPixels):
      raise ValueError, "Number of components is invalid."
	
    #self.findRoutine(algoName)
   
    # Determine what to use for initial matrix mu
    if self.initMatrices == 'Random':
      muInit = None
    elif self.initMatrices == 'Cluster':
      muInit = self.PCAAnalz.clusterspectra.T
 
    print "run %12s" % algoName
    sys.stdout.flush()
    startTime = time.time()
    print("startTime = ", startTime)

    #self.mu, self.t, DRecon = self.NNMA(D, self.energies, self.kNNMA, self.maxIters)	# calculate NNMA here by invoking the NNMATemplate.__call__ function
    self.mu, self.t, DRecon = self.NNMASparse(D, self.energies, self.kNNMA, self.maxIters, mu0=muInit)	# calculate NNMA here by invoking the NNMATemplate.__call__ function

    endTime = time.time()
    timeTaken = endTime - startTime
    print("Time taken = ", timeTaken)

    tNorm = np.linalg.norm(self.t[0,:])
    print("tNorm first row = ", tNorm)
    tNorm = np.linalg.norm(self.t[kComponents-1,:])
    print("tNorm last row = ", tNorm)

    self.ODRecon = DRecon.reshape(self.nEnergies, self.nRows, self.nCols)
    figure(0) 
    image = self.ODRecon[0, :, :]
    imshow(image)

    DError = np.abs(DRecon - D) / (D + 1e-9) * 100
    np.savetxt("D.txt", D)
    np.savetxt("DRecon.txt", DRecon)
    #np.savetxt("DError.txt", DError)
    np.savetxt("energies.txt", self.energies)
    self.ODError = DError.reshape(self.nEnergies, self.nRows, self.nCols)
    image = self.ODError[0, :, :]
    figure(1)
    imshow(image)

    mu0 = self.mu[:, 0]
    figure(2)
    plot(self.energies, self.mu, 'b')
    show()

    muImage = self.mu.T		# transpose mu matrix so that energy is displayed on x axis
    figure(3)
    maxEnergyIndex = self.nEnergies - 1
    maxMuIndex = kComponents - 1
    imshow(muImage, cmap=cm.gray, interpolation='nearest', extent=[0, maxEnergyIndex, 0, maxMuIndex]) 
    colorbar()

#----------------------------------------------------------------------
# Match up algorithm name with routine name
#----------------------------------------------------------------------
  def findRoutine(self, algoName):
    if algoName == "Basic":
      print("Inside if clause")
      print("algoName = ", algoName)
      self.NNMA.update("Updated string; inside findRoutine()")
    elif algoName == "Smooth":
      print("algoName = ", algoName)
    elif algoName == "Sparse":
      print("algoName = ", algoName)
    else:
      print("No matching NNMA routine found for ", algoName)

    print("End of findRoutine()")

#----------------------------------------------------------------------
  def printTest(self):
 
      print("self.stack.n_ev = ", self.stack.n_ev)
      print("self.stack.ev = ", self.stack.ev)
      print("self.stack.n_cols = ", self.stack.n_cols)
      print("self.stack.n_rows = ", self.stack.n_rows)

      if self.data_struct.spectromicroscopy.normalization.white_spectrum is not None:
	print("self.stack.i0data = ", self.stack.i0data)
      else:
	print("No i0data")

      print("self.stack.absdata.shape = ", self.stack.absdata.shape)

#----------------------------------------------------------------------
