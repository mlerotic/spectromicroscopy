import h5py
import numpy as np
import scipy.interpolate
import scipy.integrate
import matplotlib
matplotlib.use('Agg')  # used to save plots without needing X display
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import time

matplotlib.rcParams['pdf.fonttype'] = 42

#*************************************************************************************************
# Read in HDF5 file and extract transmission data, energy, white spectrum
#*************************************************************************************************
f = h5py.File("/home/rmak/scratch/testSparse/sperm/spermAligned2WithI0LimE.hdf5", 'r')

if "exchange" in f:
  exchangeGrp = f["exchange"]
  transmData = (exchangeGrp["data"])[...]  # reading transmission data from file
  print("transmData.shape = ", transmData.shape)
  energies = (exchangeGrp["energy"])[...]

if "spectromicroscopy" in f:
  spectrGrp = f["spectromicroscopy"]
  if "normalization" in spectrGrp:
    normGrp = spectrGrp["normalization"]
    whiteSpect = (normGrp["white_spectrum"])[...]  # reading I_0 data from file
    whiteSpectE = (normGrp["white_spectrum_energy"])[...]

f.close()


#*************************************************************************************************
# Some data pre-processing
#*************************************************************************************************
nEnergies = len(energies)
nCols = transmData.shape[0]  # x dimension of image
nRows = transmData.shape[1]  # y dimension of image
nPixels = nRows * nCols 
D = np.zeros((nEnergies, nPixels))

# Check for negative values in transmission data; if found, set to small value
negTransmInd = np.where(transmData <= 0)
if negTransmInd: transmData[negTransmInd] = 0.01

# Interpolate I0 energy measurements to transmission data energies
energiesI0 = energies.copy()
funcI0Interp = (scipy.interpolate.interp1d(whiteSpectE, whiteSpect, kind="cubic",
         bounds_error=False, fill_value=0.0))
I0 = funcI0Interp(energiesI0)

# Flatten 3D data into 2D, then convert transmission to optical density
for n in range(nEnergies):
  D[n, :] = transmData[:, :, n].flatten()  # does it matter if we use data[:, :, n].T (transpose)) instead?
  D[n, :] = -np.log(D[n, :] / I0[n])  

# Zero out NaN and negative values in optical density matrix
nanInd = np.where(np.isfinite(D) == False)
if nanInd: D[nanInd] = 0.
negInd = np.where(D < 0)
if negInd: D[negInd] = 0.
print("len(negInd) = ", len(negInd[0]))
print("DMin = ", np.amin(D))
print("DMax = ", np.amax(D))


#*************************************************************************************************
# NNMA calculations
#*************************************************************************************************
# Initialize mu, t
kNNMA = 5   # estimate for number of spectral components
lambdaMu = 0.
lambdaMuCluster = 0.	 # regularization param for smoothness
lambdaT = 0.

#muInit = np.genfromtxt("/home/rmak/scratch/appendixNNMATests/clusterSpectraSperm.txt")
#tInit = np.genfromtxt("/home/rmak/scratch/appendixNNMATests/tLetters.txt")
muInit = np.random.rand(nEnergies, kNNMA)
tInit = np.random.rand(kNNMA, nPixels)
#tInitColSum = np.sum(tInit, axis=0)
#for p in range(nPixels):
#  tInit[:, p] = tInit[:, p] / tInitColSum[p]

muCluster = np.genfromtxt("/home/rmak/scratch/appendixNNMATests/clusterSpectraSperm.txt")
# Integrate each cluster spectrum for use in scaling later
muClusterNorm = np.zeros(kNNMA)
for k in range(kNNMA):
  muClusterNorm[k] =  scipy.integrate.trapz(muCluster[:, k], x=energies)

countMax = 10000
count = 0
costFnOld = 1e99  # set some large number for initial cost function
# Array to keep track of cost fn as well as change in cost fn
costFnArray = np.zeros((countMax+1, 4))

mu = muInit
t = tInit

dJSmoothdMu = 0.
JSmooth = 0.

#########################################################################################################
# Function to calculate smoothness cost JSmooth,
# given by sum of second derivatives of mu wrt E
#########################################################################################################
def calc_JSmooth(E, mu):
  N = len(E)  # no. of energies
  energyInd = np.arange(0, N, 1)
  d2mudE2Sq = np.zeros(mu.shape)

  for n in energyInd:
   if n == 0:
     d2mudE2Sq[0, :] = 4 * ( mu[0, :]**2 / ((E[1]-E[0])**2*(E[2]-E[0])**2) 
     				+ mu[1, :]**2 / ((E[2]-E[1])**2*(E[1]-E[0])**2)
				+ mu[2, :]**2 / ((E[2]-E[1])**2*(E[2]-E[0])**2)
				- 2*mu[0, :]*mu[1, :] / ((E[2]-E[1])*(E[1]-E[0])**2*(E[2]-E[0]))
				+ 2*mu[0, :]*mu[2, :] / ((E[2]-E[1])*(E[1]-E[0])*(E[2]-E[0])**2)
				- 2*mu[1, :]*mu[2, :] / ((E[2]-E[1])**2*(E[1]-E[0])*(E[2]-E[0])) )
   elif n == (N - 1):
     d2mudE2Sq[N-1, :] = 4 * ( mu[N-1, :]**2 / ((E[N-1]-E[N-2])**2*(E[N-1]-E[N-3])**2) 
     				+ mu[N-2, :]**2 / ((E[N-1]-E[N-2])**2*(E[N-2]-E[N-3])**2)
				+ mu[N-3, :]**2 / ((E[N-2]-E[N-3])**2*(E[N-1]-E[N-3])**2)
				- 2*mu[N-1, :]*mu[N-2, :] / ((E[N-1]-E[N-2])**2*(E[N-2]-E[N-3])*(E[N-1]-E[N-3]))
				+ 2*mu[N-1, :]*mu[N-3, :] / ((E[N-1]-E[N-2])*(E[N-2]-E[N-3])*(E[N-1]-E[N-3])**2)
				- 2*mu[N-2, :]*mu[N-3, :] / ((E[N-1]-E[N-2])*(E[N-2]-E[N-3])**2*(E[N-1]-E[N-3])) )
   else:
     d2mudE2Sq[n, :] = 4 * ( mu[n+1, :]**2 / ((E[n+1]-E[n])**2*(E[n+1]-E[n-1])**2) 
     				+ mu[n, :]**2 / ((E[n+1]-E[n])**2*(E[n]-E[n-1])**2)
				+ mu[n-1, :]**2 / ((E[n]-E[n-1])**2*(E[n+1]-E[n-1])**2)
				- 2*mu[n+1, :]*mu[n, :] / ((E[n+1]-E[n])**2*(E[n]-E[n-1])*(E[n+1]-E[n-1]))
				+ 2*mu[n+1, :]*mu[n-1, :] / ((E[n+1]-E[n])*(E[n]-E[n-1])*(E[n+1]-E[n-1])**2)
				- 2*mu[n, :]*mu[n-1, :] / ((E[n+1]-E[n])*(E[n]-E[n-1])**2*(E[n+1]-E[n-1])) )
     
   JSmooth = np.sum(d2mudE2Sq)
   return JSmooth

#########################################################################################################
# Function to calculate derivative of JSmooth wrt mu;
# used in multiplicative update algorithm.
#########################################################################################################
def calc_dJSmoothdMu(E, mu):
  N = len(E)  # no. of energies
  energyInd = np.arange(0, N, 1)
  dJSmoothdMu = np.zeros(mu.shape)
  for n in energyInd:
    if n == 0:
      dJSmoothdMu[0, :] = (16 * ( mu[0,:] / ((E[1]-E[0])**2*(E[2]-E[0])**2)
      		- mu[1,:] / ((E[2]-E[1])*(E[1]-E[0])**2*(E[2]-E[0]))
		+ mu[2,:] / ((E[2]-E[1])*(E[1]-E[0])*(E[2]-E[0])**2) ) )
    elif n == 1:
      dJSmoothdMu[1, :] = (8 * (mu[0,:] * -2 / ((E[2]-E[1])*(E[1]-E[0])**2*(E[2]-E[0]))
      		+ mu[1,:] * ( 2 / ((E[2]-E[1])**2*(E[1]-E[0])**2) + 1/((E[2]-E[1])**2*(E[3]-E[1])**2) )
		+ mu[2,:] * ( -2/((E[2]-E[1])**2*(E[1]-E[0])*(E[2]-E[0])) + -1/((E[3]-E[2])*(E[2]-E[1])**2*(E[3]-E[1])) )
		+ mu[3,:] * 1 / ((E[3]-E[2])*(E[2]-E[1])*(E[3]-E[1])**2) ) )
    elif n == 2:
      dJSmoothdMu[2, :] = (8 * (mu[0,:] * 2 / ((E[2]-E[1])*(E[1]-E[0])*(E[2]-E[0])**2)
      		+ mu[1,:] * ( -2 / ((E[2]-E[1])**2*(E[1]-E[0])*(E[2]-E[0])) + -1/((E[3]-E[2])*(E[2]-E[1])**2*(E[3]-E[1])) )  
		+ mu[2,:] * ( 2 / ((E[2]-E[1])**2*(E[2]-E[0])**2) + 1/((E[3]-E[2])**2*(E[2]-E[1])**2) + 1/((E[3]-E[2])**2*(E[4]-E[2])**2) )
		+ mu[3, :] * (-1 / ((E[3]-E[2])**2*(E[2]-E[1])*(E[3]-E[1])) + -1/((E[4]-E[3])*(E[3]-E[2])**2*(E[4]-E[2])) )
		+ mu[4,:] * 1 / ((E[4]-E[3])*(E[3]-E[2])*(E[4]-E[2])**2) ) )
    elif n == (N - 3):
      dJSmoothdMu[n, :] = (8 * ( mu[N-5,:] * 1 / ((E[N-3]-E[N-4])*(E[N-4]-E[N-5])*(E[N-3]-E[N-5])**2)
      		+ mu[N-4,:] * ( -1/((E[N-3]-E[N-4])**2*(E[N-4]-E[N-5])*(E[N-3]-E[N-5])) + -1/((E[N-2]-E[N-3])*(E[N-3]-E[N-4])**2*(E[N-2]-E[N-4])) )
		+ mu[N-3,:] * ( 1/((E[N-3]-E[N-4])**2*(E[N-3]-E[N-5])**2) + 1/((E[N-2]-E[N-3])**2*(E[N-3]-E[N-4])**2) + 2/((E[N-2]-E[N-3])**2*(E[N-1]-E[N-3])**2) )
		+ mu[N-2, :] * (-1 / ((E[N-2]-E[N-3])**2*(E[N-3]-E[N-4])*(E[N-2]-E[N-4])) + -2/((E[N-1]-E[N-2])*(E[N-2]-E[N-3])**2*(E[N-1]-E[N-3])) )
		+ mu[N-1,:] * 2 / ((E[N-1]-E[N-2])*(E[N-2]-E[N-3])*(E[N-1]-E[N-3])**2) ) )
    elif n == (N - 2):
      dJSmoothdMu[n, :] = (8 * ( mu[N-4,:] * 1 / ((E[N-2]-E[N-3])*(E[N-3]-E[N-4])*(E[N-2]-E[N-4])**2)
      		+ mu[N-3,:] * ( -1 / ((E[N-2]-E[N-3])**2*(E[N-3]-E[N-4])*(E[N-2]-E[N-4])) + -2/((E[N-1]-E[N-2])*(E[N-2]-E[N-3])**2*(E[N-1]-E[N-3])) )
		+ mu[N-2,:] * ( 1/((E[N-2]-E[N-3])**2*(E[N-2]-E[N-4])**2) + 2/((E[N-1]-E[N-2])**2*(E[N-2]-E[N-3])**2) )
		+ mu[N-1,:] * -2 / ((E[N-1]-E[N-2])**2*(E[N-2]-E[N-3])*(E[N-1]-E[N-3])) ) )
    elif n == (N - 1):
      dJSmoothdMu[n, :] = (16 * ( mu[N-3,:] / ((E[N-1]-E[N-2])*(E[N-2]-E[N-3])*(E[N-1]-E[N-3])**2)
      		- mu[N-2,:] / ((E[N-1]-E[N-2])**2*(E[N-2]-E[N-3])*(E[N-1]-E[N-3]))
		+ mu[N-1,:] / ((E[N-1]-E[N-2])**2*(E[N-1]-E[N-3])**2) ) )
    else:
      dJSmoothdMu[n, :] = (8 * ( mu[n-2,:] / ((E[n]-E[n-1])*(E[n-1]-E[n-2])*(E[n]-E[n-2])**2) 
      		+ mu[n-1,:] * ( -1/((E[n+1]-E[n])*(E[n]-E[n-1])**2*(E[n+1]-E[n-1])) + -1/((E[n]-E[n-1])**2*(E[n-1]-E[n-2])*(E[n]-E[n-2])) ) 
      		+ mu[n,:] * ( 1/((E[n+1]-E[n])**2*(E[n]-E[n-1])**2) + 1/((E[n]-E[n-1])**2*(E[n]-E[n-2])**2) + 1/((E[n+1]-E[n])**2*(E[n+2]-E[n])**2) ) 
      		+ mu[n+1,:] * ( -1/((E[n+1]-E[n])**2*(E[n]-E[n-1])*(E[n+1]-E[n-1])) + -1/((E[n+2]-E[n+1])*(E[n+1]-E[n])**2*(E[n+2]-E[n])) ) 
      		+ mu[n+2,:] / ((E[n+2]-E[n+1])*(E[n+1]-E[n])*(E[n+2]-E[n])**2) ) )

  print("dJSmoothdMu = ", dJSmoothdMu)
  return dJSmoothdMu


#########################################################################################################
# Compute current DRecon and cost function
DRecon = np.dot(muInit, tInit)
if lambdaMu > 0.:
  JSmooth = calc_JSmooth(energies, muInit)
muDiff = muInit - muCluster
costFnNew = ( 0.5 * (np.linalg.norm(D - DRecon))**2 + (lambdaMu * JSmooth) + 
	(lambdaT * np.sum(np.sum(np.abs(tInit)))) + lambdaMuCluster * (np.linalg.norm(muDiff))**2 )
costFnChange = costFnNew - costFnOld  # should be negative (decreasing cost fn)
costFnArray[0, 0] = costFnNew
costFnArray[0, 1] = costFnChange
costFnArray[0, 2] = 0.5 * (np.linalg.norm(D - DRecon))**2
costFnArray[0, 3] = lambdaT * np.sum(np.sum(np.abs(tInit)))
costFnOld = costFnNew

## Save each iteration of mu for plotting
#colours = ['b', 'r', 'g', 'k', 'c', 'm', 'y']    # colours for plotting
#plt.figure(100)
#for k in range(kNNMA):
#  plt.clf()
#  plt.plot(energies, muInit[:, k], color=colours[k], linewidth=5.)
#  #plt.savefig("mu%s_0.pdf" %k, format='pdf')
#  plt.savefig("mu%s_0.png" %k, format='png')

# Update mu, t -- regularizing 1-norm of t to penalize non-sparseness, 
# and regularizing each column of t to sum to 1 (stochasticity)
startTime = time.time()
#while ( (count < countMax) and (costFnChange < 0) ):
while ((count < countMax) and (costFnChange < 1e-3)):
  count = count + 1
  print("count = ", count)

  tUpdateFactor = np.dot(mu.T, D) / ( np.dot(mu.T, np.dot(mu, t)) + lambdaT + 1e-9 )
  tUpdated = t * tUpdateFactor

  if lambdaMu > 0.:
    dJSmoothdMu = calc_dJSmoothdMu(energies, mu)
  
  negDJInd = np.where(dJSmoothdMu < 0)
  if negDJInd:
    print("len(negDJInd[0] = ", len(negDJInd[0]))
    #print("dJSmoothdMu[negDJInd] = ", dJSmoothdMu[negDJInd])
  
  # Calculate difference between current mu and muCluster
  muDiff = mu - muCluster
  
  muUpdateFactor = np.dot(D, tUpdated.T) / ( np.dot(mu, np.dot(tUpdated, tUpdated.T)) 
   			+ lambdaMu * dJSmoothdMu + lambdaMuCluster * 2 * muDiff + 1e-9 )
  #muUpdateFactor = np.dot(D, tUpdated.T) / (np.dot(mu, np.dot(tUpdated, tUpdated.T)) + 1e-9)
  muUpdated = mu * muUpdateFactor

  mu = muUpdated
  t = tUpdated
  
  negTInd = np.where(t < 0)
  if negTInd:
    print("len(negTInd[0]) = ", len(negTInd[0]))
    print("t[negTInd] = ", t[negTInd])
    t[negTInd] = 0.
  
  norm = True
  if norm == True:
  #  # Normalize each column of t
  #  tColSum = np.sum(t, axis=0)
  #  for p in range(nPixels):
  #    t[:, p] = t[:, p] / tColSum[p]
    
    ## Normalize integral under each mu
    #for k in range(kNNMA):
    #  muNorm = scipy.integrate.trapz(mu[:, k], x=energies)
    #  mu[:, k] = mu[:, k] / muNorm

    # Scale integral under each mu to same as cluster spectra
    for k in range(kNNMA):
      muNorm = scipy.integrate.trapz(mu[:, k], x=energies)
      mu[:, k] = (mu[:, k] / muNorm) * muClusterNorm[k]
    
    ## Normalize each column of mu
    #muColSum = np.sum(mu, axis=0)
    #for k in range(kNNMA):
    #  mu[:, k] = mu[:, k] / muColSum[k]


  ## Initialize scaling matrix
  #muScaling = np.eye(kNNMA, kNNMA)
  ##indLow = np.where(energies < 291.99)[0][-1]
  ##indHigh = np.where(energies < 292.0)[0][-1]
  #ind = np.where(energies < 292.0)[0][-1]
  #for k in range(kNNMA):
  #  #avg = np.mean(mu[indLow:indHigh, k])
  #  avg = mu[ind, k]
  #  muScaling[k, k] = 1. / avg
  ##for k in range(kNNMA):
  ##  muScaling[k, k] = 1. / refMax
  #print("muScaling = ", muScaling)
  #mu = np.dot(mu, muScaling)
  ## Need to multiply t by inverse(muScaling) to preserve cost function
  #t = np.dot(np.linalg.inv(muScaling), t)

  DRecon = np.dot(mu, t)
  negDInd = np.where(DRecon < 0)
  if negDInd:
    print("len(negDInd[0]) = ", len(negDInd[0]))
  
  if lambdaMu > 0.:
    JSmooth = calc_JSmooth(energies, mu)

  costFnNew = ( 0.5 * (np.linalg.norm(D - DRecon))**2 + (lambdaMu * JSmooth) 
  		+ (lambdaT * np.sum(np.sum(np.abs(t))))
  		+ lambdaMuCluster * (np.linalg.norm(muDiff))**2 )
  costFnChange = costFnNew - costFnOld
  print("costFnNew = ", costFnNew)
  print("costFnChange = ", costFnChange)
  costFnArray[count, 0] = costFnNew
  costFnArray[count, 1] = costFnChange
  costFnArray[count, 2] = 0.5 * (np.linalg.norm(D - DRecon))**2  # only basic cost fn
  costFnArray[count, 3] = lambdaT * np.sum(np.sum(np.abs(t)))    # only sparse part of cost fn
  costFnOld = costFnNew

  #if np.mod(count, 50) == 0:
  #  for k in range(kNNMA):
  #    plt.clf()
  #    plt.figure(100)
  #    plt.plot(energies, mu[:, k], color=colours[k], linewidth=5.)
  #    plt.savefig("mu%s_%s.png" %(k, count), format='png')

endTime = time.time()
timeTaken = endTime - startTime
print("Time taken = ", timeTaken)

# Save cost function to file
fileName = "costFn_sparse1Norm_k" + str(kNNMA) + "_it" + str(countMax) + ".txt"
np.savetxt(fileName, costFnArray)
fileName = "mu.txt"
np.savetxt(fileName, mu)
fileName = "t.txt"
np.savetxt(fileName, t)

# Print out some column norms for t, check to see how close they are to 1
for i in range(5):
  tColSum = np.sum(t[:, i])
  print("t col %s sum = %s" % (i, tColSum))
for i in range(5):
  tColSum = np.sum(t[:, -i-1])
  print("t col %s sum = %s" % (-i-1, tColSum))

for i in range(5):
  tColSum = np.sum(tInit[:, i])
  print("tInit col %s sum = %s" % (i, tColSum))
for i in range(5):
  tColSum = np.sum(tInit[:, -i-1])
  print("tInit col %s sum = %s" % (-i-1, tColSum))

# Print out some column norms for mu, check to see how close they are to 1
for i in range(kNNMA):
  muNorm = scipy.integrate.trapz(mu[:, i], x=energies)
  print("mu col %s norm = %s" % (i, muNorm))

for i in range(kNNMA):
  muNorm = scipy.integrate.trapz(mu[:, i], x=energies)
  print("muInit col %s norm = %s" % (i, muNorm))

#*************************************************************************************************
# Plot cost functions, spectra, thickness maps, etc.
#*************************************************************************************************
# Plot cost function
plt.figure(0)
plt.plot(np.arange(1, countMax+1, 1), costFnArray[1:countMax+1, 0], 'b', linewidth=3.)
ax = plt.gca()
plt.setp(ax.get_xticklabels(), fontsize=16, fontweight="bold")
plt.setp(ax.get_yticklabels(), fontsize=16, fontweight="bold")
ax.set_xlabel("Iteration number", fontsize=18, fontweight="bold")
ax.set_ylabel("Cost function", fontsize=18, fontweight="bold")
plt.savefig("costFn.pdf", format="pdf")

# Plot reconstructed spectra
colours = ['k', 'r', 'b', 'g', 'c', 'm', 'y', 'k', 'k', 'k']
for k in range(kNNMA):
  plt.figure(1 + k)
  plt.plot(energies, mu[:, k], color=colours[k], linewidth=3.)
  plt.title(r"$\mu_{%s}$" % str(k + 1), fontsize=18, fontweight="bold")
  ax = plt.gca()
  plt.setp(ax.get_xticklabels(), fontsize=16, fontweight="bold")
  plt.setp(ax.get_yticklabels(), fontsize=16, fontweight="bold")
  ax.set_xlabel(r"Energy (eV)", fontsize=18, fontweight="bold")
  ax.set_ylabel(r"Absorption $\mu$ (au)", fontsize=18, fontweight="bold")
  fileName = "mu" + str(k+1) + ".pdf"
  plt.savefig(fileName, format="pdf")

## Visualize sparseness of t
#fig = plt.figure(num=100, figsize=(16, 9))
#ax = Axes3D(fig)
#colours = ['r', 'g', 'b', 'y', 'c', 'm', 'k', 'k' 'k', 'k']
#x = np.arange(kNNMA) + 1
#P = np.arange(nPixels) + 1
#ind = np.arange(0, len(P), 1)
#for i in range(kNNMA):
#  tComp = t[i, ind]
#  ax.bar(P[ind], tComp, zs=x[i], zdir='y', color=colours[i], edgecolor=colours[i], alpha=0.8)
#ax.set_xlabel("P")
#ax.set_ylabel("k")
#ax.set_zlabel("t")
#ax.set_title(r"$\bf{t}$, $\lambda = %s$, iters = $%s$" % (lambdaT, countMax))
#plt.savefig("tVisual.pdf", format="pdf")

# Plot thickness maps
t = t.reshape(kNNMA, nCols, nRows)
for k in range(kNNMA):
  plt.figure(1 + kNNMA + k)
  plt.imshow(t[k, :, :], cmap=cm.gray)
  #plt.colorbar()
  plt.title(r"$\bf{t_{%s}}$" % str(k + 1), fontsize=18, fontweight="bold")
  ax = plt.gca()
  plt.setp(ax.get_xticklabels(), fontsize=16, fontweight="bold")
  plt.setp(ax.get_yticklabels(), fontsize=16, fontweight="bold")
  ax.set_xlabel("x", fontsize=18, fontweight="bold")
  ax.set_ylabel("y", fontsize=18, fontweight="bold")
  fileName = "t" + str(k+1) + ".pdf"
  plt.savefig(fileName, format="pdf")

plt.show()
