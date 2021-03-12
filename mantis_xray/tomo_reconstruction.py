# 
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
# 
#   Copyright (C) 2015 Mirna Lerotic, 2nd Look
#   http://2ndlookconsulting.com
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
from __future__ import absolute_import

import os
import numpy as np
import scipy as sp
from time import time

import multiprocessing
from functools import partial

try:
    from skimage.transform import iradon, radon
except:
    print ('SIRT reconstruction requires skimage library.')

import warnings
warnings.simplefilter('ignore', DeprecationWarning)

from .TomoCS.forward_backward_tv import fista_tv,gfb_tv, gfb_tv_weng
from .TomoCS.projections import build_projection_operator
from .TomoCS.util import generate_synthetic_data

from .TomoCS import sirt as st

from . import Mrc
        
#----------------------------------------------------------------------
def write_mrc(stack, mrcfn):
    
    #stackf32 = stack.astype(np.float32)
    Mrc.save(stack, mrcfn,ifExists='overwrite',calcMMM=False)
    
    
#----------------------------------------------------------------------
def load_mrc( mrcfn):

    stack = Mrc.load(mrcfn)

    return stack

#----------------------------------------------------------------------
def calc_sirt(R, theta, Rmin, Rmax, nonnegconst, maxiter, nrows):

    R = (R-Rmin)/(Rmax-Rmin)

    S1 = np.sum(R)

    At = iradon(R, theta=theta, output_size=nrows)
    S2 = np.sum(At)
    At = (At/S2)*S1
    xk = At

    for k in range(maxiter):
        t = iradon(radon(xk,theta),theta)
        #normalize
        St = np.sum(t)
        t = (t/St)*S1
        #update using (At g - At A x_k)
        #new xk = xk + difference between reconstruction At_starting - t_previuous_step
        xk = xk + At - t

        if nonnegconst == 1:
            #delete values <0 aka not real!
            xk = xk.clip(min=0)

    return xk

#----------------------------------------------------------------------
def calc_cs(R, theta, Rmin, Rmax, l, beta, initx0, nonnegconst, maxiter):

    R = (R-Rmin)/(Rmax-Rmin)

    proj = R.ravel()[:, np.newaxis]
    proj_operator = build_projection_operator(l, n_dir=len(theta), angles=theta)

    res, engs = gfb_tv(proj, beta, maxiter, H=proj_operator, x0=initx0, nonnegconst=nonnegconst)

    return res[-1]

#----------------------------------------------------------------------
class Ctomo:
    def __init__(self, stkdata):

        self.stack = stkdata
        
        self.tomorec = []


#----------------------------------------------------------------------   
# Calculate tomo  reconstruction
# Algorithm = 0 : CS reconstruction for 1 dataset 
# Algorithm = 1 : SIRT reconstruction for 1 dataset 
    def calc_tomo(self, tomodata, theta, maxiter, beta, samplethickness, algorithm = 0,
                  x0=[], comp = 0, beta2 = 0, nonnegconst = 1, nprocessors=1):
        
        tomodata = tomodata
        #Algorithm 
        if algorithm == 0:
            if nprocessors <= 1:
                self.calc_tomo_cs(tomodata.astype(np.float32), theta, maxiter, beta, samplethickness,
                              nonnegconst = nonnegconst)
            else:
                self.calc_tomo_cs_multi(tomodata.astype(np.float32), theta, maxiter, beta, samplethickness,
                              nonnegconst = nonnegconst, nprocessors=nprocessors)
        elif algorithm == 1:
            self.calc_tomo_sirt(tomodata.astype(np.float32), theta, maxiter, beta, samplethickness, 
                                nonnegconst = nonnegconst, nprocessors=nprocessors)
        else:
            self.calc_tomo_cs_tveng(tomodata.astype(np.float32), theta, maxiter, beta, samplethickness,
                                    x0, comp, beta2, nonnegconst = nonnegconst)  
             
        
#----------------------------------------------------------------------   
# Calculate tomo - CS reconstruction for 1 dataset  with multiprocessing
    def calc_tomo_cs_multi(self, tomodata, theta, maxiter, beta, samplethickness, nonnegconst=1, nprocessors=6):
        

        print ('Compressed sensing TV regression')
        #print ('Angles ', theta)
        print ('TV beta ', beta)
        print ('MAX iterations ', maxiter)
        print ("Sample thickness ", samplethickness)
        
        #Check if we have negative angles, if yes convert to 0-360deg
        for i in range(len(theta)):
            if theta[i] < 0:
                theta[i] = theta[i] + 360
        #print ('Angles in 0-360 range: ', theta)
                      
        nang = len(theta)
        #print ('Number of angles ', nang)
                
        dims = tomodata.shape
        print ('Dimensions ', dims, tomodata.dtype)
    
        stack = np.swapaxes(tomodata, 0, 1)

        theta = np.deg2rad(theta)
        
        dims = stack.shape
        ncols = dims[0]
        nrows = dims[1]
        l = dims[1]

        Rmax = np.amax(stack)
        Rmin = np.amin(stack)

        recondata = []
        projections = []
        
        initx0 = np.zeros((nrows,nrows), dtype=np.float32)
        
        t1 = time()

        R = stack[int(ncols/2),:, :].T
        R = (R-Rmin)/(Rmax-Rmin)

        proj = R.ravel()[:, np.newaxis]
        l = dims[1]
        proj_operator = build_projection_operator(l, n_dir=len(theta), angles=theta)

        res, engs = gfb_tv(proj, beta, maxiter, H=proj_operator, x0=initx0, nonnegconst = nonnegconst)

        initx0 = res[-1]

        
        for j in range(ncols):

            projections.append(stack[j,:, :].T)
        
        Rmax = np.amax(stack)
        Rmin = np.amin(stack)

        partial_calc_cs = partial(calc_cs, theta=theta, Rmin=Rmin, Rmax=Rmax,
                                  l=l, beta=beta, initx0=initx0, nonnegconst=nonnegconst,
                                  maxiter=maxiter)
        pool = multiprocessing.Pool(processes=nprocessors)

        res = pool.map(partial_calc_cs, projections)
        pool.close()
        pool.join()

        recondata = res
            
        t2 = time()
        print ("reconstruction done in %f s" %(t2 - t1))
                
        #Save the 3D tomo reconstruction to a HDF5 file
        recondata=np.array(recondata)
        dims = recondata.shape
        print ('final dims', dims)
        
        #Crop the data is sample thickness is defined
        if (samplethickness > 0) and (samplethickness<dims[2]-2):
            recondata = recondata[:,:,dims[2]/2-samplethickness/2:dims[2]/2+samplethickness/2]
               
        
        recondata = np.swapaxes(recondata, 0, 1)
        
        self.tomorec = recondata
        

#         write_mrc(recondata, 'testMantistomo.mrc')
        
        return


#----------------------------------------------------------------------
# Calculate tomo - CS reconstruction for 1 dataset
    def calc_tomo_cs(self, tomodata, theta, maxiter, beta, samplethickness, nonnegconst = 1):


        print ('Compressed sensing TV regression')
        #print ('Angles ', theta)
        print ('TV beta ', beta)
        print ('MAX iterations ', maxiter)
        print ("Sample thickness ", samplethickness)

        #Check if we have negative angles, if yes convert to 0-360deg
        for i in range(len(theta)):
            if theta[i] < 0:
                theta[i] = theta[i] + 360
        #print ('Angles in 0-360 range: ', theta)

        nang = len(theta)
        #print ('Number of angles ', nang)

        dims = tomodata.shape
        print ('Dimensions ', dims, tomodata.dtype)

        stack = np.swapaxes(tomodata, 0, 1)


        theta = np.deg2rad(theta)

        dims = stack.shape
        ncols = dims[0]
        nrows = dims[1]

        Rmax = np.amax(stack)
        Rmin = np.amin(stack)

        recondata = []

        initx0 = np.zeros((nrows,nrows), dtype=np.float32)

        t1 = time()

        for j in range(ncols):
            print (j, ' of ', ncols)
        #             for j in range(1):
            #j=ncols/2
            R = stack[j,:, :].T
            R = (R-Rmin)/(Rmax-Rmin)

            proj = R.ravel()[:, np.newaxis]
            l = dims[1]
            proj_operator = build_projection_operator(l, n_dir=len(theta), angles=theta)

            # Reconstruction
            res, engs = gfb_tv(proj, beta, maxiter, H=proj_operator, x0=initx0, nonnegconst = nonnegconst)

            recondata.append(res[-1])
            initx0 = res[-1]

        t2 = time()
        print ("reconstruction done in %f s" %(t2 - t1))

        #Save the 3D tomo reconstruction to a HDF5 file
        recondata=np.array(recondata)
        dims = recondata.shape
        print ('final dims', dims)

        #Crop the data is sample thickness is defined
        if (samplethickness > 0) and (samplethickness<dims[2]-2):
            recondata = recondata[:,:,dims[2]/2-samplethickness/2:dims[2]/2+samplethickness/2]

        recondata = np.swapaxes(recondata, 0, 1)

        self.tomorec = recondata

        #write_mrc(recondata, 'CS_reconstruction.mrc')

        return


#----------------------------------------------------------------------   
# Calculate tomo - SIRT reconstruction for 1 dataset 
    def calc_tomo_sirt(self, tomodata, theta, maxiter, beta, samplethickness, nonnegconst = 1, nprocessors=6):
        
        try:
            from skimage.transform import iradon, radon
        except:
            return

        print ('SIRT')
        #print ('Angles ', theta)
        print ('MAX iterations ', maxiter)
        print ("Sample thickness ", samplethickness)
        
        #Check if we have negative angles, if yes convert to 0-360deg
        for i in range(len(theta)):
            if theta[i] < 0:
                theta[i] = theta[i] + 360
        #print ('Angles in 0-360 range: ', theta)
                      
        nang = len(theta)
        print( 'Number of angles ', nang)
    
        tomodata = np.swapaxes(tomodata, 0, 1)
                
        dims = tomodata.shape
        print ('Dimensions ', dims, tomodata.dtype)
     
        stack = tomodata.astype(np.float32)

        #theta = np.deg2rad(theta)
        
        dims = stack.shape
        ncols = dims[0]
        nrows = dims[1]

        t1 = time()

        print ('Calculating SIRT reconstruction')
        print ('Number of iterations: ', maxiter)
        #N Iterations
        error = []
        recondata = np.zeros((ncols, nrows, nrows))

        projections = []
        for j in range(ncols):
            projections.append(stack[j,:, :])

        Rmax = np.amax(stack)
        Rmin = np.amin(stack)

        partial_calc_sirt = partial(calc_sirt, theta=theta, Rmin=Rmin, Rmax=Rmax,
                                    nonnegconst=nonnegconst, maxiter=maxiter, nrows=nrows)
        pool = multiprocessing.Pool(processes=nprocessors)

        res = pool.map(partial_calc_sirt, projections)
        pool.close()
        pool.join()

        recondata = res
            
        t2 = time()
        print ("reconstruction done in %f s" %(t2 - t1))

        recondata = np.array(recondata)

        dims = recondata.shape

        #Crop the data is sample thickness is defined
        if (samplethickness > 0) and (samplethickness<dims[2]-2):
            print ('Cropping the data')
            recondata = recondata[:,:,dims[2]/2-samplethickness/2:dims[2]/2+samplethickness/2]

        recondata = np.swapaxes(recondata, 0, 1)
        recondata = np.swapaxes(recondata, 0, 2)
           
        self.tomorec = recondata
        
        return

#----------------------------------------------------------------------
# Calculate tomo - SIRT reconstruction for 1 dataset
    def calc_tomo_sirt_singleprocessor(self, tomodata, theta, maxiter, beta, samplethickness, nonnegconst = 1):

        print ('SIRT')
        #print ('Angles ', theta)
        print ('MAX iterations ', maxiter)
        print ("Sample thickness ", samplethickness)

        #Check if we have negative angles, if yes convert to 0-360deg
        for i in range(len(theta)):
            if theta[i] < 0:
                theta[i] = theta[i] + 360
        #print ('Angles in 0-360 range: ', theta

        nang = len(theta)
        print ('Number of angles ', nang)

        tomodata = np.swapaxes(tomodata, 0, 1)

        dims = tomodata.shape
        print ('Dimensions ', dims, tomodata.dtype)

        stack = tomodata.astype(np.float32)

        #theta = np.deg2rad(theta)

        dims = stack.shape
        ncols = dims[0]
        nrows = dims[1]

        t1 = time()

        print ('Calculating SIRT reconstruction')
        print ('Number of iterations: ', maxiter)
        #N Iterations
        n = maxiter
        error = []
        recondata = np.zeros((ncols, nrows, nrows))

        Rmax = np.amax(stack)
        Rmin = np.amin(stack)

        for j in range(ncols):
        #for j in range(1):
            #j=ncols/2
            print( 'processing ', j,' of ', ncols)

            R = stack[j,:, :]
            # Normalize the sinogram
            R = (R-Rmin)/(Rmax-Rmin)

            S1 = np.sum(R)

            At = iradon(R, theta=theta, output_size=nrows)
            S2 = np.sum(At)
            At = (At/S2)*S1
            xk = At

            for k in range(n):
                t = iradon(radon(xk,theta),theta)
                #normalize
                St = np.sum(t)
                t = (t/St)*S1
                #update using (At g - At A x_k)
                #new xk = xk + difference between reconstruction At_starting - t_previuous_step
                xk = xk + At - t

                if nonnegconst == 1:
                    #delete values <0 aka not real!
                    xk = xk.clip(min=0)

            recondata[j,:,:] = xk.copy()

        t2 = time()
        print( "reconstruction done in %f s" %(t2 - t1))

        dims = recondata.shape

        #Crop the data is sample thickness is defined
        if (samplethickness > 0) and (samplethickness<dims[2]-2):
            print( 'Cropping the data')
            recondata = recondata[:,:,dims[2]/2-samplethickness/2:dims[2]/2+samplethickness/2]

        recondata = np.swapaxes(recondata, 0, 1)
        recondata = np.swapaxes(recondata, 0, 2)

        self.tomorec = recondata

        return

    

#----------------------------------------------------------------------   
# Calculate tomo - CS reconstruction for 1 dataset with energy TV regularization
    def calc_tomo_cs_tveng(self, tomodata, theta, maxiter, beta, samplethickness, initrecs, comp, beta2, nonnegconst = 1):
        

        print ('Compressed sensing TV regression with Energy Regularization')
        #print ('Angles ', theta)
        print ('TV beta ', beta)
        print ('TV beta2', beta2)
        print ('MAX iterations ', maxiter)
        print ("Sample thickness ", samplethickness)
        
        #Check if we have negative angles, if yes convert to 0-360deg
        for i in range(len(theta)):
            if theta[i] < 0:
                theta[i] = theta[i] + 360
        #print ('Angles in 0-360 range: ', theta)
                      
        nang = len(theta)
        #print ('Number of angles ', nang)
                
        dims = tomodata.shape
        #print ('Dimensions ', dims, tomodata.dtype  )
    
        stack = np.swapaxes(tomodata, 0, 1)


        theta = np.deg2rad(theta)
        
        dims = stack.shape
        ncols = dims[0]
        nrows = dims[1]
    

        recondata = []
        
        t1 = time()
        
        for j in range(ncols):
            #(print j, ' of ', ncols)
        #             for j in range(1):
            #j=ncols/2
            R = stack[j,:, :].T

 
            Rmax = np.amax(R)
            Rmin = np.amin(R)
         
            R = (R-Rmin)/(Rmax-Rmin)

#             pl.imshow(R, cmap=pl.cm.Greys_r)
#             pl.show()
                           
            proj = R.ravel()[:, np.newaxis]
            l = dims[1]
            proj_operator = build_projection_operator(l, n_dir=len(theta), angles=theta)
        
            
            if comp==0:
                xb=initrecs[comp+1][j]
                xa=initrecs[comp+2][j]
            elif comp == len(initrecs)-1:
                xb=initrecs[comp-2][j]
                xa=initrecs[comp-1][j]
            else:
                xb=initrecs[comp-1][j]
                xa=initrecs[comp+1][j]
            
            # Reconstruction            
            res, engs = gfb_tv_weng(proj, beta, maxiter, H=proj_operator, 
                                        x0=np.array(initrecs[comp][j]),
                                        xb=xb, xa=xa, beta2=beta2,
                                        nonnegconst = nonnegconst)

    
        
            recondata.append(res[-1])
            
        t2 = time()
        print ("reconstruction done in %f s" %(t2 - t1))
                
        #Save the 3D tomo reconstruction to a HDF5 file
        recondata=np.array(recondata)
        dims = recondata.shape
        #print ('final dims', dims)
        
        #Crop the data is sample thickness is defined
        if (samplethickness > 0) and (samplethickness<dims[2]-2):
            recondata = recondata[:,:,dims[2]/2-samplethickness/2:dims[2]/2+samplethickness/2]
               
        
        recondata = np.swapaxes(recondata, 0, 1)
        
        self.tomorec = recondata

        self.save_mrc(recondata, 'testMantistomo.mrc')
        
        return
    
        
#----------------------------------------------------------------------   
# Save mrc
    def save_mrc(self, path, data):  
        
        data = np.array(data, dtype=np.float32)  
        
        write_mrc(data, path)
