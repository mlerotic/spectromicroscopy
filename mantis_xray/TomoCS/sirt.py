# 
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
# 
#   Copyright (C) 2016 Mirna Lerotic, 2nd Look
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

import os
import numpy as np
#from skimage.morphology._greyreconstruct import reconstruction_loop

M_PI = 3.14159265358979323846264338327


#----------------------------------------------------------------------
def preprocessing(ry, rz, num_pixels, center, gridx, gridy) :

    for i in range(ry):
        gridx[i] = -ry/2.+i
    
    for i in range(rz):
        gridy[i] = -rz/2.+i

    mov = float(num_pixels)/2.0-center
    if (mov-np.ceil(mov) < 1e-2):
        mov += 1e-2
    
    return mov, gridx, gridy


#----------------------------------------------------------------------
def calc_quadrant(theta_p) :

    if ((theta_p >= 0 and theta_p < M_PI/2) or (theta_p >= M_PI and theta_p < 3*M_PI/2)) :
        quadrant = True
    
    else :
        quadrant = False
    
    return quadrant

#----------------------------------------------------------------------
def calc_coords(ry, rz, xi, yi, sin_p, cos_p,
                gridx, gridy, coordx, coordy):


    srcx = xi*cos_p-yi*sin_p
    srcy = xi*sin_p+yi*cos_p
    detx = -xi*cos_p-yi*sin_p
    dety = -xi*sin_p+yi*cos_p    

    slope = (srcy-dety)/(srcx-detx)
    islope = 1/slope
    for n in range(ry+1): 
        coordy[n] = slope*(gridx[n]-srcx)+srcy

    for n in range(rz+1) :
        coordx[n] = islope*(gridy[n]-srcy)+srcx
    

    return coordx, coordy

#----------------------------------------------------------------------
def trim_coords( ry, rz, coordx, coordy, 
                 gridx, gridy, ax, ay, 
                 bx, by):

    asize = 0
    bsize = 0
    for n in range(rz+1): 
        if (coordx[n] > gridx[0]) :
            if (coordx[n] < gridx[ry]) :
                ax[asize] = coordx[n]
                ay[asize] = gridy[n]
                asize +=1
            
    for n in range(ry+1): 
        if (coordy[n] > gridy[0]) :
            if (coordy[n] < gridy[rz]) :
                bx[bsize] = gridx[n]
                by[bsize] = coordy[n]
                bsize +=1

    return asize, ax, ay, bsize, bx, by


#----------------------------------------------------------------------
def sort_intersections( ind_condition, asize, ax, ay, 
                        bsize, bx, by, 
                        coorx, coory) :

    i=0
    j=0
    k=0

    while (i<asize and j<bsize):
        a_ind = i if ind_condition else (asize-1-i)
        if (ax[a_ind] < bx[j]) :
        
            coorx[k] = ax[a_ind]
            coory[k] = ay[a_ind]
            i += 1
            k += 1
        
        else :
        
            coorx[k] = bx[j]
            coory[k] = by[j]
            j += 1
            k +=1
        
    
    while (i < asize) :
        a_ind = i if ind_condition else (asize-1-i)
        coorx[k] = ax[a_ind]
        coory[k] = ay[a_ind]
        i += 1
        k += 1
    
    while (j < bsize) :
        coorx[k] = bx[j]
        coory[k] = by[j]
        j += 1
        k += 1
    
    
    csize = asize+bsize
    
    return csize


#----------------------------------------------------------------------
def calc_dist( ry, rz, csize, coorx, coory, indi, dist):

    for n in range(csize-1): 
    
        diffx = coorx[n+1]-coorx[n]
        diffy = coory[n+1]-coory[n]
        dist[n] = np.sqrt(diffx*diffx+diffy*diffy)
        midx = (coorx[n+1]+coorx[n])/2
        midy = (coory[n+1]+coory[n])/2
        x1 = midx+ry/2.
        x2 = midy+rz/2.
        i1 = (int)(midx+ry/2.)
        i2 = (int)(midy+rz/2.)
        indx = i1-(i1>x1)
        indy = i2-(i2>x2)
        indi[n] = indy+(indx*rz)
    
    return indi, dist



#----------------------------------------------------------------------
def calc_simdata( p, s, c, ry, rz, num_slices, num_pixels, 
                  csize, indi, dist, model, simdata):

    index_model = s*ry*rz
    index_data = c+s*num_pixels+p*num_slices*num_pixels

    for n in range(csize-1): 
        simdata[index_data] += model[indi[n]+index_model]*dist[n]
        #simdata[index_data] += model[p,s,c]*dist[n]
    
    return simdata




    
#----------------------------------------------------------------------   
#     Calculate sirt
#     tomo : ndarray
#         3D tomographic data.
#     theta : array
#         Projection angles in radian.
#     recon : ndarray, optional
#         Initial values of the reconstruction object.
#     num_gridx, num_gridy : int, optional
#         Number of pixels along x- and y-axes in the reconstruction grid.
#     num_iter : int, optional
#         Number of algorithm iterations performed.
def calculate_sirt(tomo, theta, num_iter):
    
    dx, dy, dz = tomo.shape
    print ('Calculate SIRT reconstruction')
    center = np.ones(dy, dtype='float32') * dz / 2.
    ngridx = dz
    ngridy = dz
    #tomo = -np.log(tomo)
    #recon = 1e-6 * np.ones((dy, ngridx, ngridy), dtype='float32')
    recon = 1e-6 * np.ones((dy*ngridx*ngridy), dtype='float32')
    
    print ('tomoshape', tomo.shape)
    
        
    data = tomo
    
    gridx = np.zeros((ngridx+1), dtype='float32')
    gridy = np.zeros((ngridy+1), dtype='float32')

    coordx = np.zeros((ngridy+1), dtype='float32')
    coordy = np.zeros((ngridx+1), dtype='float32') 
    
    ax = np.zeros((ngridx+ngridy), dtype='float32')    
    ay = np.zeros((ngridx+ngridy), dtype='float32') 
    bx = np.zeros((ngridx+ngridy), dtype='float32')       
    by = np.zeros((ngridx+ngridy), dtype='float32') 

    coorx = np.zeros((ngridx+ngridy), dtype='float32') 
    coory = np.zeros((ngridx+ngridy), dtype='float32') 
    dist = np.zeros((ngridx+ngridy), dtype='float32') 
    indi = np.zeros((ngridx+ngridy), dtype='int') 

          
    for i in range(num_iter): 
        print ('Iteration ', i)
        simdata = np.zeros((dx*dy*dz), dtype='float32')

        #For each slice
        for s in range(dy): 
            print ('Slice', s)

            mov, gridx, gridy = preprocessing(ngridx, ngridy, dz, center[s], gridx, gridy) 

            sum_dist = np.zeros((ngridx*ngridy), dtype='float32') 
            update = np.zeros((ngridx*ngridy), dtype='float32') 
            
            # For each projection angle 
            for p in range(dx): 
                # Calculate the sin and cos values 
                # of the projection angle and find
                # at which quadrant on the cartesian grid.
                theta_p = np.fmod(theta[p], 2*M_PI)
                quadrant = calc_quadrant(theta_p)
                sin_p = np.sin(theta_p)
                cos_p = np.cos(theta_p)
    
                # For each detector pixel 
                for d in range(dz): 
                
                    # Calculate coordinates
                    xi = -1e6
                    yi = -(dz-1)/2.0+d+mov
                    coordx, coordy = calc_coords(
                        ngridx, ngridy, xi, yi, sin_p, cos_p, gridx, gridy, 
                        coordx, coordy)

                    # Merge the (coordx, gridy) and (gridx, coordy)
                    asize, ax, ay, bsize, bx, by = trim_coords(
                        ngridx, ngridy, coordx, coordy, gridx, gridy, 
                        ax, ay, bx, by)

                    # Sort the array of intersection points (ax, ay) and
                    # (bx, by). The new sorted intersection points are 
                    # stored in (coorx, coory). Total number of points 
                    # are csize.
                    csize = sort_intersections(
                        quadrant, asize, ax, ay, bsize, bx, by, 
                        coorx, coory)

                    # Calculate the distances (dist) between the 
                    # intersection points (coorx, coory). Find the 
                    # indices of the pixels on the reconstruction grid.
                    indi, dist = calc_dist(
                        ngridx, ngridy, csize, coorx, coory, 
                        indi, dist)

                    # Calculate simdata 
                    simdata = calc_simdata(p, s, d, ngridx, ngridy, dy, dz,
                        csize, indi, dist, recon, simdata)

                    # Calculate dist*dist
                    sum_dist2 = 0.0
                    for n in range(csize-1): 
                        sum_dist2 += dist[n]*dist[n]
                        sum_dist[indi[n]] += dist[n]
                    

                    # Update
                    if (sum_dist2 != 0.0) :
                    
                        ind_data = d+s*dz+p*dy*dz
                        upd = (data[p,s,d]-simdata[ind_data])/sum_dist2
                        for n in range(csize-1): 
                            update[indi[n]] += upd*dist[n]
                
            m = 0
            for n in range(ngridx*ngridy):
                if (sum_dist[n] != 0.0) :
                    ind_recon = s*ngridx*ngridy
                    recon[m+ind_recon] += update[m]/sum_dist[n]
                m+=1
                           

    recon = np.reshape(recon,(dy, ngridx, ngridy), order='C')
    
    print ('SIRT Reconstruction Done')
    return recon
        
            