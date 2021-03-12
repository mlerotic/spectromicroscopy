# 
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
# 
#   Copyright (C) 2015 Mirna Lerotic, 2nd Look
#   http://2ndlookconsulting.com
#   License: GNU GPL v3
#
#   Originally part of tomotv, code licensed under BSD license.
#   Website: https://github.com/emmanuelle/tomo-tv
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

from __future__ import absolute_import

import numpy as np
from scipy import sparse

from .tv_denoising import tv_denoise_fista
from .projections import back_projection, projection

# ------------------ Computing energies ---------------------------

def tv_norm(im):
    """Compute the (isotropic) TV norm of an image"""
    grad_x1 = np.diff(im, axis=0)
    grad_x2 = np.diff(im, axis=1)
    return np.sqrt(grad_x1[:, :-1]**2 + grad_x2[:-1, :]**2).sum()

def tv_norm_anisotropic(im):
    """Compute the anisotropic TV norm of an image"""
    grad_x1 = np.diff(im, axis=0)
    grad_x2 = np.diff(im, axis=1)
    return np.abs(grad_x1[:, :-1]).sum() + np.abs(grad_x2[:-1, :]).sum()

# ------------------ Proximal iterators ----------------------------

def fista_tv(y, beta, niter, H, verbose=0, mask=None):
    """
    TV regression using FISTA algorithm
    (Fast Iterative Shrinkage/Thresholding Algorithm)

    Parameters
    ----------

    y : ndarray of floats
        Measures (tomography projection). If H is given, y is a column
        vector. If H is not given, y is a 2-D array where each line
        is a projection along a different angle

    beta : float
        weight of TV norm

    niter : number of forward-backward iterations to perform

    H : sparse matrix
        tomography design matrix. Should be in csr format.

    mask : array of bools

    Returns
    -------

    res : list
        list of iterates of the reconstructed images

    energies : list
        values of the function to be minimized at the different
        iterations. Its values should be decreasing.

    Notes
    -----
    This algorithm minimizes iteratively the energy

    E(x) = 1/2 || H x - y ||^2 + beta TV(x) = f(x) + beta TV(x)

    by forward - backward iterations:

    u_n = prox_{gamma beta TV}(x_n - gamma nabla f(x_n)))
    t_{n+1} = 1/2 * (1 + sqrt(1 + 4 t_n^2))
    x_{n+1} = u_n + (t_n - 1)/t_{n+1} * (u_n - u_{n-1})

    References
    ----------

    - A. Beck and M. Teboulle (2009). A fast iterative
      shrinkage-thresholding algorithm for linear inverse problems.
      SIAM J. Imaging Sci., 2(1):183-202.

    - Nelly Pustelnik's thesis (in French),
      http://tel.archives-ouvertes.fr/tel-00559126_v4/
      Paragraph 3.3.1-c p. 69 , FISTA

    """
    n_meas, n_pix = H.shape
    if mask is not None:
        l = len(mask)
    else:
        l = int(np.sqrt(n_pix))
    n_angles = n_meas / l
    Ht = sparse.csr_matrix(H.transpose())
    x0 = np.zeros(n_pix)[:, np.newaxis]
    res, energies = [], []
    gamma = .9/ (l * n_angles)
    x = x0
    u_old = np.zeros((l, l))
    t_old = 1
    for i in range(niter):
        if verbose:
            print(i)
        eps = 1.e-4
        err = H * x - y
        back_proj = Ht * err
        tmp = x - gamma * back_proj
        if mask is not None:
            tmp2d = np.zeros((l, l))
            tmp2d[mask] = tmp.ravel()
        else:
            tmp2d = tmp.reshape((l, l))
        u_n = tv_denoise_fista(tmp2d,
                weight=beta*gamma, eps=eps)
        t_new = (1 + np.sqrt(1 + 4 * t_old**2))/2.
        t_old = t_new
        x = u_n + (t_old - 1)/t_new * (u_n - u_old)
        u_old = u_n
        res.append(x)
        data_fidelity_err = 1./2 * (err**2).sum()
        tv_value = beta * tv_norm(x)
        energy = data_fidelity_err + tv_value
        energies.append(energy)
        if mask is not None:
            x = x[mask][:, np.newaxis]
        else:
            x = x.ravel()[:, np.newaxis]
    return res, energies


def ista_tv(y, beta, niter, H=None):
    """
    TV regression using ISTA algorithm
    (Iterative Shrinkage/Thresholding Algorithm)

    Parameters
    ----------

    y : ndarray of floats
        Measures (tomography projection). If H is given, y is a column
        vector. If H is not given, y is a 2-D array where each line
        is a projection along a different angle

    beta : float
        weight of TV norm

    niter : number of forward-backward iterations to perform

    H : sparse matrix or None
        tomography design matrix. Should be in csr format. If H is none,
        the projections as well as the back-projection are computed by
        a direct method, without writing explicitely the design matrix.

    Returns
    -------

    res : list
        list of iterates of the reconstructed images

    energies : list
        values of the function to be minimized at the different
        iterations. Its values should be decreasing.

    Notes
    -----

    This algorithm minimizes iteratively the energy

    E(x) = 1/2 || H x - y ||^2 + beta TV(x)

    by simple forward - backward iterations:

    x_{n + 1} = prox_{gamma beta TV(.)} (x_n - gamma nabla f(x_n))

    where f(x) = 1/2 || H x - y ||^2

    References
    ----------

    - Proximal Splitting Methods in Signal Processing, P. Combettes
      and J.-C. Pesquet, Fixed-Point Algorithms for Inverse Problems
      in Science and Engineering, p. 185 (2011). Algorithm 10.3 with
      lamba_n = 1.

    - Nelly Pustelnik's thesis (in French),
      http://tel.archives-ouvertes.fr/tel-00559126_v4/
      Paragraph 3.3.1-c p. 68 , ISTA

    """
    if H is None:
        method = 'direct'
        n_angles, l = y.shape
        n_pix = l ** 2
    else:
        method = 'matrix'
        n_angles, l = y.shape
        n_meas, n_pix = H.shape
        l = int(np.sqrt(n_pix))
        n_angles = n_meas / l
    if method == 'matrix':
        Ht = sparse.csr_matrix(H.transpose())
    x0 = np.zeros((l, l))
    res, energies = [], []
    # l * n_angles is the Lipschitz constant of Ht H
    gamma = .9/ (l * n_angles)
    x = x0
    for i in range(niter):
        eps = 1.e-4
        # Forward part
        if method == 'matrix':
            x = x.ravel()[:, np.newaxis]
            err = H * x - y
            back_proj = Ht * err
            tmp = x - gamma * back_proj
            tmp = tmp.reshape((l, l))
        else:
            err = projection(x, n_angles) - y
            back_proj = back_projection(err)
            tmp = x - gamma * back_proj
        # backward: TV prox
        x = tv_denoise_fista(tmp, weight=beta*gamma, eps=eps)
        res.append(x)
        # compute the energy
        data_fidelity_err = 1./2 * (err**2).sum()
        tv_value = beta * tv_norm(x)
        energy = data_fidelity_err + tv_value 
        energies.append(energy)
    return res, energies

def gfb_tv(y, beta, niter, H=None, val_min=0, val_max=1, x0=None,
           stop_tol=1.e-4, nonnegconst=1):
    """
    TV regression + interval constraint using the generalized
    forward backward splitting (GFB).

    Parameters
    ----------

    y : ndarray of floats
        Measures (tomography projection). If H is given, y is a column
        vector. If H is not given, y is a 2-D array where each line
        is a projection along a different angle

    beta : float
        weight of TV norm

    niter : number of forward-backward iterations to perform

    H : sparse matrix or None
        tomography design matrix. Should be in csr format. If H is none,
        the projections as well as the back-projection are computed by
        a direct method, without writing explicitely the design matrix.

    val_min, val_max: floats
        We impose that the image values are in [val_min, val_max]

    x0 : ndarray of floats, optional (default is None)
        Initial guess

    Returns
    -------

    res : list
        list of iterates of the reconstructed images

    energies : list
        values of the function to be minimized at the different
        iterations. Its values should be decreasing.

    Notes
    -----

    This algorithm minimizes iteratively the energy

    E(x) = 1/2 || H x - y ||^2 + beta TV(x) + i_C(x)

    where TV(.) is the total variation pseudo-norm and
    i_C is the indicator function of the convex set [val_min, val_max].

    The algorithm used the generalized forward-backward scheme

    z1_{n + 1} = z1_n - x_n +
                prox_{2 gamma beta TV(.)} (2*x_n - z1_n - gamma nabla f(x_n))
    z2_{n+1} =  z1_n - x_n +
                prox_{i_C(.)}(2*x_n - z2_n - gamma nabla f(x_n)

    where f(x) = 1/2 || H x - y ||^2

    This method can in fact be used for other sums of non-smooth functions
    for which the prox operator is known.

    References
    ----------
    Hugo Raguet, Jalal M. Fadili and Gabriel Peyre, Generalized
    Forward-Backward Splitting Algorithm, preprint arXiv:1108.4404v2, 2011.

    See also
    http://www.ceremade.dauphine.fr/~peyre/numerical-tour/tours/inverse_9b_gfb/
    """
    n_angles, l = y.shape
    n_meas, n_pix = H.shape
    l = int(np.sqrt(n_pix))
    n_angles = n_meas / l
    Ht = sparse.csr_matrix(H.transpose())
    if x0 is None:
        x0 = np.zeros((l, l))
    z_1 = np.zeros((l**2, 1))
    z_2 = np.zeros((l**2, 1))
    res, energies = [], []
    # l * n_angles is the Lipschitz constant of Ht H
    gamma = 2 * .9/ (l * n_angles)
    x = x0
    energy = np.inf
    for i in range(niter):
        eps = 1.e-4
        # Forward part
        x = x.ravel()[:, np.newaxis]
        err = H * x - y
        back_proj = Ht * err
        # backward: TV and i_c proxs
        # TV part
        tmp_z_1 = 2 * x - z_1 - gamma * back_proj
        tmp_z_1 = tmp_z_1.reshape((l, l))
        z_1 = z_1 + tv_denoise_fista(tmp_z_1, weight=2 * beta * gamma,
                eps=eps).ravel()[:, np.newaxis] - x
        # Projection on the interval
        tmp_z_2 = 2 * x - z_2 - gamma * back_proj
        if nonnegconst == 1: 
            tmp_z_2[tmp_z_2 < val_min] = val_min
            #print 'Non negative constraint imposed.'
        tmp_z_2[tmp_z_2 > val_max] = val_max
        z_2 = z_2 - x + tmp_z_2
        # update x: average of z_i
        x = (0.5 * (z_1 + z_2)).reshape(l, l)
        res.append(x)
        # compute the energy
        data_fidelity_err = 1./2 * (err**2).sum()
        tv_value = beta * tv_norm(x)
        energy = data_fidelity_err + tv_value
        energies.append(energy)
        # stop criterion
        if i>2 and np.abs(energy - energies[-2]) < stop_tol*energies[1]:
            break
    return res, energies


def gfb_tv_weng(y, beta, niter, H=None, val_min=0, val_max=1, x0=None,
           stop_tol=1.e-4, xb=None, xa=None, beta2=None, nonnegconst=1):
    """
    TV regression + interval constraint using the generalized
    forward backward splitting (GFB) with adjacent energies TV regularization 

    Parameters
    ----------

    y : ndarray of floats
        Measures (tomography projection). If H is given, y is a column
        vector. If H is not given, y is a 2-D array where each line
        is a projection along a different angle

    beta : float
        weight of TV norm

    niter : number of forward-backward iterations to perform

    H : sparse matrix or None
        tomography design matrix. Should be in csr format. If H is none,
        the projections as well as the back-projection are computed by
        a direct method, without writing explicitely the design matrix.

    val_min, val_max: floats
        We impose that the image values are in [val_min, val_max]

    x0 : ndarray of floats, optional (default is None)
        Initial guess at energy Em
        
    xb : ndarray of floats, optional (default is None)
        Initial guess at adjacent energy Em-1
        
    xa : ndarray of floats, optional (default is None)
        Initial guess at adjacent energy Em+1

    Returns
    -------

    res : list
        list of iterates of the reconstructed images

    CSenergies : list
        values of the function to be minimized at the different
        iterations. Its values should be decreasing.

    Notes
    -----

    This algorithm minimizes iteratively the CSenergy

    E(x) = 1/2 || H x - y ||^2 + beta TV(x) + +beta2 TV(x,xa,xb) + i_C(x)

    where TV(.) is the total variation pseudo-norm and
    i_C is the indicator function of the convex set [val_min, val_max].

    The algorithm used the generalized forward-backward scheme

    z1_{n + 1} = z1_n - x_n +
                prox_{2 gamma beta TV(.)} (2*x_n - z1_n - gamma nabla f(x_n))
    z2_{n+1} =  z1_n - x_n +
                prox_{i_C(.)}(2*x_n - z2_n - gamma nabla f(x_n))
                
    z3_{n + 1} = z1_n - x_n + 0.1 (xa+xb) + 
                prox_{2 gamma beta2 TV(.)} (2*x_n - z1_n - gamma nabla f(x_n))
                


    where f(x) = 1/2 || H x - y ||^2

    This method can in fact be used for other sums of non-smooth functions
    for which the prox operator is known.

    References
    ----------
    Hugo Raguet, Jalal M. Fadili and Gabriel Peyre, Generalized
    Forward-Backward Splitting Algorithm, preprint arXiv:1108.4404v2, 2011.

    See also
    http://www.ceremade.dauphine.fr/~peyre/numerical-tour/tours/inverse_9b_gfb/
    """
    n_angles, l = y.shape
    n_meas, n_pix = H.shape
    l = int(np.sqrt(n_pix))
    n_angles = n_meas / l
    Ht = sparse.csr_matrix(H.transpose())
    if x0 is None:
        x0 = np.zeros((l, l))
    z_1 = np.zeros((l**2, 1))
    z_2 = np.zeros((l**2, 1))
    z_3 = np.zeros((l**2, 1))
    res, energies = [], []
    xb = xb.ravel()[:, np.newaxis]
    xa = xa.ravel()[:, np.newaxis]
    if beta2 == None:
        #beta2=beta
        beta2=0
    # l * n_angles is the Lipschitz constant of Ht H
    gamma = 2 * .9/ (l * n_angles)
    x = x0
    energy = np.inf
    for i in range(niter):
        eps = 1.e-4
        # Forward part
        x = x.ravel()[:, np.newaxis]
        err = H * x - y
        back_proj = Ht * err
        # backward: TV and i_c proxs
        # TV part
        tmp_z_1 = 2 * x - z_1 - gamma * back_proj
        tmp_z_1 = tmp_z_1.reshape((l, l))
        z_1 = z_1 + tv_denoise_fista(tmp_z_1, weight=2 * beta * gamma,
                eps=eps).ravel()[:, np.newaxis] - x
        # Projection on the interval
        tmp_z_2 = 2 * x - z_2 - gamma * back_proj
        if nonnegconst:
            tmp_z_2[tmp_z_2 < val_min] = val_min
        tmp_z_2[tmp_z_2 > val_max] = val_max
        z_2 = z_2 - x + tmp_z_2
        if beta2 != 0:
            #Use images taken at energy before and after the current one
            tmp_z_3 = 2 * x - z_3 - gamma * back_proj
            tmp_z_3 = (tmp_z_3+(xa+xb)*.1).reshape((l, l))
            z_3 = z_3 + tv_denoise_fista(tmp_z_3, weight=2 * beta2 * gamma,
                    eps=eps).ravel()[:, np.newaxis] - x       
            # update x: average of z_i
            x = (0.3333 * (z_1 + z_2 + z_3)).reshape(l, l)
        else:
            x = (0.5 * (z_1 + z_2)).reshape(l, l)
        res.append(x)
        # compute the energy
        data_fidelity_err = 1./2 * (err**2).sum()
        tv_value = beta * tv_norm(x)
        energy = data_fidelity_err + tv_value
        energies.append(energy)
        # stop criterion
        if i>2 and np.abs(energy - energies[-2]) < stop_tol*energies[1]:
            break
    
    return res, energies

def gfb_tv_local(y, beta, niter, mask_pix, mask_reg, H=None,
                                val_min=0, val_max=1, x0=None):
    """
    TV regression + interval constraint using the generalized
    forward backward splitting (GFB), in local tomography mode.

    Parameters
    ----------

    y : ndarray of floats
        Measures (tomography projection). If H is given, y is a column
        vector. If H is not given, y is a 2-D array where each line
        is a projection along a different angle

    beta : float
        weight of TV norm

    niter : number of forward-backward iterations to perform

    mask_pix: ndarray of bools
        Domain where pixels are reconstructed (typically, the disk
        inside a square).

    mask_reg: ndarray of bools
        Domain where the spatial regularization is performed

    H : sparse matrix or None
        tomography design matrix. Should be in csr format. If H is none,
        the projections as well as the back-projection are computed by
        a direct method, without writing explicitely the design matrix.

    val_min, val_max: floats
        We impose that the image values are in [val_min, val_max]

    x0 : ndarray of floats, optional (default is None)
        Initial guess

    Returns
    -------

    res : list
        list of iterates of the reconstructed images

    energies : list
        values of the function to be minimized at the different
        iterations. Its values should be decreasing.

    Notes
    -----

    This algorithm minimizes iteratively the energy

    E(x) = 1/2 || H x - y ||^2 + beta TV(x) + i_C(x)

    where TV(.) is the total variation pseudo-norm and
    i_C is the indicator function of the convex set [val_min, val_max].

    The algorithm used the generalized forward-backward scheme

    z1_{n + 1} = z1_n - x_n +
                prox_{2 gamma beta TV(.)} (2*x_n - z1_n - gamma nabla f(x_n))
    z2_{n+1} =  z1_n - x_n +
                prox_{i_C(.)}(2*x_n - z2_n - gamma nabla f(x_n)

    where f(x) = 1/2 || H x - y ||^2

    This method can in fact be used for other sums of non-smooth functions
    for which the prox operator is known.

    References
    ----------
    Hugo Raguet, Jalal M. Fadili and Gabriel Peyre, Generalized
    Forward-Backward Splitting Algorithm, preprint arXiv:1108.4404v2, 2011.

    See also
    http://www.ceremade.dauphine.fr/~peyre/numerical-tour/tours/inverse_9b_gfb/
    """
    mask_reg = mask_reg[mask_pix]
    n_meas, n_pix = H.shape
    l = len(mask_pix)
    n_angles = n_meas / l
    Ht = sparse.csr_matrix(H.transpose())
    z_1 = np.zeros((n_pix, 1))
    z_2 = np.zeros((n_pix, 1))
    res, energies = [], []
    # l * n_angles is the Lipschitz constant of Ht H
    gamma = 2 * .5/ (l * n_angles)
    x0 = np.zeros(n_pix)[:, np.newaxis]
    x = x0
    for i in range(niter):
        eps = 1.e-4
        # Forward part
        err = H * x - y
        back_proj = Ht * err
        grad_descent = x - gamma * back_proj
        # backward: TV and i_c proxs
        # TV part
        tmp_z_1 = 2 * x - z_1 - gamma * back_proj
        tmp_z_1_2d = np.zeros((l, l))
        tmp_z_1_2d[mask_pix] = tmp_z_1.ravel()
        z_1 = z_1 + tv_denoise_fista(tmp_z_1_2d, weight=2 * beta * gamma,
                eps=eps)[mask_pix][:, np.newaxis] - x
        # Projection on the interval
        tmp_z_2 = 2 * x - z_2 - gamma * back_proj
        tmp_z_2[tmp_z_2 < val_min] = val_min
        tmp_z_2[tmp_z_2 > val_max] = val_max
        z_2 = z_2 - x + tmp_z_2
        # update x: average of z_i
        x = (0.5 * (z_1 + z_2))
        x[~mask_reg] = grad_descent[~mask_reg]
        tmp = np.zeros((l, l))
        tmp[mask_pix] = x.ravel()
        res.append(tmp)
    return res, energies
