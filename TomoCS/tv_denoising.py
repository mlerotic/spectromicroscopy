from __future__ import print_function
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

import numpy as np

def div(grad):
    """ Compute divergence of image gradient """
    res = np.zeros(grad.shape[1:])
    for d in range(grad.shape[0]):
        this_grad = np.rollaxis(grad[d], d)
        this_res = np.rollaxis(res, d)
        this_res[:-1] += this_grad[:-1]
        this_res[1:-1] -= this_grad[:-2]
        this_res[-1] -= this_grad[-2]
    return res

def gradient(img):
    """
    Compute gradient of an image

    Parameters
    ===========
    img: ndarray
        N-dimensional image

    Returns
    =======
    gradient: ndarray
        Gradient of the image: the i-th component along the first
        axis is the gradient along the i-th axis of the original
        array img
"""
    shape = [img.ndim, ] + list(img.shape)
    gradient = np.zeros(shape, dtype=img.dtype)
    # 'Clever' code to have a view of the gradient with dimension i stop
    # at -1
    slice_all = [0, slice(None, -1),]
    for d in range(img.ndim):
        gradient[tuple(slice_all)] = np.diff(img, axis=d)
        slice_all[0] = d + 1
        slice_all.insert(1, slice(None))
    return gradient


def _projector_on_dual(grad):
    """
    modifies in place the gradient to project it
    on the L2 unit ball
    """
    norm = np.maximum(np.sqrt(np.sum(grad**2, 0)), 1.)
    for grad_comp in grad:
        grad_comp /= norm
    return grad


def dual_gap(im, new, gap, weight):
    """
    dual gap of total variation denoising
    see "Total variation regularization for fMRI-based prediction of behavior",
    by Michel et al. (2011) for a derivation of the dual gap
    """
    im_norm = (im**2).sum()
    gx, gy = np.zeros_like(new), np.zeros_like(new)
    gx[:-1] = np.diff(new, axis=0)
    gy[:, :-1] = np.diff(new, axis=1)
    if im.ndim == 3:
        gz = np.zeros_like(new)
        gz[..., :-1] = np.diff(new, axis=2)
        tv_new = 2 * weight * np.sqrt(gx**2 + gy**2 + gz**2).sum()
    else:
        tv_new = 2 * weight * np.sqrt(gx**2 + gy**2).sum()
    dual_gap = (gap**2).sum() + tv_new - im_norm + (new**2).sum()
    return 0.5 / im_norm * dual_gap

def tv_denoise_fista(im, weight=50, eps=5.e-5, n_iter_max=200,
                                check_gap_frequency=3):

    """
    Perform total-variation denoising on 2-d and 3-d images

    Find the argmin `res` of
        1/2 * ||im - res||^2 + weight * TV(res),

    where TV is the isotropic l1 norm of the gradient.

    Parameters
    ----------
    im: ndarray of floats (2-d or 3-d)
        input data to be denoised. `im` can be of any numeric type,
        but it is cast into an ndarray of floats for the computation
        of the denoised image.

    weight: float, optional
        denoising weight. The greater ``weight``, the more denoising (at
        the expense of fidelity to ``input``)

    eps: float, optional
        precision required. The distance to the exact solution is computed
        by the dual gap of the optimization problem and rescaled by the l2
        norm of the image (for contrast invariance).

    n_iter_max: int, optional
        maximal number of iterations used for the optimization.

    Returns
    -------
    out: ndarray
        denoised array

    Notes
    -----
    The principle of total variation denoising is explained in
    http://en.wikipedia.org/wiki/Total_variation_denoising

    The principle of total variation denoising is to minimize the
    total variation of the image, which can be roughly described as
    the integral of the norm of the image gradient. Total variation
    denoising tends to produce "cartoon-like" images, that is,
    piecewise-constant images.

    This function implements the FISTA (Fast Iterative Shrinkage
    Thresholding Algorithm) algorithm of Beck et Teboulle, adapted to
    total variation denoising in "Fast gradient-based algorithms for
    constrained total variation image denoising and deblurring problems"
    (2009).
    """
    if not im.dtype.kind == 'f':
        im = im.astype(np.float)
    shape = [im.ndim, ] + list(im.shape)
    grad_im = np.zeros(shape)
    grad_aux = np.zeros(shape)
    t = 1.
    i = 0
    while i < n_iter_max:
        error = weight * div(grad_aux) - im
        grad_tmp = gradient(error)
        grad_tmp *= 1./ (8 * weight)
        grad_aux += grad_tmp
        grad_tmp = _projector_on_dual(grad_aux)
        t_new = 1. / 2 * (1 + np.sqrt(1 + 4 * t**2))
        t_factor = (t - 1) / t_new
        grad_aux = (1 + t_factor) * grad_tmp - t_factor * grad_im
        grad_im = grad_tmp
        t = t_new
        if (i % check_gap_frequency) == 0:
            gap = weight * div(grad_im)
            new = im - gap
            dgap = dual_gap(im, new, gap, weight)
            if dgap < eps:
                break
        i += 1
    return new


if __name__ == '__main__':
    from scipy.misc import face
    import matplotlib.pyplot as plt
    from time import time
    l = face().astype(np.float)
    # normalize image between 0 and 1
    l /= l.max()
    l += 0.1 * l.std() * np.random.randn(*l.shape)
    t0 = time()
    res = tv_denoise_fista(l, weight=0.05, eps=5.e-5)
    t1 = time()
    print(t1 - t0)
    plt.figure()
    plt.subplot(121)
    plt.imshow(l, cmap='gray')
    plt.subplot(122)
    plt.imshow(res, cmap='gray')
    plt.show()
