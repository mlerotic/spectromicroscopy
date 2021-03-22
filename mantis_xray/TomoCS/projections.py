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
from scipy import ndimage
from scipy import fftpack
from ._rank_order import rank_order

# --------------- Tomo projection operator  --------------------

def build_projection_operator(l_x, n_dir=None, angles=None, l_det=None, subpix=1,
                                offset=0, pixels_mask=None):
    """
    Compute the tomography design matrix.

    Parameters
    ----------

    l_x : int
        linear size of image array

    n_dir : int, default l_x
        number of angles at which projections are acquired. n_dir
        projection angles are regularly spaced between 0 and 180.

    l_det : int, default is l_x
        number of pixels in the detector. If l_det is not specified,
        we suppose that l_det = l_x.

    subpix : int, default 1
        number of linear subdivisions used to compute the projection of
        one image pixel onto a detector pixel. For example, if subpix=2,
        one image pixel is divided into 2x2 subpixels that are projected
        onto the detector, and the value of the projections is computed
        from these 4 projections.

    offset : int, default 0
        width of the strip of image pixels not covered by the detector. 
        offset > 0 means that the image is acquired in local tomography
        (aka ROI) mode, with the image larger than the detector. If the
        linear size of the array is l_x, the size of the detector is 
        l_x - 2 offset.

    pixels_mask : 1-d ndarray of size l_x**2
        mask of pixels to keep in the matrix (useful if one wishes
        to remove pixels inside or outside of a circle, for example)

    Returns
    -------
    p : sparse matrix of shape (n_dir l_x, l_x**2), in csr format
        Tomography design matrix. The csr (compressed sparse row)
        allows for efficient subsequent matrix multiplication. The
        dtype of the elements is float32, in order to save memory.

    Notes
    -----
    The returned matrix is sparse, but may nevertheless require a lot
    of memory for large l_x. For example, with l_x=512 and n_dir=512,
    the operator takes around 3 Gb of memory. The memory cost is of
    the order of l_x^2 x n_dir x 8 in bytes.

    For a given angle, the center of the pixels are rotated by the
    angle, and projected onto the detector. The value of the data pixel
    is added to the two pixels of the detector in between which the
    projection is located, with weights determined by a linear
    interpolation.

    Using subpix > 1 slows down the computation of the operator, because
    a histogram in 2-D has to be computed in order to group the projections 
    of subpixels corresponding to a single image pixel.
    (this should be accelerated by a Cython function... to be written)

    Examples
    --------
    >>> # Image with 256 pixels, 128 directions 
    >>> op = build_projection_operator(256, n_dir=128)

    >>> # Image with 128 pixels (to be reconstructed), 256 detector pixels
    >>> # subpix = 2 is used for a good precision of the projection of the
    >>> # coarse image pixels
    >>> op = build_projection_operator(128, n_dir=256, l_det=256, subpix=2)
    
    >>> # Image with 256 pixels, that is twice the size of the detector that
    >>> # has 128 pixels.
    >>> op = build_projection_operator(256, n_dir=256, l_det=128, offset=64)
    
    >>> # Image with 256 pixels, that is twice the size of the detector that
    >>> # has 256 pixels. We use subpixels for better precision.
    >>> op = build_projection_operator(256, n_dir=256, l_det=256, offset=64)
    
    >>> # Using a mask: projection operator only for pixels inside a
    >>> # central circle
    >>> l_x = 128
    >>> X, Y = np.ogrid[:l_x, :l_x]
    >>> mask = (X - l_x/2)**2 + (Y - l_x/2)**2 < (l_x/2)**2
    >>> op = build_projection_operator(l_x, pixels_mask=mask)
    >>> op.shape
    (16384, 12849)
    """
    if l_det is None:
        l_det = l_x
    X, Y = _generate_center_coordinates(subpix*l_x)
    X *= 1./subpix
    Y *= 1./subpix
    Xbig, Ybig = _generate_center_coordinates(l_det)
    Xbig *= (l_x - 2*offset) / float(l_det)
    orig = Xbig.min()
    labels = None
    if subpix > 1:
        # Block-group subpixels
        Xlab, Ylab = np.mgrid[:subpix * l_x, :subpix * l_x]
        labels = (l_x * (Xlab / subpix) + Ylab / subpix).ravel()
    if n_dir is None:
        n_dir = l_x
    if angles is None:
        angles = np.linspace(0, np.pi, n_dir, endpoint=False)
    else:
        n_dir = len(angles)
    weights, data_inds, detector_inds = [], [], []
    # Indices for data pixels. For each data, one data pixel
    # will contribute to the value of two detector pixels.
    for i, angle in enumerate(angles):
        # rotate data pixels centers
        Xrot = np.cos(angle) * X - np.sin(angle) * Y
        # compute linear interpolation weights
        inds, dat_inds, w = _weights_fast(Xrot, dx=(l_x - 2*offset)/float(l_det),
                                    orig=orig, labels=labels)
        # crop projections outside the detector
        mask = np.logical_and(inds >= 0, inds < l_det)
        weights.append(w[mask])
        detector_inds.append((inds[mask] + i * l_det).astype(np.int32))
        data_inds.append(dat_inds[mask])
    weights = np.concatenate(weights)
    weights /= subpix**2
    detector_inds = np.concatenate(detector_inds)
    data_inds = np.concatenate(data_inds)
    if pixels_mask is not None:
        if pixels_mask.ndim > 1:
            pixels_mask = pixels_mask.ravel()
        mask = pixels_mask[data_inds]
        data_inds = data_inds[mask]
        data_inds = rank_order(data_inds)[0]
        detector_inds = detector_inds[mask]
        weights = weights[mask]
    proj_operator = sparse.coo_matrix((weights, (detector_inds, data_inds)))
    return sparse.csr_matrix(proj_operator)

def _weights_fast(x, dx=1, orig=0, ravel=True, labels=None):
    """
    Compute linear interpolation weights for projection array `x`
    and regularly spaced detector pixels separated by `dx` and
    starting at `orig`.
    """
    if ravel:
        x = np.ravel(x)
    floor_x = np.floor((x - orig) / dx).astype(np.int32)
    alpha = ((x - orig - floor_x * dx) / dx).astype(np.float32)
    inds = np.hstack((floor_x, floor_x + 1))
    weights = np.hstack((1 - alpha, alpha))
    data_inds = np.arange(x.size, dtype=np.int32)
    data_inds = np.hstack((data_inds, data_inds))
    if labels is not None:
        data_inds = np.hstack((labels, labels))
        order = np.argsort(inds)
        inds, data_inds, weights = inds[order], data_inds[order], weights[order]
        steps = np.nonzero(np.diff(inds) > 0)[0] + 1
        steps = np.concatenate(([0], steps))
        inds_s, data_inds_s, weights_s = [], [], []
        for i in range(len(steps) - 1):
            d, w = data_inds[steps[i]:steps[i+1]], \
                      weights[steps[i]:steps[i+1]]
            count = np.bincount(d, weights=w)
            mask = count>0
            w = count[mask]
            weights_s.append(w)
            datind = np.arange(len(mask))[mask] 
            data_inds_s.append(datind)
            detind = inds[steps[i]]*np.ones(mask.sum()) 
            inds_s.append(detind)
            #stop
        inds = np.concatenate(inds_s)
        data_inds = np.concatenate(data_inds_s)
        weights = np.concatenate(weights_s)
    return inds, data_inds, weights

def _weights(x, dx=1, orig=0, ravel=True, labels=None):
    """
    Compute linear interpolation weights for projection array `x`
    and regularly spaced detector pixels separated by `dx` and
    starting at `orig`.
    """
    if ravel:
        x = np.ravel(x)
    floor_x = np.floor((x - orig) / dx).astype(np.int32)
    alpha = ((x - orig - floor_x * dx) / dx).astype(np.float32)
    inds = np.hstack((floor_x, floor_x + 1))
    weights = np.hstack((1 - alpha, alpha))
    data_inds = np.arange(x.size, dtype=np.int32)
    data_inds = np.hstack((data_inds, data_inds))
    if labels is not None:
        data_inds = np.hstack((labels, labels))
        w = np.histogram2d(data_inds, inds,
             bins=(np.arange(data_inds.max()+1.5), np.arange(inds.max()+1.5)),
             weights=weights)[0]
        data_inds, inds = np.argwhere(w>0).T
        weights = w[w>0]
    return inds, data_inds, weights

def _weights_nn(x, dx=1, orig=0, ravel=True):
    """
    Nearest-neighbour interpolation
    """
    if ravel:
        x = np.ravel(x)
    floor_x = np.floor(x - orig)
    return floor_x.astype(np.float32)



def _generate_center_coordinates(l_x):
    """
    Compute the coordinates of pixels centers for an image of
    linear size l_x
    """
    l_x = float(l_x)
    X, Y = np.mgrid[:l_x, :l_x]
    center = l_x / 2.
    X += 0.5 - center
    Y += 0.5 - center
    return X, Y



# ----------------- Direct projection method -------------------------
# (without computing explicitely the design matrix)


def back_projection(projections):
    """
    Back-projection (without filtering)

    Parameters
    ----------
    projections: ndarray of floats, of shape n_dir x l_x
        Each line of projections is the projection of a data image
        acquired at a different angle. The projections angles are
        supposed to be regularly spaced between 0 and 180.

    Returns
    -------
    recons: ndarray of shape l_x x l_x
        Reconstructed array

    Notes
    -----
    A linear interpolation is used when rotating the back-projection.
    This function uses ``scipy.ndimage.rotate`` for the rotation.
    """
    n_dir, l_x = projections.shape
    recons = np.zeros((l_x, l_x), dtype=np.float)
    angles = np.linspace(0, 180, n_dir, endpoint=False)
    for angle, line in zip(angles, projections):
        # BP: repeat the detector line along the direction of the beam
        tmp = np.tile(line[:, np.newaxis], (1, l_x))
        # Rotate the back-projection of the detector line, and add
        # it to the reconstructed image
        recons += ndimage.rotate(tmp, -angle, order=1, \
                            reshape=False)
    return recons


def projection(im, n_dir=None, interpolation='nearest'):
    """
    Tomography projection of an image along n_dir directions.

    Parameters
    ----------
    im : ndarray of square shape l_x x l_x
        Image to be projected

    n_dir : int
        Number of projection angles. Projection angles are regularly spaced
        between 0 and 180.

    interpolation : str, {'interpolation', 'nearest'}
        Interpolation method used during the projection. Default is
        'nearest'.

    Returns
    -------
    projections: ndarray of shape n_dir x l_x
        Array of projections.

    Notes
    -----
    The centers of the data pixels are projected onto the detector, then
    the contribution of a data pixel to detector pixels is computed
    by nearest neighbor or linear interpolation. The function 
    ``np.bincount`` is used to compute the projection, with weights
    corresponding to the values of data pixels, multiplied by interpolation
    weights in the case of linear interpolation.
    """
    l_x = len(im)
    if n_dir is None:
        n_dir = l_x
    im = im.ravel()
    projections = np.empty((n_dir, l_x))
    X, Y = _generate_center_coordinates(l_x)
    angles = np.linspace(0, np.pi, n_dir, endpoint=False)
    for i, angle in enumerate(angles):
        Xrot = np.cos(angle) * X - np.sin(angle) * Y 
        if interpolation == 'nearest':
            inds = _weights_nn(Xrot, dx=1, orig=X.min())
            mask = inds>= 0
            w = im[mask]
        elif interpolation == 'linear':
            inds, _, w = _weights(Xrot, dx=1, orig=X.min())
            w[:l_x**2] *= im
            w[l_x**2:] *= im
            mask = inds >= 0
            w = w[mask]
        projections[i] = np.bincount(inds[mask].astype(np.int), \
                                            weights=w)[:l_x]
    return projections

# -----------------Filtered back-projection----------------------

def filter_projections(proj_set, reg=False):
    """
    Ramp filter used in the filtered back projection.
    We use zero padding.

    Parameters
    ----------
    proj_set: 2-d ndarray
        each line is one projection (1 line of the detector) to be filtered

    Returns
    -------

    res: 2-d ndarray
        filtered projections

    Notes
    -----

    We use zero padding. However, we do not use any filtering (hanning, etc.)
    in the FFT yet.
    """
    nb_angles, l_x = proj_set.shape
    # Assum l is even for now
    ramp = 1./l_x * np.hstack((np.arange(l_x), np.arange(l_x, 0, -1)))
    return fftpack.ifft(ramp * fftpack.fft(proj_set, 2*l_x, axis=1), axis=1)[:,:l_x]

