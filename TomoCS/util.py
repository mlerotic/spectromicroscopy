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
from scipy import ndimage


def generate_synthetic_data(l_x=128, seed=None, crop=True, n_pts=25):
    """
    Generate synthetic binary data looking like phase separation

    Parameters
    ----------

    l_x: int, default 128
        Linear size of the returned image

    seed: int, default 0
        seed with which to initialize the random number generator.

    crop: bool, default True
        If True, non-zero data are found only within a central circle
        of radius l_x / 2

    n_pts: int, default 25
        number of seeds used to generate the structures. The larger n_pts,
        the finer will be the structures.

    Returns
    -------

    res: ndarray of float32, of shape lxl
        Output binary image

    Examples
    --------
    >>> im = generate_synthetic_data(l_x=256, seed=2, n_pts=25)
    >>> # Finer structures
    >>> im = generate_synthetic_data(l_x=256, n_pts=100)
    """
    if seed is None:
        seed = 0
    # Fix the seed for reproducible results
    rs = np.random.RandomState(seed)
    x, y = np.ogrid[:l_x, :l_x]
    mask = np.zeros((l_x, l_x))
    points = l_x * rs.rand(2, n_pts)
    mask[(points[0]).astype(np.int), (points[1]).astype(np.int)] = 1
    mask = ndimage.gaussian_filter(mask, sigma=l_x / (4. * np.sqrt(n_pts)))
    # Limit the non-zero data to a central circle
    if crop:
        mask_outer = (x - l_x / 2) ** 2 + (y - l_x / 2) ** 2 < (l_x / 2) ** 2
        mask = np.logical_and(mask > mask.mean(), mask_outer)
    else:
        mask = mask > mask.mean()
    return mask.astype(np.float32)

def tv_l0_norm(im):
    """Compute the (isotropic) TV norm of an image"""
    grad_x1 = np.diff(im, axis=0)
    grad_x2 = np.diff(im, axis=1)
    return (grad_x1[:, :-1]**2 + grad_x2[:-1, :]**2 > 0).mean()


def compute_sparsity(im):
    l_x = len(im)
    X, Y = np.ogrid[:l_x, :l_x]
    mask = ((X - l_x/2)**2 + (Y - l_x/2)**2 <= (l_x/2)**2)
    grad1 = ndimage.morphological_gradient(im, footprint=np.ones((3, 3)))
    grad2 = ndimage.morphological_gradient(im, footprint=ndimage.generate_binary_structure(2, 1))
    return (grad1[mask] > 0).mean(), (grad2[mask] > 0).mean() 
