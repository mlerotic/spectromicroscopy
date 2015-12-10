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

import os
import numpy as np
import scipy as sp


import warnings
warnings.simplefilter('ignore', DeprecationWarning)

from TomoCS.forward_backward_tv import fista_tv,gfb_tv
from TomoCS.projections import build_projection_operator
from TomoCS.util import generate_synthetic_data


#----------------------------------------------------------------------
class Ctomo:
    def __init__(self, stkdata):

        self.stack = stkdata
        
        
        
#----------------------------------------------------------------------   
# Calculate tomo - CS reconstruction for 1 dataset 
    def calc_tomo1(self,tomodata, angles):
        
        pass