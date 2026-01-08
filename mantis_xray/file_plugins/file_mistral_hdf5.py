# -*- coding: utf-8 -*-
#
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
#
#   Copyright (C) 2025 ALBA-MISTRAL Beamline
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

"""
MISTRAL ALBA Beamline HDF5 File Plugin

This plugin reads spectromicroscopy data from ALBA's MISTRAL beamline.
The data format contains aligned and normalized spectroscopic stacks with
metadata including energy points, pixel sizes, currents, and exposure times.

Data structure (fixed format):
    /SpecNormalized/
        spectroscopy_normalized_aligned  - Main data stack (n_energy, y, x) [float32]
        energy                           - Energy points in eV (n_energy,) [float64]
        x_pixel_size                     - Pixel size in x (1,) [float64] (micrometers)
        y_pixel_size                     - Pixel size in y (1,) [float64] (micrometers)
        Currents                         - Ring current values (n_energy,) [float64]
        ExpTimes                         - Exposure times (n_energy,) [float64]
        rotation_angle                   - Sample rotation angles (n_energy,) [float64]
        spectroscopy_normalized_trans    - Transformation matrix (n_energy, 3, 3) [float64]
"""

from __future__ import print_function
import sys, os, numpy, h5py
from collections import OrderedDict

title = "MISTRAL-ALBA"
extension = ["*.hdf5", "*.hdf", "*.h5"]
read_types = ["stack"]
write_types = []


def identify(filename):
    """
    Identify if this is a MISTRAL ALBA HDF5 file.

    Checks for the presence of 'SpecNormalized' group with the expected datasets.

    Args:
        filename: Path to the HDF5 file

    Returns:
        Boolean indicating if this plugin can read the file
    """
    try:
        f = h5py.File(filename, "r")
        # Check for SpecNormalized group with required datasets
        if "SpecNormalized" in f:
            grp = f["SpecNormalized"]
            has_data = "spectroscopy_normalized_aligned" in grp
            has_energy = "energy" in grp
            has_pixel_size = "x_pixel_size" in grp and "y_pixel_size" in grp
            f.close()
            return has_data and has_energy and has_pixel_size
        f.close()
        return False
    except:
        return False


def GetFileStructure(filename):
    """
    Get the internal structure of the MISTRAL HDF5 file.

    Args:
        filename: Path to the HDF5 file

    Returns:
        OrderedDict containing file structure information
    """
    try:
        f = h5py.File(filename, "r")

        # Check for SpecNormalized group
        if "SpecNormalized" not in f:
            f.close()
            return None

        grp = f["SpecNormalized"]
        D = OrderedDict()
        D["SpecNormalized"] = OrderedDict()
        D["SpecNormalized"]["type"] = "MISTRAL-ALBA STXM Stack"

        # Get data shape and dtype
        if "spectroscopy_normalized_aligned" in grp:
            data = grp["spectroscopy_normalized_aligned"]
            D["SpecNormalized"]["data_shape"] = data.shape
            D["SpecNormalized"]["data_dtype"] = str(data.dtype)

        # Get energy information
        if "energy" in grp:
            energy = grp["energy"]
            D["SpecNormalized"]["n_energies"] = len(energy)
            D["SpecNormalized"]["energy_range"] = (float(energy[0]), float(energy[-1]))

        # Check for optional datasets
        D["SpecNormalized"]["has_currents"] = "Currents" in grp
        D["SpecNormalized"]["has_exptimes"] = "ExpTimes" in grp
        D["SpecNormalized"]["has_rotation"] = "rotation_angle" in grp

        f.close()
        return D
    except:
        return None


def read(filename, stack_object, selection=None, *args, **kwargs):
    """
    Read MISTRAL ALBA HDF5 file into a MANTiS stack object.

    Args:
        filename: Path to the HDF5 file
        stack_object: MANTiS data stack object to populate
        selection: Tuple for selection (not used for MISTRAL)

    Returns:
        None (following MANTiS plugin convention)
    """
    # Open HDF5 file
    f = h5py.File(filename, "r")

    # Check for SpecNormalized group
    if "SpecNormalized" not in f:
        print("Error: SpecNormalized group not found in file")
        f.close()
        return

    grp = f["SpecNormalized"]

    # Read the main spectroscopy data
    # MISTRAL format: (n_energy, y, x)
    # MANTiS expects: (x, y, n_energy)
    data = numpy.array(grp["spectroscopy_normalized_aligned"]).astype("float32")
    stack_object.absdata = numpy.transpose(data, (2, 1, 0))

    # Read energy values
    stack_object.ev = numpy.array(grp["energy"])
    stack_object.n_ev = len(stack_object.ev)

    # Read pixel sizes (stored as arrays with shape (1,))
    x_pixel_size = float(grp["x_pixel_size"][0])
    y_pixel_size = float(grp["y_pixel_size"][0])

    # Get image dimensions
    stack_object.n_cols = stack_object.absdata.shape[0]
    stack_object.n_rows = stack_object.absdata.shape[1]

    # Create spatial coordinate arrays (in micrometers)
    stack_object.x_dist = numpy.linspace(
        0, stack_object.n_cols * x_pixel_size, num=stack_object.n_cols
    )
    stack_object.y_dist = numpy.linspace(
        0, stack_object.n_rows * y_pixel_size, num=stack_object.n_rows
    )

    # Read exposure times (dwell times)
    if "ExpTimes" in grp:
        stack_object.data_dwell = numpy.array(grp["ExpTimes"])
    else:
        stack_object.data_dwell = numpy.ones(stack_object.n_ev)

    # Optional: Store ring current data
    if "Currents" in grp:
        currents = numpy.array(grp["Currents"])
        # Data is already normalized, currents stored for reference

    # Optional: Store rotation angle
    if "rotation_angle" in grp:
        rotation_angle = numpy.array(grp["rotation_angle"])
        # Can be stored in metadata if needed

    f.close()

    # Fill the HDF5 structure for MANTiS
    stack_object.fill_h5_struct_from_stk()

    if verbose:
        print("Successfully loaded MISTRAL data:")
        print("  Data shape: %s" % (stack_object.absdata.shape,))
        print(
            "  Energy range: %.2f - %.2f eV" % (stack_object.ev[0], stack_object.ev[-1])
        )
        print(
            "  Spatial size: %d x %d pixels"
            % (stack_object.n_cols, stack_object.n_rows)
        )
        print("  Pixel size: %.4f x %.4f µm" % (x_pixel_size, y_pixel_size))

    return


# Module-level verbose flag (consistent with other plugins)
verbose = 0
