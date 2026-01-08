# MANTiS Development Guide for AI Coding Agents

This guide provides essential information for AI coding agents working on the MANTiS (Multivariate ANalysis Tool for x-ray Spectromicroscopy) codebase.

## Project Overview

- **Project**: MANTiS - Scientific GUI application for X-ray spectromicroscopy data analysis
- **Language**: Python 3 (>=3.0)
- **Version**: 3.2.10
- **GUI Framework**: PyQt5
- **License**: GNU GPL v3
- **Main Package**: `mantis_xray`

## Repository Structure

```
mantis_xray/                    # Main package directory
├── __init__.py                 # Package version
├── __main__.py                 # Entry point
├── mantis_qt.py                # Qt GUI application (main)
├── mantis.py                   # Core batch mode logic
├── data_stack.py               # Data stack management
├── data_struct.py              # Data structures
├── analyze.py                  # Analysis functions (PCA, clustering)
├── nnma.py                     # Non-negative matrix analysis
├── tomo_reconstruction.py      # Tomography reconstruction
├── henke.py                    # Henke database handling
├── helpers.py                  # Helper functions
├── Mrc.py                      # MRC file handling
├── logos.py                    # Logo data
├── *.ui                        # Qt Designer UI files
├── stylesheet_global.qss       # Qt stylesheet
├── file_plugins/               # File format plugins
│   ├── __init__.py             # Plugin loader system
│   ├── file_*.py               # Individual file format handlers
└── TomoCS/                     # Tomography compressed sensing
    └── *.py                    # Tomography algorithms
```

## Build, Install, and Run Commands

### Installation
```bash
# Install in development mode
pip install -e .

# Install with all dependencies
pip install .

# Install with optional NetCDF support
pip install .[netCDF]

# Upgrade to latest version
pip install mantis-xray -U --upgrade-strategy "eager"
```

### Running the Application
```bash
# Run GUI (after installation)
mantis-xray

# Alternative method
python3 -m mantis_xray
```

### Building Installers
```bash
# Windows installer
pyinstaller installer_files/mantis-win.spec

# macOS installer
pyinstaller installer_files/mantis-mac.spec

# Linux installer
pyinstaller installer_files/mantis-linux.spec
```

### Documentation
```bash
# Build documentation (requires mkdocs)
mkdocs build

# Serve documentation locally
mkdocs serve
```

## Testing

**NOTE**: This project currently has NO automated test suite.

When adding new features:
- Manually test GUI functionality
- Test with real spectromicroscopy data files
- Verify cross-platform compatibility (Windows, macOS, Linux)
- Test file plugins with various data formats

## Code Style Guidelines

### File Headers
Every Python file should include the GPL v3 header:
```python
# -*- coding: utf-8 -*-
#
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
#
#   Copyright (C) 2011-2013 Mirna Lerotic, 2nd Look
#   http://2ndlookconsulting.com
#   License: GNU GPL v3
#
#   [Standard GPL v3 disclaimer...]
```

### Imports
- **Order**: Standard library → Third-party → Local imports
- **Standard library**: `os, sys, re, urllib, platform`, etc.
- **Third-party**: `numpy, scipy, PyQt5, matplotlib, h5py, PIL, lxml, pyqtgraph, scikit-image`
- **Local**: Use relative imports: `from . import module_name` or `from .. import module_name`

Example:
```python
import os
import sys
from urllib.error import URLError

import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets
import matplotlib.pyplot as plt

from . import data_stack
from . import helpers
```

### Naming Conventions
- **Classes**: PascalCase (e.g., `DataStack`, `FigureCanvas`)
- **Functions**: snake_case (e.g., `save_keyeng`, `check_for_updates`, `resource_path`)
- **Variables**: snake_case (e.g., `file_path`, `current_version`, `data_array`)
- **Constants**: UPPER_CASE or PascalCase (e.g., `PlotH`, `PlotW`, `MAX_ATTEMPT`)
- **Private methods**: Prefix with underscore (e.g., `_internal_method`)

### Formatting
- **Indentation**: 4 spaces (no tabs)
- **Line endings**: CRLF (`\r\n`) for Windows compatibility
- **Max line length**: No strict limit, but keep readable (~100-120 characters preferred)
- **Docstrings**: Triple quotes for module/class/function documentation
- **Comments**: Use `#` for inline comments

### Type Hints
- Type hints are NOT currently used in this codebase
- Focus on clear variable names and docstrings instead

### Error Handling
```python
# Use try-except blocks for I/O and external calls
try:
    with urllib.request.urlopen(url, timeout=timeout) as response:
        data = response.read()
except URLError:
    print("Connection failed")
except Exception as e:
    print(f"Unexpected error: {e}")
```

### Docstrings
```python
def function_name(param1, param2):
    """ Brief description of what the function does """
    # Implementation
    return result
```

## File Plugin Architecture

When creating a new file format plugin in `mantis_xray/file_plugins/`:

### Required Attributes
```python
title = 'Plugin Name'                    # Short descriptive name
extension = ['*.ext', '*.ext2']          # List of file extensions
read_types = ['image', 'stack']          # Data types plugin can read
write_types = ['image', 'stack']         # Data types plugin can write (optional)
```

### Required Functions
```python
def identify(filename):
    """Returns boolean indicating if plugin can read this file"""
    # Check file signature, extension, or contents
    return True/False

def GetFileStructure(filename):
    """Returns structure describing internal file organization"""
    # Return dict with data sets available
    return structure_dict

def read(filename, self, selection, *args, **kwargs):
    """Loads data from filename into stack_object (self)"""
    # Load data into self.data_struct
    # Set metadata in self.data_struct
    return success_boolean
```

### Data Types
Common data types: `'spectrum'`, `'image'`, `'stack'`, `'results'`

## Dependencies

### Core Dependencies
```python
PyQt5 >= 5.15.9
numpy
scipy >= 1.11.4
matplotlib >= 3.6.0
h5py
Pillow
lxml
pyqtgraph >= 0.13.7
scikit-image >= 0.19.1
xdrlib3  # For Python 3.13+
```

### Optional Dependencies
```python
netcdf4-python  # For NetCDF file support
```

## Common Patterns

### NumPy Arrays
- Use NumPy arrays for all image and spectral data
- Prefer vectorized operations over loops

### PyQt5 GUI
- UI files (*.ui) designed in Qt Designer
- Load UI files at runtime with `uic.loadUi()`
- Connect signals/slots for event handling

### Resource Paths (for PyInstaller compatibility)
```python
from .helpers import resource_path

# Get path to bundled resource
path = resource_path('relative/path/to/resource')
```

## Version Management

- Version is stored in `mantis_xray/__init__.py`
- Update version string when making releases: `__version__ = 'X.Y.Z'`

## Best Practices

1. **Compatibility**: Maintain cross-platform compatibility (Windows, macOS, Linux)
2. **PyInstaller**: Ensure code works both in development and when frozen with PyInstaller
3. **Scientific Data**: Handle large arrays efficiently with NumPy
4. **File I/O**: Use context managers (`with` statements) for file operations
5. **Error Messages**: Provide clear, user-friendly error messages in GUI dialogs
6. **Plugin System**: Use the plugin architecture for new file format support
7. **Documentation**: Update `docs/index.md` when adding user-facing features

## Common Tasks

### Adding a New File Format
1. Create `mantis_xray/file_plugins/file_newformat.py`
2. Implement required plugin interface (title, extension, identify, GetFileStructure, read)
3. Plugin will be auto-discovered on next import

### Modifying the GUI
1. Edit *.ui files in Qt Designer
2. UI files are loaded at runtime - no compilation needed
3. Update corresponding Python code in `mantis_qt.py`

### Adding Analysis Functions
1. Add new methods to `mantis_xray/analyze.py`
2. Integrate with GUI in `mantis_qt.py`
3. Update user documentation in `docs/index.md`

## Resources

- **Website**: https://spectromicroscopy.com
- **Documentation**: https://docs.spectromicroscopy.com
- **Repository**: https://github.com/mlerotic/spectromicroscopy
- **PyPI**: https://pypi.org/project/mantis-xray/
- **Citation**: Lerotic M, et al. J. Synchrotron Rad. 2014; 21(5); 1206–1212

---
*Generated for AI coding agents working on MANTiS*
