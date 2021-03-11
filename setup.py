import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mantis-xray", 
    version="3.0.03",
    author="Mirna Lerotic",
    author_email="mirna@2ndlookconsulting.com",
    description="MANTiS is a Multivariate ANalysis Tool for x-ray Spectromicroscopy",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://spectromicroscopy.com/",
    project_urls={
        "Code": "https://github.com/mlerotic/spectromicroscopy",
        "Documentation": "https://docs.spectromicroscopy.com",
    },
    install_requires=['PyQt5','numpy', 'scipy','matplotlib','h5py','Pillow','lxml','pyqtgraph'],
    extras_require={
        "netCDF":  "netcdf4-python",
        "SIRT":    "scikit-image"
    },
    entry_points={
        "console_scripts": "mantis-xray = mantis_xray.mantis_qt:main"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
    ],
    packages=setuptools.find_packages(),
    python_requires=">=3",
)
