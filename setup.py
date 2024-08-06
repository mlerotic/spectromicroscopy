import setuptools
import mantis_xray

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mantis-xray", 
    version=mantis_xray.__version__,
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
    install_requires=['PyQt5>=5.15.9','numpy', 'scipy>=1.14.0', 'matplotlib>=3.6.0', 'h5py', 'Pillow', 'lxml', 'pyqtgraph>=0.13.7', "scikit-image>=0.19.1"],
    extras_require={
        "netCDF":  "netcdf4-python"},
    entry_points={
        "gui_scripts": "mantis-xray = mantis_xray.mantis_qt:main"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
    ],
    packages=setuptools.find_packages(),
    include_package_data=True,
    python_requires=">=3",
)
