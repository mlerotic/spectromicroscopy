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
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
    ],
    packages=setuptools.find_packages(),
    python_requires=">=3.0",
)
