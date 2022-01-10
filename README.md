# Spectromicroscopy #
[Spectromicroscopy](http://spectromicroscopy.com) combines spectral data with microscopy,
where typical datasets consist of a stack of microscopic images
taken across an energy range. Due to the data complexity, manual analysis 
can be time consuming and inefficient, whereas multivariate analysis tools 
not only reduce the time needed but also can uncover hidden trends in the data.

# Mantis #
[MANTiS](http://spectromicroscopy.com) is Multivariate ANalysis Tool for Spectromicroscopy developed in Python by [2nd Look Consulting](http://2ndlookconsulting.com). It uses principal component analysis and cluster analysis to classify pixels according to spectral similarity.

## Download ##
Mantis package and binaries can be downloaded from 
[spectromicroscopy.com](http://spectromicroscopy.com).
Alternatively, you can install [Python](https://www.python.org/downloads/) and then run the command: `python3 -m pip install mantis-xray`

## Update ##
You can upgrade to the latest package release with the command: `pip3 install mantis-xray -U`.
It is recommended that you also upgrade the dependencies with: `pip3 install mantis-xray -U --upgrade-strategy "eager"`

## Run ##
Installation via pip provides the `mantis-xray` command (alternatively `python3 -m mantis_xray`) to start the Mantis GUI.

## User Guide ##
Mantis User Guide can be found on the project Wiki pages [Home](https://github.com/mlerotic/spectromicroscopy/wiki).

## References ##

Please use the following reference when quoting Mantis

Lerotic M, Mak R, Wirick S, Meirer F, Jacobsen C. MANTiS: a program for the analysis of X-ray spectromicroscopy data. J. Synchrotron Rad. 2014 Sep; 21(5); 1206–1212 [http://dx.doi.org/10.1107/S1600577514013964]
