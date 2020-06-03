![genetIC](./genetic.svg)[![DOI](https://zenodo.org/badge/81329348.svg)](https://zenodo.org/badge/latestdoi/81329348) [![Build Status](https://travis-ci.com/pynbody/genetIC.svg?token=Kwgna3AKWpdHTHRrmaYX&branch=master)](https://travis-ci.com/pynbody/genetIC)

GenetIC
=======

GenetIC is a code to generate initial conditions for cosmological simulations, especially for zoom simulations of galaxies. It provides support for 'genetic modifications' as described by e.g. [Roth et al 2015](https://arxiv.org/abs/1504.07250), [Rey & Pontzen 2018](https://arxiv.org/abs/1706.04615).
* The main code is written in C++ and can be found (with build instructions and a test suite) in the subfolder `genetIC`. 
* A set of python tools for development and testing is provded in `tools`.
* If you wish to use the code for science, we recommend adopting a particular version. The current release is v1.0.2 and can be downloaded at https://github.com/pynbody/genetIC/releases/tag/v1.0.2
* Documentation can be found in the form of a PDF user manual, also at https://github.com/pynbody/genetIC/releases/tag/v1.0.2
* The technical underpinnings of the code are described in a release paper at https://github.com/pynbody/genetIC/releases/tag/v1.0.2 (this will be replaced with the arXiv link when it is available)

License and acknowledgement policy
----------------------------------

Legally, genetIC is released under the GNU public license v3; see LICENSE.md. However, for scientific integrity, if 
you use genetIC in any stage of preparing a publication, you must at a minimum cite the code release paper (see link 
above). We would recommend additionally citing the specific version of the code you used, using an appropriate Zenodo 
doi citation. These will be added to the repository when a first release is made.

If you receive assistance from any of the genetIC authors, we  also ask that you consider offering an acknowledgement 
or, in cases where the assistance is material to the development of a project, coauthorship to those who are helping you. 
