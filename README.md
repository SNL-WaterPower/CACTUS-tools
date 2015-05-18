CACTUS-tools
============
Tools for pre-processing CACTUS inputs and post-processing CACTUS results.

## pyCactus
Contributors: Phillip Chiu, pchiu@sandia.gov

A Python module containing post-processing classes and scripts.

### Installation
pyCactus can be installed as a module via setuptools.
Download the pyCactus directory, unzip, and run:

    python setup.py install
    
If setuptools is correctly configured, this will download and install the PyEVTK module, which is a dependency for writing VTK files.

## CreateGeom
Contributors: Matthew Barone, Jonathan C. Murray

A MATLAB/Octave library for creating CACTUS-formatted geometry files.