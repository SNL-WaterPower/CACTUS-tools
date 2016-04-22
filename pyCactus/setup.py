"""CACTUS-tools setup file."""

from setuptools import setup

setup(
    name='pyCactus',
    version='0.1',
    description='Post-processing tools for CACTUS output files.',
    author='Phillip Chiu',
    author_email='pchiu@sandia.gov',
    url='https://github.com/SNL-WaterPower/CACTUS-tools',
    py_modules=['CactusRun',
                'CactusGeom',
                'CactusWakeElems',
                'CactusField',
                'CactusProbes',
                'CactusBladeElem'
                'common_utils'],
    install_requires=['PyEVTK', 'f90nml'],
    dependency_links=['https://bitbucket.org/pauloh/pyevtk',
                      'https://pypi.python.org/pypi/f90nml'],
    scripts=['./scripts/pyCactusCloneCase.py',
             './scripts/pyCactusCsvToVtk.py',
             './scripts/pyCactusPvdSubset.py',
             './scripts/pyCactusSplitFileByTime.py',
             './scripts/pyCactusWPDataToVtk.py']
)
