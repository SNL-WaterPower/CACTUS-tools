from setuptools import setup

setup(
    name = 'pyCactus',
    version = '0.1',
    description = 'Post-processing tools for CACTUS output files.',
    author = 'Phillip Chiu',
    author_email = 'pchiu@sandia.gov',
    url = 'https://github.com/SNL-WaterPower/CACTUS-tools',
    py_modules = ['pyCactus', 'pyCactusGeom', 'pyCactusWake'],
    install_requires = ['pyevtk'],
    dependency_links=['https://bitbucket.org/pauloh/pyevtk']
)
