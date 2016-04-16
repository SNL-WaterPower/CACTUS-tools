import os
import f90nml
import warnings

class CactusInput(object):
    """Class for reading and writing CACTUS input namelist files.

    Attributes
    ----------
    namelist : f90nml.namelist.Namelist
        A f90nml Namelist object with the namelist contents.
    """

    def __init__(self, filename):
        """Reads in the CACTUS input namelist."""
        if os.path.exists(filename):
            # read in the input file namelist
            self.namelist = f90nml.read(filename)
        else:
            warnings.warn("Input file %s does not exist." % filename,
                          RuntimeWarning)
