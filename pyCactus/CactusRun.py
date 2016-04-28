# CactusRun.py
"""A module for parsing CACTUS run (input and output) files."""

import os
import time as pytime
import math

from CactusGeom import CactusGeom
from CactusWakeElems import CactusWakeElems
from CactusField import CactusField
from CactusProbes import CactusProbes
from CactusInput import CactusInput
from CactusBladeElem import CactusBladeElem

import warnings

from common_utils import load_data, df_subset_time_index, recursive_glob


class CactusRun(object):
    """Class for interrogating a CACTUS input deck.

    Attributes
    ----------
    input : CactusInput class
        Input file class.
    geom : CactusGeom class
        Geometry data class.
    param_data : Pandas DataFrame
        Parameter data.
    rev_data: Pandas DataFrame
        Revolution-averaged data.
    time_data : Pandas DataFrame
        Time data.
    bladeelem_data : CactusBladeElem class
        Blade element data.
    wakeelems : CactusWakeElems class
        Wake element data class.
    field : CactusField class
        Field data class.
    probes : CactusProbes class.
        Probe class.
    nti : int
        Number of timesteps per iteration.
    tsr : float
        Tip speed ratio (non-dimensional).
    dt : float
        Normalized timestep length (non-dimensional).
    input_fname : str
        Input data filename.
    geom_fname : str
        Geometry data filename.
    param_fname : str
        Parameter data filename.
    rev_fname : str
        Revolution-averaged data filename.
    elem_fname : str
        Blade element data filename.
    time_fname : str
        Time data filename.
    wake_filenames : list
        List of filenames containing wake element data.
    field_filenames : list
        List of filenames containing field data.
    probe_filenames : list
        List of filenames containing probe data.
    """

    def __init__(self, run_directory, case_name,
                 input_fname='',
                 geom_fname='',
                 load_field_output=True,
                 load_wakeelem_output=True,
                 load_probe_output=True,
                 wakeelem_fnames_pattern='*WakeElemData_*.csv',
                 field_fnames_pattern='*FieldData_*.csv',
                 probe_fnames_pattern='probe_*.csv*',
                 quiet=False):
        """Initialize the class, reading some data to memory.

        This method relies on recursive searches within the specified run
        directory to find the appropriate CACTUS output files. Therefore, each
        run directory should only contain one set of output files (or else the
        behavior cannot be guaranteed).

        Parameters
        ----------
        run_directory : str
            Path to the directory containing the CACTUS run.
        case_name : str
            'case name' which precedes all input and output files.
        input_fname : Optional[str]
            Input filename (default `./[case_name].in`).
        geom_fname : Optional[str]
            Geometry filename (default `./[case_name].geom`)
        load_field_output : bool
            True (default) to load field data, False otherwise.
        load_wakeelem_output : bool
            True (default) to load wake element data, False otherwise.
        load_probe_output : bool
            True (default) to load probe data, False otherwise.
        wakeelem_fnames_pattern : Optional[str]
            Glob pattern for wake element data filenames (default is
            `*WakeElemData_*.csv`)
        field_fnames_pattern : Optional[str]
            Glob pattern for field data filenames (default is
            `*FieldData_*.csv`)
        probe_fnames_pattern : Optional[str]
            Glob pattern for probe filenames (default is `probe_*.csv`)
        quiet : Optional[bool]
            Set True to hide print statements (default is False).
        """
        # if an input file is specified, use that
        if input_fname:
            self.input_fname = os.path.abspath(os.path.join(run_directory,
                                                            input_fname))
        else:
            # otherwise, look for one using [case_name].in as a glob pattern
            self.input_fname = self.__find_single_file(run_directory,
                                                       case_name + '.in')

        # if a geom file is specified, use that
        if geom_fname:
            self.geom_fname = os.path.abspath(os.path.join(run_directory,
                                                           geom_fname))
        else:
            # otherwise, look for one using [case_name].geom as a glob pattern
            self.geom_fname = self.__find_single_file(run_directory,
                                                      case_name + '.geom')

        # assemble filename patterns
        bladeelem_fname_pattern = case_name + '_ElementData.csv'
        param_fname_pattern     = case_name + '_Param.csv'
        rev_fname_pattern       = case_name + '_RevData.csv'
        time_fname_pattern      = case_name + '_TimeData.csv'

        # Load the input, geometry, blade element, rev-averaged, parameter,
        # and time data. Only one of each file should be expected. The function
        # find_single_file is used to warn if multiple files (or none) are
        # found.

        # load the input namelist
        if self.input_fname:
            tic = pytime.time()
            self.input = CactusInput(self.input_fname)
            if not quiet:
                print 'Read input namelist in %2.2f s' % (pytime.time() - tic)
        else:
            warnings.warn("Input file not loaded.")

        # load geometry data
        if self.geom_fname:
            tic = pytime.time()
            # load the geometry data
            self.geom = CactusGeom(self.geom_fname)
            if not quiet:
                print 'Read geometry file in %2.2f s' % (pytime.time() - tic)
        else:
            warnings.warn("Geometry file not loaded.")

        # load parameter data
        self.param_fname = self.__find_single_file(
            run_directory, param_fname_pattern)
        if self.param_fname:
            tic = pytime.time()
            self.param_data  = load_data(self.param_fname)
            if not quiet:
                print 'Read parameter data in %2.2f s' % (pytime.time() - tic)
        else:
            warnings.warn("Parameter data file not loaded.")

        # load revolution-averaged data
        self.rev_fname = self.__find_single_file(
            run_directory, rev_fname_pattern)
        if self.rev_fname:
            tic = pytime.time()
            self.rev_data  = load_data(self.rev_fname)
            if not quiet:
                print 'Read revolution-averaged data in %2.2f s' %\
                    (pytime.time() - tic)

        else:
            warnings.warn("Revolution-averaged data file not loaded.")

        # load blade element data
        self.bladeelem_fname = self.__find_single_file(
            run_directory, bladeelem_fname_pattern)
        if self.bladeelem_fname:
            tic = pytime.time()
            self.bladeelem_data  = CactusBladeElem(self.bladeelem_fname)
            if not quiet:
                print 'Read blade element data in %2.2f s' % (pytime.time() -
                                                              tic)
        else:
            warnings.warn("Blade element data file not loaded.")

        # time data
        self.time_fname = self.__find_single_file(
            run_directory, time_fname_pattern)
        if self.time_fname:
            tic = pytime.time()
            self.time_data  = load_data(self.time_fname)
            if not quiet:
                print 'Read time data in %2.2f s' % (pytime.time() - tic)
        else:
            warnings.warn("Time data file not loaded.")

        # The following sections initialize the CactusWakeElems, CactusField,
        # and CactusProbes classes. Initializing these classes will search for
        # files in the run_directory and parse the first line of each. This may
        # be slow, depending on the number of files

        # search for wake element, field files, and probe files anywhere in
        # the run directory
        if load_wakeelem_output:
            self.wake_filenames = sorted(recursive_glob(run_directory,
                                                        wakeelem_fnames_pattern))
            if self.wake_filenames:
                self.wakeelems = CactusWakeElems(self.wake_filenames)
            else:
                if not quiet:
                    print 'Warning: Could not find any wake element data files \
in the work directory matching %s.' % \
                          (wakeelem_fnames_pattern)

        if load_field_output:
            self.field_filenames = sorted(recursive_glob(run_directory,
                                                         field_fnames_pattern))
            if self.field_filenames:
                self.field = CactusField(self.field_filenames)
            else:
                if not quiet:
                    print 'Warning: Could not find any field data files in \
the work directory matching %s.' % \
                          (field_fnames_pattern)

        if load_probe_output:
            self.probe_filenames = sorted(recursive_glob(run_directory,
                                                         probe_fnames_pattern))
            if self.probe_filenames:
                self.probes = CactusProbes(self.probe_filenames)
            else:
                if not quiet:
                    print 'Warning: Could not find any probe data files in \
the work directory matching %s.' % \
                          (probe_fnames_pattern)

        if not quiet:
            print 'Loaded case `%s` from path `%s`\n' % (case_name,
                                                         run_directory)

    #####################################
    #         Private Functions         #
    #####################################
    def __find_single_file(self, directory, pattern):
        # Look for a glob pattern in a specified directory and return the
        # first file, throwing a warning if multiple files are found.
        # Return None if no files are found.
        results = recursive_glob(directory, pattern)

        # warn if we found too many files or none at all
        if results:
            if len(results) > 1:
                warnings.warn("Found multiple files matching %s in %s, \
                              using %s" %
                              (pattern, directory, results[0]), RuntimeWarning)
        else:
            warnings.warn("Warning: Could not find file %s in %s" %
                          (pattern, directory), RuntimeWarning)

        if results:
            return results[0]
        else:
            return None

    ####################################
    #         Public Functions         #
    ####################################
    def rotor_data_at_time_index(self, time_index):
        """Extract a single time instance from the time dataframe.

        Returns the time corresponding to the given time_index, and a dataframe
        containing the appropriate subset of data.

        Parameters
        ----------
            time_index : int
                An integer of the time_index

        Returns
        -------
            time : float
                The time corresponding to the given time index
            df : pandas.DataFrame
                The dataframe containing the time data at that particular
                instance.

        """
        # get the data series
        df = self.time_data

        # set column names
        time_col_name  = 'Normalized Time (-)'

        # extract data subset
        df, time = df_subset_time_index(df, time_index, time_col_name)

        return time, df

    def rev_to_time(self, rev):
        """Compute the normalized time from a revolution."""
        timestep_number = rev * self.nti
        return timestep_number * self.dt

    def time_to_rev(self, time):
        """Compute the fractional revolution from a normalized time."""
        return time * self.nti / self.dt

    def rev_to_timestep(self, rev):
        """Compute the normalized time from a fractional revolution."""
        tol = 1e-5
        timestep = rev * self.nti
        if (timestep % 1) > tol:
            warnings.warn("The computed timestep is not within %e of an \
integer, but it was rounded to an integer value." % tol)
        return int(timestep)

    @property
    def nti(self):
        """Simulation nti parameter (timesteps per iteration).

        Extracted from the input namelist file."""
        return self.input.namelist['configinputs']['nti']

    @property
    def tsr(self):
        """Simulation tsr parameter (tip speed ratio)."""
        return self.param_data['TSR (-)'].values[0]

    @property
    def dt(self):
        """The simulation timestep (non-dimensional)."""
        return 2 * math.pi / (self.tsr * self.nti)

    @property
    def period(self):
        """The simulation time for one revolution period (non-dimensional)."""
        return self.dt * self.nti
