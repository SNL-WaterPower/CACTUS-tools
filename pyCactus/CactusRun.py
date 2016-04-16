# CactusRun.py
"""A module for parsing CACTUS run (input and output) files."""

import os
import time as pytime

import numpy as np
import pandas as pd

from CactusGeom import CactusGeom
from CactusWakeElems import CactusWakeElems
from CactusField import CactusField
from CactusProbes import CactusProbes
from CactusInput import CactusInput

import warnings

from recursive_glob import recursive_glob


class CactusRun():
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
    elem_data : Pandas DataFrame
        Element data.
    wakeelems : CactusWakeElems class
        Wake element data class.
    field : CactusField class
        Field data class.
    probes : CactusProbes class.
        Probe class.
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
                 probe_fnames_pattern='probe_*.csv*'):
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
        elem_fname_pattern      = case_name + '_ElementData.csv'
        param_fname_pattern     = case_name + '_Param.csv'
        rev_fname_pattern       = case_name + '_RevData.csv'
        time_fname_pattern      = case_name + '_TimeData.csv'

        # search for wake element, field files, and probe files anywhere in
        # the run directory
        if load_wakeelem_output:
            self.wake_filenames = sorted(recursive_glob(run_directory,
                                                        wakeelem_fnames_pattern))
            if not self.wake_filenames:
                print 'Warning: Could not find any wake element data files in\
                       the work directory matching %s.' %\
                    (wakeelem_fnames_pattern)

        if load_field_output:
            self.field_filenames = sorted(recursive_glob(run_directory,
                                                         field_fnames_pattern))
            if not self.field_filenames:
                print 'Warning: Could not find any field data files in\
                       the work directory matching %s.' %\
                    (field_fnames_pattern)

        if load_probe_output:
            self.probe_filenames = sorted(recursive_glob(run_directory,
                                                         probe_fnames_pattern))
            if not self.probe_filenames:
                print 'Warning: Could not find any probe data files in\
                       the work directory matching %s.' %\
                    (probe_fnames_pattern)

        # Load the input, geometry, blade element, rev-averaged, parameter,
        # and time data. Only one of each file should be expected. The function
        # find_single_file is used to warn if multiple files (or none) are
        # found.

        # load the input namelist
        if self.input_fname:
            tic = pytime.time()
            self.input = CactusInput(self.input_fname)
            print 'Read input namelist in %2.2f s' % (pytime.time() - tic)
        else:
            warnings.warn("Input file not loaded.")

        # load geometry data
        if self.geom_fname:
            tic = pytime.time()
            # load the geometry data
            self.geom = CactusGeom(self.geom_fname)
            print 'Read geometry file in %2.2f s' % (pytime.time() - tic)
        else:
            warnings.warn("Geometry file not loaded.")

        # load parameter data
        self.param_fname = self.__find_single_file(
            run_directory, param_fname_pattern)
        if self.param_fname:
            tic = pytime.time()
            self.param_data  = self.__load_data(self.param_fname)
            print 'Read parameter data in %2.2f s' % (pytime.time() - tic)
        else:
            warnings.warn("Parameter data file not loaded.")

        # load revolution-averaged data
        self.rev_fname = self.__find_single_file(
            run_directory, rev_fname_pattern)
        if self.rev_fname:
            tic = pytime.time()
            self.rev_data  = self.__load_data(self.rev_fname)
            print 'Read revolution-averaged data in %2.2f s' %\
                (pytime.time() - tic)

        else:
            warnings.warn("Revolution-averaged data file not loaded.")

        # load blade element data
        self.elem_fname = self.__find_single_file(
            run_directory, elem_fname_pattern)
        if self.elem_fname:
            tic = pytime.time()
            self.elem_data  = self.__load_data(self.elem_fname)
            print 'Read blade element data in %2.2f s' % (pytime.time() - tic)
        else:
            warnings.warn("Blade element data file not loaded.")

        # time data
        self.time_fname = self.__find_single_file(
            run_directory, time_fname_pattern)
        if self.time_fname:
            tic = pytime.time()
            self.time_data  = self.__load_data(self.time_fname)
            print 'Read time data in %2.2f s' % (pytime.time() - tic)
        else:
            warnings.warn("Time data file not loaded.")

        # The following sections initialize the CactusWakeElems, CactusField,
        # and CactusProbes classes. Initializing these classes will search for
        # files in the run_directory and parse the first line of each. This may
        # be slow, depending on the number of files

        # wake element data
        self.wakeelems = CactusWakeElems(self.wake_filenames)

        # field data
        self.field = CactusField(self.field_filenames)

        # probe data
        self.probes = CactusProbes(self.probe_filenames)

        print ''
        print 'Success: Loaded case `%s` from path `%s`\n' % (case_name,
                                                              run_directory)

    #####################################
    #         Private Functions         #
    #####################################
    def __load_data(self, data_filename):
        # read a CSV file using pandas and returns a pandas dataframe
        reader = pd.read_csv(data_filename, iterator=True, chunksize=1000)
        df = pd.concat(reader, ignore_index=True)

        # strip whitespace from colnames
        df.rename(columns=lambda x: x.strip(), inplace=True)
        return df

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
    def blade_data_at_time_index(self, time_index):
        """Get blade data a specified time.

        Extract a subset dataframe of the "Element Data" dataframe by time
        index. Returns the time corresponding to the given time_index, and a
        list of dataframes containing the data, with one dataframe per blade.

        Parameters
        ----------
            time_index : int
                The time index.

        Returns
        -------
            time : float
                The non-dimensionalized time.
            dfs_blade : list
                A list of dataframes, each holding blade data.

        """
        # set column names
        time_col_name  = 'Normalized Time (-)'
        blade_col_name = 'Blade'
        elem_col_name  = 'Element'

        # get the data series
        df = self.elem_data

        # extract data subset
        df, time = self.df_subset_time_index(df, time_index, time_col_name)

        # get number of blades
        num_blades = self.geom.globalvars['NBlade'][0]
        # num_elems  = [blade['NElem'][0] for blade in self.geom.blades]

        # organize into a list of dataframes by blade
        dfs_blade = []
        for blade in range(1, num_blades + 1):
            # extract dataframe for particular blade
            df_blade = df[df[blade_col_name] == blade]

            # pop off the elemnts
            elements = df_blade.pop(elem_col_name)

            # set the elements as the df index
            df_blade.index = elements

            # append the df to the list of dfs
            dfs_blade.append(df_blade)

        return time, dfs_blade

    def blade_qty_avg(self, qty_name, timesteps_from_end, blade_num):
        """Compute atime-averaged a blade quantity distribution.

        Compute a time-averaged blade quantity distribution over a number of
        timesteps from the end. This may be used to, for example, compute the
        revolution-averaged blade quantities, such as circulation distribution,
        angle of attack, and local blade relative velocity.

        Parameters
        ----------
            qty_name : str
                The column name for the desired blade quantity
            timesteps_from_end : int
                Number of timesteps from the end
            blade_num : int
                Index of the blade number

        Returns
        -------
            qty_avg : numpy.array
                The array of averaged blade data.

        """
        # set up time indices over which we would like average
        time_indices = range(-1, -timesteps_from_end, -1)

        # initialize empty array to store circulation at final timestep
        qty_avg = np.zeros(self.geom.blades[blade_num]['NElem'],)

        # compute average
        for time_index in time_indices:
            # get time (float), and a tuple containing Pandas dataframes of
            # length num_blades.
            time, elem_data = self.blade_data_at_time_index(time_index)

            # get the circulation distribution at the final timestep
            qty = elem_data[0][qty_name]
            qty_avg = qty_avg + qty.values

        qty_avg = qty_avg / ((len(time_indices)))

        # return the averaged quantity
        return qty_avg

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
        df, time = self.df_subset_time_index(df, time_index, time_col_name)

        return time, df

    def df_subset_time_index(self, df, time_index, time_col_name):
        """Extract a subset dataframe of a given dataframe by time index.

        Parameters
        ----------
            df : pandas.DataFrame
                The dataframe to extract a subset of.
            time_index : int
                An integer of the time index.
            time_col_name : str
                The name of the column which contains time data.

        Returns
        -------
            df : pandas.DataFrame
                The dataframe subset and the time corresponding to the given
                time_index.
            time : float
                The time corresponding to the given time index.

        """
        # get unique times
        times = df.loc[:, time_col_name].unique()

        # time at which we wish to extract data
        time = times[time_index]

        # extract the subset of data corresponding to the desired time_index
        df = df[df[time_col_name] == times[time_index]]

        return df, time
