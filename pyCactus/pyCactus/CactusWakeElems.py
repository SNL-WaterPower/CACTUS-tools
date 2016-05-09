
import os
import time as pytime

import numpy as np
import pandas as pd

from warnings import *

from common_utils import load_data, get_file_time


class CactusWakeElems(object):
    """Class for manipulating wake element data from CACTUS.

    Attributes
    ----------
    filenames : list
        Filenames (str) of wake element data files.
    num_times : int
        Number of timesteps of wake element data.
    times : list
        List of times (float) of wake element data.
    fdict : dict
        A dictionary of `{time : fname}`.
    """

    def __init__(self, filenames,
                 read_headers=True):
        """Initializes the instance."""
        
        self.filenames = filenames
        self.num_times = len(filenames)
        self.times = []
        self.fdict = {}

        # read the file headers for times
        if read_headers:
            tic = pytime.time()
            self.read_file_headers()
            print 'Read %d wake element data headers in %2.2f s' %\
                        (len(self.times),
                         pytime.time() - tic)

    def read_file_headers(self):
        """Gets the times from the headers of each data file in the instance
        `filenames` attribute amd stores to instance variables `self.times` and
        `self.fdict`

        Returns
        -------
        times : list
            List of times (float) of wake element data.
        fdict : dict
            A dictionary of `{time : fname}`.
        """
        # the name of the column containing time info
        time_col_name = 'Normalized Time (-)'
    
        self.fdict = {}
        self.times = []

        for fname in self.filenames:
            time = get_file_time(fname, time_col_name)
            self.times.append(time)
            self.fdict[time] = fname

        # sort the times
        self.times = np.sort(np.array(self.times))

        return self.times, self.fdict

    def get_df_inst(self, time=None, fname=None):
        """Gets the data from a specified time or filename.

        Either the time or the filename must be specified.

        Parameters
        ----------
        time : Optional[float]
            The time at which to extract the dataframe.
        fname : Optional[str]
            The filename to read (defaults to self.fdict[time]).

        Returns
        -------
        df_inst : pandas.DataFrame
            DataFrame of the time.
        """

        if (time is None) and (fname is None):
             print 'Error: must specify either the time or filename of the desired data.'

        if time is not None:
            # if the time is specified, get the filename
            fname = self.fdict[time]
        else:
            # otherwise, use the filename given
            pass

        # read the CSV data file
        df_inst = load_data(fname)

        return df_inst


    def wakedata_from_df(self, df):
        """Takes a dataframe containing wake node data at a single timestep and
        returns a dictionary of the data as np.arrays.

        NumPy arrays are keyed by a descriptive variable name.

        Parameters
        ----------
        df : pandas.DataFrame

        Returns
        -------
        data_arrays : dict
            Dictionary of np.arrays containing the data.
        has_node_ids : bool
            True if wake element data has 'Node ID' column, False if not.
        """

        # column names
        id_col_name   = 'Node ID'
        elem_col_name = 'Origin Node'
        x_col_name    = 'X/R (-)'
        y_col_name    = 'Y/R (-)'
        z_col_name    = 'Z/R (-)'
        u_col_name    = 'U/Uinf (-)'
        v_col_name    = 'V/Uinf (-)'
        w_col_name    = 'W/Uinf (-)'

        # check if data has node IDs
        if id_col_name in df:
            has_node_ids = True
        else:
            has_node_ids = False
            print 'Warning: Wake element data does not include node IDs.'
        
        # extract columns
        elems = df.loc[:,elem_col_name].values
        x     = df.loc[:,x_col_name].values
        y     = df.loc[:,y_col_name].values
        z     = df.loc[:,z_col_name].values
        u     = df.loc[:,u_col_name].values
        v     = df.loc[:,v_col_name].values
        w     = df.loc[:,w_col_name].values

        # store data as a list of np.arrays
        data_arrays = {'elems' : elems,
                       'x' : x,
                       'y' : y,
                       'z' : z,
                       'u' : u,
                       'v' : v,
                       'w' : w}

        # if the data has node ids, append this data to the data list
        if has_node_ids:
            node_ids = df.loc[:,id_col_name].values
            data_arrays['node_ids'] =  node_ids

        return data_arrays, has_node_ids


    def write_vtk_series(self, path, name,
                         id_flag=False,
                         num_blade_elems=[],
                         print_status=True):
        """Writes the wake element data to a time series of VTK files

        Data is written as VTK unstructured data (.vtu). ID data is a scalar
        (integer). Velocity data is a vector.

        Also writes a Paraview collection file (.pvd) which contains the
        normalized times of each timestep.
        
        Parameters
        ----------
        path : str
            The path to which to write the VTK files.
        name : str
            The prefix of the VTK filenames.
        id_flag : Optional[bool]
            If True, each particle is assigned a unique ID so that it can be
            tracked. If the data contains a node_ids column already, then these
            will be written to the VTK file regardless of id_flag.
        num_blade_elems : Optional[list]
            The number of elements for each blade.
        print_status: bool
            True to print the status of VTK conversion, False to suppress.

        Returns
        -------
        data_filenames : list
            List of the VTK filenames.
        pvd_filename : str
            .pvd collection filename.           
        """


        from pyevtk.hl import pointsToVTK   # evtk module - import only if this function is called
        import xml.etree.cElementTree as ET # xml module  -                 "

        # set the collection filename
        collection_fname = name + ".pvd"

        # set up blank list of the vtk filenames
        data_filenames = []

        # set up XML tree for PVD collection file
        root = ET.Element("VTKFile")
        root.set("type", "Collection")
        collection = ET.SubElement(root, "Collection")

        # compute number of wake element nodes
        if id_flag:
            num_blade_nodes = np.array(num_blade_elems) + 1

        # loop through all the files, write the VTK files
        for i, time in enumerate(np.sort(self.times)):
            # get the system time (for elapsed time)
            t_start = pytime.time()

            # get the filename containing the data at current time
            fname = self.fdict[time]

            # base name of data file
            vtk_name = name + '_' + str(i)

            ###### Read in the data ######
            df_inst                   = self.get_df_inst(time=time)
            data_arrays, has_node_ids = self.wakedata_from_df(df_inst)

            # unpack the data from dictionary
            elems = data_arrays['elems']
            x = np.float32(data_arrays['x'])
            y = np.float32(data_arrays['y'])
            z = np.float32(data_arrays['z'])
            u = np.float32(data_arrays['u'])
            v = np.float32(data_arrays['v'])
            w = np.float32(data_arrays['w'])

            # if the data has node ids already, load the data
            if has_node_ids:
                node_ids = data_arrays['node_ids']

            # if id_flag is set, generate new node id numbers. Note that this will overwrite loaded node ids.
            if id_flag:
                ###### Generate node IDs ######
                # generate the id numbers for the elements as an array with length of num_wake_elems
                # newer elements have a higher id number.
                
                # compute the number of wake elements at this particular timestep
                num_wake_elems = len(u)

                # compute the timestep number from the number of elements
                elems_per_timestep = sum(num_blade_nodes)
                nt = num_wake_elems/elems_per_timestep

                # generate the id numbers
                multiplier = np.mod((np.arange(num_wake_elems)), nt)
                adder      = multiplier * elems_per_timestep
                node_ids   = adder + elems

            else:
                # otherwise, just use the node_ids that we loaded from the data set before
                pass

            # convert to a uint32 (fixes some problems importing into VTK/ParaView)
            node_ids = np.uint32(node_ids)


            # store vector field in a dict
            data = {'velocity' : (u,
                                  v,
                                  w),
                    'node_id' : node_ids}

            # write data
            data_filename = pointsToVTK(os.path.abspath(os.path.join(path,vtk_name)), x, y, z, data)

            # append filename to list
            data_filenames.append(data_filename)

            # add elements to XML tree for PVD collection file
            dataset = ET.SubElement(collection, "DataSet")
            dataset.set("timestep", str(time))
            dataset.set("file", os.path.basename(data_filename))

            # print status message
            elapsed_time = pytime.time() - t_start
            if print_status:
                print 'Converted: ' + fname + ' -->\n\t\t\t' + data_filename + ' in %2.2f s\n' % (elapsed_time)

        # write the collection file
        tree = ET.ElementTree(root)
        pvd_filename = os.path.abspath(os.path.join(path,collection_fname))
        tree.write(pvd_filename, xml_declaration=True)

        if print_status:
            print 'Wrote ParaView collection file: ' + pvd_filename

        return data_filenames, pvd_filename
