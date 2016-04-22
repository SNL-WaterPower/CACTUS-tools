"""CactusBladeElem module containing CactusBladeElem class."""

import numpy as np

from common_utils import load_data, df_subset_time_index


class CactusBladeElem(object):
    """Class to manage CACTUS blade element data.

    Attributes
    ----------
    data : Pandas dataframe
        Raw blade element data.
    filename : str
        Path to raw blade element data file.
    """

    def __init__(self, filename):
        """Initialize class, read in data."""
        self.filename = filename
        self.data = load_data(self.filename)

    def data_at_time_index(self, time_index):
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
        df = self.data

        # extract data subset
        df, time = df_subset_time_index(df, time_index, time_col_name)

        # get number of blades
        num_blades = len(self.data[blade_col_name].unique())

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

    def data_time_average(self, qty_name, timesteps, blade_num):
        """Compute a time-averaged a blade quantity distribution.

        Compute a time-averaged blade quantity distribution over a number of
        timesteps. This may be used to, for example, compute the revolution-
        averaged blade quantities, such as circulation distribution, angle
        of attack, and local blade relative velocity.

        Parameters
        ----------
            qty_name : str
                The column name for the desired blade quantity.
            timesteps : list
                List of timesteps (integer).
            blade_num : int
                Index of the blade number, indexed from 0).

        Returns
        -------
            qty_avg : numpy.array
                The array of averaged blade data.

        """
        # set column names
        blade_col_name = 'Blade'
        elem_col_name  = 'Element'

        # get number of elements for the specified blade
        num_elems = len((self.data[elem_col_name]
                        [self.data[blade_col_name] == blade_num + 1]).unique())

        # initialize empty array to store circulation at final timestep
        qty_avg = np.zeros(num_elems)

        # compute average
        for time_index in timesteps:
            # get time (float), and a tuple containing Pandas dataframes of
            # length num_blades.
            time, elem_data = self.data_at_time_index(time_index)

            # get the circulation distribution at the final timestep
            qty = elem_data[blade_num + 1][qty_name]
            qty_avg = qty_avg + qty.values

        qty_avg = qty_avg / ((len(timesteps)))

        # return the averaged quantity
        return qty_avg
