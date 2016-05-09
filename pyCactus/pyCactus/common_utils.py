"""Common utilities for pyCactus."""


def recursive_glob(rootdir='.', pattern='*'):
    """Search recursively for files matching a specified patterns.

    Arguments
    ---------
    rootdir : str, optional
        Path to directory to search (default is '.')
    pattern : str, optional
        Glob-style pattern (default is '*')

    Returns
    -------
    matches : list
        List of matching files.
    """
    # Adapted from http://stackoverflow.com/questions/2186525/
    import os
    import fnmatch

    matches = []
    for root, dirnames, filenames in os.walk(rootdir):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))

    return matches


def load_data(data_filename):
    """Read data from a CSV file and return as Pandas dataframe."""
    import pandas as pd

    # read a CSV file using pandas and returns a pandas dataframe
    reader = pd.read_csv(data_filename, iterator=True, chunksize=1000)
    df = pd.concat(reader, ignore_index=True)

    # strip whitespace from colnames
    df.rename(columns=lambda x: x.strip(), inplace=True)
    return df


def df_subset_time_index(df, time_index, time_col_name):
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


def get_file_time(data_filename, time_col_name):
    """Return the time of an instantaneous data set.

    Read the first two rows of a data file in order to get the time.
    Expects that the normalized time is in the first data column.

    Arguments
    ---------
    data_filename : str
        Path to data file.
    time_col_name : str
        Name of column which contains time data.

    Returns
    -------
    time : float
        Value of time in first data row.
    """
    with open(data_filename) as f:
        header = f.readline()
        row1   = f.readline()

        time_col_num = header.split(',').index(time_col_name)
        time = float(row1.split(',')[time_col_num])

    return time
