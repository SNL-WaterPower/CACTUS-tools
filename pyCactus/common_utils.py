def recursive_glob(rootdir='.', pattern='*'):
    import os
    import fnmatch

    """ A function to search recursively for files matching a specified pattern.
        Adapted from http://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python """

    matches = []
    for root, dirnames, filenames in os.walk(rootdir):
      for filename in fnmatch.filter(filenames, pattern):
          matches.append(os.path.join(root, filename))

    return matches

def load_data(data_filename):
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
