# pyCactusProbes.py
"""Module for handling probe data from a CACTUS simulation."""

import glob
import pandas as pd
import matplotlib.pyplot as plt

class Probes():
    def __init__(self):
        """Initialize Probes class."""
        self.probe_locations = {}
        self.probe_filenames = {}


    def read_probe_files(self, run_directory):
        """Find probe files in run directory and read the header."""

        # find files in run_directory which match probe*.csv
        probe_filenames = glob.glob(run_directory + '/probe*.csv')

        # read in the location of each probe (not the data)
        # add the probe location and filenames
        for probe_id,probe_filename in enumerate(sorted(probe_filenames)):
            with open(probe_filename) as f:
                # read the location of the probe
                tmp = f.readline()
                splitline = f.readline().split(',')
                
                if len(splitline) == 3:
                    x,y,z = splitline
                elif len(splitline) == 4:
                    x,y,z,tmp = splitline
                    
                # store probe location and filename to dictionary
                self.probe_locations[probe_id] = (float(x),float(y),float(z))
                self.probe_filenames[probe_id] = probe_filename


    def get_probe_data_by_id(self, probe_id):
        """Return a Pandas dataframe of the probe data with given id."""

        probe_filename = self.probe_filenames[probe_id]
        return pd.read_csv(probe_filename,skiprows=2)


    def plot_probe_data_by_id(self, probe_id, timestep=False, plot_fs=False):
        """Plot the time-series velocities a probe with given id.

        Keyword arguments:
        timestep -- True to plot x-axis as timestep instead of normalized time (default False)
        plot_fs -- True to plot freestream velocities in addition to induced velocities (default False)
        """

        df = self.get_probe_data_by_id(probe_id)
        
        plt.figure()

        # plot velocities with either timestep or normalized time for x-axis
        if timestep:
            plt.plot(df['U/Uinf (-)'], label='U/Uinf (-)')
            plt.plot(df['V/Uinf (-)'], label='V/Uinf (-)')
            plt.plot(df['W/Uinf (-)'], label='W/Uinf (-)')
            plt.xlabel('Timestep')
        else:
            plt.plot(df['Normalized Time (-)'], df['U/Uinf (-)'], label='U/Uinf (-)')
            plt.plot(df['Normalized Time (-)'], df['V/Uinf (-)'], label='V/Uinf (-)')
            plt.plot(df['Normalized Time (-)'], df['W/Uinf (-)'], label='W/Uinf (-)')
            plt.xlabel('Normalized Time (-)')

        # plot the free-stream (if requested)
        if plot_fs:
            plt.plot(df['Ufs/Uinf (-)'], label='Ufs/Uinf (-)')
            plt.plot(df['Vfs/Uinf (-)'], label='Vfs/Uinf (-)')
            plt.plot(df['Wfs/Uinf (-)'], label='Wfs/Uinf (-)')

        plt.ylabel('Normalized Velocity (-)')
        plt.title(r'Velocity at $p=(%2.2f,%2.2f,%2.2f)$' % (self.probe_locations[probe_id]))
        plt.grid(True)
        plt.legend()
        plt.show()