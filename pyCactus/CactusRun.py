# CactusRun.py
"""A module for parsing CACTUS run (input and output) files."""

import os
import glob
import time as pytime

import numpy as np
import pandas as pd
import f90nml

from CactusGeom import CactusGeom
from CactusWakeElems import CactusWakeElems
from CactusField import CactusField
from CactusProbes import CactusProbes

from recursive_glob import recursive_glob

class CactusRun():
	def __init__(self, run_directory, case_name, input_fname='', geom_fname=''):
		"""Initialize the class. Reads some data to memory.

		Parameters
		----------
		run_directory : str
			Path to the directory containing the CACTUS run.
		case_name : str
			'case name' which precedes all input and output files.
		input_fname : Optional[str]
			Input filename (default: ./[case_name].in).
		geom_fname : Optional[str]
			Geometry filename (default: be ./[case_name].geom)
		"""


		# input file
		if not input_fname:
			self.input_fname = case_name + '.in'
		else:
			self.input_fname = input_fname

		# if a geometry filename is specified, use that. otherwise, assume its [case_name].geom
		if not geom_fname:
			self.geom_fname = case_name + '.geom'
		else:
			self.geom_fname = geom_fname

		## assemble filenames
		self.elem_fname      = case_name + '_ElementData.csv'
		self.param_fname     = case_name + '_Param.csv'
		self.rev_fname       = case_name + '_RevData.csv'
		self.time_fname      = case_name + '_TimeData.csv'
		
		## search for wake data files anywhere in the directory
		self.wake_filenames     = sorted(recursive_glob(run_directory,'*WakeData_*.csv'))
		self.field_filenames    = sorted(recursive_glob(run_directory, '*WakeDefData_*csv'))
		
		## read in the input file namelist
		results = recursive_glob(run_directory, self.input_fname)
		if results:
			self.namelist = f90nml.read(results[0])
		else:
			print 'Warning: Could not find file %s in %s' % (self.input_fname, run_directory)

		## load data
		# blade element data
		results = recursive_glob(run_directory, self.elem_fname)
		if results:
			self.elem_data  = self.load_data(results[0])
		else:
			print 'Warning: Could not find file %s in %s' % (self.elem_fname, run_directory)

		# revolution-averaged data
		results = recursive_glob(run_directory, self.rev_fname)
		if results:
			self.rev_data  = self.load_data(results[0])
		else:
			print 'Warning: Could not find file %s in %s' % (self.rev_fname, run_directory)
			
		# parameter data
		results = recursive_glob(run_directory, self.param_fname)
		if results:
			self.param_data  = self.load_data(results[0])
		else:
			print 'Warning: Could not find file %s in %s' % (self.param_fname, run_directory)
			
		# time data
		results = recursive_glob(run_directory, self.time_fname)
		if results:
			self.time_data  = self.load_data(results[0])
		else:
			print 'Warning: Could not find file %s in %s' % (self.time_fname, run_directory)
			

		# wake element data
		if self.wake_filenames:
			try:
				tic = pytime.time() 
				self.wakeelems = CactusWakeElems(self.wake_filenames)
				print 'Read wake element data headers in %2.2f s' % (pytime.time() - tic)
			except:
				print 'Warning: Wake node data was found, but could not be loaded properly.'
		else:
			print 'Warning: Could not find any wake data files in the work directory matching \'*WakeData_*.csv\'.'

		# field data
		if self.field_filenames:
			try:
				tic = pytime.time()
				self.field = CactusField(self.field_filenames)
				print 'Read wake grid data headers in %2.2f s' % (pytime.time() - tic)
			except:
				print 'Warning: Wake grid data was found, but could not be loaded properly.'
		else:
			print 'Warning: Could not find any wake grid data files in the work directory matching \'*WakeGridData_*.csv\'.'


		# probe data
		tic = pytime.time() 
		self.probes = CactusProbes()
		self.probes.read_probe_files(run_directory)
		print 'Read probe data in %2.2f s' % (pytime.time() - tic)


		# geometry data
		results = recursive_glob(run_directory, self.geom_fname)
		if results:
			self.geom = CactusGeom(results[0])
		else:
			print 'Warning: Could not find file %s in %s' % (self.geom_fname, run_directory)


		print 'Success: Loaded case `%s` from path `%s`\n' % (case_name, run_directory)


	#####################################
	######### Private Functions #########
	#####################################
	def load_data(self, data_filename):
		# reads a CSV file using pandas and returns a pandas dataframe
		reader = pd.read_csv(data_filename, iterator=True, chunksize=1000)
		df = pd.concat(reader, ignore_index=True)
		
		df.rename(columns=lambda x: x.strip(), inplace=True)	# strip whitespace from colnames
		return df


	####################################
	######### Public Functions #########
	####################################
	def blade_data_at_time_index(self, time_index):
		"""Extracts a subset dataframe of the "Element Data" dataframe by time
		index.

		Returns the time corresponding to the given time_index, and a list of
		dataframes containing the data, with one dataframe per blade.

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
		num_elems  = [blade['NElem'][0] for blade in self.geom.blades]

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
		"""Computes the average value of a specified blade quantity over a
		number of timesteps from the end.

		This may be used to, for example, compute the revolution-averaged blade 
		quantities, such as circulation distribution, angle of attack, and local
		blade relative velocity.

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
		time_indices = range(-1,-timesteps_from_end,-1)
		
		# initialize empty array to store circulation at final timestep
		qty_avg = np.zeros(self.geom.blades[blade_num]['NElem'],)
		
		# compute average 
		for time_index in time_indices:
			# get time (float), and a tuple containing Pandas dataframes of length num_blades.
			time, elem_data = self.blade_data_at_time_index(time_index)

			# get the circulation distribution at the final timestep
			qty = elem_data[0][qty_name]
			qty_avg = qty_avg + qty.values
		
		qty_avg = qty_avg/((len(time_indices)))
		
		# return the averaged quantity
		return qty_avg
	
	def rotor_data_at_time_index(self, time_index):
		"""Extracts a single time instance from the time dataframe.
		
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
		"""Extracts a subset dataframe of a given dataframe by time index.

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
		times = df.loc[:,time_col_name].unique()

		# time at which we wish to extract data
		time = times[time_index]

		# extract the subset of data corresponding to the desired time_index
		df = df[df[time_col_name] == times[time_index]]

		return df, time

