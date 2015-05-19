# pyCactus.py
""" A module for parsing CACTUS run (input and output) files. """

import os
import glob
import fnmatch
import time as pytime

import numpy as np
import pandas as pd
import f90nml

from pyCactusGeom import *
from pyCactusWake import *

class CactusRun():
	def __init__(self, run_directory, case_name, input_fname='', geom_fname='', load_wake_node=True, load_wake_grid=True):
		# input file
		if not input_fname:
			self.input_fname = case_name + '.in'
		else:
			self.input_fname = input_fname

		# read in the input file namelist
		
		if os.path.isfile(self.input_fname):
			self.namelist = f90nml.read(self.input_fname)
		else:
			print 'Warning: Could not find file %s' % (self.input_fname)

		# if a geometry filename is specified, use that. otherwise, assume its [case_name].geom
		if not geom_fname:
			self.geom_fname = case_name + '.geom'
		else:
			self.geom_fname = geom_fname

		## assemble filenames
		self.geom_filename      = os.path.abspath(run_directory + '/' + self.geom_fname)
		self.elem_filename      = os.path.abspath(run_directory + '/' + case_name + '_ElementData.csv')
		self.param_filename     = os.path.abspath(run_directory + '/' + case_name + '_Param.csv')
		self.rev_filename       = os.path.abspath(run_directory + '/' + case_name + '_RevData.csv')
		self.time_filename      = os.path.abspath(run_directory + '/' + case_name + '_TimeData.csv')
		
		# search for wake data files anywhere in the directory
		self.wake_filenames     = self.recursive_glob(run_directory, '*WakeData_*.csv')
		self.wakegrid_filenames = self.recursive_glob(run_directory, '*WakeDefData_*csv')
		

		## load data
		# elem_data
		if os.path.isfile(self.elem_filename):
			self.elem_data  = self.load_data(self.elem_filename)
		else:
			print 'Warning: File ' + self.elem_filename + ' does not exist.'

		# rev_data
		if os.path.isfile(self.rev_filename):
			self.rev_data  = self.load_data(self.rev_filename)
		else:
			print 'Warning: File ' + self.rev_filename + ' does not exist.'
			
		# param_data
		if os.path.isfile(self.param_filename):
			self.elem_data  = self.load_data(self.param_filename)
		else:
			print 'Warning: File ' + self.param_filename + ' does not exist.'			
			
		# time_data
		if os.path.isfile(self.time_filename):
			self.elem_data  = self.load_data(self.time_filename)
		else:
			print 'Warning: File ' + self.time_filename + ' does not exist.'
			

		## load wake data
		if self.wake_filenames:
			try:
				tic = pytime.time() 
				self.wakeelems = CactusWakeElems(self.wake_filenames)
				print 'Loaded wake element data in %2.2f s' % (pytime.time() - tic)
			except:
				print 'Warning: Wake node data was found, but could not be loaded properly.'
		else:
			print 'Warning: Could not find any wake data files in the work directory matching \'*WakeData_*.csv\'.'

		if self.wakegrid_filenames:
			try:
				tic = pytime.time()
				self.wakegrid = CactusWakeGrid(self.wakegrid_filenames)
				print 'Loaded wake grid data in %2.2f s' % (pytime.time() - tic)
			except:
				print 'Warning: Wake grid data was found, but could not be loaded properly.'
		else:
			print 'Warning: Could not find any wake grid data files in the work directory matching \'*WakeGridData_*.csv\'.'


		## load geometry data
		if os.path.isfile(self.geom_filename):
			self.geom = CactusGeom(self.geom_filename)
		else:
			print 'Warning: File ' + self.geom_filename + ' does not exist.'


		print 'Success: Loaded case `%s` from path `%s`\n' % (case_name, run_directory)


	#####################################
	######### Private Functions #########
	#####################################
	def recursive_glob(self, rootdir='.', pattern='*'):
		""" A function to search recursively for files matching a specified pattern.
			Adapted from http://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python """

		matches = []
		for root, dirnames, filenames in os.walk(rootdir):
		  for filename in fnmatch.filter(filenames, pattern):
		      matches.append(os.path.join(root, filename))

		return matches

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
		""" Extracts a subset dataframe of the "Element Data" dataframe by time index.
			Returns the time corresponding to the given time_index, and a list of dataframes 
			containing the data, with one dataframe per blade.

			Acts on: Element Data """

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

	def rotor_data_at_time_index(self, time_index):
		""" Extracts a subset dataframe of the "Time Data" dataframe by time index.
			Returns the time corresponding to the given time_index, and a dataframe 
			containing the appropriate subset of data.

			Acts on: Time Data """

		# get the data series
		df = self.time_data

		# set column names
		time_col_name  = 'Normalized Time (-)'
		
		# extract data subset
		df, time = self.df_subset_time_index(df, time_index, time_col_name)

		return time, df

	def df_subset_time_index(self, df, time_index, time_col_name):
		""" Extracts a subset dataframe of the given dataframe by time index.
			Returns the dataframe subset and the time corresponding to the given time_index. """

		# get unique times
		times = df.loc[:,time_col_name].unique()

		# time at which we wish to extract data
		time = times[time_index]

		# extract the subset of data corresponding to the desired time_index
		df = df[df[time_col_name] == times[time_index]]

		return df, time