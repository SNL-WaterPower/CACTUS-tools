#!/usr/bin/python
# pyCactus.py : A module for parsing CACTUS run (input and output) files.

import numpy as np
import pandas as pd

from pyCactusGeom import *

class CactusRun():
	def __init__(self, run_directory, case_name, input_fname='', geom_fname=''):
		# input file
		if not input_fname:
			self.input_fname = case_name + '.in'
		else:
			self.input_fname = input_fname

		if not geom_fname:
			self.geom_fname = case_name + '.geom'
		else:
			self.geom_fname = geom_fname

		# create geometry instance (reads in geometry variables)
		self.geom = CactusGeom(run_directory + '/' + self.geom_fname)

		# output files
		self.elem_fname    = case_name + '_ElementData.csv'
		self.param_fname   = case_name + '_Param.csv'
		self.rev_fname     = case_name + '_RevData.csv'
		self.time_fname    = case_name + '_TimeData.csv'
		self.wake_fname    = case_name + '_WakeData.csv'
		self.wakedef_fname = case_name + '_WakeDefData.csv'
		
		# read data from CSV-formatted output files, incorporating try/except for possible
		# file read errors.

		# element data
		try:
			self.elem_data = self.load_data(run_directory + '/' + self.elem_fname)
		except:
			print "Could not load file : ", run_directory + '/' + self.elem_fname

		# parameters
		try:
			self.param_data = self.load_data(run_directory + '/' + self.param_fname)
		except:
			print "Could not load file : ", run_directory + '/' + self.param_fname

		# revolution-averaged data
		try:
			self.rev_data = self.load_data(run_directory + '/' + self.rev_fname)
		except:
			print "Could not load file : ", run_directory + '/' + self.rev_fname

		# time data
		try:
			self.time_data = self.load_data(run_directory + '/' + self.time_fname)
		except:
			print "Could not load file : ", run_directory + '/' + self.time_fname

		# wake data
		try:
			self.wake_data = self.load_data(run_directory + '/' + self.wake_fname)
		except:
			print "Could not load file : ", run_directory + '/' + self.wake_fname

		# wake element data
		try:
			self.wakedef_data = self.load_data(run_directory + '/' + self.wakedef_fname)
		except:
			print "Could not load file : ", run_directory + '/' + self.wakedef_fname

	#####################################
	######### Private Functions #########
	#####################################
	def load_data(self, data_filename):
		# reads a CSV file using pandas and returns a pandas dataframe
		df = pd.read_csv(data_filename)
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
		df = df[df['Normalized Time (-)'] == times[time_index]]

		return df, time