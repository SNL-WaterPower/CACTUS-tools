#!/usr/bin/python
# pyCactus.py : A module for parsing CACTUS run (input and output) files.

import numpy as np

class CactusRun():
	def __init__(self, run_directory, case_name):
		# input file
		self.input_fname   = case_name + '.in'

		# output files
		self.elem_fname    = case_name + '_ElementData.csv'
		self.param_fname   = case_name + '_Param.csv'
		self.rev_fname     = case_name + '_RevData.csv'
		self.time_fname    = case_name + '_TimeData.csv'
		self.wake_fname    = case_name + '_WakeData.csv'
		self.wakedef_fname = case_name + '_WakeDefData.csv'
		
		# read data from CSV-formatted output files (may need to change to binary files later...)
		self.elem_data     = self.__load_data(run_directory + '/' + self.elem_fname)
		self.param_data    = self.__load_data(run_directory + '/' + self.param_fname)
		self.rev_data      = self.__load_data(run_directory + '/' + self.rev_fname)
		self.time_data     = self.__load_data(run_directory + '/' + self.time_fname)
		self.wake_data     = self.__load_data(run_directory + '/' + self.wake_fname)
		self.wakedef_data  = self.__load_data(run_directory + '/' + self.wakedef_fname)

	#####################################
	######### Private Functions #########
	#####################################
	def __get_col_by_name(self, col_name, data):
		data_array = data[0]
		data_headers = data[1]

		col_num = data_headers[col_name]
		return data_array[:,col_num]


	def __get_col_by_names(self, col_names, data):
		data_columns = []
		for col_name in col_names:
			data_columns.append(self.__get_col_by_name(col_name, data))
		
		return data_columns


	def __get_col_by_number(self, col_num, data):
		data_array = data[0]
		return data_array[:,col_num]


	def __col_stats(self, data_column):
		return dict({	'mean' : np.mean(data_column),
						'std'  : np.std(data_column),
						'min'  : np.min(data_column),
						'max'  : np.max(data_column)	})


	def __load_data(self, data_filename):
		# reads data from a file. returns a dictionary pairing headers/column numbers, and 
		# a numpy array of the data
		with open(data_filename) as f:
			header_txt = f.readline()

		# read the header row
		headers = header_txt.split(',')
		headers = [header.strip() for header in headers]

		# create dictionary
		data_headers = dict([[header, col_num] for col_num, header in enumerate(headers)])

		# load the data
		data_array = np.loadtxt(data_filename,skiprows=1,delimiter=',')

		return data_array, data_headers


	####################################
	######### Public Functions #########
	####################################
	def get_named_subset(self, subset_name):
		# extracts a predefined subset of data identified by a string
		# returns the data subset (list of np.arrays) and the column headers (list of strings)

		options = {	'rev_Cp' : [ self.rev_data, ['Rev', 'Power Coeff. (-)'] ] ,
					'time_Cp': [ self.time_data, ['Normalized Time (-)','Power Coeff. (-)'] ] 
					# '': '', # add more as necessary
					# '': '',
					# '': '',
					# '': '',
					}

		data, col_names = options[subset_name]
		return self.__get_col_by_names(col_names, data), col_names