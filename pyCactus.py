#!/usr/bin/python
# pyCactus.py : A module for parsing CACTUS run (input and output) files.

import numpy as np
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
		try: self.elem_data     = self.load_data(run_directory + '/' + self.elem_fname)
		except:	print "Could not load file : ", run_directory + '/' + self.elem_fname

		# parameters
		try: self.param_data    = self.load_data(run_directory + '/' + self.param_fname)
		except: print "Could not load file : ", run_directory + '/' + self.param_fname

		# revolution-averaged data
		try: self.rev_data      = self.load_data(run_directory + '/' + self.rev_fname)
		except: print "Could not load file : ", run_directory + '/' + self.rev_fname

		# time data
		try: self.time_data     = self.load_data(run_directory + '/' + self.time_fname)
		except: print "Could not load file : ", run_directory + '/' + self.time_fname

		# wake data
		try: self.wake_data     = self.load_data(run_directory + '/' + self.wake_fname)
		except: print "Could not load file : ", run_directory + '/' + self.wake_fname

		# wake element data
		try: self.wakedef_data  = self.load_data(run_directory + '/' + self.wakedef_fname)
		except: print "Could not load file : ", run_directory + '/' + self.wakedef_fname



	#####################################
	######### Private Functions #########
	#####################################
	def get_col_by_name(self, col_name, data):
		data_array = data[0]
		data_headers = data[1]

		col_num = data_headers[col_name]
		return data_array[:,col_num]

	def get_col_by_names(self, col_names, data):
		data_columns = []
		for col_name in col_names:
			data_columns.append(self.get_col_by_name(col_name, data))
		
		return data_columns

	def get_col_num_by_name(self, col_name, data):
		data_headers = data[1]
		return data_headers[col_name]

	def get_col_by_number(self, col_num, data):
		data_array = data[0]
		return data_array[:,col_num]

	def col_stats(self, data_column):
		return dict({	'mean' : np.mean(data_column),
						'std'  : np.std(data_column),
						'min'  : np.min(data_column),
						'max'  : np.max(data_column)	})

	def load_data(self, data_filename):
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
		data_array = np.loadtxt(data_filename, skiprows=1, delimiter=',')

		return data_array, data_headers


	####################################
	######### Public Functions #########
	####################################
	def get_named_integral_subset(self, subset_name):
		""" Extracts a predefined subset of "integral" (blade or rotor-integrated) data identified
			by a string returns the data subset (list of np.arrays) and the column headers (list of
			strings) """

		options = {	'time_Cp'	: [ self.time_data, ['Normalized Time (-)', 'Power Coeff. (-)'  ] ], 
					'time_CTq'	: [ self.time_data, ['Normalized Time (-)', 'Torque Coeff. (-)' ] ],
					'time_CFx'	: [ self.time_data, ['Normalized Time (-)', 'Fx Coeff. (-)'     ] ],
					'time_CFy'	: [ self.time_data, ['Normalized Time (-)', 'Fy Coeff. (-)'     ] ],
					'time_CFz'	: [ self.time_data, ['Normalized Time (-)', 'Fz Coeff. (-)'     ] ],
					'rev_Cp' 	: [ self.rev_data,  ['Rev'                , 'Power Coeff. (-)'  ] ],
					'rev_CTq'	: [ self.rev_data,  ['Rev'                , 'Torque Coeff. (-)' ] ],
					'rev_CFx' 	: [ self.rev_data,  ['Rev'                , 'Fx Coeff. (-)'     ] ],
					'rev_CFy' 	: [ self.rev_data,  ['Rev'                , 'Fy Coeff. (-)'     ] ],
					'rev_CFz' 	: [ self.rev_data,  ['Rev'                , 'Fz Coeff. (-)'     ] ],
					'rev_power' : [ self.rev_data,  ['Rev'                , 'Power (kW)'        ] ],
					'rev_torque' : [ self.rev_data, ['Rev'                , 'Torque (ft-lbs)'   ] ],
					# '': '',
					# '': '',
					}

		data, col_names = options[subset_name]
		return self.get_col_by_names(col_names, data), col_names

	def get_named_elem_subset(self, subset_name):
		""" Extracts a predefined subset of "element data" identified by a string.
			Returns the the 1-D independent variable data (an array of time or revolution data),
			the 2-D element data (a list of np.arrays, indexed by blade number), and the column headers
			(a single list of strings). """

		options = {	'elem_time_Re'		: [ self.elem_data, ['Normalized Time (-)', 'Re (-)'  	 ] ],	# element Reynolds number	
					'elem_time_CL'		: [ self.elem_data, ['Normalized Time (-)', 'CL (-)'  	 ] ],	# element lift coefficient
					'elem_time_CD'		: [ self.elem_data, ['Normalized Time (-)', 'CD (-)'  	 ] ],	# element drag coefficient
					'elem_time_CM'		: [ self.elem_data, ['Normalized Time (-)', 'CM25 (-)'   ] ],	# element moment coefficient
					'elem_time_CN'		: [ self.elem_data, ['Normalized Time (-)', 'CN (-)'  	 ] ],	# element normal force coeff
					'elem_time_CT'		: [ self.elem_data, ['Normalized Time (-)', 'CT (-)'  	 ] ],	# element tangential force coeff
					'elem_time_CFx'		: [ self.elem_data, ['Normalized Time (-)', 'Fx (-)'  	 ] ],	# element x-dir force coeff
					'elem_time_CFy'		: [ self.elem_data, ['Normalized Time (-)', 'Fy (-)'  	 ] ],	# element y-dir force coeff
					'elem_time_CFz'		: [ self.elem_data, ['Normalized Time (-)', 'Fz (-)'  	 ] ],	# element z-dir force coeff
					'elem_time_Ur'		: [ self.elem_data, ['Normalized Time (-)', 'Ur (-)'  	 ] ],	# element local velocity
					'elem_time_Ma'		: [ self.elem_data, ['Normalized Time (-)', 'Mach (-)'   ] ],	# element Mach number
					'elem_time_CTq'		: [ self.elem_data, ['Normalized Time (-)', 'te (-)'  	 ] ],	# element torque coeff
					'elem_time_AoA_25'	: [ self.elem_data, ['Normalized Time (-)', 'AOA25 (deg)'] ],	# element angle of attack at quarter-chord
					'elem_theta_IndU'	: [ self.elem_data, ['Theta (rad)', 'IndU (-)'] ],	# induced velocity, u
					'elem_theta_IndV'	: [ self.elem_data, ['Theta (rad)', 'IndV (-)'] ],	# induced velocity, v
					'elem_theta_IndW'	: [ self.elem_data, ['Theta (rad)', 'IndW (-)'] ],	# induced velocity, w
					'elem_theta_x'	: [ self.elem_data, ['Theta (rad)', 'x/R (-)'] ],	# element position, center of quarter-chord line
					'elem_theta_y'	: [ self.elem_data, ['Theta (rad)', 'y/R (-)'] ],	# element position, center of quarter-chord line
					'elem_theta_z'	: [ self.elem_data, ['Theta (rad)', 'z/R (-)'] ],	# element position, center of quarter-chord line
					'elem_time_GB'	: [ self.elem_data, ['Normalized Time (-)', 'GB (?)'] ],	# local bound vorticity
					# '': '',
					# '': '',
					}

		data, col_names = options[subset_name]

		# get number of blades 
		num_blades, num_elems = self.get_nblade_nelem()

		# get the data columns
		columns = self.get_col_by_names(col_names, data)

		# collapse independent variable data (first column)
		ind_var = np.unique(columns[0])

		# organize into lists by blade
		total_rows      = len(columns[0])			# get the number of rows
		num_total_elems = sum(num_elems)			# total number of elements
		num_blocks      = total_rows / num_total_elems	# number of time steps represented in the data

		# initialize empty lists
		blade_data = []
		blade_indices = []
		for blade in range(num_blades):
			blade_indices.append([])
			
		# get the indices corresponding to each blade
		start_row = 0
		for block in range(num_blocks):
			for blade in range(num_blades):
				blade_indices[blade] = blade_indices[blade] + (range(start_row, start_row + num_elems[blade]))
				start_row = start_row + num_elems[blade]

		# extract the appropriate data, reshape into a 2-D array indexed by independent variable and element number (in that order)
		for blade in range(num_blades):
			blade_data.append(np.reshape((columns[1])[blade_indices[blade]], (num_blocks, num_elems[blade])))

		return ind_var, blade_data, col_names
		

	def get_nblade_nelem(self):
		""" Gets the number of blades and number of elements per blade by parsing the blade element data."""
		# get the blades and elements rows
		blades, elements = self.get_col_by_names(['Blade', 'Element'], self.elem_data)

		# get the number of blades
		num_blades = len(np.unique(blades))

		# get the indices where the blade number changes
		change_indices = [ rownum for rownum, zipped in enumerate(zip(blades[:], blades[1:])) if zipped[0] != zipped[1]]

		# do a diff to get a list containing the number of elements per blade
		num_elems = (np.diff([-1] + change_indices))[:num_blades]

		return num_blades, num_elems