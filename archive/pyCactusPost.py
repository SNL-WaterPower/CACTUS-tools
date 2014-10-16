# pyCactusPost.py

import numpy as np
from warnings import *

def rotor_quant_inst(CactusRun, name, t):
	""" rotor_quant_inst() : returns the instantaneous rotor quantity at normalized time t. """

	# dictionary for get_named_subset options
	name_to_subset = { 'Cp'  : 'time_Cp',
					   'CFx' : 'time_CFx'	}

	subset_name = name_to_subset[name]

	# get the data series
	columns, col_names = CactusRun.get_named_integral_subset(subset_name)

	# break out data
	time  = columns[0]
	quant = columns[1]

	# if t == -1, return the latest value
	if t == -1:
		quant_at_t = quant[-1]

	# otherwise, interpolate in time
	else:
		# warn if the specified time is outside of the data range
		if t < min(time) or t > max(time):
			warn('The specified time is outside of the range of the data. Data at the endpoint is returned.')
		quant_at_t = np.interp(t, time, quant)

	return quant_at_t


def rotor_quant_rev(CactusRun, name, rev):
	""" rotor_quant_inst() : returns the revolution-averaged rotor quantity at a specified revolution rev. """

	# dictionary for get_named_subset options
	name_to_subset = { 'Cp'  : 'time_Cp',
					   'CFx' : 'time_CFx'	}

	subset_name = name_to_subset[name]

	# get the data series
	columns, col_names = CactusRun.get_named_integral_subset(subset_name)

	# break out data
	revs  = columns[0]
	quant = columns[1]	

	# interpolate in revs (note that we would expect to have data for the rev we specify)
	quant_at_rev = np.interp(rev, revs, quant)

	return quant_at_rev


def blade_distribution_inst(CactusRun, name, t):
	""" blade_loads_inst() : returns the distributions of an instantaneous blade quantity at normalized time t. """

	# dictionary for get_named_subset options
	name_to_subset = { 'CFx' : 'elem_time_CFx',		# contribution of element to force coefficient CFx
					   'CTq' : 'elem_time_CTq',		# contribution of element to torque coefficient CTq
					   'Ma' : 'elem_time_Ma', 		# Mach number
					   'Ur' : 'elem_time_Ur', 		# Local velocity
					   'alpha' : 'elem_time_AoA_25', # Local angle of attack
   					   'GB' : 'elem_time_GB', 		# Local bound vorticity

   					 }	

	subset_name = name_to_subset[name]

	# get the data
	time, blade_data, col_names = CactusRun.get_named_elem_subset(subset_name)

	# interpolate in time and return the data
	return interp_blade_distribution(t, time, blade_data)


def blade_distribution_rev(CactusRun, name, rev):
	""" blade_loads_inst() : returns the distributions of a revolution-averaged blade quantity at
			revolution rev. """

	# dictionary for get_named_subset options
	name_to_subset = { 'Cp' : 'rev_Cp',			# power coefficient
					   'CTq' : 'rev_CTq',		# torque coefficient
					   'CFx' : 'rev_CFx', 		# force coefficient, x-dir
					   'CFy' : 'rev_CFy', 		# force coefficient, y-dir
					   'CFz' : 'rev_CFz', 		# force coefficient, z-dir
					   }	
	return


def interp_blade_distribution(ind_var, ind_var_data, blade_data):
	""" interp_blade_distribution() : interpolates blade data in the independent variable (typically time
		or revolution).	Returns a list of distributions (list of np.array) for each blade."""

	if ind_var < min(ind_var_data) or ind_var > max(ind_var_data):
		warn('The specified time is outside of the range of the data. Data at the endpoint is returned.')

	interpolated_data = []
	for data in blade_data:
		num_rows, num_cols = data.shape

		temp = np.zeros([0])
		# the inner loop can be replaced with scipy.interpolate.interp1d() if available
		for col_num in range(num_cols):
			temp = np.hstack( (temp, np.interp(ind_var, ind_var_data, data[:,col_num]) ) )

		interpolated_data.append(temp)

	return interpolated_data