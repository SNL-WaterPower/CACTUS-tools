# pyCactusWake.py
""" Functions for manipulating wake data from CACTUS."""

import numpy as np
import scipy.integrate

from warnings import *

#######################
# Cartesian Wake Data #
#######################
def find_nearest(array, value):
	""" find_nearest() : returns the value and index of the nearest point in a 1-d array. """
	index = (np.abs(array-value)).argmin()
	nearest_value = array[index]

	return nearest_value, index


def field_time_stats(field, times, t_start, t_end):
	""" field_time_statistics() : provides time statistics for a scalar field (ordered with time as the first index)
		from t_start to t_end. The times are specified with the 1-D array `times'. The scalar field may be of any dimension.""" 

	if t_start < times[0]:
		warn('The specified start time %f for the averaging operation is outside of the range of the data. Averaging from beginning.' % t_start)
		index_start = 0
	else:
		# find the index corresponding to t_start
		nearest_t_start, index_start = find_nearest(times, t_start)	
	
	if t_end > times[-1]:
		warn('The specified end time %f for the averaging operation is outside of the range of the data. Averaging to end.' % t_end)
		index_end = -1
	elif t_end == -1:
		index_end = -1
	else:
		# find the index corresponding to t_end
		nearest_t_end  , index_end   = find_nearest(times, t_end)

	# compute statistics in time
	mean = np.mean(field[index_start:index_end, Ellipsis], axis=0)
	stdf = np.std(field[index_start:index_end, Ellipsis], axis=0)
	minf = np.min(field[index_start:index_end, Ellipsis], axis=0)
	maxf = np.max(field[index_start:index_end, Ellipsis], axis=0)

	return mean, stdf, minf, maxf

def get_line_quantity(x_loc, X, Y, U):
	""" Extracts the value of a scalar field U stored at locations X,Y
		on a vertical line running through x_loc. """

	# find the nearest location to the desired slice
	x_grid = np.unique(X)
	x_nearest, idx_nearest = find_nearest(x_grid, x_loc)

	return U[:,idx_nearest], x_nearest, idx_nearest


def calc_thrust_coeff(x_loc, X, Y, U):
	""" Calculates the integral velocity deficit of scalar field U stored at locations X,Y,
		on a line coincident with x_loc. """
	
	U_line, x_line, x_idx_line = get_line_quantity(x_loc, X, Y, U)
	y_line = Y[:,x_idx_line]

	return scipy.integrate.trapz(y_line, U_line)