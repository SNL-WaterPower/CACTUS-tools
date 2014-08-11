# pyCactusWakeGrid.py
""" Class and functions for manipulating Cartesian wake induced velocity data from CACTUS."""

import numpy as np
import scipy.integrate

import matplotlib.pyplot as plt
from warnings import *

class WakeGridData():
	""" Class which loads WakeGridData from a pandas dataframe and creates appropriately-shaped
		Numpy arrays. """

	def __init__(self, df):
		time_col_name = 'Normalized Time (-)'
		x_col_name = 'X/R (-)'
		y_col_name = 'Y/R (-)'
		z_col_name = 'Z/R (-)'
		u_col_name = 'U/Uinf (-)'
		v_col_name = 'V/Uinf (-)'
		w_col_name = 'W/Uinf (-)'

		# get unique times
		self.times = df.loc[:,time_col_name].unique()

		# get number of times
		self.num_times = len(self.times)

		# extract the data
		x  = df.loc[:,x_col_name]
		y  = df.loc[:,y_col_name]
		z  = df.loc[:,z_col_name]
		u  = df.loc[:,u_col_name]
		v  = df.loc[:,v_col_name]
		w  = df.loc[:,w_col_name]

		self.nt = self.num_times
		self.nx = len(x.unique())
		self.ny = len(y.unique())
		self.nz = len(z.unique())

		# reshape into 4-D numpy array
		# note that in Python, the final index is the fastest changing
		self.x = np.reshape(x, [self.nt, self.nz, self.ny, self.nx])
		self.y = np.reshape(y, [self.nt, self.nz, self.ny, self.nx])
		self.z = np.reshape(z, [self.nt, self.nz, self.ny, self.nx])
		self.u = np.reshape(u, [self.nt, self.nz, self.ny, self.nx])
		self.v = np.reshape(v, [self.nt, self.nz, self.ny, self.nx])
		self.w = np.reshape(w, [self.nt, self.nz, self.ny, self.nx])


#####################################
######### Module  Functions #########
#####################################
def find_nearest(array, value):
	""" find_nearest() : returns the value and index of the nearest point in a 1-d array. """
	index = (np.abs(array-value)).argmin()
	nearest_value = array[index]

	return nearest_value, index


def field_time_stats(field, times, t_start, t_end):
	""" field_time_stats() : provides time statistics for a scalar field (ordered with time as 
			the first index) from t_start to t_end. The times are specified with the 1-D array
			`times'. The scalar field may be of any dimension.""" 

	if t_start < times[0]:
		warn('The specified start time %f for the averaging operation is outside of the range of \
				the data. Averaging from beginning.' % t_start)
		index_start = 0

	else:
		# find the index corresponding to t_start
		nearest_t_start, index_start = find_nearest(times, t_start)	
	
	if t_end > times[-1]:
		warn('The specified end time %f for the averaging operation is outside of the range of the \
				data. Averaging to end.' % t_end)
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
	""" get_line_quantity() : Extracts the value of a scalar field U stored at locations X,Y
		on a vertical line running nearest to x_loc. """

	# find the nearest location to the desired slice
	x_grid = np.unique(X)
	x_nearest, idx_nearest = find_nearest(x_grid, x_loc)

	return U[:,idx_nearest], x_nearest, idx_nearest


def calc_momentum_def(x_loc, X, Y, U):
	""" calc_momentum_def() : Calculates the integral momentum deficit of scalar field U stored at \
		locations X,Y on a vertical line that runs nearest to x_loc. """
	
	U_line, x_line, x_idx_line = get_line_quantity(x_loc, X, Y, U)
	y_line = Y[:,x_idx_line]

	return scipy.integrate.trapz(U_line*(1-U_line), y_line)


def calc_tke(U, V, W, times, t_start, t_end):
	""" calc_tke() : Calculates the turbulent kinetic energy of a velocity field given by scalar
			fields U,V,W on a plane with coordinates X,Y. Note that U,V,W are a time series, with
			time as the first index. The time statistics are computed only from t_start to t_end.

			Returns the TKE as a scalar field.

			http://en.wikipedia.org/wiki/Turbulence_kinetic_energy"""

	# calculate the mean flow of each component
	U_mean, _, _, _ = field_time_stats(U, times, t_start, t_end)
	V_mean, _, _, _ = field_time_stats(V, times, t_start, t_end)
	W_mean, _, _, _ = field_time_stats(W, times, t_start, t_end)

	# calculate the fluctation of each velocity component
	U_fluc = U_mean - U
	V_fluc = V_mean - V
	W_fluc = W_mean - W

	# compute the TKE
	tke = 0.5*(np.mean(np.power(U_fluc,2),axis=0) +  
			np.mean(np.power(V_fluc,2),axis=0) +  
			np.mean(np.power(W_fluc,2),axis=0))

	return tke