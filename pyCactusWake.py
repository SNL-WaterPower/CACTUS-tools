# pyCactusWake.py
""" Class and functions for manipulating wake velocity data from CACTUS."""

import os

import numpy as np
import scipy.integrate

import matplotlib.pyplot as plt
from warnings import *


#####################################
######### Wake Element Data #########
#####################################
class CactusWakeElems():
	""" Class which loads WakeData (element) from a pandas dataframe and creates Numpy arrays
			to hold the data. """

	def __init__(self, df):
		time_col_name = 'Normalized Time (-)'
		id_col_name = 'Node ID'
		elem_col_name = 'Origin Node'
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
		
		# generate empty lists
		elem_list = []
		x_list = []
		y_list = []
		z_list = []
		u_list = []
		v_list = []
		w_list = []
		node_id_list = []

		# check if data has node IDs
		if id_col_name in df:
			has_node_ids = True
		else:
			has_node_ids = False
			print 'Warning: Wake element data does not include node IDs.'

		# generate a set of numpy arrays for each timestep
		for time in self.times:
			# slice dataframe
			df_slice = df[df[time_col_name] == time]

			# extract columns
			elem = df_slice.loc[:,elem_col_name].values
			x    = df_slice.loc[:,x_col_name].values
			y    = df_slice.loc[:,y_col_name].values
			z    = df_slice.loc[:,z_col_name].values
			u    = df_slice.loc[:,u_col_name].values
			v    = df_slice.loc[:,v_col_name].values
			w    = df_slice.loc[:,w_col_name].values
			
			if has_node_ids:
				node_id = df_slice.loc[:,id_col_name].values
			
			# append data to lists
			elem_list.append(elem)
			x_list.append(x)
			y_list.append(y)
			z_list.append(z)
			u_list.append(u)
			v_list.append(v)
			w_list.append(w)

			if has_node_ids:
				node_id_list.append(node_id)

		# save to class variables
		self.elem_list = elem_list
		self.x_list = x_list
		self.y_list = y_list
		self.z_list = z_list
		self.u_list = u_list
		self.v_list = v_list
		self.w_list = w_list

		if has_node_ids:
			self.node_id_list = node_id_list

	def write_vtk(self, path, name, id_flag=False, num_blade_elems=[]):
		""" write_vtk(path, name) : writes the wake element data to a time series of VTK files in a
			location specified by `path`.

			A Paraview .pvd file that contains the normalized times at each timestep is also written.
			If id_flag=True, each particle is assigned an ID so that it can be tracked. The input 
			list blade_elems should contain the number of elements for each blade.
			
			ID data is a scalar (integer)
			Velocity data is a vector."""

		from evtk.hl import pointsToVTK 	# evtk module - import only if this function is called
		import xml.etree.cElementTree as ET # xml module  -                 "

		# compute number of wake element nodes
		num_blade_nodes = np.array(num_blade_elems) + 1

		# get variables
		times  = self.times
		elem_list = self.elem_list
		x_list = self.x_list
		y_list = self.y_list
		z_list = self.z_list
		u_list = self.u_list
		v_list = self.v_list
		w_list = self.w_list

		# set the collection filename
		collection_fname = name + ".pvd"

		# set up XML tree for PVD collection file
		root = ET.Element("VTKFile")
		root.set("type", "Collection")
		collection = ET.SubElement(root, "Collection")

		# write vtk unstructured (point) file
		for ti, time in enumerate(times):
			# base name of data file
			vtk_name = name + '_' + str(ti)

			if id_flag:
				# generate the id numbers for the elements as an array with length of num_wake_elems
				# newer elements have a higher id number.
				
				# get the list of elements
				elems = elem_list[ti] 

				# compute the number of wake elements at this particular timestep
				num_wake_elems = len(u_list[ti])

				# compute the timestep number from the number of elements
				elems_per_timestep = sum(num_blade_nodes)
				nt = num_wake_elems/elems_per_timestep

				# generate the id numbers
				multiplier = np.mod((np.arange(num_wake_elems)), nt)
				adder = multiplier * elems_per_timestep
				id_nums = adder + elems

				# convert to a uint32 (fixes some problems importing into VTK/ParaView)
				id_nums = np.int32(id_nums)

				# store vector field in a dict
				data = {'velocity' : (u_list[ti],
									  v_list[ti],
									  w_list[ti]),
						'node_id' : id_nums}
			
			else:
				data = {'velocity' : (u_list[ti],
									  v_list[ti],
									  w_list[ti])}

			# write data
			data_filename = pointsToVTK(path + '/' + vtk_name, x_list[ti], y_list[ti], z_list[ti], data)

			# add elements to XML tree for PVD collection file
			dataset = ET.SubElement(collection, "DataSet")
			dataset.set("timestep", str(time))
			dataset.set("file", os.path.basename(data_filename))

		# write the collection file
		tree = ET.ElementTree(root)
		tree.write(path + '/' + collection_fname, xml_declaration=True)


#####################################
########### Wake Grid Data ##########
#####################################
class CactusWakeGrid():
	""" Class which loads WakeGridData from a pandas dataframe and creates appropriately-shaped
			Numpy arrays. Grid node locations X,Y,Z are assumed to be invariant in time. """

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

		# extract the columns
		x  = df.loc[:,x_col_name]
		y  = df.loc[:,y_col_name]
		z  = df.loc[:,z_col_name]
		u  = df.loc[:,u_col_name]
		v  = df.loc[:,v_col_name]
		w  = df.loc[:,w_col_name]

		# get grid dimensions
		self.nt = self.num_times
		self.nx = len(x.unique())
		self.ny = len(y.unique())
		self.nz = len(z.unique())

		self.dx = (x.max() - x.min())/self.nx
		self.dy = (y.max() - y.min())/self.ny
		self.dz = (z.max() - z.min())/self.nz
		
		# reshape grid nodes into 3-D numpy array
		# note that in Python, the final index is the fastest changing
		X = np.reshape(x, [self.nt, self.nz, self.ny, self.nx])
		Y = np.reshape(y, [self.nt, self.nz, self.ny, self.nx])
		Z = np.reshape(z, [self.nt, self.nz, self.ny, self.nx])

		self.X = np.reshape(X[0,:,:,:], [self.nz, self.ny, self.nx])
		self.Y = np.reshape(Y[0,:,:,:], [self.nz, self.ny, self.nx])
		self.Z = np.reshape(Z[0,:,:,:], [self.nz, self.ny, self.nx])
		
		# reshape velocity fields into 4-D numpy array
		self.U = np.reshape(u, [self.nt, self.nz, self.ny, self.nx])
		self.V = np.reshape(v, [self.nt, self.nz, self.ny, self.nx])
		self.W = np.reshape(w, [self.nt, self.nz, self.ny, self.nx])

	def slice_in_space(self, axis, i, inplace=False):
		""" Returns a subset of of the data on a plane normal to the specified axis
			at the n-th index from the minimum. Return a 3-D array with the first 
			index being in time.

			If the inplace flag is set to True, the class memory will be overwritten with the sliced data.
			This is useful for reducing memory usage of wake data, if not all the wake data is required 
			in memory."""

		def slice_reshape(axis, F, i, nt, n1, n2):
			if axis == 'x':
				F = np.reshape(F[:,:,:,i], [nt, n2, n1])
			elif axis == 'y':
				F = np.reshape(F[:,:,i,:], [nt, n2, n1])
			elif axis == 'z':
				F = np.reshape(F[:,i,:,:], [nt, n2, n1])

			return F

		# set the grid dimensions depending on which axis is selected
		nt = self.nt
		if axis=='x':
			n1 = self.ny
			n2 = self.nz
		if axis=='y':
			n1 = self.nx
			n2 = self.nz
		if axis=='z':
			n1 = self.nx
			n2 = self.ny	

		# reshape the velocity data
		U = slice_reshape(axis, self.U, i, nt, n1, n2)
		V = slice_reshape(axis, self.V, i, nt, n1, n2)
		W = slice_reshape(axis, self.W, i, nt, n1, n2)

		if inplace==True:
			self.U = U
			self.V = V
			self.W = W

		return U, V, W

	def write_vtk(self, path, name):
		""" write_vtk(path, name) : writes the wake grid data to a time series of VTK files in a
				location specified by `path`. A Paraview .pvd file that contains the normalized
				times at each timestep is also written. Velocity data is written as a vector field."""

		from evtk.hl import gridToVTK 		# evtk module - import only if this function is called
		import xml.etree.cElementTree as ET # xml module  -                 "

		# load data
		times = self.times
		X = self.X
		Y = self.Y
		Z = self.Z
		U = self.U
		V = self.V
		W = self.W

		# set the collection filename
		collection_fname = name + ".pvd"

		# set up XML tree for PVD collection file
		root = ET.Element("VTKFile")
		root.set("type", "Collection")
		collection = ET.SubElement(root, "Collection")

		# write the VTK files
		for ti, time in enumerate(times):
			# base name of data file
			vtk_name = name + '_' + str(ti)
			
			# write vector fields to data file
			velocity = (U[ti,Ellipsis], V[ti,Ellipsis], W[ti,Ellipsis])
			data_filename = gridToVTK(path + '/' + vtk_name, X, Y, Z, pointData={"velocity": velocity})

			# add elements to XML tree for PVD collection file
			dataset = ET.SubElement(collection, "DataSet")
			dataset.set("timestep", str(time))
			dataset.set("file", os.path.basename(data_filename))

		# write the collection file
		tree = ET.ElementTree(root)
		tree.write(path + '/' + collection_fname, xml_declaration=True)


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
