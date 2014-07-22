#!/usr/bin/python
# pyCactusGeom.py : A module for parsing CACTUS geometry files.


import numpy as np

class CactusGeom():
	def __init__(self, geom_filename):

		# read the geometry file into public variables
		# 	globalvars, blades, and struts are dictionaries whose index values
		# 	are the variable names as strings	
		self.globalvars, self.blades, self.struts = self.__read_geom(geom_filename)

		# calculate the radial positions of blade elements
		# 	(useful for calculating torques/moments, and for plotting against radial position on
		# 	axial flow turbines)
		self.__calculate_r_elem()


	#####################################
	######### Private Functions #########
	#####################################		
	def __calculate_r_elem(self):
		""" Calculates the non-dimensionalized distance from blade elements to the rotation axis. """
		for blade in self.blades:
			# allocate space for an array
			blade['r_elem'] = np.zeros(blade['NElem'])

			# get the coordinates element centers
			pex = blade['PEx']
			pey = blade['PEy']
			pez = blade['PEz']

			# get axis of rotation and rotation axis coincident point
			rot_axis       = self.globalvars['RotN']
			rot_coincident = self.globalvars['RotP']

			# reshape, loop through the points
			for elem_num, element_center in enumerate(np.transpose(np.vstack((pex,pey,pez)))):
				(blade['r_elem'])[elem_num] = self.distance_to_rotation_axis(element_center, rot_axis, rot_coincident)


	def __read_geom(self, filename):
		""" Reads the .geom file from a CACTUS run."""

		# number of data lines for blade and strut, excluding header
		blade_num_lines = 24
		strut_num_lines = 18

		def read_block(lines):
			int_vars = ['NElem', 'FlipN']

			block_vars = {}

			for line in lines:
				# split the line at the colon
				var_name, value = (line.strip()).split(':')

				# if the variable is an integer, cast it as such
				if var_name in int_vars:
					values = map(int, (value.strip()).split())
					block_vars[var_name] = values

				# otherwise, store it as a float np.array
				else:
					values = np.array((value.strip()).split(), dtype=np.float)
					block_vars[var_name] = values

			return block_vars

		int_vars    = ['iSect', 'NBlade', 'NStrut']
		string_vars = ['Type']

		# open the file for reading
		f = open(filename, 'r')

		# initialize dict for global variables
		global_vars = {}

		# initialize empty list for blade and strut variables
		blade = []
		strut = []

		# initialize lists for starting line numbers of blades and struts
		blade_line_nums = []
		strut_line_nums = []

		# read file line by line
		with open(filename, 'r') as fid:
			lines = fid.readlines()
			for line_num, line in enumerate(lines):

				# split the line at the colon
				var_name, value = (line.strip()).split(':')

				# if the variable name is Blade or Strut, take note of the line number
				if var_name.startswith('Blade'):
					blade_line_nums.append(line_num)

				if var_name.startswith('Strut'):
					strut_line_nums.append(line_num)

			# starting at the beginning, read in the global variables
			for line in lines[:min(blade_line_nums + strut_line_nums)]:
				var_name, value = (line.strip()).split(':')

				# if the variable is an integer or string, cast it as such
				if var_name in int_vars:
					values = map(int, (value.strip()).split())
					global_vars[var_name] = values

				elif var_name in string_vars:
					values = str(value)
					global_vars[var_name] = values

				# otherwise, cast it as a float np.array
				else:
					values = np.array((value.strip()).split(), dtype=np.float)
					global_vars[var_name] = values

			# read in the blocks
			for line_start in blade_line_nums:
				blade.append(read_block(lines[line_start:line_start + blade_num_lines + 1]))

			for line_start, line_end in zip(strut_line_nums, strut_line_nums[1:]):
				strut.append(read_block(lines[line_start:line_start + strut_num_lines + 1]))
		
		return global_vars, blade, strut


	####################################
	######### Public Functions #########
	####################################
	def distance_to_rotation_axis(self, p, n, a):
		""" Computes the distance of a point p from the rotation axis with direction specified by unit
			vector rot_n with a coincident point given by a.

			http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Vector_formulation """

		return np.linalg.norm((a-p) - np.dot((a-p), n)*n)