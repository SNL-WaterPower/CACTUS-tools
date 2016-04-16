import numpy as np

class CactusGeom():
    """Class for parsing CACTUS geometry data.

    Attributes
    ----------
    globalvars : dict
        Global variables from geometry file.
    blades : list
        List of dicts of blade geometry variables.
    struts : list
        List of dicts of strut geometry variables.
    """

    def __init__(self, geom_filename):
        """Initializes the class, reads in geometry data to class attributes."""
        # read the geometry file into public variables
        #   globalvars, blades, and struts are dictionaries whose index values
        #   are the variable names as strings   
        self.globalvars, self.blades, self.struts = self.read_geom(geom_filename)

        # calculate the radial positions of blade elements
        #   (useful for calculating torques/moments, and for plotting against radial position on
        #   axial flow turbines)
        self.__calculate_r_elem()
        self.__calculate_dr_elem()


    #####################################
    ######### Private Functions #########
    #####################################       
    def __calculate_r_elem(self):
        """Calculates the non-dimensionalized distance from blade elements to\
        the rotation axis."""
        for blade_num, blade in enumerate(self.blades):
            # allocate space for an array
            r_temp = np.zeros(blade['NElem'])

            # get the coordinates element centers
            pex = blade['PEx']
            pey = blade['PEy']
            pez = blade['PEz']

            # get axis of rotation and rotation axis coincident point
            rot_axis       = self.globalvars['RotN']
            rot_coincident = self.globalvars['RotP']

            # reshape, loop through the points
            positions = np.transpose(np.vstack((pex,pey,pez)))
            for elem_num, element_center in enumerate(positions):
                r_temp[elem_num] = self.distance_to_rotation_axis(element_center, rot_axis, rot_coincident)

            self.blades[blade_num]['r_over_R_elem'] = r_temp

    def __calculate_dr_elem(self):
        """Calculates the non-dimensionalized dr distribution."""
        for blade_num, blade in enumerate(self.blades):
            # allocate space for an array
            dr_temp = np.zeros(blade['NElem'])

            # get the location of element node points (based on quarter chord lines)
            qcx = blade['QCx']
            qcy = blade['QCy']
            qcz = blade['QCz']

            # get axis of rotation and rotation axis coincident point
            rot_axis       = self.globalvars['RotN']
            rot_coincident = self.globalvars['RotP']

            # for each node location
            node_locations = np.transpose(np.vstack((qcx,qcy,qcz)))
            for elem_num, node_loc_pair in enumerate(zip(node_locations[:-1], node_locations[1:])):
                # get the endpoints of node and subsequent node
                node_A = node_loc_pair[0]
                node_B = node_loc_pair[1]

                # compute radial location of both
                r_over_R_A = self.distance_to_rotation_axis(node_A, rot_axis, rot_coincident)
                r_over_R_B = self.distance_to_rotation_axis(node_B, rot_axis, rot_coincident)

                # compute dr as the distance between radial location of node and next node
                dr_temp[elem_num] = r_over_R_B - r_over_R_A

            self.blades[blade_num]['dr_over_R'] = dr_temp


    ####################################
    ######### Public Functions #########
    ####################################
    def read_geom(self, filename,
                  blade_num_lines=24,
                  strut_num_lines=18,):
        """Reads data from a CACTUS .geom file.

        Parameters
        ----------
        filename : str
            The geometry filename.
        blade_num_lines : int
            Number of lines of data per blade.
        strut_num_lines : int
            Number of lines of data per strut.

        Returns
        -------
        globalvars : dict
            Global variables from geometry file.
        blade : list
            List of dicts of blade geometry variables.
        strut : list
            List of dicts of strut geometry variables.
        """



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
                if line.strip():
                    var_name, value = (line.strip()).split(':')
                else:
                    break

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

    def distance_to_rotation_axis(self, p, n, a):
        """Computes the distance between a point and a rotation axis.

        http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Vector_formulation

        Parameters
        ----------
        p : numpy.array
            The coordinates of the point.
        n : numpy.array
            Vector of the rotation axis.
        a : numpy.array
            Coincident point of the rotation axis.
        """

        return np.linalg.norm((a-p) - np.dot((a-p), n)*n)
        