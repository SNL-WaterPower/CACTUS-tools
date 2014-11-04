def slice_in_space(self, axis, i):
    """ Returns a subset of of the data on a plane normal to the specified axis
    at the i-th index from the minimum. Return a 2-D array.
    """

    def slice_reshape(axis, F, i, n1, n2):
        if axis == 'x':
            F = np.reshape(F[:,:,i], [n2, n1])
        elif axis == 'y':
            F = np.reshape(F[:,i,:], [n2, n1])
        elif axis == 'z':
            F = np.reshape(F[i,:,:], [n2, n1])

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

    return U, V, W

def find_nearest(array, value):
    """ find_nearest() : returns the value and index of the nearest point in a 1-d array. """
    index = (np.abs(array-value)).argmin()
    nearest_value = array[index]

    return nearest_value, index


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