import numpy as np
from numpy import ndarray

# The function inputs are a matrix A and vector b. Note that vector b is a list, but not a nested list.
# The function output is a 1D ndarray.
def gauss_eli(A, b):
    # Forming the upper triangular matrix

    # Make input matrix into a nparray object
    A = np.array(A, dtype=float)

    # Check that the determinant is greater than approximately zero
    if abs(np.linalg.det(A)) > 10 ** (-14):
        np.array(b, dtype=float)

        # Make an augmented matrix with a given b
        augment = np.column_stack((A, b))

        # Find the dimension of the matrix
        global dim
        dim = np.shape(augment)

        # Loop through the ith row of the augmented matrix depending on its dim
        for i in range(1, dim[0]):

            # Begin a loop if the ith row does not have zeros in position we want
            while str(np.array_equal(augment[i][0:i], np.zeros(i))) == 'False':

                # Loop through the ath entries of the ith row
                for a in range(0, i):
                    """If the ath entry of a row is not zero, 
                    add the ath row of the augmented matrix to the ith row 
                    we are in to cancel the ath term in the ith row"""
                    if augment[a][a] != 0 and a != i:
                        pivot = augment[a]
                        k = augment[i][a] / pivot[a]
                        augment[i] = np.copy(augment[i]) + pivot * -k
                        """If the ath row of the augmented matrix has zero in the ath position, 
                        interchange ith and ath row in the augmented matrix"""
                    else:
                        permute = np.copy(augment[i])
                        augment[i] = np.copy(augment[a])
                        augment[a] = np.copy(permute)
        # Back Substitution
        # Create a ndarray with entries as zeros
        z = np.zeros(dim[0])

        # Loop through the rows of the upper triangular matrix from bottom to top
        for i in range(dim[0] - 1, -1, -1):
            # Create sum_term to account for the terms that have to be subtracted in each sys. of equations
            sum_term = 0
            for c in range(0, dim[0] - 1 - i):
                sum_term = sum_term + (z[dim[0] - 1 - c] * augment[i][dim[0] - 1 - c])
            # Make the solution to ith row of the matrix as the ith entry of z
            z[i] = (augment[i][dim[1] - 1] - sum_term) / (augment[i][i])
        print(z)
    else:
        print("The matrix is not invertible.")