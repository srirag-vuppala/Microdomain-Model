import numpy as np


def create_A_83():
    """ The coeff matrix of the LHS in 83"""
    return 0

def create_A(delta_x, number_of_nodes):
    """ A = [L_2 0]"""
    L_2 = create_L_2(delta_x, number_of_nodes) 
    


def create_L_2(delta_x, number_of_nodes ):
    """ Create a laplacian matrix for the 3 point stencil that corresponds to equation 83"""
    constant = 1/(delta_x**2)

    L = np.diagflat([constant for _ in range(number_of_nodes-1)], -1) + np.diagflat([-2*constant for _ in range(number_of_nodes)]) + np.diagflat([constant for _ in range(number_of_nodes-1)], 1)
    # Defining the boundary conditions
    # TODO Make sure this is right
    L[0][0] = 1+constant
    L[-1][-1] = 1+constant
    return L

def create_L_3(delta_x, number_of_nodes):
    """ Create a laplacian matrix for the 5 point stencil that corresponds to equation 84"""
    constant = 1/((delta_x**3)*2)

    L = np.diagflat([-1*constant for _ in range(number_of_nodes-2)],-2) + np.diagflat([constant*2 for _ in range(number_of_nodes-1)], -1) + np.diagflat([-2*constant for _ in range(number_of_nodes-1)], 1) + np.diagflat([constant for _ in range(number_of_nodes-2)], 2)

    # Defining the boundary conditions
    # TODO Make sure this is right | WHAT IS THIS
    L[0][0] = 1+constant
    L[-1][-1] = 1+constant
    return L
