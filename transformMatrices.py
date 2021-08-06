import numpy as np
from utilities import *


def create_A_84(delta_x, n_nodes, J, S, delta_t):
    Z = create_Z(delta_x, n_nodes)
    B = create_B(n_nodes)
    return (J*Z) - (S*B/delta_t) 

def create_B_84(delta_x, n_nodes, J, S, delta_t):
    Z = create_Z(delta_x, n_nodes)
    B = create_B(n_nodes)
    return (-S*B/delta_t) - (J*Z)

def create_A_83(delta_x, n_nodes, S, delta_t, integral):
    """ The coeff matrix of the LHS in 83"""
    term_1 = create_A(delta_x, n_nodes)/2
    term_2 = S*create_B(n_nodes)/delta_t
    term_3 = create_G(n_nodes, integral)/2
    return term_1 - term_2 - term_3

def create_B_83(delta_x, n_nodes, S, delta_t , integral):
    """ The coeff matrix of the RHS in 83"""
    term_1 = create_A(delta_x, n_nodes)/2
    term_2 = S*create_B(n_nodes)/delta_t
    term_3 = create_G(n_nodes, integral)/2
    return term_3 - term_1 - term_2

def create_Z(delta_x, n_nodes):
    L_3 = create_L_3(delta_x, n_nodes)
    zero = np.zeros([n_nodes, n_nodes])
    return np.concatenate((zero, L_3), axis=1)

def create_G(n_nodes, integral):
    E = create_E(n_nodes)
    zero = np.zeros([n_nodes, n_nodes])
    return integral*np.concatenate((E, zero), axis=1)

def create_E(n_nodes):
    E = np.diagflat([2 for _ in range(n_nodes)], 0) + np.diagflat([-1 for _ in range(n_nodes-1)], -1) + np.diagflat([-1 for _ in range(n_nodes-1)], 1)
    E[0][0] = 1
    E[-1][-1] = 1
    return E


def create_A(delta_x, n_nodes):
    """ A = [L_2 0]"""
    L_2 = create_L_2(delta_x, n_nodes) 
    zero = np.zeros([n_nodes, n_nodes]) 
    return np.concatenate((L_2, zero), axis=1) 

def create_B(n_nodes):
    """ B = [I -I] """
    return np.concatenate((np.identity(n_nodes), np.negative(np.identity(n_nodes))), axis=1)
    

def create_L_2(delta_x, n_nodes ):
    """ Create a laplacian matrix for the 3 point stencil that corresponds to equation 83"""
    constant = 1/(delta_x**2)

    L = np.diagflat([constant for _ in range(n_nodes-1)], -1) + np.diagflat([-2*constant for _ in range(n_nodes)], 0) + np.diagflat([constant for _ in range(n_nodes-1)], 1)
    # Defining the boundary conditions
    # TODO Make sure this is right
    L[0][0] = -constant
    L[-1][-1] = -constant
    return L

def create_L_3(delta_x, n_nodes):
    """ Create a laplacian matrix for the 5 point stencil that corresponds to equation 84"""
    constant = 1/((delta_x**3)*2)

    L = np.diagflat([-1*constant for _ in range(n_nodes-2)],-2) + np.diagflat([constant*2 for _ in range(n_nodes-1)], -1) + np.diagflat([-2*constant for _ in range(n_nodes-1)], 1) + np.diagflat([constant for _ in range(n_nodes-2)], 2)

    # Defining the boundary conditions
    # TODO Make sure this is right 
    L[0][0] = constant
    L[1][1] = constant
    L[n_nodes-2][n_nodes-2] = constant
    L[-1][-1] = constant
    return L

def check_laplace_matrix(L):
    # I might be checking this wrong
    # Each row's elements sum should be = 0
    # Each column's elements sum should be = 0
    # Thus total sum of all elements should be 0 too 
    tsum = 0
    for arr in L:
        for i in arr:
            tsum += i
    if tsum == 0:
        print("yay it works")    
    else:
        print("Try again")

def main():
    L2 = create_L_2(0.01, 10)
    L3 = create_L_3(0.01, 10)
    matprint(L2)
    matprint(L3)
    
if __name__ == '__main__':
    main()