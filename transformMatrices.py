import numpy as np


def create_A_83(delta_x, n_nodes, S, delta_t):
    """ The coeff matrix of the LHS in 83"""
    A = create_A(delta_x, n_nodes)
    lhs = A/(2*S)
    B = create_B(n_nodes)
    rhs = B/delta_t
    return lhs - rhs

def create_B_83(delta_x, n_nodes, S, delta_t , integral):
    """ The coeff matrix of the RHS in 83"""
    B = create_B(n_nodes)
    term_1 = B/delta_t
    A = create_A(delta_x, n_nodes)
    term_2 = A/(2*S)
    C = create_C(n_nodes)
    term_3 = (integral*C)/S
    return term_3 - term_1 - term_2

def create_Z(delta_x, n_nodes):
    L_3 = create_L_3(delta_x, n_nodes)
    zero = np.zeros([n_nodes, n_nodes])
    '''matrix should look like [ 0 L3], pretty sure its incorrect in the coding prep since it should be third 
    derivative of phi_e '''
    return np.concatenate((zero, L_3), axis=1)

def create_C(n_nodes):
    E = create_E(n_nodes)
    zero = np.zeros([n_nodes, n_nodes])
    return np.concatenate((E, zero), axis=1)

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
    '''second identity matrix should be negative'''
    return np.concatenate((np.identity(n_nodes), np.negative(np.identity(n_nodes))), axis=1)
    

def create_L_2(delta_x, n_nodes ):
    """ Create a laplacian matrix for the 3 point stencil that corresponds to equation 83"""
    constant = 1/(delta_x**2)

    L = np.diagflat([constant for _ in range(n_nodes-1)], -1) + np.diagflat([-2*constant for _ in range(n_nodes)]) + np.diagflat([constant for _ in range(n_nodes-1)], 1)
    # Defining the boundary conditions
    # TODO Make sure this is right
    L[0][0] = 1+constant
    L[-1][-1] = 1+constant
    return L

def create_L_3(delta_x, n_nodes):
    """ Create a laplacian matrix for the 5 point stencil that corresponds to equation 84"""
    constant = 1/((delta_x**3)*2)

    L = np.diagflat([-1*constant for _ in range(n_nodes-2)],-2) + np.diagflat([constant*2 for _ in range(n_nodes-1)], -1) + np.diagflat([-2*constant for _ in range(n_nodes-1)], 1) + np.diagflat([constant for _ in range(n_nodes-2)], 2)

    # Defining the boundary conditions
    # TODO Make sure this is right | WHAT IS THIS
    L[0][0] = 1+constant
    L[-1][-1] = 1+constant
    return L

def main():
    print(create_B(5))
    print(create_A(5, 5))

if __name__ == '__main__':
    main()