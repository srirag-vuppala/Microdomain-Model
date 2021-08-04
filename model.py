"""
**IF THERE IS A FUNCTION YOU CAN'T FIND THE FUNCTION ITS PROBABLY PLACED IN THE UTILITIES**
"""
# Imports 
import numpy as np

from integral import generate_integral
from utilities import *
from transformMatrices import *
import os
import hh

# Prints the arrays properly with elements as X.YYY
np.set_printoptions(precision=3)


def create_strand(n_cells):
    # Creating a strand representing n cells
    # The extracellular space is initially at a potential of   0 mV 
    # The strand is [phi_i phi_e]
    strand = np.zeros(n_cells * 2)
    # The intracellular space is initially at a potential of -70 mV
    strand[:n_cells] += -70
    return strand 


def generate_ionic_current(V, A, delta_t):
    V_send = np.matmul(A, V)
    V_send = np.add(V_send, 70)
    I_ion = hh.HodgkinHuxley().main(flat(V_send))
    return 1000 * delta_t * np.asarray(I_ion)


def generate_A_matrix(n_cells, delta_x):
    # Creates "A" matrix from Microdomain Coding Prep
    # "A" matrix approximates second derivative of phi intra, A is n_cells x 2*n_cells
    A = np.zeros([n_cells, n_cells])
    A += np.diagflat(np.ones(n_cells) * -2)  # add main diagonal of -2's
    A += np.diagflat(np.ones(n_cells - 1), 1)  # add super diagonal of 1's
    A += np.diagflat(np.ones(n_cells - 1), -1)  # add sub diagonal of 1's
    zero = np.zeros([n_cells, n_cells])  # create zero matrix to be appended to A
    A = np.concatenate((A, zero), axis=1)  # horizontally merge A and zero
    A *= (1 / delta_x ** 2)  # multiply A by 1/delta_x^2
    return A


def generate_B_matrix(n_cells):
    # Creates "B" matrix from Microdomain Coding Prep
    # "B" matrix transforms phi_n into transmembrane potential, B is n_cells x 2*n_cells
    B = np.eye(n_cells)
    B_merge = np.zeros([n_cells, n_cells])
    B_merge += np.diagflat(np.ones(n_cells) * -1)  # created -I component of B this way to avoid -0's
    B = np.concatenate((B, B_merge), axis=1)  # horizontally merge B and B_merge
    return B


def generate_E_matrix(n_cells):
    # creates "E" matrix from Microdomain coding prep
    # "E" matrix transforms phi_n into phi_i - phi_neighbors, E is n_cells x 2*n_cells
    main_diag = np.ones(n_cells) * 2  # create main diagonal of 2's
    main_diag[0] = 1  # change first entry to 1
    main_diag[-1] = 1  #  change last entry to 1
    # first and last entries changed to 1 to account for boundary cases
    E = np.zeros([n_cells, n_cells])
    E += np.diagflat(main_diag)  # add main_diag to diagonal of E
    E += np.diagflat(np.ones(n_cells - 1) * -1, 1)  # create super diagonal of -1's
    E += np.diagflat(np.ones(n_cells - 1) * -1, -1)  # create sub diagonal of -1's
    zero = np.zeros([n_cells, n_cells])
    E = np.concatenate((E, zero), axis=1)  # horizontally concatenate E and zero matrix
    return E

def generate_Z_matrix(n_cells, delta_x):
    # Creates "Z" matrix from Microdomain Coding Prep
    # "Z" matrix approximates third derivative of phi extra, Z is n_cells x 2*n_cells
    Z = np.zeros([n_cells, n_cells])
    Z += np.diagflat(np.ones(n_cells - 1) * -2, 1)  # create super diagonal of -2's
    Z += np.diagflat(np.ones(n_cells - 1) * 2, -1)  # create sub diagonal of 2's
    Z += np.diagflat(np.ones(n_cells - 2), 2)  # create diagonal above super diagonal of 1's
    Z += np.diagflat(np.ones(n_cells - 2) * -1, -2)  # create diagonal below sub diagonal of -1's
    zero = np.zeros([n_cells, n_cells])
    Z = np.concatenate((zero, Z), axis=1)  # horizontally concatenate zero matrix and Z
    Z *= 1/(2*delta_x**3)
    return Z

def simulate(A, B, strand, S ):
    phi_now = strand
    for i in range(1000):
        # stimulus
        if i < 5:
            phi_now[0] = -40
            phi_now[-1] = -40
        left_term = np.matmul(B, phi_now)
        right_term = generate_ionic_current(phi_now, ...)/S
        soln_term = left_term + right_term
        phi_next = np.linalg.solve(A, soln_term)

        # set up for next iteration 
        phi_now = phi_next





def main():
    # Prep work 
    # Create the strand first
    n_cells = 10
    strand = create_strand(n_cells)
    delta_x = 0.01
    delta_t = 0.01
    J = 0.1
    S = 0.1
    v_resting = -70

    #TODO : make sure the integral is proper
    integral = generate_integral(v_resting)

    A_84 = create_A_84(delta_x, n_cells, J, S, delta_t)
    B_84 = create_B_84(delta_x, n_cells, J, S, delta_t)
    A_83 = create_A_83(delta_x, n_cells, S, delta_t, integral)
    B_83 = create_B_83(delta_x, n_cells, S, delta_t, integral)

    A_comb = np.concatenate(A_83, A_84)
    B_comb = np.concatenate(B_83, B_84)

    combined_strand = np.concatenate(strand, strand)

    # Simulate
    simulate(A_comb, B_comb, combined_strand, S)

    # os.system("ffmpeg -y -i 'foo%03d.jpg' bidomain.mp4")
    # os.system("rm -f *.jpg")


if __name__ == '__main__':
    main()
