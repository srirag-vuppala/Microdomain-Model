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


def create_strand(n_cells, phi_resting):
    # Creating a strand representing n cells
    # The extracellular space is initially at a potential of   0 mV 
    # The strand is [phi_i phi_e]
    strand = np.zeros(n_cells * 2)
    # The intracellular space is initially at a potential of -70 mV
    strand[:n_cells] += phi_resting 
    return strand 


def generate_ionic_current(V, A, delta_t):
    V_send = np.matmul(A, V)
    V_send = np.add(V_send, 70)
    I_ion = hh.HodgkinHuxley().main(flat(V_send))
    return 1000 * delta_t * np.asarray(I_ion)


def simulate(A, B, strand, n_cells):
    phi_now = strand
    identity_B = create_B(n_cells)
    for i in range(1000):
        # stimulus
        if i < 5:
            phi_now[0] = -40
            phi_now[-1] = -40
        left_term = np.matmul(B, phi_now)
        right_term = generate_ionic_current(phi_now, ...)
        soln_term = left_term + right_term
        phi_next = np.linalg.solve(A, soln_term)

        # set up for next iteration 
        phi_now = phi_next

def main():
    # Prep work 
    # Create the strand first
    n_cells = 10
    delta_x = 0.01
    delta_t = 0.01
    J = 0.1
    S = 0.1
    phi_resting = -70
    strand = create_strand(n_cells, phi_resting)

    #TODO : make sure the integral is proper
    integral = generate_integral(phi_resting)

    A_84 = create_A_84(delta_x, n_cells, J, S, delta_t)
    B_84 = create_B_84(delta_x, n_cells, J, S, delta_t)
    A_83 = create_A_83(delta_x, n_cells, S, delta_t, integral)
    B_83 = create_B_83(delta_x, n_cells, S, delta_t, integral)

    A_comb = np.concatenate((A_83, A_84))
    B_comb = np.concatenate((B_83, B_84))

    # Simulate
    simulate(A_comb, B_comb, strand, n_cells)

    # os.system("ffmpeg -y -i 'foo%03d.jpg' bidomain.mp4")
    # os.system("rm -f *.jpg")


if __name__ == '__main__':
    main()
