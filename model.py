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
import matplotlib.pyplot as plt

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

def generate_ionic_dummy(v):
    return np.asarray(hh.HodgkinHuxley().main(v))


def simulate(A, B, strand, n_cells):
    phi_now = strand
    transmembrane_transform = create_B(n_cells)
     

    for i in range(1000):
        # stimulus
        if i < 5:
            # phi_now[0] = -70
            # phi_now[n_cells-1] = -70
            phi_now[n_cells] = -40
            phi_now[-1] = 40
        if i < 10:
            plt.figure(i%100 + 1)
            # plt.scatter(np.arange(1, 11, 1), np.matmul(create_B(n_cells), phi_now))
            phi_trans = np.matmul(transmembrane_transform, phi_now)
            plt.plot(np.arange(1, n_cells+1, 1), phi_trans)
            plt.xlabel('Node')
            plt.ylabel('Transmembrane Potential')
            plt.title('time = ' + str(i/100))
            plt.savefig('timestep'+str(i)+'.jpg')
            plt.clf()
            print(phi_trans)
        RHS_left_term = np.matmul(B, phi_now)
        if i == 0:
            # phi_trans[-1] = 1 | setting a control to prevent blowup ?
            # left_term[n_cells-1] = 1
            RHS_left_term[-1] = 1
            A[-1] = 1
        #right_term = generate_ionic_current(phi_now, ...)
        #right_term = generate_ionic_dummy(phi_now)
        #soln_term = left_term + right_term
        soln_term = RHS_left_term
        phi_next = np.linalg.solve(A, soln_term)

        # set up for next iteration 
        phi_now = phi_next

def main():
    # Prep work 
    # Create the strand first
    n_cells = 10
    delta_x = 0.01
    delta_t = 0.01
    J = 100 # TODO find this 
    S = 0.1
    phi_resting = -70
    strand = create_strand(n_cells, phi_resting)
    print(strand)

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

    os.system("ffmpeg -y -i 'foo%03d.jpg' microdomain.mp4")
    # os.system("rm -f *.jpg")


if __name__ == '__main__':
    main()
