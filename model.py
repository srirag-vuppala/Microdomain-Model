"""
**IF THERE IS A FUNCTION YOU CAN'T FIND THE FUNCTION ITS PROBABLY PLACED IN THE UTILITIES**
"""
# Imports 
import numpy as np

from integral import *
from utilities import *
from transformMatrices import *
import os
# import hh
import matplotlib.pyplot as plt
from ionic_current_cubic import *

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



def simulate(A, B, strand, n_cells, delta_v, phi_resting, a):
    phi_now = strand

    # row control to make A full rank
    # replacement_row = np.ones(n_cells*2)
    # replacement_row[n_cells+1:] = 0
    # replacement_row_2 = np.ones(n_cells*2)
    # replacement_row_2[:n_cells+1] = 0
    # A[0] = replacement_row
    # A[-1] = replacement_row_2

    # print("Rank of A: " + str(np.linalg.matrix_rank(A)))

    # setting an arbitrary number of 1000 
    # can be tweaked if wanted to run for longer or shorter 
    for i in range(1000):

        intra, extra = split_list(phi_now) # split the strand to make sure we are modifying intra and extra properly

        # stimulus
        # Setting a potential difference by setting the extremeties(first and last points) of the strand to opposite voltages
        if i < 5:
            extra[0] += -100
            extra[-1] += 100
        
        phi_trans = intra - extra
        phi_now = join_list(intra, extra)

        if i < 200 or i % 100 == 0 or i == 999:
            display_plot(n_cells, phi_trans, i)
            print("i: " + str(i))
            print(phi_trans)

        # check handwritten doc for variable naming
        RHS_left_term = np.matmul(B, phi_now)
        # if i == 0:
            # phi_trans[-1] = 1 | setting a control to prevent blowup ?
            # RHS_left_term[n_cells-1] = 1
            # RHS_left_term[-1] = 1
            # A[-1] = 1
        #RHS_right_term = generate_ionic_current(phi_now, ...)
        #RHS_right_term = generate_ionic_dummy(phi_now*delta_v)
        #RHS_right_term = generate_cubic_integral(phi_resting, a)
        #RHS_soln_term = RHS_left_term + RHS_right_term
        soln_term = RHS_left_term
        # soln_term[0] = 1
        #soln_term[-1] = 1
        phi_next = np.linalg.solve(A, soln_term)
        #phi_next = np.linalg.lstsq(A, soln_term)

        # set up for next iteration 
        phi_now = phi_next

def main():
    # Prep work 
    # Create the strand first
    a = .5
    delta_v = 120
    phi_resting = -70/delta_v
    n_cells = 50
    delta_x = 0.1
    delta_t = 1
    J = -3.62*(10**-5) 
    cell_len = 1
    S = 4*cell_len**2 + (2*cell_len**2)/36

    strand = create_strand(n_cells, phi_resting)
    print(strand)

    #TODO : make sure the integral is proper
    integral = generate_integral(phi_resting * delta_v)
    #integral = generate_cubic_integral(phi_resting, a)

    A_84 = create_A_84(delta_x, n_cells, J, S, delta_t)
    B_84 = create_B_84(delta_x, n_cells, J, S, delta_t)
    A_83 = create_A_83(delta_x, n_cells, S, delta_t, integral, cell_len)
    B_83 = create_B_83(delta_x, n_cells, S, delta_t, integral, cell_len)

    A_comb = np.concatenate((A_83, A_84))
    A_comb = np.where(A_comb == -0,0, A_comb )
    B_comb = np.concatenate((B_83, B_84))

    # A_comb = np.identity(2*n_cells) 

    matprint(A_comb)
    matprint(B_comb)
    # Simulate
    simulate(A_comb, B_comb, strand, n_cells, delta_v, phi_resting, a)

    os.system("ffmpeg -y -i 'timestep%d.jpg' microdomain.mp4")
    # to delete the images run this command "rm -f *.jpg" in the terminal. Make sure you are in the right directory/folder.
    # os.system("rm -f *.jpg")


if __name__ == '__main__':
    main()
