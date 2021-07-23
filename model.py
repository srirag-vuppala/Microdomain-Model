"""
**IF THERE IS A FUNCTION YOU CAN'T FIND THE FUNCTION ITS PROBABLY PLACED IN THE UTILITIES**
"""
# Imports 
import numpy as np
from utilities import *
import os
import hh

# Prints the arrays properly with elements as X.YYY
np.set_printoptions(precision=3)


def create_strand(n_cells):
    # Creating a strand representing n cells
    # The extracellular space is initially at a potential of   0 mV 
    # The strand is [phi_i phi_e]
    strand = np.zeros(n_cells*2)
    # The intracellular space is initially at a potential of -70 mV
    strand[:n_cells] += -70
    return strand 
     
def create_laplace_matrix(V, c):
    L = []
    return L

def generate_ionic_current(V, A, delta_t):
    V_send = np.matmul(A, V)
    V_send = np.add(V_send, 70)
    I_ion = hh.HodgkinHuxley().main(flat(V_send))
    return 1000*delta_t*np.asarray(I_ion) 


def simulate(strand, L):
    i = 0
 

def main():
    # Create the strand first
    n_cells = 10
    strand = create_strand(n_cells)


    # create our laplacian matrices

    # Merge the Laplacian matrices 
    # We do this to to unify all of our variables to make it easy to keep track off of 

    # Simulate

    os.system("ffmpeg -y -i 'foo%03d.jpg' bidomain.mp4")
    os.system("rm -f *.jpg")

if __name__ == '__main__':
    main()