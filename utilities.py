import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


# TODO figure out to way to set precision with matprint
np.set_printoptions(precision=3)


def matprint(mat, fmt="g"):
    """
    USE: This function prints out the matrix in a nice and clean way.
         It also represents the matrix as we view it conceptually not just how its stored in numpy.
    
    If an argument likeTheory is given it'll print out the array like how we conceptualize it in our theory. 
    i.e our sheet is a array of arrays of columns
    """
    mat = np.asarray(mat)
    msg2 = "Actual data representation"
    print('-'*len(msg2))
    print(msg2)
    print('-'*len(msg2))

    col_maxes = [max([len(("{:"+fmt+"}").format(x)) for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="   ")
        print("")
        print("")
    
    print('-'*len(msg2))
        

def check_laplace_matrix(L):
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
    
def display_plot(n_cells, phi_trans, i):
    # plt.scatter(np.arange(1, 11, 1), np.matmul(create_B(n_cells), phi_now))
    plt.plot(np.arange(1, n_cells+1, 1), phi_trans)
    plt.xlabel('Node')
    plt.ylabel('Transmembrane Potential')
    plt.title('time = ' + str(i/100))
    plt.savefig('timestep'+str(i)+'.jpg')
    plt.clf()

def split_list(a_list):
    half = len(a_list)//2
    return a_list[:half], a_list[half:]

def join_list(A, B):
    return np.concatenate((A, B))



def create_block_diag_matrix(A, B):
    # block diagonal form 
    # Create the big laplace matrix with Vi's L and Ve's L
    final = []
    zero = np.zeros(len(A))
    for ele in A:
        final.append(np.concatenate((ele, zero), axis=0))
    for ele in B:
        final.append(np.concatenate((zero, ele), axis=0))
    return np.asarray(final)


def main():
    arr = np.ones(4)
    print(arr)
    print(join_list(arr, arr))


if __name__ == '__main__':
    main()