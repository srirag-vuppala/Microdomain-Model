import numpy as np 
import matplotlib.pyplot as plt 

# Define global variables
J = 3.62*(10 ** (-8))
tstep = .00025
cLen = 1.5*(10 ** (-2))
cWid = 1.3*(10 ** (-3))
intraConduct = (1.0/18)*(cLen ** 3)
sArea = 2*(cWid ** 2) + 4*(cWid*cLen) 
cVol = cLen*(cWid ** 2)

# defines derivative matrix, inputs vector length and 
# step size for the derivative

def derivMatrix(n, l):
	mat = np.zeros((n,n))
	for i in range(n-1):
		mat[i, i] = -1
		mat[i, i+1] = 1
	mat[n-1, n-2] = 1
	mat[n-1, n-1] = -1
	mat = (1.0/l)*mat
	return mat

# defines the membrane operator, inputs the number of 
# total nodes, should be even since same number of int/ext
# nodes. Returns phi_i - phi_e 

def membraneOp(n):
	k = int(n/2.0)
	mat = np.zeros((k, n))
	for i in range(k):
		mat[i, i] = 1
		mat[i, i+k] = -1
	return mat

#defines projection onto first n coordinates. Allows for 
#us to apply operations solely to phi_i 

def intProj(n):
	k = int(n/2)
	mat = np.zeros((k, n))
	for i in range(k):
		mat[i,i] = 1
	return mat

#defines projection for second n coordinates, similar to above

def extProj(n):
	k = int(n/2)
	mat = np.zeros((k, n))
	for i in range(k):
		mat[i, i + k] = 1
	return mat

#defines gap junctional operator. Inputs constant for surface 
#average of g_{in}(y).

def gapJuncOp(n, g):
	k = int(n/2)
	mat = np.zeros((k, n))
	mat[0,0] = 1
	mat[0,1] = -1
	mat[k-1, k-2] = -1
	mat[k-1, k-1] = 1
	for i in range(1, k-1):
		mat[i, i-1] = -1
		mat[i, i] = 2
		mat[i, i+1] = -1
	mat = (g/sArea)*mat
	return mat

#Creates a matrix for taking the second derivative, l is 
#delta_x

def secDerivMat(n, l):
	mat = np.zeros((n,n))
	mat[0,0] = -1
	mat[0,1] = 1
	mat[n-1, n-2] = 1
	mat[n-1, n-1] = -1
	for i in range(1, n-1):
		mat[i, i-1] = -1
		mat[i, i] = 2
		mat[i, i+1] = -1
	mat = ((1.0/l)**2)*mat
	return mat



#inputs a float and a number a that determines the root of the cubic 

def pointIon(x, a):
	return x*(a-x)*(x-1)

#to compute how the nodes change with a timestep we consider the process 
#as three steps. Given nodes we transform them with a matrix, then add
#the ionic current, then transform them with an inverse matrix.
#this function defines the matrix used in the first node transformation. 

def firstStep(nodes, xstep, gap):
	n = len(nodes)
	k = int(n/2)
	timeStepMat = (1/tstep)*membraneOp(n)
	intracellOp = (0.5)*(1.0/sArea)*(intraConduct)*(secDerivMat(k, xstep).dot(intProj(n)))
	extracellOp = (0.5)*(1.0/sArea)*(J)*(secDerivMat(k, xstep).dot(derivMatrix(k, xstep))).dot(extProj(n))
	upper = timeStepMat + intracellOp - (0.5)*gapJuncOp(n, gap)
	lower = timeStepMat + extracellOp 
	mat = np.concatenate((upper, lower))
	return mat


#inverseStep creates the matrix that we need to invert. We however will not invert the 
#matrix due to stability issues.

def inverseStep(nodes, xstep):
	n = len(nodes)
	k = int(n/2)
	intracellOp = (0.5)*(1/sArea)*(intraConduct)*(secDerivMat(k, xstep).dot(intProj(n)))
	extracellOp = (0.5)*(1/sArea)*(J)*(secDerivMat(k, xstep).dot(derivMatrix(k,xstep))).dot(extProj(n))
	upperMat = (1/tstep)*membraneOp(n) - intracellOp + (0.5)*gapJuncOp(n, 259)
	lowerMat = (1/tstep)*membraneOp(n) -  extracellOp
	operator = np.concatenate((upperMat,lowerMat))
	if np.linalg.matrix_rank(operator) >= n-1:
		for i in range(n):
			operator[n-1, i] = 1
		return operator
	else:
		print('Operator of Insufficient Rank')

#defines the nodes to be .05, -.05 at either end. #TODO : Verify why it says .05  
#buids the diffusive piece. Defines inv as the 
#matrix to be solved. In the for-loop we iterate over 20 time steps
#solves equation inv(nodes) = firstStep(nodes) + ionic(nodes)
def main():

	length = 100
	nodes = -.08*np.ones(length) # array of 1s with length 100
	# nodes = -.08*nodes # array of -0.08
	nodes[0] = .2
	nodes[int(length/2)] = -.2
	totalCurr = np.sum(nodes)
	diffusivePiece = firstStep(nodes, cLen, 259)
	inv = inverseStep(nodes, 1)
	vecIon = np.vectorize(pointIon)
	print(nodes)
	plt.plot(nodes)
	for _ in range(20):
		ionic = vecIon(nodes, .6)
		ionic[length-1] = 0
		step = diffusivePiece.dot(nodes) 
		#step[length-1] = totalCurr
		nodes = np.linalg.solve(inv, step)
		print(nodes)
		print(np.sum(nodes))
		plt.plot(nodes)
	plt.show()


if __name__ == '__main__':
	main()