import numpy as np 

import matplotlib.pyplot as plt 

#Define the global variables J_int, J_ext, sigma_int

J_int = 1

J_ext = 1

sigma_int = np.ones((3,3))





#accepts an m x n x k matrix, a coordinate for the derivative, and a step size. 

#xstep moves to the next x value. 

def x_partial(vec, x, step):

	xstep = x

	xstep[0] += 1

	deriv = (vec[xstep] - vec[x])/step

	return deriv



def y_partial(vec, x, step):

	ystep = x

	ystep[1] += 1

	deriv = (vec[ystep] - vec[x])/step

	return deriv



def z_partial(vec, x, step):

	zstep = x

	zstep[2] += 1

	deriv = (vec[zstep] - vec[x])/step

	return deriv



def grad(vec, x, step):

	gradient = [x_partial(vec, x, step), y_partial(vec, x, step), z_partial(vec, x, step)]

	return gradient



#takes the second derivative of function by taking the derivative of the 

#derivative. Centers at x. 



def x_2partial(vec, x, step):

	bstep = x

	bstep[0] -= 1

	2deriv = (x_partial(vec, x, step) - x_partial(vec, bstep, step))/step

	return 2deriv



def y_2partial(vec, x, step):

	bstep = x

	bstep[1] -= 1

	2deriv = (x_partial(vec, x, step) - x_partial(vec, bstep, step))/step

	return 2deriv



def z_2partial(vec, x, step):

	bstep = x

	bstep[2] -= 1

	2deriv = (x_partial(vec, x, step) - x_partial(vec, bstep, step))/step

	return 2deriv



def laplacian_2d(vec, x, step):

	return x_2partial(vec, x, step) + y_2partial(vec, x, step)



def laplacian_3d(vec, x, step):

	return x_2partial(vec, x, step) + y_2partial(vec, x, step) + z_2(vec, x, step)





#The Model



#Code up the ionic model. 




