import sys
sys.path.insert(0, '../tools/')
from tools import *
L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points = readParameters()
# setting up boundary conditions

def BC_1D_Dirichlet(T, x, mask):
	Tbl = 0 # temperature at the left boundary
	Tbr = 0 # temperature at the right boundary
	for i in range(len(mask)):
		mask_i = mask[i]
		mask_im1 = mask[i-1]
		mask_ip1 = mask[i+1]
		if mask_i==0: # at the ghost points
			if x[i] < 0:
				T[i] = Tbl
			elif x[i] > L:
				T[i] = Tbr
			else:
				print("Error: shouldnt be here.")
		else: # within the domain
			if mask_im1==0: # left boundary of the domain
				T[i] = Tbl
			if mask_ip1==0: # right boundary of the domain
				T[i] = Tbr
	return T
