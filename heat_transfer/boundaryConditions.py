import sys
sys.path.insert(0, '../../')
from tools.tools import *
parameters_directory="../../heat_transfer/parameters.txt"
L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points, num_of_timeSteps_for_plotting = readParameters(parameters_directory)
# setting up boundary conditions

def BC_1D_Dirichlet(T, x, mask):
	Tbl = 500 # temperature at the left boundary
	Tbr = 0 # temperature at the right boundary
	for i in range(len(mask)):
		mask_i = mask[i]
		
		if mask_i==0: # at the ghost points
			if x[i] < 0:
				T[i] = Tbl
			elif x[i] > L:
				T[i] = Tbr
			else:
				print("Error: shouldnt be here.")
		else: # within the domain
			mask_im1 = mask[i-1]
			mask_ip1 = mask[i+1]
			if mask_im1==0: # left boundary of the domain
				T[i] = Tbl
			if mask_ip1==0: # right boundary of the domain
				T[i] = Tbr
	return T


def BC_1D_Mixed(T, x, _lambda, mask):
	Tbl = 0 # temperature at the left boundary
	heat_flux_density = 0 # because 1D, heat flux density (W/m2) = dT/dx * k = (T[i]-T[i-1])/dx * k[i], where k is thermal conductivity
	for i in range(len(mask)):
		mask_i = mask[i]
		
		if mask_i==0: # at the ghost points
			if x[i] < 0:
				T[i] = Tbl
			elif x[i] > L:
				continue
				#T[i] = T[i-1] # heat_flux_density/_lambda[2][i] * dx + T[i-1]
			else:
				print("Error: shouldnt be here.")
		else: # within the domain
			mask_im1 = mask[i-1]
			mask_ip1 = mask[i+1]
			if mask_im1==0: # left boundary of the domain
				T[i] = Tbl

	return T



def BC_1D_Dirichlet_Tbl1_Tbr0(T, x, mask):
	Tbl = 1  # temperature at the left boundary
	Tbr = 0  # temperature at the right boundary
	for i in range(len(mask)):
		mask_i = mask[i]

		if mask_i == 0:  # at the ghost points
			if x[i] < 0:
				T[i] = Tbl
			elif x[i] > L:
				T[i] = Tbr
			else:
				print("Error: shouldnt be here.")
		else:  # within the domain
			mask_im1 = mask[i - 1]
			mask_ip1 = mask[i + 1]
			if mask_im1 == 0:  # left boundary of the domain
				T[i] = Tbl
			if mask_ip1 == 0:  # right boundary of the domain
				T[i] = Tbr
	return T


def BC_1D_Dirichlet_Tbl500_Convection(T, x, _lambda, mask): # _lambda= [ [rho1, rho2,...,rhon], [cp1, cp2, ..., cpn], [k1, k2, ..., kn] ....[lambda1, lambda2, ..., lambdan]]    ]
	Tbl = 500  # temperature at the left boundary
	Tbr = 150  # temperature at the right ghost point as the reference/free-stream temperature
	# heat being convected away at the right boundary, according to the heat transfer coefficient h
	# -K * (dT/dx) = qx = h * dT, h = -k/dx, k=-h*dx
	h_br = 20
	
	for i in range(len(mask)):
		mask_i = mask[i]

		if mask_i == 0:  # at the ghost points
			if x[i] < 0:
				T[i] = Tbl
			elif x[i] > L:
				T[i] = Tbr # assigning free stream temperature here, e.g. temperature of air
			else:
				print("Error: shouldnt be here.")
				break
		else:  # within the domain
			mask_im1 = mask[i - 1]
			mask_ip1 = mask[i + 1]
			if mask_im1 == 0:  # left boundary of the domain
				T[i] = Tbl
			if mask_ip1 == 0:  # right boundary of the domain
				
				rho_i = _lambda[0][i] 
				cp_i = _lambda[1][i] 
				dx = x[i+1] - x[i]
				k_br = h_br*dx
				_lambda[2][i] = k_br
				_lambda[3][i] = k_br/(cp_i*rho_i)
	return T, _lambda


def BC_2D_Dirichlet(T, x, mask):
	Tbl = 500 # temperature at the left boundary
	Tbr = 0 # temperature at the right boundary




    # TODO!

	for i in range(len(mask)):
		mask_i = mask[i]
		
		if mask_i==0: # at the ghost points
			if x[i] < 0:
				T[i] = Tbl
			elif x[i] > L:
				T[i] = Tbr
			else:
				print("Error: shouldnt be here.")
		else: # within the domain
			mask_im1 = mask[i-1]
			mask_ip1 = mask[i+1]
			if mask_im1==0: # left boundary of the domain
				T[i] = Tbl
			if mask_ip1==0: # right boundary of the domain
				T[i] = Tbr
	return T