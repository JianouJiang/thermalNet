import numpy as np
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

# Solve 2D Transient Heat Conduction Problem with Convection BCs using FTCS Finite Difference Method: https://www.youtube.com/watch?v=Ip47nsJOQqs
# make the energy equation balance, as in how much energy in should be equal to how much out...in a sense
'''                                      |     |----->y-axis, j
(0,0)------------------------------------*-----*(i-1,j)
   |                                  dx |  ---|             h=20
   |------------------------------(i,j-1)*--+--*(i,j)	     Tinf = 150
   |                                  dx |  L--| 
   --------------------------------------*-----*(i+1,j)
   v                                     |  dy |
   i 
  x-axis                            '''
# consider a finite volume with the centre at T[i][j], then this volume is a rectangle with height dx and width dy/2
# then the net heat flux from the left hand side dx, the upper wall dy/2, the lower wall dy/2, and the right hand side dx
# should be equal to
# the heat flux change over the period of dt of this rectangle
# k*A*dT/dy + k*A'*dT/dx + k*A'*dT/dx + h*A*(Tinf-T[i][j]) = rho*Cp*A'*A*dT/dt
# discretisation:
# k*dx*(T[i][j]-T[i][j-1])/dy + k*dy/2*(T[i][j]-T[i-1][j])/dx + k*dy/2*(T[i+1][j]-T[i][j])/dx + h*dx*(Tinf-T[i][j])=rho*Cp*dy/2*dx*(Tnew[i][j]-T[i][j])/dt
# obviously, we are doing a 1D situation, so we can assume  that the net heat flux in the x-direction is zero
# k*dT/dy + h*(Tinf-T[i][j]) = rho*Cp*A'*dT/dt
def BC_1D_Dirichlet_Tbl500_Convection(T, x, _lambda, mask): # _lambda= [ [rho1, rho2,...,rhon], [cp1, cp2, ..., cpn], [k1, k2, ..., kn] ....[lambda1, lambda2, ..., lambdan]]    ]
	Tbl = 500  # temperature at the left boundary
	Tbr = 150  # temperature at the right ghost point as the reference/free-stream temperature
	# heat being convected away at the right boundary, according to the heat transfer coefficient h
	# -K * (dT/dx) = qx = h * dT, h = -k/dx, k=-h*dx
	h_br = 20
	T_bc =  np.zeros(len(T))
	for i in range(len(mask)):
		mask_i = mask[i]

		if mask_i == 0:  # at the ghost points
			if x[i] < 0:
				T_bc[i] = Tbl
			elif x[i] > L:
				T_bc[i] = Tbr # assigning free stream temperature here, e.g. temperature of air
			else:
				print("Error: shouldnt be here.")
				break
		else:  # within the domain
			mask_im1 = mask[i - 1]
			mask_ip1 = mask[i + 1]
			if mask_im1 == 0:  # left boundary of the domain
				T_bc[i] = Tbl
			if mask_ip1 == 0:  # right boundary of the domain
				
				rho_i = _lambda[0][i] 
				cp_i = _lambda[1][i] 
				k_i = _lambda[2][i] 
				dx = x[i+1] - x[i]

				a_i = 1+ ( -2*k_i*dt/(rho_i*cp_i*dx*dx) - 2*h_br*dt/(rho_i*cp_i*dx)  ) # is it minus here?
				b_i = (2*h_br*dt/(rho_i*cp_i*dx)) # is it plus here?
				c_i = 2*k_i*dt/(rho_i*cp_i*dx*dx)
				print(str(a_i*T[i])+" "+str(c_i*T[i-1])+" "+str(b_i*Tbr))

				T_bc[i] = a_i*T[i]+ c_i*T[i-1] + b_i*Tbr
				print(T_bc[i])

	return T_bc


'''                           periodic  --> j, y-axis
(0,0)------------------------------------------
   | |                                         |
   | |periodic  zero degree initially  periodic| 0.33L
   v |                                         |
   i --------------------------------------(0.33L,L)
  x-axis          periodic '''
def BC_2D_Periodic_Dirichlet(T,  x,  mask):
	T1 = 3.5  # 500 # temperature at the left boundary
	x1 = 0.2*L
	y1 = 0.2*L
	T2 = 4  # temperature at the right boundary
	x2 = 0.8*L
	y2 = 0.2*L
	T3 = 6
	x3 = 0.2*L
	y3 = 0.8*L
	T4 = 8
	x4 = 0.8*L
	y4 = 0.8*L
	for i in range(len(mask)):
		for j in range(len(mask[0])):

			xi = x[i][j][0]
			yi = x[i][j][1]
			mask_ij = mask[i][j]

			len_T = len(T)-1

			if mask_ij ==0: # at ghost points, apply periodic B.C

				if xi<0: # at the upper ghost points
					i_index = int(len_T - 2*number_of_ghost_points)
					T[i][j] = T[i_index][j]
				elif xi>L: # at the bottom ghost points
					i_index = int(2*number_of_ghost_points)
					T[i][j] = T[i_index][j]
				elif yi<0: # at the left ghost points
					j_index = int(len_T - 2*number_of_ghost_points)
					T[i][j] = T[i][j_index]
				elif yi>L: # at the right ghost points
					j_index = int(2 * number_of_ghost_points)
					T[i][j] = T[i][j_index]
				else: # should not be here
					continue


			else: # in the domain, apply the four fixed temperature from T1 to T4
				if x1-dx/2<xi<x1+dx/2 and y1-dx/2<yi<y1+dx/2:
					T[i][j] = T1
				if x2 - dx/2 < xi < x2 +dx/2 and y2 -dx/2< yi < y2 +dx/2:
					T[i][j] = T2
				if x3-dx/2<xi<x3+dx/2 and y3-dx/2<yi<y3+dx/2:
					T[i][j] = T3
				if x4-dx/2<xi<x4+dx/2 and y4-dx/2<yi<y4+dx/2:
					T[i][j] = T4
				else:
					continue
	return T


def BC_2D_Dirichlet(T,x,  mask):  # TODO! check later
	Tbl = 1  # 500 # temperature at the left boundary
	Tbr = 0  # temperature at the right boundary
	for i in range(len(mask)):
		for j in range(len(mask[0])):

			xi = x[i][j][0]
			yi = x[i][j][1]

			mask_i = mask[i][j]

			if yi <= 0 + 10e-9:
				T[i][j] = Tbl
			elif yi >= L - 10e-9:
				T[i][j] = Tbr
			elif xi < 0:  # we are at the upper and bottom boundary where we have zero flux
				# so the value of the ghost points depends on the value of the interface point inside domain
				T[i][j] = T[i + 1][j]
			elif xi > L:
				T[i][j] = T[i - 1][j]
			else:
				continue

	return T

'''	T=30					       T=10000 if y>=0.75L      --> j, y-axis
(0,0)-----------------------------------------------------
   | |                                     					|
   | |Tbl=30  thirty degree initially  Tbr=reflective BC | L
   v |                                     					|
   i ----------------------------------------------------|(L,L)
  x-axis   T=30 '''

def BC_2D_Crater(T,x,  mask):  # TODO! check later
	Tamb = 303  # ambient temperature = 303 kelvin : ref
	T_crater = 10000  # temperature at the upper boundary if y>=0.75L : ref
	for i in range(len(mask)):
		for j in range(len(mask[0])):

			xi = x[i][j][0]
			yi = x[i][j][1]

			mask_i = mask[i][j]

			if yi <= 0 + 10e-9:
				T[i][j] = Tamb
			elif yi > L:
				T[i][j] = T[i][j-2]
			elif xi <= 0+10e-9:  
				if yi<0.75*L:
					T[i][j]=Tamb
				else:
					T[i][j]=T_crater
			elif xi >= L-10e-9:
				T[i][j] = Tamb
			else:
				continue

	return T

'''    insulation(zero flux)  --> j, y-axis
(0,0)--------------------------------------
   | |                                     |
   | |Tbl=500  zero degree initially  Tbr=0| 0.33L
   v |                                     |
   i --------------------------------------(0.33L,L)
  x-axis   insulation(zero flux) '''
def BC_2D_Dirichlet_2Layers(T, T_fine, x, x_fine, mask):# TODO! check later
	Tbl = 1#500 # temperature at the left boundary
	Tbr = 0 # temperature at the right boundary
	for i in range(len(mask)):
		for j in range(len(mask[0])):

			xi = x[i][j][0]
			yi = x[i][j][1]

			mask_i = mask[i][j]
		
			if yi <= 0+10e-9:
				T[i][j] = Tbl
			elif yi >= L-10e-9:
				T[i][j] = Tbr
			elif xi<0: # we are at the upper and bottom boundary where we have zero flux
				  # so the value of the ghost points depends on the value of the interface point inside domain
				T[i][j] = T[i+1][j]
			elif xi>L:
				T[i][j] = T[i-1][j]
			else:
				continue

	for i in range(len(T_fine)):
		for j in range(len(T_fine[0])):
			xi = x_fine[i][j][0]
			yi = x_fine[i][j][1]

			if yi <= 0 + 10e-9:
				T_fine[i][j] = Tbl
			elif yi >= L - 10e-9:
				T_fine[i][j] = Tbr
			elif xi < 0:  # we are at the upper and bottom boundary where we have zero flux
				# so the value of the ghost points depends on the value of the interface point inside domain
				T_fine[i][j] = T_fine[i + 1][j]
			elif xi > L:
				T_fine[i][j] = T_fine[i - 1][j]
			else:
				continue
	return T, T_fine
