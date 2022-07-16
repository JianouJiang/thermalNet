# importing libs
import numpy as np
import sys
sys.path.insert(0, '../../')
from tools.tools import *
from heat_transfer.boundaryConditions import *
from heat_transfer.initialConditions import *

parameters_directory="../../heat_transfer/parameters.txt"
L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points, num_of_timeSteps_for_plotting = readParameters(parameters_directory)

# Explicit forward time centred space (FTCS), first-order in time and second order convergence in space
def FTCS_Dirichlet(T, mask, _lambda, dx, dt): # if T=[Tbl, T1, T2, T3, Tbr] then mask=[0, 1, 1, 1, 0]
  print("FTCS_Dirichlet()")

  length_T = len(T)
  length_mask = len(mask)

  Tnew = np.zeros(len(T))
  
  if length_T!=length_mask:
    print("Error: len(T) must be the same as len(mask))")
    return None
  
  for i in range(len(T)):
    mask_i = mask[i]
    _lambda_i = _lambda[3][i]
    
    if mask_i == 0:
      # at the ghost cell
      Tnew[i] = T[i]

    else: # inside the domain
      mask_ip1 = mask[i+1]
      mask_im1 = mask[i-1]
      if mask_ip1==0 or mask_im1==0:
        Tnew[i] = T[i]
      else:
        gamma_i = _lambda_i * dt/(dx*dx)
        Tnew[i] = gamma_i * (- 2*T[i] + T[i-1] + T[i+1]) + T[i]
  
  return Tnew


def FTCS_Dirichlet_TwoMaterials(T, mask, _lambda, dx, dt):  # if T=[Tbl, T1, T2, T3, Tbr] then mask=[0, 1, 1, 1, 0]
  print("FTCS_Dirichlet_TwoMaterials()")

  Tnew = np.zeros(len(T))

  for i in range(len(T)):
    at_interface = 0
    if i == 52:
      at_interface = 1
    mask_i = mask[i]
    _lambda_i = _lambda[3][i]
    if mask_i==0:
      continue
    else: # within the domain
      if at_interface:
        k_ip1 =0.1# _lambda[2][i+1]  # k[i+i]
        k_im1 =1# _lambda[2][i-1]  # k[i]
        _lambda_im1 =1# _lambda[3][i-1]
        _lambda_ip1 =0.1# _lambda[3][i+1]  # _lambda[i+1]

        bi_star = (-4 * k_ip1 * _lambda_im1 + 2 * (3 * k_im1 + k_ip1) * _lambda_ip1)
        ci_star = (2 * (k_im1 + 3 * k_ip1) * _lambda_im1 - 4 * k_im1 * _lambda_ip1)
        di_star = (k_ip1 * _lambda_im1 - (3 * k_im1 + 2 * k_ip1) * _lambda_ip1)
        ei_star = (-(2 * k_im1 + 3 * k_ip1) * _lambda_im1 + k_im1 * _lambda_ip1)

        gamma_i = - dt / (12 * (k_ip1 + k_im1) * dx * dx)
        Tnew[i] = gamma_i * (ci_star * T[i - 1] + bi_star * T[i + 1] + ei_star * T[i - 2] + di_star * T[i + 2]) + T[i]
      else:
        gamma_i = _lambda_i * dt / (dx * dx)
        Tnew[i] = gamma_i * (- 2 * T[i] + T[i - 1] + T[i + 1]) + T[i]
  return Tnew


# Explicit forward time centred space (FTCS), first-order in time and second order convergence in space
def FTCS_Neumann(T, mask, _lambda, dx, dt): # if T=[Tbl, T1, T2, T3, Tbr] then mask=[0, 1, 1, 1, 0]
  print("FTCS_Neumann()") # zero heat flux at both sides

  length_T = len(T)
  length_mask = len(mask)
  
  Tnew = np.zeros(len(T))
  
  if length_T!=length_mask:
    print("Error: len(T) must be the same as len(mask))")
    return None
  
  for i in range(len(T)):
    mask_i = mask[i]
    _lambda_i = _lambda[3][i]
    
    if mask_i == 0:
      # at the ghost points
      Tnew[i] = T[i]

      if i==0:    # at the left ghost point
        Tnew[i] = T[i+1]
      else:    # at the right ghost point
        Tnew[i] = T[i-1]

    else: # inside the domain
      gamma_i = _lambda_i * dt/(dx*dx)
      Tnew[i] = gamma_i * (- 2*T[i] + T[i-1] + T[i+1]) + T[i]
  
  return Tnew

# Explicit forward time centred space (FTCS), first-order in time and second order convergence in space
def FTCS_Mixed(T, mask, _lambda, dx, dt): # if T=[Tbl, T1, T2, T3, Tbr] then mask=[0, 1, 1, 1, 0]
  print("FTCS_Mixed()") # left dirichlet right neumann (zero flux)

  length_T = len(T)
  length_mask = len(mask)
  
  Tnew = np.zeros(len(T))
  
  if length_T!=length_mask:
    print("Error: len(T) must be the same as len(mask))")
    return None
  
  for i in range(len(T)):
    mask_i = mask[i]
    _lambda_i = _lambda[3][i]
    
    if mask_i == 0:
      # at the ghost cell
      Tnew[i] = T[i]

    else: # inside the domain
      mask_ip1 = mask[i+1]
      mask_im1 = mask[i-1]
      if mask_im1==0:
        Tnew[i] = T[i]
      else:
        gamma_i = _lambda_i * dt/(dx*dx)
        Tnew[i] = gamma_i * (- 2*T[i] + T[i-1] + T[i+1]) + T[i]
  
  return Tnew


# Explicit forward time centred space (FTCS), first-order in time and second order convergence in space
def FTCS_Dirichlet_2D(T, mask, _lambda, x, dt):  # if T=[Tbl, T1, T2, T3, Tbr] then mask=[0, 1, 1, 1, 0]
  print("FTCS_Dirichlet_2D()")  # left dirichlet right neumann (zero flux)
  '''   
   insulation(zero flux)  --> j, y-axis
(0,0)--------------------------------------
 | |                                     |
 | |Tbl=500  zero degree initially  Tbr=0| 0.33L
 v |                                     |
 i --------------------------------------(0.33L,L)
x-axis   insulation(zero flux)
 '''

  Tnew = np.zeros((len(T), len(T[0])))  # Tnew = np.zeros(len(T))?

  for i in range(len(T)):
    for j in range(len(T[0])):
      mask_ij = mask[i][j]

      _lambda_ij = _lambda[3][i][j]
      xij = x[i][j][0]
      yij = x[i][j][1]

      if yij <= 0 + 10e-9:  # at the left ghost points and left interface
        Tnew[i][j] = T[i][j]
      elif yij >= (L - 10e-9):  # at the right ghost points and right interface
        Tnew[i][j] = T[i][j]
      elif xij < 0:  # at the upper ghost points
        Tnew[i][j]= T[i+1][j]
      elif xij > L:  # at the lower ghost points
        Tnew[i][j]= T[i-1][j]

      else:  # in the domain
        gamma_ij = _lambda_ij * dt / (dx * dx)
        Tnew[i][j] = gamma_ij * (- 4 * T[i][j] + T[i - 1][j] + T[i + 1][j] + T[i][j-1] + T[i][j + 1] ) + T[i][j]


  return Tnew


