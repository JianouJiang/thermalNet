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
  
  Tnew = T
  
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

# Explicit forward time centred space (FTCS), first-order in time and second order convergence in space
def FTCS_Mixed(T, mask, _lambda, dx, dt): # if T=[Tbl, T1, T2, T3, Tbr] then mask=[0, 1, 1, 1, 0]
  print("FTCS_Mixed()") # left dirichlet right neumann (zero flux)

  length_T = len(T)
  length_mask = len(mask)
  
  Tnew = T
  
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

