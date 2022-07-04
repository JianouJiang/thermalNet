# importing libs
import numpy as np
import sys
sys.path.insert(0, '../../')
from tools.tools import *
from heat_transfer.boundaryConditions import *
from heat_transfer.initialConditions import *

L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points, num_of_timeSteps_for_plotting = readParameters()
# crank-Nicolson function
# https://people.sc.fsu.edu/~jpeterson/5-CrankNicolson.pdf

# crankNicolson1D_Dirichlet
def crankNicolson1D_Dirichlet(T, mask, _lambda, dx, dt): # if T=[Tbl, T1, T2, T3, Tbr] then mask=[0, 1, 1, 1, 0]
  # print("crankNicolson1D_Dirichlet()")

  length_T = len(T)
  length_mask = len(mask)
  
  A = np.zeros((length_T, length_T))
  b = np.zeros(length_T)
  
  if length_T!=length_mask:
    print("Error: len(T) must be the same as len(mask))")
    return None
  
  for i in range(len(T)):
    mask_i = mask[i]

    _lambda_i = _lambda[i]
    
    if mask_i == 0:
      # at the ghost cell
      A[i][i] = 1
      b[i] = T[i] 

    else:
      # in the domain
      mask_ip1 = mask[i+1]
      mask_im1 = mask[i-1]
      if mask_im1==0: # at the left boundary, but in the domain
        A[i][i] = 1
        b[i] = T[i] 

      elif mask_ip1==0: # at the right boundary, but in the domain
        A[i][i] = 1
        b[i] = T[i] 

      else: # inside the domain
        ai = 1/dt + _lambda_i / (dx*dx)
        bi = - _lambda_i/(2*dx*dx)
        ci = - _lambda_i/(2*dx*dx)
        fi = 1/dt - _lambda_i/(dx*dx)
        A[i][i] = ai
        A[i][i+1] = bi
        A[i][i-1] = ci
        b[i] = fi * T[i] - ci * T[i-1] - bi * T[i+1]
  
  return A, b




# crankNicolson1D_Neumann

def crankNicolson1D_Neumann(T, mask): # if T=[Tbl, T1, T2, T3, Tbr] then mask=[0, 1, 1, 1, 0]
  # print("crankNicolson1D_Dirichlet()")
  # TODO
  length_T = len(T)
  length_mask = len(mask)
  
  A = np.zeros((length_T, length_T))
  b = np.zeros(length_T)
  
  if length_T!=length_mask:
    print("Error: len(T) must be the same as len(mask))")
    return None
  
  for i in range(len(T)):
    mask_i = mask[i]
    mask_ip1 = mask[i+1]
    mask_im1 = mask[i-1]
    
    if mask_i == 0:
      # at the ghost cell
      pass
    else:
      # in the domain
      if mask_ip1==0 or mask_im1==0:
        # at the boundary, but in the domain
        pass
      else:
        # inside the domain
        pass
    
  
  return A, b


# crankNicolson1D_Mixed

def crankNicolson1D_Mixed(T, mask): # if T=[Tbl, T1, T2, T3, Tbr] then mask=[0, 1, 1, 1, 0]
  # print("crankNicolson1D_Mixed()")
  # TODO
  length_T = len(T)
  length_mask = len(mask)
  
  A = np.zeros((length_T, length_T))
  b = np.zeros(length_T)
  
  if length_T!=length_mask:
    print("Error: len(T) must be the same as len(mask))")
    return None
  
  for i in range(len(T)):
    mask_i = mask[i]
    mask_ip1 = mask[i+1]
    mask_im1 = mask[i-1]
    
    if mask_i == 0:
      # at the ghost cell
      pass
    else:
      # in the domain
      if mask_ip1==0 or mask_im1==0:
        # at the boundary, but in the domain
        pass
      else:
        # inside the domain
        pass
    
  
  return A, b






def crankNicolson2D():
  # print("crankNicolson2D()")
  # TODO
  
  
  
  return
