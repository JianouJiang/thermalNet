# importing libs
import numpy as np
import sys
sys.path.insert(0, '../../')
from tools.tools import *
from heat_transfer.boundaryConditions import *
from heat_transfer.initialConditions import *

parameters_directory="../../heat_transfer/parameters.txt"
L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points, num_of_timeSteps_for_plotting = readParameters(parameters_directory)
# crank-Nicolson function
# https://people.sc.fsu.edu/~jpeterson/5-CrankNicolson.pdf

# crankNicolson1D_Dirichlet
def crankNicolson1D_Dirichlet(T, mask, _lambda, dx, dt): # if T=[Tbl, T1, T2, T3, Tbr] then mask=[0, 1, 1, 1, 0]
  print("crankNicolson1D_Dirichlet()")

  length_T = len(T)
  length_mask = len(mask)
  
  A = np.zeros((length_T, length_T))
  b = np.zeros(length_T)
  
  if length_T!=length_mask:
    print("Error: len(T) must be the same as len(mask))")
    return None
  
  for i in range(len(T)):
    mask_i = mask[i]

    _lambda_i = _lambda[3][i]
    
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


def crankNicolson1D_Dirichlet_Convec(T, mask, _lambda, dx, dt): # if T=[Tbl, T1, T2, T3, Tbr] then mask=[0, 1, 1, 1, 0]
  print("crankNicolson1D_Dirichlet_Convec()") # dirichlet B.C at the left and Convection at the right with heat transfer coefficient h = ? and free stream temp of air=?

  length_T = len(T)
  length_mask = len(mask)
  
  A = np.zeros((length_T, length_T))
  b = np.zeros(length_T)
  
  if length_T!=length_mask:
    print("Error: len(T) must be the same as len(mask))")
    return None
  
  for i in range(len(T)):
    mask_i = mask[i]

    _lambda_i = _lambda[3][i]
    
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


# crankNicolson1D_Dirichlet_TwoMaterials
def crankNicolson1D_Dirichlet_TwoMaterials(T, mask, _lambda, dx, dt): # if T=[Tbl, T1, T2, T3, Tbr] then mask=[0, 1, 1, 1, 0]
  # print("crankNicolson1D_Dirichlet()")

  length_T = len(T)
  length_mask = len(mask)
  
  A = np.zeros((length_T, length_T))
  b = np.zeros(length_T)
  
  if length_T!=length_mask:
    print("Error: len(T) must be the same as len(mask)")
    return None
  
  for i in range(len(T)):

    mask_i = mask[i]
    # we can also get the material interface's location from the mask
    at_interface = 0
    if i == 52:
      at_interface = 1

    _lambda_i = _lambda[3][i]

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

      elif at_interface: # at the interface between two materials
        
        k_ip1 = 0.1 #k[i+i]
        k_i = 1 #k[i]
        _lambda_i = 1
        _lambda_ip1 = 0.1# _lambda[i+1]

        ai_star = 1/dt
        bi_star = (-4*k_ip1 * _lambda_i + 2*(3*k_i+k_ip1)*_lambda_ip1) / (12*(k_i+k_ip1)*dx*dx)
        ci_star = (2*(k_i+3*k_ip1)*_lambda_i-4*k_i*_lambda_ip1)/(12*(k_i+k_ip1)*dx*dx)
        di_star = (k_ip1*_lambda_i - (3*k_i+2*k_ip1)*_lambda_ip1 )/(12*(k_i+k_ip1)*dx*dx)
        ei_star = (-(2*k_i + 3*k_ip1)*_lambda_i + k_i*_lambda_ip1)/(12*(k_i+k_ip1)*dx*dx)
        fi_star = 1/dt

        A[i][i] = ai_star
        A[i][i+1] = bi_star
        A[i][i-1] = ci_star
        A[i][i+2] = di_star
        A[i][i-2] = ei_star
        b[i] = fi_star * T[i] - ci_star * T[i-1] - bi_star * T[i+1] - ei_star * T[i-2] - di_star * T[i+2]

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
def crankNicolson1D_Mixed(T, mask, _lambda, dx, dt): # if T=[Tbl, T1, T2, T3, Tbr] then mask=[0, 1, 1, 1, 0]
  print("crankNicolson1D_Mixed()") # left BC dirichelet T=0, right BC neumann with zero heat flux
  Tbl = 0
  heat_flux = 0

  length_T = len(T)
  length_mask = len(mask)
  
  A = np.zeros((length_T, length_T))
  b = np.zeros(length_T)
  
  if length_T!=length_mask:
    print("Error: len(T) must be the same as len(mask))")
    return None
  
  for i in range(len(T)):
    mask_i = mask[i]

    k_i = _lambda[2][i]
    _lambda_i = _lambda[3][i]
    
    if mask_i == 0:
      # at the left ghost point
      if i == 0:
        A[i][i] = 1
        b[i] = T[i]
      else: # at the right ghost point
        '''
        ai = 1/dt + _lambda_i / (dx*dx)
        ci = - _lambda_i/(2*dx*dx)
        fi = 1/dt - _lambda_i/(dx*dx)
        A[i][i] = ai
        A[i][i-1] = ci
        b[i] = fi * T[i] - ci * T[i-1] - 2*(_lambda_i /(2*dx*dx) * heat_flux/k_i *dx)
        '''
        A[i][i] = 1
        A[i][i-1] = -1
        b[i] = -T[i] + T[i-1]
    else:
      # in the domain
      mask_ip1 = mask[i+1]
      mask_im1 = mask[i-1]
      if mask_im1==0: # at the left boundary, but in the domain
        A[i][i] = 1
        b[i] = Tbl

      else: # inside the domain
        ai = 1/dt + _lambda_i / (dx*dx)
        bi = - _lambda_i/(2*dx*dx)
        ci = - _lambda_i/(2*dx*dx)
        fi = 1/dt - _lambda_i/(dx*dx)
        A[i][i] = ai
        A[i][i+1] = bi
        A[i][i-1] = ci
        b[i] = fi * T[i] - ci * T[i-1] - bi * T[i+1]
        '''
      elif mask_ip1==0: # at the right boundary, but in the domain
        ai = 1/dt + _lambda_i / (dx*dx)
        bi = - _lambda_i/(2*dx*dx)
        fi = 1/dt - _lambda_i/(dx*dx)
        A[i][i] = ai
        A[i][i+1] = bi
        b[i] = fi * T[i] - bi * T[i+1] - 2*(_lambda_i /(2*dx*dx) * heat_flux/k_i *dx)
      '''

  return A, b







def crankNicolson2D():
  # print("crankNicolson2D()")
  # TODO
  
  
  
  return
