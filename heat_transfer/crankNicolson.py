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

      _lambda_ip1 = _lambda[3][i+1]
      _lambda_im1 = _lambda[3][i-1]

      if mask_im1==0: # at the left boundary, but in the domain
        A[i][i] = 1
        b[i] = T[i] 

      elif mask_ip1==0: # at the right boundary, but in the domain
        A[i][i] = 1
        b[i] = T[i] 

      else: # inside the domain
        ai = 1/dt + _lambda_i / (dx*dx)
        bi = - _lambda_ip1/(2*dx*dx)
        ci = - _lambda_im1/(2*dx*dx)
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

      _lambda_ip1 = _lambda[3][i+1]
      _lambda_im1 = _lambda[3][i-1]

      if mask_im1==0: # at the left boundary, but in the domain
        A[i][i] = 1
        b[i] = T[i] 

      else: # inside the domain
        ai = 1/dt + _lambda_i / (dx*dx)
        bi = - _lambda_ip1/(2*dx*dx)
        ci = - _lambda_im1/(2*dx*dx)
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
      
        A[i][i] = 1
        A[i][i-2] = -1
        b[i] = -T[i] + T[i-2]
    else:
      # in the domain
      mask_ip1 = mask[i+1]
      mask_im1 = mask[i-1]
      _lambda_ip1 = _lambda[3][i+1]
      _lambda_im1 = _lambda[3][i-1]
      if mask_im1==0: # at the left boundary, but in the domain
        A[i][i] = 1
        b[i] = Tbl

      else: # inside the domain
        ai = 1/dt + _lambda_i / (dx*dx)
        bi = - _lambda_ip1/(2*dx*dx)
        ci = - _lambda_im1/(2*dx*dx)
        fi = 1/dt - _lambda_i/(dx*dx)
        A[i][i] = ai
        A[i][i+1] = bi
        A[i][i-1] = ci
        b[i] = fi * T[i] - ci * T[i-1] - bi * T[i+1]


  return A, b






'''    insulation(zero flux)  --> j, y-axis
(0,0)--------------------------------------
   | |                                     |
   | |Tbl=500  zero degree initially  Tbr=0| 0.33L
   v |                                     |
   i --------------------------------------(0.33L,L)
  x-axis   insulation(zero flux)
  '''
def crankNicolson2D_Dirichlet(T, mask, _lambda, x, dt):
  print("crankNicolson2D_Dirichlet(x)")

  length_Tx = len(T)
  length_Ty = len(T[0])
  
  A = np.zeros((length_Tx*length_Ty, length_Tx*length_Ty))
  b = np.zeros(length_Tx*length_Ty)
  
  index = 0
  for i in range(len(T)):
    for j in range(len(T[0])):
      mask_ij = mask[i][j]

      _lambda_ij = _lambda[3][i][j]
      xij = x[i][j][0]
      yij = x[i][j][1]

      if mask_ij == 0:
      # at the ghost points
        if yij<0 or yij>L: # at the left or right ghost points, should be fixed temp, dirichlet
          A[index][index] = 1
          b[index] = T[i][j] 
        else: # at the top or bottom ghost points, should be zero flux
          if xij < 0: # at the upper ghost points
            A[index][index] = 1
            A[index+2][index] = -1
            b[index] = -T[i][j] + T[i+2][j]
          else: # at the lower ghost points
            A[index][index] = 1
            A[index-2][index] = -1
            b[index] = -T[i][j] + T[i-2][j]

      else:
        # in the domain
        mask_ijp1 = mask[i][j+1]
        mask_ijm1 = mask[i][j-1]
        mask_ip1j = mask[i+1][j]
        mask_im1j = mask[i-1][j]

        _lambda_ijp1 = _lambda[3][i][j+1]
        _lambda_ijm1 = _lambda[3][i][j-1]
        _lambda_ip1j = _lambda[3][i+1][j]
        _lambda_im1j = _lambda[3][i-1][j]
        
        if mask_ijm1==0: # at the left boundary, but in the domain
          A[index][index] = 1
          b[index] = T[i][j] 

        elif mask_ijp1==0: # at the right boundary, but in the domain
          A[index][index] = 1
          b[index] = T[i][j] 

        else: # inside the domain
          ai = 1/dt + 2 * _lambda_ij / (dx*dx)
          bi = - _lambda_ijp1/(2*dx*dx)
          ci = - _lambda_ijm1/(2*dx*dx)
          di = - _lambda_ip1j/(2*dx*dx)
          ei = - _lambda_im1j/(2*dx*dx)
          fi = 1/dt - 2*_lambda_ij/(dx*dx)
          A[index][index] = ai
          A[index][index+1] = bi
          A[index][index-1] = ci
          A[index][index+length_Ty] = di
          A[index][index-length_Ty] = ei
          
          b[index] = fi * T[i][j] - bi * T[i][j+1] - ci * T[i][j-1] - di * T[i+1][j] - ei * T[i-1][j]
      index = index + 1    
  return A, b