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
        Tnew[i] = T[i] # this is because we set the values in the B.C function
      else:
        gamma_i = _lambda_i * dt/(dx*dx)
        Tnew[i] = gamma_i * (- 2*T[i] + T[i-1] + T[i+1]) + T[i]
  
  return Tnew


def FTCS_Dirichlet_TwoMaterials(T, mask, _lambda, x, dt, x_interface=0.5*L):  # if T=[Tbl, T1, T2, T3, Tbr] then mask=[0, 1, 1, 1, 0]
  print("FTCS_Dirichlet_TwoMaterials()")

  #Tnew = np.zeros(len(T))
  Tnew = BC_1D_Dirichlet_Tbl500_Convection(T, x, _lambda, mask)

  for i in range(len(T)):
    mask_i = mask[i]
    
    _lambda_i = _lambda[3][i]
    if mask_i==0:
      continue
    else: # within the domain
      mask_ip1 = mask[i+1]
      mask_im1 = mask[i-1]

      if mask_ip1==0 or mask_im1==0: # we are at the left or right boundary, do not calculate here 
        continue                     # because they were obtained in the B.C function
      else:
        at_interface = 0
        mask_ip1 = mask[i+1]
        if mask_i!=mask_ip1 and mask_ip1!=0:
          at_interface = 1
        if at_interface:
          k_ip1 =_lambda[2][i+1]  # k[i+i] # 0.1
          k_im1 = _lambda[2][i-1]  # k[i] # 1
          _lambda_im1 = _lambda[3][i-1] # 1
          _lambda_ip1 = _lambda[3][i+1]  # _lambda[i+1] # 0.1

          bi_star = (-4 * k_ip1 * _lambda_im1 + 2 * (3 * k_im1 + k_ip1) * _lambda_ip1)
          ci_star = (2 * (k_im1 + 3 * k_ip1) * _lambda_im1 - 4 * k_im1 * _lambda_ip1)
          di_star = (k_ip1 * _lambda_im1 - (3 * k_im1 + 2 * k_ip1) * _lambda_ip1)
          ei_star = (-(2 * k_im1 + 3 * k_ip1) * _lambda_im1 + k_im1 * _lambda_ip1)

          gamma_i = - dt*2 / (12 * (k_ip1 + k_im1) * dx * dx)
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

      if mask_ij ==0:
        Tnew[i][j] = T[i][j]
      else:  # in the domain
        gamma_ij = _lambda_ij * dt / (dx * dx)
        Tnew[i][j] = gamma_ij * (- 4 * T[i][j] + T[i - 1][j] + T[i + 1][j] + T[i][j-1] + T[i][j + 1] ) + T[i][j]
      
        if Tnew[i][j] > 10e9 or Tnew[i][j] < -10e3:
          return None

  return Tnew


'''    insulation(zero flux)  --> j, y-axis
  (0,0)--------------------------------------
     | |  Inconel800HT |         Alu         |
     | |Tbl=500  zero degree initially  Tbr=0| L
     v |               |                     |
     i --------------------------------------(L,L)
    x-axis   insulation(zero flux)
'''
# Explicit forward time centred space (FTCS), first-order in time and second order convergence in space
def FTCS_Dirichlet_2D_TwoMaterials(T, T_fine, mask, _lambda, _lambda_fine, x, x_fine, dt):  # if T=[Tbl, T1, T2, T3, Tbr] then mask=[0, 1, 1, 1, 0]
  print("FTCS_Dirichlet_2D_TwoMaterials()")  # left dirichlet right neumann (zero flux)

  Tnew = np.zeros((len(T), len(T[0])))
  Tnew_fine = np.zeros(( int(len(T)*2-1), int(len(T[0])*2-1) )) # 2 times finer of the resolution
  Tamp=150
  for i_f in range(len(Tnew_fine)):
    for j_f in range(len(Tnew_fine[0])):
      # iterating at the fine mesh
      xij = x_fine[i_f][j_f][0]
      yij = x_fine[i_f][j_f][1]
      i = int(i_f / 2)
      j = int(j_f / 2)
      if yij <= 0 :  # at the left ghost points and left interface
        Tnew_fine[i_f][j_f] = T_fine[i_f][j_f+1]
      elif yij >= L :  # at the right ghost points and right interface
        Tnew_fine[i_f][j_f] = T_fine[i_f][j_f-1]
      elif xij < 0:  # at the upper ghost points
        Tnew_fine[i_f][j_f] = Tamp
      elif xij > L:  # at the lower ghost points
        Tnew_fine[i_f][j_f] = 500
      else:  # in the domain
        # maybe at the interface:
        if j_f % 2 != 0:
          i = int(i_f / 2)
          j = int(j_f / 2)
          mask_ij = mask[i][j]
          mask_ijp1 = mask[i][j + 1]

          if mask_ij != 0 and mask_ijp1 != 0 and mask_ij != mask_ijp1:  # right at the point of the vertical interface

            k_jp1 =  _lambda_fine[2][i_f][j_f + 1]  # k[i+i]
            k_jm1 = _lambda_fine[2][i_f][j_f -1]  # k[i]
            _lambda_jm1 = _lambda_fine[3][i_f][j_f-1]
            _lambda_jp1 = _lambda_fine[3][i_f][j_f + 1]  # _lambda[i+1]

            T_ij = T_fine[i_f][j_f]
            T_ijm1 = T_fine[i_f][j_f-1]  # T[i][j]
            T_ijp1 = T_fine[i_f][j_f + 1]  # T[i ][j+ 1]
            T_ijm2 = T_fine[i_f][j_f - 2]  # (T[i][j] + T[i][j - 1]) / 2
            T_ijp2 = T_fine[i_f][j_f + 2]  # (T[i ][j+ 1] + T[i ][j+ 2]) / 2

            bi_star = (-4 * k_jp1 * _lambda_jm1 + 2 * (3 * k_jm1 + k_jp1) * _lambda_jp1)
            ci_star = (2 * (k_jm1 + 3 * k_jp1) * _lambda_jm1 - 4 * k_jm1 * _lambda_jp1)
            di_star = (k_jp1 * _lambda_jm1 - (3 * k_jm1 + 2 * k_jp1) * _lambda_jp1)
            ei_star = (-(2 * k_jm1 + 3 * k_jp1) * _lambda_jm1 + k_jm1 * _lambda_jp1)

            dx_fine = dx / 2  # because we are operating on the fine mesh

            gamma_i = - dt * 2 / (12 * (k_jp1 + k_jm1) * dx_fine * dx_fine)
            Tnew_fine_ij = gamma_i * (ci_star * T_ijm1 + bi_star * T_ijp1 + ei_star * T_ijm2 + di_star * T_ijp2) + \
                           T_ij
            Tnew_fine[i_f][j_f] = Tnew_fine_ij
          else: # not at the interface
            _lambda_ij = _lambda_fine[3][i_f][j_f]
            dx_fine = dx / 2
            gamma_ij = _lambda_ij * dt / (dx_fine * dx_fine)  # because dx=dx/2 here in the fine mesh
            T_fine_ij = gamma_ij * (
                      - 4 * T_fine[i_f][j_f] + T_fine[i_f - 1][j_f] + T_fine[i_f + 1][j_f] + T_fine[i_f][j_f - 1] +
                      T_fine[i_f][j_f + 1]) + T_fine[i_f][j_f]
            Tnew_fine[i_f][j_f] = T_fine_ij

        # maybe at the interface
        elif i_f % 2 != 0:
          i = int(i_f / 2)
          j = int(j_f / 2)
          mask_ij = mask[i][j]
          mask_ip1j = mask[i+1][j]

          if mask_ij != 0 and mask_ip1j != 0 and mask_ij != mask_ip1j: # right at the hori interface
            k_ip1 = _lambda_fine[2][i_f + 1][j_f]  # k[i+i]
            k_im1 = _lambda_fine[2][i_f-1][j_f]  # k[i]
            _lambda_im1 = _lambda_fine[3][i_f-1][j_f]
            _lambda_ip1 = _lambda_fine[3][i_f + 1][j_f]  # _lambda[i+1]

            T_ij = T_fine[i_f][j_f]
            T_im1j = T_fine[i_f-1][j_f]  # T[i][j]
            T_ip1j = T_fine[i_f + 1][j_f]  # T[i+1][j]
            T_im2j = T_fine[i_f - 2][j_f]  # (T[i][j]+T[i-1][j])/2
            T_ip2j = T_fine[i_f + 2][j_f]  # (T[i+1][j] + T[i +2][j]) / 2

            bi_star = (-4 * k_ip1 * _lambda_im1 + 2 * (3 * k_im1 + k_ip1) * _lambda_ip1)
            ci_star = (2 * (k_im1 + 3 * k_ip1) * _lambda_im1 - 4 * k_im1 * _lambda_ip1)
            di_star = (k_ip1 * _lambda_im1 - (3 * k_im1 + 2 * k_ip1) * _lambda_ip1)
            ei_star = (-(2 * k_im1 + 3 * k_ip1) * _lambda_im1 + k_im1 * _lambda_ip1)

            dx_fine = dx / 2

            gamma_i = - dt * 2 / (12 * (k_ip1 + k_im1) * dx_fine * dx_fine)
            Tnew_fine_ij = gamma_i * (ci_star * T_im1j + bi_star * T_ip1j + ei_star * T_im2j + di_star * T_ip2j) + \
                           T_ij
            Tnew_fine[i_f][j_f] = Tnew_fine_ij
          else: # not at the interface
            _lambda_ij = _lambda_fine[3][i_f][j_f]
            dx_fine = dx / 2
            gamma_ij = _lambda_ij * dt / (dx_fine * dx_fine)  # because dx=dx/2 here in the fine mesh
            T_fine_ij = gamma_ij * (
                      - 4 * T_fine[i_f][j_f] + T_fine[i_f - 1][j_f] + T_fine[i_f + 1][j_f] + T_fine[i_f][j_f - 1] +
                      T_fine[i_f][j_f + 1]) + T_fine[i_f][j_f]
            Tnew_fine[i_f][j_f] = T_fine_ij


        else: # definitely not at the interface, and should be at the coarse point, calculate as usual
          _lambda_ij = _lambda_fine[3][i_f][j_f]
          dx_fine = dx/2
          gamma_ij = _lambda_ij * dt / (dx_fine * dx_fine) # because dx=dx/2 here in the fine mesh
          T_fine_ij = gamma_ij * (- 4 * T_fine[i_f][j_f] + T_fine[i_f - 1][j_f] + T_fine[i_f + 1][j_f] + T_fine[i_f][j_f - 1] + T_fine[i_f][j_f + 1]) + T_fine[i_f][j_f]
          Tnew_fine[i_f][j_f] = T_fine_ij
          Tnew[i][j] = T_fine_ij


  #for i in range(len(Tnew_fine[5])):
    #print(Tnew_fine[5][i])
  #print(np.array_str(Tnew, precision=3, suppress_small=True))
  return Tnew, Tnew_fine


'''    Converctive B.C h=20 W/m2/k --> j, y-axis
(0,0)--------------------------------------------------------
   | |               Inconel 800HT                          | L=3mm
   ----------------------------------------------------------
   | |zero heat flux   zero degree initially  zero heat flux| L_total= 0.01 m
   v |                Coke                                  |                     
                             -------------------------------        y=ax+b   (  a=L/(Lr-L), b=-L^2/(Lr-L)  )
   i -------------------------           T=500 degree        (L,L)
  x-axis   T=500 degree
'''
def FTCS_Dirichlet_2D_CokeInconel800HT(T, T_fine, mask, _lambda, _lambda_fine, x, x_fine, dt):  # if T=[Tbl, T1, T2, T3, Tbr] then mask=[0, 1, 1, 1, 0]
  print("FTCS_Dirichlet_2D_CokeInconel800HT()")
  h=20
  Tamb=150
  Tnew = np.zeros((len(T), len(T[0])))
  Tnew_fine = np.zeros(( int(len(T)*2-1), int(len(T[0])*2-1) )) # 2 times finer of the resolution

  for i_f in range(len(Tnew_fine)):
    for j_f in range(len(Tnew_fine[0])):
      # iterating at the fine mesh
      xij = x_fine[i_f][j_f][0]
      yij = x_fine[i_f][j_f][1]
      i = int(i_f / 2)
      j = int(j_f / 2)
      if yij <= 0 + 10e-9:  # at the left ghost points and left interface
        Tnew_fine[i_f][j_f] = T_fine[i_f][j_f]
      elif yij >= (L - 10e-9):  # at the right ghost points and right interface
        Tnew_fine[i_f][j_f] = T_fine[i_f][j_f]
      elif xij < 0:  # at the upper ghost points
        Tnew_fine[i_f][j_f] = T_fine[i_f + 1][j_f]
      elif xij > L:  # at the lower ghost points
        Tnew_fine[i_f][j_f] = T_fine[i_f - 1][j_f]
      else:  # in the domain
        # maybe at the interface:
        if j_f % 2 != 0:
          i = int(i_f / 2)
          j = int(j_f / 2)
          mask_ij = mask[i][j]
          mask_ijp1 = mask[i][j + 1]

          if mask_ij != 0 and mask_ijp1 != 0 and mask_ij != mask_ijp1:  # right at the point of the vertical interface

            k_jp1 =  _lambda_fine[2][i_f][j_f + 1]  # k[i+i]
            k_jm1 = _lambda_fine[2][i_f][j_f -1]  # k[i]
            _lambda_jm1 = _lambda_fine[3][i_f][j_f-1]
            _lambda_jp1 = _lambda_fine[3][i_f][j_f + 1]  # _lambda[i+1]

            T_ij = T_fine[i_f][j_f]
            T_ijm1 = T_fine[i_f][j_f-1]  # T[i][j]
            T_ijp1 = T_fine[i_f][j_f + 1]  # T[i ][j+ 1]
            T_ijm2 = T_fine[i_f][j_f - 2]  # (T[i][j] + T[i][j - 1]) / 2
            T_ijp2 = T_fine[i_f][j_f + 2]  # (T[i ][j+ 1] + T[i ][j+ 2]) / 2

            bi_star = (-4 * k_jp1 * _lambda_jm1 + 2 * (3 * k_jm1 + k_jp1) * _lambda_jp1)
            ci_star = (2 * (k_jm1 + 3 * k_jp1) * _lambda_jm1 - 4 * k_jm1 * _lambda_jp1)
            di_star = (k_jp1 * _lambda_jm1 - (3 * k_jm1 + 2 * k_jp1) * _lambda_jp1)
            ei_star = (-(2 * k_jm1 + 3 * k_jp1) * _lambda_jm1 + k_jm1 * _lambda_jp1)

            dx_fine = dx / 2  # because we are operating on the fine mesh

            gamma_i = - dt * 2 / (12 * (k_jp1 + k_jm1) * dx_fine * dx_fine)
            Tnew_fine_ij = gamma_i * (ci_star * T_ijm1 + bi_star * T_ijp1 + ei_star * T_ijm2 + di_star * T_ijp2) + \
                           T_ij
            Tnew_fine[i_f][j_f] = Tnew_fine_ij+dt*dx_fine*h*2/_lambda_ij*(Tamb-T_ij)
          else: # not at the interface
            _lambda_ij = _lambda_fine[3][i_f][j_f]
            dx_fine = dx / 2
            gamma_ij = _lambda_ij * dt / (dx_fine * dx_fine)  # because dx=dx/2 here in the fine mesh
            T_fine_ij = gamma_ij * (
                      - 4 * T_fine[i_f][j_f] + T_fine[i_f - 1][j_f] + T_fine[i_f + 1][j_f] + T_fine[i_f][j_f - 1] +
                      T_fine[i_f][j_f + 1]) + T_fine[i_f][j_f]
            Tnew_fine[i_f][j_f] = T_fine_ij

        # maybe at the interface
        elif i_f % 2 != 0:
          i = int(i_f / 2)
          j = int(j_f / 2)
          mask_ij = mask[i][j]
          mask_ip1j = mask[i+1][j]

          if mask_ij != 0 and mask_ip1j != 0 and mask_ij != mask_ip1j: # right at the hori interface
            k_ip1 = _lambda_fine[2][i_f + 1][j_f]  # k[i+i]
            k_im1 = _lambda_fine[2][i_f-1][j_f]  # k[i]
            _lambda_im1 = _lambda_fine[3][i_f-1][j_f]
            _lambda_ip1 = _lambda_fine[3][i_f + 1][j_f]  # _lambda[i+1]

            T_ij = T_fine[i_f][j_f]
            T_im1j = T_fine[i_f-1][j_f]  # T[i][j]
            T_ip1j = T_fine[i_f + 1][j_f]  # T[i+1][j]
            T_im2j = T_fine[i_f - 2][j_f]  # (T[i][j]+T[i-1][j])/2
            T_ip2j = T_fine[i_f + 2][j_f]  # (T[i+1][j] + T[i +2][j]) / 2

            bi_star = (-4 * k_ip1 * _lambda_im1 + 2 * (3 * k_im1 + k_ip1) * _lambda_ip1)
            ci_star = (2 * (k_im1 + 3 * k_ip1) * _lambda_im1 - 4 * k_im1 * _lambda_ip1)
            di_star = (k_ip1 * _lambda_im1 - (3 * k_im1 + 2 * k_ip1) * _lambda_ip1)
            ei_star = (-(2 * k_im1 + 3 * k_ip1) * _lambda_im1 + k_im1 * _lambda_ip1)

            dx_fine = dx / 2

            gamma_i = - dt * 2 / (12 * (k_ip1 + k_im1) * dx_fine * dx_fine)
            Tnew_fine_ij = gamma_i * (ci_star * T_im1j + bi_star * T_ip1j + ei_star * T_im2j + di_star * T_ip2j) + \
                           T_ij
            Tnew_fine[i_f][j_f] = Tnew_fine_ij
          else: # not at the interface
            _lambda_ij = _lambda_fine[3][i_f][j_f]
            dx_fine = dx / 2
            gamma_ij = _lambda_ij * dt / (dx_fine * dx_fine)  # because dx=dx/2 here in the fine mesh
            T_fine_ij = gamma_ij * (
                      - 4 * T_fine[i_f][j_f] + T_fine[i_f - 1][j_f] + T_fine[i_f + 1][j_f] + T_fine[i_f][j_f - 1] +
                      T_fine[i_f][j_f + 1]) + T_fine[i_f][j_f]
            Tnew_fine[i_f][j_f] = T_fine_ij


        else: # definitely not at the interface, and should be at the coarse point, calculate as usual
          _lambda_ij = _lambda_fine[3][i_f][j_f]
          dx_fine = dx/2
          gamma_ij = _lambda_ij * dt / (dx_fine * dx_fine) # because dx=dx/2 here in the fine mesh
          T_fine_ij = gamma_ij * (- 4 * T_fine[i_f][j_f] + T_fine[i_f - 1][j_f] + T_fine[i_f + 1][j_f] + T_fine[i_f][j_f - 1] + T_fine[i_f][j_f + 1]) + T_fine[i_f][j_f]
          Tnew_fine[i_f][j_f] = T_fine_ij
          Tnew[i][j] = T_fine_ij


  #for i in range(len(Tnew_fine[5])):
    #print(Tnew_fine[5][i])
  #print(np.array_str(Tnew, precision=3, suppress_small=True))
  return Tnew, Tnew_fine

