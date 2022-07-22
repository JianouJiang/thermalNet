# import libs 
print("import libs")
import sys
sys.path.insert(0, '../../')
import numpy as np
from heat_transfer.initialConditions import *
from heat_transfer.boundaryConditions import *
from heat_transfer.linearSystemSolvers import *
from heat_transfer.crankNicolson import *
from heat_transfer.explicitDiscretisation import *
from heat_transfer.materialProperties import *

from tools.tools import *
parameters_directory="../parameters.txt"
L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points, num_of_timeSteps_for_plotting = readParameters(parameters_directory)
dt_for_plotting = t_max / num_of_timeSteps_for_plotting
plot_times = np.arange(0.0,t_max,dt_for_plotting)
plot_times = np.array([0.0, 0.001,0.05,1.0])
# evolve temperature
print("started evolveT_2D_Dirichlet_Linear0_CokeInconel800HT_Continuity_FTCS.")
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
def evolveT_2D_Dirichlet_Linear0_CokeInconel800HT_Continuity_FTCS():
  st = time.time()
  directory = "../../data/FTCS_2D_linear0_Dirichlet_Convec_CokeInconel800HT_T.txt"
  # if the file (i.e. directory) exists, delete it first, then we add new data to a brand new .txt file
  deleteFile(directory)
  # importing initial conditions
  t, x, x_fine, T, T_fine, mask, _lambda, _lambda_fine = IC_2D_Linear0_CokeInconel800HT()

  for i in range(len(t)):
    ti = t[i]
    print("t=" + str(ti) + "s; t_max=" + str(t_max))
    # making boundary conditions
    T, T_fine = BC_2D_Dirichlet_Convec_2Layers(T, T_fine, x, x_fine, mask)
    #print(np.array_str(T, precision=2, suppress_small=True))

    # saving Temperature at t=n to .txt under /data
    if ti in plot_times:
      writeData2D(directory, ti, x, T, _lambda)

    Tn, Tn_fine = FTCS_Dirichlet_2D_CokeInconel800HT(T, T_fine, mask, _lambda, _lambda_fine, x, x_fine, dt)

    # giving the new temperature to the old temperature for the next iteration
    T = Tn
    T_fine = Tn_fine
    valid = 1
    # getting the lambda based on the newly obtained temperature for the coarse mesh
    for i in range(len(T)):
      for j in range(len(T[0])):
        Tij = T[i][j]
        # check if temp seems valid
        if Tij>10e6 or Tij<-10e3:
          valid=0
        mask_ij = mask[i][j]
        if mask_ij == 1:
          _lambda[2][i][j] = k_Aluminium(Tij)
          _lambda[3][i][j] = lambda_Aluminium(Tij)
        elif mask_ij==2:
          _lambda[2][i][j] = k_Inconel800HT(Tij)
          _lambda[3][i][j] = lambda_Inconel800HT(Tij)
        else:
          continue # this is at the ghost point, only matters if there is convective B.C
    if valid==0:
      print("temperature unstable, breaking the simulation...")
      break

    # getting the lambda based on the newly obtained temperature for the fine mesh
    for i_f in range(len(T_fine)):
      for j_f in range(len(T_fine[0])):
        if i_f%2==0 and j_f%2==0: # at the coarse point
          # getting the indexes for the coarse mesh
          i=int(i_f/2)
          j=int(j_f/2)
          _lambda_fine[2][i_f][j_f] = _lambda[2][i][j]
          _lambda_fine[3][i_f][j_f] = _lambda[3][i][j]
        elif i_f%2==0 and j_f%2!=0: # point right in between T[i][j] and T[i][j+1]
          i = int(i_f / 2)
          j = int(j_f / 2)
          _lambda_fine[2][i_f][j_f] = (_lambda[2][i][j-1] + _lambda[2][i][j])/2
          _lambda_fine[3][i_f][j_f] = (_lambda[3][i][j-1] + _lambda[3][i][j])/2
        elif i_f%2!=0 and j_f%2==0: # point right in between T[i][j] and T[i+1][j]
          i = int(i_f / 2)
          j = int(j_f / 2)
          _lambda_fine[2][i_f][j_f] = (_lambda[2][i][j] + _lambda[2][i-1][j])/2
          _lambda_fine[3][i_f][j_f] = (_lambda[3][i][j] + _lambda[3][i-1][j])/2
        else: # this point is at the centre of four coarse points, so should be the avg of those four
          i = int(i_f / 2)
          j = int(j_f / 2)
          _lambda_fine[2][i_f][j_f] = (_lambda[2][i][j] + _lambda[2][i-1][j] +_lambda[2][i][j-1] + _lambda[2][i-1][j-1] ) / 4
          _lambda_fine[3][i_f][j_f] = (_lambda[3][i][j] + _lambda[3][i-1][j] +_lambda[3][i][j-1] + _lambda[3][i-1][j-1] ) / 4

  et = time.time()
  duration = et - st
  return duration

duration = evolveT_2D_Dirichlet_Linear0_CokeInconel800HT_Continuity_FTCS()

print("finished evolveT_2D_Dirichlet_Linear0_CokeInconel800HT_Continuity_FTCS() in " + str(duration) +"s.")

