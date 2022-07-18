# Chennuo: EvolveT 1D Mixed Sines Aluminium

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
# evolve temperature
print("started evolveT_2D_Dirichlet_Linear0_Aluminium_TwoMaterials_Continuity_FTCS().")
'''    insulation(zero flux)  --> j, y-axis
(0,0)--------------------------------------
   | |  Inconel800HT |         Alu         |
   | |Tbl=500  zero degree initially  Tbr=0| L
   v |               |                     |
   i --------------------------------------(L,L)
  x-axis   insulation(zero flux)
'''
def evolveT_2D_Dirichlet_Linear0_Aluminium_TwoMaterials_Continuity_FTCS():
  st = time.time()
  directory = "../../data/FTCS_2D_linear0_Dirichlet500_0_TwoMaterials_T.txt"
  # if the file (i.e. directory) exists, delete it first, then we add new data to a brand new .txt file
  deleteFile(directory)
  # importing initial conditions
  t, x, T, mask, _lambda = IC_2D_Linear0_TwoMaterials()

  for i in range(len(t)):
    ti = t[i]
    print("t=" + str(ti) + "s; t_max=" + str(t_max))
    # making boundary conditions
    T = BC_2D_Dirichlet(T, x, mask)
    
    # saving Temperature at t=n to .txt under /data
    if ti in plot_times:
      writeData2D(directory, ti, x, T, _lambda)

    Tn = FTCS_Dirichlet_2D_TwoMaterials(T, mask, _lambda, x, dt)
    
    # giving the new temperature to the old temperature for the next iteration
    T = Tn
    # getting the lambda based on the newly obtained temperature
    for i in range(len(T)):
      for j in range(len(T[0])):
        Tij = T[i][j]
        mask_ij = mask[i][j]
        if mask_ij == 1:
          _lambda[2][i][j] = k_Aluminium(Tij)
          _lambda[3][i][j] = lambda_Aluminium(Tij)
        elif mask_ij==2:
          _lambda[2][i][j] = k_Inconel800HT(Tij)
          _lambda[3][i][j] = lambda_Inconel800HT(Tij)
        else:
          continue # this is at the ghost point, only matters if there is convective B.C

  et = time.time()
  duration = et - st
  return duration

duration = evolveT_2D_Dirichlet_Linear0_Aluminium_FTCS()

print("finished evolveT_2D_Dirichlet_Linear0_Aluminium_FTCS() in " + str(duration) +"s.")

