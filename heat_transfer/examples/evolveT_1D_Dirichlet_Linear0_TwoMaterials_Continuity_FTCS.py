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
print("started evolveT_1D_Dirichlet_Linear0_TwoMaterials_Continuity_FTCS().")
def evolveT_1D_Dirichlet_Linear0_TwoMaterials_Continuity_FTCS():
  directory = "../../data/FTCS_linear0_Dirichlet_TwoMaterials_Continuity_T.txt"
  # if the file (i.e. directory) exists, delete it first, then we add new data to a brand new .txt file
  deleteFile(directory)
  # importing initial conditions
  t, x, T, mask, _lambda = IC_1D_Linear0_TwoMaterials()
  for i in range(len(t)):
    ti = t[i]
    print("t=" + str(ti) + "s; t_max=" + str(t_max))
    # making boundary conditions
    T = BC_1D_Dirichlet_Tbl1_Tbr0(T, x, mask)

    # saving Temperature at t=n to .txt under /data
    if ti in plot_times:
      writeData(directory, ti, x, T, _lambda)
    
    Tn = FTCS_Dirichlet_TwoMaterials(T, mask, _lambda, dx, dt)

    # giving the new temperature to the old temperature for the next iteration
    T = Tn
    # getting the lambda based on the newly obtained temperature
    for i in range(len(T)):
      Ti = T[i]
      _lambda[3][i] = 1 #_lambda_Aluminium(Ti)
  
  return

evolveT_1D_Dirichlet_Linear0_TwoMaterials_Continuity_FTCS()
print("finished evolveT_1D_Dirichlet_Linear0_TwoMaterials_Continuity_FTCS().")
