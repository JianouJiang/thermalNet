# import libs
print("import libs")
import sys
sys.path.insert(0, '../../')
import numpy as np
from heat_transfer.initialConditions import *
from heat_transfer.boundaryConditions import *
from heat_transfer.linearSystemSolvers import *
from heat_transfer.crankNicolson import *
from heat_transfer.materialProperties import *

from tools.tools import *
parameters_directory="../parameters.txt"
L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points, num_of_timeSteps_for_plotting = readParameters(parameters_directory)

# evolve temperature
print("started crankNicolson_linear0_Dirichlet_TwoMaterials_Continuity_T().")
def evolveT_1D_Dirichlet_Linear0_TwoMaterials_Continuity():
  directory = "../../data/crankNicolson_linear0_Dirichlet_TwoMaterials_Continuity_T.txt"
  # if the file (i.e. directory) exists, delete it first, then we add new data to a brand new .txt file
  deleteFile(directory)
  # importing initial conditions
  t, x, T, mask, _lambda = IC_1D_Linear0_TwoMaterials()
  for i in range(len(t)):
    ti = t[i]
    # making boundary conditions
    T = BC_1D_Dirichlet_Tbl1_Tbr0(T, x, mask)
    
    # saving Temperature at t=n to .txt under /data
    writeData(directory, ti, x, T, _lambda)
    
    A, b = crankNicolson1D_Dirichlet_TwoMaterials(T, mask, _lambda, dx, dt)
    Tn, residual = LU_Decomposition(A,b,N=100, r=10**-6) 
    #print(Tn)
    
    # giving the new temperature to the old temperature for the next iteration
    T = Tn
    # getting the lambda based on the newly obtained temperature
    for i in range(len(T)):
      Ti = T[i]
      _lambda[3][i] = _lambda[i] #1 #_lambda_Aluminium(Ti)
  
  return

evolveT_1D_Dirichlet_Linear0_TwoMaterials_Continuity()
print("finished crankNicolson_linear0_Dirichlet_TwoMaterials_Continuity_T().")
