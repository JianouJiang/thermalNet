# import libs
print("import libs")
import time
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
print("started evolveT_1D_Dirichlet_Linear0_Aluminium().")
def evolveT_1D_Dirichlet_Linear0_Aluminium():
  st = time.time()
  directory = "../../data/crankNicolson_linear0_Dirichlet500_0_Aluminium_T.txt"
  # importing initial conditions
  t, x, T, mask, _lambda = IC_1D_Linear0_Aluminium()
  for i in range(len(t)):
    ti = t[i]
    print("t=" + str(ti) + "s; t_max="+str(t_max))
    # making boundary conditions
    T = BC_1D_Dirichlet(T, x, mask)
    
    # saving Temperature at t=n to .txt under /data
    writeData(directory, ti, x, T, _lambda)
    
    A, b = crankNicolson1D_Dirichlet(T, mask, _lambda, dx, dt)
    Tn, residual_list_CG = Conjugate_Gradient(A, b, x0=None, N=64, reltol=1e-6, verbose=True) # 2.08s
    

    # giving the new temperature to the old temperature for the next iteration
    T = Tn
    # getting the lambda based on the newly obtained temperature
    for i in range(len(T)):
      Ti = T[i]
      _lambda[3][i] = lambda_Aluminium(Ti)
  et = time.time()
  duration = et-st
  return duration

duration = evolveT_1D_Dirichlet_Linear0_Aluminium()
print("finished evolveT_1D_Dirichlet_Linear0_Aluminium() in " + str(duration) +"s.")
