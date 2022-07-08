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
print("started evolveT_1D__Dirichlet_UnitPulse_Aluminium().")
def evolveT_1D_Dirichlet_UnitPulse_Aluminium():
  directory = "../../data/crankNicolson_unitPulse_Dirichlet_Aluminium_T.txt"
  # importing initial conditions
  t, x, T, mask, _lambda = IC_1D_UnitPulse_Aluminium()
  for i in range(len(t)):
    ti = t[i]
    # making boundary conditions
    T = BC_1D_Dirichlet(T, x, mask)
    
    # saving Temperature at t=n to .txt under /data
    writeData(directory, ti, x, T, _lambda)
    
    A, b = crankNicolson1D_Dirichlet(T, mask, _lambda, dx, dt)
    Tn, residual = SOR(A, b, x0=T, N=100, r=10 ** -9, w=1.1)# Conjugate_Gradient(A, b, x0=T, N=1024, reltol=1e-12, verbose=True) 
    #print(Tn) # Jacobi(A, b, x0=T, N=100, r=10 **-12)
    
    # giving the new temperature to the old temperature for the next iteration
    T = Tn
    # getting the lambda based on the newly obtained temperature
    for i in range(len(T)):
      Ti = T[i]
      _lambda[3][i] = 1 #_lambda_Aluminium(Ti)
  
  return

evolveT_1D_Dirichlet_UnitPulse_Aluminium()
print("finished evolveT_1D__Dirichlet_UnitPulse_Aluminium().")
