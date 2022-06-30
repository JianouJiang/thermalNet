# import libs
print("import libs")

import numpy as np
from initialConditions import IC_1D_UnitPulse_Aluminium
from boundaryConditions import BC_1D_Dirichlet
from linearSystemSolvers import Jacobi, SOR
from crankNicolson import crankNicolson1D_Dirichlet
from materialProperties import _lambda_Aluminium
import sys
sys.path.insert(0, '../tools/')
from tools import *
L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points = readParameters()

# evolve temperature
print("started evolveT_1D__Dirichlet_UnitPulse_Aluminium().")
def evolveT_1D_Dirichlet_UnitPulse_Aluminium():
  directory = "../data/crankNicolson_unitPulse_Dirichlet_Aluminium_T.txt"
  # importing initial conditions
  t, x, T, mask, _lambda = IC_1D_UnitPulse_Aluminium()
  for i in range(len(t)):
    ti = t[i]
    # making boundary conditions
    T = BC_1D_Dirichlet(T, x, mask)
    
    # saving Temperature at t=n to .txt under /data
    writeData(directory, ti, T, _lambda)
    
    A, b = crankNicolson1D_Dirichlet(T, mask, _lambda, dx, dt)
    Tn, residual = SOR(A, b, x0=None, N=100, r=10 ** -6, w=1.5)
    print(Tn)
    
    # giving the new temperature to the old temperature for the next iteration
    T = Tn
    # getting the lambda based on the newly obtained temperature
    for i in range(len(T)):
      Ti = T[i]
      _lambda[i] = _lambda_Aluminium(Ti)
  
  return

evolveT_1D_Dirichlet_UnitPulse_Aluminium()
print("finished evolveT_1D__Dirichlet_UnitPulse_Aluminium().")
