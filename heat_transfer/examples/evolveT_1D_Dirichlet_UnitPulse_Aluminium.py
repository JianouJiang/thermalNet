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

    Tn, residual_list_LU = LU_Decomposition(A,b,x=T,N=1024, r=10**-12)
    Tn, residual_list_Jacobi = Jacobi(A, b, x0=T, N=1024, r=10 **-12)
    w = 1.1
    
    Tn, residual_list_SOR = SOR(A, b, x0=T, N=1024, r=10 ** -12, w=w) # Conjugate_Gradient(A, b, x0=T, N=1024, reltol=1e-12, verbose=True) 
    Tn, residual_list_CG = Conjugate_Gradient(A, b, x0=T, N=1024, reltol=1e-12, verbose=True) 
    Tn, residual_list_GS = Gauss_Seidel(A, b, x0=T, N=1024, r=10 **-12)
    

    # plotting residual to show convergence
    '''
    plt.plot(residual_list_LU,"ob",label="LU-D")
    plt.plot(residual_list_Jacobi,"ok",label="Jacobi")
    plt.plot(residual_list_GS,"or",label="G-S")
    plt.plot(residual_list_SOR,"og",label="SOR (w={})".format(w))
    plt.plot(residual_list_CG,"oy",label="C-G")
    plt.xlabel('iterations')
    plt.ylabel('Residual')
    plt.yscale('log')
    plt.legend()
    plt.savefig("../../img/linearSystemSolverConvergence.png")
    plt.show()
    '''

    # giving the new temperature to the old temperature for the next iteration
    T = Tn
    # getting the lambda based on the newly obtained temperature
    for i in range(len(T)):
      Ti = T[i]
      _lambda[3][i] = 1 #_lambda_Aluminium(Ti)
  
  return

evolveT_1D_Dirichlet_UnitPulse_Aluminium()
print("finished evolveT_1D__Dirichlet_UnitPulse_Aluminium().")
