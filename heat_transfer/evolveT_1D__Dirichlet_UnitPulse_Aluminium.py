# import libs
from initialConditions import *
from boundaryConditions import *
from linearSystemSolvers import *
from AnalyticalHeatTransferSolution1D import *
from materialProperties import *



# evolve temperature
def evolveT_1D__Dirichlet_UnitPulse_Aluminium():
  directory = "../data/crankNicolson_unitPulse_Dirichlet_Aluminium_T.txt"
  # importing initial conditions
  t, x, T, mask, _lambda = IC_1D_UnitPulse_Aluminium()
  for ti in range(len(t)):
    # making boundary conditions
    T = BC_1D_Dirichlet(T, x, mask)
    A, b = crankNicolson1D_Dirichlet(T, mask, _lambda, dx, dt)
    Tn, residual = SOR(A, b, x0=None, N=100, r=10 ** -6, w=1.5)

    # saving Temperature at t=n to .txt under /data
    writeData(directory, ti, T, _lambda)

    # giving the new temperature to the old temperature for the next iteration
    T = Tn
    # getting the lambda based on the newly obtained temperature
    _lambda = _lambda_Aluminium(T)

    return

evolveT_1D__Dirichlet_UnitPulse_Aluminium()
print("finished evolveTemperature().")
