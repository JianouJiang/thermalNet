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
plot_times = np.array([0.0, 5, 15,300])

# evolve temperature
print("started evolveT_1D_Dirichlet_Convec_Linear0_CokeInconel800HT_Continuity_FTCS().")
def evolveT_1D_Dirichlet_Linear0_TwoMaterials_Continuity_FTCS():
  st = time.time()
  directory = "../../data/FTCS_linear0_Dirichlet_Convec_CokeInconel800HT_Continuity_T.txt"
  # if the file (i.e. directory) exists, delete it first, then we add new data to a brand new .txt file
  deleteFile(directory)
  x_interface=0.7*L
  # importing initial conditions
  t, x, T, mask, _lambda = IC_1D_Linear0_TwoMaterials(x_interface)

  for i in range(len(t)):
    ti = t[i]
    print("t=" + str(ti) + "s; t_max=" + str(t_max))
    # making boundary conditions
    # put it in the FTCS_Dirichlet_TwoMaterials(T, mask, _lambda, x, dt, x_interface)

    #print(T)

    # saving Temperature at t=n to .txt under /data
    if ti in plot_times:
      print("writing data at "+str(ti))
      writeData(directory, ti, x, T, _lambda)
    
    Tn = FTCS_Dirichlet_TwoMaterials(T, mask, _lambda, x, dt, x_interface)


    # giving the new temperature to the old temperature for the next iteration
    T = Tn
    # getting the lambda based on the newly obtained temperature
    x_interface = 0.5*L
    valid=1
    for i in range(len(T)):
      Ti = T[i]
      if Ti>10e9 or Ti<-10e3:
        valid=0
      xi = x[i]
      if xi<x_interface:
        _lambda[2][i] = k_EthaneCoke(Ti) #  1.0 
        _lambda[3][i] = lambda_EthaneCoke(Ti) #  1.0 
      else:
        _lambda[2][i] = k_Inconel800HT(Ti) #  0.1
        _lambda[3][i] = lambda_Inconel800HT(Ti) # 0.1
    
    if valid!=1:
      print("invalid temperature, breaking the simulation...")
      break

  et = time.time()
  duration = et-st
  return duration

duration = evolveT_1D_Dirichlet_Linear0_TwoMaterials_Continuity_FTCS()
print("finished evolveT_1D_Dirichlet_Linear0_TwoMaterials_Continuity_FTCS() in " +str(duration)+" s.")
