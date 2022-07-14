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
dt_for_plotting = t_max / num_of_timeSteps_for_plotting
plot_times = np.arange(0.0,t_max,dt_for_plotting)
# evolve temperature
print("started evolveT_2D_Dirichlet_Linear0_Aluminium().")
def evolveT_2D_Dirichlet_Linear0_Aluminium():
    '''    insulation(zero flux)  --> j, y-axis
(0,0)--------------------------------------
   | |                                     |
   | |Tbl=500  zero degree initially  Tbr=0| 0.33L
   v |                                     |
   i --------------------------------------(0.33L,L)
  x-axis   insulation(zero flux)
  '''

    st = time.time()
    directory = "../../data/crankNicolson2D_linear0_Dirichlet500_0_Aluminium_T.txt"
    # if the file (i.e. directory) exists, delete it first, then we add new data to a brand new .txt file
    deleteFile(directory)
    # importing initial conditions
    t, x, T, mask, _lambda = IC_2D_Linear0_Aluminium()
    for i in range(len(t)):
        ti = t[i]
        print("t=" + str(ti) + "s; t_max="+str(t_max))
        # making boundary conditions
        T = BC_2D_Dirichlet(T, x, mask)
        #print(np.array_str(T, precision=1, suppress_small=True))
        # saving Temperature at t=n to .txt under /data
        if ti in plot_times:
            writeData2D(directory, ti, x, T, _lambda)
        
        A, b = crankNicolson2D_Dirichlet(T, mask, _lambda, x, dt)
        
        # flatting the matrix
        # using array.flatten() method
        T_flat = T.flatten()
        #print(np.array_str(T_flat, precision=1, suppress_small=True))

        #Tn, residual_list_CG = Conjugate_Gradient(A, b, x0=T_flat, N=1024, reltol=1e-6, verbose=True) 
        Tn, residual_list_Jacobi = Jacobi(A, b, x0=T_flat, N=128, r=1e-6)
        #Tn, residual_list_GS = Gauss_Seidel(A, b, x0=T_flat, N=512, r=1e-6) # 4.58s
        #print(np.array_str(Tn, precision=1, suppress_small=True))
        # giving the new temperature to the old temperature for the next iteration
        # also, converting Tn in a vector to a 2d.array matrix T
        index = 0
        for i in range(len(T)):
            for j in range(len(T[0])):
                T[i][j] = Tn[index]
                index = index + 1
        print(np.array_str(T, precision=1, suppress_small=True))
        # getting the lambda based on the newly obtained temperature
        for i in range(len(T)):
            for j in range(len(T[0])):
                Tij = T[i][j]
                _lambda[3][i][j] = lambda_Aluminium(Tij)

    et = time.time()
    duration = et-st
    return duration

duration = evolveT_2D_Dirichlet_Linear0_Aluminium()
print("finished evolveT_2D_Dirichlet_Linear0_Aluminium() in " + str(duration) +"s.")
