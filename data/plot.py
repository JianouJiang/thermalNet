print("importing libs")
import sys
sys.path.insert(0, '../tools/')
import numpy as np
from math import *
import matplotlib.pyplot as plt

L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points = readParameters()

# plot function
# exact solution, numerical solution, ansys solution, neural network solution
def plot(plot_times, fileName1, fileName2=None, fileName3=None, fileName3=None):  
    # appending data from fileName1 to x, T and _lambda at time in plot_times
    x_exact = []
    T_exact = []
    _lambda_exact = []
    with open(fileName1) as fp:
        Lines = fp.readlines()
        for i in range(len(Lines)):
            line = Lines[i]
            ti = line[0]
            xi = line[1]
            Ti = line[2]
            _lambda_i = line[3]
            
            index = 0
            ti_plot = plot_times[index]
            if ti==ti_plot:
                x_exact.append(xi)
                T_exact.append(Ti)
                _lambda_exact.append(_lambda_i)
                index = index + 1
      
    x_numerical = []
    T_numerical = []
    _lambda_numerical = []
    if (fileName2!=None):
        with open(fileName2) as fp:
        Lines = fp.readlines()
        for i in range(len(Lines)):
            line = Lines[i]
            ti = line[0]
            xi = line[1]
            Ti = line[2]
            _lambda_i = line[3]
            
            index = 0
            ti_plot = plot_times[index]
            if ti==ti_plot:
                x_numerical.append(xi)
                T_numerical.append(Ti)
                _lambda_numerical.append(_lambda_i)
                index = index + 1
      
    x_ansys = []
    T_ansys = []
    _lambda_ansys = []    
    if (fileName3!=None):
        with open(fileName3) as fp:
        Lines = fp.readlines()
        for i in range(len(Lines)):
            line = Lines[i]
            ti = line[0]
            xi = line[1]
            Ti = line[2]
            _lambda_i = line[3]
            
            index = 0
            ti_plot = plot_times[index]
            if ti==ti_plot:
                x_ansys.append(xi)
                T_ansys.append(Ti)
                _lambda_ansys.append(_lambda_i)
                index = index + 1
                
    x_PINN = []
    T_PINN = []
    _lambda_PINN = []    
    if (fileName4!=None):
        with open(fileName3) as fp:
        Lines = fp.readlines()
        for i in range(len(Lines)):
            line = Lines[i]
            ti = line[0]
            xi = line[1]
            Ti = line[2]
            _lambda_i = line[3]
            
            index = 0
            ti_plot = plot_times[index]
            if ti==ti_plot:
                x_PINN.append(xi)
                T_PINN.append(Ti)
                _lambda_PINN.append(_lambda_i)
                index = index + 1
    
                
    # plot the figure
    plt.figure(figsize=(7,5))
    color_list = ['k','r','b','g','y']
    index = 0
    for ti in plot_times:
        dot_color = 'o'+ color_list[index]
        star_color = "*" + color_list[index]
        plt.plot(x_exact,T_exact,dot_color,label='analytic at t={}s'.format(ti),markersize=3)
        if len(x_numerical)!=0:
            plt.plot(x_numerical,T_numerical,colori,label='numerical at t={}s'.format(ti),markersize=3)
        if len(x_ansys)!=0:
            plt.plot(x_numerical,T_numerical,colori,label='ansys at t={}s'.format(ti),markersize=3)
        if len(x_PINN)!=0:
            plt.plot(x_numerical,T_numerical,colori,label='PINN at t={}s'.format(ti),markersize=3)
        plt.legend(fontsize=12) 
        index = index + 1
   
    index = index + 1
    plt.xlabel('x (m)',fontsize=12)
    plt.ylabel('T (k)',fontsize=12)
    plt.title('Temperature Evolution in 1D Aluminium Rod'.format(_lambda))
    plt.savefig('../img/Temperature_Evolution_in_1D_Aluminium_Rod.png')
    plt.show()
    return


plot_times = np.arange(0.0,t_max,dt)
fileName1 =  "exact_unitPulse_Neumann_T.txt"
plot(plot_times, fileName1)
print("finished plotting...")
