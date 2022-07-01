print("importing libs")
import sys
sys.path.insert(0, '../tools/')
import numpy as np
from math import *
from tools import *
import matplotlib.pyplot as plt

L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points = readParameters()

# plot function
# exact solution, numerical solution, ansys solution, neural network solution
def plot_1D(plot_times, fileName1, fileName2=None, fileName3=None, fileName4=None):
    # appending data from fileName1 to x, T and _lambda at time in plot_times
    len_t_plot = len(plot_times)
    x_exact = [[]] * len_t_plot
    T_exact = [[]] * len_t_plot
    _lambda_exact = [[]] * len_t_plot
    r_plotting = 10*-6

    with open(fileName1, 'r+') as fp:
        lines = fp.read().splitlines()
        index = 0
        for i in range(len(lines)):
            line = lines[i]
            values = line.split()

            ti = float(values[0])
            xi = float(values[1])
            Ti = float(values[2])
            _lambda_i = float(values[3])


            ti_plot = float(plot_times[index])

            if ti==ti_plot:

                x_exact[index].append(xi)
                T_exact[index].append(Ti)
                _lambda_exact[index].append(_lambda_i)
            else:
                index = index + 1
                x_exact[index].append(xi)
                T_exact[index].append(Ti)
                _lambda_exact[index].append(_lambda_i)

    x_numerical = [[]] * len_t_plot
    T_numerical = [[]] * len_t_plot
    _lambda_numerical = [[]] * len_t_plot
    if (fileName2!=None):
        with open(fileName2, 'r+') as fp:
            lines = fp.read().splitlines()
            index = 0
            for i in range(len(lines)):
                line = lines[i]
                values = line.split()
                ti = float(values[0])
                xi = float(values[1])
                Ti = float(values[2])
                _lambda_i = float(values[3])

                ti_plot = float(plot_times[index])
                if ti==ti_plot:
                    x_numerical[index].append(xi)
                    T_numerical[index].append(Ti)
                    _lambda_numerical[index].append(_lambda_i)
                else:
                    index = index + 1
                    x_numerical[index].append(xi)
                    T_numerical[index].append(Ti)
                    _lambda_numerical[index].append(_lambda_i)

    x_ansys = [[]] * len_t_plot
    T_ansys = [[]] * len_t_plot
    _lambda_ansys = [[]] * len_t_plot
    if (fileName3!=None):
        with open(fileName3, 'r+') as fp:
            lines = fp.read().splitlines()
            index = 0
            for i in range(len(lines)):
                line = lines[i]
                values = line.split()
                ti = float(values[0])
                xi = float(values[1])
                Ti = float(values[2])
                _lambda_i = float(values[3])

                ti_plot = float(plot_times[index])
                if ti==ti_plot:
                    x_ansys[index].append(xi)
                    x_ansys[index].append(Ti)
                    _lambda_ansys[index].append(_lambda_i)
                else:
                    index = index + 1
                    x_ansys[index].append(xi)
                    x_ansys[index].append(Ti)
                    _lambda_ansys[index].append(_lambda_i)
                
    x_PINN = [[]] * len_t_plot
    T_PINN = [[]] * len_t_plot
    _lambda_PINN = [[]] * len_t_plot
    if (fileName4!=None):
        with open(fileName4, 'r+') as fp:
            lines = fp.read().splitlines()
            index = 0
            for i in range(len(lines)):
                line = lines[i]
                values = line.split()
                ti = float(values[0])
                xi = float(values[1])
                Ti = float(values[2])
                _lambda_i = float(values[3])

                ti_plot = float(plot_times[index])
                if ti==ti_plot:
                    x_PINN[index].append(xi)
                    T_PINN[index].append(Ti)
                    _lambda_PINN[index].append(_lambda_i)
                else:
                    index = index + 1
                    x_PINN[index].append(xi)
                    T_PINN[index].append(Ti)
                    _lambda_PINN[index].append(_lambda_i)
    
                
    # plot the figure
    plt.figure(figsize=(7,5))
    color_list = ['k','r','b','g','y']
    index = 0
    for ti in plot_times:
        dot_color = '.'+ color_list[index]
        star_color = "*" + color_list[index]
        cross_color = "x" + color_list[index]
        circle_color = "o" + color_list[index]
        plt.plot(x_exact[index],T_exact[index],dot_color,label='analytic at t={}s'.format(ti),markersize=3)
        if len(x_numerical[index])!=0:
            plt.plot(x_numerical[index],T_numerical[index],star_color,label='numerical at t={}s'.format(ti),markersize=3)
        if len(x_ansys[index])!=0:
            plt.plot(x_ansys[index],T_ansys[index],cross_color,label='ansys at t={}s'.format(ti),markersize=3)
        if len(x_PINN[index])!=0:
            plt.plot(x_PINN[index],T_PINN[index],circle_color,label='PINN at t={}s'.format(ti),markersize=3)
        plt.legend(fontsize=5)
        index = index + 1
        
    plt.xlabel('x (m)',fontsize=12)
    plt.ylabel('T (k)',fontsize=12)
    plt.title('Temperature Evolution in 1D Aluminium Rod')
    plt.savefig('../img/Temperature_Evolution_in_1D_Aluminium_Rod.png')
    plt.show()
    return


plot_times = np.arange(0.0,t_max,dt)
fileName1 = "exact_unitPulse_Dirichlet_T.txt"
fileName2 = "crankNicolson_unitPulse_Dirichlet_Aluminium_T.txt"
plot_1D(plot_times, fileName1, fileName2)
print("finished plotting...")
