print("importing libs")
import sys
sys.path.insert(0, '../')
import numpy as np
from math import *
from tools.tools import *
import matplotlib.pyplot as plt

parameters_directory="../heat_transfer/parameters.txt"
L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points, num_of_timeSteps_for_plotting = readParameters(parameters_directory)
dt_for_plotting = t_max / num_of_timeSteps_for_plotting
# plot function
# exact solution, numerical solution, ansys solution, neural network solution
def plot_1D(plot_times, figureName, fileName1, fileName2=None, fileName3=None, fileName4=None):
    # appending data from fileName1 to x, T and _lambda at time in plot_times
    len_t_plot = len(plot_times)
    x_exact = [[]] * len_t_plot
    T_exact = [[]] * len_t_plot
    _lambda_exact = [[]] * len_t_plot
    r_plotting = 10*-6

    with open(fileName1, 'r+') as fp:
        lines = fp.read().splitlines()
        for i in range(len(lines)):
            line = lines[i]
            values = line.split()

            ti = float(values[0])
            xi = float(values[1])
            Ti = float(values[2])
            _lambda_i = float(values[3])



            if ti in plot_times:

                index = plot_times.tolist().index(ti)
                x_exact[index].append(xi)
                T_exact[index].append(Ti)
                _lambda_exact[index].append(_lambda_i)
                

    x_numerical = [[]] * len_t_plot
    T_numerical = [[]] * len_t_plot
    _lambda_numerical = [[]] * len_t_plot
    if (fileName2!=None):
        with open(fileName2, 'r+') as fp:
            lines = fp.read().splitlines()
            for i in range(len(lines)):
                line = lines[i]
                values = line.split()
                ti = float(values[0])
                xi = float(values[1])
                Ti = float(values[2])
                _lambda_i = float(values[3])

                if ti in plot_times:
                    index = plot_times.tolist().index(ti)
                    x_numerical[index].append(xi)
                    T_numerical[index].append(Ti)
                    _lambda_numerical[index].append(_lambda_i)
                

    x_ansys = [[]] * len_t_plot
    T_ansys = [[]] * len_t_plot
    _lambda_ansys = [[]] * len_t_plot
    if (fileName3!=None):
        with open(fileName3, 'r+') as fp:
            lines = fp.read().splitlines()
            for i in range(len(lines)):
                line = lines[i]
                values = line.split()
                ti = float(values[0])
                xi = float(values[1])
                Ti = float(values[2])
                _lambda_i = float(values[3])

                if ti in plot_times:
                    index = plot_times.tolist().index(ti)
                    x_ansys[index].append(xi)
                    x_ansys[index].append(Ti)
                    _lambda_ansys[index].append(_lambda_i)
               
                
    x_PINN = [[]] * len_t_plot
    T_PINN = [[]] * len_t_plot
    _lambda_PINN = [[]] * len_t_plot
    if (fileName4!=None):
        with open(fileName4, 'r+') as fp:
            lines = fp.read().splitlines()
            for i in range(len(lines)):
                line = lines[i]
                values = line.split()
                ti = float(values[0])
                xi = float(values[1])
                Ti = float(values[2])
                _lambda_i = float(values[3])

                if ti==plot_times:
                    index = plot_times.tolist().index(ti)
                    x_PINN[index].append(xi)
                    T_PINN[index].append(Ti)
                    _lambda_PINN[index].append(_lambda_i)
                
                
    # plot the figure
    plt.figure(figsize=(7,5))
    color_list = ['k','r','b','y','g']
    index = 0
    for ti in plot_times:
        dot_color = '.'+ "k"#color_list[index]
        star_color = "*" + "r"#color_list[index]
        cross_color = "x" + "b"#color_list[index]
        circle_color = "o" + "y"#color_list[index]
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
    plt.title(figureName)
    save_directory = "../img/" + figureName + ".png"
    plt.savefig(save_directory)
    plt.show()
    return


plot_times = np.arange(0.0,t_max,dt_for_plotting)
figureName = 'Temperature Evolution in Two Materials (1D)'
fileName1 = "crankNicolson_linear0_Dirichlet_TwoMaterials_T.txt"
fileName2 = "crankNicolson_linear0_Dirichlet_TwoMaterials_Continuity_T.txt"
plot_1D(plot_times, figureName, fileName1, fileName2)
print("finished plotting...")
