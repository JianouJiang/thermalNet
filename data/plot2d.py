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
def plot_2D(plot_times, figureName, fileName):
    # appending data from fileName1 to x, T and _lambda at time in plot_times
    len_t_plot = len(plot_times)
    # plot the figure
    # generate 2 2d grids for the x & y bounds
    num_x = int(L/dx + 1)
    dy = dx
    num_y = int(L/dy + 1)
    x, y = np.mgrid[slice(-dx*number_of_ghost_points, L + dx*number_of_ghost_points + dx - 10e-9, dx),
                slice(-dy*number_of_ghost_points, L + dy*number_of_ghost_points + dy - 10e-9, dy)]
    T = x + y
    # x and y are bounds, so z should be the value *inside* those bounds.
    # Therefore, remove the last value from the z array.
    Ts = [[]] * len_t_plot # that contains len_t_plot of T
    r_plotting = 10*-6
    with open(fileName, 'r+') as fp:
        lines = fp.read().splitlines()
        x_index = 0
        y_index = 0
        for i in range(len(lines)):
            try:
                line = lines[i]
                nextline = lines[i+1]
                values = line.split()
                nextvalues = nextline.split()
                xi_next = float(nextvalues[1])

                ti = float(values[0])
                xi = float(values[1])
                yi = float(values[2])
                Ti = float(values[3])
                _lambda_i = float(values[4])

                x_temp = xi
                ti_index_temp = 0
                if ti in plot_times:
                    t_index = plot_times.tolist().index(ti)
                    if ti_index_temp != t_index:
                        x_index = 0
                        ti_index_temp = ti_index_temp + 1
                    
                    # appending data

                    if x_temp == xi_next:
                        y[x_index][y_index] = yi
                        x[x_index][y_index] = xi
                        T[x_index][y_index] = Ti

                        y_index = y_index + 1

                    else:
                        y[x_index][y_index] = yi
                        x[x_index][y_index] = xi
                        T[x_index][y_index] = Ti

                        x_index = x_index + 1  
                        y_index = 0

                    Ts[t_index] = T
            except IndexError: # reaching the last line
                continue     
        
    fig, axes = plt.subplots(nrows=len_t_plot, ncols=1)
     
    index = 0
    for i in range(len_t_plot):
        ti = plot_times[i]
        fig.set_figheight(15)
        fig.set_figwidth(5)
        plt.subplot(len_t_plot, 1, index+1)
        T = Ts[index]
        T_min, T_max = np.abs(T).min(), np.abs(T).max() # TODO
        plt.pcolor(x, y, T, cmap='RdBu_r', vmin=T_min, vmax=T_max)
        plt.title('T at t={}s'.format(ti))
        # set the limits of the plot to the limits of the data
        plt.axis([x.min(), x.max(), y.min(), y.max()])
        plt.colorbar()
        fig.tight_layout()
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        index = index + 1
        
    plt.title(figureName)
    save_directory = "../img/" + figureName + ".png"
    plt.savefig(save_directory)
    plt.show()

    return


plot_times = np.arange(0.0,t_max,dt_for_plotting)
figureName = 'Heat Transfer in 2D Aluminium'
fileName = "crankNicolson2D_linear0_Dirichlet500_0_Aluminium_T.txt"
plot_2D(plot_times, figureName, fileName)
print("finished plotting 2d...")
