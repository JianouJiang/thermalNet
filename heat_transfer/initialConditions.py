# setting up initial conditions:
# importing libs:
print("importing libs for initialConditions.py")
import numpy as np
from math import *
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import os
import sys
sys.path.insert(0, '../../')
from tools.tools import *
from heat_transfer.materialProperties import *

# setting parameters:
print("setting parameters")
parameters_directory="../../heat_transfer/parameters.txt"
L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points, num_of_timeSteps_for_plotting = readParameters(parameters_directory)# https://en.wikipedia.org/wiki/Thermal_diffusivity, lambda = k/(cp*rho) with the unit m2/s




# defining initial conditions:
print("defining initial conditions for temperature")

# unit pulse function:
m = 500 # magnitude of the unit pulse function
def unitPulse(xi,yi=None): 
    if 0.25*L<=xi<=0.75*L:
        T0 = m
    else:
        T0 = 0
    return T0 



# sines function: 
def sines(xi,yi=None): # input x is a np array
    m1 = 1 # magnitude of the 1st sine function
    m2 = 0.5 # 
    f1 = 1 # frequency of the 1st sine function
    f2 = 5
    sin1 = m1*np.sin(f1*np.pi*xi/L)
    sin2 = m2*np.sin(f2*np.pi*xi/L) 
    T0 = sin1 + sin2
    return T0 

# Linear function:
def linear(xi,yi=None): # input x is a np array
    
    xi = x[i]
    T0 = m*xi
    return T0 

def linear0(xi,yi=None): # input x is a np array
    T0 = 0.0
    return T0 

# generate meshes first:
def IC_1D_UnitPulse_Aluminium(): 

    x = np.arange(-number_of_ghost_points*dx,L+dx + number_of_ghost_points*dx-10e-9,dx) 

    mask = np.array([1 if 0<=xi<=L else 0 for xi in x])
    t = np.arange(0,t_max+dt,dt)

    T = np.array([0 for xi in x])
    for xi in x:
        T[i] = unitPulse(xi) # temperature

    rho = np.array([1.0 for i in range(len(T))])
    Cp = np.array([1.0 for i in range(len(T))]) 
    k = np.array([1.0 for i in range(len(T))])
    _lambda = k/(Cp*rho)
    for i in range(len(T)):
        Ti=T[i]
        rho[i] = rho_Aluminium(Ti) # np.array([1.0 for i in range(len(T))]) # # density 
        Cp[i] = Cp_Aluminium(Ti)#np.array([1.0 for i in range(len(T))]) #  # specific heat capacity
        k[i] = k_Aluminium(Ti)#np.array([1.0 for i in range(len(T))]) #  # thermal conductivity
        _lambda[i] = k[i] /(Cp[i] *rho[i])
    #print(_lambda)
    return t, x, T, mask, np.array([rho, Cp, k, _lambda])

def IC_1D_Linear0_Aluminium(): 

    x = np.arange(-number_of_ghost_points*dx,L+dx + number_of_ghost_points*dx-10e-9,dx) 

    mask = np.array([1 if 0<=xi<=L else 0 for xi in x])
    t = np.arange(0,t_max+dt,dt)

    T = np.array([0 for xi in x])
    for xi in x:
        T[i] = linear0(xi) # temperature
    rho = np.array([1.0 for i in range(len(T))])
    Cp = np.array([1.0 for i in range(len(T))]) 
    k = np.array([1.0 for i in range(len(T))])
    _lambda = k/(Cp*rho)
    for i in range(len(T)):
        Ti=T[i]
        rho[i] = rho_Aluminium(Ti) # np.array([1.0 for i in range(len(T))]) # # density 
        Cp[i] = Cp_Aluminium(Ti)#np.array([1.0 for i in range(len(T))]) #  # specific heat capacity
        k[i] = k_Aluminium(Ti)#np.array([1.0 for i in range(len(T))]) #  # thermal conductivity
        _lambda[i] = k[i] /(Cp[i] *rho[i])
    #print(_lambda)
    return t, x, T, mask, np.array([rho, Cp, k, _lambda])

def IC_1D_Linear0_Convect_Aluminium(): 
    h_br = 20 # heat transfer coefficient at the right hand side
    x = np.arange(-number_of_ghost_points*dx,L+dx + number_of_ghost_points*dx-10e-9,dx) 
    mask = np.array([1 if 0<=xi<=L else 0 for xi in x])
    t = np.arange(0,t_max+dt,dt)

    T = np.array([0 for xi in x])
    for xi in x:
        T[i] = linear0(xi) # temperature

    rho = np.array([1.0 for i in range(len(T))])
    Cp = np.array([1.0 for i in range(len(T))]) 
    k = np.array([1.0 for i in range(len(T))])
    _lambda = k/(Cp*rho)
    for i in range(len(T)):
        xi = x[i]
        Ti=T[i]
        rho[i] = rho_Aluminium(Ti) # np.array([1.0 for i in range(len(T))]) # # density 
        Cp[i] = Cp_Aluminium(Ti)#np.array([1.0 for i in range(len(T))]) #  # specific heat capacity
        k[i] = k_Aluminium(Ti)#np.array([1.0 for i in range(len(T))]) #  # thermal conductivity
        if xi<L:
            k[i] = k_Aluminium(Ti)
        else:
            k[i] = dx * h_br
        _lambda[i] = k[i] /(Cp[i] *rho[i])
    return t, x, T, mask, np.array([rho, Cp, k, _lambda])


def IC_1D_Sines_Aluminium(): 

    x = np.arange(-number_of_ghost_points*dx,L+dx + number_of_ghost_points*dx-10e-9,dx) 

    mask = np.array([1 if 0<=xi<=L else 0 for xi in x])
    t = np.arange(0,t_max+dt,dt)

    T = np.array([0 for xi in x])
    for xi in x:
        T[i] = sines(xi) # temperature
    
    rho = np.array([1.0 for i in range(len(T))])
    Cp = np.array([1.0 for i in range(len(T))]) 
    k = np.array([1.0 for i in range(len(T))])
    _lambda = k/(Cp*rho)
    for i in range(len(T)):
        Ti=T[i]
        rho[i] = rho_Aluminium(Ti) # np.array([1.0 for i in range(len(T))]) # # density 
        Cp[i] = Cp_Aluminium(Ti)#np.array([1.0 for i in range(len(T))]) #  # specific heat capacity
        k[i] = k_Aluminium(Ti)#np.array([1.0 for i in range(len(T))]) #  # thermal conductivity
        _lambda[i] = k[i] /(Cp[i] *rho[i])

    return t, x, T, mask, np.array([rho, Cp, k, _lambda])

def IC_1D_Linear0_TwoMaterials():
    x_interface = 0.5*L

    x = np.arange(-number_of_ghost_points * dx, L + dx + number_of_ghost_points * dx-10e-9, dx)

    mask = np.array([1 if 0 <= xi <= L else 0 for xi in x])
    t = np.arange(0, t_max + dt, dt)

    T = np.array([0 for xi in x])
    for xi in x:
        T[i] = linear0(xi) # temperature

    rho = np.array([1.0 for i in range(len(T))])  # rho_Aluminium(T) # density
    Cp = np.array([1.0 for i in range(len(T))])  # Cp_Aluminium(T) # specific heat capacity
    k = np.array([1.0 for i in range(len(T))])  # k_Aluminium(T) # thermal conductivity
    _lambda = np.array([1.0 if xi < x_interface else 0.1 for xi in x])

    return t, x, T, mask, np.array([rho, Cp, k, _lambda])

def IC_2D_Linear0_Aluminium(): 
    t = np.arange(0,t_max+dt,dt)
    num_points_x = L / dx + 1
    dy = dx
    num_points_y = L / dy + 1
    x = np.ones(( int(num_points_x + 2*number_of_ghost_points), int(num_points_y + 2*number_of_ghost_points)))

    mask = np.ones(( int(num_points_x + 2*number_of_ghost_points), int(num_points_y + 2*number_of_ghost_points)))
    T = np.ones(( int(num_points_x + 2*number_of_ghost_points), int(num_points_y + 2*number_of_ghost_points)))
    rho = np.ones(( int(num_points_x + 2*number_of_ghost_points), int(num_points_y + 2*number_of_ghost_points)))
    Cp = np.ones(( int(num_points_x + 2*number_of_ghost_points), int(num_points_y + 2*number_of_ghost_points)))
    k =np.ones(( int(num_points_x + 2*number_of_ghost_points), int(num_points_y + 2*number_of_ghost_points)))
    _lambda =np.ones(( int(num_points_x + 2*number_of_ghost_points), int(num_points_y + 2*number_of_ghost_points)))
    x = np.array(x ,dtype = object)
    
    
    for i in range(len(T)):
        xij = -number_of_ghost_points * dx + dx * i
        for j in range(len(T[0])):
            yij = -number_of_ghost_points * dy + dy * j

            if (xij<0 or xij>L):
                mask[i][j] = 0
            elif (yij<0 or yij>L):
                mask[i][j] = 0
            else:
                mask[i][j] = 1
            x[i][j] = [xij,yij]
            Tij = linear0(xij,yij) # temperature
            T[i][j] = Tij
            rho[i][j] = rho_Aluminium(Tij) # np.array([1.0 for i in range(len(T))]) # # density 
            Cp[i][j] = Cp_Aluminium(Tij)#np.array([1.0 for i in range(len(T))]) #  # specific heat capacity
            k[i][j] = k_Aluminium(Tij)#np.array([1.0 for i in range(len(T))]) #  # thermal conductivity
            _lambda[i][j] = k[i][j] /(Cp[i][j] *rho[i][j])
    return t, x, T, mask, np.array([rho, Cp, k, _lambda])