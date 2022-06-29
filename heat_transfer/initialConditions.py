# setting up initial conditions:
# importing libs:
print("importing libs")
import numpy as np
from math import *
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import os
import sys
sys.path.insert(0, '../tools/')
from tools import *

# setting parameters:
print("setting parameters")
L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points = readParameters()# https://en.wikipedia.org/wiki/Thermal_diffusivity, lambda = k/(cp*rho) with the unit m2/s




# defining initial conditions:
print("defining initial conditions for temperature")

# unit pulse function:
m = 1 # magnitude of the unit pulse function
def unitPulse(x): # input x is a np array
    T0 = np.arange(0-number_of_ghost_points*dx,L+dx + number_of_ghost_points*dx,dx) 
    for i in range(len(x)):
        xi = x[i]
        if 0.25*L<=xi<=0.75*L:
            T0[i] = m
        else:
            T0[i] = 0
    return T0 # output T0 is a np array



# sines function: 
def sines(x): # input x is a np array
    T0 = np.arange(0,L+dx,dx) 
    m1 = 1 # magnitude of the 1st sine function
    m2 = 0.5 # 
    f1 = 1 # frequency of the 1st sine function
    f2 = 5
    for i in range(len(x)):
        xi = x[i]
        sin1 = m1*np.sin(f1*np.pi*xi/L)
        sin2 = m2*np.sin(f2*np.pi*xi/L) 
        T0[i] = sin1 + sin2
    return T0 # output T0 is a np array

# Linear function:
def linear(x): # input x is a np array
    T0 = np.arange(0,L+dx,dx) 
    for i in range(len(x)):
        xi = x[i]
        T0[i] = m*xi
    return T0 # output T0 is a np array    


# generate meshes first:
def IC_1D_UnitPulse_Aluminium(): 

    x = np.arange(-number_of_ghost_points*dx,L+dx + number_of_ghost_points*dx,dx) 
    mask = np.array([1 if 0<=xi<=L else 0 for xi in x])
    t = np.arange(0,t_max+dt,dt)

    T = unitPulse(x) # temperature
    
    rho = np.array([1.0 for i in range(len(T))])  # rho_Aluminium(T) # density
    Cp = np.array([1.0 for i in range(len(T))]) # Cp_Aluminium(T) # specific heat capacity
    k = np.array([1.0 for i in range(len(T))]) # k_Aluminium(T) # thermal conductivity
    _lambda = k/(Cp*rho)

    return t, x, T, mask, _lambda

