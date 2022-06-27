# setting up initial conditions:
# importing libs:
print("importing libs")
import numpy as np
from math import *
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import os

# setting parameters:
print("setting parameters")
L = 1
dx = 0.01
t_max = 0.01
dt = 0.002
x = np.arange(0,L+dx,dx) 
T = np.arange(0,L+dx,dx) 
t = np.arange(0,t_max+dt,dt)
# https://en.wikipedia.org/wiki/Thermal_diffusivity, lambda = k/(cp*rho) with the unit m2/s
_lambda1 = 1.5
_lambda2 = 0.5








# defining initial conditions:
print("defining initial conditions")

# unit pulse function:
m = 1 # magnitude of the unit pulse function
def unitPulse(x): # input x is a np array
    T0 = np.arange(0,L+dx,dx) 
    for i in range(len(x)):
        xi = x[i]
        if 0.25*L<=xi<=0.75*L:
            T0[i] = m
        else:
            T0[i] = 0
    return T0 # output T0 is a np array

#T0 = unitPulse(x)
#print(T0)


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


# but we need to generate meshes first:
T = []
rho = []
Cp = []
k = []


def initialConditions(T, rho, Cp, k): # temperature, density, specific heat capacity and thermal conductivity
  
  # TODO
  
  return

