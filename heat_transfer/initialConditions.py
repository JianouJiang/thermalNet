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

def IC_1D_Linear0_TwoMaterials(x_interface=0.5*L):
    #x_interface = 0.5*L

    x = np.arange(-number_of_ghost_points * dx, L + dx + number_of_ghost_points * dx-10e-9, dx)

    mask = np.array([1 if (0-10e-9) <= xi <= (L+10e-9) else 0 for xi in x])
    t = np.arange(0, t_max + dt, dt)

    T = np.array([0 for xi in x])
    rho = np.array([1.0 for i in range(len(T))])  # rho_Aluminium(T) # density
    Cp = np.array([1.0 for i in range(len(T))])  # Cp_Aluminium(T) # specific heat capacity
    k = np.array([1.0 if xi < x_interface else 0.1 for xi in x])  # k_Aluminium(T) # thermal conductivity
    _lambda = np.array([1.0 if xi < x_interface else 0.1 for xi in x])
    
    for i in range(len(x)):
        xi = x[i]
        T[i] = linear0(xi) # temperature
        mask_i = mask[i]
        if mask_i!=0:
            if xi<x_interface: # material 1
                mask[i] = 1
            else: # material 2
                mask[i] = 2

        if xi<x_interface: # material 1: deposited coke
            rho[i] = rho_EthaneCoke(T[i]) #1
            Cp[i] = Cp_EthaneCoke(T[i]) #1
            k[i] = k_EthaneCoke(T[i]) #1
            _lambda[i] = lambda_EthaneCoke(T[i]) #1
                
        else: # material 2: Inconel800HT wall
            rho[i] = rho_Inconel800HT(T[i]) #1
            Cp[i] = Cp_Inconel800HT(T[i]) #1
            k[i] = k_Inconel800HT(T[i]) #1
            _lambda[i] = lambda_Inconel800HT(T[i]) #0.1    



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


def IC_2D_Linear0_Inconel800HT():
    t = np.arange(0, t_max + dt, dt)
    num_points_x = L / dx + 1
    dy = dx
    num_points_y = L / dy + 1
    x = np.ones((int(num_points_x + 2 * number_of_ghost_points), int(num_points_y + 2 * number_of_ghost_points)))

    mask = np.ones((int(num_points_x + 2 * number_of_ghost_points), int(num_points_y + 2 * number_of_ghost_points)))
    T = np.ones((int(num_points_x + 2 * number_of_ghost_points), int(num_points_y + 2 * number_of_ghost_points)))
    rho = np.ones((int(num_points_x + 2 * number_of_ghost_points), int(num_points_y + 2 * number_of_ghost_points)))
    Cp = np.ones((int(num_points_x + 2 * number_of_ghost_points), int(num_points_y + 2 * number_of_ghost_points)))
    k = np.ones((int(num_points_x + 2 * number_of_ghost_points), int(num_points_y + 2 * number_of_ghost_points)))
    _lambda = np.ones((int(num_points_x + 2 * number_of_ghost_points), int(num_points_y + 2 * number_of_ghost_points)))
    x = np.array(x, dtype=object)

    for i in range(len(T)):
        xij = -number_of_ghost_points * dx + dx * i
        for j in range(len(T[0])):
            yij = -number_of_ghost_points * dy + dy * j

            if (xij < 0 or xij > L):
                mask[i][j] = 0
            elif (yij < 0 or yij > L):
                mask[i][j] = 0
            else:
                mask[i][j] = 1
            x[i][j] = [xij, yij]
            Tij = linear0(xij, yij)  # temperature
            T[i][j] = Tij
            rho[i][j] = rho_Inconel800HT(Tij)  # np.array([1.0 for i in range(len(T))]) # # density
            Cp[i][j] = Cp_Inconel800HT(Tij)  # np.array([1.0 for i in range(len(T))]) #  # specific heat capacity
            k[i][j] = k_Inconel800HT(Tij)  # np.array([1.0 for i in range(len(T))]) #  # thermal conductivity
            _lambda[i][j] = k[i][j] / (Cp[i][j] * rho[i][j])
    return t, x, T, mask, np.array([rho, Cp, k, _lambda])


def IC_2D_Linear0_Sandstone(): 
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
            Tij = 303 #linear0(xij,yij) # temperature
            T[i][j] = Tij
            rho[i][j] = rho_Sandstone(Tij) # np.array([1.0 for i in range(len(T))]) # # density 
            Cp[i][j] = Cp_Sandstone(Tij)#np.array([1.0 for i in range(len(T))]) #  # specific heat capacity
            k[i][j] = k_Sandstone(Tij)#np.array([1.0 for i in range(len(T))]) #  # thermal conductivity
            _lambda[i][j] = k[i][j] /(Cp[i][j] *rho[i][j])
    return t, x, T, mask, np.array([rho, Cp, k, _lambda])


def IC_2D_Linear0_TwoMaterials():
    t = np.arange(0, t_max + dt, dt)
    num_points_x = L / dx + 1
    dy = dx
    num_points_y = L / dy + 1
    x = np.ones((int(num_points_x + 2 * number_of_ghost_points), int(num_points_y + 2 * number_of_ghost_points)))
    x_fine = np.ones((int( (num_points_x + 2 * number_of_ghost_points)*2-1 ), int( (num_points_y + 2 * number_of_ghost_points) *2-1) ))
    mask = np.ones((int(num_points_x + 2 * number_of_ghost_points), int(num_points_y + 2 * number_of_ghost_points)))
    T = np.zeros((int(num_points_x + 2 * number_of_ghost_points), int(num_points_y + 2 * number_of_ghost_points)))
    T_fine = np.zeros((int((num_points_x + 2 * number_of_ghost_points) * 2 - 1),
                      int((num_points_y + 2 * number_of_ghost_points) * 2 - 1)))
    rho = np.ones((int(num_points_x + 2 * number_of_ghost_points), int(num_points_y + 2 * number_of_ghost_points)))
    rho_fine = np.ones((int((num_points_x + 2 * number_of_ghost_points) * 2 - 1),
                      int((num_points_y + 2 * number_of_ghost_points) * 2 - 1)))
    Cp = np.ones((int(num_points_x + 2 * number_of_ghost_points), int(num_points_y + 2 * number_of_ghost_points)))
    Cp_fine = np.ones((int((num_points_x + 2 * number_of_ghost_points) * 2 - 1),
                      int((num_points_y + 2 * number_of_ghost_points) * 2 - 1)))
    k = np.ones((int(num_points_x + 2 * number_of_ghost_points), int(num_points_y + 2 * number_of_ghost_points)))
    k_fine = np.ones((int((num_points_x + 2 * number_of_ghost_points) * 2 - 1),
                      int((num_points_y + 2 * number_of_ghost_points) * 2 - 1)))
    _lambda = np.ones((int(num_points_x + 2 * number_of_ghost_points), int(num_points_y + 2 * number_of_ghost_points)))
    _lambda_fine = np.ones((int((num_points_x + 2 * number_of_ghost_points) * 2 - 1),
                      int((num_points_y + 2 * number_of_ghost_points) * 2 - 1)))
    x = np.array(x, dtype=object)
    x_fine = np.array(x_fine, dtype=object)

    for i_f in range(len(T_fine)):
        x_f_ij = -number_of_ghost_points*2 * dx/2 + dx/2 * i_f
        for j_f in range(len(T_fine[0])):
            y_f_ij = -number_of_ghost_points*2 * dy/2 + dy/2 * j_f

            # append data into fine mesh
            if (x_f_ij < 0 or x_f_ij > L):
                x_fine[i_f][j_f] = [x_f_ij, y_f_ij]
                T_f_ij = linear0(x_f_ij, y_f_ij)  # temperature
                T_fine[i_f][j_f] = T_f_ij
            elif (y_f_ij < 0 or y_f_ij > L):
                x_fine[i_f][j_f] = [x_f_ij, y_f_ij]
                T_f_ij = linear0(x_f_ij, y_f_ij)  # temperature
                T_fine[i_f][j_f] = T_f_ij
            else:
                x_fine[i_f][j_f] = [x_f_ij, y_f_ij]
                T_f_ij = linear0(x_f_ij, y_f_ij)  # temperature
                T_fine[i_f][j_f] = T_f_ij

                if y_f_ij > 0.5*L: # aluminium
                    rho_fine[i_f][j_f] = rho_Aluminium(T_f_ij)  # density
                    Cp_fine[i_f][j_f] = Cp_Aluminium(T_f_ij)  # specific heat capacity
                    k_fine[i_f][j_f] = k_Aluminium(T_f_ij)    # thermal conductivity
                    _lambda_fine[i_f][j_f] = k_fine[i_f][j_f] / (Cp_fine[i_f][j_f] * rho_fine[i_f][j_f])
                else: # Inconel800HT
                    rho_fine[i_f][j_f] = rho_Inconel800HT(T_f_ij)  # density
                    Cp_fine[i_f][j_f] = Cp_Inconel800HT(T_f_ij)    # specific heat capacity
                    k_fine[i_f][j_f] = k_Inconel800HT(T_f_ij)    # thermal conductivity
                    _lambda_fine[i_f][j_f] = k_fine[i_f][j_f] / (Cp_fine[i_f][j_f] * rho_fine[i_f][j_f])

            if i_f % 2 == 0 and j_f % 2 == 0:  # at the coarse point
                # getting the indexes for the coarse mesh
                i = int(i_f / 2)
                j = int(j_f / 2)
                xij = x_f_ij
                yij = y_f_ij

                if (xij < 0 or xij > L):
                    mask[i][j] = 0
                    x[i][j] = [xij, yij]
                elif (yij < 0 or yij > L):
                    mask[i][j] = 0
                    x[i][j] = [xij, yij]
                else:
                    x[i][j] = [xij, yij]
                    Tij = linear0(xij, yij)  # temperature
                    T[i][j] = Tij
                    if yij > 0.5*L: # aluminium
                        mask[i][j] = 1
                        rho[i][j] = rho_Aluminium(Tij)  # np.array([1.0 for i in range(len(T))]) # # density
                        Cp[i][j] = Cp_Aluminium(Tij)  # np.array([1.0 for i in range(len(T))]) #  # specific heat capacity
                        k[i][j] = k_Aluminium(Tij)  # np.array([1.0 for i in range(len(T))]) #  # thermal conductivity
                        _lambda[i][j] = k[i][j] / (Cp[i][j] * rho[i][j])
                    else: # Inconel800HT
                        mask[i][j] = 2
                        rho[i][j] = rho_Inconel800HT(Tij)  # np.array([1.0 for i in range(len(T))]) # # density
                        Cp[i][j] = Cp_Inconel800HT(Tij)  # np.array([1.0 for i in range(len(T))]) #  # specific heat capacity
                        k[i][j] = k_Inconel800HT(Tij)  # np.array([1.0 for i in range(len(T))]) #  # thermal conductivity
                        _lambda[i][j] = k[i][j] / (Cp[i][j] * rho[i][j])
    print(_lambda_fine)
    print(mask)
    return t, x, x_fine, T, T_fine, mask, np.array([rho, Cp, k, _lambda]), np.array([rho_fine, Cp_fine, k_fine, _lambda_fine])
