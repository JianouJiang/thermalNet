# Analytical heat transfer equation in 1D under the Dirichlet, Neumann and Mixed Boundary Conditions
# with the initial conitions of unit pulse function, the double sine waves and the linear function.
# importing libs:
print("importing libs")
import sys
sys.path.insert(0, '../tools/')
from tools import *
import numpy as np
from math import *
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import os


# setting parameters:
print("setting parameters")

parameters_directory="../heat_transfer/parameters.txt"
L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points, num_of_timeSteps_for_plotting = readParameters(parameters_directory)

x = np.arange(0,L+dx,dx) 
T = np.arange(0,L+dx,dx) 
t = np.arange(0,t_max+dt,dt)
dt_for_plotting = t_max / num_of_timeSteps_for_plotting
# https://en.wikipedia.org/wiki/Thermal_diffusivity, lambda = k/(cp*rho) with the unit m2/s




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



# other functions:
# this function is part (which is inside the integral) of the unitPulse_Dirichlet_T
def sin_n_pi_x_L(n,x):
    return np.sin(n*np.pi*x/L)

def cos_n_pi_x_L(n,x):
    return np.cos(n*np.pi*x/L)

# the following is for mixed B.C
func = lambda tau : np.tan(tau*L) + (tau)
mu = []
tau = np.linspace(0, 200, 201)
for n in range(1, 101):
    tau_initial_guess = (2*n-1)*np.pi/2
    tau_solution = fsolve(func, tau_initial_guess)
    mu.append(tau_solution)
    
# this is for a two-diffusivity (lambda) situation (e.g. steel at the left and copper at the right of the 1D rod)
def _lambda(x,Ti=1): # thermal diffusivity is related to space (e.g. steel for x<0.5 and copper for x>0.5) and temperature Ti, which is 1 by default
    interface_xi = 0.5
    _lambda_list = []
    
    for i in range(len(x)):
        xi = x[i]
        if xi<interface_xi:
            _lambda_list.append(_lambda1*Ti) # the lambda as a function of T is needed to add
        else:
            _lambda_list.append(_lambda2*Ti) # the lambda as a function of T is needed to add
    return _lambda_list

_lambda_list = _lambda(x)

   
    
# deriving analytical solutions:
print("deriving analytical solutions")
# unitPulse_Dirichlet_T
def unitPulse_Dirichlet_T(xi,ti, _lambda_i): # dirichlet bc: T(0,t)=0; T(L,t)=0
    sum_n = 0
    N = 100
    for n in range(1,N):
        f = sin_n_pi_x_L(n,x) * unitPulse(x)
        #print(f)
        A_n = 2/L * integral(x, f)
        time_part = np.exp(-n**2 * np.pi**2 * _lambda_i * ti/(L**2))
        space_part = np.sin(n*np.pi*xi/L)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n


# plotting:
print("plotting")


    
    
# plotting unitPulse_Dirichlet_T:

directory = "../data/exact_unitPulse_Dirichlet_T.txt"
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt_for_plotting)
print("plotting unitPulse_Dirichlet_T at t="+str(plot_times))
color_list = ['k','r','b','g','y']
index = 0
for ti in t:
    
    #plt.plot(y,V[int(t/dt),:],'Gray',label='numerical')
    colori = 'o'+ color_list[index]
    if ti == 0:
        T = unitPulse(x)
        plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3) # plot in dots
        plt.plot(x,T,'-k',markersize=3) # also plot in line
        plt.legend(fontsize=12)
        writeData(directory, ti, x, T, _lambda_list)
        index = index + 1
    else:
        
        for i in range(len(x)):
            xi = 0 + i*dx
            _lambda_i = _lambda_list[i]
            T[i] = unitPulse_Dirichlet_T(xi,ti,_lambda_i)
        writeData(directory, ti, x, T, _lambda_list)
        try:
            if ti==plot_times[index]:
                plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
                plt.legend(fontsize=12) 
                index = index + 1
        except IndexError:
            pass
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of Unit Pulse Function with Dirichlet B.C in 1D with $\lambda$={} m2/s'.format(_lambda))
plt.savefig('../img/unitPulse_Dirichlet_T.png')
plt.show()
print("finished plotting unitPulse_Dirichlet_T")


# sines_Dirichlet_T



# deriving analytical solutions:
print("deriving analytical solutions")
# sines_Dirichlet_T
def sines_Dirichlet_T(xi,ti,_lambda_i): # dirichlet bc: T(0,t)=0; T(L,t)=0
    sum_n = 0
    N = 100
    for n in range(1,N):
        f = sin_n_pi_x_L(n,x) * sines(x)
        #print(f)
        A_n = 2/L * integral(x, f)
        time_part = np.exp(-n**2 * np.pi**2 * _lambda_i * ti/(L**2))
        space_part = np.sin(n*np.pi*xi/L)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n



# plotting sines_Dirichlet_T:

directory = "../data/exact_sines_Dirichlet_T.txt"
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt_for_plotting)
print("plotting sines_Dirichlet_T at t="+str(plot_times))
color_list = ['k','r','b','g','y']
index = 0
for ti in plot_times:
    
    #plt.plot(y,V[int(t/dt),:],'Gray',label='numerical')
    for i in range(len(x)):
        xi = 0 + i*dx
        _lambda_i = _lambda_list[i]
        T[i] = sines_Dirichlet_T(xi,ti,_lambda_i)
    colori = 'o'+ color_list[index]
    if ti == 0.0:
        plt.plot(x,sines(x),colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.plot(x,sines(x),'-k',markersize=3) # also plot in line
        plt.legend(fontsize=12)
        writeData(directory, ti, x, T, _lambda_list)
        index =index+1
    else:
                
        for i in range(len(x)):
            xi = 0 + i*dx
            _lambda_i = _lambda_list[i]
            T[i] = sines_Dirichlet_T(xi,ti,_lambda_i)
        writeData(directory, ti, x, T, _lambda_list)
        
        try:
            if ti==plot_times[index]:
                plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
                plt.legend(fontsize=12) 
                index = index + 1
        except IndexError:
            pass
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of sines function with Dirichlet B.C in 1D with $\lambda$={} m2/s'.format(_lambda))
plt.savefig('../img/sines_Dirichlet_T.png')
plt.show()
print("finished plotting sines_Dirichlet_T")

# linear_Dirichlet_T




# deriving analytical solutions:
print("deriving analytical solutions")
# linear_Dirichlet_T
def linear_Dirichlet_T(xi,ti,_lambda_i): # dirichlet bc: T(0,t)=0; T(L,t)=0
    sum_n = 0
    N = 100
    for n in range(1,N):
        f = sin_n_pi_x_L(n,x) * linear(x)
        #print(f)
        A_n = 2/L * integral(x,f)
        time_part = np.exp(-n**2 * np.pi**2 * _lambda_i * ti/(L**2))
        space_part = np.sin(n*np.pi*xi/L)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n



# plotting linear_Dirichlet_T:

directory = "../data/exact_linear_Dirichlet_T.txt"
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt_for_plotting)
print("plotting linear_Dirichlet_T at t="+str(plot_times))
color_list = ['k','r','b','g','y']
index = 0
for ti in plot_times:
    
    #plt.plot(y,V[int(t/dt),:],'Gray',label='numerical')
    
    colori = 'o'+ color_list[index]
    if ti == 0.0:
        plt.plot(x,linear(x),colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.plot(x,linear(x),'-k',markersize=3) # also plot in line
        plt.legend(fontsize=12)
        writeData(directory, ti, x, T, _lambda_list)
        index = index + 1
    else:
        
        for i in range(len(x)):
            xi = 0 + i*dx
            _lambda_i = _lambda_list[i]
            T[i] = linear_Dirichlet_T(xi,ti,_lambda_i)
        writeData(directory, ti, x, T, _lambda_list)
        
        try:
            if ti==plot_times[index]:
                plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
                plt.legend(fontsize=12) 
                index = index + 1
        except IndexError:
            pass
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of linear Function with Dirichlet B.C in 1D with $\lambda$={} m2/s'.format(_lambda))
plt.savefig('../img/linear_Dirichlet_T.png')
plt.show()
print("finished plotting linear_Dirichlet_T")


# Neumann conditions:

# unitPulse_Neumann_T:

# deriving analytical solutions:
print("deriving analytical solutions")
# unitPulse_Neumann_T
def unitPulse_Neumann_T(xi,ti,_lambda_i): # Neumann bc: T_x(0,t)=0; T_x(L,t)=0
    sum_n = 0
    N = 100
    A0 = m/2
    sum_n = sum_n + A0 
    for n in range(1,N):
        f = cos_n_pi_x_L(n,x) * unitPulse(x)
        #print(f)
        A_n = 2/L * integral(x, f)
        time_part = np.exp(-n**2 * np.pi**2 * _lambda_i * ti/(L**2))
        space_part = np.cos(n*np.pi*xi/L)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n



# plotting unitPulse_Neumann_T:
print("plotting unitPulse_Neumann_T")
directory = "../data/exact_unitPulse_Neumann_T.txt"
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt_for_plotting)
print("plotting unitPulse_Neumann_T at t="+str(plot_times))
color_list = ['k','r','b','g','y']
index = 0
for ti in plot_times:
    
    #plt.plot(y,V[int(t/dt),:],'Gray',label='numerical')
    for i in range(len(x)):
        xi = 0 + i*dx
        _lambda_i = _lambda_list[i]
        T[i] = unitPulse_Neumann_T(xi,ti,_lambda_i)
    colori = 'o'+ color_list[index]
    if ti == 0.0:
        plt.plot(x,unitPulse(x),colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.plot(x,unitPulse(x),'-k',markersize=3) # also plot in line
        plt.legend(fontsize=12)
        writeData(directory, ti, x, T, _lambda_list)
        index = index + 1
    else:
        
        for i in range(len(x)):
            xi = 0 + i*dx
            _lambda_i = _lambda_list[i]
            T[i] = unitPulse_Neumann_T(xi,ti,_lambda_i)
        writeData(directory, ti, x, T, _lambda_list)
        
        try:
            if ti==plot_times[index]:
                plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
                plt.legend(fontsize=12) 
                index = index + 1
        except IndexError:
            pass
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of Unit Pulse Function with Neumann B.C in 1D with $\lambda$={} m2/s'.format(_lambda))
plt.savefig('../img/unitPulse_Neumann_T.png')
plt.show()
print("finished plotting unitPulse_Neumann_T")


# sines_Neumann_T:

# deriving analytical solutions:
print("deriving analytical solutions")
# sines_Neumann_T
def sines_Neumann_T(xi,ti,_lambda_i): # Neumann bc: T_x(0,t)=0; T_x(L,t)=0
    sum_n = 0
    N = 100
    A0 = 2*m/np.pi
    sum_n = sum_n + A0 
    for n in range(1,N):
        f = cos_n_pi_x_L(n,x) * sines(x)
        #print(f)
        A_n = 2/L * integral(x, f)
        time_part = np.exp(-n**2 * np.pi**2 * _lambda_i * ti/(L**2))
        space_part = np.cos(n*np.pi*xi/L)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n


# plotting sines_Neumann_T:

directory = "../data/exact_sines_Neumann_T.txt"
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt_for_plotting)
print("plotting sines_Neumann_T at t="+str(plot_times))
color_list = ['k','r','b','g','y']
index = 0
for ti in plot_times:
    
    #plt.plot(y,V[int(t/dt),:],'Gray',label='numerical')
    for i in range(len(x)):
        xi = 0 + i*dx
        _lambda_i = _lambda_list[i]
        T[i] = sines_Neumann_T(xi,ti,_lambda_i)
    colori = 'o'+ color_list[index]
    if ti == 0.0:
        plt.plot(x,sines(x),colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.plot(x,sines(x),'-k',markersize=3) # also plot in line
        plt.legend(fontsize=12)
        writeData(directory, ti, x, T, _lambda_list)
        index = index + 1
    else:
        
        for i in range(len(x)):
            xi = 0 + i*dx
            _lambda_i = _lambda_list[i]
            T[i] = sines_Neumann_T(xi,ti,_lambda_i)
        writeData(directory, ti, x, T, _lambda_list)
        
        try:
            if ti==plot_times[index]:
                plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
                plt.legend(fontsize=12) 
                index = index + 1
        except IndexError:
            pass
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of sine Function with Neumann B.C in 1D with $\lambda$={} m2/s'.format(_lambda))
plt.savefig('../img/sines_Neumann_T.png')
plt.show()
print("finished plotting sines_Neumann_T")


# linear_Neumann_T


# deriving analytical solutions:
print("deriving analytical solutions")
# linear_Neumann_T
def linear_Neumann_T(xi,ti,_lambda_i): # Neumann bc: T_x(0,t)=0; T_x(L,t)=0
    sum_n = 0
    N = 100
    A0 = L/2
    sum_n = sum_n + A0 
    for n in range(1,N):
        f = cos_n_pi_x_L(n,x) * linear(x)
        #print(f)
        A_n = 2/L * integral(x,f)
        time_part = np.exp(-n**2 * np.pi**2 * _lambda_i * ti/(L**2))
        space_part = np.cos(n*np.pi*xi/L)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n


# plotting linear_Neumann_T:

directory = "../data/exact_linear_Neumann_T.txt"
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt_for_plotting)
print("plotting linear_Neumann_T at t="+str(plot_times))
color_list = ['k','r','b','g','y']
index = 0
for ti in plot_times:
    
    #plt.plot(y,V[int(t/dt),:],'Gray',label='numerical')
    
    colori = 'o'+ color_list[index]
    if ti == 0.0:
        plt.plot(x,linear(x),colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.plot(x,linear(x),'-k',markersize=3) # also plot in line
        plt.legend(fontsize=12)
        writeData(directory, ti, x, T, _lambda_list)
        index = index + 1
    else:
        
        for i in range(len(x)):
            xi = 0 + i*dx
            _lambda_i = _lambda_list[i]
            T[i] = linear_Neumann_T(xi,ti,_lambda_i)
        writeData(directory, ti, x, T, _lambda_list)
        
        try:
            if ti==plot_times[index]:
                plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
                plt.legend(fontsize=12) 
                index = index + 1
        except IndexError:
            pass
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of linear Function with Neumann B.C in 1D with $\lambda$={} m2/s'.format(_lambda))
plt.savefig('../img/linear_Neumann_T.png')
plt.show()
print("finished plotting linear_Neumann_T")


# Mixed conditions:


# unitPulse_Mixed_T:
# deriving analytical solutions:
print("deriving analytical solutions")
# unitPulse_Mixed_T
def unitPulse_Mixed_T(xi,ti,_lambda_i): # dirichlet bc: T(0,t)=0; T(L,t)=0
    sum_n = 0
    N = 100
    for n in range(1,N):
        f1 = np.sin(mu[n-1]*x) * unitPulse(x)
        f2 = (np.sin(mu[n-1]*x))**2
        A_n = integral(x,f1) / integral(x,f2)
        time_part = np.exp((-1)*_lambda_i* ((mu[n-1])**2) * ti)
        space_part = np.sin(mu[n-1]*xi)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n



# plotting unitPulse_Mixed_T:

directory = "../data/exact_unitPulse_Mixed_T.txt"
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt_for_plotting)
print("plotting unitPulse_Mixed_T at t="+str(plot_times))
color_list = ['k','r','b','g','y']
index = 0
for ti in plot_times:
    
    #plt.plot(y,V[int(t/dt),:],'Gray',label='numerical')
    for i in range(len(x)):
        xi = 0 + i*dx
        _lambda_i = _lambda_list[i]
        T[i] = unitPulse_Mixed_T(xi,ti,_lambda_i)
    colori = 'o'+ color_list[index]
    if ti == 0.0:
        plt.plot(x,unitPulse(x),colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.plot(x,unitPulse(x),'-k',markersize=3) # also plot in line
        plt.legend(fontsize=12)
        writeData(directory, ti, x, T, _lambda_list)
        index = index + 1
    else:
        
        for i in range(len(x)):
            xi = 0 + i*dx
            _lambda_i = _lambda_list[i]
            T[i] = unitPulse_Mixed_T(xi,ti,_lambda_i)
        writeData(directory, ti, x, T, _lambda_list)
        
        try:
            if ti==plot_times[index]:
                plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
                plt.legend(fontsize=12) 
                index = index + 1
        except IndexError:
            pass
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of Unit Pulse Function with Mixed B.C in 1D with $\lambda$={} m2/s'.format(_lambda))
plt.savefig('../img/unitPulse_Mixed_T.png')
plt.show()
print("finished plotting unitPulse_Mixed_T")


# Sines_Mixed_T:

# deriving analytical solutions:
print("deriving analytical solutions")
# sines_Mixed_T
def sines_Mixed_T(xi,ti,_lambda_i): # dirichlet bc: T(0,t)=0; T(L,t)=0
    sum_n = 0
    N = 100
    for n in range(1,N):
        f1 = np.sin(mu[n-1]*x) * sines(x)
        f2 = (np.sin(mu[n-1]*x))**2
        A_n = integral(x,f1) / integral(x,f2)
        time_part = np.exp((-1)*_lambda_i* ((mu[n-1])**2) * ti)
        space_part = np.sin(mu[n-1]*xi)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n



# plotting sines_Mixed_T:

directory = "../data/exact_sines_Mixed_T.txt"
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt_for_plotting)
print("plotting sines_Mixed_T at t="+str(plot_times))
color_list = ['k','r','b','g','y']
index = 0
for ti in plot_times:
    
    #plt.plot(y,V[int(t/dt),:],'Gray',label='numerical')
    for i in range(len(x)):
        xi = 0 + i*dx
        _lambda_i = _lambda_list[i]
        T[i] = sines_Mixed_T(xi,ti,_lambda_i)
    colori = 'o'+ color_list[index]
    if ti == 0.0:
        plt.plot(x,sines(x),colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.plot(x,sines(x),'-k',markersize=3) # also plot in line
        plt.legend(fontsize=12)
        writeData(directory, ti, x, T, _lambda_list)
        index = index + 1
    else:
        
        for i in range(len(x)):
            xi = 0 + i*dx
            _lambda_i = _lambda_list[i]
            T[i] = sines_Mixed_T(xi,ti,_lambda_i)
        writeData(directory, ti, x, T, _lambda_list)
        
        try:
            if ti==plot_times[index]:
                plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
                plt.legend(fontsize=12) 
                index = index + 1
        except IndexError:
            pass
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of sines Function with Mixed B.C in 1D with $\lambda$={} m2/s'.format(_lambda))
plt.savefig('../img/sines_Mixed_T.png')
plt.show()
print("finished plotting sines_Mixed_T")



# Linear_Mixed_T:

# deriving analytical solutions:
print("deriving analytical solutions")
# linear_Mixed_T
def linear_Mixed_T(xi,ti,_lambda_i): # dirichlet bc: T(0,t)=0; T(L,t)=0
    sum_n = 0
    N = 100
    for n in range(1,N):
        f1 = np.sin(mu[n-1]*x) * linear(x)
        f2 = (np.sin(mu[n-1]*x))**2
        A_n = integral(x,f1) / integral(x,f2)
        time_part = np.exp((-1)*_lambda_i* ((mu[n-1])**2) * ti)
        space_part = np.sin(mu[n-1]*xi)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n


# plotting linear_Mixed_T:

directory = "../data/exact_linear_Mixed_T.txt"
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt_for_plotting)
print("plotting linear_Mixed_T at t="+str(plot_times))
color_list = ['k','r','b','g','y']
index = 0
for ti in plot_times:
    
    #plt.plot(y,V[int(t/dt),:],'Gray',label='numerical')
    
    colori = 'o'+ color_list[index]
    if ti == 0.0:
        plt.plot(x,linear(x),colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.plot(x,linear(x),'-k',markersize=3) # also plot in line
        plt.legend(fontsize=12)
        writeData(directory, ti, x, T, _lambda_list)
        index = index + 1
    else:
        
        for i in range(len(x)):
            xi = 0 + i*dx
            _lambda_i = _lambda_list[i]
            T[i] = linear_Mixed_T(xi,ti,_lambda_i)
        writeData(directory, ti, x, T, _lambda_list)
        
        try:
            if ti==plot_times[index]:
                plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
                plt.legend(fontsize=12) 
                index = index + 1
        except IndexError:
            pass
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of linear Function with Mixed B.C in 1D with $\lambda$={} m2/s'.format(_lambda))
plt.savefig('../img/linear_Mixed_T.png')
plt.show()
print("finished plotting linear_Mixed_T")





# TODO: @zhengyang
# two layer heat transfer

# deriving analytical solutions:

# Dirichlet boundary condition:
import numpy as np
import matplotlib.pyplot as plt
k1 = 1
k2 = 1
d1 = 1
d2 = 1
l1 = 1
l2 = 1
rho1 = 1
rho2 = 1
c1 = 1
c2 = 1
D1 = 1
D2 = 1
H = 0.5
a = 1
b = 0
d = 0
x_0 = 0
x_1 = 0.5
x_2 = 1

# initial conditions:
def f(i, x):
    if i ==1:
        return 0
    else:
        return 0

# steady state:
def w(i, x):
    if i==1:
        return 1 - ((k2*H*(a-d)*x)/(b*k1*H + a*k2*l1*H + a*k1*l2*H + a*k1*k2))
    else:
        return 1 - (((a - d)*(k1*H*(x-x_1) + k2*(l1*H + k1)))/(b*k1*H + a*k2*l1*H + a*k1*l2*H + a*k1*k2))


# find \mu_m:
def J2m(mu_m):
    return k1/k2*d2/d1*np.cos(mu_m*l1/d1)

def K2m(mu_m):
    return np.sin(mu_m*l1/d1) + k1/d1 * mu_m/H*np.cos(mu_m*l1/d1)

def h(mu_m):
    return J2m(mu_m)*(a*np.sin(mu_m*l2/d2) + mu_m*b/d2 * np.cos(mu_m*l2/d2)) + K2m(mu_m)*(-mu_m*b/d2*np.sin(mu_m*l2/d2) + a*np.cos(mu_m*l2/d2))


mu_start = 0
mu_end = 1000000
dmu = 0.1
    
def findMusWhenFunctionEqualZero(mu_start, mu_end, dmu):   
    mu_m = np.arange(mu_start,mu_end,dmu)
    mu_solution = []
    func = h(mu_m)
    for i in range((len(func)-1)):
        mu_i = mu_m[i]
        fi = func[i]
        fip1 = func[i+1]
        if fi*fip1<0:
            # we are at the cross point, get y=ax+b based on the two points, then get x when y=0
            _a = (fip1 - fi) / dmu
            _b = fi - _a*mu_i
            mu_sol_i = - _b/_a
            mu_solution.append(mu_sol_i)
            
    print("number of mus: "+str(len(mu_solution)))        
    return mu_solution

mu_solution=findMusWhenFunctionEqualZero(mu_start,mu_end,dmu)


mu_m = np.arange(mu_start,mu_end,dmu)
print(mu_solution)

func = h(mu_m)
plt.plot(mu_m, func)

plt.xlabel('$\mu$')
plt.show()

# space function X_i_m(x):
def J_1(m):
    return 1

def J_2(m):
    return k1/k2*d2/d1*np.cos(mu_solution[m-1]*l1/d1)

def K_1(m):
    return 0

def K_2(m):
    return np.sin(mu_solution[m-1]*l1/d1) + k1/d1 * mu_solution[m-1]/H*np.cos(mu_solution[m-1]*l1/d1)

def X_1(m, x):
    return J_1(m)*np.sin(mu_solution[m-1]*(x- x_0)) + K_1(m)*np.cos(mu_solution[m-1]*(x- x_0))

def X_2(m, x):
    return J_2(m)*np.sin(mu_solution[m-1]*(x- x_0)) + K_2(m)*np.cos(mu_solution[m-1]*(x- x_0))

# coefficients C_m:
def N_1(x):
    return (f(1, x) - w(1, x))*X_1(m, x)

def N_2(x):
    return (f(2, x) - w(2, x))*X_2(m, x)

def P_1(x):
    return (X_1(m, x))**2

def P_2(x):
    return (X_2(m, x))**2

def C(m):
    numerator = rho1*c1*integral(x, N_1, x_0, x_1) + rho2*c2*integral(x, N_2, x_1, x_2)
    denominator = rho1*c1*integral(x, P_1, x_0, x_1) + rho1*c1*integral(x, P_2, x_1, x_2)
    return numerator/denominator


# deriving analytical solutions:
print("deriving analytical solutions")
# two layer example:
def twolayer_transient_T(xi,ti):
    sum_n = 0
    N = 80
    if xi<=x_1:
        for n in range(1,N):
            C_n = C(n)
            time_part = np.exp((-1)*((mu[n-1])**2) * ti)
            space_part = X_1(n, xi)
            T_n = C_n * time_part * space_part
            sum_n = sum_n + T_n
        return sum_n
    else:
        for n in range(1,N):
            C_n = C(n)
            time_part = np.exp((-1)*((mu[n-1])**2) * ti)
            space_part = X_2(n, xi)
            T_n = C_n * time_part * space_part
            sum_n = sum_n + T_n
        return sum_n

# function of temperature:
def U(xi, ti):
    if xi<=x_1:
        return w(1, xi) + twolayer_transient_T(xi, ti)
    else:
        return w(2, xi) + twolayer_transient_T(xi, ti)


# plotting:
print("plotting")

def f_initial(x):
    return 0

# plotting U:
print("plotting U")
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt_for_plotting)
color_list = ['k','r','b','g','y']
index = 0
for ti in plot_times:
    
    #plt.plot(y,V[int(t/dt),:],'Gray',label='numerical')
    
    colori = 'o'+ color_list[index]
    if ti == 0.0:
        plt.plot(x,f_initial(x),colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.plot(x,f_initial(x),'-k',markersize=3) # also plot in line
        plt.legend(fontsize=12)
        writeData(directory, ti, x, T, _lambda_list)
    else:
        for i in range(len(x)):
            xi = 0 + i*dx
            _lambda_i = _lambda_list[i]
            T[i] = U(xi,ti)
        plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.legend(fontsize=12)
        writeData(directory, ti, x, T, _lambda_list)
    index = index + 1
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of two layer example')
plt.savefig('../img/two_layer_example.png')
plt.show()
print("finished plotting two layer example")


# plotting linear_Mixed_T:
print("plotting twoLayer_linear_Dirichlet_T")
directory = "data/exact_twoLayer_linear_Dirichlet_T.txt"






print("finished plotting twoLayer_linear_Dirichlet_T")
