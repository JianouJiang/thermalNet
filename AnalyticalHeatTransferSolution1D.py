# Analytical heat transfer equation in 1D under the Dirichlet, Neumann and Mixed Boundary Conditions
# with the initial conitions of unit pulse function, the double sine waves and the linear function.
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
# https://en.wikipedia.org/wiki/Thermal_diffusivity, lambda = k/(cp*rho) with the unit m2/s
_lambda1 = 1.5
_lambda2 = 0.5
_lambda_list = _lambda(x)
dx = 0.01
t_max = 0.01
dt = 0.002
x = np.arange(0,L+dx,dx) 
T = np.arange(0,L+dx,dx) 
t = np.arange(0,t_max+dt,dt)
func = lambda tau : np.tan(tau*L) + (tau)
mu = []
tau = np.linspace(0, 200, 201)
for n in range(1, 101):
    tau_initial_guess = (2*n-1)*np.pi/2
    tau_solution = fsolve(func, tau_initial_guess)
    mu.append(tau_solution)





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

# TODO: @zhengyang
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

# Liner function:
def linear(x): # input x is a np array
    T0 = np.arange(0,L+dx,dx) 
    for i in range(len(x)):
        xi = x[i]
        T0[i] = xi
    return T0 # output T0 is a np array    



# other functions:
# this function is part (which is inside the integral) of the unitPulse_Dirichlet_T
def sin_n_pi_x_L(n,x):
    return np.sin(n*np.pi*x/L)

def cos_n_pi_x_L(n,x):
    return np.cos(n*np.pi*x/L)


def _lambda(xi, Ti=1): # thermal diffusivity is related to space (e.g. steel for x<0.5 and copper for x>0.5) and temperature Ti, which is 1 by default
    interface_xi = 0.5
    if xi<interface_xi:
        return _lambda1*Ti # the lambda as a function of T is needed to add
    else:
        return _lambda2*Ti # the lambda as a function of T is needed to add


# tool box:
# integrating from start to end
def integral(f): 
    
    sum_fx = 0
    for i in range(0,(len(f)-1)):
        fi = (f[i]+f[i+1] )/2 # central differencing
        sum_fx = sum_fx + fi*dx
    return sum_fx
    
# open file in write mode and write data
def writeData(directory, ti, T, _lambda): # note T and lambda here are lists, whereas time t is single value
    isFile = os.path.isfile(directory) 
    with open(directory, 'a') as fp:
        index = 0
        for Ti in T:
            xi = 0 + dx * index
            _lambda_i = _lambda[index]
            # write each ti, xi, Ti, lambda on a new line
            line = str(ti) + " " + str(xi) + " " + str(Ti) + " " + str(_lambda_i) + "\n"
            fp.write(line)
            index = index + 1
    return
    
    
    
# deriving analytical solutions:
print("deriving analytical solutions")
# unitPulse_Dirichlet_T
def unitPulse_Dirichlet_T(xi,ti, _lambda_i): # dirichlet bc: T(0,t)=0; T(L,t)=0
    sum_n = 0
    N = 100
    for n in range(1,N):
        f = sin_n_pi_x_L(n,x) * unitPulse(x)
        #print(f)
        A_n = 2/L * integral(f)
        time_part = np.exp(-n**2 * np.pi**2 * _lambda_i * ti/(L**2))
        space_part = np.sin(n*np.pi*xi/L)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n


# plotting:
print("plotting")


    
    
# plotting unitPulse_Dirichlet_T:
print("plotting unitPulse_Dirichlet_T")
directory = "data/unitPulse_Dirichlet_T.txt"
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt)
color_list = ['k','r','b','g','y']
index = 0
for ti in plot_times:
    
    #plt.plot(y,V[int(t/dt),:],'Gray',label='numerical')
    colori = 'o'+ color_list[index]
    if ti == 0.0:
        T = unitPulse(x)
        plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3) # plot in dots
        plt.plot(x,T,'-k',markersize=3) # also plot in line
        plt.legend(fontsize=12)
        writeData(directory, ti, T, _lambda_list)
        
    else:
        
        for i in range(len(x)):
            xi = 0 + i*dx
            _lambda_i = _lambda_list[i]
            T[i] = unitPulse_Dirichlet_T(xi,ti,_lambda_i)
        plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.legend(fontsize=12) 
        writeData(directory, ti, T, _lambda_list)
    index = index + 1
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of Unit Pulse Function with Dirichlet B.C in 1D with $\lambda$={} m2/s'.format(_lambda))
plt.savefig('img/unitPulse_Dirichlet_T.png')  
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
        A_n = 2/L * integral(f)
        time_part = np.exp(-n**2 * np.pi**2 * _lambda_i * ti/(L**2))
        space_part = np.sin(n*np.pi*xi/L)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n


# plotting:
print("plotting")


# plotting sines_Dirichlet_T:
print("plotting sines_Dirichlet_T")
directory = "data/sines_Dirichlet_T.txt"
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt)
color_list = ['k','r','b','g','y']
index = 0
for ti in plot_times:
    
    #plt.plot(y,V[int(t/dt),:],'Gray',label='numerical')
    for i in range(len(x)):
        xi = 0 + i*dx
        T[i] = sines_Dirichlet_T(xi,ti)
    colori = 'o'+ color_list[index]
    if ti == 0.0:
        plt.plot(x,sines(x),colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.plot(x,sines(x),'-k',markersize=3) # also plot in line
        plt.legend(fontsize=12)
        writeData(directory, ti, T, _lambda_list)
    else:
        plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.legend(fontsize=12)
        _lambda_i = _lambda_list[i]
        writeData(directory, ti, T, _lambda_list)
    index = index + 1
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of sines function with Dirichlet B.C in 1D with $\lambda$={} m2/s'.format(_lambda))
plt.savefig('img/sines_Dirichlet_T.png')  
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
        A_n = 2/L * integral(f)
        time_part = np.exp(-n**2 * np.pi**2 * _lambda_i * ti/(L**2))
        space_part = np.sin(n*np.pi*xi/L)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n


# plotting:
print("plotting")


# plotting linear_Dirichlet_T:
print("plotting linear_Dirichlet_T")
directory = "data/linear_Dirichlet_T.txt"
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt)
color_list = ['k','r','b','g','y']
index = 0
for ti in plot_times:
    
    #plt.plot(y,V[int(t/dt),:],'Gray',label='numerical')
    for i in range(len(x)):
        xi = 0 + i*dx
        _lambda_i = _lambda_list[i]
        T[i] = linear_Dirichlet_T(xi,ti,_lambda_i)
    colori = 'o'+ color_list[index]
    if ti == 0.0:
        plt.plot(x,linear(x),colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.plot(x,linear(x),'-k',markersize=3) # also plot in line
        plt.legend(fontsize=12)
        writeData(directory, ti, T, _lambda_list)
    else:
        plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.legend(fontsize=12)
        writeData(directory, ti, T, _lambda_list)
    index = index + 1
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of linear Function with Dirichlet B.C in 1D with $\lambda$={} m2/s'.format(_lambda))
plt.savefig('img/linear_Dirichlet_T.png')  
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
        A_n = 2/L * integral(f)
        time_part = np.exp(-n**2 * np.pi**2 * _lambda_i * ti/(L**2))
        space_part = np.cos(n*np.pi*xi/L)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n


# plotting:
print("plotting")


# plotting unitPulse_Neumann_T:
print("plotting unitPulse_Neumann_T")
directory = "data/unitPulse_Neumann_T.txt"
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt)
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
        writeData(directory, ti, T, _lambda_list)
    else:
        plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.legend(fontsize=12)
        _lambda_i = _lambda
        writeData(directory, ti, T, _lambda_list)
    index = index + 1
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of Unit Pulse Function with Neumann B.C in 1D with $\lambda$={} m2/s'.format(_lambda))
plt.savefig('img/unitPulse_Neumann_T.png')  
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
        A_n = 2/L * integral(f)
        time_part = np.exp(-n**2 * np.pi**2 * _lambda_i * ti/(L**2))
        space_part = np.cos(n*np.pi*xi/L)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n


# plotting:
print("plotting")


# plotting sines_Neumann_T:
print("plotting sines_Neumann_T")
directory = "data/sines_Neumann_T.txt"
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt)
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
        writeData(directory, ti, T, _lambda_list)
    else:
        plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.legend(fontsize=12)
        writeData(directory, ti, T, _lambda_list)
    index = index + 1
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of sine Function with Neumann B.C in 1D with $\lambda$={} m2/s'.format(_lambda))
plt.savefig('img/sines_Neumann_T.png')  
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
        A_n = 2/L * integral(f)
        time_part = np.exp(-n**2 * np.pi**2 * _lambda_i * ti/(L**2))
        space_part = np.cos(n*np.pi*xi/L)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n


# plotting:
print("plotting")


# plotting linear_Neumann_T:
print("plotting linear_Neumann_T")
directory = "data/linear_Neumann_T.txt"
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt)
color_list = ['k','r','b','g','y']
index = 0
for ti in plot_times:
    
    #plt.plot(y,V[int(t/dt),:],'Gray',label='numerical')
    for i in range(len(x)):
        xi = 0 + i*dx
        _lambda_i = _lambda_list[i]
        T[i] = linear_Neumann_T(xi,ti,_lambda_i)
    colori = 'o'+ color_list[index]
    if ti == 0.0:
        plt.plot(x,linear(x),colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.plot(x,linear(x),'-k',markersize=3) # also plot in line
        plt.legend(fontsize=12)
        writeData(directory, ti, T, _lambda_list)
    else:
        plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.legend(fontsize=12)
        writeData(directory, ti, T, _lambda_list)
    index = index + 1
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of linear Function with Neumann B.C in 1D with $\lambda$={} m2/s'.format(_lambda))
plt.savefig('img/linear_Neumann_T.png')  
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
        A_n = integral(f1) / integral(f2)
        time_part = np.exp((-1)*_lambda_i* ((mu[n-1])**2) * ti)
        space_part = np.sin(mu[n-1]*xi)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n


# plotting:
print("plotting")


# plotting unitPulse_Mixed_T:
print("plotting unitPulse_Mixed_T")
directory = "data/unitPulse_Mixed_T.txt"
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt)
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
        writeData(directory, ti, T, _lambda_list)
    else:
        plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.legend(fontsize=12)
        writeData(directory, ti, T, _lambda_list)
    index = index + 1
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of Unit Pulse Function with Mixed B.C in 1D with $\lambda$={} m2/s'.format(_lambda))
plt.savefig('img/unitPulse_Mixed_T.png')  
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
        A_n = integral(f1) / integral(f2)
        time_part = np.exp((-1)*_lambda_i* ((mu[n-1])**2) * ti)
        space_part = np.sin(mu[n-1]*xi)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n


# plotting:
print("plotting")


# plotting sines_Mixed_T:
print("plotting sines_Mixed_T")
directory = "data/sines_Mixed_T.txt"
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt)
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
        _lambda_i = _lambda
        writeData(directory, ti, T, _lambda_list)
    else:
        plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.legend(fontsize=12)
        _lambda_i = _lambda
        writeData(directory, ti, T, _lambda_list)
    index = index + 1
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of sines Function with Mixed B.C in 1D with $\lambda$={} m2/s'.format(_lambda))
plt.savefig('img/sines_Mixed_T.png')  
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
        A_n = integral(f1) / integral(f2)
        time_part = np.exp((-1)*_lambda_i* ((mu[n-1])**2) * ti)
        space_part = np.sin(mu[n-1]*xi)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n


# plotting:
print("plotting")


# plotting linear_Mixed_T:
print("plotting linear_Mixed_T")
directory = "data/linear_Mixed_T.txt"
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt)
color_list = ['k','r','b','g','y']
index = 0
for ti in plot_times:
    
    #plt.plot(y,V[int(t/dt),:],'Gray',label='numerical')
    for i in range(len(x)):
        xi = 0 + i*dx
        _lambda_i = _lambda_list[i]
        T[i] = linear_Mixed_T(xi,ti,_lambda_i)
    colori = 'o'+ color_list[index]
    if ti == 0.0:
        plt.plot(x,linear(x),colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.plot(x,linear(x),'-k',markersize=3) # also plot in line
        plt.legend(fontsize=12)
        _lambda_i = _lambda
        writeData(directory, ti, T, _lambda_list)
    else:
        plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
        plt.legend(fontsize=12)
        _lambda_i = _lambda
        writeData(directory, ti, T, _lambda_list)
    index = index + 1
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution of linear Function with Mixed B.C in 1D with $\lambda$={} m2/s'.format(_lambda))
plt.savefig('img/linear_Mixed_T.png')  
plt.show()
print("finished plotting linear_Mixed_T")
