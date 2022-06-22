# Analytical heat transfer equation in 1D under the Dirichlet, Neumann and Mixed Boundary Conditions
# with the initial conitions of unit pulse function, the double sine waves and the linear function.
# importing libs:
print("importing libs")
import numpy as np

# setting parameters:
print("setting parameters")
L = 1
_lambda = 1
dx = 0.01
t_max = 0.01
dt = 0.002
x = np.arange(0,L+dx,dx) 
T = np.arange(0,L+dx,dx) 
t = np.arange(0,t_max+dt,dt)





# defining initial conditions:
print("defining initial conditions")

# unit pulse function:
def unitPulse(x): # input x is a np array
    T0 = np.arange(0,L+dx,dx) 
    m = 1 # magnitude of the unit pulse function
    for i in range(len(x)):
        xi = x[i]
        if 0.25<=xi<=0.75:
            T0[i] = m
        else:
            T0[i] = 0
    return T0 # output T0 is a np array

T0 = unitPulse(x)
#print(T0)

# TODO:
# sines function:


# linear function:




# other functions:
# this function is part (which is inside the integral) of the unitPulse_Dirichlet_T
def sin_n_pi_x_L(n,x):
    return np.sin(n*np.pi*x/L)







# tool box:
# integrating from start to end
def integral(f): 
    
    sum_fx = 0
    for i in range(len(f)):
        fi = f[i]
        sum_fx = sum_fx + fi*dx
    return sum_fx
    
    
    
    
    
    
# deriving analytical solutions:
print("deriving analytical solutions")
# unitPulse_Dirichlet_T
def unitPulse_Dirichlet_T(xi,ti): # dirichlet bc: T(0,t)=0; T(L,t)=0
    sum_n = 0
    N = 100
    for n in range(1,N):
        f = sin_n_pi_x_L(n,x) * unitPulse(x)
        #print(f)
        A_n = 2/L * integral(f)
        time_part = np.exp(-n**2 * np.pi**2 * _lambda * ti/(L**2))
        space_part = np.sin(n*np.pi*xi/L)
        T_n = A_n * time_part * space_part
        sum_n = sum_n + T_n
    return sum_n


# sines_Dirichlet_T


# linear_Dirichlet_T


# plotting:
print("plotting")


# plotting unitPulse_Dirichlet_T:
print("plotting unitPulse_Dirichlet_T")
plt.figure(figsize=(7,5))
plot_times = np.arange(0.0,t_max,dt)
color_list = ['k','r','b','g','y']
index = 0
for ti in plot_times:
    
    #plt.plot(y,V[int(t/dt),:],'Gray',label='numerical')
    for i in range(len(x)):
        xi = 0 + i*dx
        T[i] = unitPulse_Dirichlet_T(xi,ti)
    colori = 'o'+ color_list[index]
    plt.plot(x,T,colori,label='analytic at t={}s'.format(ti),markersize=3)
    #print(u)
    #if ti==dt:
    plt.legend(fontsize=12)
    index = index + 1
plt.xlabel('x (m)',fontsize=12)
plt.ylabel('T (k)',fontsize=12)
plt.title('Analytic Solution in 1D')
plt.savefig('unitPulse_Dirichlet_T.jpeg')  
plt.show()
print("finished plotting unitPulse_Dirichlet_T")

# TODO:
# plotting sines_Dirichlet_T:




