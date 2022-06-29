# different solvers for system of equations
#importing libraries
import numpy as np
from numpy import array, zeros, diag, diagflat, dot
import matplotlib.pyplot as plt


# defining linear system solver functions


# xiao yang: Jacobi
# x0 = our initial guess, N = number of iterations, r = tolerance

def Jacobi(A, b,x= None, N=100, r=10**-6): #Ax=b,we are finding x, N = number of iterations, r is the residual
    residual_list = []
    if (x is None) :
        x=zeros(len(A[0]))
    D= diag(A)
    R= A-diagflat(D)
    residual =1000
    xn_minus1=x
    #iterate for N times
    for Ni in range(N):
        
        x= (b-dot(R,x))/D
        difference= xn_minus1-x
        magnitude= np.linalg.norm(difference)
        residual =magnitude/len(x) 
        residual_list.append(residual)
        if residual <r:
            print("Jacobi: The number of iterations is: {}".format(Ni))
            break
        xn_minus1=x
        
    return x, residual_list




# chen nuo: Gauss-Seidel method
# import libraries
import numpy as np
from scipy.linalg import solve

def GaussSeidel(A,b,x=None,N=100, r=10**-6):
    residual_list = []
    if (x is None) :
        x=zeros(len(A[0]))
    xn_minus1 = x
    L = np.tril(A) # define the lower triangular matrix (L) for A
    U = A - L # define the upper triangular matrix (U)
    residual = 1000
    for i in range (N): # create a "for" loop
        # define x as the dot product of inverse L and B mins the dot product of U and x
        x = np.dot(np.linalg.inv(L),b-np.dot(U,x))
        difference = xn_minus1 - x
        magnitude = np.linalg.norm(difference)
        residual = magnitude / len(x)
        residual_list.append(residual)
        if residual < r:
            print("GaussSeidel: The number of iterations is: {}".format(Ni))
            break
    return x, residual_list





# kaiyi: Successive over-relaxation (SOR)
# x0 = our initial guess, N = number of iterations, r = tolerance
# w = relaxation factor, 1<w<2. If w=1, it's same as Gauss-Seidel Method

def SOR(A, b, x0=None, N=100, r=10**-6, w=1.5): 
    residual_list = []
    n = b.shape
    if (x0 is None):
        x0=zeros(len(A[0]))
    x = x0
    for Ni in range (1, N): 
        for i in range(n[0]): 
            new_values_sum = dot(A[i, :i], x[:i])
            old_values_sum = dot(A[i, i+1 :], x0[ i+1: ]) 
            x[i] = (b[i] - (old_values_sum + new_values_sum)) / A[i, i] 
            x[i] = dot(x[i], w) + dot(x0[i], (1 - w))  
 
        residual = np.linalg.norm(dot(A, x)-b )
        residual_list.append(residual)
        if (residual < r):
            print("SOR: The number of iterations is: {}".format(Ni))
            break 
        x0 = x 
        
    return x, residual_list



# xiao yang: Conjugate Gradient








# linear system:  5x1-x2+2x3=12
#                 3x1+8x2-2x3=-25
#                 x1+x2+4x3=6
#       x1=1, x2=-3, x3= 2

b=array([12,-25,6])
A=array([[5,-1,2],[3,8,-2],[1,1,4]])
x0=array([1000.0,1000.0,1000.0])
N= 100
r= 10**-9
w=1.1

solution_Jacobi, residual_list_Jacobi = Jacobi(A, b,x0, N, r)
print("Solution by Jacobi:")
print(solution_Jacobi)
solution_SOR, residual_list_SOR = SOR(A, b,x0, N, r)
print("Solution by SOR:")
print(solution_SOR)
solution_GaussSeidel, residual_list_GaussSeidel = GaussSeidel(A,b,x0,N, r)
print("Solution by GaussSeidel:")
print(solution_GaussSeidel)



# plotting residual to show convergence
plt.plot(residual_list_Jacobi,"ok",label="Jacobi")
plt.plot(residual_list_GaussSeidel,"or",label="Gauss-Seidel")
plt.plot(residual_list_SOR,"og",label="SOR (w={})".format(w))
plt.xlabel('iterations')
plt.ylabel('Residual')
plt.yscale('log')
plt.legend()
plt.savefig("../img/linearSystemSolverConvergence.png")
plt.show()
