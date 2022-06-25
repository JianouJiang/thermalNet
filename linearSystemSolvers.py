# different solvers for system of equations
# chennuo: complicated tested function

#importing libraries
import numpy as np
from numpy import array, zeros, diag, diagflat, dot


def Jacobi(x): #Ax=b,we are finding x
    y = 2 * x
    return y

y = Jacobi(2)
#print(y)

# xiao yang:
def Jacobi(A, b, N, x= None , r=10**-6): #Ax=b,we are finding x, N = number of iterations, r is the residual
    if (x is None) :
        x=zeros(len(A[0]))
    D= diag(A)
    R= A-diagflat(D)
    residual =1000
    xn_minus1=x
    #iterate for N times
    for i in range(N) :
        if residual <r:
            print(residual)
            break
        x= (b-dot(R,x))/D
        difference= xn_minus1-x
        magnitude= np.linalg.norm(difference)
        residual =magnitude/len(x) # TODO: check the grammar
        xn_minus1=x
        print (residual)
    return x


# kaiyi, Successive over-relaxation
# x0 = our initial guess, N = number of iterations, T = tolerance
# w = relaxation factor, 1<w<2. If w=1, it's same as Gauss-Seidel Method

def SOR(A, b, x0, N, T, w): 
    n = b.shape
    x = x0 
    for step in range (1, N): 
        for i in range(n[0]): 
            new_values_sum = dot(A[i, :i], x[:i])
            old_values_sum = dot(A[i, i+1 :], x0[ i+1: ]) 
            x[i] = (b[i] - (old_values_sum + new_values_sum)) / A[i, i] 
            x[i] = dot(x[i], w) + dot(x0[i], (1 - w))  
 
        if (np.linalg.norm(dot(A, x)-b ) < T):
            print(step) 
            break 
        x0 = x
        
    print("X = {}".format(x)) 
    print("The number of iterations is: {}".format(step))
    return x

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
T=10**-7
solution = Jacobi(A, b, N,x0, r)
print(solution)