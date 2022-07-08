# different solvers for system of equations
#importing libraries
import numpy as np
from numpy import array, zeros, diag, diagflat, dot
import matplotlib.pyplot as plt
from scipy.linalg import solve
import time
# defining linear system solver functions


# xiao yang: Jacobi
# x0 = our initial guess, N = number of iterations, r = tolerance

def Jacobi(A, b,x0= None, N=1024, r=10**-6): #Ax=b,we are finding x, N = number of iterations, r is the residual
    st = time.time()
    residual_list = []
    x = x0
    if (x0 is None) :
        x=zeros(len(A[0]))
    D= diag(A)
    R= A-diagflat(D)
    residual =1000
    xn_minus1=x
    #iterate for N times
    Ni=0
    while (Ni<N):
        Ni = Ni+1
        x= (b-dot(R,xn_minus1))/D
        residual = np.linalg.norm(dot(A, x)-b)
        residual_list.append(residual)
        if residual <r:
            et = time.time()
            duration = et-st
            print("Jacobi: The final residual is: {} with {} iterations in {}s.".format(residual, Ni,duration))
            return x, residual_list
        xn_minus1=x
    et = time.time()
    duration = et-st
    print("Jacobi: The final residual is: {} with {} iterations in {}s.".format(residual_list[-1], Ni,duration))       
    return x, residual_list




# chen nuo: LU Decomposition method


def LU_Decomposition(A,b,x=None,N=1024, r=10**-6):
    st = time.time()
    residual_list = []
    if (x is None) :
        x=zeros(len(A[0]))
    xn_minus1 = x
    L = np.tril(A) # define the lower triangular matrix (L) for A
    U = A - L # define the upper triangular matrix (U)
    residual = 1000
    Ni=0
    while (Ni<N):
        Ni = Ni+1
        # define x as the dot product of inverse L and B mins the dot product of U and x
        x = np.dot(np.linalg.inv(L),b-np.dot(U,xn_minus1))
        
        residual = np.linalg.norm(dot(A, x)-b)
        residual_list.append(residual)
        if residual < r:
            et = time.time()
            duration = et-st
            print("LU_Decomposition: The final residual is: {} with {} iterations in {}s.".format(residual, Ni,duration))
            return x, residual_list
        xn_minus1 = x
    et = time.time()
    duration = et-st
    print("LU_Decomposition: The final residual is: {} with {} iterations in {}.".format(residual_list[-1], Ni,duration))   
    return x, residual_list




# chennuo: Gauss-Seidel method
# By using Gauss-Seidel method, we have 3 elements - matrix a, matrix b and solution X
def Gauss_Seidel(A,b,x0=None,N=1024,r=10**-6):
    st = time.time()
    residual_list = []
    if (x0 is None) :
        x0=zeros(len(A[0]))
    x = x0
    n = len(A) # find length of matrix a(3)
    
    Ni=0
    while (Ni<N):
        Ni = Ni+1
        xn_minus1  = x.copy()
        #Loop over rows
        for i in range(A.shape[0]):
            x[i] = (b[i] - np.dot(A[i,:i], x[:i]) - np.dot(A[i,(i+1):], xn_minus1[(i+1):])) / A[i ,i]
            
        #Stop condition 
        residual = np.linalg.norm(dot(A, x)-b)
        residual_list.append(residual)
        if  residual < r:
            et = time.time()
            duration = et-st
            print("Gauss_Seidel: The final residual is: {} with {} iterations in {}s.".format(residual, Ni,duration))
            return x, residual_list
    et = time.time()
    duration = et-st
    print("Gauss_Seidel: The final residual is: {} with {} iterations in {}s.".format(residual, Ni,duration))        
    return x, residual_list



# kaiyi: Successive over-relaxation (SOR)
# x0 = our initial guess, N = number of iterations, r = tolerance
# w = relaxation factor, 1<w<2. If w=1, it's same as Gauss-Seidel Method

def SOR(A, b, x0=None, N=1024, r=10**-6, w=1.5): 
    st = time.time()
    residual_list = []
    n = b.shape
    if (x0 is None):
        x0=zeros(len(A[0]))
    x = x0
    Ni=0
    while (Ni<N):
        Ni = Ni+1
        for i in range(n[0]): 
            new_values_sum = dot(A[i, :i], x[:i])
            old_values_sum = dot(A[i, i+1 :], x0[ i+1: ]) 
            x[i] = (b[i] - (old_values_sum + new_values_sum)) / A[i, i] 
            x[i] = dot(x[i], w) + dot(x0[i], (1 - w))  
 
        residual = np.linalg.norm(dot(A, x)-b)
        residual_list.append(residual)
        if (residual < r):
            et = time.time()
            duration = et-st
            print("SOR: The final residual is: {} with {} iterations in {}s.".format(residual, Ni,duration))
            return x, residual_list
        x0 = x 
    et = time.time()
    duration = et-st
    print("SOR: The final residual is: {} with {} iterations in {}s.".format(residual_list[-1], Ni,duration))   
    return x, residual_list



# xiao yang: Conjugate Gradient
def Conjugate_Gradient(A, b, x0=None, N=1024, reltol=1e-6, verbose=True):
    """
    Implements conjugate gradient method to solve Ax=b for a large matrix A that is not
    computed explicitly, but given by the linear function A. 
    """
    st = time.time()
    residual_list = []
    if x0 is None:
        x0=zeros(len(A[0]))
    x = x0
    # cg standard
    r=b-np.dot(A,x)
    d=r
    rsnew=np.sum(r.conj()*r).real
    rs0=rsnew
    Ni=0
    while (Ni<N):
        residual = np.linalg.norm(dot(A, x)-b)
        residual_list.append(residual)
        Ni=Ni+1
        Ad=np.dot(A,d)
        alpha=rsnew/(np.sum(d.conj()*Ad))
        x=x+alpha*d
        if Ni%50==0:
            #every now and then compute exact residual to mitigate
            # round-off errors
            r=b-np.dot(A,x)
            d=r
        else:
            r=r-alpha*Ad
        rsold=rsnew
        rsnew=np.sum(r.conj()*r).real
        d=r+rsnew/rsold*d

        if (residual<reltol):
            et = time.time()
            duration = et-st
            print("Conjugate Gradient: The final residual is: {} with {} iterations in {}s.".format(residual,Ni,duration))
            return x, residual_list

    et = time.time()
    duration = et-st
    if verbose:
        print("Conjugate Gradient: The final residual is: {} with {} iterations in {}s.".format(residual,Ni,duration))
    return x, residual_list



'''

# linear system:  5x1-x2+2x3=12
#                 3x1+8x2-2x3=-25
#                 x1+x2+4x3=6
#       x1=1, x2=-3, x3= 2

b=array([12,-25,6])
A=array([[5,-1,2],[3,8,-2],[1,1,4]])
x0=array([10000.0,10000.0,10000.0])
N= 100
r= 10**-15
w=1.1

solution_Jacobi, residual_list_Jacobi = Jacobi(A, b,x0, N, r)
print("Solution by Jacobi:")
print(solution_Jacobi)
solution_SOR, residual_list_SOR = SOR(A, b,x0, N, r)
print("Solution by SOR:")
print(solution_SOR)
solution_LU_Decomposition, residual_list_LU_Decomposition = LU_Decomposition(A,b,x0,N, r)




# plotting residual to show convergence
plt.plot(residual_list_Jacobi,"ok",label="Jacobi")
plt.plot(residual_list_SOR,"og",label="SOR (w={})".format(w))
plt.xlabel('iterations')
plt.ylabel('Residual')
plt.yscale('log')
plt.legend()
# plt.savefig("../img/linearSystemSolverConvergence.png")
plt.show()

'''