# different solvers for system of equations
# chennuo: complicated tested function

#importing libraries
from numpy import array, zeros, diag, diagflat, dot


def Jacobi(x): #Ax=b,we are finding x
    y = 2 * x
    return y

y = Jacobi(2)
print(y)

# xiao yang:
def Jacobi(A, b, N): #Ax=b,we are finding x, N = number of iterations
    x=zeros(len(A[0]))
    D= diag(A)
    r= A-diagflat(D)
    #iterate for N times
    for i in range(N) :
        x= (b-dot(r,x))/D
    return x


# linear system:  5x1-x2+2x3=12
#                 3x1+8x2-2x3=-25
#                 x1+x2+4x3=6

b=array([12,-25,6])
A=array([[5,-1,2],[3,8,-2],[1,1,4]])
x=array([0,0,0])
N= 10
solution = Jacobi(A, b, N)
print(solution)