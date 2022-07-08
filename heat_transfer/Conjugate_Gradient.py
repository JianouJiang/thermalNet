import numpy as np
from numpy import array, zeros, diag, diagflat, dot
#b=array([24,24,0])
#A=array([[12,-6,0],[-6,12,-6],[0,-6,6]])


b=array([12,-25,6])
A=array([[5,-1,2],[3,8,-2],[1,1,4]])

#A=array([[0.7444,-0.5055,-0.0851],[-0.5055 , 3.4858 , 0.0572],[-0.0851 , 0.0572 , 0.4738]])
#b=array([-0.0043 , 2.2501,  0.2798])
def CG(A,b,N,x=None):
    if (x is None) :   #x0,r0,p0
        x=zeros(len(A[0]))
    r=b-dot(A,x)
    p=r
    for i in range(N):
        alpha=np.divide(dot(r.transpose(),r),dot(p.transpose(),dot(A,p)))
        x=x+alpha*p  #x n+1
        print(x)
        r_n1= r-alpha* dot(A,p) #r n+1
        beta=np.divide(dot(r_n1.transpose(),r_n1),dot(r.transpose(),r)) #beta n
        p=r_n1+beta*p
        r=r_n1
    return (x)

def LinearCG(A, b, x0, tol=1e-5):
    xk = x0
    rk = np.dot(A, xk) - b
    pk = -rk
    rk_norm = np.linalg.norm(rk)
    
    num_iter = 0
    curve_x = [xk]
    while rk_norm > tol:
        apk = np.dot(A, pk)
        rkrk = np.dot(rk, rk)
        
        alpha = rkrk / np.dot(pk, apk)
        xk = xk + alpha * pk
        rk = rk + alpha * apk
        beta = np.dot(rk, rk) / rkrk
        pk = -rk + beta * pk
        
        num_iter += 1
        curve_x.append(xk)
        rk_norm = np.linalg.norm(rk)
        print('Iteration: {} \t x = {} \t residual = {:.4f}'.
              format(num_iter, xk, rk_norm))
    
    print('\nSolution: \t x = {}'.format(xk))
        
    return np.array(curve_x)


#solution = LinearCG(A, b,[1.1,-2.9,2.1])
#print(solution)
# linear system:  5x1-x2+2x3=12
#                 3x1+8x2-2x3=-25
#                 x1+x2+4x3=6
#       x1=1, x2=-3, x3= 2
def conjugate_gradient(A, b, x=None, max_iter=512, reltol=1e-6, verbose=True):
    """
    Implements conjugate gradient method to solve Ax=b for a large matrix A that is not
    computed explicitly, but given by the linear function A. 
    """
    if verbose:
        print("Starting conjugate gradient...")
    if x is None:
        x=np.zeros_like(b)
    # cg standard
    r=b-np.dot(A,x)
    d=r
    rsnew=np.sum(r.conj()*r).real
    rs0=rsnew
    if verbose:
        print("initial residual: {}".format(rsnew))
    ii=0
    while ((ii<max_iter) and (rsnew>(reltol**2*rs0))):
        ii=ii+1
        Ad=np.dot(A,d)
        alpha=rsnew/(np.sum(d.conj()*Ad))
        x=x+alpha*d
        if ii%50==0:
            #every now and then compute exact residual to mitigate
            # round-off errors
            r=b-np.dot(A,x)
            d=r
        else:
            r=r-alpha*Ad
        rsold=rsnew
        rsnew=np.sum(r.conj()*r).real
        d=r+rsnew/rsold*d
        if verbose:
            print("{}, residual: {}".format(ii, rsnew))
    return x


x0=[10,-30,20]
solution=conjugate_gradient(A, b, x=x0, max_iter=5120, reltol=1e-12)
print(solution)