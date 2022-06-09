"""
@author: Jianou Jiang
"""
print("starting to import libs")
#export PYTHONPATH=c:\users\sann7537\appdata\local\programs\python\python35\lib\site-packages:$PYTHONPATH.
import sys

#!$sys.executable -m pip install --upgrade https://storage.googleapis.com/tensorflow/windows/cpu/tensorflow-0.12.0rc0-cp35-cp35m-win_amd64.whl

import tensorflow as tf

#print(tf.version)
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from scipy.interpolate import griddata
import time
from itertools import product, combinations
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
#from plotting import newfig, savefig
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

np.random.seed(1234)
tf.set_random_seed(1234)

import subprocess
subprocess.check_output(['latex', '--version'])
#!/usr/bin/env python3
# -*- coding: utf-8 -*-


print("importing plotting libs")

import numpy as np
import matplotlib as mpl
#mpl.use('pgf')

def figsize(scale, nplots = 1):
    fig_width_pt = 390.0                          # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = nplots*fig_width*golden_mean              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 10,               # LaTeX default is 10pt font.
    "font.size": 10,
    "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.figsize": figsize(1.0),     # default fig size of 0.9 textwidth
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        ]
    }
mpl.rcParams.update(pgf_with_latex)

import matplotlib.pyplot as plt

# I make my own newfig and savefig functions
def newfig(width, nplots = 1):
    fig = plt.figure(figsize=figsize(width, nplots))
    ax = fig.add_subplot(111)
    return fig, ax

def savefig(filename, crop = True):
    if crop == True:
#        plt.savefig('{}.pgf'.format(filename), bbox_inches='tight', pad_inches=0)
        plt.savefig('{}.pdf'.format(filename), bbox_inches='tight', pad_inches=0)
        plt.savefig('{}.eps'.format(filename), bbox_inches='tight', pad_inches=0)
    else:
#        plt.savefig('{}.pgf'.format(filename))
        plt.savefig('{}.pdf'.format(filename))
        plt.savefig('{}.eps'.format(filename))

## Simple plot
#fig, ax  = newfig(1.0)
#
#def ema(y, a):
#    s = []
#    s.append(y[0])
#    for t in range(1, len(y)):
#        s.append(a * y[t] + (1-a) * s[t-1])
#    return np.array(s)
#    
#y = [0]*200
#y.extend([20]*(1000-len(y)))
#s = ema(y, 0.01)
#
#ax.plot(s)
#ax.set_xlabel('X Label')
#ax.set_ylabel('EMA')
#
#savefig('ema')
print("finished importing libs")


class PhysicsInformedNN:
    # Initialize the class
    def __init__(self, x, y, t, T, phi1, phi2, layers):
        
        X = np.concatenate([x, y, t], 1)
        
        self.lb = X.min(0)
        self.ub = X.max(0)
                
        self.X = X
        
        self.x = X[:,0:1]
        self.y = X[:,1:2]
        self.phi1 = X[:,1:2] # this is a scalar field, to distinguish steel wall from coke and others
        
        self.phi2 = X[:,1:2] # this is a scalar field, from -2 to 2, can be a level-set function to differentiate material1(solid:steel wall+coke) from material2(air)
        self.t = X[:,2:3]
        
        self.T = T
        #self.v = v
        
        self.layers = layers
        
        # Initialize NN
        self.weights, self.biases = self.initialize_NN(layers)        
        
        # Initialize parameters
        self.lambda_1 = tf.Variable([0.0], dtype=tf.float32)
        self.lambda_2 = tf.Variable([0.0], dtype=tf.float32)
        
        # tf placeholders and graph
        self.sess = tf.Session(config=tf.ConfigProto(allow_soft_placement=True,
                                                     log_device_placement=True))
        
        self.x_tf = tf.placeholder(tf.float32, shape=[None, self.x.shape[1]])
        self.y_tf = tf.placeholder(tf.float32, shape=[None, self.y.shape[1]])
        self.phi1_tf = tf.placeholder(tf.float32, shape=[None, self.y.shape[1]]) # this is a scalar field, steel wall
        self.phi2_tf = tf.placeholder(tf.float32, shape=[None, self.y.shape[1]]) # this is a scalar field
        self.t_tf = tf.placeholder(tf.float32, shape=[None, self.t.shape[1]])
        
        self.T_tf = tf.placeholder(tf.float32, shape=[None, self.T.shape[1]])
        #self.v_tf = tf.placeholder(tf.float32, shape=[None, self.v.shape[1]])
        
        self.T_pred, self.phi2_pred, self.f_T_pred, self.f_i_pred, self.f_phi2_pred = self.net_HT(self.x_tf, self.y_tf, self.t_tf) 
        
        """ f_i_pred: reinforce the interface heat flux """
        self.loss = tf.reduce_sum(tf.square(self.T_tf - self.T_pred)) + \
                    tf.reduce_sum(tf.square(self.f_T_pred)) + \
                    tf.reduce_sum(tf.square(self.f_i_pred)) + \
                    tf.reduce_sum(tf.square(self.f_phi2_pred))
        """f_phi2_pred: this is to make sure that the interface between coke and air stays normal...iterative reinitialisation method!"""
        """maxiter: 500000"""
        self.optimizer = tf.contrib.opt.ScipyOptimizerInterface(self.loss, 
                                                                method = 'L-BFGS-B', 
                                                                options = {'maxiter': 50,
                                                                           'maxfun': 50000,
                                                                           'maxcor': 50,
                                                                           'maxls': 50,
                                                                           'ftol' : 1.0 * np.finfo(float).eps})        
        
        self.optimizer_Adam = tf.train.AdamOptimizer()
        self.train_op_Adam = self.optimizer_Adam.minimize(self.loss)                    
        
        init = tf.global_variables_initializer()
        self.sess.run(init)

    def initialize_NN(self, layers):        
        weights = []
        biases = []
        num_layers = len(layers) 
        for l in range(0,num_layers-1):
            W = self.xavier_init(size=[layers[l], layers[l+1]])
            b = tf.Variable(tf.zeros([1,layers[l+1]], dtype=tf.float32), dtype=tf.float32)
            weights.append(W)
            biases.append(b)        
        return weights, biases
        
    def xavier_init(self, size):
        in_dim = size[0]
        out_dim = size[1]        
        xavier_stddev = np.sqrt(2/(in_dim + out_dim))
        return tf.Variable(tf.truncated_normal([in_dim, out_dim], stddev=xavier_stddev), dtype=tf.float32)
    
    def neural_net(self, X, weights, biases):
        num_layers = len(weights) + 1
        
        H = 2.0*(X - self.lb)/(self.ub - self.lb) - 1.0
        for l in range(0,num_layers-2):
            W = weights[l]
            b = biases[l]
            H = tf.tanh(tf.add(tf.matmul(H, W), b))
        W = weights[-1]
        b = biases[-1]
        Y = tf.add(tf.matmul(H, W), b)
        return Y
        
    def net_HT(self, x, y, t):
        print("net_HT()")
        lambda_1 = self.lambda_1 # heat conductivity of material 1
        lambda_2 = self.lambda_2 # heat conductivity of material 2
        
        
        
        phi1 = self.phi1
        # now, we need to get the individual phi1 out for the following PDE, but maybe we dont have to...
        interface_phi1 = []
        dx = 0.0816
        e = 1e-4
        for i in range(len(phi1)):
            phi1i = phi1[i]
            #print(phi1i)
            if -dx/2-e<=phi1i<dx/2: # at the interface between steel and coke
                interface_phi1.append([1])
            else:
                interface_phi1.append([0])
                
        T_and_phi2 = self.neural_net(tf.concat(1,[x,y,t]), self.weights, self.biases)
        T = T_and_phi2[:,0:1]
        
        
        phi2 = T_and_phi2[:,1:2]
      
        interface_phi2 = []
        with tf.Session() as sess:   
            tf.initialize_all_variables().run() # need to initialize all variables

            print("-----shape--------")
            print(self.phi2.shape) # (5000,1)
            print(self.phi2.shape[0]) # 5000
            for i in range(self.phi2.shape[0]):
                j=0
                phi2i = self.phi2[i,j] 
                #print(phi2i)
                if -dx/2 - e <= phi2i < dx/2:
                    #print("here!!!!!!!!!")
                    interface_phi2.append([1.0])
                else:
                    interface_phi2.append([0.0])
           
        print("interface_phi2")
        #print(interface_phi2)
        
        
        _lambda = []
        for i in range(len(self.phi1)):
            phi1i = phi1[i]
            if phi1i>0: # on the steel side
                _lambda.append([lambda_1])
            else: # on the coke side
                _lambda.append([lambda_2])
        
        T_t = tf.gradients(T, t)[0]
        print("T_t:")
        print(T_t)
        T_x = tf.gradients(T, x)[0]
        print("T_x:")
        print(T_x)
        T_y = tf.gradients(T, y)[0]
        T_xx = tf.gradients(T_x, x)[0]
        T_yy = tf.gradients(T_y, y)[0]

        
        phi2_x = tf.gradients(phi2, x)[0]
        phi2_y = tf.gradients(phi2, y)[0]

        # TODO! the _lambda here is based on phi1 and phi2, also varies with temperature T.
        heatSource = 800.0 # K
        print(type(T_xx))
        f_T = T_t - _lambda*(T_xx + T_yy) - heatSource * tf.Variable(interface_phi2) # the heat source is only applied at the interface between solid and air...e.g. the inner layer
       
       # TODO! reinforce the interface (short version: i) between material 1 and material 2, as in, make sure that the heat fluxes are the same
        f_i = 0.0 # v_t + lambda_1*(T*v_x + v*v_y) + phi2_y - lambda_2*(v_xx + v_yy)
        
        # this is to reinitialise the level-set function so that the interface stays normal --> gradient of 1
        f_phi2 = phi2_x*phi2_x + phi2_y*phi2_y - 1.0

        
        return T, phi2, f_T, f_i, f_phi2
    
    def callback(self, loss, lambda_1, lambda_2):
        print('Loss: %.3e, l1: %.3f, l2: %.5f' % (loss, lambda_1, lambda_2))
      
    def train(self, nIter): 
        print("started training...") 
        tf_dict = {self.x_tf: self.x, self.y_tf: self.y, self.t_tf: self.t, self.phi2_tf: self.phi2,
                   self.T_tf: self.T}
        
        start_time = time.time()
        for it in range(nIter): # also stops when reads 0 from stopit.txt file
            self.sess.run(self.train_op_Adam, tf_dict)
            
            # Print
            if it % 10 == 0:
                elapsed = time.time() - start_time
                loss_value = self.sess.run(self.loss, tf_dict)
                lambda_1_value = self.sess.run(self.lambda_1)
                lambda_2_value = self.sess.run(self.lambda_2)
                print('It: %d, Loss: %.3e, l1: %.3f, l2: %.5f, Time: %.2f' % 
                      (it, loss_value, lambda_1_value, lambda_2_value, elapsed))
                start_time = time.time()
            
            # breaks when reads - from stopit.txt file
            f = open("stopit.txt", "r")
            stopit = int(f.read())            
            if stopit == 1: # breaks the iteration
                break
            
        self.optimizer.minimize(self.sess,
                                feed_dict = tf_dict,
                                fetches = [self.loss, self.lambda_1, self.lambda_2],
                                loss_callback = self.callback)
        print("finished training...") 
    
    def predict(self, x_star, y_star, t_star):
        
        tf_dict = {self.x_tf: x_star, self.y_tf: y_star, self.t_tf: t_star}
        
        # i think the sess.run() is only taking out the value e.g. u_pred based on x,y,t in tf_dict
        T_star = self.sess.run(self.T_pred, tf_dict) 
        #v_star = self.sess.run(self.v_pred, tf_dict)
        #phi1_star = self.sess.run(self.phi1_pred, tf_dict)
        phi2_star = self.sess.run(self.phi2_pred, tf_dict)
        
        return T_star, phi2_star

def plot_solution(X_star, T_star, index):
    print("plot_solution...")
    lb = X_star.min(0)
    ub = X_star.max(0)
    nn = 200
    x = np.linspace(lb[0], ub[0], nn)
    y = np.linspace(lb[1], ub[1], nn)
    X, Y = np.meshgrid(x,y)
    
    T_star = griddata(X_star, T_star.flatten(), (X, Y), method='cubic')
    #print(U_star)
    plt.figure(index)
    plt.pcolor(X,Y,T_star, cmap = 'jet')
    plt.colorbar()
    #plt.savefig("/plot_solution.pdf")
    plt.show()
    
    return
    
    
def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/4
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)
        
        
if __name__ == "__main__": 
    print("main...")
    N_train = 5000
    
    layers = [3, 20, 20, 20, 20, 20, 20, 20, 20, 2]
    
    # Load Data
    data = scipy.io.loadmat('cylinder_nektar_wake.mat')
           
    T_star = data['U_star'] # N x 2 x T
    phi2_star = data['p_star'] # N x T, exact phi2 field that distinguishes solid and air, we will use this to compare with our predicted phi2 field
    t_star = data['t'] # T x 1
    t_star_new = []
    for i in range(len(t_star)):   
        a = np.empty(1)
        a.fill(t_star[i][0])
        
        t_star_new.append(a)
    t_star = np.array(t_star_new)  
    
    X_star = data['X_star'] # N x 2
    
    N = X_star.shape[0] # N=5000, X_star.shape = (5000,2) # 2 is for x and y, 5000=100*50 data points in total for 4*8 dimension box.
    
    time_steps = len(t_star) # t_star.shape[0] # 200 time steps
    print("time t="+str(type(time_steps)))
    
    # Rearrange Data 
    XX = np.tile(X_star[:,0:1], (1,time_steps)) # N x t
    YY = np.tile(X_star[:,1:2], (1,time_steps)) # N x t
    
    tt = []
    for i in range(N_train):   
        a = np.linspace(0,19.9,200)  
        tt.append(a)
    tt = np.array(tt)  
    #tt = np.tile(t_star, (1,N)).t # N x t
    print(tt)
    
    TT = T_star[:,0,:] # Temperature N x t
    #VV = T_star[:,1,:] # N x t
    Phi2Phi2 = phi2_star # N x t
    
    x = XX.flatten()[:,None] # NT x 1
    #for i in range(len(x)):
    #    print(x[i])
    y = YY.flatten()[:,None] # NT x 1
    t = tt.flatten()[:,None] # time, NT x 1
    print(t)
    
    T = TT.flatten()[:,None] # NT x 1
    #v = VV.flatten()[:,None] # NT x 1
    phi2 = Phi2Phi2.flatten()[:,None] # NT x 1
    #print(p)
    
    ######################################################################
    ######################## Noiseles Data ###############################
    ######################################################################
    # Training Data   
    print(time_steps)
    idx = np.random.choice(N*time_steps, N_train, replace=False)
    print(time_steps)
    x_train = x[idx,:]
    y_train = y[idx,:]
    phi1_train = y[idx,:] # for level-set function
    phi2_train = phi2[idx,:] # for level-set function
    t_train = t[idx,:]
    T_train = T[idx,:]
    #v_train = v[idx,:]
    #print(t_train)
    #print(v_train)

    # Training
    model = PhysicsInformedNN(x_train, y_train, t_train, T_train,  phi1_train, phi2_train, layers)
    model.train(200) # 200000
    
    # Test Data
    snap = np.array([100])
    x_star = X_star[:,0:1]
    y_star = X_star[:,1:2]
    t_star = TT[:,snap]
    #print(T_star)
    t_star = T_star[:,0,snap]
    #v_star = T_star[:,1,snap]
    phi2_star = phi2_star[:,snap]
    
    # Prediction
    T_pred, phi2_pred = model.predict(x_star, y_star, t_star)
    lambda_1_value = model.sess.run(model.lambda_1)
    lambda_2_value = model.sess.run(model.lambda_2)
    
    # Error
    error_T = np.linalg.norm(t_star-T_pred,2)/np.linalg.norm(t_star,2)
    #error_v = np.linalg.norm(v_star-v_pred,2)/np.linalg.norm(v_star,2)
    error_phi2 = np.linalg.norm(phi2_star-phi2_pred,2)/np.linalg.norm(phi2_star,2)

    error_lambda_1 = np.abs(lambda_1_value - 1.0)*100
    error_lambda_2 = np.abs(lambda_2_value - 0.01)/0.01 * 100
    
    print('Error u: %e' % (error_T))    
    
    print('Error p: %e' % (error_phi2))    
    print('Error l1: %.5f%%' % (error_lambda_1))                             
    print('Error l2: %.5f%%' % (error_lambda_2))                  
    
    # Plot Results
    print(len(X_star))
    print(X_star)
    
    plot_solution(X_star, T_pred, 1)
    plot_solution(X_star, phi2_pred, 1)
#    plot_solution(X_star, v_pred, 2)
#    plot_solution(X_star, p_pred, 3)    
#    plot_solution(X_star, p_star, 4)
#    plot_solution(X_star, p_star - p_pred, 5)
    
    # Predict for plotting
    lb = X_star.min(0)
    ub = X_star.max(0)
    nn = 200
    x = np.linspace(lb[0], ub[0], nn)
    y = np.linspace(lb[1], ub[1], nn)
    X, Y = np.meshgrid(x,y)
    
    TT_star = griddata(X_star, T_pred.flatten(), (X, Y), method='cubic')
    #VV_star = griddata(X_star, v_pred.flatten(), (X, Y), method='cubic')
    Phi2Phi2_star = griddata(X_star, phi2_pred.flatten(), (X, Y), method='cubic')
    Phi2_exact = griddata(X_star, phi2_star.flatten(), (X, Y), method='cubic')
    
    
    ######################################################################
    ########################### Noisy Data ###############################
    ######################################################################
    print("started noisy data...")
    noise = 0.01        
    T_train = T_train + noise*np.std(T_train)*np.random.randn(T_train.shape[0], T_train.shape[1])
    #v_train = v_train + noise*np.std(v_train)*np.random.randn(v_train.shape[0], v_train.shape[1])    

    # Training
    model = PhysicsInformedNN(x_train, y_train, t_train, T_train, phi1_train, phi2_train, layers)
    model.train(200) # 200000
        
    lambda_1_value_noisy = model.sess.run(model.lambda_1)
    lambda_2_value_noisy = model.sess.run(model.lambda_2)
      
    error_lambda_1_noisy = np.abs(lambda_1_value_noisy - 1.0)*100
    error_lambda_2_noisy = np.abs(lambda_2_value_noisy - 0.01)/0.01 * 100
        
    print('Error l1: %.5f%%' % (error_lambda_1_noisy))                             
    print('Error l2: %.5f%%' % (error_lambda_2_noisy))     

             
    
    ######################################################################
    ############################# Plotting ###############################
    ######################################################################    
    print("started plotting")
    # Load Data
    data_vort = scipy.io.loadmat('cylinder_nektar_t0_vorticity.mat')
           
    x_vort = data_vort['x'] 
    y_vort = data_vort['y'] 
    w_vort = data_vort['w'] 
    modes = np.asscalar(data_vort['modes'])
    nel = np.asscalar(data_vort['nel'])    
    
    xx_vort = np.reshape(x_vort, (modes+1,modes+1,nel), order = 'F')
    yy_vort = np.reshape(y_vort, (modes+1,modes+1,nel), order = 'F')
    ww_vort = np.reshape(w_vort, (modes+1,modes+1,nel), order = 'F')
    
    box_lb = np.array([1.0, -2.0])
    box_ub = np.array([8.0, 2.0])
    
    fig, ax = newfig(1.0, 1.2)
    ax.axis('off')
    
    ####### Row 0: Vorticity ##################    
    gs0 = gridspec.GridSpec(1, 2)
    gs0.update(top=1-0.06, bottom=1-2/4 + 0.12, left=0.0, right=1.0, wspace=0)
    ax = plt.subplot(gs0[:, :])
    
    for i in range(0, nel):
        h = ax.pcolormesh(xx_vort[:,:,i], yy_vort[:,:,i], ww_vort[:,:,i], cmap='seismic',shading='gouraud',  vmin=-3, vmax=3) 
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(h, cax=cax)
    
    ax.plot([box_lb[0],box_lb[0]],[box_lb[1],box_ub[1]],'k',linewidth = 1)
    ax.plot([box_ub[0],box_ub[0]],[box_lb[1],box_ub[1]],'k',linewidth = 1)
    ax.plot([box_lb[0],box_ub[0]],[box_lb[1],box_lb[1]],'k',linewidth = 1)
    ax.plot([box_lb[0],box_ub[0]],[box_ub[1],box_ub[1]],'k',linewidth = 1)
    
    ax.set_aspect('equal', 'box')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_title('Vorticity', fontsize = 10)
    
    
    ####### Row 1: Training data ##################
    ########      u(t,x,y)     ###################        
    print("u(t,x,y)")
    gs1 = gridspec.GridSpec(1, 2)
    gs1.update(top=1-2/4, bottom=0.0, left=0.01, right=0.99, wspace=0)
    ax = plt.subplot(gs1[:, 0],  projection='3d')
    ax.axis('off')

    r1 = [x_star.min(), x_star.max()]
    r2 = [data['t'].min(), data['t'].max()]       
    r3 = [y_star.min(), y_star.max()]
    
    for s, e in combinations(np.array(list(product(r1,r2,r3))), 2):
        if np.sum(np.abs(s-e)) == r1[1]-r1[0] or np.sum(np.abs(s-e)) == r2[1]-r2[0] or np.sum(np.abs(s-e)) == r3[1]-r3[0]:
            ax.plot3D(*zip(s,e), color="k", linewidth = 0.5)   

    ax.scatter(x_train, t_train, y_train, s = 0.1)
    ax.contourf(X,TT_star,Y, zdir = 'y', offset = t_star.mean(), cmap='rainbow', alpha = 0.8)
              
    ax.text(x_star.mean(), data['t'].min() - 1, y_star.min() - 1, '$x$')
    ax.text(x_star.max()+1, data['t'].mean(), y_star.min() - 1, '$t$')
    ax.text(x_star.min()-1, data['t'].min() - 0.5, y_star.mean(), '$y$')
    ax.text(x_star.min()-3, data['t'].mean(), y_star.max() + 1, '$u(t,x,y)$')    
    ax.set_xlim3d(r1)
    ax.set_ylim3d(r2)
    ax.set_zlim3d(r3)
    axisEqual3D(ax)
    
    ########      v(t,x,y)     ###################   
    print("v(t,x,y)")
    ax = plt.subplot(gs1[:, 1],  projection='3d')
    ax.axis('off')
    
    r1 = [x_star.min(), x_star.max()]
    r2 = [data['t'].min(), data['t'].max()]       
    r3 = [y_star.min(), y_star.max()]
    
    for s, e in combinations(np.array(list(product(r1,r2,r3))), 2):
        if np.sum(np.abs(s-e)) == r1[1]-r1[0] or np.sum(np.abs(s-e)) == r2[1]-r2[0] or np.sum(np.abs(s-e)) == r3[1]-r3[0]:
            ax.plot3D(*zip(s,e), color="k", linewidth = 0.5)   

    ax.scatter(x_train, t_train, y_train, s = 0.1)
    #ax.contourf(X,VV_star,Y, zdir = 'y', offset = t_star.mean(), cmap='rainbow', alpha = 0.8)
              
    ax.text(x_star.mean(), data['t'].min() - 1, y_star.min() - 1, '$x$')
    ax.text(x_star.max()+1, data['t'].mean(), y_star.min() - 1, '$t$')
    ax.text(x_star.min()-1, data['t'].min() - 0.5, y_star.mean(), '$y$')
    ax.text(x_star.min()-3, data['t'].mean(), y_star.max() + 1, '$v(t,x,y)$')    
    ax.set_xlim3d(r1)
    ax.set_ylim3d(r2)
    ax.set_zlim3d(r3)
    axisEqual3D(ax)
    
    savefig('NavierStokes_data') 

    
    fig, ax = newfig(1.015, 0.8)
    ax.axis('off')
    
    ######## Row 2: Pressure #######################
    ########      Predicted p(t,x,y)     ########### 
    print("Predicted p(t,x,y)")
    gs2 = gridspec.GridSpec(1, 2)
    gs2.update(top=1, bottom=1-1/2, left=0.1, right=0.9, wspace=0.5)
    ax = plt.subplot(gs2[:, 0])
    h = ax.imshow(Phi2Phi2_star, interpolation='nearest', cmap='rainbow', 
                  extent=[x_star.min(), x_star.max(), y_star.min(), y_star.max()], 
                  origin='lower', aspect='auto')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    fig.colorbar(h, cax=cax)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_aspect('equal', 'box')
    ax.set_title('Predicted pressure', fontsize = 10)
    
    ########     Exact p(t,x,y)     ########### 
    print("Exact p(t,x,y)")
    ax = plt.subplot(gs2[:, 1])
    h = ax.imshow(Phi2_exact, interpolation='nearest', cmap='rainbow', 
                  extent=[x_star.min(), x_star.max(), y_star.min(), y_star.max()], 
                  origin='lower', aspect='auto')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    fig.colorbar(h, cax=cax)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_aspect('equal', 'box')
    ax.set_title('Exact pressure', fontsize = 10)
    
    
    ######## Row 3: Table #######################
    gs3 = gridspec.GridSpec(1, 2)
    gs3.update(top=1-1/2, bottom=0.0, left=0.0, right=1.0, wspace=0)
    ax = plt.subplot(gs3[:, :])
    ax.axis('off')
    
    s = r'$\begin{tabular}{|c|c|}';
    s = s + r' \hline'
    s = s + r' Correct PDE & $\begin{array}{c}'
    s = s + r' u_t + (u u_x + v u_y) = -p_x + 0.01 (u_{xx} + u_{yy})\\'
    s = s + r' v_t + (u v_x + v v_y) = -p_y + 0.01 (v_{xx} + v_{yy})'
    s = s + r' \end{array}$ \\ '
    s = s + r' \hline'
    s = s + r' Identified PDE (clean data) & $\begin{array}{c}'
    s = s + r' u_t + %.3f (u u_x + v u_y) = -p_x + %.5f (u_{xx} + u_{yy})' % (lambda_1_value, lambda_2_value)
    s = s + r' \\'
    s = s + r' v_t + %.3f (u v_x + v v_y) = -p_y + %.5f (v_{xx} + v_{yy})' % (lambda_1_value, lambda_2_value)
    s = s + r' \end{array}$ \\ '
    s = s + r' \hline'
    s = s + r' Identified PDE (1\% noise) & $\begin{array}{c}'
    s = s + r' u_t + %.3f (u u_x + v u_y) = -p_x + %.5f (u_{xx} + u_{yy})' % (lambda_1_value_noisy, lambda_2_value_noisy)
    s = s + r' \\'
    s = s + r' v_t + %.3f (u v_x + v v_y) = -p_y + %.5f (v_{xx} + v_{yy})' % (lambda_1_value_noisy, lambda_2_value_noisy)
    s = s + r' \end{array}$ \\ '
    s = s + r' \hline'
    s = s + r' \end{tabular}$'
 
    ax.text(0.015,0.0,s)
    
    savefig('NavierStokes_prediction') 
    print("finished the whole simulation...")