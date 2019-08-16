# -*- coding: utf-8 -*-
"""
Created on Thu May  9 11:39:01 2019

@author: WS1
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import psd

length = 25. #nondimensional length of window
height = 0.2 #window height
N = 10000 #number of vortices
gammas = 1. #vortex strength RMS (normal distribution)
rscale = 0.1 #vortex size scale (rayleigh distribution parameter)
t0 = -10.#start time for observation of convection
t1 = 10.#end time
ts = 0.001 # time step
v0 = 5 #convection speed


nn = 20 # number of injection lines
#vortx = (np.random.uniform(low=-length/2,high=length/2,size=N/nn))[np.newaxis].T
#vortx = np.hstack(([vortx]*nn))
#vortx[:,0] = np.random.uniform(low=-length/2,high=length/2,size=(N/nn,20))
vortx = np.random.uniform(low=-length/2,high=length/2,size=(N/nn,nn))
vorty = np.ones_like(vortx)
K = np.arange(-0.1,0.1,0.01)
vorty = K * vorty
vortx = vortx.T
vorty = vorty.T
vortx = vortx.flatten()
vorty = vorty.flatten()
#vortx = np.concatenate((vortx1, vortx2, vortx3, vortx4,vortx5,vortx6), axis=0)
#vorty = np.concatenate((vorty1, vorty2, vorty3, vorty4, vorty5, vorty6), axis=0)

#vorty = np.random.uniform(low=-height/2,high=height/2,size=N)
#vortx = np.random.uniform(low=-length/2,high=length/2,size=N)
# vorty = np.full((N), 0.05)

vortX = np.vstack((vortx,vorty))
gamma = np.random.normal(scale=gammas,size=N)
rho = np.random.rayleigh(scale=rscale,size=N)
plt.figure(1)
plt.scatter(vortx,vorty,s=0.1)
# plt.hist(rho)
#plt.scatter(vortx1,vorty1,s=0.1)
#plt.scatter(vortx2,vorty2,s=0.1)
#plt.scatter(vortx3,vorty3,s=0.1)
#plt.scatter(vortx4,vorty4,s=0.1)