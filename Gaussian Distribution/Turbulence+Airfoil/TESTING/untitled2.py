# -*- coding: utf-8 -*-
"""
Created on Fri May 10 09:31:35 2019

@author: WS1
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import psd
from steapy import *
import os
import numpy
from numpy import *
import math
from scipy import integrate, linalg
from matplotlib import pyplot
from pylab import *
from sympy.abc import x, y
import essentials
import h5py

q = essentials.naca_airfoil(0012, 101)
q_mid = (q[1:] + q[:-1])/2
length = 100. #nondimensional length of window
height = 0.1 #window height
N = 10000 #number of vortices
gammas = 1. #vortex strength RMS (normal distribution)
rscale = 0.1 #vortex size scale (rayleigh distribution parameter)
t0 = -10.#start time for observation of convection
t1 = 10.#end time
ts = 0.01 # time step
v0 = 5 #convection speed
t = np.arange(t0,t1,ts) # number of time-steps
nn = 20 # number of injection lines

vortx = np.random.uniform(low=-length/2,high=length/2,size=(N/nn,nn))
vorty = np.ones_like(vortx)
K = np.arange(-0.1,0.1,0.01)
vorty = K * vorty
vortx = vortx.T
vorty = vorty.T
vortx = vortx.flatten()
vorty = vorty.flatten()
vortX = np.vstack((vortx,vorty))
gamma = np.random.normal(scale=gammas,size=N)
rho = np.random.rayleigh(scale=rscale,size=N)

for i in range(len(t)):
    i = i * ts * v0
    obsx = q_mid[:,0] + (length/2)
    obsx = (obsx - i )
    obsy = q_mid[:,1]
    obsX = np.vstack((obsx,obsy))
    dist = obsX[:,:,np.newaxis]-vortX[:,np.newaxis,:] # dim 2 x timesteps x N
    r = np.sqrt((dist*dist).sum(0)) # dim timesteps x N