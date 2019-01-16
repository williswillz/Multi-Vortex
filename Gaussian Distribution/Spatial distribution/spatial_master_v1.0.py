# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 15:17:27 2019

@author: WS1
"""

from pysces import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import scipy.fftpack
import random
import time
from scipy.stats import norm

#mu, sigma = 0, 0.1 # mean and standard deviation
#s = np.random.normal(mu, sigma, 100)

total_num = 100
panels = Vortices()
panels.core_radius = 0.01

xvort1 = np.zeros((total_num,2))
vel1 = np.zeros((total_num,2))
q = np.zeros((total_num,2))
a = 0

for i in range (0, total_num):
    i = i * 0.1
    q[a,0] = np.random.uniform(0,50,1)
    q[a,1] = np.random.normal(0,1,1)
    xvort1[a,0] = i
    xvort1[a,1] = 0
    vel1[:,:] = panels.induced_velocity_single(q, xvort1, random.uniform(-1,1))
    a = a + 1

print xvort1  
#print vel1
print q

plt.scatter (q[:,0],q[:,1])
plt.show()

plt.hist(q[:,1])
#plt.hist(q[:,0])
plt.show()

t = np.linspace(0,50,total_num)
vel_tot_mag = np.zeros((total_num,1))
vel_tot_mag = (vel1[:,0]**2 + vel1[:,1]**2)**0.5
print vel_tot_mag

fig, ax = plt.subplots()
u1 = plt.plot(t, vel1[:,0], c='g', label='u')
v1 = plt.plot(t, vel1[:,1], c='r', label='v')
vel = plt.plot(t, vel_tot_mag[:], c='k', label='U')


plt.xlabel('time')
plt.ylabel('induced_velocity u')
plt.legend(loc='upper right')
#savefig(str(number_of_vortices)+'_vortices_u_velocity_rand_5000.pdf')
plt.show()

# Number of samplepoints
N = total_num
# sample spacing
T = 1.0 / (total_num/2)
x = np.linspace(0.0, N*T, N)
y = vel_tot_mag
yf = scipy.fftpack.fft(y)
xf = np.linspace(0.0, 1.0/(2.0*T), N/2)

fig, ax = plt.subplots()
ax.plot(xf, 2.0/N * np.abs(yf[:N//2]))
ax.set_xscale('log')
plt.xlabel('frequency (Hz)')
plt.ylabel('induced_velocity')
#savefig('frequency_(5,1)_rand_5000'+str(N)+str(number_of_vortices)+'.pdf')
plt.show()