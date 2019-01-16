# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 15:17:27 2019

@author: WS1
**using pysces
"""

from pysces import *
import numpy as np
import matplotlib.pyplot as plt
#from pylab import *
import scipy.fftpack
import random
import time
from scipy.stats import norm
import matplotlib.mlab as mlab

#mu, sigma = 0, 0.1 # mean and standard deviation
#s = np.random.normal(mu, sigma, 100)

total_num = 500
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
#    q[a,1] = np.random.rayleigh(0,1)
    xvort1[a,0] = i
    xvort1[a,1] = 0
    vel1[:,:] = panels.induced_velocity_single(q, xvort1, random.uniform(-1,1))
    a = a + 1

#print xvort1  
#print vel1
#print q
#
#plt.scatter (q[:,0],q[:,1])
#plt.show()

#plt.hist(q[:,1])
##plt.hist(q[:,0])
#plt.show()

t = np.linspace(0,50,total_num)
vel_tot_mag = np.zeros((total_num,1))
vel_tot_mag = (vel1[:,0]**2 + vel1[:,1]**2)**0.5
#print vel_tot_mag

#fig, ax = plt.subplots()
#u1 = plt.plot(t, vel1[:,0], c='g', label='u')
#v1 = plt.plot(t, vel1[:,1], c='r', label='v')
#vel = plt.plot(t, vel_tot_mag[:], c='k', label='U')


#plt.xlabel('time')
#plt.ylabel('induced_velocity u')
#plt.legend(loc='upper right')
##savefig(str(number_of_vortices)+'_vortices_u_velocity_rand_5000.pdf')
#plt.show()

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


import matplotlib.mlab as mlab

def L_p(x):
	return 20*np.log10(np.abs(x)/2.e-5)

blocksize=100
j=0

(werte,freq)=plt.psd(vel_tot_mag[:], NFFT=blocksize, Fs=100, detrend=mlab.detrend_none,window=mlab.window_hanning, noverlap=4, pad_to=None,sides='default', scale_by_freq='True')
#xx=30


pegel=L_p(werte)
#freq=freq[xx:]
plt.plot()

	
#Sr=freq*durchmesser/U
#solldrucklist=[1000,1200,1500,1700]
alphalist=[0.25,0.5,0.75,1.0]
stylelist=['--','-.',':','-']

plt.figure(1,figsize=(8.4/2.54,4.0/2.54))
plt.semilogx(freq,pegel,linestyle=stylelist[j],linewidth=0.6,alpha=alphalist[j])
j+=1
#plt.savefig('SPL_10000_windows.pdf')
plt.show()

fig = plt.figure(figsize=(8,40))

ax0 = fig.add_subplot(511).set_xlabel('X')
ax0 = fig.add_subplot(511).set_ylabel('Y')
ax0 = fig.add_subplot(511).set_title('Vortices injection from X=0 using Normal distribution')
ax0 = plt.scatter(q[:,0],q[:,1], label='vortex')
ax0 = fig.add_subplot(511).legend(loc='best')

ax1 = fig.add_subplot(512)
ax1 = fig.add_subplot(512).set_xlabel('X')
ax1 = fig.add_subplot(512).set_ylabel('Y')
ax1 = fig.add_subplot(512).set_title('Histogram for vortices injection from X=0')
ax1 = plt.hist(q[:,1], orientation='horizontal')

ax2 = fig.add_subplot(513)
ax2 = fig.add_subplot(513).set_xlabel('time')
ax2 = fig.add_subplot(513).set_ylabel('induced velocity')
ax2 = fig.add_subplot(513).set_title('induced velocity with respect to time')
ax2 = plt.plot(t, vel1[:,0], c='g', label='u (x-component of induced velocity)')
ax2 = plt.plot(t, vel1[:,1], c='r', label='v (y-component of induced velocity)')
ax2 = plt.plot(t, vel_tot_mag[:], c='k', label='U (total induced velocity)')
ax2 = fig.add_subplot(513).legend(loc='best')

ax3 = fig.add_subplot(514)
ax3 = fig.add_subplot(514).set_ylabel('induced velocity')
ax3 = fig.add_subplot(514).set_xlabel('Frequency (Hz)')
ax3 = fig.add_subplot(514).set_title('FFT of induced velocity')
ax3 = plt.semilogx(xf, 2.0/N * np.abs(yf[:N//2]))

ax4 = fig.add_subplot(515)
ax4 = fig.add_subplot(515).set_ylabel('Power Spectral Density (dB)')
ax4 = fig.add_subplot(515).set_xlabel('Frequency (Hz)')
ax4 = fig.add_subplot(515).set_title('Power Spectral Density (PSD) of induced velocity')
ax4 = plt.semilogx(freq,pegel,linewidth=0.6)

fig.savefig('multipleplots.pdf')
plt.show()