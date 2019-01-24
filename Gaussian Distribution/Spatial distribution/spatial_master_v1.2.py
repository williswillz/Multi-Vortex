# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 10:14:07 2019

@author: WS1
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 08:26:17 2019

@author: WS1

independent of the main code
addition of the variable
** core radius (size of the vortex)
** gamma (strength of the vortex)
"""

import numpy as np
#from pysces import *
#import numpy as np
import matplotlib.pyplot as plt
#from pylab import *
import scipy.fftpack
#import random
import time
#from scipy.stats import norm
import matplotlib.mlab as mlab
from matplotlib.lines import Line2D


start = time.time()

def induced_velocity_single(x, xvort, gam, core_radius):
        """Compute velocity induced at points x by a single vortex

        Parameters
        ----------
        x : 2d array
            Locations at which to compute induced velocity.  Expressed as
            column vectors (i.e., shape should be (n,2))
        xvort : 1d array
            Location of vortex (shape should be (2,))
        gam : float
            Strength of vortex

        Notes
        -----
        Induced velocity is

        .. math:: u_\theta = \frac{\Gamma}{2 \pi r}

        where r is the distance between the point and the vortex.  If this
        distance is less than :class:`core_radius` :math:`r_0`, the velocity is
        regularized as solid-body rotation, with

        .. math:: u_\theta = \frac{\Gamma r}{2\pi r_0^2}
        """
        r = np.array(x, ndmin=2) - np.array(xvort)
        rsq = np.maximum(np.sum(r * r, 1), core_radius**2)
        # alternative regularization (Krasny, Eldredge)
        # rsq = np.sum(r * r, 1) + core_radius**2
        vel = np.transpose(np.array([-r[:,1], r[:,0]]))
        vel = gam / (2 * np.pi) * vel / rsq[:,np.newaxis]
        return np.squeeze(vel)
    

#mu, sigma = 0, 0.1 # mean and standard deviation
#s = np.random.normal(mu, sigma, 100)

total_num = 500

xvort1 = np.zeros((total_num,2))
vel1 = np.zeros((total_num,2))
a = 0

q = np.zeros((total_num,2))
q[:,0] = np.random.uniform(0,50,total_num)
q[:,1] = np.random.uniform(-0.5,0.5,total_num)
# q[:,1] = np.random.normal(0,1,total_num)

core_radius = np.zeros((total_num,1))
core_radius[:,0] = np.random.uniform(0,1,total_num)
#print core_radius

gam = np.zeros((total_num,1))
gam[:,0] = np.random.uniform(-1,1,total_num)
#print gam

for i in range (0, total_num):
    i = i * 0.1
    xvort1[a,0] = i
    xvort1[a,1] = 0
    vel1[:,:] = induced_velocity_single(q, xvort1, gam, core_radius[a])
    # vel1[:,:] = induced_velocity_single(q, xvort1, gam, core_radius)
#    print vel1
    a = a + 1
    print a

#savetxt('core_radius'+'.csv',np.column_stack((core_radius[:,0])), fmt='%5s', delimiter=',')

#print xvort1  
#print vel1
#plt.scatter (q[:,0],q[:,1])
#plt.show()

#plt.hist(q[:,1])
##plt.hist(q[:,0])
#plt.show()

t = np.linspace(0,50,total_num)
vel_tot_mag = np.zeros((total_num,1))
vel_tot_mag = (vel1[:,0]**2 + vel1[:,1]**2)**0.5
#print vel_tot_mag

end = time.time()
print 'It took just '+ str(end - start) + ' seconds!'
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
col = np.zeros((total_num))
col = np.array(np.where(gam<0,'r','g'))
#print col
ax0 = fig.add_subplot(611).set_xlabel('X')
ax0 = fig.add_subplot(611).set_ylabel('Y')
ax0 = fig.add_subplot(611).set_title('Vortices injection from X=0 using Normal distribution (not to scale, vortices size is 25 times the original)')
ax0 = plt.scatter(q[:,0],q[:,1], c=col[:,0],s=25*core_radius, label='vortex')
legend_elements = [Line2D([0], [0], marker='o', color='w', label='+ Vortex',
                          markerfacecolor='g', markersize=7),
                   Line2D([0], [0], marker='o', color='w', label='- Vortex',
                          markerfacecolor='r', markersize=7)]
ax0 = fig.add_subplot(611).legend(handles=legend_elements,loc='upper right')

ax1 = fig.add_subplot(612)
ax1 = fig.add_subplot(612).set_xlabel('X')
ax1 = fig.add_subplot(612).set_ylabel('Y')
ax1 = fig.add_subplot(612).set_title('Histogram for vortices injection from X=0')
ax1 = plt.hist(q[:,1], orientation='horizontal')

ax2 = fig.add_subplot(613)
ax2 = fig.add_subplot(613).set_xlabel('time')
ax2 = fig.add_subplot(613).set_ylabel('induced velocity')
ax2 = fig.add_subplot(613).set_title('induced velocity with respect to time')
ax2 = plt.plot(t, vel1[:,0], c='g', label='u (x-component of induced velocity)')
ax2 = plt.plot(t, vel1[:,1], c='r', label='v (y-component of induced velocity)')
ax2 = plt.plot(t, vel_tot_mag[:], c='k', label='U (total induced velocity)')
ax2 = fig.add_subplot(613).legend(loc='best')

ax3 = fig.add_subplot(614)
ax3 = fig.add_subplot(614).set_ylabel('induced velocity')
ax3 = fig.add_subplot(614).set_xlabel('Frequency (Hz)')
ax3 = fig.add_subplot(614).set_title('FFT of induced velocity')
ax3 = plt.semilogx(xf, 2.0/N * np.abs(yf[:N//2]))

ax4 = fig.add_subplot(615)
ax4 = fig.add_subplot(615).set_ylabel('Power Spectral Density (dB)')
ax4 = fig.add_subplot(615).set_xlabel('Frequency (Hz)')
ax4 = fig.add_subplot(615).set_title('Power Spectral Density (PSD) of induced velocity')
ax4 = plt.semilogx(freq,pegel,linewidth=0.6)

ax5 = fig.add_subplot(616)
ax5 = fig.add_subplot(616).set_xlabel('X')
ax5 = fig.add_subplot(616).set_ylabel('Y')
ax5 = fig.add_subplot(616).set_title('Histogram for radius size of the vortices injected from X=0')
ax5 = plt.hist(core_radius[:], orientation='vertical')
fig.savefig('spatial_master_v1.2.pdf')
plt.show()

