from pysces import *
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from pylab import *
import time

start = time.time()

angle = 0
airfoil = naca_airfoil("0012", 501, zero_thick_te=True)  
airfoil = TransformedBody(airfoil)
airfoil = TransformedBody(airfoil, angle)
bound = BoundVortices(airfoil)

num_steps = 500
Uinfty = (1,0)
dt = 0.01
Vortices.core_radius = 0.005

#free vortices

#v1	=	(	1.724339281	,	3.1	)
v1 = (5.0,-0.1)
v2 = (5.0,0.1)
#s1	=	-0.00810137
s1 = 0.24
s2 = 0.24

vort = Vortices([v1,v2],[s1,s2])

flow = ExplicitEuler(dt, Uinfty, bound, wake=vort, need_force='wake_impulse')

for i in range(1,num_steps):
    flow.advance()

#plot force
f = flow.force
steps = np.arange(0, num_steps*dt, dt)
expected_Cl = (np.pi * np.pi / 90) * angle
expected = np.array([expected_Cl, expected_Cl])
steps2 = np.array([0, num_steps*dt])

#saving data
data = (steps, 2*f[:,0])  
savetxt('single-vortex_'+str(v1)+str(v2)+str(num_steps)+str(Uinfty)+str(dt)+str(s1)+str(s2)+'_'+'.csv',np.column_stack((steps, 2*f[:,1])), fmt='%5s', delimiter=',')
end = time.time()
print end - start

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(steps, 2*f[:,0], c='b', edgecolors='none', label='Cd')
ax1.scatter(steps, 2*f[:,1], c='r', edgecolors='none', label='Cl')
ax1.plot(steps2, expected, c='g', label='expected Cl')
plt.legend();
plt.xlabel('time')
#plt.grid(True)
plt.savefig('cl_cd_'+str(v1)+str(v2)+str(num_steps)+str(Uinfty)+str(dt)+str(s1)+str(s2)+'.pdf')
plt.show()

average_thrust = np.average(flow.force[:,0])
#print average_thrust

average_lift = np.average(flow.force[:,1])
#print average_lift

def Curles_loadingNoise(y_int,c_sound,r_dist,L,dt,Velo):
	p_acoustic = ((y_int*L)/(4*np.pi*dt*c_sound*(r_dist**2)))*(0.5*1.225*pow(Velo,2))
	return p_acoustic

noise = Curles_loadingNoise(1, 343,1, 2*f[:,1],dt,25.16)
#print (noise)

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
import math

def L_p(x):
	return 20*np.log10(np.abs(x)/2.e-5)

blocksize=100
j=0

(werte,freq)=plt.psd(noise, NFFT=blocksize, Fs=100, detrend=mlab.detrend_none,window=mlab.window_hanning, noverlap=4, pad_to=None,sides='default', scale_by_freq='True')
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
plt.savefig('SPL_single-vortex_'+str(v1)+str(v2)+str(num_steps)+str(Uinfty)+str(dt)+str(s1)+str(s2)+'.pdf')
plt.show()