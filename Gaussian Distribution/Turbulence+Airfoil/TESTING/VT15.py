# -*- coding: utf-8 -*-
"""
Created on Fri May  3 03:51:17 2019

@author: WS1


# ### Tasks

- discretize the vortex parameters (that means you will allow not any, but only certain - a finite number of different - values for y0 and size)
- verify that we still get the right spectra
- compute the farfield sound pressure for any combination of these
- put together the farfield sound pressures for a certain vortex distribution using proper delays and scaling on the individual farfield sound pressures (this is easily done for 10000 vortices if needed)
- compare the results to reference (Amiet, Gershfeld, HiFi methods)
"""
# ### Libraries
# In[2]:


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

# ### Steady background flow around the airfoil
# #### Solved using the steady vortex particle panel code 'steapy'.

# In[3]:


# Parameters to define a NACA 4 digits airfoil
x=1
m=0.03
p=0.4
t=0.12 # maximum thickeness
c=1

coordinates_symmetrical(x, t, c)
#coordinate = coordinates_cambered(x, m, p, t, c)

'''
Comment or Incomment the file path for the specific type of airfoil

'''
# load geometry from data file
naca_filepath = os.path.abspath('NACA'+'_'+'00'+str(int(t*100))+'.txt') #SYMMETRICAL
#naca_filepath = os.path.abspath('NACA'+'_'+str(int(m*100))+str(int(p*10))+str(int(t*100))+'.txt') #CAMBERED
with open(naca_filepath, 'r') as infile:
    x, y = numpy.loadtxt(infile, dtype=float, unpack=True)

# discretize geoemetry into panels
panels = define_panels(x, y, N=40)

'''DEFINE VELOCITY AND ANGLE OF ATTACK'''
# define freestream conditions
freestream = Freestream(u_inf=1.0, alpha=0.0)

#solve the source contribution
A_source = source_contribution_normal(panels)

#solve the vortex contribution
B_vortex = vortex_contribution_normal(panels)

#build LHS and RHS
A = build_singularity_matrix(A_source, B_vortex)
b = build_freestream_rhs(panels, freestream)

# solve for singularity strengths
strengths = numpy.linalg.solve(A, b)

# store source strength on each panel
for i , panel in enumerate(panels):
    panel.sigma = strengths[i]
    
# store circulation density
gamma = strengths[-1]


# ### Plot the velocity field 
# #### (just for visualisation and to prove the methodology is right)

# In[4]:
"""
PLOTS THE STREAMLINE AROUND THE AIRFOIL OF INTEREST
"""

#define a mesh grid
#nx, ny = 50, 50  # number of points in the x and y directions
#x_start, x_end = -1.0, 2.0
#y_start, y_end = -0.5, 0.5
#X, Y = numpy.meshgrid(numpy.linspace(x_start, x_end, nx),
#                      numpy.linspace(y_start, y_end, ny))
#
##compute the velocity field on the mesh grid
#u, v = get_velocity_field(panels, freestream, X, Y)
#
##plot the velocity field
#width = 10
#
#pyplot.figure(figsize=(width, width))
##pyplot.figure(figsize=(20,20))
#pyplot.xlabel('x', fontsize=16)
#pyplot.ylabel('y', fontsize=16)
#pyplot.streamplot(X, Y, u, v,
#                 density=1, linewidth=1, arrowsize=1, arrowstyle='->')
#pyplot.fill([panel.xc for panel in panels],
#           [panel.yc for panel in panels],
#           color='k', linestyle='solid', linewidth=2, zorder=3)
#pyplot.axis('scaled', adjustable='box')
#pyplot.xlim(x_start, x_end)
#pyplot.ylim(y_start, y_end)
##pyplot.ylim(1,-1)
##pyplot.title('Streamlines around a NACA 0012 airfoil (AoA = ${}^o$)'.format(alpha), fontsize=16);
## savefig('Velo_NACA'+'_'+str(int(m*100))+str(int(p*10))+str(int(t*100))+'.pdf')
#show()


# # Tracer particle coordinates
# ### Declaring a variable to store cooridnates of airfoil

# In[6]:


q = essentials.naca_airfoil(0012, 101)


# ### Collacation point (obervations point on the airfoil where the vortex influence will be captured)

# In[7]:


q_mid = (q[1:] + q[:-1])/2


# ### Defining the computational domain where the vortices are distributed

# <img src="tracer.JPG" alt="Drawing" style="width: 700px;"/>

# ### Injecting vortex from x = -1.0, y = random (y_start, y_end)
# origin at window center

# In[8]:


import numpy as np
x_start = -1.0
y_start, y_end = -0.1, 0.1
ylim = y_start, y_end

#function to give a random number between y_lim
def random_y(ylim):
    yrange = np.diff(ylim)
    return yrange * np.random.rand(1)[0] + ylim[0]

random_y(ylim)
# print (random_y)

#pts_y = np.arange(-1,1,0.1)
pts_y = np.arange(-1,1,0.1)
pts_x = (np.ones_like(pts_y))*-1
#pts_y = [0.05]
#pts_x = -1
pts =  np.vstack((pts_x,pts_y))
pts = pts[np.newaxis]
#pts = [] # empty list
#pts = list(pts)
#pts=((x_start, random_y(ylim)[0]))
##pts = (-1, 0.05)
## pts=((x_start, 0))
## print pts
#pts = np.asarray(pts)
x = np.zeros((len(pts_y)))
y = np.zeros((len(pts_y)))
x[:], y[:] = np.asarray(pts[0,:,:].T).transpose()
print (x,y)


# In[8]:

t0 = x_start  #start time for observation of convection
t1 = x_start + 3   #end time
dt = 0.01 # time step
t = np.arange(t0,t1,dt) # number of time-steps
a = np.zeros((len(pts_y),len(t)))
b = np.zeros((len(pts_y),len(t)))
#xx = np.zeros((len(pts_y),len(t)))
#yy = np.zeros((len(pts_y),len(t)))

for i in range(0,len(t)):
    vel = (((np.asarray(get_velocity_field(panels, freestream, x[:], y[:]))))[np.newaxis]).T
#    print (vel)
    pts[0,:,:] = pts[0,:,:] + ((vel*dt).T)
#    pts[:,:,0] = pts[:,:,0] + (vel*dt)
    x[:], y[:] = np.asarray(pts[0,:,:].T).transpose()
    
#    x, y = np.asarray(pts).transpose()
#     print (pts)
    a[:,i] = x[:]
    b[:,i] = y[:]
    
#vortX = np.vstack((a,b)))
vortX = a
vortX = vortX[np.newaxis,:,:]
vortY = b
vortY = vortY[np.newaxis,:,:]
vortZ = np.vstack((vortX,vortY))
#np.append(vortX[:,:,:], b, axis=0) 
#vortX = vortX[:,np.newaxis,:]
#print vortX
'''
plt.figure(figsize=(15, 5))
plt.plot(vortX[11,:],vortX[31,:],label='streamline')
#plt.plot(a[9,:],b[9,:],label='streamline')
pyplot.fill([panel.xc for panel in panels],
           [panel.yc for panel in panels],
           color='k', linestyle='solid', linewidth=2, zorder=3)
plt.title('tracer particle',fontsize=16)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
pyplot.axis('scaled', adjustable='box')
#pyplot.xlim(-1,2)
pyplot.ylim(-1.2,1.2)
plt.savefig('tracer_.pdf')
plt.show
'''

# ### Calculating distance between the vortex and each panel on airfoil surface

# In[20]:


obsx = q_mid[:,0]
# obsx = (obsx - i )
obsy = q_mid[:,1]
obsX = np.vstack((obsx,obsy))
# print obsX
#dist = obsX[:,:,np.newaxis]-vortX[:,np.newaxis,:] # dim 2 x timesteps x N
obsX = obsX[:,np.newaxis,:]
obsX = np.hstack(([obsX]*len(pts_y)))
obsX = obsX[:,:,:,np.newaxis]
vortZ = vortZ[:,:,np.newaxis,:]
dist = obsX-vortZ# dim 2 x timesteps x N
r = np.sqrt((dist*dist).sum(0)) # dim timesteps x N
print r
with h5py.File('r.h5', 'w') as hf:
    hf.create_dataset("dataset",  data=r)
#with h5py.File('aaa.h5', 'w') as hf:
#    hf.create_dataset("dataset",  data=aaa)


# ### Vortex strength and size

# In[10]:


gammas = 1. #vortex strength RMS (normal distribution)
rscale = 0.5 #vortex size scale (rayleigh distribution parameter)
gamma = np.random.normal(scale=gammas)
rho = np.random.rayleigh(scale=rscale, size = 10)
print gamma, rho


# ### Velocity calculation

# In[11]:

#gamma = 1.
#gamma = 1.035260581
#gamma = 10
#rho = np.array([0.9])
#r = aaa**0.5
r = r[np.newaxis]
r = np.vstack(([r]*len(rho)))
rho = (rho[:,np.newaxis]).T
savetxt('rho.csv',np.column_stack(rho), fmt='%0.5f', delimiter=',')

utheta = 16 * gamma * (rho[:,0]**(-3)) * np.exp(-8*(rho[:,0]**(-4)) * r[:,:,:,:]**2) * (3-(16 * (rho[:,0]**(-4)) * r[:,:,:,:]**2)) * r[:,:,:,:]   # Mexican-hat shape
# utheta = format(utheta)
# print utheta
with h5py.File('utheta.h5', 'w') as hf:
    hf.create_dataset("dataset",  data=utheta)
#savetxt('utheta.txt',np.column_stack(utheta), fmt='%0.5f', delimiter=',')

dist = dist[:,np.newaxis,:,:,:]
dist = np.vstack(([dist]*len(rho[0,:])))
uind = utheta * dist # dim 2 x timesteps x N
# print uind
uind[0,:,:,:,:] *= -1 # change sign for ux (to get correct rotation)
ux = uind[0,:,:,:,:].T
# print uind
# print uind.T
#with h5py.File('uind.h5', 'w') as hf:
#    hf.create_dataset("dataset",  data=uind)
with h5py.File('ux.h5', 'w') as hf:
    hf.create_dataset("dataset",  data=ux)
#savetxt('utheta_mag_i_j.txt',np.column_stack((utheta,uind[0],uind[1])), fmt='%0.5f', delimiter=',')
# utot = uind.sum(2) # dim 2 x timesteps
# print utot


# ### Pressure calculation using Bernoulli's equation

# In[22]:

'''update velocity here to to get correct bernoulli's pressure'''
v0 = 1 
q_newx = q[:,0]
q_newy = q[:,1]
q = np.array([q_newx, q_newy]).T
dq = np.diff(q, axis=0)
numpanels = dq.shape[0]
lengths = np.linalg.norm(dq, axis=1) 
# print lengths
normals = np.transpose(np.array([dq[:,1], -dq[:,0]]) / lengths)
tangents = -np.transpose(np.array([dq[:,0], dq[:,1]]) / lengths)
# print tangents
#utot_tangent = uind.T * tangents #old
#utheta_ = utheta[np.newaxis] #new
#utot_tangent = utheta_.T * tangents #new
# utot_tangent = uind.T
# print utot_tangent
# print utot_tangent[:,:,0]
#utot_tangent_magnitude = pow((pow(utot_tangent[:,:,0],2) + pow(utot_tangent[:,:,1],2)), 0.5)
# print utot_tangent_magnitude
p_ref =0
#p = p_ref + (0.5 * 1.225 * (1**2 - utot_tangent_magnitude**2))
'''taking uind[0] which is x component velocity'''
p = p_ref + (0.5 * 1.225 * (v0**2 - ux**2))
#with h5py.File('p.h5', 'w') as hf:
#    hf.create_dataset("dataset",  data=p)
#savetxt('p.csv',np.column_stack(p), fmt='%0.5f', delimiter=',')
dp_ = np.diff(p, axis=0)
lengths = lengths[:,np.newaxis]
dp = dp_.T * lengths
#dp = dp*lengths
# print len(p)
#F = dp_ * lengths
# print F
#G = F.sum(axis=1) # dim 2 x timesteps
# print len(G)
#cp = (1 - ((utot_tangent_magnitude**2)/v0**2))
#cp = cp * normals[:,1]
# print cp


# In[1]:
'''
TIME DOMAIN
'''
obs_x = 0
# obsx = (obsx - i )
obs_y = 1
obs_X = np.vstack((obs_x,obs_y))
# print obsX
vort_X = np.vstack((q_mid[:,0],q_mid[:,1]))
#print vort_X
dist_ = obs_X[:,:,np.newaxis]-vort_X[:,np.newaxis,:] # dim 2 x timesteps x N
r_ = np.sqrt((dist_*dist_).sum(0)) # dim timesteps x N
#print r_
r_cap = dist_[:,0,:] / r_
tau = dt - (r_/340)

#p_acoustic = (((obs_X*(normals.T)).sum(0))/(4*np.pi*340)) * ((dp/(tau*(r_**2)))[np.newaxis]) wrong r_cap
#normalst = normals.T
dp = dp.T
nr = ((dist_[:,0,:]*normals.T).sum(0))
nr = nr.T
nr = nr[np.newaxis]
nr = np.vstack(([nr]*len(pts_y)))
nr = nr[np.newaxis]
#nr = np.vstack(([nr]*len(rho)))
nr = np.vstack(([nr]*len(rho[0,:])))
nr = nr.T

tau_ = (tau*(r_**2))
tau_ = tau_.T
tau_ = tau_[np.newaxis]
tau_ = np.vstack(([tau_]*len(pts_y)))
tau_ = tau_[np.newaxis]
tau_ = np.vstack(([tau_]*len(rho[0,:])))
tau_ = tau_.T
#p_acoustic = (((dist_[:,0,:]*normals.T).sum(0))/(4*np.pi*340)) * ((dp/(tau*(r_**2)))[np.newaxis])
p_acoustic = ((nr)/(4*np.pi*340)) * ((dp/tau_)[np.newaxis])
with h5py.File('p_acoustic.h5', 'w') as hf:
    hf.create_dataset("dataset",  data=p_acoustic)
#var_1 = (dp/(tau*(r_**2)))[np.newaxis]
p_acoustic_1 = p_acoustic[0,:,:,:,:].sum(axis=1)
with h5py.File('p_acoustic_1.h5', 'w') as hf:
    hf.create_dataset("dataset",  data=p_acoustic_1)
savetxt('p_acoustic_1.csv', np.column_stack((p_acoustic_1)), fmt='%0.5f', delimiter=',')


p_acoustic_2 = ((np.abs(p_acoustic_1))/(2.e-5))**2
with h5py.File('p_acoustic_2.h5', 'w') as hf:
    hf.create_dataset("dataset",  data=p_acoustic_2)
savetxt('p_acoustic_2.csv', np.column_stack((p_acoustic_2)), fmt='%0.5f', delimiter=',')
    
    # In[10]:
(val_NC, freq_NC) = psd(H_NC,Fs=1/dt,detrend='mean')
SPL = 10*np.log10(val_NC)
#savetxt('SPL.csv', np.column_stack((freq,SPL)), fmt='%0.5f', delimiter=',')
# print SPL
#plt.savefig('SPL.pdf')
#print freq
plt.semilogx(freq_NC,10*np.log10(val_NC),label='SPL')
plt.title('Sound Pressure Level (dB)')
plt.xlabel('Frequency')
plt.ylabel('SPL, dB',labelpad=1.5)
grid(color='0.5',linestyle=':',linewidth=0.2)
plt.legend()
#plt.savefig('SPL1.pdf')
plt.show

# In[14]:
'''
FREQUENCY DOMAIN
'''
#sp = numpy.fft.fft(dp_)
#p1 = p[:,np.newaxis,:]
#p2 = p1[:,:,:] * r_cap
#p2 = p1 * r_cap
p3 = ((dist_[:,0,:]*normals.T).sum(0))
#p4 = p*p3.T
p4 = p*nr
'''
important
'''
sp = numpy.fft.fft(p4) 
sp = numpy.delete(sp, (0), axis=0)
lengths_fft = lengths[np.newaxis]
freqfft = np.fft.fftfreq(len(t)-1)



mat = (freqfft * r_.T).T
mat1 = (2*np.pi/340)*(-1J)
mat2 = np.exp(mat1*mat)
mat2 = mat2.T #exp(-jkr)
mat2 = mat2[np.newaxis]
mat2 = np.vstack(([mat2]*len(pts_y)))
mat2 = mat2[np.newaxis]
mat2 = np.vstack(([mat2]*len(rho[0,:])))
mat2 = mat2.T
#
#mat3 = np.multiply(sp,mat2)
#mat4 = mat3.sum(axis=1)
#
#mat5 = ((freqfft / r_.T).T) / (2*340)
#mat6 = np.multiply(mat5,mat3)
#mat7 = mat6.sum(axis=1)
#mat8 = mat7 * (-1J)
mat9 = ((1J*freqfft) / (((r_.T)**2) * 340 * 2))
mat9 = mat9[np.newaxis]
mat9 = np.vstack(([mat9]*len(pts_y)))
mat9 = mat9[np.newaxis]
mat9 = np.vstack(([mat9]*len(rho[0,:])))
mat9 = mat9.T


mat10 = np.multiply(mat9,mat2,sp)

lengths_fft = lengths_fft[:,:,:,np.newaxis] #cis it possile to achieve this wwithout multplying by len(rho) and len(pts_y)?
mat14 = mat10*lengths_fft
mat11 = mat14.sum(axis=1)

mat12 = np.fft.ifft(mat11)
L = ((np.abs(mat12.real))/(2.e-5))**2

(val_L, freq_L) = psd(L,Fs=1/dt,detrend='mean')
plt.semilogx(freq_L,10*np.log10(val_L),label='SPL')
# In[14]:


print G
print len(G)


# ### Curle's analogy to obtain p' at a distance of 1 unit at theta=90 deg

# In[15]:


def Curles_loadingNoise(y_int,c_sound,r_dist,L,dt,Velo):
#     p_acoustic = (((y_int*L)/(4*np.pi*dt*c_sound*(r_dist**2)))*(0.5*1.225*pow(Velo,2)))
    p_acoustic = (((y_int*L)/(4*np.pi*dt*c_sound*(r_dist**2))))
    return p_acoustic

noise = Curles_loadingNoise(1,343,1,G,dt,v0)
# print noise
H = ((noise)/(2.e-5))**2
# print H
(val, freq) = psd(H,Fs=1/dt,detrend='mean')
SPL = 10*np.log10(val)
savetxt('SPL.csv', np.column_stack((freq,SPL)), fmt='%0.5f', delimiter=',')
# print SPL
plt.savefig('SPL.pdf')
print freq


# In[16]:


plt.semilogx(freq,10*np.log10(val),label='SPL')
plt.title('Sound Pressure Level (dB)')
plt.xlabel('Frequency')
plt.ylabel('SPL, dB')
plt.legend()
plt.savefig('SPL1.pdf')
plt.show


# In[221]:


plt.plot(-G/(0.5*1.225*v0**2*lengths.sum()))
plt.title('Coefficient of lift')
plt.xlabel('n')
plt.ylabel('cl')
plt.legend()
# plt.savefig('SPL1.pdf')
plt.show


# In[ ]:

plt.figure(2)
#plt.subplot(1,2,1)
plt.semilogx(freq_NC,10*np.log10(val_NC),label='SPL_time domain',color='green') # compact source
plt.semilogx(freq,10*np.log10(val),label='SPL_time domain-compact source') # non-compact source
plt.semilogx(freq_L,10*np.log10(val_L),label='SPL_frequency domain')
plt.legend()
savefig('acoustic.pdf')
#plt.subplot(1,2,2)
#(valu,freq) = psd(utot[0],Fs=1/ts,detrend='mean')
#(valv,freq) = psd(utot[1],Fs=1/ts,detrend='mean')
#plt.loglog(freq[1:],valu[1:],label='u')
#plt.loglog(freq[1:],valv[1:],label='v')
#plt.loglog(freq[1:],valu[1:]+valv[1:],label='tot')
#plt.legend()


