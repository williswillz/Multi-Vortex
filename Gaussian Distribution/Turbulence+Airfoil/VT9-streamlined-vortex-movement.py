#!/usr/bin/env python
# coding: utf-8

# ### Tasks

# - compute the steady background flow around an airfoil using the panel method you have implemented earlier
# - compute the pathlines a massless particle realeased well upstream the airfoil at certain y0 will follow (these are the streamlines) until well downstream the airfoil. I am not totally sure what would be the best approach to do so, but has to be repeated for all relevant y0
# - create an acoustic 1d mesh at the airfoil surface (maybe you could use just the panels if the dimensions make sense)
# - define a sampling frequency and move a vortex along the pathline/streamline, compute the time-depended pressure from the vortex at all surface elements
# - integrate over the surface elements to get the farfield sound pressure (Curle) at any generic location (field point x_p,y_p). This involves also to take into account the retarded time instants due to different distances between the panels and the field point. One option would be to interpolate, another option would be to take the FFT of p at each surface element and multiply by exp(-jkr). The latter is more easy to implement (frequency domain Curle's formulation). Finally you would have the farfield sound pressure due to 1 vortex (strength,size) that started at y0.
# - repeat everything for different y0 and also for same y0, but different vortex size
# - verify that the results are different and do not just scale with a factor
# - repeat everything for same y0 and different strength
# - verify that the results are the same, except for a certain factor

# ### Libraries

# In[ ]:


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
import sympy
from sympy.abc import x, y


# ### Steady background flow around the airfoil
# #### Solved using the steady vortex particle panel code 'steapy'.

# In[47]:


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

# In[48]:


#define a mesh grid
nx, ny = 50, 50  # number of points in the x and y directions
x_start, x_end = -5.0, 2.0
y_start, y_end = -0.5, 0.5
X, Y = numpy.meshgrid(numpy.linspace(x_start, x_end, nx),
                      numpy.linspace(y_start, y_end, ny))

#compute the velocity field on the mesh grid
u, v = get_velocity_field(panels, freestream, X, Y)

#plot the velocity field
width = 10

pyplot.figure(figsize=(width, width))
#pyplot.figure(figsize=(20,20))
pyplot.xlabel('x', fontsize=16)
pyplot.ylabel('y', fontsize=16)
pyplot.streamplot(X, Y, u, v,
                 density=1, linewidth=1, arrowsize=1, arrowstyle='->')
pyplot.fill([panel.xc for panel in panels],
           [panel.yc for panel in panels],
           color='k', linestyle='solid', linewidth=2, zorder=3)
pyplot.axis('scaled', adjustable='box')
pyplot.xlim(x_start, x_end)
pyplot.ylim(y_start, y_end)
pyplot.ylim(1,-1)
#pyplot.title('Streamlines around a NACA 0012 airfoil (AoA = ${}^o$)'.format(alpha), fontsize=16);
# savefig('Velo_NACA'+'_'+str(int(m*100))+str(int(p*10))+str(int(t*100))+'.pdf')
show()


# # Tracer particle coordinates

# ### NACA 4 digit airfoil generator

# In[4]:


def naca_airfoil(code, num_points, zero_thick_te=False, uniform=False):
    """Return a NACA 4-digit series airfoil"""
    # extract parameters from 4-digit code
    code_str = "%04d" % int(code)
    if len(code_str) != 4:
        raise ValueError("NACA designation is more than 4 digits")
    max_camber = 0.01 * int(code_str[0])
    p = 0.1 * int(code_str[1])  # location of max camber
    thickness = 0.01 * int(code_str[2:])
    if uniform:
        x = np.linspace(0, 1, num_points)
    else:
        # closer spacing near leading edge
        theta = np.linspace(0, 0.5 * np.pi, num_points)
        x = 1 - np.cos(theta)

    # thickness
    coefs = [-0.1015, 0.2843, -0.3516, -0.1260, 0, 0.2969]
    if zero_thick_te:
        coefs[0] = -0.1036
    y_thick = 5 * thickness * (np.polyval(coefs[:5], x) +
                               coefs[5] * np.sqrt(x))

    # camber
    front = np.where(x <= p)
    back = np.where(x > p)
    y_camber = np.zeros_like(x)
    if p:
        y_camber[front] = max_camber * x[front] / p**2 * (2 * p - x[front])
        y_camber[back] = max_camber * ((1. - x[back])/(1. - p)**2 *
                                       (1 + x[back] - 2 * p))
    x = np.hstack([x[-1:0:-1], x])
    y = np.hstack([y_camber[-1:0:-1] + y_thick[-1:0:-1],
                   y_camber - y_thick])
    return (np.array([x, y]).T)


# ### Declaring a variable to store cooridnates of airfoil

# In[5]:


q = naca_airfoil(0012, 11)


# ### Collacation point (obervations point on the airfoil where the vortex influence will be captured)

# In[6]:


q_mid = (q[1:] + q[:-1])/2


# ### Defining the computational domain where the vortices are distributed

# <img src="tracer.JPG" alt="Drawing" style="width: 700px;"/>

# ### Injecting vortex from x = -1.0, y = random (y_start, y_end)
# origin at window center

# In[51]:


import numpy as np
x_start, x_end = -1.0, 2.0
y_start, y_end = -0.5, 0.5
ylim = y_start, y_end

#function to give a random number between y_lim
def random_y(ylim):
    yrange = np.diff(ylim)
    return yrange * np.random.rand(1)[0] + ylim[0]

random_y(ylim)
# print (random_y)


pts = [] # empty list
pts = list(pts)
pts=((x_start, random_y(ylim)[0]))
# print pts
pts = np.asarray(pts)
x, y = np.asarray(pts).transpose()
# print (x,y)

n = 300
a = np.zeros((1,n))
b = np.zeros((1,n))

for i in range(0,n):

    dt = 0.01
#    u, v = get_velocity_field(panels, freestream, x, y)
#    print (u,v)

    vel = (np.asarray(get_velocity_field(panels, freestream, x, y))).T
#    print (vel)
    pts = pts + (vel*dt)
    x, y = np.asarray(pts).transpose()
#     print (pts)
    a[0,i] = x
    b[0,i] = y
    
vortX = np.vstack((a,b))
print vortX


# ### Calculating distance between the vortex and each panel on airfoil surface

# In[56]:


obsx = q_mid[:,0]
# obsx = (obsx - i )
obsy = q_mid[:,1]
obsX = np.vstack((obsx,obsy))
# print obsX
dist = obsX[:,:,np.newaxis]-vortX[:,np.newaxis,:] # dim 2 x timesteps x N
r = np.sqrt((dist*dist).sum(0)) # dim timesteps x N
# print r


# ### Vortex strength and size

# In[71]:


gammas = 1. #vortex strength RMS (normal distribution)
rscale = 0.1 #vortex size scale (rayleigh distribution parameter)
gamma = np.random.normal(scale=gammas)
rho = np.random.rayleigh(scale=rscale)
print gamma, rho


# ### Velocity calculation

# In[72]:


# gamma = 100
# rho = 40
utheta = 16 * gamma * (rho**(-3)) * np.exp(-8*(rho**(-4)) * r**2) * (3-(16 * (rho**(-4)) * r**2)) * r   # Mexican-hat shape
# utheta = format(utheta)
# print utheta
savetxt('utheta.txt',np.column_stack(utheta), fmt='%0.5f', delimiter=',')
uind = utheta * dist[::-1] # dim 2 x timesteps x N
# print uind
uind[0] *= -1 # change sign for ux (to get correct rotation)
uindT = uind.T
# print uind
print uind.T
savetxt('utheta_mag_i_j.txt',np.column_stack((utheta,uind[0],uind[1])), fmt='%0.5f', delimiter=',')
# utot = uind.sum(2) # dim 2 x timesteps
# print utot


# ### Pressure calculation using Bernoulli's equation

# In[73]:


v0 = 1 # empty array to store cl value after each iteration
q_newx = q[:,0]
q_newy = q[:,1]
q = np.array([q_newx, q_newy]).T
dq = np.diff(q, axis=0)
numpanels = dq.shape[0]
lengths = np.linalg.norm(dq, axis=1) 
normals = np.transpose(np.array([dq[:,1], -dq[:,0]]) / lengths)
tangents = -np.transpose(np.array([dq[:,0], dq[:,1]]) / lengths)
# print tangents
utot_tangent = uind.T * tangents
# print utot_tangent
# print utot_tangent[:,:,0]
utot_tangent_magnitude = pow((pow(utot_tangent[:,:,0],2) + pow(utot_tangent[:,:,1],2)), 0.5)
# print utot_tangent_magnitude
p_ref =0
p = p_ref + (0.5 * 1.225 * (1**2 - utot_tangent_magnitude**2))
# print p
F = p * lengths
# print F
G = F.sum(axis=1) # dim 2 x timesteps
# print G
cp = (1 - ((utot_tangent_magnitude**2)/v0**2))
cp = cp * normals[:,1]
# print cp


# ### Curle's analogy to obtain p' at a distance of 1 unit at theta=90 deg

# In[74]:


def Curles_loadingNoise(y_int,c_sound,r_dist,L,dt,Velo):
    p_acoustic = (((y_int*L)/(4*np.pi*dt*c_sound*(r_dist**2)))*(0.5*1.225*pow(Velo,2)))
    return p_acoustic

noise = Curles_loadingNoise(1,343,1,G,dt,v0)
# print noise
H = ((noise)/(2.e-5))**2
# print H
(val, freq) = psd(H, NFFT=512,Fs=1/dt,detrend='mean')
SPL = 10*np.log10(val)
# print SPL


# ### Plot unsteady lift on the airfoil and the SPL in far-field 

# In[75]:


cl = G / (0.5 * 1.225 * v0**2)
plt.figure(1)
plt.plot(G,label='cl')
plt.title('Unsteady coefficient of lift')
plt.xlabel('time, s')
plt.ylabel('cl')
plt.legend()
plt.savefig('cl_SPL.pdf')
plt.show


# In[ ]:




