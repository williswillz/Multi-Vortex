#!/usr/bin/env python
# coding: utf-8

# In[72]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import psd
import essentials


# ### NACA 4 digit airfoil generator

# In[73]:


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

# In[74]:


q = naca_airfoil(0012, 101)


# ### Collacation point (obervations point on the airfoil where the vortex influence will be captured)

# In[75]:


q_mid = (q[1:] + q[:-1])/2


# ### Defining the computational domain where the vortices are distributed

# In[76]:


length = 100. #nondimensional length of window
height = 0.1 #window height
N = 10000 #number of vortices
gammas = 1. #vortex strength RMS (normal distribution)
rscale = 0.1 #vortex size scale (rayleigh distribution parameter)
t0 = -10-1.#start time for observation of convection
t1 = 10+1.#end time
ts = 0.02 # time step
v0 = 5 #convection speed


# ### Set random distribution for vortex location, size and strength
# origin at window center

# In[77]:


vortx = np.random.uniform(low=-length/2,high=length/2,size=N)
vorty = np.random.uniform(low=-height/2,high=height/2,size=N)
vortX = np.vstack((vortx,vorty))
gamma = np.random.normal(scale=gammas,size=N)
rho = np.random.rayleigh(scale=rscale,size=N)


# In[78]:


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
rho = essentials.round_up(rho,1)


# ### Solver/time-stepper

# In[79]:


t = np.arange(t0,t1,ts) # number of time-steps
# print t
print len(t)


# In[80]:


A = [] # empty array to store cl value after each iteration
for i in range(len(t)):
    print "this is the iteration {} / {}".format(i, len(t))
    i = i * ts * v0
    obsx = q_mid[:,0] + ((length+1)/2)
    obsx = (obsx - i )
    obsy = q_mid[:,1]
    obsX = np.vstack((obsx,obsy))
    dist = obsX[:,:,np.newaxis]-vortX[:,np.newaxis,:] # dim 2 x timesteps x N
    r = np.sqrt((dist*dist).sum(0)) # dim timesteps x N
    utheta = 16 * gamma * (rho**(-3)) * np.exp(-8*(rho**(-4)) * r**2) * (3-(16 * (rho**(-4)) * r**2)) * r   # Mexican-hat shape
    uind = utheta * dist[::-1] # dim 2 x timesteps x N
    uind[0] *= -1 # change sign for ux (to get correct rotation)
    ux = uind[0].T
    utot = uind.sum(2) # dim 2 x timesteps
    q_newx = q[:,0] + (length/2)
    q_newx = q_newx - i
    q_newy = q[:,1]
    q = np.array([q_newx, q_newy]).T
    dq = np.diff(q, axis=0)
    numpanels = dq.shape[0]
    lengths = np.linalg.norm(dq, axis=1) 
    normals = np.transpose(np.array([dq[:,1], -dq[:,0]]) / lengths)
    tangents = -np.transpose(np.array([dq[:,0], dq[:,1]]) / lengths)
    utot_tangent = utot.T * tangents
    utot_tangent_magnitude = pow((pow(utot_tangent[:,0],2) + pow(utot_tangent[:,1],2)), 0.5)
    p_ref =0
    p = p_ref + (0.5 * 1.225 * (v0**2 - ux**2))
    cp = (1 - ((utot_tangent_magnitude**2)/v0**2))
#     cp = (1 - ((ux**2)/v0**2))
    cp = cp * normals[:,1]
    cl = cp.sum(0)
    A.append(cl)


# #### Unsteady lift

# In[81]:


print p
dp_ = np.diff(p, axis=0)
dp = dp_ * lengths
print dp
#np.savetxt('dp.csv', np.column_stack((dp)), fmt='%0.5f', delimiter=',')
dt = ts


# In[82]:


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

tau = dt - (r_/340)

p_acoustic = (((obs_X*(normals.T)).sum(0))/(4*np.pi*340)) * ((dp/(tau*(r_**2)))[np.newaxis])
#var_1 = (dp/(tau*(r_**2)))[np.newaxis]
p_acoustic_1 = p_acoustic[0,:,:].sum(axis=1)

H_NC = ((np.abs(p_acoustic_1))/(2.e-5))**2
# print H
(val_NC, freq_NC) = psd(H_NC,Fs=1/dt,detrend='mean')
SPL = 10*np.log10(val_NC)
#savetxt('SPL.csv', np.column_stack((freq,SPL)), fmt='%0.5f', delimiter=',')
# print SPL
#plt.savefig('SPL.pdf')
#print freq
plt.figure(1)
plt.semilogx(freq_NC,10*np.log10(val_NC),linestyle='-',alpha=1.0,label='SPL')
plt.title('Sound Pressure Level (dB)')
plt.xlabel('Frequency')
plt.ylabel('SPL, dB',labelpad=1.5)
plt.grid(color='0.5',linestyle=':',linewidth=0.2)
plt.legend()
plt.savefig('SPL_check.pdf')
plt.show


# In[83]:


A = np.array(A)


# ### Curle's analogy to obtain p' at a distance of 1 unit at theta=90 deg

# In[84]:


def Curles_loadingNoise(y_int,c_sound,r_dist,L,dt,Velo):
#     p_acoustic = (((y_int*L)/(4*np.pi*dt*c_sound*(r_dist**2)))*(0.5*1.225*pow(Velo,2)))
    p_acoustic = (((y_int*L)/(4*np.pi*dt*c_sound*(r_dist**2))))

    return p_acoustic

noise = Curles_loadingNoise(1,343,1,np.abs(A),ts,v0)
H = ((noise)/(2.e-5))**2


# ### Plot unsteady lift on the airfoil and the SPL in far-field 

# In[85]:


plt.figure(2)
plt.subplot(1,2,1)
plt.plot(t,A,label='cl')
plt.title('Unsteady coefficient of lift')
plt.xlabel('time, s')
plt.ylabel('cl')
plt.legend()
plt.subplot(1,2,2)
(val, freq) = psd(H,Fs=1/ts,detrend='mean')
plt.semilogx(freq,10*np.log10(val),label='SPL')
plt.title('Sound Pressure Level (dB)')
plt.xlabel('Frequency')
plt.ylabel('SPL, dB')
plt.legend()
plt.savefig('cl_SPL.pdf')
plt.show


# In[ ]:
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 19:28:10 2019

@author: WS1
"""

#Amiet formulation
#Amiet formulation with Gershfeld's correction

from pylab import *
from numpy import *

# Beschriftungseigenschaften
ts=8
rc('font', family='sans-serif',size=ts)
rc('axes', labelsize=ts)
rc('xtick',labelsize=ts)
rc('ytick',labelsize=ts)
rc('legend',fontsize=ts,handlelength=1.6,labelspacing=0.3,handletextpad=0.5,borderpad=0.3)
width=8.
height=6.
ymin,ymax=-20,120

collist=['b','g','r','c','m']



###################
##     Berechnung LEN
###################

# Amiet
def calcAmiet(bx,Rx,Ux,ux,Lambdax,tfreqx,cx):
	Mx=Ux/cx
	erg_ami=zeros(len(tfreqx))
	zaehler=0
	for fx in tfreqx:
		Kx=8.*pi*fx*Lambdax/(3.*Ux)
		erg_ami[zaehler]=10.*log10((bx*Lambdax*Mx**5.*ux**2.*Kx**3.)/(Rx**2.*Ux**2.*(1+Kx**2.)**(7./3.)))+181.3
		zaehler=zaehler+1
	return(erg_ami)

# Amiet + Gershfeld
def calcAmietGershfeld(bx,dx,Rx,Ux,ux,Lambdax,tfreqx,cx):
	Mx=Ux/cx
	erg_ami=zeros(len(tfreqx))
	zaehler=0
	for fx in tfreqx:
		Kx=8.*pi*fx*Lambdax/(3.*Ux)
		erg_ami[zaehler]=10.*log10(exp(-2.*pi*fx*dx/(2.*Ux))*(bx*Lambdax*Mx**5.*ux**2.*Kx**3.)/(Rx**2.*Ux**2.*(1+Kx**2.)**(7./3.)))+181.3
		zaehler=zaehler+1
	return (erg_ami)


#####################
##        Tu
#####################
chord=1
span=5
U=5.0
Lambda=0.0058
freq=array((10,12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,3150,4000,5000,6300,8000,10000,12500,16000,20000))

Tulist=[3.8]

labellist=[r'$Tu$ = 3.8%']


figure(2,figsize=(width/2.54,height/2.54))
for Tu in Tulist:
	j=Tulist.index(Tu)
	
	b=span/2.
	d=0.12*chord
	pegelamiet=calcAmietGershfeld(b,d,1.0,U,Tu*U,Lambda,freq,343.)
	pegelamiet2=calcAmiet(b,1.0,U,Tu*U,Lambda,freq,343.)
	
	print max(pegelamiet), Tu

	semilogx(freq,pegelamiet,'-',linewidth=1.0,color=collist[j],label=labellist[j])
	semilogx(freq,pegelamiet2,'--',linewidth=1.0,color=collist[j])
semilogx(freq_NC,10*np.log10(val_NC),linestyle='-',alpha=1.0,label='SPL',color='g')

grid(color='0.5',linestyle=':',linewidth=0.2)
# savetxt('Tu_038_Amiet+Gershfeld.csv',np.column_stack((pegelamiet, pegelamiet2)), fmt='%5s', delimiter=',')
# savetxt('Tu_038_Amiet.csv',np.column_stack((freq, pegelamiet, pegelamiet2)), fmt='%5s', delimiter=',')
xticks((20,50,100,200,500,1000,2000,5000,10000,20000),('0,02','0,05','0,1','0,2','0,5','1','2','5','10','20'))
xlim(20,20000)
ylim(ymin,ymax)
legend(loc='best')
xlabel('$f_m$ in kHz',labelpad=1.0)
ylabel(r'$L_{p}$ in dB',labelpad=1.0)
gca().set_position([0.135,0.13,0.81,0.84] )
# savefig('Flatplate_LE_Tu038_.pdf',dpi=600)
savefig('naca0012_U5.pdf',dpi=600)
plt.show()



