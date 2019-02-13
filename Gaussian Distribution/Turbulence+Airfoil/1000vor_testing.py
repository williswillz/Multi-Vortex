# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 14:39:50 2019

@author: WS1
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 11:19:52 2019

@author: WS1
"""


# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import psd


# ## define parameters

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

q = naca_airfoil(0012, 101)
#print q
#print q[:,0]
#print q[:,1]
        
q_mid = (q[1:] + q[:-1])/2
#print q_mid

#dq = np.diff(q, axis=0)
#numpanels = dq.shape[0]
#lengths = np.linalg.norm(dq, axis=1) 
#normals = np.transpose(np.array([dq[:,1], -dq[:,0]]) / lengths)
#tangents = -np.transpose(np.array([dq[:,0], dq[:,1]]) / lengths)

length = 10. #nondimensional length of window
height = 0.2 #window height
v0 = 5 #convection speed
N = 10000 #number of vortices
gammas = 1. #vortex strength RMS (normal distribution)
rscale = 0.1 #vortex size scale (rayleigh distribution parameter)
t0 = -(length / v0 / 2)  #start time for observation of convection
t1 = (length / v0 / 2)   #end time
ts = 0.0001 # time step


# ## set random distribution for vortex location, size and strength
# origin at window center



vortx = np.random.uniform(low=-length/2,high=length/2,size=N)
vorty = np.random.uniform(low=-height/2,high=height/2,size=N)
vortX = np.vstack((vortx,vorty))
gamma = np.random.normal(scale=gammas,size=N)
rho = np.random.rayleigh(scale=rscale,size=N)


# ## set relative locations for observation
# vortex window moves to the right
# t=0 means observation in the center of the window



t = np.arange(t0,t1,ts)
delta = -v0*t
print len(t)
#obsy = np.zeros_like(obsx)
#obsX = np.vstack((obsx,obsy))

#for i in range(len(t)):

A = []
B = []
cp = []
for i in range(len(t)):
#    print i
    print "this is the iteration {} / {}".format(i, len(t))
    i = i * ts * v0
#    print i
    obsx = q_mid[:,0] + (length/2)
    obsx = (obsx - i )
#    print len(obsx)
    obsy = q_mid[:,1]
    obsX = np.vstack((obsx,obsy))
#    print obsX
    dist = obsX[:,:,np.newaxis]-vortX[:,np.newaxis,:] # dim 2 x timesteps x N
    r = np.sqrt((dist*dist).sum(0)) # dim timesteps x N
#    print r
    utheta = 16 * gamma * (rho**(-3)) * np.exp(-8*(rho**(-4)) * r**2) * (3-(16 * (rho**(-4)) * r**2)) * r   # Mexican-hat shape

    uind = utheta * dist[::-1] # dim 2 x timesteps x N
    uind[0] *= -1 # change sign for ux (to get correct rotation)
    utot = uind.sum(2) # dim 2 x timesteps
#    print utot.T
    q_newx = q[:,0] + (length/2)
    q_newx = q_newx - i
#    print len(q_newx)
    q_newy = q[:,1]
#    q = np.vstack((q_newx, q_newy))
    q = np.array([q_newx, q_newy]).T
#    print q
    dq = np.diff(q, axis=0)
#    print dq
    numpanels = dq.shape[0]
    lengths = np.linalg.norm(dq, axis=1) 
    normals = np.transpose(np.array([dq[:,1], -dq[:,0]]) / lengths)
#    print normals[:,0]
    tangents = -np.transpose(np.array([dq[:,0], dq[:,1]]) / lengths)
#    print tangents
    utot_tangent = utot.T * tangents
#    print utot_tangent
    utot_tangent_magnitude = pow((pow(utot_tangent[:,0],2) + pow(utot_tangent[:,1],2)), 0.5)
#    print utot_tangent_magnitude
    cp = (1 - ((utot_tangent_magnitude**2)/v0**2))
#    print cp
    cp = cp * normals[:,1]
#    print cp
    cl = cp.sum(0)
    print cl
    A.append(cl)
#    A.append(utot_tangent[:,0].sum())
#    B.append(utot_tangent[:,1].sum())
#    pref = 100000 #Pa
#    p = -(utot_tangent_magnitude**2)/2 + (v0**2)/2 + pref
#    print p
#    l = p * normals
#    print l
    
A = np.array(A)
#B = np.array(B)
#C = np.array((A**2 + B**2)**0.5)
#print C

plt.figure(1)
#plt.subplot(1,2,1)
#plt.plot(t,utot_tangent[:,0],label='u')
#plt.plot(t,utot_tangent[:,1],label='v')
plt.plot(t,A,label='cl')
#plt.plot(t,B,label='v')
#plt.plot(t,C,label='V')
plt.legend()
#plt.subplot(1,2,2)
#(valu,freq) = psd(utot[0],Fs=1/ts,detrend='mean')
#(valv,freq) = psd(utot[1],Fs=1/ts,detrend='mean')
#plt.loglog(freq[1:],valu[1:],label='u')
#plt.loglog(freq[1:],valv[1:],label='v')
#plt.loglog(freq[1:],valu[1:]+valv[1:],label='tot')
#plt.legend()
plt.show
#
pref = 100000 #Pa
#p = -(C**2)/2 + (v0**2)/2 + pref
#print p
def Curles_loadingNoise(y_int,c_sound,r_dist,L,dt,Velo):
	p_acoustic = (((y_int*L)/(4*np.pi*dt*c_sound*(r_dist**2)))*(0.5*1.225*pow(Velo,2)))
	return p_acoustic

noise = Curles_loadingNoise(10,343,10,np.abs(A),ts,v0)
print noise 
H = ((noise)/(2.e-5))**2
print H
#def L_p(x):
#	return 20*np.log10(np.abs(x)/2.e-5)
#spl = L_p(noise)
#print spl

plt.figure(2)
plt.subplot()
#(val,freq) = psd(noise,Fs=1/ts,detrend='mean')
(val, freq) = psd(H, NFFT=512,Fs=1000,detrend='mean')
#plt.psd(noise,Fs=1/ts,detrend='mean')
#plt.psd(H, NFFT=512,Fs=1000,detrend='mean')
plt.semilogx(freq,10*np.log10(val),label='SPL')
#print freq
#plt.loglog(freq,val,label='SPL')
plt.title('PSD: power spectral density')
plt.xlabel('Frequency')
plt.ylabel('PSD')
plt.tight_layout()
plt.legend()
plt.show

#def L_p(x):
#	return 20*np.log10(np.abs(x)/2.e-5)
#spl = L_p(noise)
#
#plt.figure(3)
#plt.subplot()
#(val,freq) = psd(noise,Fs=1/ts,detrend='mean')
##plt.semilogx(freq,spl,label='SPL')
#plt.semilogx(freq,20*np.log10(val/2.e-5),label='SPL')
#plt.title('PSD: power spectral density')
#plt.xlabel('Frequency')
#plt.ylabel('PSD')
#plt.tight_layout()
#plt.legend()
#plt.show


'''
def compute_force_pressure(bound_old, wake_old, bound_new, wake_new, pref=0, Uinfty=(1,0)): 
    """Finds the forces on the body by integrating the pressure (from unsteady Bernoulli equation)"""
    #q is the local velocity at collocation point 
    vel = _compute_velocity(bound_new, Uinfty, wake_new)
    q = np.linalg.norm(vel, axis=1)
    #vref is the kinematic velocity
    vref = _compute_kin_velocity(bound_new, Uinfty)
    if bound_old and wake_old:
        dPhi = _vortex_potential(bound_new, wake_new) - _vortex_potential(bound_old, wake_old) 
        dPhi = np.sum(dPhi)
    else:
        dPhi = 0 
    p = -(q**2)/2 + (vref**2)/2 - dPhi + pref
    Cp = 1 - ((2*q)**2) #factor of two due to singularity in velocity caluclation
    #note: look at Saffman ch. 2 for better calculation of velocity at surface?
    f = p * bound_new.lengths 
    motion = bound_new._body.get_motion()    
    if motion:
        norm_inertial = motion.map_vector(bound_new.normals)
    else:
        norm_inertial = bound_new.normals
    f = f[:,np.newaxis] * norm_inertial
    force = -np.sum(f, axis=0)
    return force, Cp 

obsx = q_mid[:,0]
obsy = q_mid[:,1]
#print obsx
#print obsy
obsX = np.vstack((obsx,obsy))
#print obsX


# ## Vortex models
# #### Gaussian shape vortex: $u_\theta = 18\Gamma  \rho^{-3} e^{-9 \rho^{-4} {\boldsymbol{{r}}}^{2}}r$
# #### Point vortex (singularity at the center): $u_\theta = \frac{\Gamma}{2 \pi r}$
# #### Zero velocity at the absolute center: $u_\theta = \frac{\Gamma r}{2\pi \rho^2}$
# #### Combination of point vortex and zero velocity at the center vortex : $minimum of$ $u_\theta = \frac{\Gamma}{2 \pi r}$ $and$ $u_\theta = \frac{\Gamma r}{2\pi \rho^2}$
# #### Lamb-Oseen vortex: $u_\theta = \frac{\Gamma}{2{\pi}r}\left ( 1-exp\left ( -\frac{r^2}{\rho^2} \right ) \right )$
# #### Mexican-Hat shape vortex: $u_\theta = 16 \Gamma r \rho^{-3}exp\left ( -\left ( \sqrt{8}\rho^{-2}r \right )^{2} \right )\left ( 3-\left ( 4\rho^{-2}r \right )^{2} \right )$



dist = obsX[:,:,np.newaxis]-vortX[:,np.newaxis,:] # dim 2 x timesteps x N
r = np.sqrt((dist*dist).sum(0)) # dim timesteps x N

#comment out one of the two following lines to get alternative vortex models:

#utheta = 18 * gamma * (rho**(-3)) * np.exp((-9*r**2) / (rho**4)) * r   # Gaussian shape function vortex

#utheta = (0.5/np.pi)*(1/r)                                             # Point vortex with singularity at the center

#utheta = (0.5/np.pi)*gamma*np.minimum(1/r,r/rho)                       # Wrong because r/rho is dimensionaless :dim timesteps x N:

#utheta = (0.5/np.pi)*gamma*np.minimum(1/r,r/rho**2)                    # Gaussian shape 2

#utheta = (0.5/np.pi)*gamma/r/r                                         # Wrong because it's divinding the function twice by 'r'

#utheta = gamma*rho**(1.5)*np.exp(-9*rho*rho*r*r)                       # Wrong

#utheta = (0.5/np.pi)*(1/r)*(1-np.exp(-(r**2)/rho**2))                  # Lamb-Oseen vortex

utheta = 16 * gamma * (rho**(-3)) * np.exp(-8*(rho**(-4)) * r**2) * (3-(16 * (rho**(-4)) * r**2)) * r   # Mexican-hat shape

# into cartesian coords

uind = utheta * dist[::-1] # dim 2 x timesteps x N
uind[0] *= -1 # change sign for ux (to get correct rotation)
# sum over vortices
utot = uind.sum(2) # dim 2 x timesteps


# ## plot time histories and psd for induced velocity



plt.figure(2)
plt.subplot(1,2,1)
plt.plot(t,utot[0],label='u')
plt.plot(t,utot[1],label='v')
plt.legend()
plt.subplot(1,2,2)
(valu,freq) = psd(utot[0],Fs=1/ts,detrend='mean')
(valv,freq) = psd(utot[1],Fs=1/ts,detrend='mean')
plt.loglog(freq[1:],valu[1:],label='u')
plt.loglog(freq[1:],valv[1:],label='v')
plt.loglog(freq[1:],valu[1:]+valv[1:],label='tot')
plt.legend()

#plt.show()

 # Find Liepmann spectrum that fits the Guu (E11) best and plot it.


def E11(f,uu,tt):
    return 4*uu*uu*tt/(1+(2*np.pi*f*tt)**2)
def E22(f,uu,tt):
    return 4*uu*uu*tt*(2*(2*np.pi*f*tt)**2)/(1+(2*np.pi*f*tt)**2)**2
from scipy import optimize
opti,_ = optimize.curve_fit(E11,freq[1:],valu[1:],method='lm') # it is probably better to use the log of Guu and of valu here
uuopt,ttopt = opti
plt.loglog(freq[1:],E11(freq[1:],uuopt,ttopt),label='E11')
plt.loglog(freq[1:],E22(freq[1:],uuopt,ttopt),label='E22')
plt.legend()
#plt.savefig('spectra.pdf')
#plt.show()

opti



'''

