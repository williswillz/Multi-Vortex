import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import psd
from pysces import *


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


angle = 0
q = naca_airfoil("0012", 5, zero_thick_te=True)  			#airfoil input
qx = q[:,0]
qy = q[:,1]
q_mid = (q[1:] + q[:-1])/2									#collacation points or mid-point on the panels
print q_mid
print q


#num_pan = 11
##


obsx = q_mid[:,0]
obsy = q_mid[:,1]
#print obsx
#print obsy
obsX = np.vstack((obsx,obsy))
#print obsX

length = np.sqrt((qx[1:] - qx[:-1])**2 + (qy[1:] - qy[:-1])**2)

#for i in range (0, (num_pan/2 * 4)):
#    if q[i+1,0] - q[i,0] <=0:
#        beta = np.arccos((q[0,i+1] - q[0,i]) / length[i])

print length 


beta = q[1:] - q[:-1]
print beta
#print beta[:,1]
slope = np.arctan(beta[:,1] / beta[:,0])
print slope    

if np.any(beta[:,0] <= 0):
    slope = np.arctan(beta[:,1] / beta[:,0])
elif np.any(beta[:,0] > 0):
    slope = np.pi + np.arctan(beta[:,1] / beta[:,0])
print slope

#print qx[1:] - qx[:-1]

#if np.any(qx[1:] - qx[:-1] <= 0.0):
#    beta = np.arccos((qy[1:] - qy[:-1]) / length)
#else:
#    beta = np.pi + np.arccos((qy[1:] - qy[:-1]) / length)
#
#print beta    
    
#beta = ((-qy[1:] + qy[:-1]) / (-qx[1:] + qx[:-1]))
#print beta

'''
   self.xa, self.ya = xa, ya  # panel starting-point
        self.xb, self.yb = xb, yb  # panel ending-point
        
        self.xc, self.yc = (xa + xb) / 2, (ya + yb) / 2  # panel center
        self.length = numpy.sqrt((xb - xa)**2 + (yb - ya)**2)  # panel length
        
        # orientation of panel (angle between x-axis and panel's normal)
        if xb - xa <= 0.0:
            self.beta = numpy.arccos((yb - ya) / self.length)
        elif xb - xa > 0.0:
            self.beta = numpy.pi + numpy.arccos(-(yb - ya) / self.length)
'''
length = 10. #nondimensional length of window
height = 5 #window height
N = 400 #number of vortices
gammas = 1. #vortex strength RMS (normal distribution)
rscale = 0.1 #vortex size scale (rayleigh distribution parameter)
t0 = -1.#start time for observation of convection
t1 = 1.#end time
ts = 0.001 # time step
v0 = 5 #convection speed


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
obsx = -v0*t
obsy = np.zeros_like(obsx)
obsX = np.vstack((obsx,obsy))