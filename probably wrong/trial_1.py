from pysces import *
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from pylab import *
import time

start = time.time()

panels = Vortices()
panels.core_radius = 0.01
num_steps = 500
Uinfty = (1.421,0)


x = 0
y = 0
q = np.array([x,y]).T
xvort = np.array([-5,0])
vel = panels.induced_velocity_single(q, xvort, 1)
plt.plot(x,vel[:,0], '-x', label="u")
plt.plot(x,vel[:,1], '-+', label="v")
plt.legend()
plt.show()


#free vortices

#v1	=	(	1.724339281	,	3.1	)
v1 = (-5.0,-0.1)
v2 = (-5.0,-0.2)
#s1	=	-0.00810137
s1 = 0.24
s2 = -0.24

vort = Vortices([v1,v2],[s1,s2])

flow = ExplicitEuler(dt, Uinfty, panels, wake=vort, need_force='wake_impulse')

for i in range(1,num_steps):
    flow.advance()

vel = panels.induced_velocity_single(q, vort, 0.24)
plt.plot(x,vel[:,0], '-x', label="u")
plt.plot(x,vel[:,1], '-+', label="v")
plt.legend()
plt.show()
