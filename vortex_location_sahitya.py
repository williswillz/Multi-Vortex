from pysces import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *


panels = Vortices()
panels.core_radius = 0.01


#x = 15
#y = 1
#q = np.array([x,y]).T

xvort1 = np.zeros((500,2))
vel1 = np.zeros((500,2))
q = np.zeros((500,2))
a = 100


for i in range (100, 500):
    
    i = i * 0.1
    q[a,0] = 15
    q[a,1] = 1
    xvort1[a,0] = i
    xvort1[a,1] = 0
    vel1[:,:] = panels.induced_velocity_single(q, xvort1, 0.14)
    a = a + 1


print xvort1  
print vel1

xvort2 = np.zeros((500,2))
vel2 = np.zeros((500,2))
q = np.zeros((500,2))
a = 0


for i in range (0, 500):
    
    i = i * 0.1
    q[a,0] = 15
    q[a,1] = 1
    xvort2[a,0] = i
    xvort2[a,1] = 0
    vel2[:,:] = panels.induced_velocity_single(q, xvort2, -0.4)
    a = a + 1


print xvort2 
print vel2

# x = 15
# y = 1
# q = np.array([x,y]).T

# xvort2 = np.zeros((500,2))
# vel2 = np.zeros((500,2))
# a = 0


# for i in range (0, 500):
#     i = i * 0.1
#     xvort2[a,0] = i
#     xvort2[a,1] = 0
#     vel2[:,:] = panels.induced_velocity_single(q, xvort2, -0.4)
#     a = a + 1


# print xvort2  
# print vel2

vel3 = np.zeros((500,2))
vel3 = vel1 + vel2
print vel3

savetxt('single-vortex_array_sp'+str()+'.csv',np.column_stack((vel1[:],vel2[:],vel3[:])), fmt='%5s', delimiter=',')

t = np.linspace(0,50,500)

vel3_mag = np.zeros((500,1))
vel3_mag = (vel3[:,0]**2 + vel3[:,1]**2)**0.5
print vel3_mag

# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# ax1.scatter(t, vel1[:,0], c='b', edgecolors='none', label='-0.24, u')
# ax1.scatter(t, vel1[:,1], c='r', edgecolors='none', label='-0.24, v')
# ax1.scatter(t, vel2[:,0], c='k', edgecolors='none', label='-0.48, u')
# ax1.scatter(t, vel2[:,1], c='g', edgecolors='none', label='-0.48, v')
# plt.xlabel('time')
# plt.legend(loc='upper left');
# plt.grid(True)
# plt.show()

fig, ax = plt.subplots()
u1 = plt.plot(t, vel1[:,0], c='b', label='-0.24, u')
v1 = plt.plot(t, vel1[:,1], c='r', label='-0.24, v')
u2 = plt.plot(t, vel2[:,0], c='k', label='-0.48, u')
v2 = plt.plot(t, vel2[:,1], c='g', label='-0.48, v')
plt.legend(loc='upper left')
savefig('induced_velocity__location_array_sp.pdf')
plt.show()

fig, ax = plt.subplots()
vx3 = plt.plot(t, vel3[:,0], c='r', label='x_tangential')
vy3 = plt.plot(t, vel3[:,1], c='g', label='y_tangential')
v3_mag = plt.plot(t, vel3_mag[:], c='b', label='induced_velocity')
plt.legend(loc='best')
savefig('induced_velocity_location_tot_array_sp.pdf')
plt.show()


# print vel  

# plt.plot(x,vel[:,0], '-x', label="u")
# plt.plot(x,vel[:,1], '-+', label="v")
# plt.legend()
# plt.show()