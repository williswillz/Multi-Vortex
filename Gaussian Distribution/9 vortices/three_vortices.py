from pysces import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *

# The test is conducted while keeping NACA0012 in mind
# observation point = (25,0)
# strength of vortex_1 = 0.5

panels = Vortices()
panels.core_radius = 0.01

# vortex 1
# starts from (0,0) to (50,0)

xvort1 = np.zeros((500,2))
vel1 = np.zeros((500,2))
q = np.zeros((500,2))
a = 0

for i in range (0, 500):
    i = i * 0.1
    q[a,0] = 25 #observation point
    q[a,1] = 0
    xvort1[a,0] = i
    xvort1[a,1] = 0.2
    vel1[:,:] = panels.induced_velocity_single(q, xvort1, 0.5)
    a = a + 1

print xvort1  
print vel1

# vortex 2
# starts from (0.1,0) to (50,0)

xvort2 = np.zeros((500,2))
vel2 = np.zeros((500,2))
a = 0


for i in range (0, 500):
    i = i * 0.1
    q[a,0] = 25
    q[a,1] = 0
    xvort2[a,0] = i
    xvort2[a,1] = 0
    vel2[:,:] = panels.induced_velocity_single(q, xvort2, 0.5)
    a = a + 1

print xvort2
print vel2

xvort3 = np.zeros((500,2))
vel3 = np.zeros((500,2))
a = 0


for i in range (0, 500):
    i = i * 0.1
    q[a,0] = 25
    q[a,1] = 0
    xvort3[a,0] = i
    xvort3[a,1] = -0.1
    vel3[:,:] = panels.induced_velocity_single(q, xvort3, 0.5)
    a = a + 1


print xvort3  
print vel3

vel4 = np.zeros((500,2))
vel4 = vel1 + vel2 + vel3
print vel4

t = np.linspace(0,50,500)

vel4_mag = np.zeros((500,1))
vel4_mag = (vel4[:,0]**2 + vel4[:,1]**2)**0.5
print vel4_mag

savetxt('region_1'+str()+'.csv',np.column_stack((vel1[:],vel2[:],vel3[:],vel4[:])), fmt='%5s', delimiter=',')

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
u1 = plt.plot(t, vel1[:,0], c='b', label='0.5, u')
v1 = plt.plot(t, vel1[:,1], c='r', label='0.5, v')
u2 = plt.plot(t, vel2[:,0], c='k', label='0.5, u')
v2 = plt.plot(t, vel2[:,1], c='g', label='0.5, v')
u3 = plt.plot(t, vel3[:,0], c='y', label='0.5, u')
v3 = plt.plot(t, vel3[:,1], c='m', label='0.5, v')
plt.legend(loc='best')
savefig('induced_velocity_individual_three_vortices.pdf')
plt.show()

fig, ax = plt.subplots()
vx4 = plt.plot(t, vel4[:,0], c='r', label='x_tangential')
vy4 = plt.plot(t, vel4[:,1], c='g', label='y_tangential')
v4_mag = plt.plot(t, vel4_mag[:], c='b', label='induced_velocity')
plt.legend(loc='best')
savefig('induced_velocity_three_vortices.pdf')
plt.show()


# print vel  

# plt.plot(x,vel[:,0], '-x', label="u")
# plt.plot(x,vel[:,1], '-+', label="v")
# plt.legend()
# plt.show()