from pysces import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *

# The test is conducted while keeping NACA0012 in mind
# observation point = (25,0)
# strength of vortex_1 = 0.5

panels = Vortices()
panels.core_radius = 0.01

number_of_vortices = 10

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
    xvort1[a,1] = 0
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
    xvort2[a,1] = 0.1
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
    xvort3[a,1] = 0.2
    vel3[:,:] = panels.induced_velocity_single(q, xvort3, 0.5)
    a = a + 1

print xvort3  
print vel3

xvort4 = np.zeros((500,2))
vel4 = np.zeros((500,2))
a = 0

for i in range (0, 500):
    i = i * 0.1
    q[a,0] = 25
    q[a,1] = 0
    xvort4[a,0] = i
    xvort4[a,1] = 0.3
    vel4[:,:] = panels.induced_velocity_single(q, xvort4, 0.5)
    a = a + 1

print xvort4
print vel4

xvort5 = np.zeros((500,2))
vel5 = np.zeros((500,2))
a = 0

for i in range (0, 500):
    i = i * 0.1
    q[a,0] = 25
    q[a,1] = 0
    xvort5[a,0] = i
    xvort5[a,1] = 0.4
    vel5[:,:] = panels.induced_velocity_single(q, xvort5, 0.5)
    a = a + 1

print xvort5
print vel5

xvort6 = np.zeros((500,2))
vel6 = np.zeros((500,2))
a = 0


for i in range (0, 500):
    i = i * 0.1
    q[a,0] = 25
    q[a,1] = 0
    xvort6[a,0] = i
    xvort6[a,1] = -0.1
    vel6[:,:] = panels.induced_velocity_single(q, xvort6, 0.5)
    a = a + 1

print xvort6
print vel6

xvort7 = np.zeros((500,2))
vel7 = np.zeros((500,2))
a = 0

for i in range (0, 500):
    i = i * 0.1
    q[a,0] = 25
    q[a,1] = 0
    xvort7[a,0] = i
    xvort7[a,1] = -0.2
    vel7[:,:] = panels.induced_velocity_single(q, xvort7, 0.5)
    a = a + 1

print xvort7 
print vel7

xvort8 = np.zeros((500,2))
vel8 = np.zeros((500,2))
a = 0

for i in range (0, 500):
    i = i * 0.1
    q[a,0] = 25
    q[a,1] = 0
    xvort8[a,0] = i
    xvort8[a,1] = -0.3
    vel8[:,:] = panels.induced_velocity_single(q, xvort8, 0.5)
    a = a + 1

print xvort8
print vel8

xvort9 = np.zeros((500,2))
vel9 = np.zeros((500,2))
a = 0

for i in range (0, 500):
    i = i * 0.1
    q[a,0] = 25
    q[a,1] = 0
    xvort9[a,0] = i
    xvort9[a,1] = -0.4
    vel9[:,:] = panels.induced_velocity_single(q, xvort9, 0.5)
    a = a + 1

print xvort9
print vel9

vel10 = np.zeros((500,2))
vel10 = vel1 + vel2 + vel3 + vel4 + vel5 + vel6 + vel7 + vel8 + vel9
print vel10

t = np.linspace(0,50,500)
vel10_mag = np.zeros((500,1))
vel10_mag = (vel10[:,0]**2 + vel10[:,1]**2)**0.5
print vel10_mag

savetxt('vertical_var'+str()+'.csv',np.column_stack((t,vel10[:])), fmt='%5s', delimiter=',')

# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# ax1.scatter(t, vel1[:,0], c='b', edgecolors='none', label='-0.24, u')
# ax1.scatter(t, vel1[:,1], c='r', edgecolors='none', label='-0.24, v')
# ax1.scatter(t, vel2[:,0], c='k', edgecolors='none', label='0.50, u')
# ax1.scatter(t, vel2[:,1], c='g', edgecolors='none', label='0.50, v')
# plt.xlabel('time')
# plt.legend(loc='upper left');
# plt.grid(True)
# plt.show()

fig, ax = plt.subplots()
u1 = plt.plot(t, vel1[:,0], c='b', label='0.5, u')
v1 = plt.plot(t, vel1[:,1], c='r', label='0.5, v')
u2 = plt.plot(t, vel2[:,0], c='k', label='0.50, u')
v2 = plt.plot(t, vel2[:,1], c='g', label='0.50, v')
u3 = plt.plot(t, vel3[:,0], c='c', label='0.50, u')
v3 = plt.plot(t, vel3[:,1], c='m', label='0.50, v')
u4 = plt.plot(t, vel4[:,0], c='y', label='0.50, u')
v4 = plt.plot(t, vel4[:,1], c='b', label='0.50, v')
u5 = plt.plot(t, vel5[:,0], c='k', label='0.50, u')
v5 = plt.plot(t, vel5[:,1], c='g', label='0.50, v')
u6 = plt.plot(t, vel6[:,0], c='c', label='0.50, u')
v6 = plt.plot(t, vel6[:,1], c='m', label='0.50, v')
u7 = plt.plot(t, vel7[:,0], c='y', label='0.50, u')
v7 = plt.plot(t, vel7[:,1], c='b', label='0.50, v')
u8 = plt.plot(t, vel8[:,0], c='k', label='0.50, u')
v8 = plt.plot(t, vel8[:,1], c='g', label='0.50, v')
u9 = plt.plot(t, vel9[:,0], c='m', label='0.50, u')
v9 = plt.plot(t, vel9[:,1], c='c', label='0.50, v')
plt.legend(loc='upper right')
# ax.axhline(y=0, color='k')
# ax.axvline(x=25, color='k')
savefig(str(number_of_vortices)+'_vortices_uv_velocity.pdf')
plt.show()

fig, ax = plt.subplots()
u1 = plt.plot(t, vel1[:,0], c='b', label='0.5, u')
#v1 = plt.plot(t, vel1[:,1], c='r', label='0.5, v')
u2 = plt.plot(t, vel2[:,0], c='k', label='0.50, u')
#v2 = plt.plot(t, vel2[:,1], c='g', label='0.50, v')
u3 = plt.plot(t, vel3[:,0], c='c', label='0.50, u')
#v3 = plt.plot(t, vel3[:,1], c='m', label='0.50, v')
u4 = plt.plot(t, vel4[:,0], c='y', label='0.50, u')
#v4 = plt.plot(t, vel4[:,1], c='b', label='0.50, v')
u5 = plt.plot(t, vel5[:,0], c='k', label='0.50, u')
#v5 = plt.plot(t, vel5[:,1], c='g', label='0.50, v')
u6 = plt.plot(t, vel6[:,0], c='c', label='0.50, u')
#v6 = plt.plot(t, vel6[:,1], c='m', label='0.50, v')
u7 = plt.plot(t, vel7[:,0], c='y', label='0.50, u')
#v7 = plt.plot(t, vel7[:,1], c='b', label='0.50, v')
u8 = plt.plot(t, vel8[:,0], c='k', label='0.50, u')
#v8 = plt.plot(t, vel8[:,1], c='g', label='0.50, v')
u9 = plt.plot(t, vel9[:,0], c='m', label='0.50, u')
#v9 = plt.plot(t, vel9[:,1], c='c', label='0.50, v')
plt.legend(loc='upper right')
savefig(str(number_of_vortices)+'_vortices_u_velocity.pdf')
plt.show()

fig, ax = plt.subplots()
#u1 = plt.plot(t, vel1[:,0], c='b', label='0.5, u')
v1 = plt.plot(t, vel1[:,1], c='r', label='0.5, v')
#u2 = plt.plot(t, vel2[:,0], c='k', label='0.50, u')
v2 = plt.plot(t, vel2[:,1], c='g', label='0.50, v')
#u3 = plt.plot(t, vel3[:,0], c='c', label='0.50, u')
v3 = plt.plot(t, vel3[:,1], c='m', label='0.50, v')
#u4 = plt.plot(t, vel4[:,0], c='y', label='0.50, u')
v4 = plt.plot(t, vel4[:,1], c='b', label='0.50, v')
#u5 = plt.plot(t, vel5[:,0], c='k', label='0.50, u')
v5 = plt.plot(t, vel5[:,1], c='g', label='0.50, v')
#u6 = plt.plot(t, vel6[:,0], c='c', label='0.50, u')
v6 = plt.plot(t, vel6[:,1], c='m', label='0.50, v')
#u7 = plt.plot(t, vel7[:,0], c='y', label='0.50, u')
v7 = plt.plot(t, vel7[:,1], c='b', label='0.50, v')
#u8 = plt.plot(t, vel8[:,0], c='k', label='0.50, u')
v8 = plt.plot(t, vel8[:,1], c='g', label='0.50, v')
#u9 = plt.plot(t, vel9[:,0], c='m', label='0.50, u')
v9 = plt.plot(t, vel9[:,1], c='c', label='0.50, v')
plt.legend(loc='upper right')
savefig(str(number_of_vortices)+'_vortices_v_velocity.pdf')
plt.show()

fig, ax = plt.subplots()
vx10 = plt.plot(t, vel10[:,0], c='r', label='x_tangential')
vy10 = plt.plot(t, vel10[:,1], c='g', label='y_tangential')
v10_mag = plt.plot(t, vel10_mag[:], c='b', label='induced_velocity')
plt.legend(loc='best')
savefig(str(number_of_vortices)+'_vortices_total_velocity.pdf')
plt.show()


# print vel  

# plt.plot(x,vel[:,0], '-x', label="u")
# plt.plot(x,vel[:,1], '-+', label="v")
# plt.legend()
# plt.show()