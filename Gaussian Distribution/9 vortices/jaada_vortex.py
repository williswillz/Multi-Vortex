from pysces import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import scipy.fftpack


panels = Vortices()
panels.core_radius = 0.01
q = np.zeros((500,2))
number_of_vortices = 10

xvort1 = np.zeros((500,2))
vel1 = np.zeros((500,2))
a = 85

for i in range (85, 500):
    i = i * 0.1
    q[a,0] = 10
    q[a,1] = 1
    xvort1[a,0] = i
    xvort1[a,1] = 0
    vel1[:,:] = panels.induced_velocity_single(q, xvort1, 0.8147)
    a = a + 1


print xvort1  
print vel1

xvort2 = np.zeros((500,2))
vel2 = np.zeros((500,2))
a = 41

for i in range (41, 500):
    i = i * 0.1
    q[a,0] = 10
    q[a,1] = 1
    xvort2[a,0] = i
    xvort2[a,1] = 0
    vel2[:,:] = panels.induced_velocity_single(q, xvort2, 0.1576)
    a = a + 1


print xvort2  
print vel2

xvort3 = np.zeros((500,2))
vel3 = np.zeros((500,2))
a = 78


for i in range (78, 500):
    i = i * 0.1
    q[a,0] = 10
    q[a,1] = 1
    xvort3[a,0] = i
    xvort3[a,1] = 0
    vel3[:,:] = panels.induced_velocity_single(q, xvort3, 0.6557)
    a = a + 1

print xvort3
print vel3

xvort4 = np.zeros((500,2))
vel4 = np.zeros((500,2))
a = 62

for i in range (62, 500):
    i = i * 0.1
    q[a,0] = 10
    q[a,1] = 1
    xvort4[a,0] = i
    xvort4[a,1] = 0
    vel4[:,:] = panels.induced_velocity_single(q, xvort4, 0.7060)
    a = a + 1


print xvort4
print vel4

xvort5 = np.zeros((500,2))
vel5 = np.zeros((500,2))
a = 23


for i in range (23, 500):
    i = i * 0.1
    q[a,0] = 10
    q[a,1] = 1
    xvort5[a,0] = i
    xvort5[a,1] = 0
    vel5[:,:] = panels.induced_velocity_single(q, xvort5, 0.9058)
    a = a + 1


print xvort5
print vel5

xvort6 = np.zeros((500,2))
vel6 = np.zeros((500,2))
a = 5


for i in range (5, 500):
    i = i * 0.1
    q[a,0] = 10
    q[a,1] = 1
    xvort6[a,0] = i
    xvort6[a,1] = 0
    vel6[:,:] = panels.induced_velocity_single(q, xvort6, 0.9706)
    a = a + 1


print xvort6
print vel6

xvort7 = np.zeros((500,2))
vel7 = np.zeros((500,2))
a = 38


for i in range (38, 500):
    i = i * 0.1
    q[a,0] = 10
    q[a,1] = 1
    xvort7[a,0] = i
    xvort7[a,1] = 0
    vel7[:,:] = panels.induced_velocity_single(q, xvort7, 0.0357)
    a = a + 1


print xvort7
print vel7

xvort8 = np.zeros((500,2))
vel8 = np.zeros((500,2))
a = 35


for i in range (35, 500):
    i = i * 0.1
    q[a,0] = 10
    q[a,1] = 1
    xvort8[a,0] = i
    xvort8[a,1] = 0
    vel8[:,:] = panels.induced_velocity_single(q, xvort8, 0.0318)
    a = a + 1


print xvort8
print vel8

xvort9 = np.zeros((500,2))
vel9 = np.zeros((500,2))
a = 35


for i in range (35, 500):
    i = i * 0.1
    q[a,0] = 10
    q[a,1] = 1
    xvort9[a,0] = i
    xvort9[a,1] = 0
    vel9[:,:] = panels.induced_velocity_single(q, xvort9, 0.1270)
    a = a + 1


print xvort9  
print vel9

xvort10 = np.zeros((500,2))
vel10 = np.zeros((500,2))
a = 90


for i in range (90, 500):
    i = i * 0.1
    q[a,0] = 10
    q[a,1] = 1
    xvort10[a,0] = i
    xvort10[a,1] = 0
    vel10[:,:] = panels.induced_velocity_single(q, xvort10, 1.9572)
    a = a + 1


print xvort10  
print vel10


vel_tot = np.zeros((500,2))
vel_tot = vel1 + vel2 + vel3 + vel4 + vel5 + vel6 + vel7 + vel8 +vel9 + vel10
print vel_tot


t = np.linspace(0,50,500)

vel_tot_mag = np.zeros((500,1))
vel_tot_mag = (vel_tot[:,0]**2 + vel_tot[:,1]**2)**0.5
print vel_tot_mag

savetxt('single-vortex_10_10_sahitya'+str()+'.csv',np.column_stack((t, vel_tot_mag[:])), fmt='%5s', delimiter=',')

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
u3 = plt.plot(t, vel3[:,0], c='c', label='-0.48, u')
v3 = plt.plot(t, vel3[:,1], c='m', label='-0.48, v')
u4 = plt.plot(t, vel4[:,0], c='y', label='-0.48, u')
v4 = plt.plot(t, vel4[:,1], c='b', label='-0.48, v')
u5 = plt.plot(t, vel5[:,0], c='k', label='-0.48, u')
v5 = plt.plot(t, vel5[:,1], c='g', label='-0.48, v')
u6 = plt.plot(t, vel6[:,0], c='c', label='-0.48, u')
v6 = plt.plot(t, vel6[:,1], c='m', label='-0.48, v')
u7 = plt.plot(t, vel7[:,0], c='y', label='-0.48, u')
v7 = plt.plot(t, vel7[:,1], c='b', label='-0.48, v')
u8 = plt.plot(t, vel8[:,0], c='k', label='-0.48, u')
v8 = plt.plot(t, vel8[:,1], c='g', label='-0.48, v')
u9 = plt.plot(t, vel9[:,0], c='m', label='-0.48, u')
v9 = plt.plot(t, vel9[:,1], c='c', label='-0.48, v')
u10 = plt.plot(t, vel10[:,0], c='y', label='-0.48, u')
v10 = plt.plot(t, vel10[:,1], c='k', label='-0.48, v')
plt.legend(loc='best')
savefig('induced_velocity_individual_10_10_sahitya.pdf')
plt.show()

fig, ax = plt.subplots()
vxtot = plt.plot(t, vel_tot[:,0], c='r', label='x_tangential')
vytot = plt.plot(t, vel_tot[:,1], c='g', label='y_tangential')
v3_mag = plt.plot(t, vel_tot_mag[:], c='b', label='induced_velocity')
plt.legend(loc='best')
savefig('induced_velocity_10_10_sahitya.pdf')
plt.show()


# Number of samplepoints
N = 500
# sample spacing
T = 1.0 / 250
x = np.linspace(0.0, N*T, N)
y = vel_tot_mag
yf = scipy.fftpack.fft(y)
xf = np.linspace(0.0, 1.0/(2.0*T), N/2)

fig, ax = plt.subplots()
ax.plot(xf, 2.0/N * np.abs(yf[:N//2]))
ax.set_xscale('log')
savefig('frequency_'+str(N)+str(number_of_vortices)+'.pdf')
plt.show()