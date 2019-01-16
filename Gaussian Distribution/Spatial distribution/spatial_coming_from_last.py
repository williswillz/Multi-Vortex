# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 10:46:52 2019

@author: WS1
"""

from pysces import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import scipy.fftpack
import random
import time

start = time.time()


panels = Vortices()

# The test is conducted while keeping NACA0012 in mind
# observation point = (25,0)
# strength of vortex_1 = 0.5

panels = Vortices()
panels.core_radius = 0.01

number_of_vortices = 22

# Region-1 = 16 vortices
# vortex 1
# starts from (0,0) to (50,0)

xvort1 = np.zeros((5000,2))
vel1 = np.zeros((5000,2))
q = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250 #observation point
    q[a,1] = 0
    xvort1[a,0] = i
    xvort1[a,1] = random.uniform(-1,0)
    vel1[:,:] = panels.induced_velocity_single(q, xvort1, random.uniform(-1,1))
    a = a + 1

print xvort1  
print vel1

# vortex 2
# starts from (0.1,0) to (50,0)

xvort2 = np.zeros((5000,2))
vel2 = np.zeros((5000,2))
a = 0


for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort2[a,0] = i
    xvort2[a,1] = random.uniform(-1,0)
    vel2[:,:] = panels.induced_velocity_single(q, xvort2, random.uniform(-1,1))
    a = a + 1

print xvort2
print vel2

xvort3 = np.zeros((5000,2))
vel3 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort3[a,0] = i
    xvort3[a,1] = random.uniform(-1,0)
    vel3[:,:] = panels.induced_velocity_single(q, xvort3, random.uniform(-1,1))
    a = a + 1

print xvort3  
print vel3

xvort4 = np.zeros((5000,2))
vel4 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort4[a,0] = i
    xvort4[a,1] = random.uniform(-1,0)
    vel4[:,:] = panels.induced_velocity_single(q, xvort4, random.uniform(-1,1))
    a = a + 1

print xvort4
print vel4

xvort5 = np.zeros((5000,2))
vel5 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort5[a,0] = i
    xvort5[a,1] = random.uniform(-1,0)
    vel5[:,:] = panels.induced_velocity_single(q, xvort5, random.uniform(-1,1))
    a = a + 1

print xvort5
print vel5

xvort6 = np.zeros((5000,2))
vel6 = np.zeros((5000,2))
a = 0


for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort6[a,0] = i
    xvort6[a,1] = random.uniform(-1,0)
    vel6[:,:] = panels.induced_velocity_single(q, xvort6, random.uniform(-1,1))
    a = a + 1

print xvort6
print vel6

xvort7 = np.zeros((5000,2))
vel7 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort7[a,0] = i
    xvort7[a,1] = random.uniform(-1,0)
    vel7[:,:] = panels.induced_velocity_single(q, xvort7, random.uniform(-1,1))
    a = a + 1

print xvort7 
print vel7

xvort8 = np.zeros((5000,2))
vel8 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort8[a,0] = i
    xvort8[a,1] = random.uniform(-1,0)
    vel8[:,:] = panels.induced_velocity_single(q, xvort8, random.uniform(-1,1))
    a = a + 1

print xvort8
print vel8

xvort9 = np.zeros((5000,2))
vel9 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort9[a,0] = i
    xvort9[a,1] = random.uniform(0,1)
    vel9[:,:] = panels.induced_velocity_single(q, xvort9, random.uniform(-1,1))
    a = a + 1

print xvort9
print vel9

xvort10 = np.zeros((5000,2))
vel10 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort10[a,0] = i
    xvort10[a,1] = random.uniform(0,1)
    vel10[:,:] = panels.induced_velocity_single(q, xvort10, random.uniform(-1,1))
    a = a + 1

print xvort10
print vel10

xvort11 = np.zeros((5000,2))
vel11 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort11[a,0] = i
    xvort11[a,1] = random.uniform(0,1)
    vel11[:,:] = panels.induced_velocity_single(q, xvort11, random.uniform(-1,1))
    a = a + 1

print xvort11
print vel11

xvort12 = np.zeros((5000,2))
vel12 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort12[a,0] = i
    xvort12[a,1] = random.uniform(0,1)
    vel12[:,:] = panels.induced_velocity_single(q, xvort12, random.uniform(-1,1))
    a = a + 1

print xvort12
print vel12

xvort13 = np.zeros((5000,2))
vel13 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort13[a,0] = i
    xvort13[a,1] = random.uniform(0,1)
    vel13[:,:] = panels.induced_velocity_single(q, xvort13, random.uniform(-1,1))
    a = a + 1

print xvort13
print vel13

xvort14 = np.zeros((5000,2))
vel14 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort14[a,0] = i
    xvort14[a,1] = random.uniform(0,1)
    vel14[:,:] = panels.induced_velocity_single(q, xvort14, random.uniform(-1,1))
    a = a + 1

print xvort14
print vel14

xvort15 = np.zeros((5000,2))
vel15 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort15[a,0] = i
    xvort15[a,1] = random.uniform(0,1)
    vel15[:,:] = panels.induced_velocity_single(q, xvort15, random.uniform(-1,1))
    a = a + 1

print xvort15
print vel15

xvort16 = np.zeros((5000,2))
vel16 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort16[a,0] = i
    xvort16[a,1] = random.uniform(0,1)
    vel16[:,:] = panels.induced_velocity_single(q, xvort16, random.uniform(-1,1))
    a = a + 1

print xvort16
print vel16

vel_tot_reg1 = np.zeros((5000,2))
vel_tot_reg1 = vel1 + vel2 + vel3 + vel4 + vel5 + vel6 + vel7 + vel8 + vel9 + vel10 + vel11 + vel12 + vel13 + vel14 + vel15 + vel16
print vel_tot_reg1

# Region-2 = 6 vortices

xvort17 = np.zeros((5000,2))
vel17 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort17[a,0] = i
    xvort17[a,1] = random.uniform(1,2)
    vel17[:,:] = panels.induced_velocity_single(q, xvort17, random.uniform(-1,1))
    a = a + 1

print xvort17
print vel17

xvort18 = np.zeros((5000,2))
vel18 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort18[a,0] = i
    xvort18[a,1] = random.uniform(1,2)
    vel18[:,:] = panels.induced_velocity_single(q, xvort18, random.uniform(-1,1))
    a = a + 1

print xvort18
print vel18

xvort19 = np.zeros((5000,2))
vel19 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort19[a,0] = i
    xvort19[a,1] = random.uniform(1,2)
    vel19[:,:] = panels.induced_velocity_single(q, xvort19, random.uniform(-1,1))
    a = a + 1

print xvort19
print vel19

xvort20 = np.zeros((5000,2))
vel20 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort20[a,0] = i
    xvort20[a,1] = random.uniform(-1,-2)
    vel20[:,:] = panels.induced_velocity_single(q, xvort20, random.uniform(-1,1))
    a = a + 1

print xvort20
print vel20

xvort21 = np.zeros((5000,2))
vel21 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort21[a,0] = i
    xvort21[a,1] = random.uniform(-1,-2)
    vel21[:,:] = panels.induced_velocity_single(q, xvort21, random.uniform(-1,1))
    a = a + 1

print xvort21
print vel21

xvort22 = np.zeros((5000,2))
vel22 = np.zeros((5000,2))
a = 0

for i in range (0, 5000):
    i = i * 0.1
    q[a,0] = 250
    q[a,1] = 0
    xvort22[a,0] = i
    xvort22[a,1] = random.uniform(-1,-2)
    vel22[:,:] = panels.induced_velocity_single(q, xvort22, random.uniform(-1,1))
    a = a + 1

print xvort22
print vel22

vel_tot_reg2 = np.zeros((5000,2))
vel_tot_reg2 = vel17 + vel18 + vel19 + vel20 + vel21 + vel22
print vel_tot_reg2

vel_tot = vel_tot_reg1 + vel_tot_reg2

t = np.linspace(0,50,5000)
vel_tot_mag_reg1 = np.zeros((5000,1))
vel_tot_mag_reg1 = (vel_tot_reg1[:,0]**2 + vel_tot_reg1[:,1]**2)**0.5
vel_tot_mag_reg2 = np.zeros((5000,1))
vel_tot_mag_reg2 = (vel_tot_reg2[:,0]**2 + vel_tot_reg2[:,1]**2)**0.5
vel_tot_mag = vel_tot_mag_reg1 + vel_tot_mag_reg2
print vel_tot_mag_reg1
print vel_tot_mag_reg2

end = time.time()
print end - start

#savetxt('25vertical_var_rand_5000'+str()+'.csv',np.column_stack((t,vel_tot[:])), fmt='%5s', delimiter=',')

# fig = plt.figure()
# ax1 = 6666fig.add_subplot(111)
# ax1.scatter(t, vel1[:,0], c='b', edgecolors='none', label='-0.24, u')
# ax1.scatter(t, vel1[:,1], c='r', edgecolors='none', label='-0.24, v')
# ax1.scatter(t, vel2[:,0], c='k', edgecolors='none', label='vortex ')
# ax1.scatter(t, vel2[:,1], c='g', edgecolors='none', label='vortex ')
# plt.xlabel('time')
# plt.legend(loc='upper left');
# plt.grid(True)
# plt.show()

# fig, ax = plt.subplots()
# u1 = plt.plot(t, vel1[:,0], c='b', label='0.5, u')
# v1 = plt.plot(t, vel1[:,1], c='r', label='vortex ')
# u2 = plt.plot(t, vel2[:,0], c='k', label='vortex ')
# v2 = plt.plot(t, vel2[:,1], c='g', label='vortex ')
# u3 = plt.plot(t, vel3[:,0], c='c', label='vortex ')
# v3 = plt.plot(t, vel3[:,1], c='m', label='vortex ')
# u4 = plt.plot(t, vel4[:,0], c='y', label='vortex ')
# v4 = plt.plot(t, vel4[:,1], c='b', label='vortex ')
# u5 = plt.plot(t, vel5[:,0], c='k', label='vortex ')
# v5 = plt.plot(t, vel5[:,1], c='g', label='vortex ')
# u6 = plt.plot(t, vel6[:,0], c='c', label='vortex ')
# v6 = plt.plot(t, vel6[:,1], c='m', label='vortex ')
# u7 = plt.plot(t, vel7[:,0], c='y', label='vortex ')
# v7 = plt.plot(t, vel7[:,1], c='b', label='vortex ')
# u8 = plt.plot(t, vel8[:,0], c='k', label='vortex ')
# v8 = plt.plot(t, vel8[:,1], c='g', label='vortex ')
# u9 = plt.plot(t, vel9[:,0], c='m', label='vortex ')
# v9 = plt.plot(t, vel9[:,1], c='c', label='vortex ')
# plt.legend(loc='upper right')
# # ax.axhline(y=0, color='k')
# # ax.axvline(x=25, color='k')
# savefig(str(number_of_vortices)+'_vortices_uv_velocity.pdf')
# plt.show()

# fig, ax = plt.subplots()
# u1 = plt.plot(t, vel1[:,0], c='gray', label='vortex 1')
# #v1 = plt.plot(t, vel1[:,1], c='r', label='vortex ')
# u2 = plt.plot(t, vel2[:,0], c='firebrick', label='vortex 2')
# #v2 = plt.plot(t, vel2[:,1], c='g', label='vortex ')
# u3 = plt.plot(t, vel3[:,0], c='sandybrown', label='vortex 3')
# #v3 = plt.plot(t, vel3[:,1], c='m', label='vortex ')
# u4 = plt.plot(t, vel4[:,0], c='gold', label='vortex 4')
# #v4 = plt.plot(t, vel4[:,1], c='b', label='vortex ')
# u5 = plt.plot(t, vel5[:,0], c='chartreuse', label='vortex 5')
# #v5 = plt.plot(t, vel5[:,1], c='g', label='vortex ')
# u6 = plt.plot(t, vel6[:,0], c='mediumspringgreen', label='vortex 6')
# #v6 = plt.plot(t, vel6[:,1], c='m', label='vortex ')
# u7 = plt.plot(t, vel7[:,0], c='seagreen', label='vortex 7')
# #v7 = plt.plot(t, vel7[:,1], c='b', label='vortex ')
# u8 = plt.plot(t, vel8[:,0], c='deepskyblue', label='vortex 8')
# #v8 = plt.plot(t, vel8[:,1], c='g', label='vortex ')
# u9 = plt.plot(t, vel9[:,0], c='royalblue', label='vortex 9')
# #v9 = plt.plot(t, vel9[:,1], c='c', label='vortex ')
# u10 = plt.plot(t, vel10[:,0], c='mediumpurple', label='vortex 10')
# u11 = plt.plot(t, vel11[:,0], c='lightcoral', label='vortex 11')
# u12 = plt.plot(t, vel12[:,0], c='violet', label='vortex 12')
# u13 = plt.plot(t, vel13[:,0], c='crimson', label='vortex 13')
# u14 = plt.plot(t, vel14[:,0], c='orange', label='vortex 14')
# u15 = plt.plot(t, vel15[:,0], c='yellow', label='vortex 15')
# u16 = plt.plot(t, vel16[:,0], c='r', label='vortex 16')
# plt.xlabel('time')
# plt.ylabel('induced_velocity u')
# plt.legend(loc='upper right')
# savefig(str(number_of_vortices)+'_vortices_u_velocity_rand_5000.pdf')
# plt.show()

fig, ax = plt.subplots()
u1 = plt.plot(t, vel1[:,0], c='g', label='vortex 1')
u2 = plt.plot(t, vel2[:,0], c='g', label='vortex 2')
u3 = plt.plot(t, vel3[:,0], c='g', label='vortex 3')
u4 = plt.plot(t, vel4[:,0], c='g', label='vortex 4')
u5 = plt.plot(t, vel5[:,0], c='g', label='vortex 5')
u6 = plt.plot(t, vel6[:,0], c='g', label='vortex 6')
u7 = plt.plot(t, vel7[:,0], c='g', label='vortex 7')
u8 = plt.plot(t, vel8[:,0], c='g', label='vortex 8')
u9 = plt.plot(t, vel9[:,0], c='g', label='vortex 9')
u10 = plt.plot(t, vel10[:,0], c='g', label='vortex 10')
u11 = plt.plot(t, vel11[:,0], c='g', label='vortex 11')
u12 = plt.plot(t, vel12[:,0], c='g', label='vortex 12')
u13 = plt.plot(t, vel13[:,0], c='g', label='vortex 13')
u14 = plt.plot(t, vel14[:,0], c='g', label='vortex 14')
u15 = plt.plot(t, vel15[:,0], c='g', label='vortex 15')
u16 = plt.plot(t, vel16[:,0], c='r', label='vortex 16')
u17 = plt.plot(t, vel17[:,0], c='r', label='vortex 17')
u18 = plt.plot(t, vel18[:,0], c='r', label='vortex 18')
u19 = plt.plot(t, vel19[:,0], c='r', label='vortex 19')
u20 = plt.plot(t, vel20[:,0], c='r', label='vortex 20')
u21 = plt.plot(t, vel21[:,0], c='r', label='vortex 21')
u22 = plt.plot(t, vel22[:,0], c='r', label='vortex 22')
plt.xlabel('time')
plt.ylabel('induced_velocity u')
plt.legend(loc='upper right')
#savefig(str(number_of_vortices)+'_vortices_u_velocity_rand_5000.pdf')
plt.show()

# fig, ax = plt.subplots()
# v1 = plt.plot(t, vel1[:,1], c='gray', label='vortex 1')
# #v1 = plt.plot(t, vel1[:,1], c='r', label='vortex ')
# v2 = plt.plot(t, vel2[:,1], c='firebrick', label='vortex 2')
# #v2 = plt.plot(t, vel2[:,1], c='g', label='vortex ')
# v3 = plt.plot(t, vel3[:,1], c='sandybrown', label='vortex 3')
# #v3 = plt.plot(t, vel3[:,1], c='m', label='vortex ')
# v4 = plt.plot(t, vel4[:,1], c='gold', label='vortex 4')
# #v4 = plt.plot(t, vel4[:,1], c='b', label='vortex ')
# v5 = plt.plot(t, vel5[:,1], c='chartreuse', label='vortex 5')
# #v5 = plt.plot(t, vel5[:,1], c='g', label='vortex ')
# v6 = plt.plot(t, vel6[:,1], c='mediumspringgreen', label='vortex 6')
# #v6 = plt.plot(t, vel6[:,1], c='m', label='vortex ')
# v7 = plt.plot(t, vel7[:,1], c='seagreen', label='vortex 7')
# #v7 = plt.plot(t, vel7[:,1], c='b', label='vortex ')
# v8 = plt.plot(t, vel8[:,1], c='deepskyblue', label='vortex 8')
# #v8 = plt.plot(t, vel8[:,1], c='g', label='vortex ')
# v9 = plt.plot(t, vel9[:,1], c='royalblue', label='vortex 9')
# #v9 = plt.plot(t, vel9[:,1], c='c', label='vortex ')
# v10 = plt.plot(t, vel10[:,1], c='mediumpurple', label='vortex 10')
# v11 = plt.plot(t, vel11[:,1], c='lightcoral', label='vortex 11')
# v12 = plt.plot(t, vel12[:,1], c='violet', label='vortex 12')
# v13 = plt.plot(t, vel13[:,1], c='crimson', label='vortex 13')
# v14 = plt.plot(t, vel14[:,1], c='orange', label='vortex 14')
# v15 = plt.plot(t, vel15[:,1], c='yellow', label='vortex 15')
# v16 = plt.plot(t, vel16[:,1], c='r', label='vortex 16')
# plt.xlabel('time')
# plt.ylabel('induced_velocity v')
# plt.legend(loc='upper right')
# savefig(str(number_of_vortices)+'_vortices_v_velocity_rand_5000.pdf')
# plt.show()

fig, ax = plt.subplots()
v1 = plt.plot(t, vel1[:,1], c='g', label='vortex 1')
v2 = plt.plot(t, vel2[:,1], c='g', label='vortex 2')
v3 = plt.plot(t, vel3[:,1], c='g', label='vortex 3')
v4 = plt.plot(t, vel4[:,1], c='g', label='vortex 4')
v5 = plt.plot(t, vel5[:,1], c='g', label='vortex 5')
v6 = plt.plot(t, vel6[:,1], c='g', label='vortex 6')
v7 = plt.plot(t, vel7[:,1], c='g', label='vortex 7')
v8 = plt.plot(t, vel8[:,1], c='g', label='vortex 8')
v9 = plt.plot(t, vel9[:,1], c='g', label='vortex 9')
v10 = plt.plot(t, vel10[:,1], c='g', label='vortex 10')
v11 = plt.plot(t, vel11[:,1], c='g', label='vortex 11')
v12 = plt.plot(t, vel12[:,1], c='g', label='vortex 12')
v13 = plt.plot(t, vel13[:,1], c='g', label='vortex 13')
v14 = plt.plot(t, vel14[:,1], c='g', label='vortex 14')
v15 = plt.plot(t, vel15[:,1], c='g', label='vortex 15')
v16 = plt.plot(t, vel16[:,1], c='r', label='vortex 16')
v17 = plt.plot(t, vel17[:,1], c='r', label='vortex 17')
v18 = plt.plot(t, vel18[:,1], c='r', label='vortex 18')
v19 = plt.plot(t, vel19[:,1], c='r', label='vortex 19')
v20 = plt.plot(t, vel20[:,1], c='r', label='vortex 20')
v21 = plt.plot(t, vel21[:,1], c='r', label='vortex 21')
v22 = plt.plot(t, vel22[:,1], c='r', label='vortex 22')
plt.xlabel('time')
plt.ylabel('induced_velocity v')
plt.legend(loc='upper right')
#savefig(str(number_of_vortices)+'_vortices_v_velocity_rand_5000.pdf')
plt.show()

fig, ax = plt.subplots()
vx10 = plt.plot(t, vel_tot_reg1[:,0], c='r', label='x_tangential')
vy10 = plt.plot(t, vel_tot_reg1[:,1], c='g', label='y_tangential')
v10_mag = plt.plot(t, vel_tot_mag_reg1[:], c='b', label='induced_velocity')
plt.xlabel('time')
plt.ylabel('induced_velocity V_ind_tot')
plt.legend(loc='best')
#savefig(str(number_of_vortices)+'_vortices_total_velocity_rand_5000_region1.pdf')
plt.show()

fig, ax = plt.subplots()
vx10 = plt.plot(t, vel_tot_reg2[:,0], c='r', label='x_tangential')
vy10 = plt.plot(t, vel_tot_reg2[:,1], c='g', label='y_tangential')
v10_mag = plt.plot(t, vel_tot_mag_reg2[:], c='b', label='induced_velocity')
plt.xlabel('time')
plt.ylabel('induced_velocity V_ind_tot')
plt.legend(loc='best')
#savefig(str(number_of_vortices)+'_vortices_total_velocity_rand_5000_region2.pdf')
plt.show()

fig, ax = plt.subplots()
vx10 = plt.plot(t, vel_tot[:,0], c='r', label='x_tangential')
vy10 = plt.plot(t, vel_tot[:,1], c='g', label='y_tangential')
v10_mag = plt.plot(t, vel_tot_mag[:], c='b', label='induced_velocity')
plt.xlabel('time')
plt.ylabel('induced_velocity V_ind_tot')
plt.legend(loc='best')
#savefig(str(number_of_vortices)+'_vortices_total_velocity_rand_5000_region1_2.pdf')
plt.show()

# Number of samplepoints
N = 5000
# sample spacing
T = 1.0 / 2500
x = np.linspace(0.0, N*T, N)
y = vel_tot_mag
yf = scipy.fftpack.fft(y)
xf = np.linspace(0.0, 1.0/(2.0*T), N/2)

fig, ax = plt.subplots()
ax.plot(xf, 2.0/N * np.abs(yf[:N//2]))
ax.set_xscale('log')
plt.xlabel('frequency (Hz)')
plt.ylabel('induced_velocity')
#savefig('frequency_(5,1)_rand_5000'+str(N)+str(number_of_vortices)+'.pdf')
plt.show()

# Number of samplepoints
N = 5000
# sample spacing
T = 1.0 / 2500
x = np.linspace(0.0, N*T, N)
y1 = vel_tot_mag
y2 = vel_tot[:,0]
y3 = vel_tot[:,1]
yf1 = scipy.fftpack.fft(y1)
yf2 = scipy.fftpack.fft(y2)
yf3 = scipy.fftpack.fft(y3)
xf = np.linspace(0.0, 1.0/(2.0*T), N/2)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(xf, 2.0/N * np.abs(yf1[:N//2]), c='b', label='V_ind_tot')
ax1.plot(xf, 2.0/N * np.abs(yf2[:N//2]), c='r', label='V_ind_u')
ax1.plot(xf, 2.0/N * np.abs(yf3[:N//2]), c='g', label='V_ind_u')
plt.legend()
ax1.set_xscale('log')
plt.xlabel('frequency (Hz)')
plt.ylabel('induced_velocity')
#savefig('frequency_(5,1)_rand_5000_All'+str(N)+str(number_of_vortices)+'.pdf')
plt.show()

# Number of samplepoints
N = 5000
# sample spacing
T = 1.0 / 2500
x = np.linspace(0.0, N*T, N)
y1 = vel_tot_mag**2
y2 = vel_tot[:,0]**2
y3 = vel_tot[:,1]**2
yf1 = scipy.fftpack.fft(y1)
yf2 = scipy.fftpack.fft(y2)
yf3 = scipy.fftpack.fft(y3)
xf = np.linspace(0.0, 1.0/(2.0*T), N/2)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(xf, 2.0/N * np.abs(yf1[:N//2]), c='b', label='V_ind_tot')
ax1.plot(xf, 2.0/N * np.abs(yf2[:N//2]), c='r', label='V_ind_u')
ax1.plot(xf, 2.0/N * np.abs(yf3[:N//2]), c='g', label='V_ind_u')
plt.legend()
ax1.set_xscale('log')
plt.xlabel('frequency (Hz)')
plt.ylabel('induced_velocity^2')
#savefig('frequency_vsquare(5,1)_rand_5000'+str(N)+str(number_of_vortices)+'.pdf')
plt.show()