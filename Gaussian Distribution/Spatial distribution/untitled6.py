# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 09:03:04 2019

@author: WS1
"""
import numpy as np
total_num = 10
core_radius = np.zeros((total_num,1))
a = 0

for i in range (0, total_num):
#    q[a,0] = np.random.uniform(0,50,1)
#    q[a,1] = np.random.normal(0,1,1)
#    q[a,1] = np.random.rayleigh(0,1)
#    xvort1[a,0] = i
#    xvort1[a,1] = 0
#    vel1[:,:] = panels.induced_velocity_single(q, xvort1, random.uniform(-1,1))
    core_radius[a] = np.random.uniform(0,1,1)
    a = a + 1

#print xvort1  
#print vel1


print core_radius[:]