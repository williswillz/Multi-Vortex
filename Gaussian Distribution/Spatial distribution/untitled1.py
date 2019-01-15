# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 15:17:27 2019

@author: WS1
"""

from pysces import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import scipy.fftpack
import random
import time
from scipy.stats import norm

#mu, sigma = 0, 0.1 # mean and standard deviation
#s = np.random.normal(mu, sigma, 100)

total_num = 100
#x, y = np.random.normal(0, 1, (2, total_num))
#print x
#print y
#xy = np.vstack([x, y])
#x = s
#y = s
#plt.scatter(x,y)
#plt.show()

x = np.random.uniform(0, 5, (1, total_num))
y = np.random.normal(0, 1, (1, total_num))
z = np.vstack([x, y])
print z[0,:]
print z[1,:]
#x = s
#y = s
plt.scatter(x,y)
plt.show()