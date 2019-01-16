# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 11:21:42 2019

@author: WS1
"""

import numpy as np
#from pysces import *
#import numpy as np
import matplotlib.pyplot as plt
#from pylab import *
#import scipy.fftpack
#import random
#import time
#from scipy.stats import norm
#import matplotlib.mlab as mlab
#from matplotlib.patches import Patch
from matplotlib.lines import Line2D


total_num = 500

q = np.zeros((total_num,2))
q[:,0] = np.random.uniform(0,50,total_num)
q[:,1] = np.random.normal(0,1,total_num)

core_radius = np.zeros((total_num,1))
core_radius[:,0] = np.random.uniform(0,1,total_num)

gam = np.zeros((total_num,1))
gam[:,0] = np.random.uniform(-1,1,total_num)

col = np.zeros((total_num))
col = np.array(np.where(gam<0,'r','g'))
#print col
#print core_radius

fig = plt.figure()

legend_elements = [Line2D([0], [0], marker='o', color='w', label='+ Vortex',
                          markerfacecolor='g', markersize=7),
                   Line2D([0], [0], marker='o', color='w', label='- Vortex',
                          markerfacecolor='r', markersize=7)]

# Create the figure
#fig, ax = plt.subplots()
#ax.legend(handles=legend_elements, loc='center')

ax0 = fig.add_subplot(111).set_xlabel('X')
ax0 = fig.add_subplot(111).set_ylabel('Y')
ax0 = fig.add_subplot(111).set_title('Vortices injection from X=0 using Normal distribution')
ax0 = plt.scatter(q[:,0],q[:,1], c=col[:,0],s=25*core_radius)
ax0 = fig.add_subplot(111).legend(handles=legend_elements,loc='upper right')


plt.show()   