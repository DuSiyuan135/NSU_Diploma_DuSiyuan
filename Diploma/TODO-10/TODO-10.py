#!/usr/bin/env python
# coding: utf-8

# In[2]:


import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline


fig,ax=plt.subplots()
# plt.xticks([])
# plt.yticks([])
plt.xlim([-30, 32])
plt.ylim([-30, 32])
my_x_ticks = np.arange(-30, 32, 2)
my_y_ticks = np.arange(-30, 32, 2)
plt.xticks(my_x_ticks)
plt.yticks(my_y_ticks)
plt.grid()

data = np.loadtxt('Data\\Cell Concentration\\k=0.txt') 

extent = (-100, 100, -100, 100 )

plt.imshow(data, extent=extent, cmap='Blues', vmin=-0, vmax=1)
plt.colorbar()


# In[3]:


import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline


fig,ax=plt.subplots()
# plt.xticks([])
# plt.yticks([])
plt.xlim([-30, 32])
plt.ylim([-30, 32])
my_x_ticks = np.arange(-30, 32, 2)
my_y_ticks = np.arange(-30, 32, 2)
plt.xticks(my_x_ticks)
plt.yticks(my_y_ticks)
plt.grid()

data=np.loadtxt('Data\\Cell Concentration\\k=100.txt') 

extent = (-100, 100, -100, 100 )

plt.imshow(data, extent=extent, cmap='Blues', vmin=-0, vmax=1)
plt.colorbar()

