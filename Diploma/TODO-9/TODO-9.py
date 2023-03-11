#!/usr/bin/env python
# coding: utf-8

# In[6]:


import numpy as np
x = np.loadtxt('Data\Cell Center Coordinates X.txt')
c0 = np.loadtxt('Data\Cell Concentration 0.txt')
c1 = np.loadtxt('Data\Cell Concentration 1.txt')
print(x)
print(c0)
print(c1)


# In[7]:


import matplotlib.pyplot as plt
import numpy as np

plt.plot(x, c0, label='t=0')
plt.plot(x, c1, label='t=100')

plt.ylabel('concentrarion')
plt.xlabel('x')
plt.ylim([-1, 2])
plt.grid()
plt.legend()
plt.title('Concentration (y=0)')
plt.show()

