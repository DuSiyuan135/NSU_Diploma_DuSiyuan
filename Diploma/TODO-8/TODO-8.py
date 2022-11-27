#!/usr/bin/env python
# coding: utf-8

# In[4]:


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

concentration = np.zeros((20,20))
concentration[0:20,5:15]=1

key_list = np.linspace(-19, 19, 20)
ax = sns.heatmap(pd.DataFrame(concentration, columns=key_list, index=key_list), annot=True, vmin=0, vmax=1)
ax.set_title('concentration', fontsize=18)
plt.show()


# In[5]:


concentration = np.zeros((50,50))
concentration[0:50,20:30]=1

key_list = np.linspace(-49, 49, 50)
ax = sns.heatmap(pd.DataFrame(concentration, columns=key_list, index=key_list), annot=True, vmin=0, vmax=1)
ax.set_title('concentration', fontsize=18)
plt.show()


# In[2]:


concentration = np.zeros((100,100))
concentration[0:100,45:55]=1

key_list = np.linspace(-99, 99, 100)
ax = sns.heatmap(pd.DataFrame(concentration, columns=key_list, index=key_list), annot=True, vmin=0, vmax=1)
ax.set_title('concentration', fontsize=18)
plt.show()


# In[ ]:




