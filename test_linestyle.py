#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 13:48:09 2019

@author: pst019
"""

import matplotlib.pyplot as plt
import numpy as np

plt.figure(0)
plt.clf()
marker = ['-', '--', ':', '-.', (0, (2,2)), (0, (3, 2, 1, 2, 1, 2)), (0, (4, 2, 4, 2, 1, 2))]

#from cycler import cycler
#mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')

#colorlist= ['b','g','r','c','m','y','k']

colorlist= ['k','#17becf','#1f77b4','#2ca02c','#d62728','#ff7f0e','#8c564b']
colorlist= [plt.cm.Spectral(h) for h in [0, 40, 80]] + ['k'] + [plt.cm.Spectral(h) for h in [170, 210, 250]] 

x= np.linspace(0, np.pi, 100)

for i in range(7):
    plt.plot(x, np.sin(x+ i/(np.pi)), color= colorlist[i]) # ls=marker[i],