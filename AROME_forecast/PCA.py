#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 18:05:34 2019

@author: pst019
"""

import time
start = time.time()


import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from f_useful import *
from f_STARS import *

fignr= 1


S = import_STARS()
S = merge_ERA5(S)


S_ind= STARS_individual_systems(S)


S= remove_morphs(S, 'Morphology', count= 50, how='replace', excludelist=['-', 'U'])
#S= remove_morphs(S, 'Morphology', count= 100, how='remove', excludelist=['-', 'U'])


S_ind= remove_morphs(S_ind, 'Morphology', count= 1, how='remove', excludelist=[''])

S_ind= S_ind.dropna()
droplist= ['Season', 'Obs', 'Month', 'PL_Morph_full', 'PL_Morph', 'Morphology', 'ID']

#S_np= S_ind.drop(['Season', 'Obs', 'Month', 'PL_Morph_full', 'PL_Morph', 'Morphology', 'ID'], axis= 1).values
#S_np= np.nan_to_num(S_np)

from sklearn.decomposition import PCA, KernelPCA

pca = KernelPCA(n_components=2, kernel='rbf', gamma=.5)


#pca.fit(S_np)
#X = pca.transform(S_np)

#X = pca.fit_transform(S_np)


X = pca.fit_transform(S_ind.drop(droplist, axis= 1))
#Y = pca.fit_transform(S_ind[S_ind['Morphology']== 'C'].drop(droplist, axis= 1))


plt.clf()
plt.scatter(X[:,0], X[:,1])
plt.scatter(X[S_ind['Morphology']== 'C'][:,0], X[S_ind['Morphology']== 'C'][:,1])

plt.show()