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

if user=='pst019':
    homedir= '/home/'+user+'/home/'
    Mediadir= '/media/'+user+'/1692A00D929FEF8B/'

else:
    Mediadir= '/run/media/pst019/1692A00D929FEF8B/'
'/home/'+user+'/home/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from f_useful import *
from f_STARS import *

fignr= 2


S = import_STARS()
S = merge_ERA5(S)


S_ind= STARS_individual_systems(S)


S= remove_morphs(S, 'Morphology', count= 50, how='replace', excludelist=['-', 'U'])
#S= remove_morphs(S, 'Morphology', count= 100, how='remove', excludelist=['-', 'U'])


S_ind= remove_morphs(S_ind, 'Morphology', count= 1, how='remove', excludelist=[''])

S_ind= S_ind.dropna()

S_ind['vo850_max_200']*= 1E5
S_ind['barotropic850_max_200'] = 1/S_ind['barotropic850_max_200'] * 1/(60**2)


"""correlation statistics"""
plt.figure(fignr, figsize=(7,7))
plt.clf()


#S_ind= S_ind.drop(['Obs', 'Month', 'PL_Morph_full', 'Morphology', 'Morphology', 'ID'], axis= 1)

S_ind.columns= [index.split('_')[0] for index in S_ind.columns.tolist()]
S_ind= S_ind.rename(columns={"grad":"grad_t850", "vert": "vert_shear_angle"})

corr= S_ind.corr()
#corr.style.background_gradient(cmap='coolwarm')
#plt.matshow(corr)
#corr.style.background_gradient()




def heatmap(x, y, value, fignr= 1):
    fig, ax = plt.subplots(num= fignr)
    
    # Mapping from column names to integer coordinates
    x_labels = [v for v in sorted(x.unique())]
    y_labels = [v for v in sorted(y.unique())][::-1]
    x_to_num = {p[1]:p[0] for p in enumerate(x_labels)} 
    y_to_num = {p[1]:p[0] for p in enumerate(y_labels)} 
    
    size_scale = 150
    sc= ax.scatter(
        x=x.map(x_to_num), # Use mapping for x
        y=y.map(y_to_num), # Use mapping for y
        s=value.abs() * size_scale, # Vector of square sizes, proportional to size parameter
        c= value, vmin= -1, vmax= 1, cmap= plt.cm.get_cmap('RdBu'),
        marker='s' # Use square as scatterplot marker
    )
    
    # Show column labels on the axes
    ax.set_xticks([x_to_num[v] for v in x_labels])
    ax.set_xticklabels(x_labels, rotation=45, horizontalalignment='right')
    ax.set_yticks([y_to_num[v] for v in y_labels])
    ax.set_yticklabels(y_labels)

    plt.colorbar(sc, ticks= [-1, -.5, 0, .5, 1], label= 'Correlation coeffient')

    

corr = pd.melt(corr.reset_index(), id_vars='index') # Unpivot the dataframe, so we can get pair of arrays for x and y
corr.columns = ['x', 'y', 'value']

heatmap(x=corr['x'], y=corr['y'], value=corr['value'], fignr= fignr)



plt.tight_layout()
fignr+=1



"""scatter plots"""

import seaborn as sns
#S_ind_red= S_ind.drop(['lat', 'lon', 'cp', 'msl', 'vert_shear_angle', 'sf'], axis= 1)

S_ind_red= S_ind[['Diameter', 'vo850', 'barotropic850', 'cape', 'Morphology']]
#cape is wrong
#baroctropic calculation should be checked

#S_ind_red= S_ind[['Diameter', 'U', 'blh', 'grad_t850', 'Morphology']]

S_ind_red= S_ind[['vo850', 'tp', 'slhf', 'sshf',  'blh', 'Morphology']]


plt.close(fignr)

g= sns.pairplot(S_ind_red, kind='scatter', hue='Morphology', height= 2, aspect= 1)
g.fig.set_tight_layout(1)
#g.fig.set_size_inches(5,5)


fignr+=1

"""PCA"""
##S_np= S_ind.drop(['Season', 'Obs', 'Month', 'PL_Morph_full', 'Morphology', 'ID'], axis= 1).values
##S_np= np.nan_to_num(S_np)
#
#from sklearn.decomposition import PCA, KernelPCA
#
#pca = KernelPCA(n_components=2, kernel='rbf', gamma=.5)
#
#
##pca.fit(S_np)
##X = pca.transform(S_np)
#
##X = pca.fit_transform(S_np)
#
#droplist= ['Season', 'Obs', 'Month', 'PL_Morph_full',  'Morphology', 'ID']
#
#droplist=['ID', 'PL_Morph_full', 'Morphology']
#X = pca.fit_transform(S_ind.drop(droplist, axis= 1))
##Y = pca.fit_transform(S_ind[S_ind['Morphology']== 'C'].drop(droplist, axis= 1))
#
#
#plt.clf()
#plt.scatter(X[:,0], X[:,1])
#plt.scatter(X[S_ind['Morphology']== 'C'][:,0], X[S_ind['Morphology']== 'C'][:,1])
#
#plt.show()