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
    Mediadir= '/media/'+user+'/PatsOrange/'

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

save= False

S = import_STARS()
S = merge_ERA5(S)

#S= remove_morphs(S, 'Morphology', count= 50, how='replace', excludelist=['-', 'U'])
S= remove_morphs(S, 'Morphology', count= 100, how='remove', excludelist=['-', 'U'])

S_ind= STARS_individual_systems(S)


#S_ind= remove_morphs(S_ind, 'Morphology', count= 1, how='remove', excludelist=[''])
#S_ind= S_ind.dropna()


"""aditional columns"""
#Duration
starttime= S.loc[S['Obs'] == 1].groupby(['ID']).first()['time']
S_ind= pd.concat([S_ind, starttime], axis= 1, sort=True)
S_ind= S_ind.rename(columns={"time":"starttime"})

#endtime
endobs= S_ind['Obs']*2 -1
endobs.loc[(endobs - endobs.astype(int)) > .01]= np.nan

endtime= S[np.logical_and(S['ID'] == endobs.index[0], S['Obs'] == endobs[0])][['ID', 'time']]
for i in range(1, len(endobs)):
    a= S[np.logical_and(S['ID'] == endobs.index[i], S['Obs'] == endobs[i])][['ID', 'time']]
    endtime= pd.concat([endtime, a])

endtime= endtime.set_index('ID')

S_ind= pd.concat([S_ind, endtime], axis= 1, sort=True)
S_ind= S_ind.rename(columns={"time":"endtime"})   

S_ind['Duration']= np.round((S_ind['endtime'] - S_ind['starttime'])/np.timedelta64(1,'h'))

#Distance
#S_ind= S_ind.rename(columns={"lat":"avg_lat", "lon":"avg_lon"})
#
#startloc= S.loc[S['Obs'] == 1].groupby(['ID']).first()[['lat', 'lon']]
#S_ind= pd.concat([S_ind, startloc], axis= 1, sort=True)
#S_ind= S_ind.rename(columns={"lat":"start_lat", "lon":"start_lon"})
#
#endobs= S_ind['Obs']*2 -1
#endobs.loc[(endobs - endobs.astype(int)) > .01]= np.nan
#
#endloc= S[np.logical_and(S['ID'] == endobs.index[0], S['Obs'] == endobs[0])][['ID', 'lat', 'lon']]
#for i in range(1, len(endobs)):
#    a= S[np.logical_and(S['ID'] == endobs.index[i], S['Obs'] == endobs[i])][['ID', 'lat', 'lon']]
#    endloc= pd.concat([endloc, a])
#
##endtime.rename(columns={"time":"endtime"})
#endloc= endloc.set_index('ID')
#
#S_ind= pd.concat([S_ind, endloc], axis= 1, sort=True)
#S_ind= S_ind.rename(columns={"lat":"end_lat", "lon":"end_lon"})
#
#S_ind['Propagation distance']= 110* np.sqrt( (S_ind['start_lat']- S_ind['end_lat'])**2+ (np.cos(np.deg2rad(0.5*(S_ind['start_lat']+ S_ind['end_lat'])))* (S_ind['start_lon']- S_ind['end_lon']))**2)



"""remove columns"""
S_ind= S_ind.drop(columns=['Obs', 'Month', 'starttime', 'endtime', #'start_lat', 'end_lat', 'avg_lat', 'avg_lon', 'start_lon', 'end_lon', 
                           'PL_Morph_full', 
                           'barotropic850_max_200','barotropic_efolding_filter10_850_min_200'])


S_ind['vo850_max_200']*= 1E5
#S_ind['barotropic850_max_200'] = 1/S_ind['barotropic850_max_200'] * 1/(60**2)


    
    
    
"""correlation statistics"""
#plt.figure(fignr, figsize=(9,7))
#plt.clf()


#Sx, Stype= S, 'step'
Sx, Stype= S_ind, 'system'

#S_ind= S_ind.drop(['Obs', , 'PL_Morph_full', 'Morphology', 'ID'], axis= 1)

#S_ind.columns= [index.split('_')[0] for index in S_ind.columns.tolist()]
    
Sx= Sx.rename(columns={
        'U_400': '10m Wind Speed (max)',
       'vo850_max_200': 'Vorticity$_{850}$ (max)',
       'msl_min_50': 'Sea-level pressure (min)',
       'blh_max_250': 'Boundary layer height (max)',
       'cape_max_250': 'Convective avail. pot. en. (max)', 
        'skt_med_200': 'Skin temperature (med)',
       'skt-t500_max_200': 'SST- T$_{500}$ (max)', 
       'skt-t700_max_200': 'SST -T$_{700}$ (max)',
       'tp_mean_300': 'Total precipitation (mean) ',
       'cp_mean_300': 'Convective precipitation (mean)',
       'sf_mean_300': 'Snow fall (mean)',
       'sshf_mean_300': 'Sensible heat flux (mean)',
       'slhf_mean_300': 'Latent heat flux (mean)',
       'grad_t850_max_300': 'Grad T$_{850}$ (max)',
       'baroclinic_efolding_filter4_925-700_min_200': 'Baroclinic efolding time (min)' ,
       'barotropic_efolding_filter4_850_min_200': 'Barotropic efolding time (min)',
       'vert_shear_angle925-700_mean_200': 'Vertical shear angle (mean)'
       })    

if Stype== 'system':     
    Sx= Sx[ ['Morphology', 'lat', 'Diameter', 'Duration', '10m Wind Speed (max)', 'Vorticity$_{850}$ (max)',
           'Sea-level pressure (min)', 'Boundary layer height (max)', 'Convective avail. pot. en. (max)', 
            'Skin temperature (med)', 'SST- T$_{500}$ (max)', 'SST -T$_{700}$ (max)',
           'Total precipitation (mean) ', 'Convective precipitation (mean)', 'Snow fall (mean)',
           'Sensible heat flux (mean)', 'Latent heat flux (mean)',
           'Grad T$_{850}$ (max)', 'Baroclinic efolding time (min)' , 'Barotropic efolding time (min)', 'Vertical shear angle (mean)'] ]

if Stype== 'step': #no Duration
    Sx= Sx[ ['Morphology', 'lat', 'Diameter', '10m Wind Speed (max)', 'Vorticity$_{850}$ (max)',
           'Sea-level pressure (min)', 'Boundary layer height (max)', 'Convective avail. pot. en. (max)', 
            'Skin temperature (med)', 'SST- T$_{500}$ (max)', 'SST -T$_{700}$ (max)',
           'Total precipitation (mean) ', 'Convective precipitation (mean)', 'Snow fall (mean)',
           'Sensible heat flux (mean)', 'Latent heat flux (mean)',
           'Grad T$_{850}$ (max)', 'Baroclinic efolding time (min)' , 'Barotropic efolding time (min)', 'Vertical shear angle (mean)'] ]

    
#corr= Sx.corr()





"""why not making this with the standard EOF package???"""


"""PCA"""
S_np= S_ind.drop(['Morphology', 'ID'], axis= 1).values
S_np= np.nan_to_num(S_np)

from sklearn.preprocessing import StandardScaler

S_np= StandardScaler().fit_transform(S_np)

from sklearn.decomposition import PCA, KernelPCA

#pca = KernelPCA(n_components=2, kernel='rbf', gamma=.5)
pca = PCA()

#pca.fit(S_np)
#X = pca.transform(S_np)

#X = pca.fit_transform(S_np)

#droplist= ['Season', 'Obs', 'Month', 'PL_Morph_full', 'Morphology', 'ID']

#droplist=['ID', 'PL_Morph_full', 'Morphology']

#droplist= ['Morphology', 'ID']
#X = pca.fit_transform(S_ind.drop(droplist, axis= 1))
X = pca.fit_transform(S_np)
#X = pca.fit(S_np)


print('Explained Variance Ratio: ', pca.explained_variance_ratio_)

#Y = pca.fit_transform(S_ind[S_ind['Morphology']== 'C'].drop(droplist, axis= 1))

plt.figure(fignr)
fignr+= 1
plt.clf()
plt.scatter(X[:,0], X[:,1])
#plt.scatter(X[S_ind['Morphology']== 'C'][:,0], X[S_ind['Morphology']== 'C'][:,1])

plt.show()