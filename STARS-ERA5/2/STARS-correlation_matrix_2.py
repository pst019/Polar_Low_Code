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
plt.figure(fignr, figsize=(9,7))
plt.clf()


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

    
corr= Sx.corr()




def heatmap(x, y, value, fignr= 1, sortalph= True, xlabel= 'top'):
    """ sortalph= [True, False] - True -the columns should be ordered alphabetically, False - the column order of x, y is taken
    xlabel= ['top', 'bootom'] - the position of the xticks and labels """

#    if xlabel =='top':
#        plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
#        plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True

    fig, ax = plt.subplots(num= fignr)
    
    if xlabel =='top':

        ax.tick_params(bottom=False, top=True, labelbottom=False, labeltop=True)
    
#     Mapping from column names to integer coordinates
    if sortalph:
        x_labels = [v for v in sorted(x.unique())]
        y_labels = [v for v in sorted(y.unique())][::-1]
    else:
        x_labels = [v for v in x.unique()]
        y_labels = [v for v in y.unique()][::-1]    
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
    if xlabel =='top': ax.set_xticklabels(x_labels, rotation=45, horizontalalignment='left')
    if xlabel =='bottom': ax.set_xticklabels(x_labels, rotation=45, horizontalalignment='right')
    
    ax.set_yticks([y_to_num[v] for v in y_labels])
    ax.set_yticklabels(y_labels)

    plt.colorbar(sc, ticks= [-1, -.5, 0, .5, 1], label= 'Correlation coeffient')

    

corr = pd.melt(corr.reset_index(), id_vars='index') # Unpivot the dataframe, so we can get pair of arrays for x and y
corr.columns = ['x', 'y', 'value']
#
heatmap(x=corr['x'], y=corr['y'], value=corr['value'], fignr= fignr, sortalph= False)



plt.tight_layout()
fignr+=1


if save:
#    if Stype== 'system': Stype += '_'+system_char
    
    savedir= homedir+ 'Polar_Low/STARS-ana/Figs/2/'+'Correlation-Matrix_'+Stype
    print(savedir)
    plt.savefig(savedir, bbox_inches='tight')



"""scatter plots"""

import seaborn as sns
#S_ind_red= S_ind.drop(['lat', 'lon', 'cp', 'msl', 'vert_shear_angle', 'sf'], axis= 1)

#S_ind_red= S_ind[['Diameter', 'vo850', 'barotropic850', 'cape', 'Morphology']]
#cape is wrong
#baroctropic calculation should be checked

#S_ind_red= S_ind[['Diameter', 'U', 'blh', 'grad_t850', 'Morphology']]

#S_ind_red= S_ind[['vo850', 'tp', 'slhf', 'sshf',  'blh', 'Morphology']]

Sx_red= Sx[['Morphology',
        #'lat', 'Diameter', '10m Wind Speed (max)',
            'Vorticity$_{850}$ (max)',
#           'Sea-level pressure (min)',
           'Boundary layer height (max)',
#           'Convective avail. pot. en. (max)', 
#            'Skin temperature (med)', 'SST- T$_{500}$ (max)', 'SST -T$_{700}$ (max)',
#           'Total precipitation (mean) ', 'Convective precipitation (mean)', 'Snow fall (mean)',
           'Sensible heat flux (mean)',
#           'Latent heat flux (mean)',
#           'Grad T$_{850}$ (max)',
           'Baroclinic efolding time (min)' , 
           'Barotropic efolding time (min)',
#           'Vertical shear angle (mean)'
           ] ]


Sx_red= Sx_red.rename(columns={
       'Vorticity$_{850}$ (max)': '$\zeta_{850} [10^{-5}$1/s]',
       'Boundary layer height (max)': 'BLH [m]',
#       'cape_max_250': 'Convective avail. pot. en. (max)', 
#        'skt_med_200': 'Skin temperature (med)',
#       'skt-t500_max_200': 'SST- T$_{500}$ (max)', 
#       'skt-t700_max_200': 'SST -T$_{700}$ (max)',
#       'tp_mean_300': 'Total precipitation (mean) ',
#       'cp_mean_300': 'Convective precipitation (mean)',
#       'sf_mean_300': 'Snow fall (mean)',
       'Sensible heat flux (mean)': 'Sens FLX [W/m$^2$]',
#       'slhf_mean_300': 'Latent heat flux (mean)',
#       'grad_t850_max_300': 'Grad T$_{850}$ (max)',
        'Baroclinic efolding time (min)': '$t_{baroclin}$ [h]',
        'Barotropic efolding time (min)': '$t_{barotrop} [h]$'
#       'vert_shear_angle925-700_mean_200': 'Vertical shear angle (mean)'
       })    


plt.close(fignr)

g= sns.pairplot(Sx_red, kind='scatter', hue='Morphology', height= 2, aspect= 1)
#g.fig.set_tight_layout(1)
#g.fig.set_size_inches(5,5)

#handles = g._legend_data.values()
#labels = g._legend_data.keys()
#
#g.fig.legend(handles=handles, labels=labels, loc='upper center', ncol=5)
#g.fig.subplots_adjust(top=0.92, bottom=0.08)

if save:
#    if Stype== 'system': Stype += '_'+system_char
    
    savedir= homedir+ 'Polar_Low/STARS-ana/Figs/2/'+'Scatterplot-Matrix_'+Stype
    print(savedir)
    plt.savefig(savedir, bbox_inches='tight')



fignr+=1




"""PCA"""
##S_np= S_ind.drop(['Season', 'Obs', 'Month', 'PL_Morph_full',  'Morphology', 'ID'], axis= 1).values
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