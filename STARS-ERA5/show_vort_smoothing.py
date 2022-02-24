#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 10:01:45 2019

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
    Mediadir= '/run/media/pst019/PatsOrange/'

homedir=Mediadir+'home/'


import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import xarray as xr
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from windrose import WindroseAxes
import matplotlib as mpl
from scipy import stats

plt.rcParams.update({'font.size': 12})


import numpy as np
from f_useful import *
from f_STARS import *

#
save= True
#save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/vort_tend/'

#
fignr= 1


ID= '3_1'
#ID= '40_1'



lifelim= 6


"""import Stoll list"""
Stoll_imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


Stoll= pd.read_csv(Stoll_imp_dir + 'Stoll_list_noRojo.csv')
Stoll= Stoll.set_index(['ID', 'time'])

Stoll_dir_speed= pd.read_csv(Stoll_imp_dir + 'Stoll_list_dir-speed_smth1E-3.csv')
Stoll_dir_speed= Stoll_dir_speed.drop(columns='Obs')
Stoll_dir_speed= Stoll_dir_speed.set_index(['ID', 'time'])

#S_nodes= pd.read_csv(Stoll_imp_dir + 'Stoll_nodes_x'+str(xSOM)+'_y'+str(ySOM)+'.csv')
#S_nodes= S_nodes.set_index(['ID', 'time'])

#S_ERA5= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5.csv')
#S_ERA5= S_ERA5.set_index(['ID', 'time'])
#
#S_ERA5_2= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5_2.csv')
#S_ERA5_2= S_ERA5_2.set_index(['ID', 'time'])

Stoll= pd.concat([Stoll, Stoll_dir_speed], axis=1)


Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps







def Stoll_param_tend_2(Stoll, parameter='Vorticity$_{850}$ (centre)', outname= 'Vort_tend'):
    """Stoll calculate the tendency in the parameter (first derivative with Obs - like time derivative - if all no gaps would be present
    this is just the tendency from one to the next timestep - inaccurate if gaps occur
    """
    import scipy.signal

    Stoll[outname]= ""
    
    for ID in remove_dublicate(Stoll['ID']):
        vo= Stoll.loc[Stoll['ID'] == ID, parameter].values
        if len(vo)>= 11: window= 11
        else: window= len(vo) -(len(vo)+1)%2 #to get the largest uneven number smaller than the length of vorticity
        
        vodf= scipy.signal.savgol_filter(vo, window_length=window, polyorder=2, deriv=1)
        vo_smooth= scipy.signal.savgol_filter(vo, window_length=window, polyorder=2, deriv=0)

        Stoll.loc[Stoll['ID'] == ID, outname]= vodf
        Stoll.loc[Stoll['ID'] == ID, 'vo_smooth']= vo_smooth
        
    return Stoll





def calc_Obs_mature_2(S, S_ind, intensity_var='vo'):
    """calculate the observation number of the mature phase, defined as maximum in the intensity variable """

    if 'ID' not in S.columns:
        reset_index= True
        S= S.reset_index()
    else: reset_index= False

    Obs_mature_list= []
    S['Obs_mature']= 0
    
    for ID_now in S_ind.index.values:
#        print(ID_now)
        S_now= S.loc[S['ID'] == ID_now]
        obs_mature= S_now.iloc[np.argmax(S_now[intensity_var].values)].Obs
        
#        S[np.logical_and(S['ID']== ID_now, S['Obs'] == obs_mature)]['Obs_mature']= 1
        Obs_mature_list += [obs_mature]
        
#    S_ind.reindex([remove_dublicate(S.ID)])
    S_ind['Obs_mature']= Obs_mature_list
    
    for i in range(len(S_ind)):
        S.loc[np.logical_and(S['ID']== S_ind['Obs_mature'].index[i], S['Obs']== S_ind['Obs_mature'][i]), 'Obs_mature']= 1
    
    if reset_index:
        S= S.set_index(['ID', 'time'])
    
    return S, S_ind



Stoll_ind= Stoll.groupby(['ID']).mean()

Stoll, Stoll_ind= calc_Obs_mature_2(Stoll, Stoll_ind, intensity_var='vo')



#
fig= plt.figure(fignr, figsize=(7, 8))
plt.clf()

ax1 = plt.subplot(311)
ax2 = plt.subplot(312, sharex = ax1)
ax3 = plt.subplot(313, sharex = ax1)


S_ID= Stoll.loc[ID]

ax1.plot(S_ID['Obs'].values -1, S_ID.vo.values*1E5, label='Spatially-filtered vorticity')

ax2.plot(S_ID['Obs'].values -1, S_ID.vo_tendency* 1E5, c='r', label='Spatially + time-filtered vorticity tendency')


ax1.scatter(S_ID[S_ID['Obs_mature']== 1]['Obs']-1, S_ID[S_ID['Obs_mature']== 1]['vo']*1E5, label= 'Mature stage')



ax3.plot(S_ID['Obs'].values -1, 24*S_ID.vo_tendency/S_ID.vo.values, c='r', label='Spatially + time-filtered vorticity tendency')




Stoll= Stoll.reset_index()
Stoll= Stoll_param_tend_2(Stoll, parameter= 'vo')

S_ID= Stoll[Stoll.ID== ID]

ax1.plot(S_ID['Obs'].values-1, S_ID.vo_smooth*1E5, label='Spatially + time-filtered vorticity', c= 'r')
ax1.set_xlim(0, S_ID['Obs'].max()-1 )
#
ax1.legend()
ax2.legend(loc= 3)

ax1.set_ylabel('Vorticity [10$^{-5}$s$^{-1}$]')
ax2.set_ylabel('Vorticity tendency \n [10$^{-5}$s$^{-1}$ per h]')
ax2.plot(S_ID['Obs']-1, len(S_ID)*[0], c='k')

ax3.set_ylabel('Growth rate \n [day$^{-1}$]')
ax3.plot(S_ID['Obs']-1, len(S_ID)*[0], c='k')


ax3.set_xlabel('Time since genesis of the polar low [h]')



plt.tight_layout()

if save:
    save_name=savedir+ ID
    print(save_name)
    plt.savefig(save_name , bbox_inches='tight', dpi= 150)
#
