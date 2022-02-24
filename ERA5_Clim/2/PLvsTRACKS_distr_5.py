#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 16:20:33 2021

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/PatsOrange/'

homedir=Mediadir+'home/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import scipy.signal as signal


fignr = 1

version= "fram_run3"
ending= 'atl'

durlim= 6

start_year= 2004
end_year= 2019


matchdist= 150


percentile= 10
PLs_larger= False

PL_name= 'Rojo'
# PL_name = 'Noer'
# PL_name = 'Smirnova'
# PL_name = 'Yanase' #it has only one location per PL
PL_name_list= [PL_name]

# PL_name_list= ['Noer', 'Rojo', 'Smirnova', 'Yanase']
# PL_name_list= ['Noer', 'Smirnova', 'Yanase']
# PL_name_list= ['Noer', 'Rojo']


var_char_ind = [['theta_diff_925-500_mean250', 'min'],
            ['theta_trop_mean250', 'min'] ]


# Stype = 'system'
# plotchar= 'mean'
# plotchar= 'max'
# plotchar='min'
#plotchar='med'
#plotchar='initial'
#plotchar='final'



bins= ''

#var, bins='Season', np.arange(1999, 2020,1)-.5
#var, bins='Month', np.arange(.5, 12.6, 1)
# var, bins= 'multiple', np.arange(10)+.5
#var, bins= 'Duration', np.arange(0, 116, 5)
#var, bins= 'Distance', np.arange(0,2101, 100)  #"""travelled distance - start-end"""
#var, bins= 'Prop_speed', np.arange(0, 100.1, 4)
#var, bins= 'Diameter', np.arange(0, 851, 50)


# var, bins, PLs_larger= 'vo', np.arange(0.10, 1.0, 0.05), True
# var= 'area'
# var= 'diameter'
# var, Stype, system_char= 'step', 'system', 'final'  # Duration - this gets the last step

plotvar, plotchar='theta_trop_mean250', 'min'
plotvar, plotchar='theta_diff_925-500_mean250', 'min'
# plotvar, plotchar, bins='month', False, np.arange(0.5, 13, 1)

# plotvar, plotchar, PLs_larger='vo', 'max', True
# plotvar, plotchar, PLs_larger='diameter', 'max', False
# plotvar, plotchar= 'duration', False
plotvar, plotchar, PLs_larger= 'vo_tendency', False, True 

"""the plot variable"""
if plotchar == False: plotvarchar= plotvar
else:
    plotvarchar= plotvar+'_'+plotchar
    if [plotvar, plotchar] not in var_char_ind: var_char_ind += [[plotvar, plotchar]]



"""import lists"""
csv_name= 'merged_tracks_'+version+'_'+ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv')
print('now it is only compared against atlantic tracks')

systems= [tracks] #list that gets tracks and PLs
systems_name_list= ['tracks'] + PL_name_list
systems_ind= []

for PL_name in PL_name_list:
    PLs=pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/matchPLlist/matchPLlist_"+ PL_name +f"_dist{matchdist}.csv")
    systems += [PLs]


for system in systems:

    system_ind= system.groupby(['ID']).max()
    system_ind= system_ind.rename(columns={'step':'duration', 'vo': 'vo_max', 'area': 'area_max'})
    

    """get individual params for systems (tracks + PLs) - calculate the system_ind variable """
    system_params= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'_params.csv').set_index(['ID', 'step'])
    
    for vc in var_char_ind:
        importvar, importchar= vc
        if importchar== 'max':
            system_ind_param=  system_params[importvar].groupby(['ID']).max().rename( importvar+'_'+importchar)
        if importchar== 'min':
            system_ind_param=  system_params[importvar].groupby(['ID']).min().rename( importvar+'_'+importchar)
        if importchar== 'mean':
            system_ind_param=  system_params[importvar].groupby(['ID']).mean().rename( importvar+'_'+importchar)
        if importchar== 'median':
            system_ind_param=  system_params[importvar].groupby(['ID']).median().rename( importvar+'_'+importchar)
        if importchar== 'initial':
            system_ind_param=  system_params[importvar].groupby(['ID']).first().rename( importvar+'_'+importchar)
        if importchar== 'final':
            system_ind_param=  system_params[importvar].groupby(['ID']).last().rename( importvar+'_'+importchar)  
            
        system_ind= system_ind.merge(system_ind_param, left_index=True, right_index=True)


    """create additional variables"""
    # if plotvarchar == 'diameter_max':
    system_ind['diameter_max']= np.sqrt(system_ind['area_max']) *4/np.pi


    if plotvar== 'month':
        system[plotvar]= pd.to_datetime(system['time']).dt.month
        system_ind[plotvar]= system.groupby(['ID']).median()[plotvar]

    if plotvar == 'vo_tendency':
        system[plotvar]= np.nan
        for index in system_ind.index:
            track_now=system[system.ID == index]
            
            if len(track_now)>= 11:  window= 11
            else: window= len(track_now) -(len(track_now)+1)%2 #to get the largest uneven number smaller than the length of vorticity
        
            # track_now= track_now.assign(vo_tendency= signal.savgol_filter(track_now['vo'].values, window_length=window, polyorder=2, deriv= 1) )
            vo_tend= signal.savgol_filter(track_now['vo'].values, window_length=window, polyorder=2, deriv= 1)

            system.loc[system.ID == index, plotvar]= vo_tend
            
        system_ind[plotvar]=  system.set_index(['ID', 'step'])[plotvar].groupby(['ID']).max()


    systems_ind += [system_ind]







"""tracks after crit"""
tracks_ind_crit = systems_ind[0]
tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['theta_diff_925-500_mean250_min'] < 17] #Smirnova
tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['theta_trop_mean250_min'] < 297.5] #Gunnar

tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['vo_max'] > 0.255] #Gunnar
tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['diameter_max'] < 450] #betw Gunnar and Smir

systems_ind += [tracks_ind_crit]
systems_name_list += ['tracks after crit'] 


print('Tracks remaining with criteria:', np.round(100*len(tracks_ind_crit)/len(systems_ind[0]), 1), '%' , len(tracks_ind_crit), 'of', len(systems_ind[0]) ,'\n')













system_var_list= [system[plotvarchar] for system in systems_ind]
   
if bins == '':
    bins= np.linspace(system_var_list[0].min(), system_var_list[0].max(), 20)


"""plot histogram"""
plt.figure(fignr)
fignr+=1
plt.clf() 

for system_var, system_name in zip(system_var_list, systems_name_list):
    plt.hist(system_var, density= True, bins= bins, alpha= 0.6, label= system_name)


plt.legend()
plt.xlabel(plotvarchar)
plt.ylabel('Density')


"""plot density distribution"""
plt.figure(fignr)
fignr+=1
plt.clf() 


xplot= np.linspace(bins[0], bins[-1], 100)
 

for system_var, system_name in zip(system_var_list, systems_name_list):
    kernel= stats.gaussian_kde(system_var) # .dropna())
    p= plt.plot(xplot, kernel(xplot), label= system_name )

    if system_name in systems_name_list[1:-1]:
        if PLs_larger:
            crit= np.percentile( system_var, percentile)
            plt.scatter(crit, 0, marker='^', c= p[0].get_color(), s= 200 )
        
            #if PL value is larger than track value
            # print(PL_name, np.round(crit, 3), 'excludes: ', np.round(len(np.where(tracks_var < crit)[0])/len(tracks_var), 2))
            print(system_name, 'crit:', np.round(crit, 3), 'excludes:' , 100*np.round(len(np.where(system_var_list[0] < crit)[0])/len(system_var_list[0]), 3), '% (', len(np.where(system_var_list[0] < crit)[0]), 'of', len(system_var_list[0]), ')' )
            print('Of remaining tracks excludes:' , 100*np.round(len(np.where(system_var_list[-1] < crit)[0])/len(system_var_list[-1]), 3), '% (', len(np.where(system_var_list[-1] < crit)[0]), 'of', len(system_var_list[-1]), ')' )
    
    
        else:
            crit= np.percentile( system_var, 100- percentile)
            plt.scatter(crit, 0, marker='^', c= p[0].get_color(), s= 200 )
        
            #if PL value is larger than track value
            print(system_name, 'crit:', np.round(crit, 3), 'excludes:' , 100*np.round(len(np.where(system_var_list[0] > crit)[0])/len(system_var_list[0]), 3), '% (', len(np.where(system_var_list[0] > crit)[0]), 'of', len(system_var_list[0]), ')' )
            print('Of remaining tracks excludes:' , 100*np.round(len(np.where(system_var_list[-1] > crit)[0])/len(system_var_list[-1]), 3), '% (', len(np.where(system_var_list[-1] > crit)[0]), 'of', len(system_var_list[-1]), ')' )
    

plt.ylim(ymin=0)
plt.legend()
plt.xlabel(plotvarchar)   
plt.ylabel('Normalized density')



