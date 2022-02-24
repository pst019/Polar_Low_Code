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


fignr = 1

version= "fram_run3"
ending= 'atl'

durlim= 6

start_year= 2004
end_year= 2019


matchdist= 150


percentile= 10
PLs_larger= False

# PL_list_name= 'Rojo'
PL_name = 'Noer'
# PL_name = 'Smirnova'
# PL_name = 'Yanase' #it has only one location per PL

# PL_name_list= ['Noer', 'Rojo', 'Smirnova', 'Yanase']
PL_name_list= ['Noer', 'Smirnova', 'Yanase']
# PL_name_list= ['Noer', 'Rojo']


Stype= 'step'

Stype = 'system'
# system_char= 'mean'
# system_char= 'max'
system_char='min'
#system_char='med'
#system_char='initial'
#system_char='final'



bins= ''

#var, bins='Season', np.arange(1999, 2020,1)-.5
#var, bins='Month', np.arange(.5, 12.6, 1)
# var, bins= 'multiple', np.arange(10)+.5
#var, bins= 'Duration', np.arange(0, 116, 5)
#var, bins= 'Distance', np.arange(0,2101, 100)  #"""travelled distance - start-end"""
#var, bins= 'Prop_speed', np.arange(0, 100.1, 4)
#var, bins= 'Diameter', np.arange(0, 851, 50)


#var_param=False #from the original track list
# var, bins, PLs_larger= 'vo', np.arange(0.10, 1.0, 0.05), True
# var= 'area'
# var= 'diameter'
# var, Stype, system_char= 'step', 'system', 'final'  # Duration - this gets the last step

# var_param= True
# var='theta_diff_925-500_mean250' #quite good at system minimum
var='theta_trop_mean250'

"""import lists"""
csv_name= 'merged_tracks_'+version+'_'+ending+f'_{start_year}-{end_year}'


if durlim: csv_name += 'durlim'+str(durlim)
tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv')
print('now it is only compared against atlantic tracks')

# if var_param:
tracks= tracks.set_index(['ID', 'step'])

tracks_params= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'_params.csv').set_index(['ID', 'step'])

tracks= tracks.merge(tracks_params, left_index=True, right_index=True)
# tracks= pd.concat([tracks, tracks_params], axis=1)
tracks= tracks.reset_index()


"""tracks after crit"""
tracks_crit= tracks[tracks['theta_diff_925-500_mean250'] < 15.7]
tracks_crit= tracks_crit[tracks_crit['theta_trop_mean250'] < 296.7]

#system vorticity > 0.3

# alltracks['time']= pd.to_datetime(alltracks['time'])
# tracks= tracks.set_index(['ID', 'step'])


# PLs= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/matchPLlist_"+ PL_name +f"_dist{matchdist}.csv")
# PLs= PLs.set_index(['ID', 'step'])

PLs_list=[]
for PL_name in PL_name_list:
    PLs=pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/matchPLlist/matchPLlist_"+ PL_name +f"_dist{matchdist}.csv").set_index(['ID', 'step'])
    
    PLs_params= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/matchPLlist/matchPLlist_"+PL_name+f'_dist{matchdist}_params.csv').set_index(['ID', 'step'])

    PLs= PLs.merge(PLs_params, left_index=True, right_index=True)

    PLs= PLs.reset_index()
    PLs_list += [PLs]

"""create variables"""
if var == 'diameter':
    tracks['diameter']= np.sqrt(tracks['area']) *4/np.pi
    tracks_crit['diameter']= np.sqrt(tracks_crit['area']) *4/np.pi
    for PLs in PLs_list:
        PLs['diameter']= np.sqrt(PLs['area']) *4/np.pi


"""make the system variable"""
if Stype == 'system':
    if system_char== 'max':
        tracks_ind= tracks.groupby(['ID']).max()[var]
        tracks_crit_ind= tracks_crit.groupby(['ID']).max()[var]
        PLs_ind_list= [PLs.groupby(['ID']).max()[var] for PLs in PLs_list]
    if system_char== 'min':
        tracks_ind= tracks.groupby(['ID']).min()[var]
        tracks_crit_ind= tracks_crit.groupby(['ID']).min()[var]                
        PLs_ind_list= [PLs.groupby(['ID']).min()[var] for PLs in PLs_list]
    if system_char== 'mean':
        tracks_ind= tracks.groupby(['ID']).mean()[var]
        tracks_crit_ind= tracks_crit.groupby(['ID']).mean()[var]        
        PLs_ind_list= [PLs.groupby(['ID']).mean()[var] for PLs in PLs_list]        
    if system_char== 'final':
        tracks_ind= tracks.groupby(['ID']).last()[var]
        tracks_crit_ind= tracks_crit.groupby(['ID']).last()[var]               
        PLs_ind_list= [PLs.groupby(['ID']).last()[var] for PLs in PLs_list]               
    # elif system_char== 'median': Sx[var]= S.groupby(['ID']).median()[var]
    # elif system_char== 'initial': Sx[var]=S.groupby(['ID']).first()[var]







"""the plot variable"""
if Stype == 'step':
    tracks_var = tracks[var]
    tracks_crit_var = tracks_crit[var]
    
    PLs_var_list= [PLs[var] for PLs in PLs_list]
    
if Stype == 'system':    
    tracks_var = tracks_ind
    tracks_crit_var = tracks_crit_ind

    PLs_var_list= PLs_ind_list


if bins == '':
    bins= np.linspace(tracks_var.min(), tracks_var.max(), 20)


"""plot histogram"""
plt.figure(fignr)
fignr+=1
plt.clf() 



plt.hist(tracks_var, density= True, bins= bins, alpha= 0.6, label= 'tracks')
plt.hist(tracks_crit_var, density= True, bins= bins, alpha= 0.6, label= 'tracks after crit')

for PLs_var, PL_name in zip(PLs_var_list, PL_name_list):
    plt.hist(PLs_var, density= True, bins= bins, alpha= 0.6, label= PL_name)


plt.legend()
plt.xlabel(var)
plt.ylabel('Density')


"""plot density distribution"""
plt.figure(fignr)
fignr+=1
plt.clf() 


xplot= np.linspace(bins[0], bins[-1], 100)
 

kernel= stats.gaussian_kde(tracks_var) # .dropna())
plt.plot(xplot, kernel(xplot), label='tracks')

kernel= stats.gaussian_kde(tracks_crit_var) # .dropna())
plt.plot(xplot, kernel(xplot), label='tracks after crit')

for PLs_var, PL_name in zip(PLs_var_list, PL_name_list):
    kernel= stats.gaussian_kde(PLs_var) # .dropna())
    p= plt.plot(xplot, kernel(xplot), label= PL_name )

    if PLs_larger:
        crit= np.percentile( PLs_var, percentile)
        plt.scatter(crit, 0, marker='^', c= p[0].get_color(), s= 200 )
    
        #if PL value is larger than track value
        print(PL_name, np.round(crit, 3), 'excludes: ', np.round(len(np.where(tracks_var < crit)[0])/len(tracks_var), 2))

    else:
        crit= np.percentile( PLs_var, 100- percentile)
        plt.scatter(crit, 0, marker='^', c= p[0].get_color(), s= 200 )
    
        #if PL value is larger than track value
        print(PL_name, 'crit:', np.round(crit, 3), 'excludes:' , 100*np.round(len(np.where(tracks_var > crit)[0])/len(tracks_var), 1), '% (', len(np.where(tracks_var > crit)[0]), 'of', len(tracks_var), ')' )


plt.ylim(ymin=0)
plt.legend()
plt.xlabel(var)   
plt.ylabel('Normalized density')



