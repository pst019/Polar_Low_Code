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
# PL_name_list= ['Noer', 'Smirnova', 'Yanase']
# PL_name_list= ['Noer', 'Rojo']


var_char_ind = [['theta_diff_925-500_mean250', 'min'],
            ['theta_trop_mean250', 'min'] ]


# Stype = 'system'
# system_char= 'mean'
# system_char= 'max'
# system_char='min'
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


# var, bins, PLs_larger= 'vo', np.arange(0.10, 1.0, 0.05), True
# var= 'area'
# var= 'diameter'
# var, Stype, system_char= 'step', 'system', 'final'  # Duration - this gets the last step

plotvar, plotchar='theta_trop_mean250', 'min'
plotvar, plotchar='theta_diff_925-500_mean250', 'min'

# plotvar, plotchar, PLs_larger='vo', 'max', True
plotvar, plotchar, PLs_larger='diameter', 'max', False
# plotvar= 'duration'


"""import lists"""
csv_name= 'merged_tracks_'+version+'_'+ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv')
print('now it is only compared against atlantic tracks')

tracks_ind= tracks.groupby(['ID']).max()
tracks_ind= tracks_ind.rename(columns={'step':'duration', 'vo': 'vo_max', 'area': 'area_max'})

"""get individual params for tracks"""
tracks_params= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'_params.csv').set_index(['ID', 'step'])

for vc in var_char_ind:
    var, char= vc
    if char== 'min':
        tracks_ind_param=  tracks_params[var].groupby(['ID']).min().rename( var+'_'+char)
    
    tracks_ind= tracks_ind.merge(tracks_ind_param, left_index=True, right_index=True)




PLs_ind_list=[]
for PL_name in PL_name_list:
    PLs=pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/matchPLlist/matchPLlist_"+ PL_name +f"_dist{matchdist}.csv")
    
    PLs_ind= PLs.groupby(['ID']).max()
    PLs_ind= PLs_ind.rename(columns={'step':'duration', 'vo': 'vo_max', 'area': 'area_max'})

    
    PLs_params= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/matchPLlist/matchPLlist_"+PL_name+f'_dist{matchdist}_params.csv').set_index(['ID', 'step'])

    for vc in var_char_ind:
        var, char= vc
        if char== 'min':
            PLs_ind_param=  PLs_params[var].groupby(['ID']).min().rename( var+'_'+char)
        
        PLs_ind= PLs_ind.merge(PLs_ind_param, left_index=True, right_index=True)

    PLs_ind_list += [PLs_ind]

"""create variables"""
if plotvar == 'diameter':
    tracks_ind['diameter_max']= np.sqrt(tracks_ind['area_max']) *4/np.pi
    for PLs in PLs_ind_list:
        PLs['diameter_max']= np.sqrt(PLs['area_max']) *4/np.pi




"""tracks after crit"""
tracks_ind_crit= tracks_ind[tracks_ind['theta_diff_925-500_mean250_min'] < 17] #Smirnova
tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['theta_trop_mean250_min'] < 297.5] #Gunnar

tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['vo_max'] > 0.255] #Gunnar
tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['diameter_max'] < 450] #betw Gunnar and Smir

print('Tracks remaining with criteria:', np.round(100*len(tracks_ind_crit)/len(tracks_ind), 1), '%' , len(tracks_ind_crit), 'of', len(tracks_ind) )

#system vorticity > 0.3

# alltracks['time']= pd.to_datetime(alltracks['time'])
# tracks= tracks.set_index(['ID', 'step'])


"""make the system variable"""
# if system_char== 'max':
#     tracks_ind= tracks.groupby(['ID']).max()[var]
#     tracks_crit_ind= tracks_crit.groupby(['ID']).max()[var]
#     PLs_ind_list= [PLs.groupby(['ID']).max()[var] for PLs in PLs_list]
# if system_char== 'min':
#     tracks_ind= tracks.groupby(['ID']).min()[var]
#     tracks_crit_ind= tracks_crit.groupby(['ID']).min()[var]                
#     PLs_ind_list= [PLs.groupby(['ID']).min()[var] for PLs in PLs_list]
# if system_char== 'mean':
#     tracks_ind= tracks.groupby(['ID']).mean()[var]
#     tracks_crit_ind= tracks_crit.groupby(['ID']).mean()[var]        
#     PLs_ind_list= [PLs.groupby(['ID']).mean()[var] for PLs in PLs_list]        
# if system_char== 'final':
#     tracks_ind= tracks.groupby(['ID']).last()[var]
#     tracks_crit_ind= tracks_crit.groupby(['ID']).last()[var]               
#     PLs_ind_list= [PLs.groupby(['ID']).last()[var] for PLs in PLs_list]               
# # elif system_char== 'median': Sx[var]= S.groupby(['ID']).median()[var]
# # elif system_char== 'initial': Sx[var]=S.groupby(['ID']).first()[var]







"""the plot variable"""
if plotvar == 'duration': plotvarchar= plotvar
else: plotvarchar= plotvar+'_'+plotchar


# if Stype == 'step':
tracks_var = tracks_ind[plotvarchar]
tracks_crit_var = tracks_ind_crit[plotvarchar]

PLs_var_list= [PLs[plotvarchar] for PLs in PLs_ind_list]
    
# if Stype == 'system':    
# tracks_var = tracks_ind
# tracks_crit_var = tracks_crit_ind

# PLs_var_list= PLs_ind_list


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
plt.xlabel(plotvarchar)
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
        # print(PL_name, np.round(crit, 3), 'excludes: ', np.round(len(np.where(tracks_var < crit)[0])/len(tracks_var), 2))
        print(PL_name, 'crit:', np.round(crit, 3), 'excludes:' , 100*np.round(len(np.where(tracks_var < crit)[0])/len(tracks_var), 3), '% (', len(np.where(tracks_var < crit)[0]), 'of', len(tracks_var), ')' )
        print('Of remaining tracks excludes:' , 100*np.round(len(np.where(tracks_crit_var < crit)[0])/len(tracks_crit_var), 3), '% (', len(np.where(tracks_crit_var < crit)[0]), 'of', len(tracks_crit_var), ')' )


    else:
        crit= np.percentile( PLs_var, 100- percentile)
        plt.scatter(crit, 0, marker='^', c= p[0].get_color(), s= 200 )
    
        #if PL value is larger than track value
        print(PL_name, 'crit:', np.round(crit, 3), 'excludes:' , 100*np.round(len(np.where(tracks_var > crit)[0])/len(tracks_var), 3), '% (', len(np.where(tracks_var > crit)[0]), 'of', len(tracks_var), ')' )
        print('Of remaining tracks excludes:' , 100*np.round(len(np.where(tracks_crit_var > crit)[0])/len(tracks_crit_var), 3), '% (', len(np.where(tracks_crit_var > crit)[0]), 'of', len(tracks_crit_var), ')' )


plt.ylim(ymin=0)
plt.legend()
plt.xlabel(plotvarchar)   
plt.ylabel('Normalized density')



