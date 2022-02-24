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


fignr = 1

version= "fram_run3"
ending= 'atl'

durlim= 6

start_year= 2004
end_year= 2019


matchdist= 150


# PL_list_name, minObs= 'Rojo', 4
PL_list_name, minObs = 'Gunnar', 4
# PL_list_name, minObs = 'Smirnova', 2
# PL_list_name, minObs, ending = 'Yanase', 1, 'pac' #it has only one location per PL



Stype= 'step'

Stype = 'system'
# system_char= 'mean'
system_char= 'max'
#system_char='min'
#system_char='med'
#system_char='initial'
#system_char='final'




#var, bins='Season', np.arange(1999, 2020,1)-.5
#var, bins='Month', np.arange(.5, 12.6, 1)
# var, bins= 'multiple', np.arange(10)+.5
#var, bins= 'Duration', np.arange(0, 116, 5)
#var, bins= 'Distance', np.arange(0,2101, 100)  #"""travelled distance - start-end"""
#var, bins= 'Prop_speed', np.arange(0, 100.1, 4)
#var, bins= 'Diameter', np.arange(0, 851, 50)

var, bins= 'vo', np.arange(0.15, 1.5, 0.05)
# var= 'area'
# var= 'Vort_diameter'


"""import lists"""
csv_name= 'merged_tracks_'+version+'_'+ending+f'_{start_year}-{end_year}'


if durlim: csv_name += 'durlim'+str(durlim)
tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv')
# alltracks['time']= pd.to_datetime(alltracks['time'])
# tracks= tracks.set_index(['ID', 'step'])


PLs= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/matchPLlist_"+ PL_list_name +f"_dist{matchdist}.csv")
# PLs= PLs.set_index(['ID', 'step'])







if Stype == 'system':
# #    Sx.reset_index(inplace=True, drop=True)
    if system_char== 'max':
        PLs_ind= PLs.groupby(['ID']).max()[var]
        tracks_ind= tracks.groupby(['ID']).max()[var]
    # elif system_char== 'min': Sx[var]= S.groupby(['ID']).min()[var]
    # elif system_char== 'median': Sx[var]= S.groupby(['ID']).median()[var]
    # elif system_char== 'initial': Sx[var]=S.groupby(['ID']).first()[var]
    # elif system_char== 'final': Sx[var]=S.groupby(['ID']).last()[var]



"""create variables"""
if var == 'Vort_diameter':
    PLs['Vort_diameter']= np.sqrt(PLs['area']) *4/np.pi
    tracks['Vort_diameter']= np.sqrt(tracks['area']) *4/np.pi




"""plot"""
plt.figure(fignr)
fignr+=1
plt.clf() 

if Stype == 'step':
    plt.hist(tracks[var], density= True, bins= bins, alpha= 0.6, label= 'tracks')
    plt.hist(PLs[var], density= True, bins= bins, alpha= 0.6, label= 'PLs')

if Stype == 'system':
    plt.hist(tracks_ind, density= True, bins= bins, alpha= 0.6, label= 'tracks')
    plt.hist(PLs_ind, density= True, bins= bins, alpha= 0.6, label= 'PLs')


plt.legend()




