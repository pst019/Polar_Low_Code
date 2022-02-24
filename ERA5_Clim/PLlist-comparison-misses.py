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

from f_useful import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import scipy.signal as signal
import matplotlib as mpl
from cycler import cycler

save=True
# save= False
savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/PLvsTRACKS_distr_2/'

fignr = 4

version= "" #"fram_run3_"
ending= 'atl-pac'

start_year= 1979
end_year= 2020

durlim= 6
matchdist= 150

# time_diff= False
time_diff= 24 #does not include the matched track if the matched track starts more than time_diff earlier or ends more than time_diff later

more_info= True
post_process= True


terrain_dist, terrain_thresh= 2, 50
# terrain_thresh= False



# PL_name_list= ['Rojo', 'Smirnova']
# PL_name_list= ['Rojo', 'Noer']
# PL_name_list= ['Rojo', 'Golubkin']
# PL_name_list= ['Golubkin', 'Smirnova'] #have no time period in common


# PL_name = 'Noer'

PL_lists= []

"""import lists"""
print('two different track lists')
csv_name= 'merged_tracks_'+version + ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
if more_info: csv_name += '_moreinfo'
if post_process: csv_name += '_post-process'
if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
# if split: csv_name += f'_split{exceed_large_to_small}'
print(csv_name)



for PL_name in PL_name_list:
    PLfile= Mediadir+"/data/ERA5_Clim/track_lists/matchPLlist2/matchPLlist_"+ PL_name +f"_dist{matchdist}"
    if time_diff and PL_name != 'Yanase': PLfile += f'_timediff{time_diff}'
    PLfile += f'_{csv_name}.csv'
    
    PLs=pd.read_csv(PLfile)
    PLs.time= pd.to_datetime(PLs.time)
    PL_lists += [PLs]

PL_lists[0]= PL_lists[0][(PL_lists[0]['time'] >= PL_lists[1].time.min()) & (PL_lists[0]['time'] <= PL_lists[1].time.max())]
PL_lists[1]= PL_lists[1][(PL_lists[1]['time'] >= PL_lists[0].time.min()) & (PL_lists[1]['time'] <= PL_lists[0].time.max())]

# if PL_name_list== ['Rojo', 'Smirnova']:
    # PL_lists[0]= PL_lists[0][(PL_lists[0]['lat'] >= 70) & (PL_lists[0]['lat'] <= PL_lists[1].lat.max())]
    # PL_lists[1]= PL_lists[1][(PL_lists[1]['lat'] >= 70) & (PL_lists[1]['lat'] <= PL_lists[0].lat.max())]
# 
# else:
PL_lists[0]= PL_lists[0][(PL_lists[0]['lat'] >= PL_lists[1].lat.min()) & (PL_lists[0]['lat'] <= PL_lists[1].lat.max())]
PL_lists[1]= PL_lists[1][(PL_lists[1]['lat'] >= PL_lists[0].lat.min()) & (PL_lists[1]['lat'] <= PL_lists[0].lat.max())]

PL_lists[0]= PL_lists[0][(PL_lists[0]['lon'] >= PL_lists[1].lon.min()) & (PL_lists[0]['lon'] <= PL_lists[1].lon.max())]
PL_lists[1]= PL_lists[1][(PL_lists[1]['lon'] >= PL_lists[0].lon.min()) & (PL_lists[1]['lon'] <= PL_lists[0].lon.max())]

ID0= set(remove_dublicate(PL_lists[0].ID))
ID1= set(remove_dublicate(PL_lists[1].ID))

print(PL_name_list[0], len(ID0))
print(PL_name_list[1], len(ID1))

print('In both', len(ID0 & ID1), '. Of', PL_name_list[1], np.round(100*len(ID0 & ID1)/len(ID1)),'%') 
print(PL_name_list[1],' False Positives: ', len(ID1- ID0), np.round(100*len(ID1- ID0)/len(ID1)),'%' )
print(PL_name_list[1],' Misses: ', len(ID0- ID1), np.round(100*len(ID0- ID1)/len(ID0)),'%' )