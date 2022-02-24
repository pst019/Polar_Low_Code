#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 15:46:14 2020

@author: pst019
"""

import time
start = time.time()


import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    Mediadir= '/media/'+user+'/PatsOrange/'
else:
    Mediadir= '/run/media/pst019/PatsOrange/'

homedir=Mediadir+'home/'


import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import xarray as xr
import pandas as pd 
import matplotlib.pyplot as plt
#import matplotlib.colors as colors
from scipy import ndimage


import numpy as np
from f_useful import *
# from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs

#
save= True
# save= False
savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/Investigate_situation/'

fignr= 2

#track_dir= Mediadir+"/data/ERA5_Clim/track_lists/fram/"

hem= 'NH'

if hem == 'NH':
    version= "fram_run3"
    ending= 'atl-pac'
    start_year= 1980
    end_year= 2019
    terrain_thresh= 50
    terrain_dist = 4
    more_info= True

if hem == 'SH':
    version= "fram_run1"
    ending= 'SH'
    start_year= 1988
    end_year= 2020
    terrain_thresh= 50
    terrain_dist=4
    more_info= True



durlim= 6




#matchdist= 150


"""import track list"""
csv_name= 'merged_tracks_'+version+'_'+ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
if more_info: csv_name += '_moreinfo'

if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'

tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv')




tracks['boundary']= 0


if hem == 'SH':
    """180E exclusion"""
    tracks.loc[np.abs(tracks.lon +.1) > 179.8, 'boundary']= 1
    

if hem == 'NH':
    """80N exclusion"""
    tracks.loc[np.abs(tracks.lat ) == 80, 'boundary']= 1
    



"""method by removing whole tracks that exceed a fraction of time at the boundary"""
tracks_ind= tracks[['ID', 'boundary']].groupby(['ID']).mean() 

thresh_boundary= 30
tracks= tracks.set_index('ID').loc[tracks_ind[tracks_ind['boundary'] <= thresh_boundary/100].index].reset_index()
tracks = tracks.drop(columns= 'boundary')

outfile= Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+f'_boundary-{thresh_boundary}.csv'
print('write:', outfile)
tracks.set_index('ID').to_csv(outfile)


"""method by removing part of the track when at least 5 consequtive time steps are along boundary"""
