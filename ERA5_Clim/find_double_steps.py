#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 16:29:52 2019

@author: pst019
only in python 3.6
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
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
from f_carto import *
from f_useful import *
from f_ERA5_clim import *
# import scipy
# from scipy.ndimage import uniform_filter
from scipy import ndimage


fignr= 1
savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/Investigate_PL/'
# save= True
save= False


hem= 'NH'
version= ''

time_start= False #in this case the first and last index of alltracks will be used


if hem == 'NH':
    # version= "fram_run3_"
    ending= 'atl-pac'
    start_year= 1979 #2004
    end_year= 2020
    terrain_thresh= False #50
    terrain_dist = 2 # it is actually 4
    more_info= True
    split, exceed_large_to_small= True, 40 
    time_start= '2007-1-1' #the vorticity and slp is needed in this time
    time_end= '2008-1-1'
    
    Plot_fields= True
if hem == 'SH':
    # version= "fram_run1_"
    ending= 'SH'
    start_year= 1988 #2004
    end_year= 2020

    Plot_fields= True
    # time_start= '2016-7-1' #the vorticity and slp is needed in this time
    # time_end= '2016-7-31'

durlim= 6



"""get all tracks"""
# csv_dir= "derived_PLs/"
# csv_name= 'PLs-from_merged_tracks_'+version+ending+f'_{start_year}-{end_year}'

csv_dir= "mergedtracks/"
csv_name= 'merged_tracks_'+version+ending+f'_{start_year}-{end_year}'

if durlim: csv_name += 'durlim'+str(durlim)
if more_info: csv_name += '_moreinfo'

if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
if split: csv_name += f'_split{exceed_large_to_small}'

# if time_start != False: csv_name += f'_{time_start}-{time_end}'

print(csv_name)

tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/"+csv_dir+csv_name+'.csv')
tracks['time']= pd.to_datetime(tracks['time'])




if time_start != False:
    tracks.time= pd.to_datetime(tracks.time)
    tracks= tracks[(tracks['time'] >= time_start) & (tracks['time'] < time_end)]

tracks= tracks.drop(columns= 'Duration')


track_ID_list= remove_dublicate(tracks['ID'])

red_tracks = pd.DataFrame(columns=tracks.columns) #make an empty dataframe to append the tracks
red_tracks = red_tracks.set_index(['ID', 'step'])

tracks= tracks.set_index(['ID', 'step'])


"""find double steps"""
# for index in track_ID_list:
#     track=tracks.loc[index].copy()

#     # track= track.set_index('step')
#     # track['same_location']= 0 #same location as previous time step
#     # for s in range(len(track)):
#     for s in remove_dublicate(track.index):

#         track.loc[s]
#         if type(track.loc[s]) == pd.core.frame.DataFrame:
#             print(index, s)

# index= 119850307390 #has double steps
# # 119920110490 0 

# index= 119981203160 #0-20

# index= 120001206140 #0-13
# index= 120080210130 
# index= 120130100430 0
# index= 120151209450 0
# index= 220070302570 0

# tracks= tracks.reset_index()
# tracks[(tracks.ID== index) & (tracks.step == 0)] #get the row_number

# tracks= tracks.drop(row_number)


"""find double times of a track"""
for index in track_ID_list:
    track=tracks.loc[index].copy()

    track= track.drop_duplicates(subset='time').sort_values(by= 'time')
    track['ID']= index

    red_tracks= pd.concat([red_tracks, track])


red_tracks['ID']= red_tracks['ID'].astype(int)
red_tracks= S_step(red_tracks)

red_tracks= red_tracks.set_index(['ID', 'Obs'])
