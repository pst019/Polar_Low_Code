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
import time
from f_carto import *
from f_useful import *
from f_ERA5_clim import *
# import scipy
# from scipy.ndimage import uniform_filter
from scipy import ndimage

start = time.perf_counter()


fignr= 1

write= True


hem= 'NH'
version= ''

time_start= False #in this case the first and last index of alltracks will be used
track_year= False

if hem == 'NH':
    # version= "fram_run3_"
    ending= 'atl-pac'
    start_year= 1979 #2004
    end_year= 2020
    terrain_thresh= False #50
    terrain_dist = 2 # it is actually 4
    more_info= True
    split, exceed_large_to_small= False, 40 
    # time_start= '2007-1-1' #the vorticity and slp is needed in this time
    # time_end= '2007-2-1'
    # track_year= 2006    
if hem == 'SH':
    # version= "fram_run1_"
    ending= 'SH'
    start_year= 1988 #2004
    end_year= 2020

    # time_start= '2016-7-1' #the vorticity and slp is needed in this time
    # time_end= '2016-7-31'

durlim= 6

import scipy.signal as signal
threshold_vo_max_min= .25 #1E-3 1/s
exceed_large_to_small= 40 #for a split the smaller local max has to be at least this percentage larger than the local minima
window= 11

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

if track_year:
    tracks= tracks[(tracks.ID%10000000000)//1E6 == track_year]
    # tracks= tracks[(tracks.ID%10000000000)//1E4 == 200601]

tracks= tracks.drop(columns= 'Duration')


track_ID_list= remove_dublicate(tracks['ID'])

red_tracks = pd.DataFrame(columns=tracks.columns) #make an empty dataframe to append the tracks
red_tracks = red_tracks.set_index(['ID', 'step'])

tracks= tracks.set_index(['ID'])
tracks=tracks.drop(columns= 'step')

before_nr_timesteps= len(tracks)
remain_double_ts= 0
remain_location_repeaters= 0
remain_boundary_repeaters = 0
total_splits= 0

print(f'Before track loop: {round(time.perf_counter() - start,2)} seconds')





for index in track_ID_list:
    track=tracks.loc[index].copy()
    track= track.reset_index() #.set_index('step')

    """exclude double time steps"""
    track= track.drop_duplicates(subset='time').sort_values(by= 'time')
    remain_double_ts += len(track)

    # track['step']= np.arange(len(track))
    if len(track) < 6: continue
    track= track.set_index(np.arange(len(track)))
    
    """exclude location repeaters at beginning and end"""
    track['same_location']= 0 #same location as previous time step
    for s in range(1, len(track)):
        if (track.loc[s, 'lat'] == track.loc[ s-1, 'lat']) & (track.loc[s, 'lon'] == track.loc[ s-1, 'lon']): #checks for similar latitude and longitude
            track.loc[s, 'same_location'] = track.loc[s-1, 'same_location'] +1
        
    end_repeater = track.loc[len(track)-1, 'same_location']
    if end_repeater > 3:  track= track[:-end_repeater]

    begin_repeater= np.max(np.where(track['same_location'] == track.index)[0]) #np.argmax(track['same_location'] == track.index)
    if begin_repeater > 3: track= track[begin_repeater:]
    
    track= track.drop(columns= 'same_location')
    remain_location_repeaters += len(track)

    if len(track) < 6: continue
    track= track.set_index(np.arange(len(track)))

    """exclude boundary repeaters at beginning and end"""
    track['boundary']= 0

    if hem == 'SH':
        track.loc[np.abs(track.lon +.1) > 179.8, 'boundary']= 1
        track.loc[np.abs(track.lat ) == -30, 'boundary']= 1    
    
    elif hem == 'NH':
        track.loc[np.abs(track.lat ) == 80, 'boundary']= 1    
        track.loc[np.abs(track.lat ) == 30, 'boundary']= 1    
        track.loc[np.abs(track.lon ) == 80, 'boundary']= 1    
    
    for s in range(1, len(track)):
        if track.loc[s, 'boundary'] != 0:
            track.loc[s, 'boundary'] += track.loc[s-1, 'boundary']
    
    end_repeater = track.loc[len(track)-1, 'boundary']
    if end_repeater > 5:  track= track[:-end_repeater]

    if len(track) < 6: continue
    if track.loc[0, 'boundary'] == 1:
        begin_repeater= np.max(np.where(track['boundary']-1 == track.index)[0])
        if begin_repeater > 5: track= track[begin_repeater:]
    
    track= track.drop(columns= 'boundary')
    remain_boundary_repeaters += len(track)
    track.set_index(np.arange(len(track)))

  

    """split track"""
    track.ID *= 10
    index *= 10
    split_this_track= 0
    
    if len(track)>= window:
    # else: window= len(vo) -(len(vo)+1)%2 #to get the largest uneven number smaller than the length of vorticity
    
        track= track.assign(vo_smooth= signal.savgol_filter(track['vo'].values, window_length=window, polyorder=2) )
        
        local_max_0= signal.argrelextrema(track['vo_smooth'].where(track['vo_smooth'] > threshold_vo_max_min).values, np.greater)[0]
        
        if len(local_max_0) >= 2:
        
            # local_min= signal.argrelextrema(track['vo_smooth'].where(track['vo_smooth'] < threshold_vo_max_min).values, np.less)[0]
            local_min= signal.argrelextrema(track['vo_smooth'].values, np.less)[0]
            
            if len(local_min) >= 1:
                #to exclude "Wendepunkte" - small local max and min
                local_max= np.array([lmax for lmax in local_max_0 if np.min(np.abs(lmax - local_min)) > 3])
                local_min= np.array([lmin for lmin in local_min if np.min(np.abs(lmin - local_max_0)) > 3])
    
                if len(local_max) >= 2:
                    #the local min has to between local max
                    local_min= local_min[(local_min > local_max[0]) &(local_min < local_max[-1])]
        
                    #for each local min - start from the end since then the time steps can be replaced without influencing the next ones
                    nsplits=0
                    lmin_split=[]
                    for lmin in local_min:
                    # lmin= local_min[0]
                        vo_min= track.iloc[lmin]['vo_smooth']
                    
                        loc_max_around= [local_max[local_max < lmin][-1] , local_max[local_max > lmin][0] ]
                        vo_smaller_max= track.iloc[loc_max_around]['vo_smooth'].min()
                    
                        if vo_smaller_max/vo_min > (1 + exceed_large_to_small/100):
                            print('fraction large to small', vo_smaller_max/vo_min)
                            print('split track', index)
                            split_this_track += 1
                            
                            track.loc[(track.ID == index+nsplits) & (track.index >= lmin), 'ID'] += 1 
                            
                            nsplits += 1
                            total_splits += 1
                            lmin_split += [lmin]

    """make the steps for the split track  -have to check this, test case: 120060114270"""
    if split_this_track > 0:
        for split_track_index in remove_dublicate(track.ID):
            split_track= track.loc[track.ID== split_track_index ].copy()

            if len(split_track) >= 6:
                split_track['step']= np.arange(len(split_track))
                split_track= split_track.set_index(['ID', 'step'])
                
                red_tracks= red_tracks.append(split_track)

                
    else:            
        track['step']= np.arange(len(track))
        track= track.set_index(['ID', 'step'])

        red_tracks= red_tracks.append(track)


print(f'After track loop: {round(time.perf_counter() - start,2)} seconds')

# red_tracks= S_step(red_tracks)

# red_tracks['ID']= red_tracks['ID'].astype(int)
# red_tracks['step']= red_tracks['step'].astype(int)

# red_tracks= red_tracks.set_index(['ID', 'step'])


if write:
    outfile= Mediadir+"/data/ERA5_Clim/track_lists/"+csv_dir+csv_name+'_post-process.csv'
    print('write:', outfile)
    red_tracks.to_csv(outfile) 
    
print('\n Nr time steps prior:', before_nr_timesteps)
print('Remain after double timestep exclusion:', remain_double_ts)
print('remain after location repeater exclusion:', remain_location_repeaters)
print('remain after boundary repeater exclusion:', remain_boundary_repeaters)
print('remain in list:', len(red_tracks))
    
print(f'\n Finished in {round(time.perf_counter() - start,2)} seconds')