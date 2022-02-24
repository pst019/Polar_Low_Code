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
import scipy.signal as signal

# index = 2008011597 *10

fignr= 1
# savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/Investigate_PL/'
# save= True

# 

version= "fram_run3"
ending= 'atl-pac'

start_year= 1980
end_year= 2019

durlim= 6

terrain_dist, terrain_thresh= 4, 50
# terrain_thresh= False

Plot= True
# Plot= False

write= True
# write= False

threshold_vo_max_min= .25 #1E-3 1/s
exceed_large_to_small= 50 #for a split the smaller local max has to be at least this percentage larger than the local minima

# # PL_list_name, minObs= 'Rojo', 4
# PL_list_name, minObs = 'Gunnar', 4
# # PL_list_name, minObs = 'Smirnova', 2
# # PL_list_name, minObs, ending = 'Yanase', 1, 'pac' #it has only one location per PL


# index= 77 #PLnr

# time_sel= 'orig_middle' #middle time step of the original PL for the PL list
# time_sel= 'orig_start'
# time_sel= 'orig_end'

# # time_sel= 'matched_start'
# # time_sel= 'matched_end'
# # time_sel= 'matched_middle'
# time_sel= 'matched_2_3'


time_start= False #in this case the first and last index of alltracks will be used
# time_start= '2008-1-1' #the vorticity and slp is needed in this time
# time_end= '2009-2-1'





"""get all tracks"""
csv_name= 'merged_tracks_'+version+'_'+ending+f'_{start_year}-{end_year}'

if durlim: csv_name += 'durlim'+str(durlim)
if more_info: csv_name += '_moreinfo'

if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'

print(csv_name)

tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv')
tracks['time']= pd.to_datetime(tracks['time'])

tracks.ID *= 10
# tracks.vo *= 1E2

# tracks= tracks.set_index(['ID', 'step'])

if time_start != False:
    tracks= tracks[(tracks['time'] >= time_start) & (tracks['time'] < time_end)]


total_splits= 0
track_ID_list= remove_dublicate(tracks['ID'])

# track_ID_list= [12008011597 *10]


for index in track_ID_list:
    track_now=tracks[tracks.ID == index]
    
    if len(track_now)>= 11:
        window= 11
    # else: window= len(vo) -(len(vo)+1)%2 #to get the largest uneven number smaller than the length of vorticity

        track_now= track_now.assign(vo_smooth= signal.savgol_filter(track_now['vo'].values, window_length=window, polyorder=2) )
        
        
        
        local_max= signal.argrelextrema(track_now['vo_smooth'].where(track_now['vo_smooth'] > threshold_vo_max_min).values, np.greater)[0]
        
        if len(local_max) >= 2:
        
            local_min= signal.argrelextrema(track_now['vo_smooth'].where(track_now['vo_smooth'] < threshold_vo_max_min).values, np.less)[0]
            
            local_min= local_min[(local_min > local_max[0]) &(local_min < local_max[-1])]
            if Plot:
                fig = plt.figure(fignr)
                plt.clf()
                plt.plot(track_now.time, track_now['vo'].values)
                plt.plot(track_now.time, track_now['vo_smooth'].values)
                plt.plot(track_now.iloc[local_max]['time'], track_now.iloc[local_max]['vo_smooth'], 'o')
                plt.plot(track_now.iloc[local_min]['time'], track_now.iloc[local_min]['vo_smooth'], 'v')
            
            
            #for each local min - start from the end since then the time steps can be replaced without influencing the next ones
            nsplits=0
            lmin_split=[]
            for lmin in local_min:
            # lmin= local_min[0]
                vo_min= track_now.iloc[lmin]['vo_smooth']
            
                loc_max_around= [local_max[local_max < lmin][-1] , local_max[local_max > lmin][0] ]
                vo_smaller_max= track_now.iloc[loc_max_around]['vo_smooth'].min()
            
                if vo_smaller_max/vo_min > (1 + exceed_large_to_small/100):
                    print('fraction large to small', vo_smaller_max/vo_min)
                    print('split track', index)
            
                    tracks.loc[(tracks.ID == index+nsplits) & (tracks.step >= lmin), 'ID'] += 1 
                    
                    nsplits += 1
                    total_splits += 1
                    lmin_split += [lmin]
                    # index += 1
                    # tracks.loc[tracks.ID == index+1, 'step'] -= lmin
            
                    # local_min -= lmin
            
            for li, lmin in enumerate(lmin_split):
                tracks.loc[(tracks.ID == index+li+1), 'step'] -= lmin 
            
            
            if Plot:
                tracks_after_split= remove_dublicate(tracks.loc[np.round(tracks.ID, -1) == index, 'ID'])
                
                for track_index in tracks_after_split:
                    track_plot= tracks.loc[tracks.ID == track_index]
                    plt.plot(track_plot.time, track_plot['vo'].values, 'x', label= track_index)
                
                
                plt.legend()


print('splits', total_splits, 'of' , len(track_ID_list))

if Plot == False:
    if write:
        print('write')
        outfile= Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+f'_split{exceed_large_to_small}.csv'
        print(outfile)
        tracks.set_index('ID').to_csv(outfile)
            
