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

write= True


version= "" #"fram_run3_"
# ending= 'atl-pac'
ending='SH180'

start_year= 1979
end_year= 2020

durlim= 6

more_info= True
post_process= True

#split, exceed_large_to_small= True, 40
# split= False

terrain_dist, terrain_thresh= 2, 50
# terrain_thresh= False



"""csv name """
csv_name= 'merged_tracks_'+version + ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
if more_info: csv_name += '_moreinfo'
if post_process: csv_name += '_post-process'
if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
print(csv_name)

# tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv')

# time_start= False #in this case the first and last index of alltracks will be used
# time_start= '2008-1-1' #the vorticity and slp is needed in this time
# time_end= '2009-1-1'

# if time_start != False:
#     tracks.time= pd.to_datetime(tracks.time)
#     tracks= tracks[(tracks['time'] >= time_start) & (tracks['time'] < time_end)]



for yi, year in enumerate(np.arange(start_year, end_year+1)):
    time_start= f'{year}-1-1' #the vorticity and slp is needed in this time
    time_end= f'{year+1}-1-1'
    tracks_params_year= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks_tomerge/"+csv_name+f'_params_{time_start}-{time_end}.csv').set_index(['ID', 'step'])
    if yi == 0:
        tracks_params_merged= tracks_params_year
    else:
        tracks_params_merged= pd.concat([tracks_params_merged, tracks_params_year])

# tracks_params_merged= tracks_params_merged.drop(columns=['theta_trop_mean250.1', 'theta_diff_500-925_mean250',
       # 'U_trop_polew_mean250', 'theta_trop_mean500', 'theta_diff_500-skt_mean250', 'U_surface_max250', ])
if ending == 'atl-pac': tracks_params_merged= tracks_params_merged.drop(columns=[
       'U_surface_max250', 'U_trop_polew_False250', 'theta_trop_mean250.1',  'SST-T_500_mean250.1', 'SST-T_500_mean250.1.1',  'U_surface_max250.1', 'U_trop_polew_False250.1'])

elif ending == 'SH': tracks_params_merged= tracks_params_merged.drop(columns=['theta_diff_500-sst_mean250.1', 'U_trop_polew_False250',
       'U_surface_max250', 'SST-T_500_mean250'])

elif ending == 'SH180': tracks_params_merged= tracks_params_merged.drop(columns=['land_dist.1', 'SST-T_500_mean250.1'])

# tracks_params_merged[tracks_params_merged['theta_diff_500-sst_mean250'] == -1000] = np.nan
tracks_params_merged.loc[tracks_params_merged['theta_diff_500-sst_mean250'] == -1000, 'theta_diff_500-sst_mean250'] = np.nan

if write:
    outfile= Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks2/"+csv_name+f'_params'
    outfile += '.csv'
    print(outfile)
    tracks_params_merged.to_csv(outfile)
