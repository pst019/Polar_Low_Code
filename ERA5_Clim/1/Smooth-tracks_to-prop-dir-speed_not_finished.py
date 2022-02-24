#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 10:01:45 2019

@author: pst019
"""


import time
start = time.time()


import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    homedir= '/home/'+user+'/home/'
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/1692A00D929FEF8B/'
'/home/'+user+'/home/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import csaps #pip install csaps==0.10.0 had to add "np.pad(u, pad_width, mode='constant'), axis=0) " in _sspumv.py

from f_useful import *
from f_meteo import *
from f_STARS import *
from datetime import datetime
import great_circle_calculator.great_circle_calculator as gcc # pip install great-circle-calculator





smooth_parameter, string_smooth_param= 0.001, '1E-3'



version= "fram_run3"
ending= 'atl-pac'

start_year= 1980
end_year= 2019

durlim= 6
matchdist= 150

split, exceed_large_to_small= True, 50
split= False

terrain_excl, terrain_dist, terrain_thresh= True, 2, 100
# terrain_excl= False

"""import lists"""
csv_name= 'merged_tracks_'+version+'_'+ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
if terrain_excl: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
if split: csv_name += f'_split{exceed_large_to_small}'
print(csv_name)

tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv')

time_start= False #in this case the first and last index of alltracks will be used
time_start= '2008-1-1' #the vorticity and slp is needed in this time
time_end= '2009-2-1'

if time_start != False:
    tracks= tracks[(tracks['time'] >= time_start) & (tracks['time'] < time_end)]




tracks['prop_dir']= np.nan
tracks['prop_speed']= np.nan


for ID in remove_dublicate(tracks['ID']):
#    print(i) #, ID)

#    S_now= S.iloc[i]
    track = tracks[tracks['ID'] == ID]
    
    # S_t = pd.to_datetime(track.time) - pd.to_datetime(track.time.iloc[0])
    # S_t /= np.timedelta64(1, 'h')
    # S_t = S_t.values
    """get the center position, the steering vector"""

    #obtain the orientation of the smoothed track
    track_data= np.array([track.lon.values, track.lat.values])
    t= np.arange(len(track_data))
    
    
    sp = csaps.CubicSmoothingSpline(t, track_data, smooth= smooth_parameter)  #include mode= 'constant' into np.pad
    Sm_lon, Sm_lat= sp(t)
#        ax.scatter(Sm_lon, Sm_lat, transform=ccrs.PlateCarree(), s= 50, marker="x", color= 'k', zorder= 2)
#        ax.plot(Sm_lon, Sm_lat, transform=ccrs.PlateCarree(), color= 'k')
#        ax.scatter(Sm_lon[Obs-1], Sm_lat[Obs-1], transform=ccrs.PlateCarree(), s= 50, marker= "s", color= 'k', zorder=2)
    
    for iObs, Obs in enumerate(track.Obs):
    
        if iObs == 0:
            beering= gcc.bearing_at_p1( (Sm_lon[iObs], Sm_lat[iObs]), (Sm_lon[iObs +1], Sm_lat[iObs +1]) )
            prop_dist= gcc.distance_between_points( (Sm_lon[iObs], Sm_lat[iObs]), (Sm_lon[iObs +1], Sm_lat[iObs +1]))
            prop_speed = prop_dist/(t[iObs+1] - t[iObs] )/3600

            Stoll.loc[np.logical_and(Stoll['ID'] == ID, Stoll['Obs'] == Obs), 'prop_dir']= beering
            Stoll.loc[np.logical_and(Stoll['ID'] == ID, Stoll['Obs'] == Obs), 'prop_speed']= prop_speed

        elif iObs == len(track.Obs) -1:
            beering= gcc.bearing_at_p2( (Sm_lon[iObs-1], Sm_lat[iObs-1]), (Sm_lon[iObs], Sm_lat[iObs]) )
            prop_dist= gcc.distance_between_points( (Sm_lon[iObs-1], Sm_lat[iObs-1]), (Sm_lon[iObs], Sm_lat[iObs])) #in meters
            prop_speed = prop_dist/(t[iObs] - t[iObs-1] )/3600

            Stoll.loc[np.logical_and(Stoll['ID'] == ID, Stoll['Obs'] == Obs), 'prop_dir']= beering
            Stoll.loc[np.logical_and(Stoll['ID'] == ID, Stoll['Obs'] == Obs), 'prop_speed']= prop_speed


        else:
            beering= gcc.bearing_at_p2( (Sm_lon[iObs-1], Sm_lat[iObs-1]), (Sm_lon[iObs+1], Sm_lat[iObs+1]) )
            prop_dist= gcc.distance_between_points( (Sm_lon[iObs-1], Sm_lat[iObs-1]), (Sm_lon[iObs+1], Sm_lat[iObs+1])) #in meters
            prop_speed = prop_dist/(t[iObs+1] - t[iObs-1] )/3600

            Stoll.loc[np.logical_and(Stoll['ID'] == ID, Stoll['Obs'] == Obs), 'prop_dir']= beering
            Stoll.loc[np.logical_and(Stoll['ID'] == ID, Stoll['Obs'] == Obs), 'prop_speed']= prop_speed
         
        
 
Stoll= Stoll.set_index(['ID', 'time'])
Stoll= Stoll[ ['Obs', 'prop_dir', 'prop_speed'] ]
Stoll.to_csv(imp_dir + 'Stoll_list_dir-speed_smth'+string_smooth_param+'.csv')



print('Time: ', time.time() - start) 