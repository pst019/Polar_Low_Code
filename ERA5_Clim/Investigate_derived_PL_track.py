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
import scipy
from scipy.ndimage import uniform_filter


fignr= 1
savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/Investigate_PL/'
save= True

# 
hem= 'NH'
version= ''

if hem == 'NH':
    version= "fram_run3_"
    ending= 'atl-pac'
    start_year= 1980 #2004
    end_year= 2019
    terrain_thresh= 50
    terrain_dist = 2 # it is actually 4
    more_info= True
    split, exceed_large_to_small= True, 40 
    time_start= '2008-1-1' #the vorticity and slp is needed in this time
    time_end= '2009-1-1'
    
    Plot_fields= True
if hem == 'SH':
    # version= "fram_run1_"
    ending= 'SH'
    start_year= 1988 #2004
    end_year= 2020

    Plot_fields= True
    time_start= False #in this case the first and last index of alltracks will be used
    time_start= '2016-7-1' #the vorticity and slp is needed in this time
    time_end= '2016-7-31'

durlim= 6



matchdist= 150

#tracks[(tracks.lat > 78) & (tracks.lon > 0) & (tracks.lon < 20)]
index= 120080106290
index= 120080110990
index= 120080116560
index= 120080212510
index= 120081205470

#tracks[(tracks.lat > 60) & (tracks.lon > -45) & (tracks.lon < -35)]
index= 120080100390
index= 120080107200
index= 120080202090

index= 120080205270
# index= 120080210030
# index= 120080210760 #this should really be excluded
# index= 120080211790
# index= 120080301330
# index= 120080304700
# index= 120080304870 #small movement
# index= 120080306300
# index= 120080306240
# index= 120080403960
# index= 120081009240
# index= 120081010190








# index= 120081213900


"""get all tracks"""
csv_name= 'PLs-from_merged_tracks_'+version+ending+f'_{start_year}-{end_year}'

if durlim: csv_name += 'durlim'+str(durlim)
# if more_info: csv_name += '_moreinfo'

if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
if split: csv_name += f'_split{exceed_large_to_small}'

if time_start != False: csv_name += f'_{time_start}-{time_end}'

tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/derived_PLs/"+csv_name+'.csv')
tracks['time']= pd.to_datetime(tracks['time'])



if time_start != False:
    tracks.time= pd.to_datetime(tracks.time)
    tracks= tracks[(tracks['time'] >= time_start) & (tracks['time'] < time_end)]

# tracks= tracks.set_index(['ID', 'step'])

track= tracks[tracks.ID== index]

"""alltracks - for a test"""
# tracks= pd.read_csv(Mediadir+"data/ERA5_Clim/track_lists/my_PC/tracks_"+test+".csv").set_index(['ID', 'step'])
# tracks['time']= pd.to_datetime(tracks['time'])



        
"""map"""
fig = plt.figure(fignr, figsize= (10,8))
fignr += 1
plt.clf()
if hem == 'NH':
    # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-50, 80, 40, 80], subplot= (1,1,1))
    # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, 30, 80], subplot= (1,1,1))
    ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, 30, 90], scalebar=False, subplot= (1,1,1), hem='NH', circular= True)

if hem == 'SH':
    ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, -80, -30], subplot= (1,1,1))   
    # ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, -90, -30], scalebar=False, subplot= (1,1,1), hem='SH', circular= True)
   


"""plot track"""
ax.plot(track.lon.values, track.lat.values, c= 'g',  transform= ccrs.PlateCarree())
ax.plot(track.lon.values, track.lat.values, 'x', c= 'g',  transform= ccrs.PlateCarree())
ax.plot(track.lon.values[0], track.lat.values[0], c= 'g', marker= 'o',  transform= ccrs.PlateCarree())



# if hem == 'NH':
#     lsm= xr.open_dataset(Mediadir + "data/ERA5_Clim/ERA5_data/lsm_shift_update.nc")

# if hem == 'SH':
#     lsm= xr.open_dataset(Mediadir + "data/ERA5_Clim/ERA5_data/lsm_SH.nc")
    
# lsm['lsm']= lsm['lsm'].where(lsm['lsm'] == 0, 1) #replaces all values that are not 0 with 1

# lsm= lsm.assign(extend = (('lat', 'lon'), uniform_filter(lsm['lsm'], size= (4, 4), mode= 'constant' ) ))
# lsm['extend']= lsm['extend'].where(lsm['extend'] == 0, 1) #replaces all values that are not 0 with 1

# cf= ax.contour(lsm.lon , lsm.lat, lsm.extend.values, transform= ccrs.PlateCarree(), colors= 'k', levels=[0.5]) #, vmin= -vextr, vmax= vextr)

# track['near_land']= [int(lsm['extend'].sel(lat= tracklat, lon= tracklon) ) for tracklat, tracklon in zip(track.lat, track.lon)]

# track['near_land'].groupby(['ID']).mean()

# print( track[['ID', 'near_land']].groupby(['ID']).mean() )

# print(track)
print(track[['lat', 'lon']])