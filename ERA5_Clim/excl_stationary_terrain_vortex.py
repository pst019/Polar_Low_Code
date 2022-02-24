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
from matplotlib import colors
from scipy.ndimage import gaussian_filter, uniform_filter

# index = 2008011597 *10

fignr= 6
# savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/Investigate_PL/'
# save= True

# 
hem= 'SH'

if hem== 'NH':
    version= "" #"fram_run3_"
    ending= 'atl-pac'
    start_year= 1979
    end_year= 2020    

if hem== 'SH':
    version= "" #"fram_run1_"
    ending= 'SH180'
    start_year= 1979
    end_year= 2020    

durlim= 6




Plot= True
# Plot= False

write= True
# write = False

# threshold_vo_max_min= 25
# frac_large_to_small= 1.5

terrain_thresh= 50
terrain_dist=2

more_info= True
post_process= True


split, exceed_large_to_small= True, 40


time_start= False #in this case the first and last index of alltracks will be used
# time_start= '2008-1-1' #the vorticity and slp is needed in this time
# time_end= '2009-2-1'

# time_start= '1980-1-1' #the vorticity and slp is needed in this time
# time_end= '1980-6-1'



"""get all tracks"""
csv_name= 'merged_tracks_'+version+ending+f'_{start_year}-{end_year}'

if durlim: csv_name += 'durlim'+str(durlim)
if more_info: csv_name += '_moreinfo'
if post_process: csv_name += '_post-process'
elif split: csv_name += f'_split{exceed_large_to_small}'

tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks2/"+csv_name+'.csv')
tracks['time']= pd.to_datetime(tracks['time'])

tracks.loc[tracks.lon >= 180, 'lon'] -= 360
# tracks.ID *= 10
# tracks.vo *= 1E2

# tracks= tracks.set_index(['ID', 'step'])

if time_start != False:
    tracks= tracks[(tracks['time'] >= time_start) & (tracks['time'] < time_end)]


# total_splits= 0
# track_ID_list= remove_dublicate(tracks['ID'])


# fig = plt.figure(fignr)#, figsize= (13,10) )
# fignr+=1
# plt.clf()


# tracks= tracks[tracks.time.dt.year == 2008]

# tlat, tlon= 64.75, -39.75
# tracks[(tracks.lat== tlat) & (tracks.lon== tlon)]


#idea: extend landmask and remove all systems that are within the landmask for their whole lifetime


fig = plt.figure(fignr, figsize= (13,10) )
fignr+=1
plt.clf()

# ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-100, 80, 30, 80], subplot= (1,1,1))

"""lsm"""
if hem == 'NH':
    ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, 30, 80], subplot= (1,1,1))
    lsm= xr.open_dataset(Mediadir + "data/ERA5_Clim/ERA5_data/lsm_shift_update.nc")
    
    # lsm= xr.open_dataset(Mediadir + "data/ERA5_Clim/ERA5_data/lsm_shift_era5_2014_01.nc")
    # lsm= lsm.isel(time= 0)
if hem == 'SH':
    ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, -80, -30], subplot= (1,1,1))
    # lsm= xr.open_dataset(Mediadir + "data/ERA5_Clim/ERA5_data/lsm_SH.nc")

    lsm= xr.open_dataset(Mediadir + "data/ERA5_Clim/ERA5_data/lsm_SH_era5_2014_01.nc")
    lsm= lsm.isel(time= 0)
    
lsm['lsm']= lsm['lsm'].where(lsm['lsm'] == 0, 1) #replaces all values that are not 0 with 1

lsm= lsm.assign(extend = (('lat', 'lon'), uniform_filter(lsm['lsm'], size= (terrain_dist, terrain_dist), mode= 'constant' ) ))
lsm['extend']= lsm['extend'].where(lsm['extend'] == 0, 1) #replaces all values that are not 0 with 1

# cmap = colors.ListedColormap(['white', 'blue']) # ['blue'] + [(plt.cm.Reds(i)) for i in range(1,256)] 
# cmap = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=256)
# cf= ax.pcolormesh(lsm.lon , lsm.lat, lsm.lsm.values, transform= ccrs.PlateCarree(), cmap= cmap) #, vmin= -vextr, vmax= vextr)
# cb= plt.colorbar(cf, ax= ax, shrink=0.5, orientation= 'horizontal')

cf= ax.contour(lsm.lon , lsm.lat, lsm.extend.values, transform= ccrs.PlateCarree(), colors= 'k', levels=[0.5]) #, vmin= -vextr, vmax= vextr)


# tracks['near_land']= np.nan

# tracks= tracks.set_index('ID')

# for track_ID in remove_dublicate(tracks.index):
#     print(track_ID)
#     track= tracks.loc[track_ID] #[tracks.ID== track_ID]
#     tracks.loc[track.index[0], 'near_land']= [int(lsm['extend'].sel(lat= tracklat, lon= tracklon) ) for tracklat, tracklon in zip(track.lat, track.lon)]
    


tracks['near_land']= [int(lsm['extend'].sel(lat= tracklat, lon= tracklon) ) for tracklat, tracklon in zip(tracks.lat, tracks.lon)]

tracks_ind= tracks[['ID', 'near_land']].groupby(['ID']).mean() 

tracks= tracks.set_index('ID').loc[tracks_ind[tracks_ind['near_land'] <= terrain_thresh/100].index].reset_index()

"""simple density map"""

loc_freq= tracks.groupby(['lat', 'lon']).size()
loc_freq= loc_freq.to_xarray().fillna(0)

# loc_freq.values= gaussian_filter(loc_freq, sigma=(3, 3))

cf= ax.pcolormesh(loc_freq.lon, loc_freq.lat, loc_freq.values, transform= ccrs.PlateCarree(), alpha= 0.5, cmap= 'Reds' ) #, levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)
cb= fig.colorbar(cf, ax= ax, shrink=0.7, orientation= 'horizontal')

figname= homedir + f'Polar_Low/ERA5_PLclim/Figs/excl_terrain/Map_frequency_terrain{terrain_dist}-excl{terrain_thresh}'
print(figname)
fig.savefig(figname, bbox_inches = 'tight')

if write:
    outfile= Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+f'_terrain{terrain_dist}-excl{terrain_thresh}.csv'
    print('write:', outfile)
    tracks.set_index('ID').to_csv(outfile)        
     
