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
import matplotlib.colors as colors
from scipy import ndimage


import numpy as np
from f_useful import *
# from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs

#
save= True
save= False
savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/Investigate_situation/'

fignr= 4

#track_dir= Mediadir+"/data/ERA5_Clim/track_lists/fram/"

hem= 'SH'

if hem == 'NH':
    version= ""
    ending= 'atl-pac'


if hem == 'SH':
    version= "" #"fram_run1_"
    ending= 'SH180'



start_year= 1979
end_year= 2020
more_info= True
post_process= True
terrain_thresh= 50
terrain_dist = 2

durlim= 6




#matchdist= 150


"""import track list"""
csv_name= 'merged_tracks_'+version+ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
if more_info: csv_name += '_moreinfo'

if post_process: csv_name += '_post-process'

if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'

tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks2/"+csv_name+'.csv')







"""make frequency distribution"""
loc_freq= tracks.groupby(['lat', 'lon']).size()
loc_freq= loc_freq.to_xarray().fillna(0)

loc_freq /= (end_year -start_year +1) #per year


"""simple density map per pixel"""
fig = plt.figure(fignr, figsize= (13,10) )
fignr+=1
plt.clf()

if hem == 'NH':
    # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-50, 80, 40, 80], subplot= (1,1,1))
    # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, 30, 80], subplot= (1,1,1))
    ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, 30, 90], scalebar=False, subplot= (1,1,1), hem='NH', circular= True)
    
if hem == 'SH':
    # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, -80, -30], subplot= (1,1,1))   
    ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, -90, -30], scalebar=False, subplot= (1,1,1), hem='SH', circular= True)


# cf= ax.pcolormesh(loc_freq.lon, loc_freq.lat, loc_freq.values, transform= ccrs.PlateCarree(), cmap='Reds' ) #, levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)
#logarithmic
cf= ax.pcolormesh(loc_freq.lon, loc_freq.lat, loc_freq.values, transform= ccrs.PlateCarree(), norm=colors.LogNorm(vmin=1, vmax=loc_freq.max()), cmap='Reds' ) #, levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)

cb= fig.colorbar(cf, ax= ax, shrink=0.7)#, orientation= 'horizontal')
cb.set_label('Track points per grid cell per year')




"""simple density map per box, noramized"""
per_dist = 500 #km
normalize= True


if per_dist:
    loc_freq.values= ndimage.uniform_filter(loc_freq, size=(1+int(per_dist/(110/4)), 1+ int(1/np.cos(np.deg2rad(50)) *per_dist/(110/4))) )
else:
    loc_freq.values= ndimage.uniform_filter(loc_freq, size=(6, 10))

if normalize: loc_freq /=np.cos(np.deg2rad(loc_freq.lat))
if per_dist: loc_freq *= (per_dist/ (0.25*110))**2


fig = plt.figure(fignr, figsize= (13,10) )
fignr+=1
plt.clf()

if hem == 'NH':
    ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, 30, 90], scalebar=False, subplot= (1,1,1), hem='NH', circular= True)   
elif hem == 'SH':
    ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, -90, -30], scalebar=False, subplot= (1,1,1), hem='SH', circular= True)


# cf= ax.pcolormesh(ds.lon - 0.25, ds.lat + 0.125, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)

# cf= ax.pcolormesh(loc_freq.lon, loc_freq.lat, loc_freq.values, transform= ccrs.PlateCarree(), cmap='Reds')#, vmax= 4000 ) #, levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)
# logarithmic
cf= ax.pcolormesh(loc_freq.lon, loc_freq.lat, loc_freq.values, transform= ccrs.PlateCarree(), norm=colors.LogNorm(vmin=50, vmax=loc_freq.max()), cmap='Reds' ) #, levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)

cb= fig.colorbar(cf, ax= ax, shrink=0.7)#, orientation= 'horizontal')

if per_dist:
    cb.set_label(f'Annual track points within {per_dist}km distance')
else:
    cb.set_label('Track points per grid cell')



if save:
    savefile= savedir+ f'{hem}_track_density.png'
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight', dpi= 100)    


"""plot tracks passing 180E"""
fig = plt.figure(fignr, figsize= (13,10) )
fignr+=1
plt.clf()

if hem == 'NH':
    ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, 30, 90], scalebar=False, subplot= (1,1,1), hem='NH', circular= True)   
elif hem == 'SH':
    # ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, -90, -30], scalebar=False, subplot= (1,1,1), hem='SH', circular= True)
    ax= Plot_PlateCarree(fig, central_longitude= 180, extent= [130, 230, -70, -30], subplot= (1,1,1))   



track_ID_crossing_180= remove_dublicate(tracks[tracks.lon== 180].ID)

for track_ID in track_ID_crossing_180[:10]:
    track= tracks[tracks.ID== track_ID]
    ax.plot(track.lon.values, track.lat.values,  transform= ccrs.PlateCarree())
    # ax.plot(track[track.time== track_time].lon.values, track[track.time== track_time].lat.values, marker='s', c= 'g',  transform= ccrs.PlateCarree())
    # ax.plot(track.lon.values[0], track.lat.values[0], c= 'g', marker= 'o',  transform= ccrs.PlateCarree())
  


