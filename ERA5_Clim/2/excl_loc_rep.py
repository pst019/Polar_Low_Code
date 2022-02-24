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


write= True

hem= 'NH'
version= ''

time_start= False #in this case the first and last index of alltracks will be used

excl_double_ts= True

if hem == 'NH':
    # version= "fram_run3_"
    ending= 'atl-pac'
    start_year= 1979 #2004
    end_year= 2020
    terrain_thresh= False #50
    terrain_dist = 2 # it is actually 4
    more_info= True
    split, exceed_large_to_small= True, 40 
    # time_start= '1980-1-1' #the vorticity and slp is needed in this time
    # time_end= '1981-1-1'
    
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

if excl_double_ts: csv_name= '_rem-2ts'

if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
if split: csv_name += f'_split{exceed_large_to_small}'

if time_start != False: csv_name += f'_{time_start}-{time_end}'

tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/"+csv_dir+csv_name+'.csv')
tracks['time']= pd.to_datetime(tracks['time'])




if time_start != False:
    tracks.time= pd.to_datetime(tracks.time)
    tracks= tracks[(tracks['time'] >= time_start) & (tracks['time'] < time_end)]

tracks= tracks.drop(columns= 'Duration')

"""remove location repeaters"""
track_ID_list= remove_dublicate(tracks['ID'])

red_tracks = pd.DataFrame(columns=tracks.columns) #make an empty dataframe to append the tracks
red_tracks = red_tracks.set_index(['ID', 'step'])

tracks= tracks.set_index(['ID', 'step'])

#index= 120080210760 for testing

for index in track_ID_list:
    track=tracks.loc[index].copy()

    # track= track.set_index('step')
    track['same_location']= 0 #same location as previous time step
    for s in range(1, len(track)):
        if (track.loc[s, 'lat'] == track.loc[ s-1, 'lat']) & (track.loc[s, 'lon'] == track.loc[ s-1, 'lon']): #checks for similar latitude and longitude
            track.loc[s, 'same_location'] = track.loc[s-1, 'same_location'] +1
        
    end_repeater = track.loc[len(track)-1, 'same_location']
    if end_repeater > 3:
        track= track[:-end_repeater]
        # tracks= tracks.drop( (index, :-end_repeater) )


    track= track.reset_index()
    begin_repeater= np.max(np.where(track['same_location']== track['step'])[0])
    if begin_repeater > 3: track= track[begin_repeater:]
    
    track= track.drop(columns= 'same_location')
    track['ID']= index
    track= track.set_index(['ID', 'step'])
    
    red_tracks= pd.concat([red_tracks, track])


"""write the file"""
if write:
    outfile= Mediadir+"/data/ERA5_Clim/track_lists/"+csv_dir+csv_name+'_rem-locrep.csv'
    print('write:', outfile)
    red_tracks.to_csv(outfile) 






    
"""make frequency distribution"""
loc_freq= red_tracks.groupby(['lat', 'lon']).size()
loc_freq= loc_freq.to_xarray().fillna(0)

loc_freq /= (end_year -start_year +1) #per year


"""simple density map per pixel"""
fig = plt.figure(fignr, figsize= (13,10) )
fignr+=1
plt.clf()

if hem == 'NH':
    ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, 30, 90], scalebar=False, subplot= (1,1,1), hem='NH', circular= True)
elif hem == 'SH':
    ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, -90, -30], scalebar=False, subplot= (1,1,1), hem='SH', circular= True)


cf= ax.pcolormesh(loc_freq.lon, loc_freq.lat, loc_freq.values, transform= ccrs.PlateCarree(), cmap='Reds' ) #, levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)
cb= fig.colorbar(cf, ax= ax, shrink=0.7)#, orientation= 'horizontal')
cb.set_label('Track points per grid cell per year')

plt.tight_layout()



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
# logarithmic
# cf= ax.pcolormesh(loc_freq.lon, loc_freq.lat, loc_freq.values, transform= ccrs.PlateCarree(), norm=colors.LogNorm(vmin=loc_freq.min(), vmax=loc_freq.max()) ) #, levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)

cf= ax.pcolormesh(loc_freq.lon, loc_freq.lat, loc_freq.values, transform= ccrs.PlateCarree(), cmap='Reds' ) #, levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)
cb= fig.colorbar(cf, ax= ax, shrink=0.7)#, orientation= 'horizontal')

if per_dist:
    cb.set_label(f'Annual track points within {per_dist}km distance')
else:
    cb.set_label('Track points per grid cell')

plt.tight_layout()



# print(track)
# print(track[['lat', 'lon']])