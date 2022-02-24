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
save= False

# 
hem= 'SH'

if hem == 'NH':
    version= "fram_run3"
    ending= 'atl-pac'
    start_year= 1980 #2004
    end_year= 2019
if hem == 'SH':
    version= "fram_run1"
    ending= 'SH'
    start_year= 1988 #2004
    end_year= 2020

    Plot_fields= True

durlim= 6



matchdist= 150
split= True
split= False



# index= 32018070460 #track ID
# index= 12008010247
# index= 12008010186
# index= 32016070129
# index= 32016073107
# index= 32016072100
index= 32016072651

# time_sel= 'orig_middle' #middle time step of the original PL for the PL list
time_sel= 'orig_start'
# time_sel= 'orig_end'



"""get all tracks"""
csv_name= 'merged_tracks_'+version+'_'+ending+f'_{start_year}-{end_year}'

if durlim: csv_name += 'durlim'+str(durlim)
if split: csv_name += '_split'


tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv')
tracks['time']= pd.to_datetime(tracks['time'])


time_start= False #in this case the first and last index of alltracks will be used
time_start= '2016-7-1' #the vorticity and slp is needed in this time
time_end= '2016-7-31'

if time_start != False:
    tracks.time= pd.to_datetime(tracks.time)
    tracks= tracks[(tracks['time'] >= time_start) & (tracks['time'] < time_end)]

# tracks= tracks.set_index(['ID', 'step'])

track= tracks[tracks.ID== index]

"""alltracks - for a test"""
# tracks= pd.read_csv(Mediadir+"data/ERA5_Clim/track_lists/my_PC/tracks_"+test+".csv").set_index(['ID', 'step'])
# tracks['time']= pd.to_datetime(tracks['time'])



for time_sel in [time_sel]: # ['orig_start']: #[ 'orig_middle' ]: #, 'orig_start', 'orig_end', 'matched_start', 'matched_end', 'matched_middle', 'matched_2_3']:
    # for index in PL_indexes:
        # print(index)
        
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

    "PL original"

    if time_sel== 'orig_start':
        track_time= track['time'].iloc[0]
    if time_sel== 'orig_end':
        track_time= track['time'].iloc[-1]

    

    
    """plot map"""
    if Plot_fields:
        # var='var138'#'vo'
        # ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/vorticity_era5_{track_time.year}_{track_time.month:02}.nc")
        if hem== 'NH':
            ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/vorticity_shift_era5_{track_time.year}_{track_time.month:02}.nc")
            var= 'vo'
        if hem== 'SH':
            ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/vorticity_SHm1_era5_{track_time.year}_{track_time.month:02}.nc")
            var= 'var138'

        ds[var]*= 1E5
        
        ds= ds.sel(time= track_time)
        ds= ds.isel(plev= 0)
        
        vextr= np.max([np.max(ds[var]), -np.min(ds[var])])*.7
        #    # cf= ax.pcolormesh(ds.lon - 0.25, ds.lat + 0.125, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
        cf= ax.contourf(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)
        #
        ax.contour(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), colors= 'red', levels=np.array([15, 20]) )
        cb= fig.colorbar(cf, ax= ax, shrink=0.7)
        # cb= fig.colorbar(cf, ax= ax, shrink=0.5, orientation='horizontal') #, extend='both')
        varlabel= r'Relative vorticity [10$^{-5}$ s$^{-1}$]'
        #if plevel != None: varlabel = var+#'_'+str(plevel)
        cb.set_label(varlabel , size=11)    
        # cb.set_label(ds[var].long_name + ' ['+ ds[var].units + ']', size=11)    
        
        # var= 'var151'  #'msl'
        # ds2= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/slp_era5_{track_time.year}_{track_time.month:02}.nc")
        # var='msl'
        # ds2= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/slp_shift_era5_{track_time.year}_{track_time.month:02}.nc")
        if hem== 'NH':
            ds2= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/slp_shift_era5_{track_time.year}_{track_time.month:02}.nc")
            var='msl'
        if hem== 'SH':
            ds2= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/slp_SH_era5_{track_time.year}_{track_time.month:02}.nc")
            var= 'var151'  #'msl'
        ds2= ds2.sel(time= track_time)
        
        cs= ax.contour(ds2.lon, ds2.lat, ds2[var]/100, np.arange(900, 1100, 5), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
        #if numbers_cont:
        plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')
        
    
    """plot track"""
    ax.plot(track.lon.values, track.lat.values, c= 'g',  transform= ccrs.PlateCarree())
    ax.plot(track[track.time== track_time].lon.values, track[track.time== track_time].lat.values, marker='s', c= 'g',  transform= ccrs.PlateCarree())
    ax.plot(track.lon.values[0], track.lat.values[0], c= 'g', marker= 'o',  transform= ccrs.PlateCarree())
    


if hem == 'NH':
    lsm= xr.open_dataset(Mediadir + "data/ERA5_Clim/ERA5_data/lsm_shift_update.nc")

if hem == 'SH':
    lsm= xr.open_dataset(Mediadir + "data/ERA5_Clim/ERA5_data/lsm_SH.nc")
    
lsm['lsm']= lsm['lsm'].where(lsm['lsm'] == 0, 1) #replaces all values that are not 0 with 1

lsm= lsm.assign(extend = (('lat', 'lon'), uniform_filter(lsm['lsm'], size= (4, 4), mode= 'constant' ) ))
lsm['extend']= lsm['extend'].where(lsm['extend'] == 0, 1) #replaces all values that are not 0 with 1

cf= ax.contour(lsm.lon , lsm.lat, lsm.extend.values, transform= ccrs.PlateCarree(), colors= 'k', levels=[0.5]) #, vmin= -vextr, vmax= vextr)



track['near_land']= [int(lsm['extend'].sel(lat= tracklat, lon= tracklon) ) for tracklat, tracklon in zip(track.lat, track.lon)]

# track['near_land'].groupby(['ID']).mean()

print( track[['ID', 'near_land']].groupby(['ID']).mean() )

