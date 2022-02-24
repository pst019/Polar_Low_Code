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
savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/Investigate_situation/'
save= True

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




"""get all tracks"""
more_info= True #True

terrain_thresh= 50
terrain_dist=4

#matchdist= 150


"""import track list"""
csv_name= 'merged_tracks_'+version+'_'+ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
if more_info: csv_name += '_moreinfo'

if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'

tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv')

tracks['time']= pd.to_datetime(tracks['time'])

tracks= tracks.set_index(['ID', 'step'])

# track= tracks[tracks.ID== index]

"""alltracks - for a test"""
# tracks= pd.read_csv(Mediadir+"data/ERA5_Clim/track_lists/my_PC/tracks_"+test+".csv").set_index(['ID', 'step'])
# tracks['time']= pd.to_datetime(tracks['time'])

# time_now= '2008-7-7 00:00:00'
# pd.DatetimeIndex('2008-7-7 00:00:00')

time_now= np.datetime64('2016-07-05T00:00')


"""map"""
fig = plt.figure(fignr, figsize= (10,10))
fignr += 1
plt.clf()
if hem == 'NH':
    # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-50, 80, 40, 80], subplot= (1,1,1))
    ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, 30, 80], subplot= (1,1,1))
if hem == 'SH':
    # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, -80, -30], subplot= (1,1,1))   
    ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, -90, -30], scalebar=False, subplot= (1,1,1), hem='SH', circular= True)




"""plot map"""
if Plot_fields:
    # var='var138'#'vo'
    # ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/vorticity_era5_{track_time.year}_{track_time.month:02}.nc")
    if hem== 'NH':
        ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/vorticity_shift_era5_{time_now.astype(object).year}_{time_now.astype(object).month:02}.nc")
        var= 'vo'
    if hem== 'SH':
        ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/vorticity_SHm1_era5_{time_now.astype(object).year}_{time_now.astype(object).month:02}.nc")
        var= 'var138'
    ds[var]*= 1E5
    
    ds= ds.sel(time= time_now)
    ds= ds.isel(plev= 0)
    
    # vextr= np.max([np.max(ds[var]), -np.min(ds[var])])*.7
    #    # cf= ax.pcolormesh(ds.lon - 0.25, ds.lat + 0.125, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
    # cf= ax.contourf(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)
    # cb= fig.colorbar(cf, ax= ax, shrink=0.7)
    # cb= fig.colorbar(cf, ax= ax, shrink=0.5, orientation='horizontal') #, extend='both')
    # varlabel= r'Relative vorticity [10$^{-5}$ s$^{-1}$]'
    #
    ax.contour(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), colors= 'red', levels=np.array([15, 20]) )

    #if plevel != None: varlabel = var+#'_'+str(plevel)
    # cb.set_label(varlabel , size=11)    
    # cb.set_label(ds[var].long_name + ' ['+ ds[var].units + ']', size=11)    
    
    # var= 'var151'  #'msl'
    # ds2= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/slp_era5_{track_time.year}_{track_time.month:02}.nc")
    
    if hem== 'NH':
        ds2= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/slp_shift_era5_{time_now.astype(object).year}_{time_now.astype(object).month:02}.nc")
        var='msl'
    if hem== 'SH':
        ds2= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/slp_SH_era5_{time_now.astype(object).year}_{time_now.astype(object).month:02}.nc")
        var= 'var151'  #'msl'

    ds2= ds2.sel(time= time_now)
    cs= ax.contour(ds2.lon, ds2.lat, ds2[var]/100, np.arange(900, 1100, 10), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
    #if numbers_cont:
    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')
    

tracks_now= tracks[tracks['time']== time_now]

ax.scatter(tracks_now.lon.values, tracks_now.lat.values, marker= 's', s= 50, c= 'b',  transform= ccrs.PlateCarree(), zorder= 2)

for ti,track_ID in enumerate(tracks_now.index.get_level_values('ID')):
    # ax.plot(tracks.loc[track_ID].lon.values, tracks.loc[track_ID].lat.values, c='g',  transform= ccrs.PlateCarree())
    ax.scatter(tracks.loc[track_ID].lon.values, tracks.loc[track_ID].lat.values, marker='x', s= 3,  c='b',  transform= ccrs.PlateCarree())
    # ax.plot(tracks.loc[track_ID].lon.values, tracks.loc[track_ID].lat.values, 'x',  c='g',  transform= ccrs.PlateCarree())
    if ti== 0: ax.scatter(tracks.loc[track_ID].lon.values[0], tracks.loc[track_ID].lat.values[0], marker= 'o', s=30, c='b',  transform= ccrs.PlateCarree(), label= 'all tracks')
    else: ax.scatter(tracks.loc[track_ID].lon.values[0], tracks.loc[track_ID].lat.values[0], marker= 'o', s= 30, c='b',  transform= ccrs.PlateCarree())


tracks_strong= tracks_now[tracks_now['slp'] < 970]

ax.scatter(tracks_strong.lon.values, tracks_strong.lat.values, marker= 's', s= 50, c= 'g',  transform= ccrs.PlateCarree(), zorder= 2)

for ti,track_ID in enumerate(tracks_strong.index.get_level_values('ID')):
    # ax.plot(tracks.loc[track_ID].lon.values, tracks.loc[track_ID].lat.values, c='g',  transform= ccrs.PlateCarree())
    ax.scatter(tracks.loc[track_ID].lon.values, tracks.loc[track_ID].lat.values, marker='x', s= 3,  c='g',  transform= ccrs.PlateCarree())
    # ax.plot(tracks.loc[track_ID].lon.values, tracks.loc[track_ID].lat.values, 'x',  c='g',  transform= ccrs.PlateCarree())
    if ti== 0: ax.scatter(tracks.loc[track_ID].lon.values[0], tracks.loc[track_ID].lat.values[0], marker= 'o', s=30, c='g',  transform= ccrs.PlateCarree(), label= 'tracks, slp < 970hPa')
    else: ax.scatter(tracks.loc[track_ID].lon.values[0], tracks.loc[track_ID].lat.values[0], marker= 'o', s= 30, c='g',  transform= ccrs.PlateCarree())



tracks_strong= tracks_strong[tracks_strong['vo'] > 0.3]

ax.scatter(tracks_strong.lon.values, tracks_strong.lat.values, marker= 's', s= 50, c= 'purple',  transform= ccrs.PlateCarree(), zorder= 2)

for ti,track_ID in enumerate(tracks_strong.index.get_level_values('ID')):
    # ax.plot(tracks.loc[track_ID].lon.values, tracks.loc[track_ID].lat.values, c='purple',  transform= ccrs.PlateCarree())
    ax.scatter(tracks.loc[track_ID].lon.values, tracks.loc[track_ID].lat.values, marker='x', s= 3,  c='purple',  transform= ccrs.PlateCarree())
    # ax.plot(tracks.loc[track_ID].lon.values, tracks.loc[track_ID].lat.values, 'x',  c='purple',  transform= ccrs.PlateCarree())
    if ti== 0: ax.scatter(tracks.loc[track_ID].lon.values[0], tracks.loc[track_ID].lat.values[0], marker= 'o', s=30, c='purple',  transform= ccrs.PlateCarree(), label= 'slp < 970hPa & vo > .3')
    else: ax.scatter(tracks.loc[track_ID].lon.values[0], tracks.loc[track_ID].lat.values[0], marker= 'o', s= 30, c='purple',  transform= ccrs.PlateCarree())


plt.legend()

# tracks_strong= tracks_now[tracks_now['slp'] < 980]

# ax.scatter(tracks_strong.lon.values, tracks_strong.lat.values, marker= 's', s= 50, c= 'r',  transform= ccrs.PlateCarree(), zorder= 2)

# for ti,track_ID in enumerate(tracks_strong.index.get_level_values('ID')):
#     # ax.plot(tracks.loc[track_ID].lon.values, tracks.loc[track_ID].lat.values, c='g',  transform= ccrs.PlateCarree())
#     ax.scatter(tracks.loc[track_ID].lon.values, tracks.loc[track_ID].lat.values, marker='x', s= 3,  c='r',  transform= ccrs.PlateCarree())
#     # ax.plot(tracks.loc[track_ID].lon.values, tracks.loc[track_ID].lat.values, 'x',  c='g',  transform= ccrs.PlateCarree())
#     if ti== 0: ax.scatter(tracks.loc[track_ID].lon.values[0], tracks.loc[track_ID].lat.values[0], marker= 'o', s=30, c='r',  transform= ccrs.PlateCarree(), label= 'strong tracks')
#     else: ax.scatter(tracks.loc[track_ID].lon.values[0], tracks.loc[track_ID].lat.values[0], marker= 'o', s= 30, c='r',  transform= ccrs.PlateCarree())





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
plt.tight_layout()

if save:
    savefile= savedir+ f'{hem}_{time_now}.png'
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight', dpi= 100)    