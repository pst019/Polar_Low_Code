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




# index= 32018070460 #track ID
# index= 12008010247
# index= 12008010186
# index= 32016070129
# index= 32016073107
# index= 32016072100
# index= 32016072651
# index= 120080100240
# index= 220081212620

# 120080100240,
#  120080100390,
#  120080100830,
#  120080101160,
#  120080101500,
index=  120080101800 #good
# index=  120080101960 # good
# index=  120080102070 #ok
# index=  120080102450 #also maybe rather an EC
# index= 120080102910 #rather EC - could be excluded by long lifetime
# index= 120080102970 #out of labrador sea
# index= 120080103000 #west of ireland in NWly flow, ok
# index= 120080103140 #coming out Labrador Sea, appear reasonable
# index= 120080103470 #North Sea, appears rather to be EC
# index=  120080103660 #in slipstream of EC, south of Iceland
# index=  120080103760 # likely transition - quite long, but ok
# index=  120080103880 #close to Ireland, ok
# index=  120080104220 #quite far out in the Atlantic, but why not
# index=  120080104270 #rather a bit of an EC
# index= 120080104690 #tip of Greenland - ok
# index= 120080104890 #in British Sea, but seems reasonable
# index= 120080104920 #ok
# index=  120080105010 #have to check this, SE of ireland
# index=  120080106060 #good
# index=  120080106070 #ok
# index= 120080106290 #ok
# index=  120080106981 #have to set different time due to split. Appear ok
# index=  120080107090 #Good
# index=  120080107200 #Denmark Straight - very short
# index= 120080107620 # close to France
# index=  120080107790

index= 120080100390


time_sel= 'orig_middle' #middle time step of the original PL for the PL list
# time_sel= 'orig_start'
# time_sel= 'orig_end'



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
    if time_sel== 'orig_middle':
        track_time= track['time'].iloc[track['Duration'].iloc[0]//2]

    
    
    """plot map"""
    if Plot_fields:
        # var='var138'#'vo'
        # ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/vorticity_era5_{track_time.year}_{track_time.month:02}.nc")
        if hem== 'NH':
            ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/vorticity_shift_era5_{track_time.year}_{track_time.month:02}.nc")
            # var= 'vo'
            var= 'var138'

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
            # var='msl'
            var= 'var151'  #'msl'
            
        if hem== 'SH':
            ds2= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/slp_SH_era5_{track_time.year}_{track_time.month:02}.nc")
            var= 'var151'  #'msl'
        ds2= ds2.sel(time= track_time)
        
        cs= ax.contour(ds2.lon, ds2.lat, ds2[var]/100, np.arange(900, 1100, 2), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
        #if numbers_cont:
        plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')
        
    
    """plot track"""
    ax.plot(track.lon.values, track.lat.values, c= 'g',  transform= ccrs.PlateCarree())
    ax.plot(track[track.time== track_time].lon.values, track[track.time== track_time].lat.values, marker='s', c= 'g',  transform= ccrs.PlateCarree())
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

print(track)