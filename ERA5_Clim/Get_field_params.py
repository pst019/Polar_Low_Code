492/12#!/usr/bin/env python3
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

import time
start = time.perf_counter()

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
from f_carto import *
from f_useful import *
from f_meteo import *

from f_ERA5_clim import *



fignr= 7
Plot= True
# savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/Investigate_PL/'
# save= False

# 

version= "fram_run3"
ending= 'atl'

durlim= 6

start_year= 2004
end_year= 2019


# param= 'N_925-500'
param= 'theta_diff_925-500'

mean_rad= 300

time_start= False #in this case the first and last index of alltracks will be used
time_start= '2013-1-1' #the field data is needed for this time
time_end= '2013-2-1'


var, file_name= 't', 'levels_500'



"""get all tracks"""
csv_name= 'merged_tracks_'+version+'_'+ending+f'_{start_year}-{end_year}'

if durlim: csv_name += 'durlim'+str(durlim)
tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv')
tracks['time']= pd.to_datetime(tracks['time'])

tracks= tracks.set_index(['ID', 'step'])

tracks= tracks[(tracks['time'] > time_start) & (tracks['time'] < time_end)]



# indexes= remove_dublicate(tracks.index.get_level_values(0)



for i in [1]: #range(len(tracks)):
    # print(i)

    track= tracks.iloc[i]
        
  
    

    
    """plot map"""
    # var='var138'#'vo'
    # ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/vorticity_era5_{PL_time.year}_{PL_time.month:02}.nc")

    if param== 'theta_diff_925-500':
        ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/{file_name}_era5_{track.time.year}_{track.time.month:02}.nc")
        ds= ds.sel(time= track.time)
    
        theta_1, theta_0= PotTemp(ds.t.isel(plev=1), ds.plev[1]/100), PotTemp(ds.t.isel(plev=0), ds.plev[0]/100)

        ds[param]= (('lat', 'lon'), theta_1 - theta_0)
        ds[param].attrs['units'] = 'K'
    
    
    param_loc= ds[param].sel(lat= track.lat, lon= track.lon, method= 'nearest').values
    
    # print('At track location:', param_loc)

    dist= (110* np.sqrt( (ds.lat - track.lat)**2 + (np.cos(np.deg2rad(track.lat)) *(ds.lon- track.lon) )**2 ) )
    # dist= distance((ds.lat, ds.lon), (track.lat, track.lon))
    param_mean= float(ds[param].where(dist < mean_rad).mean().values)
    # print('Mean around track location:', param_mean)

    
    
    """map"""
    if Plot:
        fig = plt.figure(fignr, figsize= (10,8))
        plt.clf()
        ax= Plot_Polar_Stereo(fig, central_longitude= 20, extent= [-10, 40, 60, 77], subplot= (1,1,1), scalebar=False)
        # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-45, 45, 50, 80], subplot= (1,1,1))
        # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, 30, 80], subplot= (1,1,1))
        #ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-50, 80, 40, 80], subplot= (1,1,1))
        
        # vextr= np.max([np.max(ds[param]), -np.min(ds[param])])
        # cf= ax.contourf(ds.lon, ds.lat, ds[param], transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)
        # cf= ax.contourf(ds.lon, ds.lat, ds[param].where(dist < 500, transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)
    
    
        vmax, vmin= np.max(ds[param]), np.min(ds[param])
    
        cf= ax.contour(ds.lon, ds.lat, ds[param], transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r', vmin= vmin, vmax= vmax)
        cf= ax.contourf(ds.lon, ds.lat, ds[param].where(dist < 500) , transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r', vmin= vmin, vmax= vmax)
    
        cb= fig.colorbar(cf, ax= ax, shrink=0.7)
        varlabel= param 
        cb.set_label(varlabel , size=11)    
    
        ax.plot(track.lon, track.lat, marker='o', c= 'g',  transform= ccrs.PlateCarree() , label= track.name )
           
        plt.legend()
    



finish= time.perf_counter()

print(f'Finished in {round(finish - start,2)} seconds')