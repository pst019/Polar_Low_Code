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
import matplotlib
from f_carto import *




fignr= 3

hem= 'SH'

# save=False
# savedir= homedir+ 'Polar_Low/STARS-ana/Figs/Track_eval/'


fig = plt.figure(fignr)
plt.clf()


if hem == 'NH':
    # ax= Plot_Polar_Stereo(fig, central_longitude= 20, extent= [-10, 40, 60, 77], subplot= (1,1,1), scalebar=False)
    #ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-45, 45, 50, 80], subplot= (1,1,1))
    ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, 30, 80], subplot= (1,1,1))
    #ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-50, 80, 40, 80], subplot= (1,1,1))
    

    var='lsm'#'vo'
    ds= xr.open_dataset(Mediadir + "data/ERA5_Clim/ERA5_data/lsm_shift_era5_2014_01.nc")
    ds= ds.isel(time= 0)
    #ds= ds.isel(plev= 0)

elif hem == 'SH':
    ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, -80, -30], subplot= (1,1,1))

    var='lsm'#'vo'
    ds= xr.open_dataset(Mediadir + "data/ERA5_Clim/ERA5_data/lsm_SH_era5_2014_01.nc")
    ds= ds.isel(time= 0)


ds[var]= ds[var].where(ds[var] == 0, 1) #replaces all values that are not 0 with 1


if hem == 'NH':
    
    """Bjørnøya"""
    min_lon, max_lon = 15, 20
    min_lat, max_lat = 72, 75
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    
    
    """Jan Mayen"""
    min_lon, max_lon = -10, 0
    min_lat, max_lat = 70, 75
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    
    """near lofoten"""
    min_lon, max_lon = 10, 13
    min_lat, max_lat = 67, 70
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    
    """between Iceland and Scotland"""
    min_lon, max_lon = -10, 0
    min_lat, max_lat = 59, 65
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    """SE Svalbard"""
    min_lon, max_lon = 20, 30
    min_lat, max_lat = 75, 77
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    """West Scottland"""
    min_lon, max_lon = -10, -8
    min_lat, max_lat = 56, 60
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    
    """West Spain"""
    min_lon, max_lon = -40, -15
    min_lat, max_lat = 30, 40
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    """Aleutan US"""
    min_lon, max_lon = -180, -165
    min_lat, max_lat = 50, 60
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    """Aleutan US north"""
    min_lon, max_lon = -180, -166.5
    min_lat, max_lat = 60, 62
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    """Chuckchi"""
    min_lon, max_lon = -177.25, -170
    min_lat, max_lat = 70, 75
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    """for shift Aleutan US"""
    min_lon, max_lon = -180+360, -165+360
    min_lat, max_lat = 50, 60
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    """for shift Aleutan US north"""
    min_lon, max_lon = -180+360, -166.5+360
    min_lat, max_lat = 60, 62
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    
    """for shift Chuckchi"""
    min_lon, max_lon = -177.25+360, -170+360
    min_lat, max_lat = 70, 75
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    
    """Aleutan Russia"""
    min_lon, max_lon = 165, 180
    min_lat, max_lat = 50, 57
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    
    """Japan Russia"""
    min_lon, max_lon = 145.5, 156.5
    min_lat, max_lat = 42, 51.25
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    """Sea of Japan"""
    min_lon, max_lon = 130, 135
    min_lat, max_lat = 36, 40
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    """SE of Japan"""
    min_lon, max_lon = 138, 141
    min_lat, max_lat = 32, 34.25
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    """Hudson Bay"""
    min_lon, max_lon = -81, -78.5
    min_lat, max_lat = 55.5, 63
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    """Hudson Bay"""
    min_lon, max_lon = -78.5, -77.25
    min_lat, max_lat = 56, 57.5
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
    """Lawrence Stream"""
    min_lon, max_lon = -61, -59
    min_lat, max_lat = 42, 44.5
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)


if hem == 'SH':
    """pac"""
    min_lon, max_lon = -180, -78
    min_lat, max_lat = -70, -30
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)


    """atl -ind"""
    min_lon, max_lon = -50, 110
    min_lat, max_lat = -65, -35
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)


    """Drake passage"""
    min_lon, max_lon = -65, -50
    min_lat, max_lat = -62.75, -53
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)


    """Drake passage"""
    min_lon, max_lon = -65, -60
    min_lat, max_lat = -63.25, -53
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)

    """South America"""
    min_lon, max_lon = -75, -65
    min_lat, max_lat = -63, -55.75
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)


    """NZ"""
    min_lon, max_lon = 157, 180
    min_lat, max_lat = -68, -47
    mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
    mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds[var]= ds[var].where(~(mask_lon & mask_lat), 0)
    
# cmap= 'Reds'
colors = ['blue'] + [(plt.cm.Reds(i)) for i in range(1,256)] 
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=256)

Lon, Lat= np.meshgrid(ds.lon, ds.lat)
cf= ax.pcolormesh(ds.lon , ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap) #, vmin= -vextr, vmax= vextr)
cb= plt.colorbar(cf, ax= ax, shrink=0.5, orientation= 'horizontal')


#cf= ax.pcolormesh(Lon , Lat , ds[var], transform= ccrs.PlateCarree(), cmap= cmap)

#    cf= ax.contourf(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), levels= 15, cmap= cmap, vmin= -vextr, vmax= vextr)

#    ax.contour(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), colors= 'red', levels=np.array([1.5, 2.0])*1E-4 )




if hem == 'NH':
    print(Mediadir + "data/ERA5_Clim/ERA5_data/lsm_shift_update.nc")    
    ds.to_netcdf(Mediadir + "data/ERA5_Clim/ERA5_data/lsm_shift_update.nc")

if hem == 'SH':
    print(Mediadir + "data/ERA5_Clim/ERA5_data/lsm_SH.nc")
    ds.to_netcdf(Mediadir + "data/ERA5_Clim/ERA5_data/lsm_SH.nc")