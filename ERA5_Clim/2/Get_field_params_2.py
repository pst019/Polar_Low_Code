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



fignr= 6
Plot= True
write= False

# Plot= False
# write= True

# savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/Investigate_PL/'
# save= False

# 

version= "fram_run3"
ending= 'atl'

durlim= 6

start_year= 2004
end_year= 2019


# param= 'N_925-500'
param_var, file_name= 'theta_diff_925-500', 'levels_500'
param_var, file_name= 'theta_trop', 'PV'
param_var, file_name= 'U_trop', 'PV'

param_rad= 250
param_type= 'mean'
# param_type= 'max'

param= param_var+'_'+param_type+ str(param_rad)

time_start= False #in this case the first and last index of alltracks will be used
time_start= '2013-1-1' #the field data is needed for this time
time_end= '2014-1-1'





"""get all tracks"""
csv_name= 'merged_tracks_'+version+'_'+ending+f'_{start_year}-{end_year}'

if durlim: csv_name += 'durlim'+str(durlim)
tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv')
tracks['time']= pd.to_datetime(tracks['time'])

tracks= tracks.set_index(['ID', 'step'])

tracks= tracks[(tracks['time'] >= time_start) & (tracks['time'] < time_end)]


tracks[param]= np.nan


# indexes= remove_dublicate(tracks.index.get_level_values(0)

time_steps = remove_dublicate(tracks.time)
if Plot: time_steps= time_steps[:1]

for time_step in time_steps:
    if time_step.day == 1 and time_step.hour == 0:
        print(time_step)
    
    """import the parameter"""
    if param_var== 'theta_diff_925-500':
        ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/{file_name}_era5_{time_step.year}_{time_step.month:02}.nc")
        ds= ds.sel(time= time_step)
    
        theta_1, theta_0= PotTemp(ds.t.isel(plev=1), ds.plev[1]/100), PotTemp(ds.t.isel(plev=0), ds.plev[0]/100)

        ds[param_var]= (('lat', 'lon'), theta_1 - theta_0)
        ds[param_var].attrs['units'] = 'K'


    if param_var== 'shear_925-500':
        #get the strength and the angle
        ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/{file_name}_era5_{time_step.year}_{time_step.month:02}.nc")
        ds= ds.sel(time= time_step)
    
        # theta_1, theta_0= PotTemp(ds.t.isel(plev=1), ds.plev[1]/100), PotTemp(ds.t.isel(plev=0), ds.plev[0]/100)

        # ds[param_var]= (('lat', 'lon'), theta_1 - theta_0)
        # ds[param_var].attrs['units'] = 'K'    
    
    if param_var== 'theta_trop':
        ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/{file_name}_era5_{time_step.year}_{time_step.month:02}.nc")
        ds= ds.sel(time= time_step)   
        ds= ds.isel(lev= 0)
        
        ds[param_var]= ds.pt

    if param_var== 'U_trop':
        ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/{file_name}_era5_{time_step.year}_{time_step.month:02}.nc")
        ds= ds.sel(time= time_step)   
        ds= ds.isel(lev= 0)
        
        ds[param_var]= np.sqrt(ds.u**2 +ds.v**2)

    
    tracks_now= tracks[tracks.time== time_step]
    
    if Plot: tracks_now= tracks_now[:1]
    
    for i in range(len(tracks_now)): 
        # print(i)
        track= tracks_now.iloc[i]
           
    
        # param_loc= ds[param].sel(lat= track.lat, lon= track.lon, method= 'nearest').values
        
        # print('At track location:', param_loc)
    
        dist= (110* np.sqrt( (ds.lat - track.lat)**2 + (np.cos(np.deg2rad(track.lat)) *(ds.lon- track.lon) )**2 ) )
        # dist= distance((ds.lat, ds.lon), (track.lat, track.lon))
        
        if param_type== 'mean':
            param_mean= float(ds[param_var].where(dist < param_rad).mean().values)
        elif param_type== 'max':
            param_mean= float(ds[param_var].where(dist < param_rad).max().values)            
        # print('Mean around track location:', param_mean)
    
        # tracks_now[param]= param_mean
        tracks.loc[track.name, param]= param_mean    
    
    
        
        """map"""
        if Plot:
            fig = plt.figure(fignr, figsize= (10,8))
            plt.clf()
            # ax= Plot_Polar_Stereo(fig, central_longitude= 20, extent= [-10, 40, 60, 77], subplot= (1,1,1), scalebar=False)
            # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-45, 45, 50, 80], subplot= (1,1,1))
            # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, 30, 80], subplot= (1,1,1))
            ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-50, 80, 30, 80], subplot= (1,1,1))
            
            # vextr= np.max([np.max(ds[param]), -np.min(ds[param])])
            # cf= ax.contourf(ds.lon, ds.lat, ds[param], transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)
            # cf= ax.contourf(ds.lon, ds.lat, ds[param].where(dist < 500, transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)
        
        
            vmax, vmin= np.max(ds[param_var]), np.min(ds[param_var])
        
            cf= ax.contour(ds.lon, ds.lat, ds[param_var], transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r', vmin= vmin, vmax= vmax)
            cf= ax.contourf(ds.lon, ds.lat, ds[param_var].where(dist < 500) , transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r', vmin= vmin, vmax= vmax)
        
            cb= fig.colorbar(cf, ax= ax, shrink=0.7)
            varlabel= param_var
            cb.set_label(varlabel , size=11)    
        
            ax.plot(track.lon, track.lat, marker='o', c= 'g',  transform= ccrs.PlateCarree() , label= track.name )
               
            plt.legend()
   
 


print(param)
    
if len(np.where(np.isnan(tracks[param]))[0]) != 0:
    print('There are some nan values')
else:
    if write:
        csv_out= Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'_params.csv'
        print('write', csv_out)
        exist= os.path.isfile(csv_out) 

        if exist == False:  #create new file
            tracks[param].to_csv(csv_out)
            # with open(csv_out,'w') as file:
            #     writer = csv.DictWriter(file, delimiter=',', fieldnames= [fieldname])
            #     writer.writeheader()
            #     for index in range(len(data_array)):
            #         writer.writerow({fieldname: str(data_array[index])})
        
        if exist: #write into existing file
            csv_input = pd.read_csv(csv_out)
            csv_input= csv_input.set_index(['ID', 'step'])
            
            csv_input= pd.concat([csv_input, tracks[param]], axis=1)
            csv_input= csv_input.reset_index()
            # csv_input[param] = tracks[param].values
            csv_input.to_csv(csv_out, index=False)
       
   
#tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'_params.csv')

    
finish= time.perf_counter()

print(f'Finished in {round(finish - start,2)} seconds')
