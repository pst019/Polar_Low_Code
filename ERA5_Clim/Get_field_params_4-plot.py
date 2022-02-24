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


workdir= Mediadir + 'data/ERA5_Clim/ERA5_data/'
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

fram= False

# savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/Investigate_PL/'
# save= False


version= "fram_run3"
ending= 'atl-pac'

start_year= 1980
end_year= 2019

durlim= 6
matchdist= 150

split, exceed_large_to_small= True, 40
# split= False

terrain_excl, terrain_dist, terrain_thresh= True, 2, 50
# terrain_excl= False



# 
to_list= 'Smirnova'
to_list= 'Rojo'

to_list='tracks' #to which list parameters should be derived



# param_var= 'N_500-925'
# param_var= 'theta_diff_500-925'
# param_var= 'theta_diff_500-skt'
# param_var= 'theta_diff_500-sst'
# param_var= 'SST-T_500'
# param_var= 'theta_trop'
# param_var= 'U_trop'
# param_var= 'shear_925-500'
# param_var= 'U_trop_polew' #not calculated with param_rad and type
# param_var= 'U_surface'

param_rad= 250
param_type= 'mean'
# param_type= 'max'

# param_var, param_type, param_rad= 'shear_925-500', 'mean', 500

param_var, param_type= 'land_dist', False


time_start= False #in this case the first and last index of alltracks will be used
# time_start= '2013-1-1' #the field data is needed for this time
# time_end= '2013-1-3'

time_start= '2008-1-7' #the field data is needed for this time
time_end= '2008-1-8'

csv_name= 'merged_tracks_'+version+'_'+ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
if terrain_excl: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
if split: csv_name += f'_split{exceed_large_to_small}'


"""get all tracks"""
if to_list=='tracks':
    tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv')


else:
    """get matchPL list"""
    track_dir= Mediadir+f"/data/ERA5_Clim/track_lists/matchPLlist/matchPLlist_{to_list}_dist{matchdist}_{csv_name}.csv"
    tracks= pd.read_csv(track_dir)



tracks['time']= pd.to_datetime(tracks['time'])
tracks= tracks.set_index(['ID', 'step'])

if time_start != False:
    tracks= tracks[(tracks['time'] >= time_start) & (tracks['time'] < time_end)]






if param_var == 'shear_925-500':
    param1= param_var+'_angle_'+param_type+ str(param_rad)
    param2= param_var+'_strength_'+param_type+ str(param_rad)
    tracks[param1]= np.nan
    tracks[param2]= np.nan
    
elif param_type:
    param= param_var+'_'+param_type+ str(param_rad)
    tracks[param]= np.nan

else:
    param= param_var
    tracks[param]= np.nan


# indexes= remove_dublicate(tracks.index.get_level_values(0)

time_steps = remove_dublicate(tracks.time)
# if Plot:
time_steps= time_steps[:1]

for time_step in time_steps:
    if time_step.day == 1 and time_step.hour == 0:
        print(time_step)
    
    """import the parameter"""
    if param_var== 'land_dist':
        filename= "lsm25_orig.nc"
        if fram == True: filename= 'lsm/'+filename
        ds= xr.open_dataset(workdir + filename )
        ds= ds.isel(time= 0)

        
    if param_var== 'theta_diff_500-925':
        file_name= 'levels_500'
        ds= xr.open_dataset(workdir + f"{file_name}_era5_{time_step.year}_{time_step.month:02}.nc")
        ds= ds.sel(time= time_step)
    
        theta_1, theta_0= PotTemp(ds.t.isel(plev=1), ds.plev[1]/100), PotTemp(ds.t.isel(plev=0), ds.plev[0]/100)

        ds[param_var]= (('lat', 'lon'), theta_1 - theta_0)
        ds[param_var].attrs['units'] = 'K'


    if param_var== 'N_500-925': #not tested
        filename= f"levels_500_era5_{time_step.year}_{time_step.month:02}.nc"
        if fram == True: filename= 'levels_500/'+filename
        ds= xr.open_dataset(workdir + filename )
        ds= ds.sel(time= time_step)
    
        # if fram:
        theta_1= PotTemp(ds['var130'].isel(plev=1), ds.plev[1]/100).values
        theta_0= PotTemp(ds['var130'].isel(plev=0), ds.plev[0]/100).values
        # else:
        #     theta_1 = PotTemp(ds.t.isel(plev=1), ds.plev[1]/100)
        #     theta_0 = PotTemp(ds.t.isel(plev=0), ds.plev[0]/100)
            
        filename= f"levels_z_era5_{time_step.year}_{time_step.month:02}.nc"
        if fram == True: filename= 'levels_z/'+filename
        ds2= xr.open_dataset(workdir + filename )
        ds2 = ds2.sel(time= time_step)

        h_diff=  (ds2['var129'].isel(plev=1) - ds2['var129'].isel(plev=0) ).values/9.81

        ds[param_var]= (('lat', 'lon'), np.sqrt(9.81/((theta_1+theta_0)/2) *  (theta_1- theta_0)/h_diff))
        ds[param_var].attrs['units'] = '1/s'
            

    if param_var== 'theta_diff_500-skt':
        filename= f"levels_500_era5_{time_step.year}_{time_step.month:02}.nc"
        if fram == True: filename= 'levels_500/'+filename
        ds= xr.open_dataset(workdir + filename )
        ds= ds.sel(time= time_step)
    
        theta_1= PotTemp(ds['var130'].isel(plev=1), ds.plev[1]/100)
        # else: theta_1, theta_0= PotTemp(ds.t.isel(plev=1), ds.plev[1]/100)
        
        """get the lower level potential temperature"""
        filename= f"skin_temp_era5_{time_step.year}_{time_step.month:02}.nc"
        if fram == True: filename= 'skin_temp/'+filename
        ds2= xr.open_dataset(workdir + filename )
        ds2= ds2.sel(time= time_step)
    
        filename= f"slp_shift_era5_{time_step.year}_{time_step.month:02}.nc"
        if fram == True: filename= 'tracking/'+filename
        ds3= xr.open_dataset(workdir + filename )
        ds3= ds3.sel(time= time_step)
        ds3.sel(lat= ds3.lat[::2], lon= ds3.lon[::2]) #due to a higher horizontal resolution
    
        theta_0= PotTemp(ds2['var235'], ds3['var151']/100)
        # else: theta_0= PotTemp(ds2['var235'], ds3['msl']/100)
        
        ds[param_var]= (('lat', 'lon'), theta_1 - theta_0)
        ds[param_var].attrs['units'] = 'K'


    if param_var== 'theta_diff_500-sst':
        """get sst by masking skt"""
        filename= f"levels_500_era5_{time_step.year}_{time_step.month:02}.nc"
        if fram == True: filename= 'levels_500/'+filename
        ds= xr.open_dataset(workdir + filename )
        ds= ds.sel(time= time_step)
    
        theta_1= PotTemp(ds['var130'].isel(plev=1), ds.plev[1]/100)
        # else: theta_1, theta_0= PotTemp(ds.t.isel(plev=1), ds.plev[1]/100)
        
        """get the lower level potential temperature"""
        filename= f"skin_temp_era5_{time_step.year}_{time_step.month:02}.nc"
        if fram == True: filename= 'skin_temp/'+filename
        ds2= xr.open_dataset(workdir + filename )
        ds2= ds2.sel(time= time_step)
    
        filename= f"slp_shift_era5_{time_step.year}_{time_step.month:02}.nc"
        if fram == True: filename= 'tracking/'+filename
        ds3= xr.open_dataset(workdir + filename )
        ds3= ds3.sel(time= time_step)
        ds3.sel(lat= ds3.lat[::2], lon= ds3.lon[::2]) #due to a higher horizontal resolution
    
        theta_0= PotTemp(ds2['var235'], ds3['var151']/100)
        # else: theta_0= PotTemp(ds2['var235'], ds3['msl']/100)

        """get the lower level potential temperature"""
        filename= "lsm_shift_update.nc"
        if fram == True: filename= 'lsm/'+filename
        ds4= xr.open_dataset(workdir + filename )
        ds4= ds4.sel(lat= ds4.lat[::2], lon= ds4.lon[::2])
        theta_0.values[ds4.lsm.values != 0]= np.nan #excludes land
        theta_0.values[ds2['var235'].values <= 273.15]= np.nan #excludes sea ice
        
        ds[param_var]= (('lat', 'lon'), theta_1 - theta_0)
        ds[param_var].attrs['units'] = 'K'


    if param_var== 'SST-T_500':
        """get sst by masking skt"""
        filename= f"levels_500_era5_{time_step.year}_{time_step.month:02}.nc"
        if fram == True: filename= 'levels_500/'+filename
        ds= xr.open_dataset(workdir + filename )
        ds= ds.sel(time= time_step)
    
        theta_1= PotTemp(ds['var130'].isel(plev=1), ds.plev[1]/100)
        # else: theta_1, theta_0= PotTemp(ds.t.isel(plev=1), ds.plev[1]/100)
        
        """get the lower level potential temperature"""
        filename= f"skin_temp_era5_{time_step.year}_{time_step.month:02}.nc"
        if fram == True: filename= 'skin_temp/'+filename
        ds2= xr.open_dataset(workdir + filename )
        ds2= ds2.sel(time= time_step)
    

        # theta_0= PotTemp(ds2['var235'], ds3['var151']/100)
        # else: theta_0= PotTemp(ds2['var235'], ds3['msl']/100)

        """get the lsm and exclude land and sea ice"""
        filename= "lsm_shift_update.nc"
        if fram == True: filename= 'lsm/'+filename
        ds4= xr.open_dataset(workdir + filename )
        ds4= ds4.sel(lat= ds4.lat[::2], lon= ds4.lon[::2])
        
        ds2['var235'].values[ds4.lsm.values != 0]= np.nan #excludes land
        ds2['var235'].values[ds2['var235'].values <= 273.15]= np.nan #excludes sea ice
        
        ds[param_var]= (('lat', 'lon'), ds2['var235'] - ds['var130'].isel(plev=1))
        ds[param_var].attrs['units'] = 'K'


    if param_var== 'shear_925-500':
        #get the strength and the angle
        
        ds= xr.open_dataset(workdir + f"levels_500_era5_{time_step.year}_{time_step.month:02}.nc")
        ds= ds.sel(time= time_step)

        file_name_z= 'levels_z'
        ds_z= xr.open_dataset(workdir + f"levels_z_era5_{time_step.year}_{time_step.month:02}.nc")
        ds_z= ds_z.sel(time= time_step)

    
    
    if param_var== 'theta_trop':
        file_name= 'PV'
        ds= xr.open_dataset(workdir + f"{file_name}_era5_{time_step.year}_{time_step.month:02}.nc")
        ds= ds.sel(time= time_step)   
        ds= ds.isel(lev= 0)
        
        ds[param_var]= ds.pt

    if param_var in ['U_trop', 'U_trop_polew']:
        # ds= xr.open_dataset(workdir + f"{file_name}/{file_name}_era5_{time_step.year}_{time_step.month:02}.nc")
        filename= f"PV_era5_{time_step.year}_{time_step.month:02}.nc"
        if fram == True: filename= 'PV/'+filename
        ds= xr.open_dataset(workdir + filename )
        ds= ds.sel(time= time_step)   
        ds= ds.isel(lev= 0)
        
        if fram: ds[param_var]= np.sqrt(ds['var131']**2 +ds['var132']**2) #np.sqrt(ds.u**2 +ds.v**2)
        else: ds[param_var]= np.sqrt(ds.u**2 +ds.v**2)

    if param_var == 'U_surface':
        filename= f"surface_wind_era5_{time_step.year}_{time_step.month:02}.nc"
        if fram == True: filename= 'surface_wind/'+filename
        ds= xr.open_dataset(workdir + filename )
        ds= ds.sel(time= time_step)   
        # ds= ds.isel(lev= 0)            
        
        ds[param_var]= np.sqrt(ds['var165']**2 +ds['var166']**2) #np.sqrt(ds.u**2 +ds.v**2)



    
    tracks_now= tracks[tracks.time== time_step]
    
    if Plot:
        tracks_now= tracks_now[10:11]
        print('Calculate only for one time step')
        fig = plt.figure(fignr, figsize= (10,8))
        plt.clf()
        # ax= Plot_Polar_Stereo(fig, central_longitude= 20, extent= [-10, 40, 60, 77], subplot= (1,1,1), scalebar=False)
        # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-45, 45, 50, 80], subplot= (1,1,1))
        ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, 30, 80], subplot= (1,1,1))
        # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-50, 80, 30, 80], subplot= (1,1,1))
        

    
    for i in range(len(tracks_now)): 
        # print(i)
        track= tracks_now.iloc[i]
          
   
        if param_var == 'shear_925-500':
            u_0= value_in_rad(ds['var131'].isel(plev= 0), ds.lat, ds.lon, track.lat, track.lon, distance= param_rad)
            v_0= value_in_rad(ds['var132'].isel(plev= 0), ds.lat, ds.lon, track.lat, track.lon, distance= param_rad)
            u_1= value_in_rad(ds['var131'].isel(plev= 1), ds.lat, ds.lon, track.lat, track.lon, distance= param_rad)
            v_1= value_in_rad(ds['var132'].isel(plev= 1), ds.lat, ds.lon, track.lat, track.lon, distance= param_rad)

            z_0= value_in_rad(ds_z['var129'].isel(plev= 0), ds_z.lat, ds_z.lon, track.lat, track.lon, distance= param_rad)
            z_1= value_in_rad(ds_z['var129'].isel(plev= 1), ds_z.lat, ds_z.lon, track.lat, track.lon, distance= param_rad)
            
            du= u_1- u_0
            dv= v_1- v_0
       
            differential_angle= UV2Direction(du, dv) #in compass direction
            mean_angle= UV2Direction(u_0 + u_1, v_0 + v_1)
            
            shear_angle= (differential_angle - mean_angle)%360
            
            dU= np.sqrt(du**2 + dv**2)
            dH= (z_1 - z_0)/9.81

            shear_strength=  dU/dH
            
            print('shear strength (10-3):', shear_strength*1E3, '\n shear angle:', shear_angle)
            print('(u_1, v_1)', np.round(u_1,1), np.round(v_1,1), '\n(u_0, v_0)', np.round(u_0,1), np.round(v_0,1)) 
            
            tracks.loc[track.name, param1]= shear_angle
            tracks.loc[track.name, param2]= shear_strength


            if Plot:
                print('param', param1)
                # vextr= np.max([np.max(ds[param]), -np.min(ds[param])])
                # cf= ax.contourf(ds.lon, ds.lat, ds[param], transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)
                cf= ax.contourf(ds.lon, ds.lat, ds['var130'].isel(plev= 0), transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r')
            
                cf= ax.contour(ds_z.lon, ds_z.lat, ds_z['var129'].isel(plev= 0)/9.81, transform= ccrs.PlateCarree(), levels= 15, colors= 'k') #, vmin= vmin, vmax= vmax)
                plt.clabel(cf, fontsize=10, inline=1, fmt='%1.0f')

                cf= ax.contour(ds_z.lon, ds_z.lat, ds_z['var129'].isel(plev= 1)/9.81, transform= ccrs.PlateCarree(), levels= 15, colors= 'purple') #, vmin= vmin, vmax= vmax)
                plt.clabel(cf, fontsize=10, inline=1, fmt='%1.0f')
                
                # cb= fig.colorbar(cf, ax= ax, shrink=0.7)
                # varlabel= param_var
                # cb.set_label(varlabel , size=11)    
     

        # param_loc= ds[param].sel(lat= track.lat, lon= track.lon, method= 'nearest').values
        
        # print('At track location:', param_loc)

        elif param_var== 'U_trop_polew':
            param_value = float(ds[param_var].where((np.abs(ds.lon-track.lon) < 0.5)& (ds.lat >= track.lat)).max())
            tracks.loc[track.name, param]= param_value
            
            if Plot:
                print('param', param_value)
                vmax, vmin= np.max(ds[param_var]), np.min(ds[param_var])           
                cf= ax.contourf(ds.lon, ds.lat, ds[param_var], transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r', vmin= vmin, vmax= vmax)
                # cf= ax.contourf(ds.lon, ds.lat, ds[param_var].where(dist < 500) , transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r', vmin= vmin, vmax= vmax)           
                cb= fig.colorbar(cf, ax= ax, shrink=0.7)
                varlabel= param_var
                cb.set_label(varlabel , size=11)    
            



        elif param_var== 'land_dist':
            # if track.lon >= 180: track.lon -= 360
            dist= (110* np.sqrt( (ds.lat - track.lat)**2 + (np.cos(np.deg2rad(track.lat)) *((ds.lon)%360- (track.lon)%360 ) )**2 ) ).values

            # tracks.loc[track.name, param]= np.min(dist[ds.lsm.values > 0])   
            tracks.loc[track.name, param]= np.min(dist[ds['var172'].values > 0])   #exclude distances over water

            if Plot:
                print('param', param_var, tracks.loc[track.name, param])

                # cf= ax.contourf(ds.lon, ds.lat, ds['lsm'], transform= ccrs.PlateCarree(), levels= 15, cmap= 'Reds')
                cf= ax.contourf(ds.lon, ds.lat, ds['var172'], transform= ccrs.PlateCarree(), levels= 15, cmap= 'Reds')

                # cf= ax.contourf(ds.lon, ds.lat, ds[param_var].where(dist < 500) , transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r', vmin= vmin, vmax= vmax)
            
                cb= fig.colorbar(cf, ax= ax, shrink=0.7)
                varlabel= param_var
                cb.set_label(varlabel , size=11)    
            
          
            
    
        else: #other variables where a map of the new variable could be plotted
            dist= (110* np.sqrt( (ds.lat - track.lat)**2 + (np.cos(np.deg2rad(track.lat)) *(ds.lon- track.lon) )**2 ) )
            # dist= distance((ds.lat, ds.lon), (track.lat, track.lon))
            
            if param_type== 'mean':
                param_mean= np.nanmean(ds[param_var].where(dist < param_rad))
            elif param_type== 'max':
                param_mean= np.nanmax(ds[param_var].where(dist < param_rad))         
            # print('Mean around track location:', param_mean)
        
            if param_var in ['theta_diff_500-sst', 'SST-T_500']:
                if np.isnan(param_mean):
                    param_mean= -1E3
            # tracks_now[param]= param_mean
            tracks.loc[track.name, param]= param_mean    
        
            
            if Plot:
                print('param', param_mean)          
                vmax, vmin= np.max(ds[param_var]), np.min(ds[param_var])           
                cf= ax.contour(ds.lon, ds.lat, ds[param_var], transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r', vmin= vmin, vmax= vmax)
                cf= ax.contourf(ds.lon, ds.lat, ds[param_var].where(dist < 500) , transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r', vmin= vmin, vmax= vmax)           
                cb= fig.colorbar(cf, ax= ax, shrink=0.7)
                varlabel= param_var
                cb.set_label(varlabel , size=11)    
        
        
        if Plot:
            ax.plot(track.lon, track.lat, marker='o', c= 'g',  transform= ccrs.PlateCarree() , label= track.name )
            plt.legend()
       
     

if param_var == 'shear_925-500':
    param= [param1, param2]
# else:
    
print(param)
    
if len(np.where(np.isnan(tracks[param]))[0]) != 0:
    print('There are some nan values')
else:
    if write:
        csv_out= track_dir+'_params'
        if time_start != False: csv_out+= f'_{time_start}-{time_end}'
        csv_out += '.csv'
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
       
   
#tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+track_dir+'_params.csv')

    
finish= time.perf_counter()

print(f'Finished in {round(finish - start,2)} seconds')
