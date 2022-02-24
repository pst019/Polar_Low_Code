#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 16:29:52 2019

@author: pst019

Run via mpirun -n 4 python Get_field_params_mpi.py
"""

import os
user = os.getcwd().split('/')[2]


import pandas as pd
import numpy as np
import xarray as xr
import time
import sys  #to import the functions from a different directory

start = time.perf_counter()



fram= False
if user=='pst019': #office
    Mediadir= '/media/pst019/PatsOrange/data/ERA5_Clim/track_lists/'
    workdir= '/media/pst019/PatsOrange/data/ERA5_Clim/ERA5_data/'
    sys.path.insert(0, '../Functions')

elif user == 'media': #home
    Mediadir= '/run/media/pst019/PatsOrange/data/ERA5_Clim/track_lists/'
    workdir= '/run/media/pst019/PatsOrange/data/ERA5_Clim/ERA5_data/'
    sys.path.insert(0, '../Functions')

else: #fram
    fram=True
    Mediadir= '/cluster/home/pst019/pmctrack/post_process/'
    workdir= '/cluster/work/users/pst019/ERA5_data/'

from f_useful import *
from f_meteo import *



write= False
csvoutname= 'params'

if fram: write= True


version= ""# "fram_run3_"
ending= 'atl-pac'

start_year= 1979 #80
end_year= 2020 #19

more_info= True
post_process= True


durlim= 6
matchdist= 150

split, exceed_large_to_small= True, 40
# split= False

terrain_excl, terrain_dist, terrain_thresh= True, 2, 50
# terrain_excl= False

time_diff= 24 #retains only the matched PL lists that do not start earlier or end later than 24 hours than the listed PL
# timediff= False

ds_name_lists= ['tracks']
#ds_name_lists= ['Yanase', 'Smirnova', 'Noer', 'Rojo', 'Golubkin', 'tracks']
#ds_name_lists= ['Smirnova', 'Noer', 'Rojo', 'Golubkin', 'tracks']
# ds_name_lists= ['Yanase']
#ds_name_lists= ['Smirnova', 'Rojo', 'tracks']
# ds_name_lists= ['Yanase','Noer', 'Golubkin']


#param_var, param_type= 'N_500-925', 'mean'
#param_var, param_type= 'theta_diff_500-925', 'mean'
#param_var, param_type= 'theta_diff_500-skt', 'mean'
# param_var, param_type= 'theta_diff_500-sst', 'mean'
# param_var, param_type= 'SST-T_500', 'mean'
#param_var, param_type= 'theta_trop', 'mean'
#param_var, param_type= 'U_trop', 'max'

# param_var= 'U_trop_polew' #not calculated with param_rad and type
#param_var, param_type='U_surface', 'max'

param_var, param_type= 'shear_925-500', 'mean'


param_rad= 250
# param_type= 'mean'
# param_type= 'max'

param= param_var+'_'+param_type+ str(param_rad)

print(param)


# param_var_list=[]
# param_type_list=[]

def get_track_for_list(ds_name, time_start, time_end, version= version, ending= ending,
                       durlim= durlim, terrain_dist= terrain_dist, terrain_thres= terrain_thresh, 
                       exceed_large_to_small= exceed_large_to_small):
    """get all tracks"""
    

    csv_name= 'merged_tracks_'+version + ending+f'_{start_year}-{end_year}'
    if durlim: csv_name += 'durlim'+str(durlim)
    if more_info: csv_name += '_moreinfo'
    if post_process: csv_name += '_post-process'
    if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
    # print(csv_name)
    
    # # tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks2/"+csv_name+'.csv')


    # csv_name= 'merged_tracks_'+version+ending+f'_{start_year}-{end_year}'
    # if durlim: csv_name += 'durlim'+str(durlim)
    # if terrain_excl: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
    # if split: csv_name += f'_split{exceed_large_to_small}'
    
    
    if ds_name=='tracks':
        track_dir= Mediadir+'/mergedtracks2/'+csv_name
        
    else:
        """get matchPL list"""
        track_dir= Mediadir+f"/matchPLlist/matchPLlist_{ds_name}_dist{matchdist}"
        if time_diff and ds_name != 'Yanase': track_dir += f"_timediff{time_diff}"
        track_dir += f"_{csv_name}"
    
    tracks= pd.read_csv(track_dir+'.csv')
    
    tracks['time']= pd.to_datetime(tracks['time'])
    tracks= tracks.set_index(['ID', 'step'])
    
    if time_start != False:
    	tracks= tracks[(tracks['time'] >= time_start) & (tracks['time'] < time_end)]

    return tracks, track_dir






def calc_param_to_ds(time_steps, tracks, param, param_var, param_type, workdir= workdir):
    """calculate the parameter for the tracks"""
    if param_var == 'shear_925-500':
        param1= param_var+'_angle_'+param_type+ str(param_rad)
        param2= param_var+'_strength_'+param_type+ str(param_rad)
        tracks[param1]= np.nan
        tracks[param2]= np.nan
        
    else:
        tracks[param]= np.nan


    for time_step in time_steps:
        if time_step.day == 1 and time_step.hour == 0:
            print(time_step)
        
        """import the parameter"""
        if param_var== 'theta_diff_500-925':
            filename= f"levels_500_era5_{time_step.year}_{time_step.month:02}.nc"
            if fram == True: filename= 'levels_500/'+filename
            ds= xr.open_dataset(workdir + filename )
            ds= ds.sel(time= time_step)
        
            theta_1, theta_0= PotTemp(ds['var130'].isel(plev=1), ds.plev[1]/100), PotTemp(ds['var130'].isel(plev=0), ds.plev[0]/100)
            # else: theta_1, theta_0= PotTemp(ds.t.isel(plev=1), ds.plev[1]/100), PotTemp(ds.t.isel(plev=0), ds.plev[0]/100)
            ds[param_var]= (('lat', 'lon'), theta_1 - theta_0)
            ds[param_var].attrs['units'] = 'K'

        if param_var== 'N_500-925': #not tested
            filename= f"levels_500_era5_{time_step.year}_{time_step.month:02}.nc"
            if fram == True: filename= 'levels_500/'+filename
            ds= xr.open_dataset(workdir + filename )
            ds= ds.sel(time= time_step)
        
            theta_1= PotTemp(ds['var130'].isel(plev=1), ds.plev[1]/100).values
            theta_0= PotTemp(ds['var130'].isel(plev=0), ds.plev[0]/100).values

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
    
            """get the lower level potential temperature"""
            filename= f"lsm_shift_update.nc"
            if fram == True: filename= 'tracking/'+filename
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
            
            """get the lower level potential temperature"""
            filename= f"skin_temp_era5_{time_step.year}_{time_step.month:02}.nc"
            if fram == True: filename= 'skin_temp/'+filename
            ds2= xr.open_dataset(workdir + filename )
            ds2= ds2.sel(time= time_step)
        
            """get the lsm and exclude land and sea ice"""
            filename= f"lsm_shift_update.nc"
            if fram == True: filename= 'tracking/'+filename
            ds4= xr.open_dataset(workdir + filename )
            ds4= ds4.sel(lat= ds4.lat[::2], lon= ds4.lon[::2])
            
            ds2['var235'].values[ds4.lsm.values != 0]= np.nan #excludes land
            ds2['var235'].values[ds2['var235'].values <= 273.15]= np.nan #excludes sea ice
            
            ds[param_var]= (('lat', 'lon'), ds2['var235'] - ds['var130'].isel(plev=1))
            ds[param_var].attrs['units'] = 'K'

        
        
        if param_var== 'theta_trop':
            # ds= xr.open_dataset(workdir + f"{file_name}/{file_name}_era5_{time_step.year}_{time_step.month:02}.nc")
            filename= f"PV_era5_{time_step.year}_{time_step.month:02}.nc"
            if fram == True: filename= 'PV/'+filename
            ds= xr.open_dataset(workdir + filename )
            ds= ds.sel(time= time_step)   
            ds= ds.isel(lev= 0)
            
            if fram: ds[param_var]= ds['var3'] #ds.pt
            else: ds[param_var]= ds.pt
    
    
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

        if param_var== 'shear_925-500':
            #get the strength and the angle
            
            ds= xr.open_dataset(workdir + f"levels_500_era5_{time_step.year}_{time_step.month:02}.nc")
            ds= ds.sel(time= time_step)
    
            file_name_z= 'levels_z'
            ds_z= xr.open_dataset(workdir + f"levels_z_era5_{time_step.year}_{time_step.month:02}.nc")
            ds_z= ds_z.sel(time= time_step)

        
        tracks_now= tracks[tracks.time== time_step]
        
        
        for i in range(len(tracks_now)): 
            # print(i)
            track= tracks_now.iloc[i]
               
            if param_var== 'U_trop_polew':
                tracks.loc[track.name, param]= float(ds[param_var].where((np.abs(ds.lon-track.lon) < 0.5)& (ds.lat >= track.lat)).max())

            elif param_var == 'shear_925-500':
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
                
                tracks.loc[track.name, param1]= shear_angle
                tracks.loc[track.name, param2]= shear_strength
                # tracks.loc[track.name, param]= [shear_angle, shear_strength]

    
            else: #other variables where a map of the new variable could be plotted    
                # param_loc= ds[param_var].sel(lat= track.lat, lon= track.lon, method= 'nearest').values
                # print('At track location:', param_loc)
            
                dist= (110* np.sqrt( (ds.lat - track.lat)**2 + (np.cos(np.deg2rad(track.lat)) *(ds.lon- track.lon) )**2 ) )
                
                if param_type== 'mean':
                    param_mean= float(ds[param_var].where(dist < param_rad).mean().values)
                elif param_type== 'max':
            	        param_mean= float(ds[param_var].where(dist < param_rad).max().values)
                elif param_type== 'min':
            	        param_mean= float(ds[param_var].where(dist < param_rad).min().values)                           
                # print('Mean around track location:', param_mean)
            
                tracks.loc[track.name, param]= param_mean    
                # param_list[track.name]=  param_mean
                # print(param_list)
                
    # return param_list
    return tracks


for ds_name in ds_name_lists:
    print(ds_name)
    
    # if fram == False or ds_name== 'tracks':
    #     time_start= '2013-1-1' #the field data is needed for this time
    #     time_end= '2013-12-31'

    if fram == False or ds_name== 'tracks':
        time_start= '2008-1-1' #the field data is needed for this time
        if fram: time_end= '2009-1-1'
        else: time_end= '2008-1-3'


    else:
        time_start= False #the field data is needed for this time
        time_end= False
        
    
    tracks, track_dir= get_track_for_list(ds_name, time_start= time_start, time_end= time_end)

    time_steps = np.array(remove_dublicate(tracks.time))

    tracks= calc_param_to_ds(time_steps, tracks, param, param_var, param_type, workdir)


   

        
    # if len(np.where(np.isnan(tracks[param]))[0]) != 0:
    #     print('There are some nan values')
    # else:
    if write and len(tracks)> 0:
        csv_out= track_dir+'_'+csvoutname
        if time_start != False: csv_out+= f'_{time_start}-{time_end}'
        csv_out += '.csv'

        print('write', csv_out)
        exist= os.path.isfile(csv_out) 

        if exist == False:  #create new file
            tracks[param].to_csv(csv_out)
        
        if exist: #write into existing file
            csv_input = pd.read_csv(csv_out)
            csv_input= csv_input.set_index(['ID', 'step'])
            
            csv_input= pd.concat([csv_input, tracks[param]], axis=1)
            csv_input= csv_input.reset_index()
            # csv_input[param] = tracks[param].values
            csv_input.to_csv(csv_out, index=False)
           
       
    
        
    finish= time.perf_counter()
    
    print(f'Time passed: {round(finish - start,2)} seconds')
