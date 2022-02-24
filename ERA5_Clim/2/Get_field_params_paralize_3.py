#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 16:29:52 2019

@author: pst019
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


to_list= 'Yanase'
to_list= 'Smirnova'
# to_list= 'Noer'
to_list= 'Rojo'
#to_list= 'Golubkin'
to_list='tracks' #to which list parameters should be derived

# param= 'N_925-500'
param_var= 'theta_diff_925-500'
param_var= 'theta_trop'
#param_var= 'U_trop'
# param_var= 'U_trop_polew' #not calculated with param_rad and type

param_rad= 250
param_type= 'mean'
# param_type= 'max'

param= param_var+'_'+param_type+ str(param_rad)

time_start= False #in this case the first and last index of alltracks will be used
time_start= '2013-1-1' #the field data is needed for this time
time_end= '2013-2-1'

print(to_list)
print(param)


csv_name= 'merged_tracks_'+version+'_'+ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
if terrain_excl: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
if split: csv_name += f'_split{exceed_large_to_small}'

"""get all tracks"""
if to_list=='tracks':
    track_dir= Mediadir+'/mergedtracks/'+csv_name
    
else:
    """get matchPL list"""
    track_dir= Mediadir+f"/matchPLlist/matchPLlist_{to_list}_dist{matchdist}_{csv_name}"

tracks= pd.read_csv(track_dir+'.csv')

tracks['time']= pd.to_datetime(tracks['time'])
tracks= tracks.set_index(['ID', 'step'])

if time_start != False:
	tracks= tracks[(tracks['time'] >= time_start) & (tracks['time'] < time_end)]


tracks[param]= np.nan



time_steps = remove_dublicate(tracks.time)


def calc_param_to_list(time_step, tracks= tracks, param= param, param_var= param_var, param_type= param_type, workdir= workdir):
    if time_step.day == 1 and time_step.hour == 0:
        print(time_step)
    
    """import the parameter"""
    if param_var== 'theta_diff_925-500':
        filename= f"levels_500_era5_{time_step.year}_{time_step.month:02}.nc"
        if fram == True: filename= 'levels_500/'+filename
        ds= xr.open_dataset(workdir + filename )
        ds= ds.sel(time= time_step)
    
        theta_1, theta_0= PotTemp(ds['var130'].isel(plev=1), ds.plev[1]/100), PotTemp(ds['var130'].isel(plev=0), ds.plev[0]/100)
        # theta_1, theta_0= PotTemp(ds.t.isel(plev=1), ds.plev[1]/100), PotTemp(ds.t.isel(plev=0), ds.plev[0]/100)
        ds[param_var]= (('lat', 'lon'), theta_1 - theta_0)
        ds[param_var].attrs['units'] = 'K'
    
    
    if param_var== 'theta_trop':
        # ds= xr.open_dataset(workdir + f"{file_name}/{file_name}_era5_{time_step.year}_{time_step.month:02}.nc")
        filename= f"PV_era5_{time_step.year}_{time_step.month:02}.nc"
        if fram == True: filename= 'PV/'+filename
        ds= xr.open_dataset(workdir + filename )
        ds= ds.sel(time= time_step)   
        ds= ds.isel(lev= 0)
        
        # ds[param_var]= ds['var3'] #ds.pt
        ds[param_var]= ds.pt


    if param_var in ['U_trop', 'U_trop_polew']:
        # ds= xr.open_dataset(workdir + f"{file_name}/{file_name}_era5_{time_step.year}_{time_step.month:02}.nc")
        filename= f"PV_era5_{time_step.year}_{time_step.month:02}.nc"
        if fram == True: filename= 'PV/'+filename
        ds= xr.open_dataset(workdir + filename )
        ds= ds.sel(time= time_step)   
        ds= ds.isel(lev= 0)
        
        # ds[param_var]= np.sqrt(ds['var131']**2 +ds['var132']**2) #np.sqrt(ds.u**2 +ds.v**2)
        ds[param_var]= np.sqrt(ds.u**2 +ds.v**2)

    
    tracks_now= tracks[tracks.time== time_step]
    
    
    for i in range(len(tracks_now)): 
        # print(i)
        track= tracks_now.iloc[i]
           
        if param_var== 'U_trop_polew':
            tracks.loc[track.name, param]= float(ds[param_var].where((np.abs(ds.lon-track.lon) < 0.5)& (ds.lat >= track.lat)).max())

        else: #other variables where a map of the new variable could be plotted    
            # param_loc= ds[param_var].sel(lat= track.lat, lon= track.lon, method= 'nearest').values
            # print('At track location:', param_loc)
        
            dist= (110* np.sqrt( (ds.lat - track.lat)**2 + (np.cos(np.deg2rad(track.lat)) *(ds.lon- track.lon) )**2 ) )
            
            if param_type== 'mean':
                param_mean= float(ds[param_var].where(dist < param_rad).mean().values)
            elif param_type== 'max':
    	        param_mean= float(ds[param_var].where(dist < param_rad).max().values)            
            print('Mean around track location:', param_mean)
        
            # tracks.loc[track.name, param]= param_mean    
            # results[track.name]=  param_mean

    return [track.name, param_mean]


param_list= dict()
def record_param_list(res):
    trackname, param_mean = res
    param_list[trackname]= param_mean


# tracks_param= pd.DataFrame(index= 'ID', columns= [param])

# for time_step in time_steps:
#     calc_param_to_list(time_step, tracks, param, param_var, param_type, workdir)

    
# import concurrent.futures
# with concurrent.futures.ProcessPoolExecutor() as executor:
#     results= executor.map(calc_param_to_list, time_steps)
    

import multiprocessing as mp
# from multiprocessing import Manager


# manager = Manager()
# results = manager.dict()

# pool= mp.Pool(4)
# jobs= []


# for time_step in time_steps:
#     job= pool.apply(calc_param_to_list , (time_step, tracks, param, param_var, param_type, workdir) )
#     jobs.append(job)

# # for job in jobs:
# #     job.get()
    
# pool.close()
# pool.join()
n_processes= mp.cpu_count()

with mp.Pool(processes= n_processes) as pool:
    for time_step in time_steps:
        pool.apply_async(calc_param_to_list, args=(time_step, tracks, param, param_var, param_type, workdir),
                          callback= record_param_list)
    pool.close()
    pool.join()


    
if len(np.where(np.isnan(tracks[param]))[0]) != 0:
    print('There are some nan values')
else:
    if write:
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

print(f'Finished in {round(finish - start,2)} seconds')
