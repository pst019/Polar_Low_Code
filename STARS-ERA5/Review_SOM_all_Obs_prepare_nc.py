#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 10:01:45 2019

@author: pst019
"""


import time
start = time.time()


import os
user = os.getcwd().split('/')[2]

if user=='pst019':
#    homedir= '/home/'+user+'/home/'
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/1692A00D929FEF8B/'

homedir=Mediadir+'home/'


import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import xarray as xr
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from f_useful import *
from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs

#
write_nc= True
#write_nc= False


"""SOM var"""
plevel=850

var= 't'
#var='z'
#var='vo'
#var='U'
#var='q'

#PLCG_type, vel= 'stearing_flow', 3
PLCG_type, smooth_param= 'track_smth', '1E-3'




"""lifetime threshold"""
lifetime=12
windthresh= 20

"""all Obs to remove the nan values and the time steps with a low propagation speed"""
if PLCG_type== 'stearing_flow':
    file=Mediadir + "ERA5_STARS/PL_centred_fields/" +var + '_'+ str(plevel) + '_allObs.nc'
    if os.path.isfile(file):
        ds= xr.open_dataset(file)
        
    
elif PLCG_type== 'track_smth':
    file=Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +var + '_'+ str(plevel) + '_allObs_track-smth-'+smooth_param+'.nc'
    if os.path.isfile(file):
        ds= xr.open_dataset(file)

    else:
        file=Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +var + '_all_levs_allObs_track-smth-'+smooth_param+'.nc'
        ds= xr.open_dataset(file)
        ds= ds.sel(plev= plevel)        

    
print('Initial: number timesteps:', len(ds.time), 'number systems:', len(ds.groupby('PLnr').last()['PLnr']))



ds= ds.dropna(dim='time') #drops all time steps where a variable is nan, this happens for example if part of the interpolated field crosses the ERA-5 boundary. reduces from 13345 to 12826
print('After exclusion where domain is not fully captured: number timesteps:', len(ds.time), 'number systems:', len(ds.groupby('PLnr').last()['PLnr']))



ds[var+ '_ano']= (('time', 'x', 'y'), ds[var].values - np.mean(ds[var].values, axis= (1,2))[:, np.newaxis, np.newaxis])


"""get some characteristics from the Rojo list in order to sort the dataset"""

imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


Stoll= pd.read_csv(imp_dir + 'Stoll_list_prior_calcs.csv')
Stoll= Stoll.set_index(['ID', 'time'])

S_ERA5_2= pd.read_csv(imp_dir + 'Stoll_ERA5_2.csv')
S_ERA5_2= S_ERA5_2.set_index(['ID', 'time'])

Stoll= pd.concat([Stoll, S_ERA5_2], axis=1)


#exclude matches to secondary Rojo events
Stoll= Stoll[Stoll['Rojo_ID'].str[-2] != '.']

#exclude low wind speeds
Stoll_ind= Stoll['U_msk_max-250'].groupby('ID').max()
Stoll_ind= Stoll_ind[Stoll_ind > windthresh]


ID_int= [int(ID.split('_')[0])*1E3 + int(ID.split('_')[1])for ID in Stoll_ind.index] #makes integer of the PL nrs

"""exclude steps"""
if PLCG_type== 'stearing_flow':
    ds= ds.where(ds.steer_vel > vel, drop=True) #reduces from 12826 to 10989


ds['duration'] = ('time', np.zeros(len(ds['PLnr'])) )
ds['primary_Rojo']= ('time', np.zeros(len(ds['PLnr']) ) )


duration_ds=ds.groupby(ds['PLnr']).last()['Obs']
for PLnow in duration_ds.PLnr.values:
    ds['duration'][ds['PLnr']== PLnow] =duration_ds.sel(PLnr= PLnow).values

    if PLnow in ID_int:  ds['primary_Rojo'][ds['PLnr']== PLnow] = 1


ds= ds.where(ds.duration >= lifetime, drop=True)    

print('After lifetime exclusion: number timesteps:', len(ds.time), 'number systems:', len(ds.groupby('PLnr').last()['PLnr']))


ds= ds.where(ds.primary_Rojo == 1, drop=True)    

print('After primary + windspeed exclusion: number timesteps:', len(ds.time), 'number systems:', len(ds.groupby('PLnr').last()['PLnr']))


if write_nc:
    if PLCG_type== 'stearing_flow':
        filename= Mediadir + "ERA5_STARS/PL_centred_fields/" +var + '_'+ str(plevel) + '_allObs_SOMprep_vel'+str(vel)+'_dur'+str(lifetime)+'.nc'
    elif PLCG_type== 'track_smth':
        filename= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +var + '_'+ str(plevel) + '_allObs_primary_SOMprep_track-smth-'+smooth_param+'_dur'+str(lifetime)+'_Umax'+str(windthresh)+'.nc'

    
    
    print(filename)
    ds.to_netcdf(filename, format='NETCDF3_CLASSIC' )


"""remove the time steps with low propagation spped"""
##Obs= 'Obsnr1'
#Obs= 'mature'
#
#ds= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields/" +var + '_'+ str(plevel) + '_'+Obs+'.nc')
#ds= ds.dropna(dim='time')
#
##ds[var+ '_ano']= (('time', 'x', 'y'), ds[var].values - np.mean(ds[var].values, axis= (1,2))[:, np.newaxis, np.newaxis])
#
#ds= ds.where(ds.steer_vel > 3, drop=True) #for Obsnr1: 397 -> 342 , for mature -> 358
#ds.to_netcdf(Mediadir + "ERA5_STARS/PL_centred_fields/" +var + '_'+ str(plevel) + '_'+Obs +'_prep.nc',
#             format='NETCDF3_CLASSIC' )