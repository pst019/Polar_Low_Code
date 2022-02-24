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
    homedir= '/home/'+user+'/home/'
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/1692A00D929FEF8B/'
'/home/'+user+'/home/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from f_useful import *
from f_meteo import *
from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs
from scipy.interpolate import griddata
import netCDF4 as nc
from datetime import datetime
import great_circle_calculator.great_circle_calculator as gcc # pip install great-circle-calculator


imp_lists=True
#imp_lists=False

imp_data=True
#imp_data=False

nc_save=True
#nc_save=False

#track_type='Rojo'
track_type='Stoll'


radius= 500 #radius of the local grid km
grid_dist= 25 #in km
ncells= radius//grid_dist
dist_grid_axis= np.arange(-ncells, ncells+1)*grid_dist

steering_rad= 200


interp='' #linear interpolation is done - this leads to some nan values at edges
#interp='_bd_nearest' #if the interpolated grid exceeds the ERA5 grid, neareast interpolation is done - however this makes field less smooth


"""specify system ID, Obsnr"""
#Obs= 1


import xarray as xr

plevel=None
vmin,vmax=None,None

stearing_file='tracking'
stearing_level= 850

"""specify field"""
#filetype='tracking'
#plevel= 850
#plevel= 'all_levs'
#var='vo'
#var='z'
#var='U'
#var='q'
#var='t'

#filetype='plev_w'
#plevel= 850
#plevel= 'all_levs'
#var='d' #divergence
#var='w'
#var='pv'

#filetype='boundary'
##var= 'cape'
#var= 'blh'
#var='mcc'
#var='lcc'
#var='hcc'
#var='sst'
#plevel=''

#filetype='surface'
#var='ci'
#var='skt'
#var='msl'
#plevel=''

filetype='forecast'
var='cp'
#var='tp'
#var='cin'
#var='slhf'
plevel=''

#filetype= 'lsm'
#var='lsm'
#plevel=''

#filetype='PV' #all variables should have a '_pv' at the end. this is necessary at the dataset merging
#var='pres_pv'
#var= 'z_pv'
#plevel=''

print('var', var)

if imp_lists:
    if track_type=='Rojo':
        S = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv")
    
    elif track_type== 'Stoll':
        test= 'version4'
        dist= 150
        
        Stoll_STARS_file="ERA5_STARS/STARS-PMC-merge/"+test+"_dist_"+str(dist)+"_handfix.csv"
        Stoll = pd.read_csv(Mediadir+ Stoll_STARS_file)
        Stoll['time']= pd.DatetimeIndex(Stoll.time)
        
        Stoll= Stoll_excl(Stoll)
        Stoll= Stoll_Obs_nr(Stoll) #get the Obs nr
    
        Stoll= Stoll.rename(columns={"Stoll nr": "ID"})
        S= Stoll


"""get the ERA-5 data"""
if imp_data:
    
#    for i, ID in enumerate(remove_dublicate(S['ID'])):
    for i in range(len(S)):
        print(i) #, ID)
    
#        S_now= S[np.logical_and(S['ID'] == ID, S['Obs'] == Obs ) ]
        S_now= S.iloc[i]
        
        if track_type=='Rojo':  S_now.time= S_now.time.dt.round("H")      #this has to be tested, should give a Timestamp
        
        filetime= str(S_now.time.year)+'_'+str(S_now.time.month).zfill(2)+'_'+str(S_now.time.day).zfill(2)
        hour= S_now.time.hour
    
        print(S_now.time)
    
        """general"""
        if filetype != 'lsm': #the normal case for not lsm
            ds0= xr.open_dataset(Mediadir + "ERA5_STARS/data/"+filetype+ "_era5_"+ filetime + '.nc')
            if filetype== 'PV':
                ds0= ds0.isel(lev=0)
                ds0= ds0.rename({'z': 'z_pv', 'u': 'u_pv', 'v': 'v_pv', 'pres': 'pres_pv', 'pt': 'pt_pv'})
                
            if stearing_file != filetype:
                ds2= xr.open_dataset(Mediadir + "ERA5_STARS/data/"+stearing_file+ "_era5_"+ filetime + '.nc')
                ds0 = xr.merge([ds0, ds2])
            ds0= ds0.isel(time= hour)


        if filetype== 'lsm': 
            ds0= xr.open_dataset(Mediadir + "ERA5_STARS/data/"+filetype+ '.nc')
#            ds0= ds0.isel(time= 0)
            ds0= ds0.squeeze().drop('time') #this removes the time coordinate and avoids problems at the merging
            ds0= ds0.rename({'latitude': 'lat', 'longitude': 'lon'})
            
            ds2= xr.open_dataset(Mediadir + "ERA5_STARS/data/"+stearing_file+ "_era5_"+ filetime + '.nc')
            ds2= ds2.isel(time= hour)

            ds0 = xr.merge([ds0, ds2])
        
        
        ds0['plev']/= 100

        
        """get the center position, the steering vector"""
        center_lon, center_lat= S_now.lon, S_now.lat
        
        u_steer= value_in_rad(ds0['u'].sel(plev= [1000, 700]).mean(axis= 0), ds0.lat, ds0.lon, center_lat, center_lon, steering_rad)
        v_steer= value_in_rad(ds0['v'].sel(plev= [1000, 700]).mean(axis= 0), ds0.lat, ds0.lon, center_lat, center_lon, steering_rad)

        steer_vel= np.sqrt(u_steer**2 + v_steer**2)
#        if steer_vel < 3:
#            print('Should be considered if the rotation works with such a steering wind of {} km/h'.format(np.round(steer_vel*3.6)) )


        beering= UV2Direction(u_steer, v_steer) #wind direction

        #the tangental axis in the strearing wind direction
        tanax= [list(gcc.point_given_start_and_bearing((center_lon, center_lat), beering, n*grid_dist*1E3)) for n in np.arange(-ncells, ncells+1)]
        tanax= np.array(tanax)
        
        #the beering along the tangential axis
        point0= gcc.point_given_start_and_bearing((center_lon, center_lat), beering, -(radius+grid_dist)*1E3) #point one before the start of the tangential axis, used for calculation of beering_axis
        beering_axis= [gcc.bearing_at_p2((point0), (tanax[m,0], tanax[m,1])) for m in range(len(tanax))] #the wind direction along the tangential axis
        
        #creation of the PL centred grid
        center_grid= [[list(gcc.point_given_start_and_bearing((tuple(tanax[m])), (beering_axis[m]-90)%360, n*grid_dist*1E3)) for m in range(len(tanax))] for n in np.arange(-ncells, ncells+1)]
        center_grid= np.array(center_grid)



        """ interpolate to local grid"""
        if type(plevel) == int: ds0= ds0.sel(plev= plevel)
        if var=='z': ds0[var]/=9.81
        if var=='vo': ds0[var] *=1E5
        if var=='q': ds0[var] *=1E3
        if var=='U': ds0[var]= np.sqrt(ds0['u']**2+ ds0['v']**2)
        
        if i == 0: longrid, latgrid= np.meshgrid(ds0.lon, ds0.lat)
        
        if type(plevel) == int or plevel == '':
            var_interp0= griddata(tuple([np.ravel(longrid), np.ravel(latgrid)]) , np.ravel(ds0[var]), (center_grid[:,:,0], center_grid[:,:,1]), method='linear')
            
            if interp== '_bd_nearest':
                #if the interpolated domain crosses the boundary of the ERA-5 domain, but less than 10% of the area is outside, nearest neighboor interpolation is done. This way the part outside the domain is not set to zero.
                if len(np.where(np.isnan(var_interp0))[0]) > 0 and len(np.where(np.isnan(var_interp0))[0]) < 0.1* (2*ncells +1)**2:
                    var_interp0= griddata(tuple([np.ravel(longrid), np.ravel(latgrid)]) , np.ravel(ds0[var]), (center_grid[:,:,0], center_grid[:,:,1]), method='nearest')
                    print('do nearest neighbor')
                elif len(np.where(np.isnan(var_interp0))[0]) >= 0.1* (2*ncells +1)**2:
                    print('too many interp gridcells outside')
        
        else: #for all plevels
            var_interp0= np.array([griddata(tuple([np.ravel(longrid), np.ravel(latgrid)]) , np.ravel(ds0[var][ip]), (center_grid[:,:,0], center_grid[:,:,1]), method='linear') for ip in range(len(ds0.plev))])
            
            if interp== '_bd_nearest':
                #if the interpolated domain crosses the boundary of the ERA-5 domain, but less than 10% of the area is outside, nearest neighboor interpolation is done. This way the part outside the domain is not set to zero.
                if len(np.where(np.isnan(var_interp0))[0]) > 0 and len(np.where(np.isnan(var_interp0))[0]) < 0.1* (2*ncells +1)**2:
                    var_interp0= np.array([griddata(tuple([np.ravel(longrid), np.ravel(latgrid)]) , np.ravel(ds0[var][ip]), (center_grid[:,:,0], center_grid[:,:,1]), method='nearest') for ip in range(len(ds0.plev))])
                    print('do nearest neighbor')
                elif len(np.where(np.isnan(var_interp0))[0]) >= 0.1* (2*ncells +1)**2:
                    print('too many interp gridcells outside')

        """save the relevant variables"""    
        if i == 0: #to create all the lists in the first iteration
            var_interp=np.zeros((0,) +var_interp0.shape )
            time_list= [S_now.time] #time vector for the output netcdf file
#            PLnr_list= [ID]
            PLnr_list_float= [int(S_now.ID.split('_')[0])*1E3 + int(S_now.ID.split('_')[1])]
            Obs_list= [S_now.Obs]
            steer_vel_list= [steer_vel]
            
#        if len(np.argwhere(np.isnan(var_interp0))) == 0:    #no nan values in data array, can be due to boundaries of domain
        var_interp= np.vstack((var_interp, var_interp0[np.newaxis]))
        if i != 0:
            time_list += [S_now.time]
#                PLnr_list += [ID]
            PLnr_list_float += [int(S_now.ID.split('_')[0])*1E3 + int(S_now.ID.split('_')[1])]
            Obs_list += [S_now.Obs]
            steer_vel_list += [steer_vel]
#        else: print('array contains nan values')

       
if nc_save:
    #to get the calendar and the time units
    if filetype != 'lsm': ncdstime= nc.Dataset(Mediadir + "ERA5_STARS/data/"+filetype+ "_era5_"+ filetime + '.nc')
    else: ncdstime= nc.Dataset(Mediadir + "ERA5_STARS/data/"+stearing_file+ "_era5_"+ filetime + '.nc')
    calendar= ncdstime['time'].calendar
    timeunits= ncdstime['time'].units
#    cdftime = nc.netcdftime.utime(timeunits,calendar=calendar)
    ncdstime.close()
    
    
    if type(plevel) == int:
        ncds = nc.Dataset(Mediadir + "ERA5_STARS/PL_centred_fields/"
                      +var + '_'+ str(plevel) +'_allObs'+interp+'.nc','w',format='NETCDF3_CLASSIC')
    elif plevel == '':#if problems with plevel= int this could be "else"
        ncds = nc.Dataset(Mediadir + "ERA5_STARS/PL_centred_fields/"
                      +var +'_allObs'+interp+'.nc','w',format='NETCDF3_CLASSIC')
    else:#for all levels
        ncds = nc.Dataset(Mediadir + "ERA5_STARS/PL_centred_fields/"
                      +var + '_'+ plevel +'_allObs'+interp+'.nc','w',format='NETCDF3_CLASSIC')            
        ncds.createDimension('plev',var_interp.shape[1])
        

        ncplev = ncds.createVariable('plev','f',('plev',))
        ncplev.setncattr('name', ds0.plev.standard_name)
        ncplev.setncattr('units','hPa')
        ncplev.setncattr('axis','Z')
        ncplev.setncattr('positive','down')        
        ncplev[:] = ds0.plev.values
        

    ncds.createDimension('time',var_interp.shape[0])        
    ncds.createDimension('x',var_interp.shape[-2])
    ncds.createDimension('y',var_interp.shape[-1])
    #ncds.createDimension('xy',field.shape[0]*field.shape[1])
#        ncds.createDimension('time',None)
    x = ncds.createVariable('x','f',('x',))
    x.setncattr('name','Distance in propagation direction')
    x.setncattr('units','km')
    x.setncattr('axis','X')
    x[:] = dist_grid_axis

    y = ncds.createVariable('y','f',('y',))
    y.setncattr('name','Distance perpendicular to propagation direction')
    y.setncattr('units','km')
    y.setncattr('axis','Y')
    y[:] = dist_grid_axis
    
    """add the datetime_vec also add the PL number"""
    nctime = ncds.createVariable('time','f',('time',))
    nctime.setncattr('name','time')
    nctime.setncattr('calendar', calendar)
    nctime.setncattr('axis','T')
    nctime.setncattr('units', timeunits)
    
    time_list= pd.to_datetime(time_list).values.astype('datetime64[h]')
    time_list= time_list.astype(datetime)
    time_list= nc.date2num(time_list, timeunits, calendar=calendar) #check why it does not work
    nctime[:]= time_list

    steer_vel = ncds.createVariable('steer_vel', 'f' ,('time',)) #some PLs have to be excluded.
    steer_vel.setncattr('name','Velocity of the steering wind [m/s]')
    steer_vel[:]= steer_vel_list

    ncPLnr = ncds.createVariable('PLnr', 'f' ,('time',)) #some PLs have to be excluded.
    ncPLnr.setncattr('name','Polar low number from Stoll list')
    ncPLnr[:]= PLnr_list_float
#    ncPLnr[:]= nc.stringtochar(np.array(PLnr_list, 'S'))

    ncObs = ncds.createVariable('Obs', 'f' ,('time',)) #some PLs have to be excluded.
    ncObs.setncattr('name','Observation number from Stoll list')
    ncObs[:]= Obs_list

    if type(plevel) == int or plevel == '': vari = ncds.createVariable(var,'f',('time','x','y',)) 
    else: vari = ncds.createVariable(var,'f',('time','plev', 'x','y',))
        
    if var != 'U':
        if plevel != '' and var != 'pv':
            vari.setncattr('standard_name', ds0[var].standard_name)
        vari.setncattr('long_name', ds0[var].long_name)
        vari.setncattr('units', ds0[var].units)
    if type(plevel) == int or plevel == '': vari[:,:,:] = var_interp

    else: vari[:,:,:,:] = var_interp
    

    ncds.setncattr('history',"Created by Patrick Stoll on %s." % \
                 (datetime.today().strftime("%Y-%m-%d") ) )
    ncds.close()        


