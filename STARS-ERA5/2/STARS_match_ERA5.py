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
from f_STARS import *



#save= True
save= False

fignr= 1


S = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv")



S['dt_round']= [dt.round("H") for dt in S['datetime']]


write=False
write=True

gausfilter=None

plevel=None


#filetype='surface'
#var= '2t'
#var= '10u'
#var= '10v'
#var= 'U'
#var= 'skt'
#var= 'msl'


#filetype='boundary'
#var= 'cape'
#var= 'blh'


#filetype='forecast'
#var='sshf'
#var='slhf'
#var= 'tp'
#var='cp'
#var= 'sf'
##var= 'lsp'
##var= 'ttr' #toa_outgoing_longwave_flux


filetype='plevels'
plevel= 850
#var='grad_t'
#var='barotropic'
var, gausfilter='barotropic_efolding_filter', 10

#
#
#filetype='vorticity'
#plevel= 850
#var='vo'

#filetype=['surface','plevels']
#var= 'skt-t'
#plevel= 700


#filetype=['plevels', 'vorticity']
#plevel= [925, 700]
#var='vert_shear_angle'
#var, gausfilter= 'baroclinic_efolding_filter', 4



dist= 200
#treat= 'max'
treat='min'
#treat='mean'
#treat='med'


fieldname= var
if gausfilter: fieldname += str(gausfilter) +'_'
if plevel == None: fieldname += '_' +treat +'_'+ str(dist)
elif type(plevel) == list: fieldname += str(plevel[0])+'-'+str(plevel[1])+ '_' +treat +'_'+ str(dist)
else: fieldname += str(plevel)+ '_' +treat +'_'+ str(dist)


import xarray as xr
import os.path
csv_name= Mediadir + "PL/STARS/STARS_ERA5.csv"
exist= os.path.isfile(csv_name) 

print(exist)

import csv



data_array= []


for index in range(len(S)):

    S_now = S.iloc[index]
    
    if np.isnan(S_now.lat):
        data_array += [np.nan]
    
    else:
        dtnow= S_now.dt_round
        year, month, day, hour= dtnow.year, dtnow.month, dtnow.day, dtnow.hour
        
        """get the ERA-5 data"""   
        filetime= str(year)+'_'+str(month).zfill(2)+'_'+str(day).zfill(2)
        
        if type(filetype) == str:
            ds= xr.open_dataset(Mediadir + "ERA5_STARS/"+filetype+ "_era5_"+ filetime + '.nc')
            ds= ds.isel(time= hour)
            if 'plev' in list(ds.coords.keys()):
                ds['plev']/= 100
                ds= ds.sel(plev= plevel)            
            
        else: #if two files have to be opened
            ds= xr.open_dataset(Mediadir + "ERA5_STARS/"+filetype[0]+ "_era5_"+ filetime + '.nc')
            ds= ds.isel(time= hour)
            ds1= xr.open_dataset(Mediadir + "ERA5_STARS/"+filetype[1]+ "_era5_"+ filetime + '.nc')
            ds1= ds1.isel(time= hour)            
            ds= xr.merge([ds, ds1])
            if 'plev' in list(ds.coords.keys()):
                ds['plev']/= 100
                ds= ds.sel(plev= plevel)
        
        if var== 'U': ds['U']= np.sqrt(ds['10u']**2 + ds['10v']**2)
        if var in ['ttr', 'sshf', 'slhf']:
            ds[var]/= -3600
            ds[var].attrs['units']= 'W m**-2'
        if var in ['sf', 'lsp', 'cp', 'tp']:
            ds[var] *= 1000
            ds[var].attrs['units']= 'mm/h' 

        if var== 'skt-t':
            ds[var]= ds['skt']- ds['t']
            
        if var=='grad_t':

            variable= np.sqrt(grad_x(ds.t, ds.lat, ds.lon)**2 + grad_y(ds.t, ds.lat, ds.lon)**2) *1E5
            ds[var]= (('lat', 'lon'),  variable)                    
            

        if var=='vert_shear_angle':
            z_diff= (ds['z'].sel(plev=700) - ds['z'].sel(plev=925) )
            f= np.sin(np.deg2rad(ds.lat)) * 2*2*np.pi/(24*60**2)
            f= np.tile(f, (len(ds.lon), 1)).T
    
            u_therm, v_therm= -1/f* grad_y(z_diff, ds.lat, ds.lon), 1/f*grad_x(z_diff, ds.lat, ds.lon)
                
            u_mean= value_in_rad(ds['u'].mean(dim='plev'), ds.lat, ds.lon, S_now.lat, S_now.lon, distance= dist)
            v_mean= value_in_rad(ds['v'].mean(dim='plev'), ds.lat, ds.lon, S_now.lat, S_now.lon, distance= dist)
            print('mean wind:', u_mean, v_mean)
                  
            u_therm_m= value_in_rad(u_therm, ds.lat, ds.lon, S_now.lat, S_now.lon, distance= dist, intype='np')
            v_therm_m= value_in_rad(v_therm, ds.lat, ds.lon, S_now.lat, S_now.lon, distance= dist, intype='np')    
            print('thermal wind:', u_therm_m, v_therm_m)
                
            alpha= angle_between((u_therm_m, v_therm_m), (u_mean, v_mean))    

            print(index, alpha)
        
            data_array += [alpha]
            

        if var=='barotropic':
            #calculates the barotropic growth rate by 0.2* dU/dy, to get the e-folding time: 1/sigma * 1/60**2
            variable= 0.2 * np.sqrt(grad_y(ds.u, ds.lat, ds.lon)**2 + grad_x(ds.v, ds.lat, ds.lon)**2)
            ds[var]= (('lat', 'lon'),  variable)  


        if var=='barotropic_efolding_filter':
            #calculates the barotropic growth rate by 0.2* dU/dy, to get the e-folding time: 1/sigma * 1/60**2
            ds['u'] = (('lat', 'lon'), filters.gaussian_filter(ds.u, sigma= gausfilter, mode='nearest', truncate= 1.) )
            ds['v'] = (('lat', 'lon'), filters.gaussian_filter(ds.v, sigma= gausfilter, mode='nearest', truncate= 1.) )
        
            variable= 0.2* np.sqrt(grad_y(ds.u, ds.lat, ds.lon)**2 + grad_x(ds.v, ds.lat, ds.lon)**2)
            variable= 1/variable* 1/(60**2)   

            ds[var]= (('lat', 'lon'),  variable)  


        if var=='baroclinic_efolding_filter':
            ds['u'] = (('plev', 'lat', 'lon'), filters.gaussian_filter(ds.u, sigma= (0, gausfilter, gausfilter) , mode='nearest', truncate= 1.) )
            ds['v'] = (('plev', 'lat', 'lon'), filters.gaussian_filter(ds.v, sigma= (0, gausfilter, gausfilter) , mode='nearest', truncate= 1.) )
            ds['z'] = (('plev', 'lat', 'lon'), filters.gaussian_filter(ds.z, sigma= (0, gausfilter, gausfilter) , mode='nearest', truncate= 1.) )
            ds['t'] = (('plev', 'lat', 'lon'), filters.gaussian_filter(ds.t, sigma= (0, gausfilter, gausfilter) , mode='nearest', truncate= 1.) )
        
            U_diff= np.sqrt(ds.u.isel(plev=1)**2 + ds.v.isel(plev=1)**2) - np.sqrt(ds.u.isel(plev=0)**2 + ds.v.isel(plev=0)**2)
            h_diff= (ds['z'].isel(plev=1) - ds['z'].isel(plev=0) )/9.81
            dUdH= np.abs(U_diff /h_diff)
            
            f= np.sin(np.deg2rad(ds.lat)) * 2*2*np.pi/(24*60**2)
            f= np.tile(f, (len(ds.lon), 1)).T 
            
            theta_1, theta_0= PotTemp(ds.t.isel(plev=1), plevel[1]), PotTemp(ds.t.isel(plev=0), plevel[0])
            N= np.sqrt(9.81/(theta_1+theta_0)/2 *  (theta_1- theta_0)/h_diff)
        
            variable= 0.31*f/N* dUdH    
            variable= 1/variable* 1/(60**2)   
            ds[var]= (('lat', 'lon'),  variable)



        """get the local max/min/mean/median"""
        if var != 'vert_shear_angle':
            data= ds[var].values
            
            ds_lon= np.tile(ds.lon.values, (len(ds.lat), 1))
            ds_lat= np.tile(ds.lat.values, (len(ds.lon), 1)).T  
            
            dist_S_now= 110* np.sqrt( (S_now.lat- ds_lat)**2+ (np.cos(np.deg2rad(S_now.lat))* (S_now.lon- ds_lon))**2)
            
            if treat =='max': local= np.max(ds[var].values[dist_S_now < dist])
            elif treat == 'min': local= np.min(ds[var].values[dist_S_now < dist])
            elif treat == 'mean': local= np.mean(ds[var].values[dist_S_now < dist])
            elif treat == 'med': local= np.median(ds[var].values[dist_S_now < dist])
            
            print(index, local)
        
            data_array += [local]
        

print(fieldname)
    
if len(data_array)== 3850: 
    if write:
        print('write')
        if exist == False:  #create new file
            with open(csv_name,'w') as file:
                writer = csv.DictWriter(file, delimiter=',', fieldnames= [fieldname])
                writer.writeheader()
                for index in range(len(data_array)):
                    writer.writerow({fieldname: str(data_array[index])})
        
        if exist: #write into existing file
            csv_input = pd.read_csv(csv_name)
            csv_input[fieldname] = data_array
            csv_input.to_csv(csv_name, index=False)
    






