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

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'

Mediadir= '/run/media/pst019/1692A00D929FEF8B/'
homedir= Mediadir+'home/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from f_useful import *
from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs


#save= True
save= False

fignr= 5



S = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv")
S_ind= STARS_individual_systems(S)

fig = plt.figure(fignr, figsize= (8,6) )
fignr+=1
plt.clf()


ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (1,1,1))
#ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-15, 55, 50, 80], subplot= (1,1,1))

scale_bar(ax, 500, location= (0.06, 0.04))




"""get the ERA-5 data"""
import xarray as xr

plevel=None
vmin,vmax=None,None

"""specify field"""
#filetype='surface'
#var= '2t'
#var= '10u'
#var= '10v'
#var= 'U'
#var= 'skt'


filetype='plevels'
plevel= 925
#var='U'
#var= 'u'
#var='Vort'
#var= 'dudy'
#var= 'dvdx'
var='t'
#var='grad_t'
#vmin,vmax= 0, 10

#var= 'uvec_mean

#var='skt-t'
#var='sst-t'
#vmin,vmax= 35, 52
#
#
#filetype='vorticity'
#plevel= 850
#var='vo'


#filetype='boundary'
#var= 'cape'
#var= 'blh'
#var= 'tcc'
#var= 'hcc'

#filetype='forecast'
##var='sshf'
##var='slhf'
#var= 'tp'
##var= 'sf'
##var= 'lsp'
##var= 'cp'
##var= 'lsp'
#var='sf'
#var= 'ttr' #toa_outgoing_longwave_flux



sym= False
cmap='jet'
#cmap= 'viridis'
#cmap= 'Greys'
#cmap= 'Greys_r'
#cmap= 'Blues'
#cmap= 'Reds'

#sym= True
#cmap= 'RdBu_r'

black_radius= 200

"""specify time"""
##year, month, day= 2019, 5, 3 
##year, month, day= 2000, 3, 9
##year, month, day= 1999, 12, 20
#year, month, day= 2002, 1, 26
#
#hour= 4
#
#filetime= str(year)+'_'+str(month).zfill(2)+'_'+str(day).zfill(2)
#
#datetime_now= np.datetime64(str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'T'+str(hour).zfill(2) )
#tdelta= np.abs((S['datetime'] - datetime_now) / np.timedelta64(1, 'h'))
#S_now= S[tdelta== tdelta.min()]




"""specify system ID, Obsnr"""
ID= '81' #95
Obs= 1

S_now= S[np.logical_and(S['ID'] == ID, S['Obs'] == Obs ) ]

datetime_now= S_now.datetime.dt.round("H")

#S_now['dt_round']= [dt.replace(hour= 0, minute= 0) for dt in S_now['datetime']]

filetime= str(datetime_now.dt.year.values[0])+'_'+str(datetime_now.dt.month.values[0]).zfill(2)+'_'+str(datetime_now.dt.day.values[0]).zfill(2)
hour= datetime_now.dt.hour.values[0]

print(datetime_now)
print('Morphology now:', S_now['Morphology'].values[0])

"""general"""
if plevel==None:
    d0= xr.open_dataset(Mediadir + "ERA5_STARS/surface_era5_"+ filetime + '.nc')
    d0= d0.isel(time= hour)

else:
    d0= xr.open_dataset(Mediadir + "ERA5_STARS/vorticity_era5_"+ filetime + '.nc')
    d0= d0.isel(time= hour)

ds= xr.open_dataset(Mediadir + "ERA5_STARS/"+filetype+ "_era5_"+ filetime + '.nc')
ds= ds.isel(time= hour)

ds= xr.merge([d0, ds])

if filetype in ['plevels', 'vorticity']:
    ds['plev']/= 100
    ds= ds.sel(plev= plevel)



if filetype=='surface':
    if var== 'U':
        ds['U']= np.sqrt(ds['10u']**2 + ds['10v']**2)
        ds[var].attrs['units']= 'm/s'
else:
    if var== 'U':
        ds['U']= np.sqrt(ds['u']**2 + ds['v']**2)
        ds[var].attrs['units']= 'm/s'


if var=='vo':
    ds[var]*= 1E5
    ds[var].attrs['units']= '10$^{-5}$ 1/s'    

if var== 'Vort':
    dy= 0.25 *110E3 
    dx= 0.5 *110E3 * np.cos(np.deg2rad(ds.lat))
    dx= np.tile(dx, (len(ds.lon), 1)).T
    vort= (np.gradient(ds.v, axis= 0)/dx - np.gradient(ds.u, dy, axis= 1))*1E5
    
#    vort= (np.gradient(ds.v, dy, axis= 1)- np.gradient(ds.u, dx, axis= 0))*1E5

    ds[var]= (('lat', 'lon'), vort)
    ds[var].attrs['units']= '10**-5 s**-1' 


if var== 'dudy':
    dy= 0.25 *110E3 
    dx= 0.5 *110E3 * np.cos(np.deg2rad(ds.lat))
#    dx= np.tile(dx, (len(ds.lon), 1)).T
    vort= -1*(np.gradient(ds.u, dy, axis= 1))*1E5
    
    ds[var]= (('lat', 'lon'), vort)
    ds[var].attrs['units']= '10**-5 s**-1'     

if var== 'dvdx':
    dy= 0.25 *110E3 
    dx= 0.5 *110E3 * np.cos(np.deg2rad(ds.lat))
#    vort= (np.gradient(ds.v.values, dx.values, axis= 0))*1E5

    dx= np.tile(dx, (len(ds.lon), 1)).T
    vort= (np.gradient(ds.v.values, axis= 0))*1E5 /dx

    
#    grad= np.gradient(ds.v)[0]
#    vort= grad/dx *1E5

    
    ds[var]= (('lat', 'lon'), vort)
    ds[var].attrs['units']= '10**-5 s**-1'   



if var=='skt-t':
    d0= xr.open_dataset(Mediadir + "ERA5_STARS/surface_era5_"+ filetime + '.nc')
    d0= d0.isel(time= hour)

    ds= xr.merge([d0, ds])    
    
    ds[var]= ds['skt']- ds['t']
    ds[var].attrs['units']= 'K'   

if var=='sst-t':
    d0= xr.open_dataset(Mediadir + "ERA5_STARS/boundary_era5_"+ filetime + '.nc')
    d0= d0.isel(time= hour)

    ds= xr.merge([d0, ds])    
    
    ds[var]= ds['sst']- ds['t']
    ds[var].attrs['units']= 'K' 


if var== 'grad_t':
    dy= 0.25 *110E3 
    dx= 0.5 *110E3 * np.cos(np.deg2rad(ds.lat))
    dx= np.tile(dx, (len(ds.lon), 1)).T

    ds[var]= (('lat', 'lon'),  np.sqrt( (np.gradient(ds.t, axis= 0, edge_order=1)/dx)**2 + np.gradient(ds.t, dy, axis= 1, edge_order=1) **2) *1E5)
#    ds[var]= (('lat', 'lon'),  np.sqrt( (np.gradient(ds.t, axis= 0, edge_order=1)/dx)**2 ) *1E5)

    ds[var].attrs['units']= 'K/100 km'



if 'msl' in list(ds.keys()):
    ds['msl']/= 100
    ds['msl'].attrs['units']= 'hPa'

if var in ['ttr', 'sshf', 'slhf']:
    ds[var]/= -3600
    ds[var].attrs['units']= 'W m**-2'
if var in ['sf', 'lsp', 'cp', 'tp']:
    ds[var] *= 1000
    ds[var].attrs['units']= 'mm/h' 



if vmax != None:
    cf= ax.pcolormesh(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= vmin, vmax=vmax)
    
elif sym== False: cf= ax.pcolormesh(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap)
elif sym== True:
    vextr= np.max([np.max(ds[var]), -np.min(ds[var])])
    cf= ax.pcolormesh(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)


cb= fig.colorbar(cf, ax= ax, shrink=0.8)
cb.set_label(var + ' ['+ ds[var].units + ']', size=14)    


if plevel==None :
    cs= ax.contour(ds.lon, ds.lat, ds['msl'], np.arange(950, 1030, 2), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')

else:
    cs= ax.contour(ds.lon, ds.lat, ds['z']/9.81, np.arange(0, 10000, 10), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')


if var== 'U':
    #ax.quiver(ds.lon, ds.lat, ds['u'].values, ds['v'].values, pivot='middle', scale= 2000, transform= ccrs.PlateCarree())
    PlotWind(ds.lon, ds.lat, ds['u'].values, ds['v'].values, ax, nx= 15, ny=30, arrowtype= 'quiver', alen=30, scale=False)

#
#"""plot radius"""
ax.scatter(S_now.Longitude, S_now.Latitude, transform=ccrs.PlateCarree(), s= 10, color= 'r')

lon, lat= S_now.Longitude.values[0], S_now.Latitude.values[0]
radius= S_now.Diameter.values[0]/2

plot_circle(ax, lon, lat, radius, edgecolor= 'r')

plot_circle(ax, lon, lat, black_radius, edgecolor= 'grey')


plt.tight_layout()




"""get the shear"""
dist= 200

ds= xr.open_dataset(Mediadir + "ERA5_STARS/plevels_era5_"+ filetime + '.nc')
ds= ds.isel(time= hour)
ds1= xr.open_dataset(Mediadir + "ERA5_STARS/vorticity_era5_"+ filetime + '.nc')
ds1= ds1.isel(time= hour)            
ds= xr.merge([ds, ds1])

ds['plev']/= 100
ds= ds.sel(plev= [925, 700])



ds_lon= np.tile(ds.lon.values, (len(ds.lat), 1))
ds_lat= np.tile(ds.lat.values, (len(ds.lon), 1)).T  

dist_S_now= 110* np.sqrt( (S_now.Latitude.values[0]- ds_lat)**2+ (np.cos(np.deg2rad(S_now.Latitude.values[0]))* (S_now.Longitude.values[0]- ds_lon))**2)



u_mean= np.mean(ds['u'].mean(dim= 'plev').values[dist_S_now < dist])
v_mean= np.mean(ds['v'].mean(dim= 'plev').values[dist_S_now < dist])

print(u_mean, v_mean)


f= np.sin(np.deg2rad(S_now.Latitude.values[0])) * 2*np.pi/(24*60**2)


thickness= ds['z'].sel(plev=700) - ds['z'].sel(plev=925)

dy= 0.25 *110E3 
dx= 0.5 *110E3 * np.cos(np.deg2rad(ds.lat))
dx= np.tile(dx, (len(ds.lon), 1)).T
u_therm, v_therm= - np.gradient(thickness, dy, axis= 1), np.gradient(thickness, axis= 0)/dx

u_therm_m= 1/f* np.mean(u_therm[dist_S_now < dist])
v_therm_m= 1/f* np.mean(v_therm[dist_S_now < dist])

print(u_therm_m, v_therm_m)


alpha= angle_between((u_therm_m, v_therm_m), (u_mean, v_mean))

print(alpha)


#fig = plt.figure(fignr, figsize= (8,6) )
#fignr+=1
#plt.clf()
#
#ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (1,1,1))

cs= ax.contour(ds.lon, ds.lat, thickness/10, np.arange(0, 5000, 10), transform= ccrs.PlateCarree(), colors='r', linewidths= 1)
plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')

#
#
#data= ds[var].values
#dist= 300
#
#ds_lon= np.tile(ds.lon.values, (len(ds.lat), 1))
#ds_lat= np.tile(ds.lat.values, (len(ds.lon), 1)).T
#
#
##ds['dist']= (dims=('lon', 'lat'), data= 110* np.sqrt( (S_now.Latitude.values[0]- ds_lat)**2+ (np.cos(np.deg2rad(S_now.Latitude.values[0]))* (S_now.Longitude.values[0]- ds_lon))**2)  )
#
#dist_S_now= 110* np.sqrt( (S_now.Latitude.values[0]- ds_lat)**2+ (np.cos(np.deg2rad(S_now.Latitude.values[0]))* (S_now.Longitude.values[0]- ds_lon))**2)
#
#localmax= np.max(ds[var].values[dist_S_now < dist])
#localmin= np.min(ds[var].values[dist_S_now < dist])
#
#print(localmax, localmin)
#
##import scipy.ndimage.filters as filters
##newmax= filters.maximum_filter(data, dist)
#
##def LocalMax(ds[var].values, )
#
#
#"""scatter plot of systems at this time"""
#
###ds[var].plot.contourf(ax= ax, transform= ccrs.PlateCarree())
#ds[var].plot(ax= ax, transform= ccrs.PlateCarree())
###
##
##
##fig = plt.figure(fignr)
##fignr+=1
##plt.clf()
##
##
##ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (1,1,1))
###newmax.plot(ax= ax, transform= ccrs.PlateCarree())
##



