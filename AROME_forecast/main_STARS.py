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
from f_meteo import *

#save= True
save= False

fignr= 8


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
setup= None

"""specify field"""
#filetype='surface'
#var= '2t'
#var= '10u'
#var= '10v'
#var= 'U'
#var= 'skt'


filetype='plevels'
plevel= 850
#var='U'
#var= 'u'
#var='Vort'
#var='t'
#var='grad_t'
#vmin,vmax= 0, 10

var='dU/dL'

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

#filetype='plevels'
#plevel= [925, 700]
#var='dU/dH'
#var='baroclinic'

"""setup"""
#filetype= 'plevels'
##setup, plevel, var, cont_var= 'geo', 850, 'U', 'z'
#setup, plevel, var, cont_var= 'shear', [925, 700], 'thickness', 'z'


"""the contour variable"""
if plevel== None: cont_var= 'msl'
cont_var='z'
#cont_var='t'


sym= False
#cmap='jet'
#cmap= 'viridis'
#cmap= 'Greys'
#cmap= 'Greys_r'
#cmap= 'Blues'
cmap= 'Reds'

#sym= True
#cmap= 'RdBu_r'

distance= 200





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
ID= '12' #95
Obs= 4

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




"""--------------create/prepare the plotted variable-------------------"""
if var== 'U':
    if filetype=='surface':
        ds['U']= np.sqrt(ds['10u']**2 + ds['10v']**2)
        ds[var].attrs['units']= 'm/s'
    else:
        ds['U']= np.sqrt(ds['u']**2 + ds['v']**2)
        ds[var].attrs['units']= 'm/s'

if var=='vo':
    ds[var]*= 1E5
    ds[var].attrs['units']= '10$^{-5}$ 1/s'    

if var== 'Vort': 
    vort= grad_x(ds.v, ds.lat, ds.lon) - grad_y(ds.u, ds.lat, ds.lon)
    vort *= 1E5

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
    variable= np.sqrt(grad_x(ds.t, ds.lat, ds.lon)**2 + grad_y(ds.t, ds.lat, ds.lon)**2) *1E5
    ds[var]= (('lat', 'lon'),  variable)

    ds[var].attrs['units']= 'K/100 km'

if var=='thickness':   
    z_diff= (ds['z'].sel(plev=700) - ds['z'].sel(plev=925) )
    ds[var]= (('lat', 'lon'),  z_diff/9.81)
    ds[var].attrs['units']= 'm'

if var== 'dU/dL':
    variable= np.sqrt(grad_y(ds.u, ds.lat, ds.lon)**2 + grad_x(ds.v, ds.lat, ds.lon)**2)
#    variable= np.sqrt(grad_x(ds.v, ds.lat, ds.lon)**2)
#    variable= np.sqrt(grad_y(ds.u, ds.lat, ds.lon)**2)
    ds[var]= (('lat', 'lon'),  variable)
    ds[var].attrs['units']= '1/s'


if var in ['dU/dH', 'baroclinic']:
    U_diff= np.sqrt(ds.u.isel(plev=1)**2 + ds.v.isel(plev=1)**2) - np.sqrt(ds.u.isel(plev=0)**2 + ds.v.isel(plev=0)**2)
    h_diff= (ds['z'].isel(plev=1) - ds['z'].isel(plev=0) )/9.81
    dUdH= np.abs(U_diff /h_diff)

    ds[var]= (('lat', 'lon'),  dUdH)
    ds[var].attrs['units']= '1/s'    


if var=='baroclinic':
    f= np.sin(np.deg2rad(ds.lat)) * 2*2*np.pi/(24*60**2)
    f= np.tile(f, (len(ds.lon), 1)).T
    
    theta_1, theta_0= PotTemp(ds.t.isel(plev=1), plevel[1]), PotTemp(ds.t.isel(plev=0), plevel[0])
    N= np.sqrt(9.81/(theta_1+theta_0)/2 *  (theta_1- theta_0)/h_diff)

    variable= 0.31*f/N* dUdH
    ds[var]= (('lat', 'lon'),  variable)
    ds[var].attrs['units']= '1/s'    
    
if 'msl' in list(ds.keys()):
    ds['msl']/= 100
    ds['msl'].attrs['units']= 'hPa'

if var in ['ttr', 'sshf', 'slhf']:
    ds[var]/= -3600
    ds[var].attrs['units']= 'W m**-2'
if var in ['sf', 'lsp', 'cp', 'tp']:
    ds[var] *= 1000
    ds[var].attrs['units']= 'mm/h' 




"""make the shading plot"""
if setup != 'shear':
    if vmax != None:
        cf= ax.pcolormesh(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= vmin, vmax=vmax)
        
    elif sym== False: cf= ax.pcolormesh(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap)
    elif sym== True:
        vextr= np.max([np.max(ds[var]), -np.min(ds[var])])
        cf= ax.pcolormesh(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
    
    
    cb= fig.colorbar(cf, ax= ax, shrink=0.8)
    varlabel= var
    if plevel != None: varlabel = var+'_'+str(plevel)
    cb.set_label(varlabel + ' ['+ ds[var].units + ']', size=14)    




"""make the contour plot"""
if cont_var== 'msl' :
    cs= ax.contour(ds.lon, ds.lat, ds['msl'], np.arange(950, 1030, 2), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')

elif cont_var== 'z':
    if len(ds['z'].shape)== 3: plotz= ds['z'].mean(dim='plev') #if there are several pressure levels, the mean is taken
    else: plotz= ds['z']
    
    cs= ax.contour(ds.lon, ds.lat, plotz/9.81, np.arange(0, 10000, 10), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')

elif cont_var=='t':
    cs= ax.contour(ds.lon, ds.lat, ds['t'], np.arange(200, 400, 1), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')



if var in ['U', 'dU/dL']: #plot the wind vector
    nx= 10
#    ax.quiver(ds.lon[::nx], ds.lat[::nx], ds['u'].values[::nx, ::nx], ds['v'].values[::nx, ::nx], pivot='middle', scale= 700, transform= ccrs.PlateCarree())
    PlotWind(ds.lon, ds.lat, ds['u'].values, ds['v'].values, ax, nx= 15, ny=30, arrowtype= 'quiver', alen=15, scale=False)






if setup=='geo': #geostrophic wind

    f= np.sin(np.deg2rad(ds.lat)) * 2*2*np.pi/(24*60**2)
    f= np.tile(f, (len(ds.lon), 1)).T
    
   
    u_g= -1/f* grad_y(ds.z, ds.lat, ds.lon)
    v_g= 1/f* grad_x(ds.z, ds.lat, ds.lon)
    
#    ax.quiver(ds.lon[::nx]+0.3, ds.lat[::nx]+.3, u_g[::nx, ::nx], v_g[::nx, ::nx], pivot='middle', color= 'r', scale= 700, transform= ccrs.PlateCarree())
    PlotWind(ds.lon+0.3, ds.lat+.3, u_g, v_g, ax, nx= 15, ny=30, arrowtype= 'quiver', alen=15, scale=False, color= 'r')



"""get the shear"""
if setup== 'shear':
    PlotWind(ds.lon, ds.lat, ds['u'].mean(dim='plev').values, ds['v'].mean(dim='plev').values, ax, nx= 15, ny=30, arrowtype= 'quiver', alen=15, scale=False)
    
    #calculate and plot the thermal wind
    f= np.sin(np.deg2rad(ds.lat)) * 2*2*np.pi/(24*60**2)
    f= np.tile(f, (len(ds.lon), 1)).T
    
    u_therm, v_therm= -1/f* grad_y(z_diff, ds.lat, ds.lon), 1/f*grad_x(z_diff, ds.lat, ds.lon)

    cs= ax.contour(ds.lon, ds.lat, z_diff/9.81, np.arange(0, 5000, 10), transform= ccrs.PlateCarree(), colors='r', linewidths= 1)
    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')
    
    PlotWind(ds.lon[10:], ds.lat[10:], u_therm[10:, 10:], v_therm[10:, 10:], ax, nx= 15, ny=30, arrowtype= 'quiver', alen=15, scale=False, color= 'r')
    
    
    #calculate the shear        
    u_mean= value_in_rad(ds['u'].mean(dim='plev'), ds.lat, ds.lon, S_now.Latitude.values[0], S_now.Longitude.values[0], distance= distance)
    v_mean= value_in_rad(ds['v'].mean(dim='plev'), ds.lat, ds.lon, S_now.Latitude.values[0], S_now.Longitude.values[0], distance= distance)
    print('mean wind:', u_mean, v_mean)
          
    u_therm_m= value_in_rad(u_therm, ds.lat, ds.lon, S_now.Latitude.values[0], S_now.Longitude.values[0], distance= distance, intype='np')
    v_therm_m= value_in_rad(v_therm, ds.lat, ds.lon, S_now.Latitude.values[0], S_now.Longitude.values[0], distance= distance, intype='np')    
    print('thermal wind:', u_therm_m, v_therm_m)
        
    alpha= angle_between((u_therm_m, v_therm_m), (u_mean, v_mean))    
    print('alpha:', alpha)



        
    




"""plot radius"""
ax.scatter(S_now.Longitude, S_now.Latitude, transform=ccrs.PlateCarree(), s= 10, color= 'g')

lon, lat= S_now.Longitude.values[0], S_now.Latitude.values[0]
radius= S_now.Diameter.values[0]/2

plot_circle(ax, lon, lat, radius, edgecolor= 'r')

plot_circle(ax, lon, lat, distance, edgecolor= 'green')


plt.tight_layout()




"""print local max, min, mean"""
localmax= value_in_rad(ds[var], ds.lat, ds.lon, S_now.Latitude.values[0], S_now.Longitude.values[0], distance= distance, type='max')
localmin= value_in_rad(ds[var], ds.lat, ds.lon, S_now.Latitude.values[0], S_now.Longitude.values[0], distance= distance, type='min')
localmean= value_in_rad(ds[var], ds.lat, ds.lon, S_now.Latitude.values[0], S_now.Longitude.values[0], distance= distance, type='mean')

print('localmax, localmin, localmean: ', localmax, localmin, localmean)



if var== 'dU/dL':
    barotropic= 0.2* localmax
    print('barotropic growth rate: ', np.round(barotropic, 6) )
    print('e-folding time [h]: ', np.round(1/barotropic * 1/(60**2), 1) )


if var== 'baroclinic':
    print('baroclinic growth rate: ', np.round(barotropic, 6) )
    print('e-folding time [h]: ', np.round(1/barotropic * 1/(60**2), 1) )
