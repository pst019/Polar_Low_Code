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
    Mediadir= '/run/media/pst019/PatsOrange/'

homedir=Mediadir+'home/'


import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import csaps #pip install csaps==0.10.0 had to add "np.pad(u, pad_width, mode='constant'), axis=0) " in _sspumv.py


from f_useful import *
from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs
from f_meteo import *
import scipy.ndimage.filters as filters
import great_circle_calculator.great_circle_calculator as gcc # pip install great-circle-calculator


save= True
#save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/Track_fields/'
save_extra= ''

fignr= 3


#track_type='Rojo'
track_type='Stoll'

imp_lists= True
#imp_lists= False

ID= '10_1' #20_1 (Obs20 for figure)'25_1' #'95' #'12'  #95
Obs= 15 #starts with 1

Plot_PLCG= True
if Plot_PLCG: save_extra+='_PLCG-smth'

Plot_smth_track= True
if Plot_smth_track: save_extra+='_track-smth'

Plot_PL_track=True
if Plot_PL_track == False: save_extra += '_noPLtrack'

Plot_circle=False
radius= 500 #radius of the green circle

Plot_all_my_tracks=False
if Plot_all_my_tracks: save_extra += '_allTRACKs'

Plot_Local_Vort_Max=False #only relevant if Plot_all_my_tracks is False
if Plot_Local_Vort_Max and Plot_all_my_tracks== False: save_extra += '_Vo-max'

Plot_Rojo_track=True
if Plot_Rojo_track: save_extra += '_Rojo'

legend=False
legend=True



"""get the ERA-5 data"""

plevel=None
clevels=''
vmin,vmax=None,None
setup= None #'geo'


uniform_filter=False
filtercells= 10
if uniform_filter: save_extra+= '_filter'+str(filtercells)
"""specify field"""
#filetype='surface'
#var= '2t'
#var= '10u'
#var= '10v'
#var= 'U'
#var= 'skt'


filetype='tracking'
plevel= 850
#var='U'
#var= 'u'
var, clevels='vo', [-3, -2, -1, -.5, 0.5, 1, 2, 3]
#var, clevels='t', np.arange(250, 270, 1)
#var, gausfilter='grad_t', 1
#var= 'z'
#vmin,vmax= 0, 10

#var, gausfilter='barotropic_efolding_filter', 4
#var, gausfilter='barotropic_gr_filter', 4

#var='skt-t'
#var='sst-t'
#plevel=500
#vmin,vmax= 35, 52
#
#


#filetype='plev_w'
#plevel= 925
#var='w'
#var, clevels= 'd', np.arange(-10, 10.1, 2)*1E-5
#var, clevels= 'd_calc', np.arange(-10, 10.1, 2)*1E-5 #filetype = 'tracking' would be sufficient, but is loaded anyway to get z

#var='pv'

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


#filetype='tracking'
#plevel= [1000, 850]
#var, gausfilter= 'N', 4
#var, gausfilter='baroclinic_dUdz_efolding_filter', 4
#var, gausfilter='baroclinic_dUdz_gr_filter', 4 #calculated by f/N dU/dz

#plevel= [1000, 925, 850]
#plevel= [925, 850, 500]
#var, gausfilter='baroclinic_dTdy_gr_filter', 10 #calculated by  g/NT dT/dy, 3 plevels: the first and last are used for N and the middle one for T



#vmin,vmax= 0, 30

"""setup"""
#filetype= 'tracking'
##setup, plevel, var, cont_var= 'geo', 850, 'U', 'z'
#setup, plevel, var, cont_var= 'shear', [925, 700], 'thickness', 'z'
#it prints

"""the contour variable"""
if plevel== None: cont_var= 'msl'
cont_var='z'
#cont_var='t'
if type(plevel) == int: plevel_contvar= plevel
else:
    plevel_contvar= 850
    print('Specify plevel of the contour variable to: ', plevel_contvar)

numbers_cont = True #specifies if there are numbers on the contours

sym= False
#cmap='jet'
#cmap= 'viridis'
#cmap= 'Greys'
#cmap= 'Greys_r'
#cmap= 'Blues'
cmap= 'Reds'
#cmap= 'Reds_r'

sym= True
cmap= 'RdBu_r'






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




if imp_lists:
    if track_type=='Rojo':
        S = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv")
        S_ind= STARS_individual_systems(S)
    
    elif track_type== 'Stoll':
        test= 'version4'
        dist= 150
        
        Stoll_STARS_file="ERA5_STARS/STARS-PMC-merge/"+test+"_dist_"+str(dist)+"_handfix.csv"
        Stoll = pd.read_csv(Mediadir+ Stoll_STARS_file)
        Stoll['time']= pd.DatetimeIndex(Stoll.time)
        
        Stoll= Stoll_excl(Stoll)
        Stoll= Stoll_Obs_nr(Stoll) #get the Obs nr
    
        Stoll= Stoll.rename(columns={"Stoll nr": "ID"})
        
        
        """Stoll_individual_systems"""
        Stoll_ind= Stoll_individual_systems(Stoll, ID_name= 'ID')    
    
        S= Stoll
        S_ind= Stoll_ind



fig = plt.figure(fignr, figsize= (10,6) )
fignr+=1
plt.clf()


#ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (1,1,1))
ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-10, 25, 68, 80], subplot= (1,1,1))
#ax= Plot_Polar_Stereo(fig, central_longitude= 20, extent= [15, 35, 72, 78], subplot= (1,1,1))
#ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-40, 65, 50, 85], subplot= (1,1,1))


scale_bar(ax, 500, location= (0.06, 0.04))




"""specify system ID, Obsnr"""
S_now= S[np.logical_and(S['ID'] == ID, S['Obs'] == Obs ) ]

datetime_now= S_now.time.dt.round("H")

#S_now['dt_round']= [dt.replace(hour= 0, minute= 0) for dt in S_now['datetime']]

filetime= str(datetime_now.dt.year.values[0])+'_'+str(datetime_now.dt.month.values[0]).zfill(2)+'_'+str(datetime_now.dt.day.values[0]).zfill(2)
hour= datetime_now.dt.hour.values[0]

print(datetime_now)
print('Morphology now:', S_now['Morphology'].values[0])

"""general"""

#else:
#    d0= xr.open_dataset(Mediadir + "ERA5_STARS/data/tracking_era5_"+ filetime + '.nc')
#    d0= d0.isel(time= hour)

ds= xr.open_dataset(Mediadir + "ERA5_STARS/data/"+filetype+ "_era5_"+ filetime + '.nc')
ds= ds.isel(time= hour)

if plevel==None: #to get the surface pressure
    d0= xr.open_dataset(Mediadir + "ERA5_STARS/data/surface_era5_"+ filetime + '.nc')
    d0= d0.isel(time= hour)

    ds= xr.merge([d0, ds])

if filetype == 'plev_w': #to get the geopotential for the contours
    d0= xr.open_dataset(Mediadir + "ERA5_STARS/data/tracking_era5_"+ filetime + '.nc')
    d0= d0.isel(time= hour)
    ds= xr.merge([d0, ds])
    
if filetype in ['plev_w', 'tracking']:
    ds['plev']/= 100
    dsorg= ds
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
    ds[var]*= 1E4
    ds[var].attrs['units']= '10$^{-4}$ s$^{-1}$'    

if var== 'Vort': 
    vort= grad_x(ds.v, ds.lat, ds.lon) - grad_y(ds.u, ds.lat, ds.lon)
    vort *= 1E5

    ds[var]= (('lat', 'lon'), vort)
    ds[var].attrs['units']= '10**-5 s**-1' 

if var=='skt-t':
    d0= xr.open_dataset(Mediadir + "ERA5_STARS/data/surface_era5_"+ filetime + '.nc')
    d0= d0.isel(time= hour)
    ds= xr.merge([d0, ds])    
    
    ds[var]= ds['skt']- ds['t']
    ds[var].attrs['units']= 'K'   

if var=='sst-t':
    d0= xr.open_dataset(Mediadir + "ERA5_STARS/data/boundary_era5_"+ filetime + '.nc')
    d0= d0.isel(time= hour)
    ds= xr.merge([d0, ds])    
    
    ds[var]= ds['sst']- ds['t']
    ds[var]= ds[var].fillna(value= 0) #fill nan values with 0
    
    ds[var].attrs['units']= 'K' 


if var== 'grad_t':
    ds['t'] = (('lat', 'lon'), filters.gaussian_filter(ds.t, sigma= gausfilter, mode='nearest', truncate= 1.) )

    variable= np.sqrt(grad_x(ds.t, ds.lat, ds.lon)**2 + grad_y(ds.t, ds.lat, ds.lon)**2) *1E5
    ds[var]= (('lat', 'lon'),  variable)

    ds[var].attrs['long_name']= 'Horizontal temperature gradient'
    ds[var].attrs['units']= 'K/100 km'

if var=='thickness':   
    z_diff= (ds['z'].isel(plev=1) - ds['z'].isel(plev=0) )
    ds[var]= (('lat', 'lon'),  z_diff/9.81)
    ds[var].attrs['units']= 'm'

if var== 'dU/dL':
    variable= np.sqrt(grad_y(ds.u, ds.lat, ds.lon)**2 + grad_x(ds.v, ds.lat, ds.lon)**2)
    ds[var]= (('lat', 'lon'),  variable)
    ds[var].attrs['units']= '1/s'


if 'barotropic' in var:
    ds['u'] = (('lat', 'lon'), filters.gaussian_filter(ds.u, sigma= gausfilter, mode='nearest', truncate= 1.) )
    ds['v'] = (('lat', 'lon'), filters.gaussian_filter(ds.v, sigma= gausfilter, mode='nearest', truncate= 1.) )
    variable= 0.2* np.sqrt(grad_y(ds.u, ds.lat, ds.lon)**2 + grad_x(ds.v, ds.lat, ds.lon)**2)    

#if var=='barotropic_efolding_filter':
#    variable= 1/variable* 1/(60**2)   
#    ds[var]= (('lat', 'lon'),  variable)
#    ds[var].attrs['units']= 'h'
#
#if var=='barotropic_gr_filter':
#    variable*=60**2*24
#    ds[var]= (('lat', 'lon'),  variable)
#    ds[var].attrs['units']= '1/day'
    
    
if 'baroclinic_dUdz' in var:
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

if 'baroclinic_dTdy' in var: #calculate by g/NT dT/dy, 3 plevels: the first and last are used for N and the middle one for T
    ds['z'] = (('plev', 'lat', 'lon'), filters.gaussian_filter(ds.z, sigma= (0, gausfilter, gausfilter) , mode='nearest', truncate= 1.) )
    ds['t'] = (('plev', 'lat', 'lon'), filters.gaussian_filter(ds.t, sigma= (0, gausfilter, gausfilter) , mode='nearest', truncate= 1.) )

    g= 9.81
    theta_2, theta_0= PotTemp(ds.t.isel(plev=2), plevel[2]), PotTemp(ds.t.isel(plev=0), plevel[0])
    h_diff= (ds['z'].isel(plev=2) - ds['z'].isel(plev=0) )/g

    N= np.sqrt(9.81/(theta_2+theta_0)/2 *  (theta_2- theta_0)/h_diff)

    dTdy= np.sqrt(grad_x(ds['t'].isel(plev=1), ds.lat, ds.lon)**2 + grad_y(ds['t'].isel(plev=1), ds.lat, ds.lon)**2)
    T= ds['t'].isel(plev=1)
    
    variable=0.31 * g/(N*T) * dTdy


if '_gr_' in var:
    variable*=60**2    
    ds[var]= (('lat', 'lon'),  variable)
    ds[var].attrs['units']= '1/h'     
    
if 'efolding' in var:
    variable= 1/variable* 1/(60**2)   
    ds[var]= (('lat', 'lon'),  variable)
    ds[var].attrs['units']= 'h'    
      
    
    
if var=='N':
    ds['z'] = (('plev', 'lat', 'lon'), filters.gaussian_filter(ds.z, sigma= (0, gausfilter, gausfilter) , mode='nearest', truncate= 1.) )
    ds['t'] = (('plev', 'lat', 'lon'), filters.gaussian_filter(ds.t, sigma= (0, gausfilter, gausfilter) , mode='nearest', truncate= 1.) )

    h_diff= (ds['z'].isel(plev=1) - ds['z'].isel(plev=0) )/9.81
    theta_1, theta_0= PotTemp(ds.t.isel(plev=1), plevel[1]), PotTemp(ds.t.isel(plev=0), plevel[0])
    N= np.sqrt(9.81/(theta_1+theta_0)/2 *  (theta_1- theta_0)/h_diff)        
    ds[var]= (('lat', 'lon'),  N)
    ds[var].attrs['units']= '1/s'         
    


if var=='d_calc':
    grad_u_vec= np.gradient(ds.u, 25E3) #two components of the gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km
    grad_v_vec= np.gradient(ds.v, 25E3) #two components of the gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km
    
    ds[var]= (('lat', 'lon'), (grad_x(ds.u, ds.lat, ds.lon) + grad_y(ds.v, ds.lat, ds.lon) ) )
#    ds[var]*= 1E5
    ds[var].attrs['units'] = '(10$^{-5}$s$^{-1}$'
    ds[var].attrs['long_name']='Div (self calc)' 

    
if 'msl' in list(ds.keys()):
    ds['msl']/= 100
    ds['msl'].attrs['units']= 'hPa'

if var in ['ttr', 'sshf', 'slhf']:
    ds[var]/= -3600
    ds[var].attrs['units']= 'W m**-2'
if var in ['sf', 'lsp', 'cp', 'tp']:
    ds[var] *= 1000
    ds[var].attrs['units']= 'mm/h' 



if uniform_filter:
    ds[var].values = filters.uniform_filter(ds[var], size= filtercells, mode='nearest')


"""-----------------------------------------------------"""
"""make the shading plot"""
if setup != 'shear':
    if clevels != '':
        cf= ax.contourf(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, levels= clevels, vmin= clevels[0], vmax=clevels[-1], extend='both')            
    elif vmax != None:
#        cf= ax.pcolormesh(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= vmin, vmax=vmax)
        cf= ax.contourf(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= vmin, vmax=vmax)
        
    elif sym== False:
#        cf= ax.pcolormesh(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap)
        cf= ax.contourf(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap)
    elif sym== True:
        vextr= np.max([np.max(ds[var]), -np.min(ds[var])])
#        cf= ax.pcolormesh(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
        cf= ax.contourf(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
    
    
    cb= fig.colorbar(cf, ax= ax, shrink=0.8)
    varlabel= var
    if plevel != None: varlabel = var+'_'+str(plevel)
#    cb.set_label(varlabel + ' ['+ ds[var].units + ']', size=11)    
    cb.set_label(ds[var].long_name + ' ['+ ds[var].units + ']', size=11)    




"""make the contour plot"""
if cont_var== 'msl' :
    cs= ax.contour(ds.lon, ds.lat, ds['msl'], np.arange(950, 1030, 2), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
    if numbers_cont: plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')

elif cont_var== 'z':   
    cs= ax.contour(ds.lon, ds.lat, dsorg['z'].sel(plev= plevel_contvar)/9.81, np.arange(0, 10000, 10), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
    if numbers_cont: plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')

elif cont_var=='t':
    cs= ax.contour(ds.lon, ds.lat, dsorg['t'].sel(plev= plevel_contvar), np.arange(200, 400, 1), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
    if numbers_cont: plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')



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
    sheer_rad= 200
    #plot the thickness
    cf= ax.pcolormesh(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap)
    
    cb= fig.colorbar(cf, ax= ax, shrink=0.8)
    varlabel= var
    if plevel != None: varlabel = var+'_'+str(plevel)
#    cb.set_label(varlabel + ' ['+ ds[var].units + ']', size=11)        
    cb.set_label(ds[var].long_name + ' ['+ ds[var].units + ']', size=11)        

    
    PlotWind(ds.lon, ds.lat, ds['u'].mean(dim='plev').values, ds['v'].mean(dim='plev').values, ax, nx= 15, ny=30, arrowtype= 'quiver', alen=15, scale=False)
    
    #calculate and plot the thermal wind
    f= np.sin(np.deg2rad(ds.lat)) * 2*2*np.pi/(24*60**2)
    f= np.tile(f, (len(ds.lon), 1)).T
    
    u_therm, v_therm= -1/f* grad_y(z_diff, ds.lat, ds.lon), 1/f*grad_x(z_diff, ds.lat, ds.lon)

    cs= ax.contour(ds.lon, ds.lat, z_diff/9.81, np.arange(0, 5000, 10), transform= ccrs.PlateCarree(), colors='r', linewidths= 1)
    if numbers_cont: plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')
    
    PlotWind(ds.lon[10:], ds.lat[10:], u_therm[10:, 10:], v_therm[10:, 10:], ax, nx= 15, ny=30, arrowtype= 'quiver', alen=15, scale=False, color= 'r')
    
    
    #calculate the shear        
    u_mean= value_in_rad(ds['u'].mean(dim='plev'), ds.lat, ds.lon, S_now.lat.values[0], S_now.lon.values[0], distance= sheer_rad)
    v_mean= value_in_rad(ds['v'].mean(dim='plev'), ds.lat, ds.lon, S_now.lat.values[0], S_now.lon.values[0], distance= sheer_rad)
    print('mean wind:', u_mean, v_mean)
    print('mean wind direction: ', UV2Direction(u_mean, v_mean) )
          
    u_therm_m= value_in_rad(u_therm, ds.lat, ds.lon, S_now.lat.values[0], S_now.lon.values[0], distance= sheer_rad, intype='np')
    v_therm_m= value_in_rad(v_therm, ds.lat, ds.lon, S_now.lat.values[0], S_now.lon.values[0], distance= sheer_rad, intype='np')    
    print('thermal wind:', u_therm_m, v_therm_m)
    print('thermal wind direction: ', UV2Direction(u_therm_m, v_therm_m) )

        
    alpha= angle_between((u_therm_m, v_therm_m), (u_mean, v_mean))    
    print('alpha:', alpha)



        
    




"""plot radius, steering wind"""
center_lon, center_lat= S_now.lon.values[0], S_now.lat.values[0]

if Plot_circle:
    #radius= S_now.Diameter.values[0]/2
    #plot_circle(ax, lon, lat, radius, edgecolor= 'r')
    plot_circle(ax, center_lon, center_lat, radius, edgecolor= 'green')


"""all my tracks"""
if Plot_all_my_tracks or Plot_Local_Vort_Max:
    test = 'version4'
    
    #get the correct tracktime_id for which datetime_now is between startdate and enddate
    
    file= "txt_files/track_list_PMC_appent.txt"
    tracktime_id, startyear, startmonth, startday, endyear, endmonth, endday= np.loadtxt(file).T.astype(int)
    startdate= np.array([datetime(startyear[i], startmonth[i], startday[i]) for i in range(len(startyear))])
    enddate= np.array([datetime(endyear[i], endmonth[i], endday[i], 23) for i in range(len(startyear))])
     
    #convert datetime_now from datetime64 to datetime.datetime:
    from datetime import datetime
    dd_datetime_now= datetime.utcfromtimestamp(datetime_now.values.astype(int)* 1E-9)
    
    tracktime_id_now= np.where(np.logical_and(dd_datetime_now > startdate, dd_datetime_now < enddate))[0]
    if len(tracktime_id_now) > 0:
        i = tracktime_id[tracktime_id_now ][0]
        print('PMC timegroup: ',str(i))
        track_dir= Mediadir+"/ERA5_STARS/tracks/"+test+"/tracks_"+str(i).zfill(3)
    
        t_dir= Path(".") / track_dir
        tr = TrackRun(t_dir, columns=['lon', 'lat', 'vo', 'time', 'area', 'vortex_type', 'slp'])
    
        ilegend= 0
        for (_, track) in tr.gb:
            ind_now= np.where((track.time == dd_datetime_now).values)[0]
            if len(ind_now)== 1:
                ilegend+=1
                if Plot_all_my_tracks:
                    if ilegend== 1: track.plot_track(ax=ax, zorder= 1, label='all tracks')
                    else: track.plot_track(ax=ax, zorder= 1, label='')
                
                tracklon, tracklat= track.lonlat[ind_now[0]]
                ax.scatter(tracklon, tracklat, transform=ccrs.PlateCarree(), s= 70, marker="s", color= 'g', zorder= 2)
    
    
    
    else: print('Not included in '+file)


if Plot_Rojo_track:
    Rojo = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv")
    Rojo['time']= Rojo.time.dt.round("H") #this can lead to a Rojo PL having two values for one time 

    #convert datetime_now to Timestamp   
    ts_datetime_now= pd.DatetimeIndex(datetime_now)[0]
        
    Rojo_nr_list= Rojo[Rojo['time']== ts_datetime_now].ID.values
    Rojo_nr_list= remove_dublicate(Rojo_nr_list)
    print('Rojo nr list: ', Rojo_nr_list)
    
    if len(Rojo_nr_list)== 0: Rojo_nr_list= [ID.split('_')[0]]

    for Rojo_nr in Rojo_nr_list:
        Rojo_now= Rojo[Rojo['ID'] == Rojo_nr] 
       
        R_tracks= prepare_tracks(Rojo_now, system_id_label="ID")
        
        for track in R_tracks:
            track.plot_track(ax=ax, color='red', label='Rojo track')   
            
            i_R_now= (np.abs(track.time - ts_datetime_now)).idxmin()
    
            ax.scatter(track.lon.values[i_R_now], track.lat.values[i_R_now], transform=ccrs.PlateCarree(), 
                       s= 150, marker= "x", linewidth=4, color= 'r', zorder= 2)
#            ax.text(track.lon.values[i_R_now], track.lat.values[i_R_now]+ 0.1, track.ID[0], transform=ccrs.PlateCarree(), color='r', fontsize=13, fontweight= 'bold')
        
            time_plot_iRojo= track.time.values[i_R_now]
            print('Rojo time: ', time_plot_iRojo)
    



if Plot_PL_track:
    ax.scatter(S_now.lon, S_now.lat, transform=ccrs.PlateCarree(), s= 150, linewidth= 4, marker= "x", color= 'b', zorder=2, label='')
#    ax.scatter(S_now.lon, S_now.lat, transform=ccrs.PlateCarree(), s= 150, marker= "x", color= 'k', zorder=2, label='')

    S_now_whole= S.loc[S['ID'] == ID]
    S_now_track= prepare_tracks(S_now_whole, system_id_label="ID")[0]
    S_now_track.plot_track(ax= ax, color='b', label='matched track')

#    ax.scatter(S_now_whole.lon, S_now_whole.lat, transform=ccrs.PlateCarree(), s= 70, marker="x", color= 'g', zorder= 2, label='')


if Plot_PLCG or Plot_smth_track:
    #calculate the beering of the track
    S_data= np.array([S_now_whole.lon.values, S_now_whole.lat.values])
    S_t= np.arange(len(S_now_whole.lon.values))
    
    
    sp = csaps.CubicSmoothingSpline(S_t, S_data, smooth= 0.001)
    Sm_lon, Sm_lat= sp(S_t)
#    ax.scatter(Sm_lon, Sm_lat, transform=ccrs.PlateCarree(), s= 70, marker="x", color= 'k', zorder= 2)
    if Plot_smth_track:
        ax.plot(Sm_lon, Sm_lat, transform=ccrs.PlateCarree(), color= 'k', label= 'smoothed track', linewidth= 2.5, zorder=1)
#        ax.scatter(Sm_lon[Obs-1], Sm_lat[Obs-1], transform=ccrs.PlateCarree(), s= 70, marker= "s", color= 'k', zorder=2)
    
    #obs is index +1
    beering= gcc.bearing_at_p1( (Sm_lon[Obs-2], Sm_lat[Obs-2]), (Sm_lon[Obs], Sm_lat[Obs]) )
    
    
    grid_dist= 25 #in km
    ncells= radius//grid_dist
    
    #the tangental axis
    tanax= [list(gcc.point_given_start_and_bearing((center_lon, center_lat), beering, n*grid_dist*1E3)) for n in np.arange(-ncells, ncells+1)]
    tanax= np.array(tanax)
    
    #the tanax
    #ax.scatter(tanax[:,0], tanax[:,1], transform=ccrs.PlateCarree(), s= 3, marker="x", color= 'r', zorder= 1)
    
    point0= gcc.point_given_start_and_bearing((center_lon, center_lat), beering, -(radius+grid_dist)*1E3) #point one before the start of the tangential axis, used for calculation of beering_axis
    
    beering_axis= [gcc.bearing_at_p2((point0), (tanax[m,0], tanax[m,1])) for m in range(len(tanax))] #the wind direction along the tangential axis
    
    center_grid= [[list(gcc.point_given_start_and_bearing((tuple(tanax[m])), (beering_axis[m]-90)%360, n*grid_dist*1E3)) for m in range(len(tanax))] for n in np.arange(-ncells, ncells+1)]
    center_grid= np.array(center_grid)
    
    if Plot_PLCG:
        ax.scatter(center_grid[:,:,0], center_grid[:,:,1], transform=ccrs.PlateCarree(), s= 1, marker="x", color= 'k', zorder= 1)
    



if legend: plt.legend()



#plt.tight_layout()


if save:
    savefile= savedir+ 'StollID'+ID+ '_Obsnr'+str(Obs)+ '_'+var+str(plevel)+save_extra
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight') 



if Plot_PLCG:
    """ interpolate to local grid"""
    fig= plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    
    from scipy.interpolate import griddata
    
    longrid, latgrid= np.meshgrid(ds.lon, ds.lat)
    
    var_interp= griddata(tuple([np.ravel(longrid), np.ravel(latgrid)]) , np.ravel(ds[var]), (center_grid[:,:,0], center_grid[:,:,1]), method='linear')
    cont_var_interp= griddata(tuple([np.ravel(longrid), np.ravel(latgrid)]) , np.ravel(dsorg[cont_var].sel(plev=plevel_contvar)), (center_grid[:,:,0], center_grid[:,:,1]), method='linear')
    
    dist_grid_axis= np.arange(-ncells, ncells+1)*grid_dist
    
    
    if clevels != None:
        cf= plt.contourf(dist_grid_axis, dist_grid_axis, var_interp, cmap= cmap, levels= clevels, vmin= clevels[0], vmax=clevels[-1], extend='both')            
    elif vmax != None:
        cf= plt.contourf(dist_grid_axis, dist_grid_axis, var_interp, cmap= cmap, vmin= vmin, vmax=vmax)        
    elif sym== False: cf= plt.contourf(dist_grid_axis, dist_grid_axis, var_interp, cmap= cmap)
    elif sym:
        cf= plt.contourf(dist_grid_axis, dist_grid_axis, var_interp, cmap=cmap, vmin= -vextr, vmax= vextr)
    
    cb= plt.colorbar(cf, shrink= .8)
    varlabel= var
    if plevel != None: varlabel = var+'_'+str(plevel)
#    cb.set_label(varlabel + ' ['+ ds[var].units + ']', size=11)   
#    cb.set_label(ds[var].long_name + ' ['+ ds[var].units + ']', size=11)
    if var == 'vo':
        cb.set_label('Relative vorticity ['+ ds[var].units + ']', size=11) 
    else:        
        cb.set_label(ds[var].long_name + ' ['+ ds[var].units + ']', size=11) 
    
    #plot the contour variable
    if cont_var== 'z':
        cont_var_interp/= 9.81
        cont_levels= np.arange(0, 10000, 10)
    else: print('have to specify cont_levels')
    contours= plt.contour(dist_grid_axis, dist_grid_axis, cont_var_interp, levels= cont_levels, colors='k', linewidths= 1)
    if numbers_cont: plt.clabel(contours, inline=True, fontsize=8, fmt='%1.0f') 
    
    plt.xlabel('Distance in propagation direction [km]', size= 12)
    plt.ylabel('Distance orthogonal to propagation [km]', size= 12)
    plt.yticks(np.arange(-250,501, 250))
    plt.xticks(np.arange(-500,501, 250))
    plt.scatter(0,0, marker= 'x', color= 'k', zorder= 5, s= 150, linewidth= 2)

#    plt.scatter(0,0,  color='k', s= 70, marker= "s")
    
    
    """make the propagation arrow"""
    arrow_ax = fig.add_axes([0.78, 0.04, 0.12, 0.08])
    import matplotlib.patches as mpatches
    
    arrow = mpatches.FancyArrowPatch((0.2, .8), (.7, .75), mutation_scale=20)
    arrow_ax.add_patch(arrow)
    arrow_ax.text(0,-0.2 , 'Propagation\ndirection', fontsize=12) #, fontweight= 'bold'
    arrow_ax.set_frame_on(False)
    arrow_ax.get_xaxis().set_visible(False)
    arrow_ax.get_yaxis().set_visible(False)
        
    if save:
        savefile= savedir+ 'StollID'+ID+ '_Obsnr'+str(Obs)+ '_'+var+str(plevel)+'_onPLCG_smth-track'
        print(savefile)
        plt.savefig(savefile , bbox_inches='tight') 


