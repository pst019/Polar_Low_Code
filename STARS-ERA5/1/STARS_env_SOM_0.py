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

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from f_useful import *
from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs


#save= True
save= False

fignr= 9

#imp_data=True
imp_data=False

maptype= 'PlateCarree'
#maptype= 'Polar_Stereo'




S = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv")
S_ind= STARS_individual_systems(S)






"""get the ERA-5 data"""
import xarray as xr

plevel=None
vmin,vmax=None,None

"""specify field"""
filetype='surface'
#var= '2t'
#var= '10u'
#var= '10v'
#var= 'U'
#var= 'skt'
var='msl'

#filetype='plevels'
#plevel= 850
#var='U'
#var= 'u'
#var='Vort'
#var='t'
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
#var='z'

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


"""the contour variable"""
if plevel== None: cont_var= 'msl'
else: cont_var='z'
    #cont_var='t'


sym= False
#cmap='jet'
#cmap= 'viridis'
#cmap= 'Greys'
#cmap= 'Greys_r'
#cmap= 'Blues'
#cmap= 'Reds'

#sym= True
cmap= 'RdBu_r'







if imp_data:
    """specify system ID, Obsnr"""
    Obs= 1
    
    for i, ID in enumerate(remove_dublicate(S['ID']))[:5]:
    
    
        S_now= S[np.logical_and(S['ID'] == ID, S['Obs'] == Obs ) ]
    
        datetime_now= S_now.datetime.dt.round("H")
        
    
        filetime= str(datetime_now.dt.year.values[0])+'_'+str(datetime_now.dt.month.values[0]).zfill(2)+'_'+str(datetime_now.dt.day.values[0]).zfill(2)
        hour= datetime_now.dt.hour.values[0]
    
        print(datetime_now)
    
    
    
        """general"""
    #    if plevel==None:
    #        d0= xr.open_dataset(Mediadir + "ERA5_STARS/surface_era5_"+ filetime + '.nc')
    #        d0= d0.isel(time= hour)
    #    
    #    else:
    #        d0= xr.open_dataset(Mediadir + "ERA5_STARS/vorticity_era5_"+ filetime + '.nc')
    #        d0= d0.isel(time= hour)
        
        ds0= xr.open_dataset(Mediadir + "ERA5_STARS/data/"+filetype+ "_era5_"+ filetime + '.nc')
        ds0= ds0.isel(time= hour)
        ds0= ds0[var]
    
        if filetype in ['plevels', 'vorticity']:
            ds0['plev']/= 100
            ds0= ds0.sel(plev= plevel)
    
        if i== 0: ds= ds0
        else: ds= xr.concat([ds, ds0], 'time')
        
       
    
    #
    #
    #
    #if filetype=='surface':
    #    ds['u']= ds['10u']
    #    ds['v']= ds['10v']
    #
    #    
    #if var== 'U':
    #    ds['U']= np.sqrt(ds['u']**2 + ds['v']**2)
    #    ds[var].attrs['units']= 'm/s'
    #
    #if var=='vo':
    #    ds[var]*= 1E5
    #    ds[var].attrs['units']= '10$^{-5}$ 1/s'    
    #
    #if var== 'Vort': 
    #    vort= grad_x(ds.v, ds.lat, ds.lon) - grad_y(ds.u, ds.lat, ds.lon)
    #    vort *= 1E5
    #
    #    ds[var]= (('lat', 'lon'), vort)
    #    ds[var].attrs['units']= '10**-5 s**-1' 
    #
    #
    #
    #
    #if var=='skt-t':
    #    d0= xr.open_dataset(Mediadir + "ERA5_STARS/surface_era5_"+ filetime + '.nc')
    #    d0= d0.isel(time= hour)
    #
    #    ds= xr.merge([d0, ds])    
    #    
    #    ds[var]= ds['skt']- ds['t']
    #    ds[var].attrs['units']= 'K'   
    #
    #if var=='sst-t':
    #    d0= xr.open_dataset(Mediadir + "ERA5_STARS/boundary_era5_"+ filetime + '.nc')
    #    d0= d0.isel(time= hour)
    #
    #    ds= xr.merge([d0, ds])    
    #    
    #    ds[var]= ds['sst']- ds['t']
    #    ds[var].attrs['units']= 'K' 
    #
    #
    #if var== 'grad_t':
    #    variable= np.sqrt(grad_x(ds.t, ds.lat, ds.lon)**2 + grad_y(ds.t, ds.lat, ds.lon)**2) *1E5
    #    ds[var]= (('lat', 'lon'),  variable)
    #
    #    ds[var].attrs['units']= 'K/100 km'
    #
    #
    #
    if var== 'msl':
        ds/= 100
        ds.attrs['units']= 'hPa'
    
    if var=='z':
        ds/=9.81
        ds.attrs['units']= 'm'
        
    #if var in ['ttr', 'sshf', 'slhf']:
    #    ds[var]/= -3600
    #    ds[var].attrs['units']= 'W m**-2'
    #if var in ['sf', 'lsp', 'cp', 'tp']:
    #    ds[var] *= 1000
    #    ds[var].attrs['units']= 'mm/h' 
    #
    #
    #
    

    dsraw= ds

#ds= ds.where(np.logical_and(ds.lat > 60, ds.lon <=60), drop=True)


from minisom import MiniSom



dsreshaped= np.reshape(ds.values, (len(ds.time), -1))


som = MiniSom(6, 6, dsreshaped.shape[1], sigma=0.3, learning_rate=0.5) # initialization of 6x6 SOM

#som.pca_weights_init(dsreshaped)
som.train_random(dsreshaped, 100)

#"""plot the mean"""
#if maptype== 'Polar_Stereo': fig = plt.figure(fignr, figsize= (8,6) )
#elif maptype== 'PlateCarree': fig = plt.figure(fignr, figsize= (8,6) )
#    
#fignr+=1
#plt.clf()
#
#if maptype== 'Polar_Stereo': ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (2,1,1))
#elif maptype== 'PlateCarree': ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-40, 85, 50, 85], subplot= (2,1,1))
#scale_bar(ax, 500, location= (0.06, 0.04))
#
#
# 
#dsmean= ds.mean(dim='time')    
#    
#if vmax != None:
#    cf= ax.pcolormesh(ds.lon, ds.lat, dsmean, transform= ccrs.PlateCarree(), cmap= cmap, vmin= vmin, vmax=vmax)
#    
#elif sym== False:
#    cf= ax.pcolormesh(ds.lon, ds.lat, dsmean, transform= ccrs.PlateCarree(), cmap= cmap)
#    cs= ax.contour(ds.lon, ds.lat, dsmean, transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
#    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')
#
#
#elif sym== True:
#    vextr= np.max([np.max(dsmean), -np.min(dsmean)])
#    cf= ax.pcolormesh(ds.lon, ds.lat, dsmean, transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
#
#
#cb= fig.colorbar(cf, ax= ax, shrink=0.7)
#varlabel= var
#if plevel != None: varlabel = var+'_'+str(plevel)
#cb.set_label(varlabel + ' ['+ ds.units + ']', size=14)  
#plt.title('Mean')
#
#
#"""plot the standard deviation"""
#if maptype== 'Polar_Stereo': ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (2,1,2))
#elif maptype== 'PlateCarree': ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-40, 85, 50, 85], subplot= (2,1,2))
#scale_bar(ax, 500, location= (0.06, 0.04))
#
#cf= ax.pcolormesh(ds.lon, ds.lat, ds.std(dim='time'), transform= ccrs.PlateCarree(), cmap= 'Reds')
#cb= fig.colorbar(cf, ax= ax, shrink=0.7)
#
#cs= ax.contour(ds.lon, ds.lat, ds.std(dim='time'), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
#plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')
#
#cb.set_label(varlabel + ' ['+ ds.units + ']', size=14)  
#
#plt.title('Standard deviation')
#
#
#
#
#"""calculate the EOFs"""
#from eofs.xarray import Eof
#
##dsano= ds - dsmean
#
#if maptype== 'Polar_Stereo':   fig = plt.figure(fignr, figsize= (8,6) )
#if maptype== 'PlateCarree': fig = plt.figure(fignr, figsize= (10,4) )
#fignr+=1
#plt.clf()
#
#
#coslat= np.cos(np.deg2rad(ds.coords['lat'].values)).clip(0., 1.)
#wgts = np.sqrt(coslat)[..., np.newaxis]
#solver = Eof(ds, weights=wgts)
#
#EOFnr= 4
#
#
#var_frac= np.round(solver.varianceFraction(neigs=EOFnr).values*100, 1)
#print('Variance fraction [%]: ', var_frac)
#
#eigenvalue_error= solver.northTest(neigs=EOFnr, vfscaled=True)
#print('Eigenvalue error [%]: ', np.round(eigenvalue_error.values*100, 1))
#
#
#
#solver.pcs(npcs= EOFnr)
#
#
#eof1 = solver.eofsAsCovariance(neofs=EOFnr, pcscaling= 1)#Empirical orthogonal functions (EOFs) expressed as the covariance between the PCs and the time series of the Eof input dataset at each grid point.  pcscaling= 1 - scaling to unit variance
##eof1 = solver.eofsAsCorrelation(neofs=EOFnr) #Empirical orthogonal functions (EOFs) expressed as the correlation between PCs and the time series of the Eof input dataset at each grid point.
##eof1 = solver.eofs(neofs=EOFnr)
#
#
#"""find which PC is most important for each PL"""
#PC= solver.pcs(npcs= EOFnr)
#PCabs= np.abs(PC)
#MostImpPCofPL= np.where( PCabs== np.max(PCabs, axis= 1))[1] #finds the most important PC for each PL
#PCSign= np.sign(PC.values[(PCabs== np.max(PCabs, axis= 1)).values]) #the sign of the most important PC
#
#for iEOFnr in range(EOFnr):
#
#    if maptype== 'Polar_Stereo': ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (2,2,iEOFnr+1))
#    if maptype== 'PlateCarree': ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-40, 85, 50, 85], subplot= (2,2,iEOFnr+1))
#    #[-40, 50, 85, 85]
#    
#    #proj = ccrs.Orthographic(central_longitude=-20, central_latitude=60)
#    #ax = plt.axes(projection=proj)
#    #ax.coastlines()
#    #ax.set_global()
#    
#    scale_bar(ax, 500, location= (0.06, 0.04))
#    
#    vextr= np.max([np.max(eof1), -np.min(eof1)])
#    cf= ax.pcolormesh(eof1.lon, eof1.lat, eof1.sel(mode= iEOFnr), transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
#    cb= fig.colorbar(cf, ax= ax, shrink=0.7)
#    
#    cs= ax.contour(eof1.lon, eof1.lat, eof1.sel(mode= iEOFnr), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
#    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')
##    varlabel= var
##    if plevel != None: varlabel = var+'_'+str(plevel)
##    #cb.set_label(varlabel + ' ['+ ds.units + ']', size=14)  
#    plt.title('EOF '+str(iEOFnr+1)+' ('+str(int(np.round(var_frac[iEOFnr])))+'%)')
#
#    #plot the points where the EOF is in positive phase in yellow
#    PLsThisPC= S[S['Obs'] == Obs ][np.logical_and(MostImpPCofPL == iEOFnr, PCSign >= 1)] #get all PLs that are connected to this PC
#    ax.scatter(PLsThisPC.lon, PLsThisPC.lat, transform=ccrs.PlateCarree(), s= 1, color='y')
#
#    PLsThisPC= S[S['Obs'] == Obs ][np.logical_and(MostImpPCofPL == iEOFnr, PCSign <= -1)] #get all PLs that are connected to this PC
#    ax.scatter(PLsThisPC.lon, PLsThisPC.lat, transform=ccrs.PlateCarree(), s= 1, color= 'g')


