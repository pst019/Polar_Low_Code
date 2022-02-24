#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

This is from plot_comparison
"""

import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'
AROMEdir = Mediadir + 'PL/AA/operational/'


from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt


import numpy as np
import time
#from datetime import date, datetime, timedelta
from scipy import stats
import scipy.ndimage.filters as filters
import xarray as xr

# import own modules
import sys  #to import the functions from a different directory
sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')

from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_imp_AROME import *  # Read in netcdf file
import f_imp_thorpex as thorpex
from f_meteo import *
from f_useful import *


plt.rcParams.update({'font.size': 12})

"""global variables"""
fignr= 5

#maptype='AA'
maptype='AA_half'
#maptype='Lambert'

if maptype=='AA': plt.figure(fignr, figsize= (6, 6))
elif maptype=='AA_half': plt.figure(fignr, figsize= (6, 4.7))
else: plt.figure(fignr, figsize= (6, 4.5))

fignr+=1
 
"""time specification: time of the plot"""
year, month, day, hour = 2018, 4, 6, 0

t= 24 #leadtime

#option= 'sfx' #no SLP otherwise good
option= 'extracted' #lots, ml,pl,zl
#option= 'pp' #most for the surface

    

tlist= [t]
#tlist=[6, 12, 18, 24]

"""c) plot characteristics"""

#var= 'PressWind'
#var= 'PressWind_advanced'
#var= 'PressWind_geop'


#var= 'Precip'

#var= 'Cloud_high' 
#var= 'Cloud_med'    
#var= 'Cloud_low'

#var, option='netrad', 'sfx'
#var, option= 'OLW', 'extracted' #top of atmosphere outgoing longwave radiation
#
var= 'Vort_lev'
#
#var= 'PBH_advanced' #CAPE as contours
#
#var= 'Thickness_advanced' #black contours: grad theta 850, green contours: theta_e 500- thetha_e SST


#var='surface_geop'

pleveldist=1
#var= 'PressWindarrow'
#var= 'PressWindVortmax'
#var= 'SLP_laplacian'
#var= 'SLP_filtered'
#var= 'MeanSLPminSLP'
#var= 'Theta_lev'
#
#var= 'Theta_e_lev'
#var= 'Front_lev'
#var= 'Grad_Theta_lev'
#var= 'Grad_Theta_e_lev'

#var='SST'
#var='T2m'
#var= 'SST-T500'
#var= 'thetaSST-theta500'
#var= 'theta_eSST-theta_e500'

#bounds= np.arange(264, 286, 1) #for pottempe850
#var='Shear_lev'
#var='Abs_shear_lev'

#var= 'Stretching_lev'
#var= 'OW_lev'
#var= 'Vort_lev_T'
#var= 'Thickness'
#var= 'test_rot'
#var= 'Geop_lev'
#var= 'Geop_wind_lev' #seems like the arrow is wrong

#var= 'Geop_laplacian_lev'
#var= 'RH_lev'
#var= 'SH_lev'
#var= 'RH_lev_fromSH'

#var= 'TH'
#var= 'TH_thor' #for comparison with thorpex
#var= 'LH'
#var= 'LH_thor'
#var= 'FLX'
#var= 'PBH'

#var= 'Precip'
#var='LH_release'
#var= 'Precip_acc'

#var= 'W_lev'
#var= 'W_lev_filter'
#var= 'PV_lev'
#var= 'Cloud_conv' # is somehow always 0.

#var= 'CAPE'
#var= 'WVTC'
#var= 'CWR'


#var, display='tropopause', 'z'
#var, display='tropopause', 'pres'
#var='tropopause_advanced'
#var='dynamical_trop'


if '_lev' in var: option= 'extracted'

filename= 'arome_arctic_'+option+'_2_5km_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'T'+str(hour).zfill(2)+'Z.nc'


""" presure level (if the variable ends with '_var' """
plev= 850 #hPa


Plotbox= False
boxcolor='k'

PlotCentre=True
IceEdge= False
    

save=False
#save= True

AAres=1 #2- every second datapoint is taken


title_extra='' #'jet' #''


"""Lambert coordinates"""
#lllon, lllat, urlon, urlat= -15, 65, 50, 75
#lat0, lon0= 75, 0 #(lat0, lon0) = center point


savedir= homedir+'Polar_Low/Arome-Arctic-operational/Figs/'            
                 
"""end global variables"""




for t in tlist:
    print('t', t)
    plt.clf()


    if maptype== 'AA': map = AA_map()
    elif maptype== 'AA_half': map = AA_map_half()
    else: map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0, latdist=2, londist= 5)    
    
    
    ds= xr.open_dataset(AROMEdir + filename)
    if var != 'OLW': ds= ds.isel(time= t)

    if '_lev' in var: ds= ds.sel(pressure= plev)
                                            
##    d= data(filename= AROMEdir + filename, res= AAres)
#    
    Lon,Lat = map(ds.longitude.values , ds.latitude.values)
#        
#        
#    d.imp_surf(tn= t)
#    
#    if var in ['FLX', 'LH', 'LH_thor', 'TH', 'TH_thor', 'PBH', 'PBH_advanced', 'CAPE', 'CIN', 'LCL', 'LFC', 'OLW','Precip','Precip_acc', 'LH_release'] or 'Cloud' in var:
#        d.imp_atm(tn= t)
#    
#    if '_lev' in var:
#        d.imp_level(pn= pn, tn= t)
#    
#    if 'tropopause' in var:
#        d.imp_grib(tn= t, calcheight= True)
#    
#    if var== 'dynamical_trop':
#        d.imp_standard_t(tn= t)
#    
#    
#    
    """MSLP, Surface temperature and surface winds"""
    if var== 'PressWind':
        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= 1)
        PlotWindVelo(Lon, Lat, np.sqrt(ds.x_wind_10m.values[0]**2+ ds.y_wind_10m.values[0]**2), map, Umax= 25)
#    
#    if var== 'PressWindarrow':
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= 2)
#    #    PlotWindVelo(Lon, Lat, np.sqrt(d.u10m**2+ d.v10m**2), map, Umax= 25)
#        PlotWind(ds.longitude.values, ds.latitude.values, d.u10m, d.v10m, map, rot=False, arrowtype='barb') #put this after PlotColorMap, otherwise it disappears
#                          
#    
#    """MSLP, Surface temperature and surface winds -advanced"""
#    if var== 'PressWind_advanced':
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= 1, numbers=False)
#    
#        U= np.sqrt(d.u10m**2+d.v10m**2)
#        PlotWindVelo(Lon, Lat, U, map, Umax= 25, color='YlBu')    
#        #find local pressure min and plot them in the map   
#        PlotLocalMax(ds.air_pressure_at_sea_level.values[0]/100, threshold=1010, distance=100/AAres, map= map, lon=ds.longitude.values, lat=ds.latitude.values,
#                     data2=U, threshold2=18, distance2=150/AAres)    
#      
#        PlotLocalMax(U, threshold=20, distance=100/AAres, map= map, lon=ds.longitude.values, lat=ds.latitude.values, typ='max',
#                     color='orange', dot=False, roundorder= 0)    
# 
    """MSLP, Surface temperature and surface winds, geopotential height at 500hPa"""
    if var== 'PressWind_geop':
        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= 1, numbers=False)

        plev= 500
        PlotContours(Lon, Lat, 1/9.81* ds.geopotential_pl.sel(pressure= plev).values, map, leveldist= 20, numbers= True, color='red')


        U= np.sqrt(ds.x_wind_10m.values[0]**2+ ds.y_wind_10m.values[0]**2)
        PlotWindVelo(Lon, Lat, U, map, Umax= 25, color='YlBu')
   
        #find local pressure min and plot them in the map   
        PlotLocalMax(ds.air_pressure_at_sea_level.values[0]/100, threshold=1010, distance=100/AAres,
                     map= map, lon=ds.longitude.values, lat=ds.latitude.values, #lonbound= [-15, 40], latbound= [65, 80],
                     data2=U, threshold2=18, distance2=150/AAres, color='k')    
#      
        PlotLocalMax(U, threshold=15, distance=200/AAres, map= map, 
                     lon=ds.longitude.values, lat=ds.latitude.values, typ='max', lonbound= [-15, 40], latbound= [65, 80],
                     color='orange', dot=False, roundorder= 0)      
#    
#    
#    if var== 'PressWindVortmax':
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= 2, numbers=False)
#    
#        U= np.sqrt(d.u10m**2+d.v10m**2)
#    #    PlotWindVelo(Lon, Lat, U, map, Umax= 25)    
#        #find local pressure min and plot them in the map   
#    #    PlotLocalMax(ds.air_pressure_at_sea_level.values[0]/100, threshold=1010, distance=100/AAres, map= map, lon=ds.longitude.values, lat=ds.latitude.values,
#    #                 data2=U, threshold2=16, distance2=80/AAres)    
#        PlotWind(ds.longitude.values, ds.latitude.values, d.u10m, d.v10m, map, rot=False, arrowtype='barb') #put this after PlotColorMap, otherwise it disappears
#    
#      
#        PlotLocalMax(U, threshold=20, distance=100/AAres, map= map, lon=ds.longitude.values, lat=ds.latitude.values, typ='max',
#                     color='r', dot=False, roundorder= 0)    
#    
#    
#        d.imp_level(pn= 4, tn= t)
#        dx= 2500* AAres
#        vort= (np.gradient(d.v, dx, axis= 1)- np.gradient(d.u, dx, axis= 0))*1E5         
#        #probably one has to rotate u and v first to proper x and y coordinates
#        gausfilterdist= 60
#        vortfilter= filters.gaussian_filter(vort, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#    #    PlotColorMap4(Lon, Lat, vortfilter, map, label= r"$\xi_{"+str(d.chosen_plevel)+"}$ [10$^{-5}$ 1/s]" #Gauss filter of std "+str(gausfilterdist)+"km"
#    #    , color='RdBu', symetric=True)
#    
#        PlotLocalMax(vortfilter, threshold=10, distance=250/(AAres*2.5), map= map, lon=ds.longitude.values, lat=ds.latitude.values, typ='max', color='y', roundorder= 0)    
#    
#    
#    """potential temperature at 700 hPa"""
#    if var== 'Theta_lev':
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist, numbers= False)
#    #    PlotContours(Lon, Lat, d.geop, map, leveldist= 10)
#        if d.chosen_plevel == 500: bounds= np.arange(278, 288, 1)
#        elif d.chosen_plevel == 700: bounds= np.arange(270, 282, 1) #for pot temp 700
#        elif d.chosen_plevel == 850: bounds= np.arange(264, 280, 1) #for pot temp 850
#        elif d.chosen_plevel >= 950: bounds= np.arange(260, 278, 1) #for pot temp 950 
#        else: bounds= None
#        
##        bounds= np.arange(267, 290.5, 1)
##        bounds= np.arange(260, 291, 2)
#        
#        theta= PotTemp(d.T, d.chosen_plevel)
#        new_map, norm =PlotColorMap4(Lon, Lat, theta, map, bounds= bounds, color='RdBu', label=r"$\theta_{"+str(d.chosen_plevel)+"}$ [K]", return_colormap= True)
#
#
#    """ equivalent potential temperature at 850hPa"""
#    if var== 'Theta_e_lev':
#        bounds= None
##        spechum= d.RH/(.263 *d.pres[pn])* np.exp(17.67*(d.T-273.15) / (d.T-29.65))
#        theta_e= EquiPotTemp(d.T, d.SH, d.chosen_plevel)
#        
#        new_map, norm =PlotColorMap4(Lon, Lat, theta_e, map, bounds= bounds, color='RdBu', label=r"$\theta_{e,"+str(d.chosen_plevel)+"}$ [K]", return_colormap= True)
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#
#
#    if var== 'SST':
#        PlotContours(Lon, Lat, d.T0m -273.15 , map, leveldist= 1)
#        PlotColorMap4(Lon, Lat, d.T0m -273.15 , map, symetric=True, maxlevel= 10, boxnr= 20, color= 'RdBu', label=r"Surface temperature [$^{\circ}$C]")
#
#    if var== 'T2m':
#        #        Tsurfbounds= np.array([-40, -20, -10, -5, -2, 0,  2,  5, 10, 20, 40])  
#
#        PlotContours(Lon, Lat, d.T2m -273.15 , map, leveldist= 1)
#        PlotColorMap4(Lon, Lat, d.T2m -273.15 , map, symetric=True, maxlevel= 10, boxnr= 20, color= 'RdBu', label=r"2m temperature [$^{\circ}$C]")
#
#
#    """stability"""
#    if var== 'thetaSST-theta500':         
#        d.imp_level(pn= 8, tn= t) # 8 = 500hPa
#        theta500= PotTemp(d.T, d.chosen_plevel)
#        thetaSST= PotTemp(d.T0m, ds.air_pressure_at_sea_level.values[0]/100)  
#        gausfilterdist= 25
#        stab_filter= filters.gaussian_filter((thetaSST- theta500), sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#
#        PlotStaticStability(Lon, Lat, stab_filter , map, difftype='thetaSST-theta500')
##
#        
#    if var== 'SST-T500':         
#        d.imp_level(pn= 8, tn= t) # 8 = 500hPa 
#        gausfilterdist= 25
#        stab_filter= filters.gaussian_filter((d.T0m- d.T), sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#
#        PlotStaticStability(Lon, Lat, stab_filter , map, difftype='SST-T500')
#
#
#    if var== 'theta_eSST-theta_e500':         
#        d.imp_level(pn= 8, tn= t) # 8 = 500hPa
#        theta_e500= EquiPotTemp(d.T, d.SH, d.chosen_plevel)
#        
#        d.imp_level(pn= 0, tn=t)
#        theta_eSST= EquiPotTemp(d.T0m, d.SH, ds.air_pressure_at_sea_level.values[0]/100)  
##        PlotStaticStability(Lon, Lat, (theta_eSST- theta_e500) , map, difftype='theta_eSST-theta_e500')
#
#        gausfilterdist= 25
#        stab_filter= filters.gaussian_filter((theta_eSST- theta_e500), sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#
#        PlotStaticStability(Lon, Lat, stab_filter , map, difftype='theta_eSST-theta_e500')
#        PlotColorMap4(Lon, Lat, stab_filter, map, symetric= True, maxlevel=25, color='RdBu', label=r"$\theta_{e,SST}-\theta_{e,500}$ [K]")
#
#
#    if var== 'Grad_Theta_lev':
#        gausfilterdist = 100 #20 in km
#        theta= PotTemp(d.T, d.chosen_plevel)
#        theta_filter= filters.gaussian_filter(theta, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#        
#        PlotContours(Lon, Lat, theta_filter, map, leveldist= 2)
#        
#        dx= 2500* AAres
#        grad_theta_filter= np.gradient(theta_filter, dx, edge_order=1)
#        abs_grad_theta_filter= np.sqrt(grad_theta_filter[0]**2+ grad_theta_filter[1]**2)
#
##        bounds=np.arange(0,np.max(abs_grad_theta_filter)*1E5 +1,1)
#        bounds=np.arange(0, 5.6, 0.5)
#        PlotColorMap4(Lon, Lat, abs_grad_theta_filter *1E3*1E2, map, color= 'red', bounds=bounds, label=r"$\nabla \theta_{"+str(d.chosen_plevel)+"}$ [K/100km]")
#
#    
#
#    if var== 'Front_lev':
#        gausfilterdist = 50 #in km
#        theta_e= EquiPotTemp(d.T, d.SH, d.chosen_plevel)
#        
#        theta_e_filter= filters.gaussian_filter(theta_e, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#        PlotContours(Lon, Lat, theta_e_filter, map, leveldist= 2)
#        
#        dx= 2500* AAres
#        grad_theta_e_filter= np.gradient(theta_e_filter, dx, edge_order=1)
#        abs_grad_theta_e_filter= np.sqrt(grad_theta_e_filter[0]**2+ grad_theta_e_filter[1]**2)
#
##        PlotColorMap4(Lon, Lat, abs_grad_theta_e_filter *1E3*1E2, map, color= 'red', bounds=np.arange(0,10.1,1), label=r"$\nabla \theta_{e,"+str(d.chosen_plevel)+"}$ [K/100km]")
#        
#        TFP= - np.sum(np.array(np.gradient(abs_grad_theta_e_filter, dx, edge_order=1)) * grad_theta_e_filter/abs_grad_theta_e_filter  , axis= 0)
#        TFP_filter= filters.gaussian_filter(TFP, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#        #the laplacian is to find the local maxima - where it is 0
#        
##        PlotColorMap4(Lon, Lat, TFP_filter, map, color= 'RdBu', symetric=True, label=r"TFP$_{"+str(d.chosen_plevel)+"}$")
#
#        front= np.zeros(np.shape(theta_e))
#        front[np.logical_and(np.abs(TFP_filter)< .7* 1E-10, abs_grad_theta_e_filter >4*1E-5)]= 1
#             
#        PlotColorMap4(Lon, Lat, front, map, color= 'red')#, label=r"$\nabla \theta_{e,"+str(d.chosen_plevel)+"}$ [K/100km]")
#
#
#
#    if var== 'Grad_Theta_e_lev':
#        gausfilterdist = 50 #20 #in km
#        theta_e= EquiPotTemp(d.T, d.SH, d.chosen_plevel)
#        
##        PlotColorMap3(Lon, Lat, theta_e, map, bounds= bounds, label=r"$\theta_{e,"+str(d.chosen_plevel)+"}$ [K]")
#        
#        theta_e_filter= filters.gaussian_filter(theta_e, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#        PlotContours(Lon, Lat, theta_e_filter, map, leveldist= 2)
#
#
#        dx= 2500* AAres
#        grad_theta_e_filter= np.gradient(theta_e_filter, dx, edge_order=1)
#        abs_grad_theta_e_filter= np.sqrt(grad_theta_e_filter[0]**2+ grad_theta_e_filter[1]**2)
##        abs_grad_theta_e_filter= filters.gaussian_filter(abs_grad_theta_e, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 4.)
#
#        PlotColorMap4(Lon, Lat, abs_grad_theta_e_filter *1E3*1E2, map, color= 'red', bounds=np.arange(0,10.1,1), label=r"$\nabla \theta_{e,"+str(d.chosen_plevel)+"}$ [K/100km]")
#
#
#        PlotLocalMax(abs_grad_theta_e_filter *1E5, threshold=4, distance=250/(AAres*2.5), map= map,
#                     lon=ds.longitude.values, lat=ds.latitude.values, typ='max', color='blue', roundorder= 1)    
#
#        
#        #SLP min
#        U= np.sqrt(d.u10m**2+d.v10m**2)
#        PlotLocalMax(ds.air_pressure_at_sea_level.values[0]/100, threshold=1010, distance=100/AAres, map= map, lon=ds.longitude.values, lat=ds.latitude.values,
#                     data2=U, threshold2=18, distance2=150/AAres) 
#
#        #vort max
#        dx= 2500* AAres
#        vort= (np.gradient(d.v, dx, axis= 1)- np.gradient(d.u, dx, axis= 0))*1E5         
#        gausfilterdist= 100# 40 #100 #60
#        vortfilter= filters.gaussian_filter(vort, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 2)
#
#        PlotLocalMax(vortfilter, threshold=10, distance=250/(AAres*2.5), map= map, lon=ds.longitude.values, lat=ds.latitude.values, typ='max', color='black', roundorder= 0)    
#
#
#
#       
#    """mean(SLP) - SLP"""
#    if var== 'MeanSLPminSLP':  # make AAres= 10 #otherwise it becomes really uneffecient
#        radius= 300
#        pixrad= radius/(2.5*AAres) #the pixelradius
#        SLP_radmean= CircularFilter_samedist(ds.air_pressure_at_sea_level.values[0]/100, radius= pixrad)
#        diffslp= ds.air_pressure_at_sea_level.values[0]/100 - SLP_radmean
#        diffslp[:int(.5*pixrad)], diffslp[:, :int(.5*pixrad)], diffslp[-int(.5*pixrad):], diffslp[:, -int(.5*pixrad):]= 0, 0, 0, 0  #cout of edges
#    
#        PlotColorMap4(Lon, Lat, diffslp, map, label=r'$\overline{SLP} - SLP$ [hPa] in '+str(radius)+'km radius', color='RdBu', symetric=True)
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        title_extra='_radius'+str(radius)+'km'
#    
#    """filtered SLP"""
#    if var== 'SLP_filtered':
#        T_low, T_up= 60, 120 
#        MSLPfourierfilter= FourierFilter2d_equaldist(ds.air_pressure_at_sea_level.values[0]/100, dist= 2.5* AAres, T_low= T_low, T_up= T_up)
#    
#        PlotColorMap4(Lon, Lat, MSLPfourierfilter, map, label= r"SLP filtered T"+str(T_low)+'-T'+str(T_up) , color='RdBu', symetric=True)
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        title_extra='_T'+str(T_low)+'-T'+str(T_up)
#    
#    
#    """filtered Laplacian of SLP"""
#    if var== 'SLP_laplacian':
#        dx= 2500* AAres
#        SLPlap= np.gradient(np.gradient(ds.air_pressure_at_sea_level.values[0]/100*1E2, dx, axis= 0, edge_order=1), dx, axis= 0, edge_order=1) + np.gradient(np.gradient(ds.air_pressure_at_sea_level.values[0]/100*1E2, dx, axis= 1, edge_order=1), dx, axis= 1, edge_order=1)
#        SLPlap[:3], SLPlap[:, :3], SLPlap[-3:], SLPlap[:, -3:]= 0, 0, 0, 0  #cout of edges
#    
#        import scipy.ndimage.filters as filters
#        gausfilterdist= 100  
#        SLPlapfilter= filters.gaussian_filter(SLPlap, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 4.)
#    
#        PlotColorMap4(Lon, Lat, SLPlapfilter, map, color= 'RdBu', symetric= True, label= r"$\nabla^2$ SLP [Pa/m$^2$] Gauss filter of std "+str(gausfilterdist)+"km")
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        title_extra= '_GaussFilter'+str(gausfilterdist)+'km'
#    
#    
    """calculate vorticity"""
    if var== 'Vort_lev':
        gausfilterdist= 60# 40 #100
        vort= CalcVort(ds.x_wind_pl.values, ds.y_wind_pl.values, dx= 2.5, filter=True, filterdist= gausfilterdist, gaustruncate= 1)
        
        PlotColorMap4(Lon, Lat, vort, map, label= r"$\xi_{"+str(plev)+"}$ [10$^{-5}$ 1/s]" #Gauss filter of std "+str(gausfilterdist)+"km"
        , color='RdBu', symetric=True)
        
        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist, numbers=False)
        U= np.sqrt(ds.x_wind_10m.values**2+ds.y_wind_10m.values**2)
    
        PlotLocalMax(vort, threshold=20, distance=250/(AAres*2.5), map= map, 
                     lon=ds.longitude.values, lat=ds.latitude.values, typ='max', color='black', roundorder= 0, latbound=[66, 83], lonbound= [-15, 20])    
    #    PlotLocalMax(ds.air_pressure_at_sea_level.values[0]/100, threshold=1010, distance=100/AAres, map= map, lon=ds.longitude.values, lat=ds.latitude.values,
    #                 data2=U, threshold2=16, distance2=80/AAres, color='g')      
        title_extra= '_GaussFilter'+str(gausfilterdist)+'km'
#    
#    
#    
#    """calculate shear deformation"""
#    if var== 'Shear_lev':
#        dx= 2500* AAres
#        shear= (np.gradient(d.v, dx, axis= 1)+ np.gradient(d.u, dx, axis= 0))*1E5 #it is not transformation invariant
#    #    shear= np.sqrt(np.gradient(d.v, dx, axis= 1)**2+ np.gradient(d.u, dx, axis= 0)**2)*1E5
#    
#    
#        import scipy.ndimage.filters as filters
#        gausfilterdist= 60
#        shearfilter= filters.gaussian_filter(shear, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#        PlotColorMap4(Lon, Lat, shearfilter, map, color='RdBu', symetric=True,
#                      label= r"$(\partial_y u + \partial_x v)_{"+str(d.chosen_plevel)+"}$ [10$^{-5}$ 1/s]") #Gauss filter of std "+str(gausfilterdist)+"km")
#    #    PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        PlotContours(Lon, Lat, d.geop, map, leveldist= 10)
#        PlotWind(ds.longitude.values, ds.latitude.values, d.u, d.v, map, rot=False)
#    
#        title_extra= '_GaussFilter'+str(gausfilterdist)+'km'
#
#    if var== 'Abs_shear_lev':
#        dx= 2500* AAres
#        shear= (np.gradient(d.v, dx, axis= 1)+ np.gradient(d.u, dx, axis= 0))*1E5 #it is not transformation invariant
#    #    absshear= np.sqrt(np.gradient(d.v, dx, axis= 1)**2+ np.gradient(d.u, dx, axis= 0)**2)*1E5
#    
#    
#        import scipy.ndimage.filters as filters
#        gausfilterdist= 20
#        shearfilter= filters.gaussian_filter(shear, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#        absshearfilter= np.abs(shearfilter)
#        PlotColorMap4(Lon, Lat, absshearfilter, map , color='viridis',
#                      label= r"$(\partial_y u + \partial_x v)_{"+str(d.chosen_plevel)+"}$ [10$^{-5}$ 1/s]")# Gauss filter of std "+str(gausfilterdist)+"km")
#    #    PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        PlotContours(Lon, Lat, d.geop, map, leveldist= 10)
#        PlotWind(ds.longitude.values, ds.latitude.values, d.u, d.v, map, rot=False)
#    
#        title_extra= '_GaussFilter'+str(gausfilterdist)+'km'
#
#
#    
#    
#    """calculate stretching deformation - this looks more like like the shear"""
#    if var== 'Stretching_lev':
#        dx= 2500* AAres
#        stretch= (np.gradient(d.u, dx, axis= 1)- np.gradient(d.v, dx, axis= 0))*1E5
#              
#    #    stretch= (np.gradient(d.u, dx, axis= 0)) *1E5 #- np.gradient(d.v, dx, axis= 0))*1E5 #to look only at the u shear
#    #    stretch= (np.gradient(d.v, dx, axis= 1)) *1E5 #- np.gradient(d.v, dx, axis= 0))*1E5
#    
#        import scipy.ndimage.filters as filters
#        gausfilterdist= 60
#        stretchfilter= filters.gaussian_filter(stretch, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#        PlotColorMap4(Lon, Lat, stretchfilter, map, label= r"$(\partial_x u - \partial_y v)_{"+str(d.chosen_plevel)+"}$ [10$^{-5}$ 1/s] Gauss filter of std "+str(gausfilterdist)+"km", color='RdBu', symetric=True)
#    #    PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        PlotContours(Lon, Lat, d.geop, map, leveldist= 10)
#        PlotWind(ds.longitude.values, ds.latitude.values, d.u, d.v, map, rot=False)
#        
#    #    PlotWind(ds.longitude.values, ds.latitude.values, d.u, np.zeros(d.v.shape), map, rot=False, nx= 50) #to look only at the u shear
#    #    PlotWind(ds.longitude.values, ds.latitude.values, np.zeros(d.v.shape), d.v, map, rot=False, ny= 50)
#           
#        title_extra= '_GaussFilter'+str(gausfilterdist)+'km'
#    
#    """calculate Okubo-Weiss parameter (vort^2 - (shear^2 + stretching^2)"""
#    if var== 'OW_lev':
#        dx= 2500* AAres
#        OW= 1E10* ( (np.gradient(d.v, dx, axis= 1)- np.gradient(d.u, dx, axis= 0))**2 - (np.gradient(d.v, dx, axis= 1)+ np.gradient(d.u, dx, axis= 0))**2 - (np.gradient(d.u, dx, axis= 1)- np.gradient(d.v, dx, axis= 0))**2 )
#        import scipy.ndimage.filters as filters
#        gausfilterdist= 60
#        OWfilter= filters.gaussian_filter(OW, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#        OWfilter[OWfilter < 0]= 0
#    #    PlotColorMap4(Lon, Lat, OWfilter, map, label= r"OW$_{"+str(d.chosen_plevel)+"}$ [10$^{-8}$ 1/s$^2$] Gauss filter of std "+str(gausfilterdist)+"km", color='RdBu', symetric=True)
#        PlotColorMap4(Lon, Lat, np.sqrt(OWfilter), map, label= r"$\sqrt{OW_{"+str(d.chosen_plevel)+"}}$ [10$^{-5}$ 1/s] Gauss filter of std "+str(gausfilterdist)+"km", color='red')
#    
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        title_extra= '_GaussFilter'+str(gausfilterdist)+'km'
#    
#        
#    """the fourier filtered vorticity"""
#    if var== 'Vort_lev_T':
#    #    d.imp_level(pn= 4, tn= t) #850hPa
#        
#        dx= 2500* AAres
#        vort= (np.gradient(d.v, dx, axis= 1)- np.gradient(d.u, dx, axis= 0))*1E5
#        
#        T_low, T_up= 40, 100          
#        vortfourierfilter= FourierFilter2d_equaldist(vort, dist= 2.5* AAres, T_low= T_low, T_up= T_up)
#    
#        PlotColorMap4(Lon, Lat, vortfourierfilter, map, label= r"$\xi_{"+str(d.chosen_plevel)+"}$ [10$^{-5}$ 1/s] filtered T"+str(T_low)+'-T'+str(T_up) , color='RdBu', symetric=True)
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        title_extra='_filtered_T'+str(T_low)+'-T'+str(T_up)
#    
#    """ different geopotential levels and mslp"""
#    if var== 'Geop_lev':
#    #    title_extra=' SLP and '+str(d.chosen_plevel)+'hPa geopotential height'
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= 2)
#        PlotContours(Lon, Lat, d.geop, map, color='r', leveldist= 10)
#    #    PlotContours(Lonm, Latm, geop500, m, color=[(plt.cm.RdBu_r(h)) for h in range(0,257,32)], leveldist= 20)
#    
#    
#    
#    if var== 'test_rot': #the rotation of the wind field
#        del_lon= (ds.longitude.values[1:, :]- ds.longitude.values[:-1, :])* np.cos(np.deg2rad(ds.latitude.values[1:, :]))
#        del_lat= ds.latitude.values[1:, :]- ds.latitude.values[:-1, :]
#        
#        rot= np.rad2deg(np.arctan2(del_lon, del_lat)) #how much to rotate the gridcell anticlockwise to be oriented N-S
#    #    PlotColorMap4(Lon[1:], Lat[1:], rot, map)
#        
#        dir_raw= np.rad2deg(np.arctan2(d.u[1:], d.v[1:])) #the wind direction in respect to the grid (0 being upwards, 90 rightwards, -90 leftw)
#    #    PlotColorMap4(Lon[1:], Lat[1:], dir_raw, map)
#    
#        dir_rot= (dir_raw + rot)%360 #the wind direction with respect to the North pole (0 Northwards, 90 E, 270 W)
#    #    PlotColorMap4(Lon[1:], Lat[1:], dir_rot, map)
#    
#        wind_vel= np.sqrt(d.u[1:]**2 +d.v[1:]**2)
#        u= wind_vel* np.sin(np.deg2rad(dir_rot))
#        v= wind_vel* np.cos(np.deg2rad(dir_rot))
#        
#        PlotWind(ds.longitude.values[1:], ds.latitude.values[1:], u, v, map, rot=True, nx= 50, ny= 50) #this should give the same as only the wind field
#    
#    
#    if var == 'Geop_wind_lev':
#        PlotContours(Lon, Lat, d.geop, map, leveldist= 10, alpha= .7, numbers= False)
#    #    PlotWindVelo(Lon, Lat, np.sqrt(d.u**2+ d.v**2), map)
#        new_map, norm= PlotColorMap4(Lon, Lat, np.sqrt(d.u**2+ d.v**2), map, color='YlBu', label='Wind velocity [m/s] at '+str(d.chosen_plevel)+'hPa', bounds= np.arange(0, 31, 3), return_colormap= True)
#
#
#        uzonal, vmeri= Wind2ZonalMeridional(ds.longitude.values, ds.latitude.values, d.u, d.v) #this rotates the wind in model coordinate direction to zonal and meridional winds
#        PlotWind2(ds.longitude.values, ds.latitude.values, uzonal, vmeri, map, everyx= 13, everyy= 13, rot=False, color='white', quiverkey= False)
#    
##        PlotWind2(ds.longitude.values, ds.latitude.values, d.u, d.v, map, rot=False, everyx= 13, everyy= 13, color='white') #this was a wrong rotation
#    
#    
#    """filtered Laplacian of the geopotential height"""
#    if var== 'Geop_laplacian_lev':
#        dx= 2500* AAres
#        geoplap= np.gradient(np.gradient(d.geop, dx, axis= 0, edge_order=1), dx, axis= 0, edge_order=1) + np.gradient(np.gradient(d.geop, dx, axis= 1, edge_order=1), dx, axis= 1, edge_order=1)
#        geoplap[:3], geoplap[:, :3], geoplap[-3:], geoplap[:, -3:]= 0, 0, 0, 0  #cout of edges
#    
#        import scipy.ndimage.filters as filters
#        gausfilterdist= 40  
#        geoplapfilter= filters.gaussian_filter(geoplap, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 4.)
#    
#        PlotColorMap4(Lon, Lat, geoplapfilter, map, color= 'RdBu', symetric= True, label= r"$\nabla^2$ Z$_{"+str(d.chosen_plevel)+"}$ [m/m$^2$] Gauss filter of std "+str(gausfilterdist)+"km")
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist, numbers=False)
#        PlotContours(Lon, Lat, d.geop, map, color='r', leveldist= 10, numbers=False)
#    
#        title_extra= '_GaussFilter'+str(gausfilterdist)+'km'
#     
#    """geopotential 500-1000 hPa thickness"""
#    if var== 'Thickness':
#        d.imp_geop1000(tn= t)
#        d.imp_level(pn= 8, tn= t)
#    
#        geop_diff= d.geop - d.geop1000
#        PlotContours(Lon, Lat, geop_diff, map)
##        PlotColorMap3(Lon, Lat, geop_diff, map, symetric=False, label='500-1000 hPa thickness')
#        
#        bounds= np.arange(4900, 5150+1, 25)
#        PlotColorMap4(Lon, Lat, geop_diff, map, bounds= bounds, color='RdBu', label='500-1000 hPa thickness')
#
#        title_extra='_500-1000hPa'
#
#
#
    """geopotential 500-1000 hPa thickness + instability + temp gradients"""
    if var== 'Thickness_advanced':
        IceEdge= True
        geop_diff=  1/9.81* (ds.geopotential_pl.sel(pressure= 500).values - ds.geopotential_pl.sel(pressure= 1000).values)
        
        bounds= np.arange(4900, 5250+1, 25)
        PlotColorMap4(Lon, Lat, geop_diff, map, color='RdBu',bounds= bounds, label='500-1000 hPa thickness [m]')

        title_extra='_500-1000hPa'

        """stability"""
        theta_e500= EquiPotTemp(ds.air_temperature_pl.sel(pressure= 500).values, ds.specific_humidity_pl.sel(pressure= 500).values, 500)
        
        theta_eSST= EquiPotTemp(ds.air_temperature_0m.values[0], ds.specific_humidity_2m.values[0], ds.surface_air_pressure.values[0]/100)  
#        PlotStaticStability(Lon, Lat, (theta_eSST- theta_e500) , map, difftype='theta_eSST-theta_e500')

        gausfilterdist= 50
        stab_filter= filters.gaussian_filter((theta_eSST- theta_e500), sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)

        PlotStaticStability(Lon, Lat, stab_filter , map, difftype='theta_eSST-theta_e500', color='g')
#        PlotColorMap4(Lon, Lat, stab_filter, map, symetric= True, maxlevel=25, color='RdBu', label=r"$\theta_{e,SST}-\theta_{e,500}$ [K]")
#
#
        """fronts"""       
        gausfilterdist = 100 #20 in km
        theta= PotTemp(ds.air_temperature_pl.sel(pressure= 850).values, 850)
        theta_filter= filters.gaussian_filter(theta, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
        
#        PlotContours(Lon, Lat, theta_filter, map, leveldist= 2)
        
        dx= 2500* AAres
        grad_theta_filter= np.gradient(theta_filter, dx, edge_order=1)
        abs_grad_theta_filter= np.sqrt(grad_theta_filter[0]**2+ grad_theta_filter[1]**2)

#        bounds=np.arange(0,np.max(abs_grad_theta_filter)*1E5 +1,1)
        levels=[2, 3, 4, 5,6]
        PlotContours(Lon, Lat, abs_grad_theta_filter *1E3*1E2, map, levels= levels)
        
        
#        PlotColorMap4(Lon, Lat, abs_grad_theta_filter *1E3*1E2, map, color= 'red', bounds=bounds, label=r"$\nabla \theta_{"+str(d.chosen_plevel)+"}$ [K/100km]")
#
#                                      
#       
    """ Precip """
    if var== 'Precip': 
        precipbounds= [0, 0.1, 0.2, 0.5, 1, 2, 5, 10]
        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level[0]/100, map, leveldist= 1)
        PlotColorMap4(Lon, Lat, ds.precipitation_amount, map, bounds= precipbounds, color='blue', label='Precipitation [mm/h]')
#
#    if var== 'Precip_acc':
##        precipbounds= [0, 0.1, 0.2, 0.5, 1, 2, 5, 10]
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        PlotColorMap4(Lon, Lat, d.prec_acc, map, color='blue', label='Precipitation in '+str(t)+'h [mm]')
#
#
#    if var== 'LH_reldase': 
##        precipbounds= [0, 0.1, 0.2, 0.5, 1, 2, 5, 10]
#        LHrelease= LatentHeat_fromSnow(d.prec)
#        
#        gausfilterdist= 100
#        LHfilter= filters.gaussian_filter(LHrelease, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        bounds= np.arange(0, 601, 100)
#        PlotColorMap4(Lon, Lat, LHfilter, map,  color='blue', bounds=bounds, label='Latent heat release [W/m$^2$]')
#
#
#    
#    """ Sensible Heat Flux"""
#    if var== 'TH':
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        
#        bounds= np.arange(0, 601, 100)
#        PlotColorMap4(Lon, Lat, d.SenH, map, color='blue', bounds= bounds, label=r"Sensible heat flux [W/m$^2$]")
#
#    if var== 'TH_thor':
#        fluxbounds= np.arange(0, 280, 2)
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        PlotColorMap4(Lon, Lat, d.SenH, map, color='jet', bounds= fluxbounds, label=r"Sensible heat flux [W/m$^2$]")
#       
#    
#    if var== 'LH':
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        
#        #        fluxbounds= [-800,  -500, -300, -200 , -100, -50, 0, 50, 100, 200, 300, 500, 800]
#        bounds= np.arange(0, 601, 100)
#
#        PlotColorMap4(Lon, Lat, d.LatH, map, color='blue', bounds= bounds, label=r"Latent heat flux [W/m$^2$]")
#
#    if var== 'LH_thor':
#        fluxbounds= np.arange(0, 280, 2)
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        PlotColorMap4(Lon, Lat, d.LatH, map, color='jet', bounds= fluxbounds, label=r"Latent heat flux [W/m$^2$]")
#       
#    
#    if var== 'FLX':
#        fluxbounds= [0,  100, 200, 300, 500, 800, 1200]
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        PlotColorMap4(Lon, Lat, d.SenH+ d.LatH, map, bounds= fluxbounds, color= 'brown', label=r"Turbulent heat flux [W/m$^2$]")
#       
#    
#    if var== 'PBH':
#        BLbounds= [0,  500, 1000, 2000, 3000, 4000, 5000]
#        BLbounds= np.arange(0, 5001, 500)       
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= 2, numbers=False)
#
#        PlotColorMap4(Lon, Lat, d.BL, map, bounds=BLbounds,  color= 'viridis_r', label=r"Boundary layer height [m]")
#
##        gausfilterdist= 60
##        BLfilter= filters.gaussian_filter(d.BL, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
##        PlotColorMap4(Lon, Lat, BLfilter, map, bounds=BLbounds,  color= 'viridis_r', label=r"Boundary layer height [m]")
#
#
    if var== 'PBH_advanced':
        BLbounds= np.array([0,  500, 1000, 2000, 3000, 4000, 5000])
        BLbounds= np.arange(0, 5001, 500)       
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= 2, numbers=False)

        PlotColorMap4(Lon, Lat, ds.atmosphere_boundary_layer_thickness[0]/1E3, map, bounds=BLbounds/1E3,  color= 'viridis_r', label=r"Boundary layer height [km]")


        gausfilterdist= 60
        CAPEfilter= filters.gaussian_filter(ds.specific_convective_available_potential_energy_old[0], sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)

        CAPEbounds= np.array([200, 400, 600]) 
        
        PlotContours(Lon, Lat, CAPEfilter, map, levels= CAPEbounds, color='white', linewidth=1.5)
#       
#    
#    if var== 'W_lev':
#    #    d.imp_level(pn= 6, tn= t)
##        wbounds= [-5, -1, -0.5, -0.2, -0.1, 0.1, 0.2, 0.5, 1, 5]
#        wbounds= [-5, -2, -1, -0.5, -0.2, 0.2, 0.5, 1, 2, 5]
#        print('max updraft: ', str(np.max(d.w)))
#        
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        PlotColorMap4(Lon, Lat, d.w, map, bounds=wbounds, color= 'RdBuWhite', label= r"vertical velocity at "+str(d.chosen_plevel)+" [m/s]")
#    
#    if var== 'W_lev_filter':
#        gausfilterdist= 60
#        wfilter= filters.gaussian_filter(d.w, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#
##        wbounds= [-5, -2, -1, -0.5, -0.2, 0.2, 0.5, 1, 2, 5]
##        print('max updraft: ', str(np.max(d.w)))
#        
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        PlotColorMap4(Lon, Lat, wfilter, map, symetric= True, color= 'RdBuWhite', label= r"vertical velocity at "+str(d.chosen_plevel)+" [m/s]")
#    
#    
#    if var== 'PV_lev':
##        PVbounds= np.array([-5, -3, -2, -1, -0.5, 0,  0.5, 1, 2, 3, 5])
#        PVbounds= np.array([0, 1, 2, 3, 5])
#
##        if d.chosen_plevel < 500: PVbounds= np.array([-10, -5, -3, -2, 0, 2, 3, 5, 10])
#        
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        PlotColorMap4(Lon, Lat, d.PV *1E6, map, color= 'red', bounds=PVbounds, label= str(d.chosen_plevel) + r"hPa PV [$10^{-6}$K m$^2$ kg$^{-1}$ s$^{-1}$]")
#              
#    
#    if var== 'RH_lev':
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        new_map, norm= PlotColorMap4(Lon, Lat, d.RH, map, symetric=False, label='Relative humidity at ' +str(d.chosen_plevel)+'hPa', color='blue', bounds= np.array([0, 0.4, 0.7, 0.8, 0.9, 0.95, 1]), return_colormap= True )
#    
#    if var== 'SH_lev':
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist, numbers= False)
#        new_map, norm= PlotColorMap4(Lon, Lat, d.SH*1E3, map, symetric=False, label='Specific humidity [g/kg] at ' +str(d.chosen_plevel)+'hPa', return_colormap= True
#                                     , bounds= np.arange(0, np.max(d.SH)*1E4, 2)*1E-1, color='jet_r')
#
#        RH= SH2RH(d.SH*1E3, d.chosen_plevel, d.T)/100
#        
#        gausfilterdist= 10
#        RH= filters.gaussian_filter(RH, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)            
#        
#        PlotContours(Lon, Lat, RH, map, levels= [0.95], numbers= False, color='white')
#
#    
#    
#    if var== 'RH_lev_fromSH':
#        import metpy.calc as mpcalc
#        from metpy.units import units
##        RH= mpcalc.relative_humidity_from_specific_humidity(d.SH, d.T* units.K, d.chosen_plevel* units.hPa)
#        RH= SH2RH(d.SH*1E3, d.chosen_plevel, d.T)/100
#  
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist, numbers= False)
#        new_map, norm= PlotColorMap4(Lon, Lat, RH, map, symetric=False, label='Relative humidity at ' +str(d.chosen_plevel)+'hPa', color='blue', bounds= np.array([0, 0.4, 0.7, 0.8, 0.9, 0.95, 1]), return_colormap= True )
#    
#    
#    if var== 'Cloud_conv':
#    #    Cloud_convbounds= np.array([0,  100,  200, 300, 500, 800])
#        Cloud_bounds= np.array([0, 0.2, 0.4, 0.6, 0.8, 1])
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        PlotColorMap4(Lon, Lat, d.cloud_conv , map, bounds=Cloud_bounds, color= 'brown', label=r"Convective Clouds")
#    
#    
    if var== 'Cloud_high':
        Cloud_bounds= np.array([0, 0.2, 0.4, 0.6, 0.8, 1])    
        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist, color='r')
        PlotColorMap4(Lon, Lat, ds.high_type_cloud_area_fraction , map, bounds=Cloud_bounds, color= 'grey', label=r"High cloud cover")

    if var== 'Cloud_med':
        Cloud_bounds= np.array([0, 0.2, 0.4, 0.6, 0.8, 1])    
        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist, color='r')
        PlotColorMap4(Lon, Lat, ds.medium_type_cloud_area_fraction , map, bounds=Cloud_bounds, color= 'grey', label=r"Medium cloud cover")        

    if var== 'Cloud_low':
        Cloud_bounds= np.array([0, 0.2, 0.4, 0.6, 0.8, 1])    
        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist, color='r')
        PlotColorMap4(Lon, Lat, ds.low_type_cloud_area_fraction , map, bounds=Cloud_bounds, color= 'grey', label=r"Low cloud cover")
#    
#    if var== 'CAPE':
#        CAPEbounds= np.array([0,  100,  200, 300, 500, 800])    
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        PlotColorMap4(Lon, Lat, d.CAPE , map, bounds=CAPEbounds, color= 'brown', label=r"CAPE [J/kg]")
#    
#    if var== 'CIN':
#        CINbounds= np.array([-100, -10, -1, -0.1, 0])    
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        PlotColorMap4(Lon, Lat, d.CIN , map, bounds=CINbounds, color= 'brown', label=r"CIN [J/kg]")
#     
#        
#    if var== 'LFC':
#        LFCbounds= np.array([0, 10, 100, 200, 500, 1000, 10000, 30000])    
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        PlotColorMap4(Lon, Lat, d.LFC , map, bounds=LFCbounds, color= 'brown', label=r"LFC [m]")
#    
#    if var== 'LCL':
#        LCLbounds= np.array([0, 10, 100, 1000, 3000, 5000, 10000, 30000])  
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        PlotColorMap4(Lon, Lat, d.LCL , map, bounds=LCLbounds, color= 'brown', label=r"LCL [m]")
#    
#    if var== 'CTT':
#        levels= np.arange(940, 1040, 3)
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, levels= levels, numbers=False, color= 'r', alpha= .6)
#        bounds= np.arange(210, 281, 2)   
#    
#        PlotColorMap4(Lon, Lat, d.CTT , map, color= 'grey', bounds= bounds, label=r"Cloud top temperature [K]")
#    #    PlotColorMap4(Lon, Lat, d.CTT , map, label=r"Cloud top temperature [K]", boxnr= 50)
#    
#    if var== 'WVTC':
#    #    Cloud_bounds= np.array([0, 0.2, 0.4, 0.6, 0.8, 1])    
#    #    PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        PlotColorMap4(Lon, Lat, d.WVTC , map, color= 'grey', label=r"Water vapor temperature with clouds [K]", boxnr= 128)
#    
#    if var== 'WVT':
#    #    Cloud_bounds= np.array([0, 0.2, 0.4, 0.6, 0.8, 1])    
#    #    PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        PlotColorMap4(Lon, Lat, d.WVT , map, color= 'grey', label=r"Water vapor temperature [K]", boxnr= 100)
#    
#    if var== 'CWR':
#    #    Cloud_bounds= np.array([0, 0.2, 0.4, 0.6, 0.8, 1])    
#    #    PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, leveldist= pleveldist)
#        PlotColorMap4(Lon, Lat, d.CWR , map, color= 'grey', label=r"Cloud water reflectivity", bounds= np.linspace(0,1., 101))
#
    if var =='OLW':
        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[t, 0]/100, map, leveldist= pleveldist*2, numbers=False, color='r')

        OLW= -1/(60**2)*(ds.integral_of_toa_outgoing_longwave_flux_wrt_time.values[t+1,0]- ds.integral_of_toa_outgoing_longwave_flux_wrt_time.values[t,0])
        PlotColorMap4(Lon, Lat, OLW , map, label= 'Outgoing longwave radiation [W/m$^2$]', color='grey', bounds= np.arange(130, 251, 10))

    if var =='netrad':
#        PlotContours(Lon, Lat, ds.air_pressure_at_sea_level.values[0]/100, map, levels= levels, numbers=False, color= 'r', alpha= .6)
        PlotColorMap4(Lon, Lat, ds.RN.values , map, label= 'Net radiation [W/m$^2$]', color='grey_r', bounds= np.arange(-150, 30, 10))

#
#    if 'tropopause' in var:
#        Lapserate= 1E3* (d.T[1:] - d.T[:-1]) /d.dz     
#        zintermediate= (d.z[1:] + d.z[:-1])/2
#                
#        """works but extremly unefficient"""
##        tropheight, troppres= np.zeros(np.shape(Lapserate)[1:]), np.zeros(np.shape(Lapserate)[1:])
##        for xi in range(Lapserate.shape[1]):
##            print(xi)
##            for yi in range(Lapserate.shape[2]):
##                
##                tropheight[xi, yi]= next(zintermediate[i,xi, yi]
##                    for i,l in enumerate(Lapserate[:,xi, yi]) if  (l >-2 
##                        and np.min(Lapserate[i: np.argmin(np.abs(zintermediate[:, xi, yi] - zintermediate[i, xi, yi] - 2000)), xi, yi]) > -2 ) )
## 
#        """this is precises"""
#        zdist= 200
#        z_int= np.arange(3000, 13001, zdist) #this defines the range in which the tropopause height can be
#        from scipy.interpolate import interp1d        
#
#        Lapserate_int= np.moveaxis( np.array([[ interp1d(zintermediate[:, x, y], Lapserate[:, x, y])(z_int) for y in range(Lapserate.shape[2]) ] for x in range(Lapserate.shape[1]) ]  ), -1, 0)
#        p_int= np.moveaxis( np.array([[ interp1d(zintermediate[:, x, y], (d.pres[1:,x,y]+ d.pres[:-1,x,y])/2)(z_int) for y in range(Lapserate.shape[2]) ] for x in range(Lapserate.shape[1]) ]  ), -1, 0)
#
#        
#        levels2km= int(2000/zdist)
#        Lapserate_2kmabove= np.array([np.min(Lapserate_int[i: i+levels2km], axis= 0) for i in range(len(z_int) -levels2km)])
#        
#        z_int_total= z_int[:, np.newaxis, np.newaxis]* np.ones(np.shape(Lapserate_int) )
#        z_int_total[Lapserate_int < -2]= np.nan
#        z_int_total[:-levels2km][Lapserate_2kmabove < -2]= np.nan
#        tropheight= np.nanmin(z_int_total, axis= 0)
#
#
#    if var == 'tropopause':
#        if display=='pres':
#            p_tropheight= np.array([[p_int[:, x,y][z_int== tropheight[x,y]][0] for y in range(Lapserate.shape[2]) ] for x in range(Lapserate.shape[1]) ]  )
#            
#            gausfilterdist= 100
#            p_tropheight= filters.gaussian_filter(p_tropheight, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)            
#            
#            bounds= np.arange(260, 361, 10)
#            PlotColorMap4(Lon, Lat, p_tropheight, map, bounds= bounds, color='viridis_r', label='Tropopause height [hPa]')
#            PlotContours(Lon, Lat, p_tropheight, map, leveldist= 10)
#
#        else: #display='z':
#            gausfilterdist= 100
#            tropheight= filters.gaussian_filter(tropheight, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#
#            bounds= np.arange(7000, 9001, 250)
#            PlotColorMap4(Lon, Lat, tropheight, map, bounds= bounds, label='Tropopause height [m]')
#            PlotContours(Lon, Lat, tropheight, map, leveldist= 250)
#
#
#    if var== 'tropopause_advanced':
#        gausfilterdist= 100
#
#        tropheight= filters.gaussian_filter(tropheight, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#
#        bounds= np.arange(7000, 9251, 200)
#        PlotColorMap4(Lon, Lat, tropheight/1E3, map, bounds= bounds/1E3, label='Tropopause height [km]')#        PlotContours(Lon, Lat, p_tropheight, map, leveldist= 10)
#
#
#        d2= data(filename= AAfilename, res= AAres)
#        d2.imp_level(pn= 10, tn= t)
#
#        PlotContours(Lon, Lat, d2.geop/1E3, map, leveldist= .02, alpha= 1, numbers= True, color='black', fmt='%1.2f')
#    #    PlotWindVelo(Lon, Lat, np.sqrt(d.u**2+ d.v**2), map)
#
#        PlotWind(d2.lon, d2.lat, d2.u, d2.v, map, rot=False, alen= 40) #put this after PlotColorMap, otherwise it disappears
#
#        
##        plt.figure(fignr+1)
##        plt.clf()
##        plt.plot(Lapserate_int[:, 0, 0], p_int[:, 0, 0])
#
#    if var== 'dynamical_trop':
#    
#        from scipy.interpolate import interp1d        
#        
#        inc= 20
#        pres_int= np.arange(100, 700, inc)
#            
#        PV_int= np.moveaxis( np.array([[ interp1d(d.pres, d.PV[:, x, y])(pres_int) for y in range(d.PV.shape[2]) ] for x in range(d.PV.shape[1]) ]  ), -1, 0)
#        
#        pres_int_total= pres_int[:, np.newaxis, np.newaxis]* np.ones(np.shape(PV_int) )
#        pres_int_total[PV_int < 2*1E-6]= np.nan
#        tropheight= np.nanmax(pres_int_total, axis= 0) +inc/2
#        
#        gausfilterdist= 20
#        tropheight= filters.gaussian_filter(tropheight, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#        
#        bounds= np.arange(200, 561, 20)
#        PlotColorMap4(Lon, Lat, tropheight, map, bounds= bounds, color='viridis_r', label='Tropopause height [hPa]')
#

    if var=='surface_geop':
        PlotColorMap4(Lon, Lat, 1/9.81* ds.surface_geopotential.values[0], map, bounds= [-20, 0, 20, 40, 60])


 
#    """Plot Thorpex points"""
#    if Thorpex== True:
#        if year == 2008 and month== 3 and day== 3 and hour<15:
#            excl= [3, 5, 12] #for the first flight       
#        elif year == 2008 and month== 3 and day== 3 and hour>=15:
##            excl=[1,9,11] #for the second flight
#            excl=[1, 2, 12, 14] #1+14 just because it is outside the domain. 11 very little data around 500hPa. no idea what happens with the dropsonde between 4 and 5, which is 5 in other papers
#        elif year == 2008 and month== 3 and day== 4:
#            if 'DA' in exp_name:
#                excl= [1, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]# for the third flight - for original AA domain
#            else:
#                excl= [1, 5, 13]# for the third flight - for AA domain moved southwards
#        
#        thor= thorpex.data(year, month, day, hour, plevels= [d.chosen_plevel], exclude=excl, excl_level= 500)
#    #    thor2= thorpex.data(year, month, day, hour, level= d.chosen_plevel)
#        
#        xpt, ypt= map(thor.lon, thor.lat)
#        #map.plot(xpt,ypt,'bo')
#        
#        for i in np.arange(len(xpt)):
##            if var=='Geop_wind_lev': scattercolor= 'r'
##            if var is not 'Geop_wind_lev':
#            plt.text(xpt[i]-10000, ypt[i]+20000, str(thor.dropnr[i]), color=scattercolor, fontsize=14, fontweight='bold')               
##                plt.text(xpt[i]+100, ypt[i]+100, str(thor.dropnr[i])+'\n'+ str(thor.datetime[i, 0])[11:16], color='r')
#                    
#            
#            if var== 'Theta_lev' or var=='Theta_newbounds_lev':
#                value= thor.theta[i]
#            elif var== 'Theta_e_lev':
#                value= thor.theta_e[i]                   
#            elif 'RH_lev' in var:
#                value= thor.RH[i]
#            elif var== 'SH_lev':
#                value= thor.SH[i]*1E3
#            elif var== 'Geop_wind_lev':
#                value= thor.U[i]                
#            else:
#                break
#
#            if np.isnan(value[0]): print('Exclude ', thor.dropnr[i])   #in case the value that should be plotted is none it exits the loop
#            else:
#               
#                map.scatter(xpt[i], ypt[i], c= value , cmap= new_map, norm= norm, s= 150, edgecolors= scattercolor)
#        
#        
#                if var== 'Geop_wind_lev':
#                    dist= np.sqrt( (thor.lat[i,0]- ds.latitude.values)**2+ (np.cos(np.deg2rad(thor.lat[i,0]))* (thor.lon[i,0]- ds.longitude.values))**2)
#                    xpos, ypos= np.where(dist== np.min(dist))
##                    plt.text(xpt[i], ypt[i]+10000, str(int(thor.alt[i,0]- d.geop[xpos[0], ypos[0]])), color='r')  #depature from the geopotential height  
#                    
#                    print(thor.dropnr[i], str(int(thor.alt[i,0]- d.geop[xpos[0], ypos[0]])))        
#                    u,v, Lon, Lat = map.rotate_vector(thor.u[i], thor.v[i] , thor.lon[i], thor.lat[i] ,returnxy=True)
#                    plt.quiver(Lon, Lat, u, v, color= 'r', scale = 500, width= .004, headwidth= 5)
#    
#    
#    
#    if Plotbox== True and '_box' not in title_extra:
#        x,y= map(boxlon, boxlat) #longitude and latitude coordinates of the corner points
#        map.plot(x, y, color=boxcolor, linewidth= 2)
#        title_extra +='_box'
#
#
    if IceEdge== True:
        Ice = ds.air_temperature_0m.values[0] -273.15
        Ice[1/9.81* ds.surface_geopotential.values[0] > 20]= 10 #mask land values
        map.contour(Lon, Lat, Ice, '-', levels= [-2], linewidths= 2, colors='white')
#
#
    if PlotCentre== True:
        threshold= 20 #10
        propdistance= 180E3 #in km   #70
        distance= int(200E3/(AAres * 2.5E3)) #100   #make distance reasonable larger than propdistance to avoid "doublematches" twice is on the safe side
        
        
        latbound= [60, 71]
        lonbound= [-5, 15]
#        tbound= [20, 30]    

#        gausfilterdist= 100 #60
#        gaustruncate= 2
#
#        #"""calculate vorticity"""
#        dx= 2500* AAres


        ds_full= xr.open_dataset(AROMEdir + filename)
        ds_full= ds_full.sel(pressure= plev)
        
#        gausfilterdist= 60# 40 #100
        vort_full= CalcVort(ds_full.x_wind_pl.values, ds_full.y_wind_pl.values, dx= 2.5, filter=True, filterdist= gausfilterdist, gaustruncate= 1)

        
#        vort_full= (np.gradient(d.v_ft, dx, axis= -1)- np.gradient(d.u_ft, dx, axis= -2))*1E5     
        
#        vortfilter_full= filters.gaussian_filter1d(vort_full, axis= -1, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= gaustruncate)
#        vortfilter_full= filters.gaussian_filter1d(vortfilter_full, axis= -2, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= gaustruncate)

        #""" a tracking algorithm"""
        vortlist, lonlist, latlist, labellist = EasyTrack(vort_full, ds.longitude.values, ds.latitude.values, distance, threshold, propdistance)
        vortlist, lonlist, latlist, labellist= Tracks_in_lat_lon(vortlist, lonlist, latlist, labellist, latbound, lonbound, tbound= None)

        #"""plot some of the data in map"""
#        ind= np.where(tcycl == t)[0]
#        if len(ind)> 0: #only if the cyclone exists for this point in time

        xpt, ypt= map(lonlist[t], latlist[t])
#        map.plot(xpt,ypt, 'o', color= 'k', markersize= 10)
        
        for i in range(len(xpt)):
            plt.text(xpt[i], ypt[i]+ 0.01*(map.ymax-map.ymin), str(int(np.round(vortlist[t][i]))), color='r', fontsize=15, fontweight= 'bold')


        for cnr in labellist[t]:
            tcycl, vortcycl, loncycl, latcycl= Data_One_Track(vortlist, lonlist, latlist, labellist, cnr)
            xpt, ypt= map(loncycl, latcycl)
            map.plot(xpt,ypt,'kx')
            
#            for i in range(len(loncycl)):
#                value= vortcycl[i]            

    
    if save== False:
        plt.title('Arome '+ ds.min_time[:10]+' + '+str(t).zfill(2))
    
#    if vertical==False and save==True:
    if save==True:
        if '_lev' in var: savevar = var+ str(plev)
        else: savevar= var
            
        savename= savedir+ 'Arome_'+ ds.min_time[:10]+'+'+str(t).zfill(2)+'_'+ savevar +title_extra
        
        plt.savefig(savename, bbox_inches= 'tight')#, dpi= 60)
        print(savename)        
    
    plt.tight_layout()
#    
#
#
#
#
#
#
#
#    
#if vertical==True:
#    """for cross section definition"""
#    
#    x_start= thor.lon[thor.dropnr==startdrop, 0]
#    x_end= thor.lon[thor.dropnr==enddrop, 0]
#    y_start=thor.lat[thor.dropnr==startdrop, 0]
#    y_end= thor.lat[thor.dropnr==enddrop, 0]
#    
#       
#    x, y= mapcross_section(ds.longitude.values, ds.latitude.values, x_start, y_start, x_end, y_end, map, extent= 15, coordinates='latlon', color=scattercolor, markersize= 1.5)
#
#
#    if save==True:
#        if '_lev' in var:
#            savevar = var+ str(d.chosen_plevel)
#        else: savevar= var
#            
#        print('exp_name:', exp_name)
#        if exp_name== 'DA_080301_cycling': save_exp_name= exp_name+'_'+str(fileday).zfill(2)+str(filehour).zfill(2)
#        else: save_exp_name= exp_name
#        
#        savename= savedir+'Arome_'+save_exp_name+'_'+str(t).zfill(2)+'_'+ savevar +title_extra
#        if vertical== True: savename += 'vertical'
#        plt.savefig(savename)
#        print(savename)
#        
#
#    if var== 'W_lev': d.imp_cross_sec(xn= x, yn= y, tn= t)
#    else: d.imp_cross_sec_grib_2(xn= x, yn= y, tn= t)
#
#
#    plt.figure(fignr)
#    plt.clf()
#    fignr += 1
#    
#    upperlevel= 300 #upper level in hPa
#    ilevel= np.where(d.pres == np.max(d.pres[d.pres < upperlevel]))[0][0] +1
#    if var== 'Theta_lev':
#
##        pressm= np.tile(d.pres, (np.shape(d.T)[1],1)).T #make a matrix with the pressure in every row
#        pressm= d.pres
#        theta= PotTemp(d.T, pressm)
#        
#        PlotCross_sec_contour(d.pres[:ilevel, 0], ds.longitude.values[x, y], theta[:ilevel], levels= np.arange(260, 290))#, color='default')
#        new_map, norm =PlotCross_sec_color(d.pres[:ilevel, 0], ds.longitude.values[x, y], theta[:ilevel], label='Potential temperature [K]', color='RdBu')#, color='default')
#
#
#    if var== 'Geop_wind_lev':
##        PCross_sec_horvel(d.pres, ds.longitude.values[y,x], d.u, d.v)
##        PCross_sec_horvel2(d.pres, ds.longitude.values[y,x], d.u, d.v)
#        U= np.sqrt(d.u**2+ d.v**2)
#        new_map, norm =PlotCross_sec_color(d.pres[:ilevel, 0], ds.longitude.values[x, y], U[:ilevel], color='YlBu', label= 'Horizontal wind [m/s]')
##        theta_e= EquiPotTemp(d.T, d.SH, d.pres)
##        PlotCross_sec_contour(d.pres[:ilevel, 0], ds.longitude.values[x, y], theta_e[:ilevel], levels= np.arange(260, 330, 2), numbers=False)
#
#        #to plot the geopotential height, but its boring, rather plot the depature, but this is not implemented yet.
##        d.imp_cross_sec(xn= x, yn= y, tn= t)
##        PlotCross_sec_contour(d.pres[:ilevel], ds.longitude.values[x, y], d.geop[:ilevel])#, color=None, cmap='RdBu_r')
##        d.geop[:ilevel] -np.mean(d.geop[:ilevel], axis= 1)
#
#
#        d.imp_cross_sec(xn= x, yn= y, tn= t)
#        everynth= 1
#        PlotCross_sec_vervel(np.delete(d.pres, [1,2,4]), ds.longitude.values[x,y][::everynth], np.delete(d.w, [1,2,4], axis= 0)[:, ::everynth], color='k')
#
#    if var== 'W_lev':
#        everynth= 3
#        PlotCross_sec_vervel(d.pres, ds.longitude.values[x,y][::everynth], d.w[:, ::everynth])
##        PlotCross_sec_color(d.pres[:ilevel], ds.longitude.values[x, y], d.w[:ilevel])#, color='default')
#
#    if var== 'Theta_e_lev':
#        print('have to transform')
#        theta_e= EquiPotTemp(d.T, d.SH, d.pres)
#        
#        PlotCross_sec_contour(d.pres[:ilevel, 0], ds.longitude.values[x, y], theta_e[:ilevel], levels= np.arange(260, 330, 2), numbers= False)#, color='default')
#        new_map, norm =PlotCross_sec_color(d.pres[:ilevel, 0], ds.longitude.values[x, y], theta_e[:ilevel], label= 'Equivalent potential temperature [K]')#, color='default')
#
#
##        RH= SH2RH(d.SH*1E3, d.pres, d.T)/100
##        PlotCross_sec_contour(d.pres[:ilevel, 0], ds.longitude.values[x, y], RH[:ilevel], levels= [0, 0.8])#, color='default')
#
#
##        d.imp_cross_sec(xn= x, yn= y, tn= t)
##        PlotCross_sec_vervel(d.pres, ds.longitude.values[x,y], d.w, color='r')
#
#
#    if var == 'RH_lev_fromSH':
#        RH= SH2RH(d.SH*1E3, d.pres, d.T)/100
#        
##        PlotCross_sec_hum(d.pres[:ilevel], ds.longitude.values[y,x], d.RH[:ilevel])
#        new_map, norm =PlotCross_sec_color(d.pres[:ilevel, 0], ds.longitude.values[x, y], RH[:ilevel], color='blue', label='Relative humidity', bounds= np.arange(0.1, 1.05, .1))
#
#        theta_e= EquiPotTemp(d.T, d.SH, d.pres)
#        PlotCross_sec_contour(d.pres[:ilevel, 0], ds.longitude.values[x, y], theta_e[:ilevel], levels= np.arange(260, 330, 2), numbers=False)
#
##        d.imp_cross_sec(xn= x, yn= y, tn= t)
##        everynth= 5
##        PlotCross_sec_vervel(d.pres, ds.longitude.values[x,y][::everynth], d.w[:, ::everynth], color='y')
#
#    plt.ylim([1000, upperlevel])
#    
#    
#    
#    thorvert= thorpex.data(year, month, day, hour, level='vertical', plevels= np.linspace(985, upperlevel, 25))
#    
#    istartdrop = np.where(thorvert.dropnr == startdrop)[0][0]
#    ienddrop = np.where(thorvert.dropnr == enddrop)[0][0]
#    
#    for idrop in np.arange(istartdrop, ienddrop+1):
#        if var== 'Theta_lev':  thorvalue= thorvert.theta[idrop]
#        elif var== 'Theta_e_lev': thorvalue=  thorvert.theta_e[idrop] 
#        elif var== 'RH_lev_fromSH':  thorvalue= thorvert.RH[idrop]
#        elif var== 'SH_lev':  thorvalue= thorvert.SH[idrop] 
#        elif var== 'Geop_wind_lev':  thorvalue=  thorvert.U[idrop]
# 
#        plt.scatter(thorvert.lon[idrop], thorvert.pres[idrop], c= thorvalue , cmap= new_map, norm= norm, s= 100, edgecolors= scattercolor)
#
#       
#        plt.text(thorvert.lon[idrop, ~np.isnan(thorvert.lon[idrop])][-1] -.2, upperlevel -10, str(thorvert.dropnr[idrop]), color=scattercolor, fontsize=14)               
#
#        plt.tight_layout()
#    if save:
#        savename= savedir+'Arome_'+save_exp_name+'_'+str(t).zfill(2)+'_'+ var + '_vertical_Drop'+str(startdrop)+'-'+str(enddrop)
#        plt.savefig(savename, pad_inches=0) #, bbox_inches= 'tight' )
#        print(savename)
    