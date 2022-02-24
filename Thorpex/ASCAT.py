#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 19:13:32 2019

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'

import numpy as np
import xarray as xr #for ASR

import sys
sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')

from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_meteo import *
from f_useful import *


"""global variables"""
fignr= 1

#maptype='AA'
maptype='AA_half'
#maptype='Lambert'
#maptype='polar'

if maptype=='AA': plt.figure(fignr, figsize= (6, 6))
elif maptype=='AA_half': plt.figure(fignr, figsize= (6, 4.7))
else: plt.figure(fignr, figsize= (6, 4.5))
plt.clf()

"""time specification: time of the plot"""
year, month = 2008, 3
#day, hour= 4, 1


var= 'Wind_advanced'

save= True
#savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/PseudoSat3/'


"""Lambert coordinates"""

#to match with Thorpex
if maptype=='Lambert':
    if year == 2008 and month== 3 and day== 3 and hour<15:
        lllon, lllat, urlon, urlat= -4, 70, 19, 75
        lat0, lon0= 75, 0 #(lat0, lon0) = center point
    
    elif year == 2008 and month== 3 and day== 3 and hour>=15:
        lllon, lllat, urlon, urlat= -6, 69.5, 13, 74.5
        lat0, lon0= 75, 0 #(lat0, lon0) = center point
    
    elif year == 2008 and month== 3 and day== 4: 
        if var=='CTT':
            lllon, lllat, urlon, urlat= -13, 69, 12, 73
            lat0, lon0= 70, -20 #(lat0, lon0) = center point 
            title_extra +='_zoom'
        elif 'DA' in exp_name: #for original AA domain
            
            lllon, lllat, urlon, urlat= -0.5, 67.5, 12.5, 68.8
            lat0, lon0= 68, -25 #(lat0, lon0) = center point 
           
        else:  #for AA domain moved south
            lllon, lllat, urlon, urlat= -3, 63.5, 15, 69.5
            lat0, lon0= 70, 0 #(lat0, lon0) = center point 


if maptype== 'AA': map = AA_map()
elif maptype== 'AA_half': map = AA_map_half()
elif maptype== 'polar': map=Polar_map(latsouth= 10, hemisp='NH')
#elif maptype== 'polar': map=Polar_map(latsouth= -10, hemisp='SH')
#else: map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0, latdist=2, londist= 5)    
#map=Lat_Lon_map(0, -70, 360, 70)




#ASCAT
file= Mediadir+ '/PL/ASCAT/1319932/'+(
#'OR1ASWC12_20080303_011622_07108_M02.nc')
#'OR1ASWC12_20080303_025745_07109_M02.nc')
#'OR1ASWC12_20080303_043905_07110_M02.nc')
#'OR1ASWC12_20080303_062028_07111_M02.nc')
#'OR1ASWC12_20080303_080148_07112_M02.nc')
#'OR1ASWC12_20080303_094311_07113_M02.nc')
#'OR1ASWC12_20080303_112431_07114_M02.nc') #ok
#'OR1ASWC12_20080303_130552_07115_M02.nc') #ok
#'OR1ASWC12_20080303_144715_07116_M02.nc') #ok
#'OR1ASWC12_20080303_162835_07117_M02.nc') #good
#'OR1ASWC12_20080303_180958_07118_M02.nc') #ok
#'OR1ASWC12_20080303_195118_07119_M02.nc') # good
#'OR1ASWC12_20080303_213241_07120_M02.nc')
#'OR1ASWC12_20080303_231401_07121_M02.nc')
#'OR1ASWC12_20080304_005524_07122_M02.nc')
#'OR1ASWC12_20080304_023645_07123_M02.nc')
#'OR1ASWC12_20080304_041807_07124_M02.nc')
#'OR1ASWC12_20080304_055928_07125_M02.nc')
#'OR1ASWC12_20080304_074050_07126_M02.nc')
'OR1ASWC12_20080304_092211_07127_M02.nc') #ok
#'OR1ASWC12_20080304_110331_07128_M02.nc') #ok
#'OR1ASWC12_20080304_124454_07129_M02.nc')
#'OR1ASWC12_20080304_142615_07130_M02.nc')
#'OR1ASWC12_20080304_160737_07131_M02.nc')
#'OR1ASWC12_20080304_174858_07132_M02.nc') #ok
#'OR1ASWC12_20080304_193020_07133_M02.nc')
#'OR1ASWC12_20080304_211141_07134_M02.nc')
#'OR1ASWC12_20080304_225303_07135_M02.nc')

SCAT= xr.open_dataset(file)
#
Lon,Lat = map(SCAT.lon.values,SCAT.lat.values)
#
#
##PlotColorMap4(Lon, Lat, SCAT.ice_age.values, map)
##PlotColorMap4(Lon, Lat, SCAT.wvc_index.values, map)
#
Wind_speed= np.copy(SCAT.wind_speed.values)
Wind_speed[np.logical_or(SCAT.wvc_index.values ==  41 , SCAT.wvc_index.values ==  42)]= np.nan
#
#PlotColorMap4(Lon, Lat, Wind_speed, map, bounds=np.arange(0,20, 0.5))
PlotWindVelo(Lon, Lat, Wind_speed, map, Umax= 25, color='YlBu')    


#if var=='Wind_advanced':
#    PlotColorMap4(Lon, Lat, SCAT.wind_speed.values, map, bounds=np.arange(0,20, 0.5))

#    PlotWindVelo(Lon, Lat, SCAT.wind_speed.values, map, Umax= 25, color='YlBu')    
#    find local pressure min and plot them in the map   
#    PlotLocalMax(d.mslp, threshold=1010, distance=100/AAres, map= map, lon=d.lon, lat=d.lat,
#                 data2=U, threshold2=18, distance2=150/AAres)    
#  
#    PlotLocalMax(U, threshold=20, distance=100/AAres, map= map, lon=d.lon, lat=d.lat, typ='max',
#                 color='orange', dot=False, roundorder= 0)    
