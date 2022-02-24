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


var= 'Wind_speed'

save= False
savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/SCAT_2/'


SCAT_type= 'QUIKSCAT'
#SCAT_type= 'ASCAT'


if maptype=='AA': plt.figure(fignr, figsize= (6, 6))
elif maptype=='AA_half': plt.figure(fignr, figsize= (6, 4.7))
else: plt.figure(fignr, figsize= (6, 4.5))
plt.clf()


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





"""QUIKSCAT"""
if SCAT_type == 'QUIKSCAT':
    filedir= Mediadir+ '/PL/SCAT/QUIKSCAT/'
    file= (
    #'OR1SWW025_20080303_004214_45328_QUIKSCAT.nc')
    #'OR1SWW025_20080303_022313_45329_QUIKSCAT.nc')
    #'OR1SWW025_20080303_040412_45330_QUIKSCAT.nc')
    #'OR1SWW025_20080303_054511_45331_QUIKSCAT.nc')
    #'OR1SWW025_20080303_072610_45332_QUIKSCAT.nc')
    #'OR1SWW025_20080303_090709_45333_QUIKSCAT.nc')
    #'OR1SWW025_20080303_104808_45334_QUIKSCAT.nc')
    #'OR1SWW025_20080303_122907_45335_QUIKSCAT.nc')
    #'OR1SWW025_20080303_141006_45336_QUIKSCAT.nc')
    #'OR1SWW025_20080303_155106_45337_QUIKSCAT.nc')
#    'OR1SWW025_20080303_173205_45338_QUIKSCAT.nc') #good
#    'OR1SWW025_20080303_191304_45339_QUIKSCAT.nc') #very good
    #'OR1SWW025_20080303_205403_45340_QUIKSCAT.nc') #ok
    #'OR1SWW025_20080303_223502_45341_QUIKSCAT.nc') #ok
    #'OR1SWW025_20080304_001601_45342_QUIKSCAT.nc') #ok
#    'OR1SWW025_20080304_015700_45343_QUIKSCAT.nc') #ok+
    'OR1SWW025_20080304_033759_45344_QUIKSCAT.nc') #good - in major phase
#    'OR1SWW025_20080304_051858_45345_QUIKSCAT.nc') #ok
#    'OR1SWW025_20080304_065957_45346_QUIKSCAT.nc')
    #'OR1SWW025_20080304_084056_45347_QUIKSCAT.nc')
    #'OR1SWW025_20080304_102155_45348_QUIKSCAT.nc')
    #'OR1SWW025_20080304_120254_45349_QUIKSCAT.nc')
    #'OR1SWW025_20080304_134354_45350_QUIKSCAT.nc')
    #'OR1SWW025_20080304_152453_45351_QUIKSCAT.nc')
#    'OR1SWW025_20080304_170552_45352_QUIKSCAT.nc') #ok - PL out of domain
    
#    time= file[10:25]

elif SCAT_type== 'ASCAT':
    filedir= Mediadir+ '/PL/SCAT/ASCAT/'
    file=(
    #'OR1ASWC12_20080303_011622_07108_M02.nc')
    #'OR1ASWC12_20080303_025745_07109_M02.nc')
    #'OR1ASWC12_20080303_043905_07110_M02.nc')
    #'OR1ASWC12_20080303_062028_07111_M02.nc')
    #'OR1ASWC12_20080303_080148_07112_M02.nc')
    #'OR1ASWC12_20080303_094311_07113_M02.nc')
#    'OR1ASWC12_20080303_112431_07114_M02.nc') #ok
#    'OR1ASWC12_20080303_130552_07115_M02.nc') # good
#    'OR1ASWC12_20080303_144715_07116_M02.nc') #ok
#    'OR1ASWC12_20080303_162835_07117_M02.nc') #good 
#    'OR1ASWC12_20080303_180958_07118_M02.nc') #ok
    'OR1ASWC12_20080303_195118_07119_M02.nc') # good - best
    #'OR1ASWC12_20080303_213241_07120_M02.nc')
    #'OR1ASWC12_20080303_231401_07121_M02.nc')
    #'OR1ASWC12_20080304_005524_07122_M02.nc')
    #'OR1ASWC12_20080304_023645_07123_M02.nc')
#    'OR1ASWC12_20080304_041807_07124_M02.nc')
#    'OR1ASWC12_20080304_055928_07125_M02.nc')
    #'OR1ASWC12_20080304_074050_07126_M02.nc')
#    'OR1ASWC12_20080304_092211_07127_M02.nc') #ok
#    'OR1ASWC12_20080304_110331_07128_M02.nc') #ok+
    #'OR1ASWC12_20080304_124454_07129_M02.nc')
    #'OR1ASWC12_20080304_142615_07130_M02.nc')
    #'OR1ASWC12_20080304_160737_07131_M02.nc')
#    'OR1ASWC12_20080304_174858_07132_M02.nc') #ok
    #'OR1ASWC12_20080304_193020_07133_M02.nc')
    #'OR1ASWC12_20080304_211141_07134_M02.nc')
    #'OR1ASWC12_20080304_225303_07135_M02.nc')

#    time= file[10:25]


SCAT= xr.open_dataset(filedir+ file)

SCAT.lon.values[SCAT.lon.values > 180] -= 360
SCAT= SCAT.where(np.logical_and(SCAT.lat > 55, SCAT.lat <83), drop=True)
SCAT= SCAT.where(np.logical_and(SCAT.lon > -20, SCAT.lon <50), drop=True)

if SCAT_type== 'ASCAT':
    SCAT= SCAT.where(np.logical_and(SCAT.wvc_index!= 41 , SCAT.wvc_index!= 42), drop=True)


Lon,Lat = map(SCAT.lon.values,SCAT.lat.values)

#
#Wind_speed= np.copy(SCAT.wind_speed.values)
#Wind_speed[np.logical_or(SCAT.wvc_index.values ==  41 , SCAT.wvc_index.values ==  42)]= np.nan
#
PlotWindVelo(Lon, Lat, SCAT.wind_speed.values, map, Umax= 25, color='YlBu')    
  
#if SCAT_type== 'ASCAT':

PlotLocalMax(SCAT.wind_speed.values, threshold=18, distance=(20,20), map= map, lon=SCAT.lon.values, lat=SCAT.lat.values, typ='max',
             color='orange', dot=False, roundorder= 0, latbound= [68, 80], lonbound= [-12, 50])    


SCAT_time= np.datetime_as_string(np.nanmin(SCAT.time), unit='m')

if save:
    savename= savedir+ SCAT_type +'_'+SCAT_time   
    plt.savefig(savename, bbox_inches= 'tight' )
    print(savename)