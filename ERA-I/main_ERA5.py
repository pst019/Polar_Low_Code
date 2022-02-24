#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 13:25:51 2017

@author: pst019
from ASR-ERA_Inv_random_3 of folder ASR
"""

import pickle
import sys
import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})


sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')
from f_imp_ERA2 import *
from f_imp_ASR import *
from f_impLists import *
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_meteo import *



"""global variables"""

fignr= 3

#maptype='AA'
maptype='AA_half'
#maptype='Lambert'
#maptype='polar'

if maptype=='AA': plt.figure(fignr, figsize= (6, 6))
elif maptype=='AA_half': plt.figure(fignr, figsize= (6, 4.7))
else: plt.figure(fignr, figsize= (6, 4.5))
fignr += 1
plt.clf()

"""Lambert coordinates"""
lllon, lllat, urlon, urlat= -15, 63, 60, 75
lat0, lon0= 75, 0 #(lat0, lon0) = center point

                   
if maptype== 'AA': map = AA_map()
elif maptype== 'AA_half': map = AA_map_half()
elif maptype=='polar': map= Polar_map(latsouth= 50) #d.lat[-1])
else: map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)

title_extra='jet'

"""time"""
year=2008
month=3
day= 3
hour= 18

t= hour
plevel= 850

#var= 'PressWind'
#var= 'PressWind_advanced'
#var= 'RH'
var='SH'
#var='OLW'



#"""c) plot characteristics"""
#pleveldist= 4
#
#"""end global variables """
#
#title_extra=''
#"""module starts"""
#   
#
save= True
savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/Fields/'
# 

varimp= 'T'
from netCDF4 import Dataset, num2date
file= Mediadir + 'ECMWF/ERA5/'+str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+varimp+'.nc'

nc= Dataset(file)
#print(nc.variables.keys())
                    
tim= nc.variables['time'][t]
datetime= num2date(tim, nc.variables['time'].units)

level= nc.variables['level'][:]

#if varimp== 'T': 
T= nc.variables['t'][t]
print('T shape ', T.shape)

#if varimp== 'SH': 
varimp= 'SH'    
file= Mediadir + 'ECMWF/ERA5/'+str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+varimp+'.nc'  
nc= Dataset(file)

tim= nc.variables['time'][t]
datetime= num2date(tim, nc.variables['time'].units)                             
SH= nc.variables['q'][t]
lat= nc.variables['latitude'][:] #it can be different lat-lon
lon= nc.variables['longitude'][:]

print('SH shape', SH.shape, 'maybe download variables new')

varimp= 'U'    
file= Mediadir + 'ECMWF/ERA5/'+str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+varimp+'.nc'  
nc= Dataset(file) 
u= nc.variables['u'][t]
v= nc.variables['v'][t]
U= np.sqrt(u**2+ v**2)

varimp= 'SLP'    
file= Mediadir + 'ECMWF/ERA5/'+str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+varimp+'.nc'  
nc= Dataset(file) 
SLP= nc.variables['msl'][t]/100

ilevel = np.where(level== plevel)[0][0]


varimp= 'OLW'
if day== 4 and hour <6: E5day= 3
else: E5day= day
file= Mediadir + 'ECMWF/ERA5/'+str(year)+'-'+str(month).zfill(2)+'-'+str(E5day).zfill(2)+'_'+varimp+'_01.nc'  
nc= Dataset(file)

datetimeOLW= num2date(nc.variables['time'][:], nc.variables['time'].units)
latOLW= nc.variables['lat'][:] #it can be different lat-lon
lonOLW= nc.variables['lon'][:]
print('Specify t such that it fits') #the data starts at 6 and 18 and has 18 hourly time steps for each - so it is 6, 7, 8,.. 24, 18, 19, ... 36
t = np.argmin(np.abs(datetimeOLW - datetime))
print(t, datetimeOLW[t])
OLW= -1* nc.variables['var179'][t]/3600 #ws/m**2 to w/m**2
#   
#    
#print(year, month, day, hour, 't: ', t)

grid= np.meshgrid(lon, lat) #pcolormesh needs corner points for projection
Lon, Lat= map(grid[0], grid[1])


if var == 'PressWind_advanced':
    PlotContours(Lon, Lat, SLP, map, leveldist= 1, numbers=False)
    
    PlotWindVelo(Lon, Lat, U[ilevel], map, Umax= 25, color='YlBu')
    PlotLocalMax(SLP, threshold=1010, distance=30, map= map, lon=lon, lat=lat,
                     data2=U[ilevel], threshold2=10, distance2= 7) 

    PlotLocalMax(U[ilevel], threshold=18, distance=20, map= map, lon=lon, lat=lat, typ='max',
                         color='orange', dot=False, roundorder= 0)
    
    
if var== 'RH':
    PlotContours(Lon, Lat, SLP, map, leveldist= 1, numbers=False)

    import metpy.calc as mpcalc
    from metpy.units import units    
    dRH= mpcalc.relative_humidity_from_specific_humidity(SH[ilevel], T[ilevel]* units.K, level[ilevel]* units.hPa)
    
    PlotColorMap4(Lon, Lat, dRH , map, label= 'Relative humidity at '+str(level[ilevel])+' hPa', color='blue', bounds= np.array([0, 0.4, 0.7, 0.8, 0.9, 0.95, 1]) )
    title_extra+= 'lev_'+str(level[ilevel])

if var== 'SH':
    PlotContours(Lon, Lat, SLP, map, leveldist= 1, numbers=False)
   
    PlotColorMap4(Lon, Lat, SH[ilevel]*1E3 , map, label= 'Specific humidity [g/kg] at '+str(level[ilevel])+' hPa',
                  color='jet_r', bounds= np.arange(0, 2.3, .2))


    import metpy.calc as mpcalc
    from metpy.units import units    
    dRH= mpcalc.relative_humidity_from_specific_humidity(SH[ilevel], T[ilevel]* units.K, level[ilevel]* units.hPa)

    PlotContours(Lon, Lat, dRH, map, levels= [0.9], numbers= False, color='white')

    
    title_extra+= 'lev_'+str(level[ilevel])



if var =='OLW':
    levels= np.arange(940, 1040, 3)
    PlotContours(Lon, Lat, SLP, map, levels= levels, numbers=False, color= 'r', alpha=.6)
    grid= np.meshgrid(lonOLW, latOLW) #pcolormesh needs corner points for projection
    Lon, Lat= map(grid[0], grid[1])
    
    PlotColorMap4(Lon, Lat, OLW , map, label= 'Outgoing longwave radiation [W/m$^2$]', color='grey', bounds= np.arange(130, 241, 10))


#
#"""ERA wind and surface pressure"""
#if var == 'PressWind':   
#    d.imp_u10()
#    U10= np.sqrt(d.u10[0]**2+d.v10[0]**2)
#
#    PlotWindVelo(Lon, Lat, U10, map, Umax= 25)
#   
#    PlotWind(d.lon, d.lat, d.u10[0], d.v10[0], map)
#    PlotContours(Lon, Lat, d.MSLP[0], map, leveldist= pleveldist)
#
#
#if var == 'PressWind_advanced':
#    PlotContours(Lon, Lat, d.MSLP[0], map, leveldist= 1, numbers=False)
#
#    d.imp_u10()
#    U10= np.sqrt(d.u10[0]**2+d.v10[0]**2)
#    PlotWindVelo(Lon, Lat, U10, map, Umax= 25, color='YlBu')
#
#    PlotLocalMax(d.MSLP[0], threshold=1010, distance=10, map= map, lon=d.lon, lat=d.lat,
#                     data2=U10, threshold2=10, distance2= 7) 
#
#    PlotLocalMax(U10, threshold=15, distance=10, map= map, lon=d.lon, lat=d.lat, typ='max',
#                 color='orange', dot=False, roundorder= 0)        
#
#
#
#
#if var == 'RH850':
#    PlotContours(Lon, Lat, d.MSLP[0], map, leveldist= pleveldist)
#    d.impvar('sHum850')
#    d.impvar('T850')
#    import metpy.calc as mpcalc
#    from metpy.units import units    
#    dRH= mpcalc.relative_humidity_from_specific_humidity(d.sHum850[0], d.T850[0]* units.K, 850* units.hPa)
#    
#    PlotColorMap4(Lon, Lat, dRH , map, 'Relative humidity at 850 hPa', color='blue', bounds= np.array([0, 0.4, 0.7, 0.8, 0.9, 0.95, 1]) )
#
#
#if save== False:
#    plt.title('ERA-I '+str(d.datetime[0])[:-3]+ title_extra)
#
#plt.tight_layout()
#
if save== True:
    print('save')
    savename= savedir + 'ERA5_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+str(hour).zfill(2)+'_'+var+title_extra
    plt.savefig(savename, pad_inches=0)
    print(savename)