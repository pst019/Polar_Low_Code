#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 12:33:32 2018

@author: patricks
"""

import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'

import sys
sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')

from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
from netCDF4 import Dataset, num2date



year, month = 2008, 3
day, hour= 4, 0

#exp_name, fileday, filehour='DA_080303_CTR', 3, 0
exp_name, fileday, filehour='DA_08030212', 2, 12
#exp_name, fileday, filehour='DA_080301_cycling', 3, 0



fileday_2, filehour_2 = fileday, filehour
#exp_name_2 ='DA_080303_m2SST'
exp_name_2 ='DA_080303_2FLX'
exp_name_2, fileday_2, filehour_2='DA_080301_cycling', 4, 0


t= (day- fileday)*24 + (hour- filehour)  # -1 
t_2= (day- fileday_2)*24 + (hour- filehour_2)  # -1 


#suffix='sfx'
#variable= 'H'

suffix='fp'
#variable='y_wind_10m'
#variable='air_temperature_2m'
variable='air_temperature_0m'


fignr= 1

#maptype='AA'
maptype='AA_half'
#maptype='Lambert'

if maptype=='AA': plt.figure(fignr, figsize= (6, 6))
elif maptype=='AA_half': plt.figure(fignr, figsize= (6, 4.7))
else: plt.figure(fignr, figsize= (6, 4.5))
plt.clf()

if maptype== 'AA': map = AA_map()
elif maptype== 'AA_half': map = AA_map_half()
else: map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)    



AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_'+suffix+'_extract.nc'
nc= Dataset(AAfilename)

var= nc.variables[variable][t, 0, :]
lon=nc.variables['longitude'][:]
lat=nc.variables['latitude'][:]
Lon,Lat = map(lon,lat)

AAfilename_2= Mediadir+'PL/AA/ec/'+exp_name_2+'_'+str(year)+str(month).zfill(2)+str(fileday_2).zfill(2)+str(filehour_2).zfill(2)+'_'+suffix+'_extract.nc'
nc_2= Dataset(AAfilename_2)
var_2= nc_2.variables[variable][t_2, 0, :]


if np.max(var - var_2)== 0:
    print('seems like both datasets are similar')
else:
    PlotColorMap4(Lon, Lat, var_2- var, map, maxlevel=1, label= "difference 2- 1", color= 'RdBuWhite', symetric= True)


