#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  4 10:49:51 2018

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

from netCDF4 import Dataset
import numpy as np
#import time
import sys

sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
#from f_imp_AROME_exp import *  # Read in netcdf file
#from f_meteo import *
#from f_useful import *


#fignr= 10
#
day= 4
hour= 3
#
#"""a) Lambert coordinates"""
#lllon, lllat, urlon, urlat= -15, 65, 50, 75
#lat0, lon0= 75, 0 #(lat0, lon0) = center point
#

#from f_useful import *
Mediadir= '/media/'+user+'/1692A00D929FEF8B/PL/Sat_AMRC/'

file= '/home/'+user+'/home/Polar_Low/Sat_image_attempt/noaa200803040132mne_Steinar.nc'
#file= '/home/'+user+'/home/Polar_Low/Sat_image_attempt/noaa200803040132mne.nc'

#file= '/home/'+user+'/home/Polar_Low/Sat_image_attempt/noaa200803040313mne.nc'

#file= '/home/'+user+'/home/Polar_Low/Sat_image_attempt/noaa201808211436mne.nc'

#for 2008 it has about a resolution of 5km
#file=Mediadir+ 'Arctic.Composite.5km.Infrared.2008.03.'+str(day).zfill(2)+'.'+str(hour).zfill(2)+'Z.nc'
nc= Dataset(file)

data= nc.variables['Band1'][:]
#data= nc.variables['data'][:][0]


plt.figure(1)
plt.clf()
#plt.plot(data)

BT = 323.15 - 0.5*data

colors= [(plt.cm.Greys(h)) for h in range(256)]
new_map= plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors) 
cs= plt.contourf(BT, cmap= new_map)
plt.colorbar(cs)
#dlon= nc.variables['lon'][:]
#dlat= nc.variables['lat'][:]
#dT= nc.variables['data'][0,:] #temperature [K]
#
#
#
#plt.figure(fignr)
#fignr+=1
#plt.clf()
#
#map = AA_map_half()
##map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
#Lon,Lat = map(dlon,dlat)
#
#PlotColorMap4(Lon, Lat, dT, map, label='Temperature [K]', color='grey', bounds= np.arange(210, 281, 1))
#
#
##plt.contourf(dT)