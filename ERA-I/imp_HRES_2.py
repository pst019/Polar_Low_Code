#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

This is for importing HRES data. It is also a test to see how gribdata actually looks like
"""

from netCDF4 import Dataset
import numpy as np
import time
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
import matplotlib.pyplot as plt

import os
user = os.getcwd().split('/')[2]

from f_useful import *

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'

#year
#month
#day
#tstart

filename= Mediadir+'ECMWF/HRES/oper_sfc20180216_0000_Sval_nice.nc'

t= 12 #lacktime


nc= Dataset(filename)

#print(nc.variables.keys())
#this one is even more interesting:
#print(nc.variables.values())

tim = nc.variables['time'][:]
lat= nc.variables['latitude'][:]
lon= nc.variables['longitude'][:]

T2= nc.variables['t2m'][:]
z= nc.variables['z'][:]/9.81 #geopotential of the surface translated to surface height

plt.figure(1)
plt.clf()                  
                    
                    
#Svalbard map
#map= Lambert_map(lllon=9, lllat=75, urlon=25, urlat=81, lat0= 70, lon0= 20, res='h', fill=False)
#maparea == "Advent":
map= Lambert_map(lllon=15, lllat=78.05, urlon=17, urlat=78.35, lat0= 75, lon0= 16, res='h', fill=False)
# maparea == 'Advent_med':
#map= Lambert_map(lllon=15.35, lllat=78.09, urlon=16.55, urlat=78.31, lat0= 75, lon0= 16, res='h', fill=False)

grid= np.meshgrid(lon, lat) 
Lon,Lat = map(grid[0],grid[1])

PlotContours(Lon, Lat, z[t], map, leveldist=None,levels=[10, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000], numbers=True, color= 'k')

PlotColorMap4(Lon, Lat, T2[t], map, bounds=np.arange(252, 268, 0.5), label='Temperature [$^{\circ}$C]')