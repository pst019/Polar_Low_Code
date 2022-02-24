#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

This is from plot_comparison
"""

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import numpy as np
import os
import time
#from datetime import date, datetime, timedelta
from scipy import stats

# import own modules
from f_imp_WRF import *  # Read in netcdf 

import sys  #to import the functions from a different directory
sys.path.insert(0, '/home/'+user+'/codeAROME/')
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)


domain= 1
year, month, day, hour= 2008, 3, 3, 0
t= 8 #which timestep of the data is regarded
l= 0 #level

#filename= (Mediadir+'PL/WRF/wrfout_d0'+str(domain)+'_'+str(year)
#+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+str(hour).zfill(2)+':00:00_INTRP' )
#
#nc= Dataset(filename)
#print(nc.variables.keys())
#u10= nc.variables['U10'][:]
#v10= nc.variables['V10'][:]
#geo= nc.variables['GHT'][:, l, :]
#T= nc.variables['TT'][:, l]


d= WRFdata(filename= Mediadir+'PL/WRF/wrfout_d0'+str(domain)+'_'+str(year)
+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+str(hour).zfill(2)+':00:00_INTRP')
d.imp_standard()
#
#
plt.figure(2)
plt.clf()
map= WRF_map(d.lat, d.lon, d.lat0, d.lon0)
Lon,Lat = map(d.lon,d.lat)


"""MSLP, Surface temperature and surface winds"""
PlotSurfTemp(Lon,Lat, d.T[t, l]-273, map)
PlotWind(d.lon, d.lat, d.u10[t], d.v10[t], map, rot= False)

PlotContours(Lon, Lat, mslp(d.geo[t,l],d.T[t,l]), map)


"""domain 2"""
domain= 2

d= WRFdata(filename= Mediadir+'PL/WRF/wrfout_d0'+str(domain)+'_'+str(year)
+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+str(hour).zfill(2)+':00:00_INTRP')
d.imp_standard()
#
#
plt.figure(3)
plt.clf()
map= WRF_map(d.lat, d.lon, d.lat0, d.lon0)
Lon,Lat = map(d.lon,d.lat)


"""MSLP, Surface temperature and surface winds"""
PlotSurfTemp(Lon,Lat, d.T[t, l]-273, map)
PlotWind(d.lon, d.lat, d.u10[t], d.v10[t], map, rot= False)

PlotContours(Lon, Lat, mslp(d.geo[t,l],d.T[t,l]), map)


"""domain 3"""
domain= 3

d= data(filename= Mediadir+'PL/WRF/wrfout_d0'+str(domain)+'_'+str(year)
+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+str(hour).zfill(2)+':00:00_INTRP')
d.imp_standard()
#
#
plt.figure(4)
plt.clf()
map= WRF_map(d.lat, d.lon, d.lat0, d.lon0)
Lon,Lat = map(d.lon,d.lat)


"""MSLP, Surface temperature and surface winds"""
PlotSurfTemp(Lon,Lat, d.T[t, l]-273, map)
PlotWind(d.lon, d.lat, d.u10[t], d.v10[t], map, rot= False)

PlotContours(Lon, Lat, mslp(d.geo[t,l],d.T[t,l]), map)

#