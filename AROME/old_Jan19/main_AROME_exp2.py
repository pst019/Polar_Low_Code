#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

This is from plot_comparison
"""


import time, calendar, datetime, numpy
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import netCDF4 as Dataset
from netCDF4 import Dataset as NetCDFFile

import numpy as np
import os
from pylab import *
from datetime import date, datetime, timedelta
from scipy import stats

# import own modules
from f_mapplot import *
from f_plot_fields import *

t= 10

# Read in netcdf file
from f_imp_AROME_exp import *

#d= data(filename= '/home/'+user+'/Data/pl/vilje/fc2015121300_ctr_fp.nc')

d= data(filename= Mediadir+'PL/ec/080303_test_3_2008030300_fp_extract.nc')
#d.imp_standard()

#ncsfx = NetCDFFile('/home/'+user+'/Data/pl/fc2015121300_ctr_sfx.nc')
#SST = ncsfx.variables['SST'][:]

"""file 2"""
#nc2 = NetCDFFile('/home/patricks/Data/fc2015121300_SST_fp.nc') #at laptop
#d2= data(filename= '/home/'+user+'/Data/pl/vilje/fc2015121300_FLX_fp.nc')


# Make the map
plt.figure(2)
plt.clf()

plt.title(str(time.localtime(d.tim[t])[:4]))
#plt.subplot(2, 2, 1)

m = AA_map()
Lon,Lat = m(d.lon,d.lat)



d.imp_level(pn= 0, tn= t)


# Draw contours for mslp
PlotContours(Lon, Lat, d.mslp, m)
#
PlotSurfTemp(Lon,Lat, d.Temp - 273.15, m)
#
#
##Umax= np.max([d.U10[t,1], d2.U10[t,1]])
##PlotWindVelo(Lon, Lat, d.U[t,-1,:,:], m)
PlotWind(d.lon, d.lat, d.u, d.v, m)


d.imp_level(pn= 10, tn= t)
#PlotContours(Lon, Lat, d.mslp, m)

PlotPV(Lon, Lat, d.PV, m)


plt.figure(3)
plt.clf()

#d.imp_cross_sec(20, tn= t)
##PlotCross_sec()
#nrlevels= 10 #np.arange(960,1080,4)
#cs= contour(d.lon[20], d.press , d.T, nrlevels, linewidths= 1. , colors= 'k' )
#plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f', color= 'black')
#plt.gca().invert_yaxis()
#plt.ylabel('Pressure [hPa]')
