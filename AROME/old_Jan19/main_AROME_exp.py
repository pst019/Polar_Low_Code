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

t= -4

# Read in netcdf file
from script_imp_AROME_exp import *

#d= data(filename= '/home/'+user+'/Data/pl/vilje/fc2015121300_ctr_fp.nc')

#ncsfx = NetCDFFile('/home/'+user+'/Data/pl/fc2015121300_ctr_sfx.nc')
#SST = ncsfx.variables['SST'][:]

"""file 2"""
#nc2 = NetCDFFile('/home/patricks/Data/fc2015121300_SST_fp.nc') #at laptop
#d2= data(filename= '/home/'+user+'/Data/pl/vilje/fc2015121300_FLX_fp.nc')


# Make the map
plt.figure(4)
plt.clf()


#plt.subplot(2, 2, 1)

m = AA_map()
Lon,Lat = m(d.lon,d.lat)



# Draw contours for mslp
PlotContours(Lon, Lat, d.mslp[t, 0]/100, m)

PlotSurfTemp(Lon,Lat, d.T[t,-1] - 273.15, m)


#Umax= np.max([d.U10[t,1], d2.U10[t,1]])
#PlotWindVelo(Lon, Lat, d.U10[t,-1,:,:], m)
PlotWind(d.lon, d.lat, d.u10[t,-1], d.v10[t,-1], m)



## Make the map
#plt.subplot(2, 2, 2)
#m = make_map(ncfile=nc2)
#
#
#
#
## Draw contours for mslp
#clevs = np.arange(960,1080,4)
#cmslp = m.contour(x,y,mslp_2[t,0,:,:]/100.,clevs,colors='black',linewidths=1.)
#
## Contour labels
#plt.clabel(cmslp, fontsize=9, inline=1, fmt = '%1.0f') 
#
## Draw filled contours for T2m
##ctemp = m.contourf(x,y,T[t,-1,:,:]-273.15,cmap=plt.cm.RdBu_r)
### Colorbar for filled contour plot
##cb = m.colorbar(ctemp,"right", size="5%", pad="2%")
##cb.set_label('Temperature at 2m (C)')
#
#
#cwind = m.pcolormesh(x,y,U10_2[t,-1,:,:],cmap=plt.cm.YlGnBu, vmax= Umax)
#cb = m.colorbar(cwind,"right", size="5%", pad="2%")
#cb.set_label('Wind speed 10 m (m/s)')
#
#
##plt.subplot(2, 2, 3)
##m = make_map(ncfile=nc2)
##x,y = m(lon,lat)
##cwind = m.contourf(x,y,SST[t,-1]-273.15,cmap=plt.cm.YlGnBu)
##cb = m.colorbar(cwind,"right", size="5%", pad="2%")
##cb.set_label('SST ($^{\circ}$C)')
##
##plt.subplot(2, 2, 4)
##m = make_map(ncfile=nc2)
##x,y = m(lon,lat)
##cb = m.contourf(x,y,SST_2[t,-1]-273.15,cmap=plt.cm.YlGnBu)
##cb = m.colorbar(cwind,"right", size="5%", pad="2%")
##cb.set_label('SST (K)')
#
#
#colors= ['#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061']
#boxnr= len(colors)
#new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N= boxnr )
#
#diffmax= np.max(np.abs(U10_2[t,-1,:,:]-U10[t,-1,:,:]))
#
#plt.subplot(2, 2, 3)
#m = make_map(ncfile=nc2)
#
#cwind = m.pcolormesh(x,y,U10_2[t,-1,:,:]-U10[t,-1,:,:],cmap= new_map, vmin= -diffmax, vmax= diffmax)
#cb = m.colorbar(cwind,"right", size="5%", pad="2%")
#cb.set_label('Wind speed 10 m (m/s)')
#
#
#diffmax= np.max(np.abs((mslp_2[t,0,:,:]-mslp[t,0,:,:])/100.))
#plt.subplot(2, 2, 4)
#m = make_map(ncfile=nc2)
#
#cwind = m.pcolormesh(x,y,(mslp_2[t,0,:,:]-mslp[t,0,:,:])/100.,cmap= new_map, vmin= -diffmax, vmax= diffmax)
#cb = m.colorbar(cwind, "right", size="5%", pad="2%")
#cb.set_label('Surface pressure difference (m/s)')
#
##title
#plt.suptitle('Time: ' +time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(tim[t])))
#
##plt.show()
## Save the figure
##fileout = "t2m.png"
##savefig(fileout, bbox_inches='tight', pad_inches=0.5)
