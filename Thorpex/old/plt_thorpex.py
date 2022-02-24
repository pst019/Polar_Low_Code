#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 12:57:58 2016

@author: pst019

Plot the thorpex data

"""

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import os

# import own modules
#from fcts import *


"""first run import_thorpex.py"""
os.system('python import_thorpex.py')


"""mape the map"""
plt.figure(1)
plt.clf()
m = Basemap(resolution='c',projection='aea',
            llcrnrlon= min(lon), llcrnrlat=min(lat)-.5,
            urcrnrlon= max(lon)+3, urcrnrlat=max(lat), 
            lon_0=(max(lon)+min(lon))/2, lat_0=(max(lat)+min(lat))/2)

parallels = np.arange(0.,90,1.)
m.drawparallels(parallels,labels=[True,False,False,False])    
meridians = np.arange(-50,361,2)
m.drawmeridians(meridians,labels=[False,False,False,True])

m.drawcoastlines(color='0.5')


"""plot the observation points and their time"""
xpt, ypt= m(lon, lat)
m.plot(xpt,ypt,'bo')

for i in range(len(UTC)):
    plt.text(xpt[i]+100, ypt[i]+100, UTC[i])


    
#from scipy.interpolate import interp2d

"""make lat lon grid"""
resol= 50
loni= np.linspace(min(lon), max(lon), resol)
lati= np.linspace(min(lat), max(lat), resol)
long, latg= np.meshgrid(loni, lati)
lonm, latm= m(long, latg) #translate it to basemap

"""interpolate the data"""
data=pres
#data= U

from matplotlib.mlab import griddata
presi= griddata(lon, lat, data, long, latg)#, interp='linear')


"""make the lines"""
#clevs = np.arange(960,1080,1)
cmslp = m.contour(lonm,latm, presi, colors='black',linewidths=1.)
plt.clabel(cmslp, fontsize=9, inline=1, fmt = '%1.0f') 

"""make the color plot"""
m.pcolormesh(lonm,latm, presi)