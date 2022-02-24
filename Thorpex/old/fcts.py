#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 12:57:58 2016

@author: pst019

Functions

"""

import numpy as np
import os
from pylab import *
from datetime import date, datetime, timedelta
import matplotlib.pyplot as plt
from scipy import stats
from mpl_toolkits.basemap import Basemap


def make_map(ncfile, ax=None, fill=False,proj='lcc',res='l',tick_incr=[5,1],domain="MetCoop"):
    """ Reads in a netcdf file and returns a map """
    
    lat = ncfile.variables['latitude'][:]
    lon = ncfile.variables['longitude'][:]
    lllon = lon[0][0]
    urlon = lon[-1][-1]
    lllat = lat[0][0]
    urlat = lat[-1][-1]
    
    # Set the projection attributes 
    lon0 = ncfile.variables['projection_lambert'].getncattr('longitude_of_central_meridian')
    lat1 = ncfile.variables['projection_lambert'].getncattr('standard_parallel')[0]
    lat2 = ncfile.variables['projection_lambert'].getncattr('standard_parallel')[1]
    lat0 = ncfile.variables['projection_lambert'].getncattr('latitude_of_projection_origin')
    earth_radius = ncfile.variables['projection_lambert'].getncattr('earth_radius')
    
    if ax is None:
        ax = plt.gca()    
    
    # Make the map
    m = Basemap(llcrnrlon=lllon, llcrnrlat=lllat, urcrnrlon=urlon, urcrnrlat=urlat, rsphere=earth_radius,\
            resolution=res,projection=proj,\
            lat_1=lat1,lat_2=lat2,lat_0=lat0,lon_0=lon0)
    
    ticklon = np.array(tick_incr)[0]
    try:
        ticklat = np.array(tick_incr)[1]
    except IndexError:
        ticklat = ticklon

    # Draw the latitudes and the longitudes
    parallels = np.arange(0.,90,5.)
    m.drawparallels(parallels,labels=[True,False,False,False])    
    meridians = np.arange(10.,361.,10.)
    m.drawmeridians(meridians,labels=[False,False,False,True])

    # Draw the coastline
    m.drawcoastlines(color='0.5')

    if fill:
        m.drawlsmask(land_color='0.8', ocean_color='w')
    return m

class Arome:
    def __init__(self, nc):
        self.tim = nc.variables['time'][:]
        self.lat = nc.variables['latitude'][:]
        self.lon = nc.variables['longitude'][:]
        self.height= nc.variables['height_above_msl'][:]
        
        # Variables
        #T = nc.variables['air_temperature_pl'][:]
        self.mslp = nc.variables['air_pressure_at_sea_level'][:]/100
        self.u = nc.variables['x_wind_pl'][:]
        self.v = nc.variables['y_wind_pl'][:]
        self.U = sqrt(self.u**2+self.v**2)
    
class thorpex:    
    def __init__(self, nc):
        self.tim = nc.variables['time'][:]
        self.lat = nc.variables['lat'][:]
        self.lon = nc.variables['lon'][:]
        self.alt= nc.variables['alt'][:]
        
        # Variables
        self.T = nc.variables['tdry'][:]
        self.pres = nc.variables['pres'][:]
        self.u = nc.variables['u_wind'][:]
        self.v = nc.variables['v_wind'][:]
        #U = sqrt(u10**2+v10**2)
        self.U = nc.variables['wspd'][:] #wind speed
        self.Udir = nc.variables['wdir'][:]#wind direction
        
        #also potential temperatures available