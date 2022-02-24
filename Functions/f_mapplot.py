#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 13:25:51 2017

@author: pst019
Plot raw maps
"""

import numpy as np


import matplotlib.pyplot as plt
from scipy import stats
from mpl_toolkits.basemap import Basemap




def make_map(ncfile, ax=None, fill=False,proj='lcc',res='l',tick_incr=[5,1]):
    """ Reads in a netcdf file and returns a map
    used in AROME/project_NWP"""
    
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

    print(lllon, lllat, urlon, urlat, earth_radius, lat1, lat2, lat0, lon0)
    # Make the map
    map = Basemap(llcrnrlon=lllon, llcrnrlat=lllat, urcrnrlon=urlon, urcrnrlat=urlat, rsphere=earth_radius,\
            resolution=res,projection=proj,\
            lat_1=lat1,lat_2=lat2,lat_0=lat0,lon_0=lon0)
    
    
    ticklon = np.array(tick_incr)[0]
    try:
        ticklat = np.array(tick_incr)[1]
    except IndexError:
        ticklat = ticklon

    # Draw the latitudes and the longitudes
    parallels = np.arange(0.,90,5.)
    map.drawparallels(parallels,labels=[True,False,False,False])    
    meridians = np.arange(10.,361.,10.)
    map.drawmeridians(meridians,labels=[False,False,False,True])

    # Draw the coastline
    map.drawcoastlines(color='0.5')

    if fill:
        map.drawlsmask(land_color='0.8', ocean_color='w')
    return map


    
def WRF_map(lat, lon, lat0, lon0, res='l', fill=False):
    """lat and lon should be 2D matrixes from which the corner points are taken
    lat0- some latitude standard parallel
    lon0- the parallel longitude"""
    lllon = lon[0][0]
    urlon = lon[-1][-1]
    lllat = lat[0][0]
    urlat = lat[-1][-1]   
    map = Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0, res, fill)
    return map


    
def ASR_map(res='c', area_thresh= None):
    """map of the ARCTIc SYSTEM REANALYSIS - just upside down (opposite with lon_0= 185)
    res= 'l' increases the resolution with area_thresh it can be further reduced (c=10000)"""
    map = Basemap(projection='npstere',boundinglat=41,lon_0=5,resolution=res, area_thresh= area_thresh)
    map.drawcoastlines()
    
    map.drawparallels(np.arange(-80.,81.,10.))
    map.drawmeridians(np.arange(-180.,181.,10.))
    return map    

    
    
def Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0, res='l', fill=False, coastline=True, latdist=5, londist=10):
    """lllon - lon lower left corner ...
    lat0 - latitude standard parallel, should be somewhere in the center of the domain
    lon0 - the parallel longitude
    http://matplotlib.org/basemap/api/basemap_api.html"""
    rsphere=(6378137.00,6356752.3142)
    map = Basemap(llcrnrlon=lllon, llcrnrlat=lllat, urcrnrlon=urlon, urcrnrlat=urlat, rsphere=rsphere,
            resolution=res,projection='lcc', lat_0=lat0,lon_0=lon0)
    
    # Draw the latitudes and the longitudes
    parallels = np.arange(0.,90,latdist)
    map.drawparallels(parallels,labels=[True,False,False,False])    
    meridians = np.arange(0,361., londist)
    map.drawmeridians(meridians,labels=[False,False,False,True])

    # Draw the coastline
    if coastline==True:
        map.drawcoastlines(color='0.5')

    if fill:
        map.drawlsmask(land_color='0.8', ocean_color='w')
    return map    


def AA_map(res='l', fill=False):
    """make a map of Arome Arctic"""
    lllon, lllat, urlon, urlat= -17.956999588, 69.2990002573, 68.8371159334, 71.153315265
    lat0, lon0= 77.5, -25.0 #(lat0, lon0) = center point
    map = Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0, res, fill)
    return map
 
def AA_South_map(res='l', fill=False):
    """make a map of Arome Arctic"""
    lllon, lllat, urlon, urlat= -19.5, 65.5, 52, 71
    lat0, lon0= 77.5, -25.0 #(lat0, lon0) = center point
    map = Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0, res, fill)
    return map    

def AA_map_half(res='l', fill=False):
    """make a map of Arome Arctic"""
    lllon, lllat, urlon, urlat= -17.6, 69.45, 45, 72
    lat0, lon0= 77.5, -25.0 #(lat0, lon0) = center point
    map = Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0, res, fill)
    return map    

def Lat_Lon_map(lllon, lllat, urlon, urlat, res='l', fill=False, proj= 'merc'):
    """quadratic map with same distance between each lon and each lat
    proj= merc - mercator (all is bad near the pole)
        cyl- equidistant cylindrical - easiest projection
        cea - cyclindrical equal area"""
    map = Basemap(llcrnrlon=lllon, llcrnrlat=lllat, urcrnrlon=urlon, urcrnrlat=urlat,
                resolution=res, lat_ts= 75 #(lllat+urlat)/2
                , projection=proj)
    
    # Draw the latitudes and the longitudes
    parallels = np.arange(0.,90,5.)
    map.drawparallels(parallels,labels=[True,False,False,False])    
    meridians = np.arange(10.,361.,10.)
    map.drawmeridians(meridians,labels=[False,False,False,True])

    # Draw the coastline
    map.drawcoastlines(color='0.5')

    if fill:
        map.drawlsmask(land_color='0.8', ocean_color='w')
    return map       


    
def Polar_map(latsouth, hemisp='NH', fillcontinents=False):
    """make a polar centered map"""
    if hemisp == 'NH':
        map= Basemap(projection= 'nplaea',resolution='c',lon_0=0,boundinglat= latsouth ,area_thresh=10000, round=1)
    elif hemisp == 'SH':
        map= Basemap(projection= 'splaea',resolution='c',lon_0=180,boundinglat= latsouth ,area_thresh=10000, round=1)

    if fillcontinents==True: map.fillcontinents(color='#ddaa66') #'0.9')#'gray')
    if fillcontinents==False: map.drawcoastlines(color='black')
    map.drawmeridians(np.arange(0,360,10))#, labels=[0,0,0,1])
    map.drawparallels(np.arange(-90,90,10))#, labels=[1,0,0,0])
#    plt.tight_layout()
    return map
    
    

