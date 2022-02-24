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
    
    print(lllon, lllat, urlon, urlat)
    print(lat1, lat2, lat0, lon0)
    
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
    

def read_stations(filepath="/home/teresav/R/work/stations_FF.txt"):
    """ Reads in station listing file from the disk. Returns latitude and longitude"""
    
    # Open the file if it exists
    if os.path.isfile(filepath) and os.stat(filepath).st_size != 0:
        f = open(filepath, 'r')

        # Get the data
        stationlist = np.genfromtxt(f, delimiter=';')
        np.delete(stationlist,0,0)
        
        # Return data
        return stationlist
    
        # Close the file
        f.close()
        
def read_txt(filepath="/home/teresav/R/work/biasFF_20130301-20130328stations_PR.txt"):
    """ Reads in text file and returns data"""
    
    # Make a check if the file exists and if it is empty
    
    # Open the file if it exists
    if os.path.isfile(filepath) and os.stat(filepath).st_size != 0:
        f = open(filepath, 'r')

        # Get the data
        data = np.genfromtxt(f, delimiter=' ')
        
        # Return data
        return data
    
        # Close the file
        f.close()    
