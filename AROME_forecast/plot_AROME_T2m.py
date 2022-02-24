
# Import modules
import time, calendar, datetime, numpy
from mpl_toolkits.basemap import Basemap
import matplotlib
import matplotlib.cm as cm
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import os
from pylab import *

#
# Map definition
def make_map(ncfile, ax=None, fill=False,proj='lcc',res='l',tick_incr=[5,1]):
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
            resolution='l',projection=proj,\
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


# Directory to save the figures 
figdir = os.path.join('./')

# Read in netcdf file
url = 'http://thredds.met.no/thredds/dodsC/aromearcticarchive/2018/04/06/arome_arctic_full_2_5km_20180406T21Z.nc'
nc = NetCDFFile(url)       

# Variables 
#
# You can check the content for example by 
# ndcump -h http://thredds.met.no/thredds/dodsC/aromearcticarchive/2018/04/06/arome_arctic_full_2_5km_20180406T21Z.nc

lat = nc.variables['latitude'][:]
lon = nc.variables['longitude'][:]
T2m = nc.variables['air_temperature_2m'][:] 
nctime = nc.variables['time'][:]

# Loop through the times
ntimes = 3
for i in range(0,ntimes,1):

    # Make the figure
    plt.figure(figsize=(20,10))
        
    # Make the map
    m = make_map(ncfile=nc)
    x,y = m(lon,lat)  
             
    # Draw the variable
    cvar = m.contourf(x,y,T2m[i,0,:,:],cmap=plt.cm.RdBu_r)

    # Add colorbar and label
    cb = m.colorbar(cvar,"right", size="5%", pad="2%")
    cb.set_label(nc.variables['air_temperature_2m'].getncattr('long_name'))

    # Make title
    plt.title(nc.title +" " +time.strftime('%Y-%m-%d %H:%M', time.gmtime(nctime[i])))

#    plt.show()

    # Save the figure
    figname = os.path.join("AA_T2m_" + time.strftime('%Y%m%d%H%M', time.gmtime(nctime[i])) + ".png")
    figout = os.path.join(figdir, figname)
    
    plf.savefig(figout, bbox_inches='tight', pad_inches=0.5)
    
    plt.clf()

