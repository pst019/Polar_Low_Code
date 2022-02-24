
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
from mapplot_def import *


t= 2

## Read in netcdf file
##nc = NetCDFFile('/home/patricks/Data/fc2015121300_fp.nc') #at laptop
#nc = NetCDFFile('/home/'+user+'/Data/pl/fc2015121300_ctr_fp.nc')
#
##print(nc.variables.keys()) #print(\"Variables in data set: {}\".format(\", \".join(nc.variables.keys())))
#
#tim = nc.variables['time'][:]
#lat = nc.variables['latitude'][:]
#lon = nc.variables['longitude'][:]
##height= nc.variables['height_above_msl'][:]
#
## Variables
##T = nc.variables['air_temperature_pl'][:]
#mslp = nc.variables['air_pressure_at_sea_level'][:] 
#u10 = nc.variables['x_wind_pl'][:]
#v10 = nc.variables['y_wind_pl'][:]
#U10 = sqrt(u10**2+v10**2)
#
#ncsfx = NetCDFFile('/home/'+user+'/Data/pl/fc2015121300_ctr_sfx.nc')
#SST = ncsfx.variables['SST'][:]
#
#"""file 2"""
#nc2 = NetCDFFile('/home/'+user+'/Data/pl/fc2015121300_SST_fp.nc')
#
## Variables
##T2 = nc2.variables['air_temperature_pl'][:]
#mslp_2 = nc2.variables['air_pressure_at_sea_level'][:] 
#u10_2 = nc2.variables['x_wind_pl'][:]   
#v10_2 = nc2.variables['y_wind_pl'][:]
#U10_2 = sqrt(u10_2**2+v10_2**2)
#
#ncsfx_2 = NetCDFFile('/home/'+user+'/Data/pl/fc2015121300_SST_sfx.nc')
#SST_2 = ncsfx_2.variables['SST'][:]


# Make the map
plt.figure(2)
plt.clf()


for t in range(2):
    plt.subplot(2, 2, t+1)
    m = make_map(ncfile=nc)
    x,y = m(lon,lat)
    
    #nc.close()
    
    
    # Draw contours for mslp
    clevs = np.arange(960,1080,4)
    cmslp = m.contour(x,y,mslp[t,0,:,:]/100.,clevs,colors='black',linewidths=1.)
    
    # Contour labels
    plt.clabel(cmslp, fontsize=9, inline=1, fmt = '%1.0f') 
    
    # Draw filled contours for T2m
    #ctemp = m.contourf(x,y,T[t,-1,:,:]-273.15,cmap=plt.cm.RdBu_r)
    ## Colorbar for filled contour plot
    #cb = m.colorbar(ctemp,"right", size="5%", pad="2%")
    #cb.set_label('Temperature at 2m (C)')
    
    
    cwind = m.contourf(x,y,U10[t,-1,:,:],cmap=plt.cm.YlGnBu)
    cb = m.colorbar(cwind,"right", size="5%", pad="2%")
    cb.set_label('Wind speed 10 m (m/s)')

    # Title 
    plt.title('Time: ' +time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(tim[t])))


#plt.show()
# Save the figure
#fileout = "t2m.png"
#savefig(fileout, bbox_inches='tight', pad_inches=0.5)
