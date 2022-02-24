
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


t= 1

# Read in netcdf file
#nc = NetCDFFile('/home/patricks/Data/fc2015121300_fp.nc') #at laptop
nc = NetCDFFile('/home/'+user+'/Data/pl/vilje/fc2015121300_ctr_fp.nc')

#print(nc.variables.keys()) #print(\"Variables in data set: {}\".format(\", \".join(nc.variables.keys())))

tim = nc.variables['time'][:]
lat = nc.variables['latitude'][:]
lon = nc.variables['longitude'][:]
height= nc.variables['height_above_msl'][:]

# Variables
#T = nc.variables['air_temperature_pl'][:]
mslp = nc.variables['air_pressure_at_sea_level'][:]/100
u10 = nc.variables['x_wind_pl'][:]
v10 = nc.variables['y_wind_pl'][:]
U10 = sqrt(u10**2+v10**2)

#ncsfx = NetCDFFile('/home/'+user+'/Data/pl/fc2015121300_ctr_sfx.nc')
#SST = ncsfx.variables['SST'][:]

"""file 2"""
#nc2 = NetCDFFile('/home/patricks/Data/fc2015121300_SST_fp.nc') #at laptop
#nc2 = NetCDFFile('/home/'+user+'/Data/pl/fc2015121300_FLX_fp.nc')
nc2 = NetCDFFile('/home/'+user+'/Data/pl/vilje/arome_arctic_extracted_2_5km_20151213T00Z.nc')

# Variablesnc2 = NetCDFFile('/home/'+user+'/Data/pl/fc2015121300_FLX_fp.nc')

#T2 = nc2.variables['air_temperature_pl'][:]

tim2 = nc2.variables['time'][:]
mslp_2 = nc2.variables['air_pressure_at_sea_level'][:]/100 
u10_2 = nc2.variables['x_wind_pl'][:]   
v10_2 = nc2.variables['y_wind_pl'][:]
U10_2 = sqrt(u10_2**2+v10_2**2)

nc2.close()
ncsfx_2 = NetCDFFile('/home/'+user+'/Data/pl/vilje/fc2015121300_SST_sfx.nc')
SST_2 = ncsfx_2.variables['SST'][:]


t_2= np.where(tim2 == tim[t])[0][0]

# Make the map
plt.figure(1)
plt.clf()


"""first plot"""
plt.subplot(2, 2, 1)
m = make_map(ncfile=nc)
x,y = m(lon,lat)



Umax= np.max([U10[t,-1], U10_2[t_2,-1]])

# Draw contours for mslp
plt.title('MPSL controlrun')
clevs = np.arange(960,1080,2)
cmslp = m.contour(x,y,mslp[t,0, :,:],clevs,colors='black',linewidths=1.)
plt.clabel(cmslp, fontsize=9, inline=1, fmt = '%1.0f') 

cwind = m.pcolormesh(x,y,U10[t,-1,:,:],cmap=plt.cm.YlGnBu, vmax= Umax)
cb = m.colorbar(cwind,"right", size="5%", pad="2%")
cb.set_label('Wind speed 10 m (m/s)')

neighborhood_size = 100
neighborhood_size2 = 50
pthreshold = 1010
Uthreshold = 15

data_max = filters.minimum_filter(mslp[t,0,:,:], neighborhood_size)
data_max_wind = filters.maximum_filter(U10[t,-1,:,:], neighborhood_size2)
maxima = (mslp[t,0, :,:] == data_max)

maxima[data_max > pthreshold]=0

maxima[data_max_wind < Uthreshold]= 0
xp, yp= m(lon[maxima], lat[maxima])
m.plot(xp, yp, 'bo')


print(mslp[t,0,:,:][maxima])

"""second plot"""
plt.subplot(2, 2, 2)
m = make_map(ncfile=nc)

plt.title('MPSL operational')
#clevs = np.arange(960,1080,2)
cmslp = m.contour(x,y,mslp_2[t_2,0,:,:],clevs,colors='black',linewidths=1.)
plt.clabel(cmslp, fontsize=9, inline=1, fmt = '%1.0f') 

cwind = m.pcolormesh(x,y,U10_2[t_2,-1,:,:],cmap=plt.cm.YlGnBu, vmax= Umax)
cb = m.colorbar(cwind,"right", size="5%", pad="2%")
cb.set_label('Wind speed 10 m (m/s)')



data_max = filters.minimum_filter(mslp_2[t_2,0,:,:], neighborhood_size)
data_max_wind = filters.maximum_filter(U10_2[t_2,-1,:,:], neighborhood_size2)
maxima_2 = (mslp_2[t_2,0,:,:] == data_max)

maxima_2[data_max > pthreshold]=0
maxima_2[data_max_wind < Uthreshold]= 0

xp_2, yp_2= m(lon[maxima_2], lat[maxima_2])
m.plot(xp_2, yp_2, 'ro')

print(mslp[t,0,:,:][maxima_2])


colors= ['#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061']
colors= colors[::-1]
boxnr= len(colors)
new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N= boxnr )

plt.suptitle('Time: ' +time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(tim[t])))


"""third plot"""
plt.figure(fignr+1)
plt.clf()

m = make_map(ncfile=nc)

diffmax= np.max(np.abs(U10_2[t_2,-1,:,:]-U10[t,-1,:,:]))

#plt.title('Control run - Operational')
cwind = m.pcolormesh(x,y,U10[t,-1,:,:]- U10_2[t_2,-1,:,:],cmap= new_map, vmin= -diffmax, vmax= diffmax)
cb = m.colorbar(cwind,"right", size="5%", pad="2%")
cb.set_label('Wind speed 10 m (m/s)')

m.plot(xp, yp, 'bo')
m.plot(xp_2, yp_2, 'ro')

#plt.savefig(outdir+"wind_ctr-oper"+str(time.localtime(tim[t])[3]), bbox_inches='tight', dpi= 60) #pad_inches=0.5,


"""fourth plot"""
plt.figure(fignr+1)
plt.clf()

diffmax= np.max(np.abs((mslp_2[t_2,0,:,:]-mslp[t,0,:,:])))
#plt.subplot(2, 2, 4)
m = make_map(ncfile=nc)

#plt.title('Control run - Operational')
cwind = m.pcolormesh(x,y,(mslp[t,0,:,:]- mslp_2[t_2,0,:,:]),cmap= new_map, vmin= -diffmax, vmax= diffmax)
cb = m.colorbar(cwind, "right", size="5%", pad="2%")
cb.set_label('Surface pressure difference (hPa)')

m.plot(xp, yp, 'bo')
m.plot(xp_2, yp_2, 'ro')
#title

#plt.savefig(outdir+"PSL_ctr-oper"+str(time.localtime(tim[t])[3]), bbox_inches='tight', dpi= 60) #pad_inches=0.5,

#plt.show()
# Save the figure
#fileout = "t2m.png"
#savefig(fileout, bbox_inches='tight', pad_inches=0.5)
