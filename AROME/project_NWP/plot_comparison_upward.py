
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
import scipy.ndimage.filters as filters


# import own modules
from mapplot_def import *

outdir= "/home/'+user+'/home/Polar_Low/AromeArcticParameterization/Graphs/"
t= 1
fignr= 1

## Read in netcdf file
##nc = NetCDFFile('/home/patricks/Data/fc2015121300_fp.nc') #at laptop
nc = NetCDFFile('/home/'+user+'/Data/pl/fc2015121300_ctr_fp.nc')
##nc = NetCDFFile('/home/'+user+'/Data/pl/fc2015121300_FLX_fp.nc')
#
##print(nc.variables.keys()) #print(\"Variables in data set: {}\".format(\", \".join(nc.variables.keys())))
#
tim = nc.variables['time'][:]
lat = nc.variables['latitude'][:]
lon = nc.variables['longitude'][:]
height= nc.variables['height_above_msl'][:]

# Variables
#T = nc.variables['air_temperature_pl'][:]
mslp = nc.variables['air_pressure_at_sea_level'][:]/100
u10 = nc.variables['x_wind_pl'][:]
v10 = nc.variables['y_wind_pl'][:]
w = nc.variables['upward_air_velocity_pl'][:]

U10 = sqrt(u10**2+v10**2)
#
##ncsfx = NetCDFFile('/home/'+user+'/Data/pl/fc2015121300_ctr_sfx.nc')
##SST = ncsfx.variables['SST'][:]
#
"""file 2"""
##nc2 = NetCDFFile('/home/patricks/Data/fc2015121300_SST_fp.nc') #at laptop
nc2 = NetCDFFile('/home/'+user+'/Data/pl/fc2015121300_FLX_fp.nc')
#nc2 = NetCDFFile('/home/'+user+'/Data/pl/arome_arctic_extracted_2_5km_20151213T00Z.nc')
#
## Variablesnc2 = NetCDFFile('/home/'+user+'/Data/pl/fc2015121300_FLX_fp.nc')
#
##T2 = nc2.variables['air_temperature_pl'][:]
#
tim_2 = nc2.variables['time'][:]
mslp_2 = nc2.variables['air_pressure_at_sea_level'][:]/100 
u10_2 = nc2.variables['x_wind_pl'][:]   
v10_2 = nc2.variables['y_wind_pl'][:]
w_2 = nc2.variables['upward_air_velocity_pl'][:]
#theta_2 = nc2.variables['relative_humidity_pl'][:]
#lh_2 = nc2.variables['integral_of_surface_downward_latent_heat_evaporation_flux_wrt_time'][:, 0]*-1 #does not include sublimation
#
U10_2 = sqrt(u10_2**2+v10_2**2)


nc2.close()
##ncsfx_2 = NetCDFFile('/home/'+user+'/Data/pl/fc2015121300_SST_sfx.nc')
##SST_2 = ncsfx_2.variables['SST'][:]


t_2= np.where(tim_2 == tim[t])[0][0]



#for t in range(len(tim)):
#    print('time', str(time.localtime(tim[t])[3]))
#
#    # Make the map
#    plt.figure(fignr)
#    fignr+=1
#    plt.clf()
#    
#    """first plot"""
#    #plt.subplot(2, 2, 1)
#    m = make_map(ncfile=nc)
#    x,y = m(lon,lat)
#    
#    t_2= np.where(tim_2 == tim[t])[0][0]
##    Umax= np.max(U10[:,-1])
#    
#    # Draw contours for mslp
#    clevs = np.arange(960,1080,3)
#    cmslp = m.contour(x,y,mslp[t,0, :,:],clevs,colors='black',linewidths=1.)
#    plt.clabel(cmslp, fontsize=9, inline=1, fmt = '%1.0f') 
#    
#    cwind = m.pcolormesh(x,y,U10[t,-1,:,:],cmap=plt.cm.YlGnBu, vmax= Umax)
#    
#    neighborhood_size = 100
#    neighborhood_size2 = 50
#    pthreshold = 1010
#    Uthreshold = 15
#    
#    data_max = filters.minimum_filter(mslp[t,0,:,:], neighborhood_size)
#    data_max_wind = filters.maximum_filter(U10[t,-1,:,:], neighborhood_size2)
#    maxima = (mslp[t,0, :,:] == data_max)
#    
#    maxima[data_max > pthreshold]=0
#    maxima[lon < -10] = 0
#    
#    maxima[data_max_wind < Uthreshold]= 0
#    xp, yp= m(lon[maxima], lat[maxima])
#    #m.plot(xp, yp, 'go')
#    
#    yoffset = 0.050*(m.ymax-m.ymin)
#    for i in range(len(xp)): #range(2): #
##        if lon[maxima][i] < -10: i += 1
#
#        m.plot(xp[i], yp[i], 'ro')
#        plt.text(xp[i], yp[i]+yoffset, str(round(mslp[t,0,:,:][maxima][i], 1)), fontsize=13, fontweight= 'bold',
#                            ha='center',va='top',color='r',
#                            )#bbox = dict(boxstyle="square",ec='None',fc=(1,1,1,0.5)))
#    
#        plwind= np.round(data_max_wind[maxima][i],1)
#        print(plwind)
#        plt.text(xp[i], yp[i]-.5* yoffset, str(plwind), fontsize=13, fontweight= 'bold',
#                            ha='center',va='top',color='black',
#                            )#bbox = dict(boxstyle="square",ec='None',fc=(1,1,1,0.5)))
#    
#    
#    plt.savefig(outdir+"PSL_FLX"+str(time.localtime(tim[t])[3]), bbox_inches='tight', dpi= 60) #pad_inches=0.5,
##    plt.title('MPSL controlrun')
#
#    cb = m.colorbar(cwind,"right", size="5%", pad="2%")
#    cb.set_label('Wind speed 10 m (m/s)')
#    
#    
#    print(np.round(mslp[t,0,:,:][maxima], 0))


colors= ['#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061']
colors= colors[::-1]
boxnr= len(colors)
new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N= boxnr )

"""first"""
plt.figure(fignr)
fignr+=1
plt.clf()

m = make_map(ncfile=nc)
x,y = m(lon,lat)

h= -3 #-3 is 850 hPa
wmax= np.max(np.abs(w[t,h,:,:]))/3
cwind = m.pcolormesh(x,y,w[t,h,:,:],cmap=new_map, vmax= wmax, vmin= -wmax)
cb = m.colorbar(cwind,"right", size="5%", pad="2%")
cb.set_label('Vertical wind speed (m/s)')

neighborhood_size = 100
neighborhood_size2 = 50
pthreshold = 1010
Uthreshold = 15

data_max = filters.minimum_filter(mslp[t,0,:,:], neighborhood_size)
data_max_wind = filters.maximum_filter(U10[t,-1,:,:], neighborhood_size2)
maxima = (mslp[t,0, :,:] == data_max)

maxima[data_max > pthreshold]=0
maxima[lon < -10] = 0

maxima[data_max_wind < Uthreshold]= 0
xp, yp= m(lon[maxima], lat[maxima])
m.plot(xp, yp, 'ro')
plt.savefig(outdir+"W_ctr"+str(time.localtime(tim[t])[3]), bbox_inches='tight', dpi= 60) 

"""second"""
plt.figure(fignr)
fignr+=1
plt.clf()

m = make_map(ncfile=nc)
x,y = m(lon,lat)

h= -3 #-3 is 850 hPa
#wmax= np.max(np.abs(w_2[t,h,:,:]))/2
cwind = m.pcolormesh(x,y,w_2[t,h,:,:],cmap=new_map, vmax= wmax, vmin= -wmax)
cb = m.colorbar(cwind,"right", size="5%", pad="2%")
cb.set_label('Vertical wind speed (m/s)')

data_max = filters.minimum_filter(mslp_2[t,0,:,:], neighborhood_size)
data_max_wind = filters.maximum_filter(U10_2[t,-1,:,:], neighborhood_size2)
maxima = (mslp_2[t,0, :,:] == data_max)

maxima[data_max > pthreshold]=0
maxima[lon < -10] = 0

maxima[data_max_wind < Uthreshold]= 0
xp, yp= m(lon[maxima], lat[maxima])
m.plot(xp, yp, 'ro')

plt.savefig(outdir+"W_sens"+str(time.localtime(tim[t])[3]), bbox_inches='tight', dpi= 60) 


"""third"""
plt.figure(fignr)
fignr+=1
plt.clf()

m = make_map(ncfile=nc)
x,y = m(lon,lat)



#m = make_map(ncfile=nc)
#x,y = m(lon,lat)
#lhmax= np.max(np.abs(lh_2[t_2,:,:]))
#cwind = m.pcolormesh(x,y,lh_2[t_2,:,:],cmap=new_map, vmax= lhmax, vmin= -lhmax)
#cb = m.colorbar(cwind,"right", size="5%", pad="2%")
#cb.set_label('Wind speed (m/s)')

#"""operational data"""
#for t in [1]: #range(len(tim)):
#    print('time', str(time.localtime(tim[t])[3]))
#    t_2= np.where(tim_2 == tim[t])[0][0]
#
#    # Make the map
#    plt.figure(fignr)
#    fignr+=1
#    plt.clf()
#    
#    """first plot"""
#    #plt.subplot(2, 2, 1)
#    m = make_map(ncfile=nc)
#    x,y = m(lon,lat)
#    
#   
#    # Draw contours for mslp
#    clevs = np.arange(960,1080,3)
#    cmslp = m.contour(x,y,mslp_2[t_2,0, :,:],clevs,colors='black',linewidths=1.)
#    plt.clabel(cmslp, fontsize=9, inline=1, fmt = '%1.0f') 
#    
#    cwind = m.pcolormesh(x,y,U10_2[t_2,-1,:,:],cmap=plt.cm.YlGnBu, vmax= Umax)
#    
#    neighborhood_size = 100
#    neighborhood_size2 = 50
#    pthreshold = 1010
#    Uthreshold = 15
#    
#    data_max = filters.minimum_filter(mslp_2[t_2,0,:,:], neighborhood_size)
#    data_max_wind = filters.maximum_filter(U10_2[t_2,-1,:,:], neighborhood_size2)
#    maxima = (mslp_2[t_2,0, :,:] == data_max)
#    
#    maxima[data_max > pthreshold]=0
#    maxima[lon < -10] = 0
#    
#    maxima[data_max_wind < Uthreshold]= 0
#    xp, yp= m(lon[maxima], lat[maxima])
#    #m.plot(xp, yp, 'go')
#    
#    yoffset = 0.050*(m.ymax-m.ymin)
#    for i in range(len(xp)): #range(2): #
##        if lon[maxima][i] < -10: i += 1
#
#        m.plot(xp[i], yp[i], 'ro')
#        plt.text(xp[i], yp[i]+yoffset, str(round(mslp_2[t_2,0,:,:][maxima][i], 1)), fontsize=13, fontweight= 'bold',
#                            ha='center',va='top',color='r',
#                            )#bbox = dict(boxstyle="square",ec='None',fc=(1,1,1,0.5)))
#    
#        plwind= np.round(data_max_wind[maxima][i],1)
#        print(plwind)
#        plt.text(xp[i], yp[i]-.5* yoffset, str(plwind), fontsize=13, fontweight= 'bold',
#                            ha='center',va='top',color='black',
#                            )#bbox = dict(boxstyle="square",ec='None',fc=(1,1,1,0.5)))
#    
#    
##    plt.savefig(outdir+"PSL_oper"+str(time.localtime(tim[t])[3]), bbox_inches='tight', dpi= 60) #pad_inches=0.5,
##    plt.title('MPSL controlrun')
#
#    cb = m.colorbar(cwind,"right", size="5%", pad="2%")
#    cb.set_label('Wind speed 10 m (m/s)')
#    
#    
#    print(np.round(mslp_2[t_2,0,:,:][maxima], 0))
#
#plt.figure(fignr+1)
#plt.clf()
#labelsize= 22
#
#cax = plt.axes([0.1, 0.2, 0.8, 0.05])
#cb= plt.colorbar(cwind, orientation= 'horizontal', cax= cax)
#cb.ax.tick_params(labelsize=labelsize-2)
#cb.set_label('Wind speed 10 m (m/s)', size= labelsize) 
##plt.savefig(outdir+'windbar', dpi= 50, bbox_inches='tight')