
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

outdir= "/home/'+user+'/home/Polar_Low/AromeArcticParameterization/Graphs/"
t= 0
fignr= 1

# Read in netcdf file
nc = NetCDFFile('/home/'+user+'/Data/pl/fc2015121300+000_EPS.nc')
#print(nc.variables.keys()) #print(\"Variables in data set: {}\".format(\", \".join(nc.variables.keys())))

tim = nc.variables['time'][:]
lat = nc.variables['latitude'][:]
lon = nc.variables['longitude'][:]
height= nc.variables['hybrid'][:]

# Variables
#T = nc.variables['air_temperature_pl'][:]
mslp = nc.variables['ga_pres_105'][:]/100
u10 = nc.variables['ga_u_105'][:]
v10 = nc.variables['ga_v_105'][:]
U10 = sqrt(u10**2+v10**2)



plt.figure(1)
plt.clf()



#plt.subplot(2, 2, 1)
m = make_map(ncfile=nc)
x,y = m(lon,lat)

#nc.close()

Umax= np.max(U10)

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


cwind = m.pcolormesh(x,y,U10[t,0,:,:],cmap=plt.cm.YlGnBu, vmax= Umax)
cb = m.colorbar(cwind,"right", size="5%", pad="2%")
cb.set_label('Wind speed 10 m (m/s)')


#
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
#
#
#
#"""operational data"""
#for t in [0]: #range(len(tim)):
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
#   
#    # Draw contours for mslp
#    clevs = np.arange(960,1080,3)
#    cmslp = m.contour(x,y,mslp_2[t_2,0, :,:],clevs,colors='black',linewidths=1.)
#    plt.clabel(cmslp, fontsize=9, inline=1, fmt = '%1.0f') 
#    
#    cwind = m.pcolormesh(x,y,U10_2[t,t,:,:],cmap=plt.cm.YlGnBu, vmax= Umax)
#    
##    neighborhood_size = 100
##    neighborhood_size2 = 50
##    pthreshold = 1010
##    Uthreshold = 15
##    
##    data_max = filters.minimum_filter(mslp_2[t_2,0,:,:], neighborhood_size)
##    data_max_wind = filters.maximum_filter(U10_2[t_2,-1,:,:], neighborhood_size2)
##    maxima = (mslp_2[t_2,0, :,:] == data_max)
##    
##    maxima[data_max > pthreshold]=0
##    maxima[lon < -10] = 0
##    
##    maxima[data_max_wind < Uthreshold]= 0
##    xp, yp= m(lon[maxima], lat[maxima])
##    #m.plot(xp, yp, 'go')
##    
##    yoffset = 0.050*(m.ymax-m.ymin)
##    for i in range(len(xp)): #range(2): #
###        if lon[maxima][i] < -10: i += 1
##
##        m.plot(xp[i], yp[i], 'ro')
##        plt.text(xp[i], yp[i]+yoffset, str(round(mslp_2[t_2,0,:,:][maxima][i], 1)), fontsize=13, fontweight= 'bold',
##                            ha='center',va='top',color='r',
##                            )#bbox = dict(boxstyle="square",ec='None',fc=(1,1,1,0.5)))
##    
##        plwind= np.round(data_max_wind[maxima][i],1)
##        print(plwind)
##        plt.text(xp[i], yp[i]-.5* yoffset, str(plwind), fontsize=13, fontweight= 'bold',
##                            ha='center',va='top',color='black',
##                            )#bbox = dict(boxstyle="square",ec='None',fc=(1,1,1,0.5)))
##    
##    
###    plt.savefig(outdir+"PSL_oper"+str(time.localtime(tim[t])[3]), bbox_inches='tight', dpi= 60) #pad_inches=0.5,
###    plt.title('MPSL controlrun')
##
##    cb = m.colorbar(cwind,"right", size="5%", pad="2%")
##    cb.set_label('Wind speed 10 m (m/s)')
##    
##    
##    print(np.round(mslp_2[t_2,0,:,:][maxima], 0))
##
##plt.figure(fignr+1)
##plt.clf()
##labelsize= 22
#
#cax = plt.axes([0.1, 0.2, 0.8, 0.05])
#cb= plt.colorbar(cwind, orientation= 'horizontal', cax= cax)
#cb.ax.tick_params(labelsize=labelsize-2)
#cb.set_label('Wind speed 10 m (m/s)', size= labelsize) 
#plt.savefig(outdir+'windbar', dpi= 50, bbox_inches='tight')