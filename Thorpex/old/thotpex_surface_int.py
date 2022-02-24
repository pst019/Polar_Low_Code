#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 12:57:58 2016

@author: pst019

Plot

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
from fcts import *


# Read in netcdf file
direct=r'/home/'+user+'/Data/pl/IPY-THOTPEX_DROPSONDES_3-4-MARCH_2008/QC_'
date='080303'
flight='a'
times=['102731', '104912', '110458', '111154', '112742', '113414', '113922', '114338', 
'115202', '120832', '121359', '121815', '122309', '123128', '123624', '124345', '125351',
'130301', '131613', '132830']

#flight='b'
#times=['153245','154424','155225','161223','162442','163020',
#'165104','170842','171817','172626','174236']

tim=np.zeros((len(times)))
lat=np.zeros((len(times)))
lon=np.zeros((len(times)))
alt=np.zeros((len(times)))
T=np.zeros((len(times)))
pres=np.zeros((len(times)))
u=np.zeros((len(times)))
v=np.zeros((len(times)))
U=np.zeros((len(times)))
theta=np.zeros((len(times)))

for t in range(len(times)):

    file=direct+date+flight+r'/D20'+date+'_'+times[t]+'QC.nc'
    nc = NetCDFFile(file)
    
    
    #print(nc.variables.keys()) #print(\"Variables in data set: {}\".format(\", \".join(nc.variables.keys())))
    #ob= thorpex(nc) #this way all the relevant information can be imported
#    altr= nc.variables['alt'][:]
#    
#    alt[t]=altr[np.where(altr > 0)][0]
#    tim[t] = nc.variables['time'][:][np.where(altr > 0)][0]
#    lat[t] = nc.variables['lat'][:][np.where(altr > 0)][0]
#    lon[t] = nc.variables['lon'][:][np.where(altr > 0)][0]
#    
#    # Variables
#    T[t] = nc.variables['tdry'][:][np.where(altr > 0)][0]
#    pres[t] = nc.variables['pres'][:][np.where(altr > 0)][0]
#    u[t] = nc.variables['u_wind'][:][np.where(altr > 0)][0]
#    v[t] = nc.variables['v_wind'][:][np.where(altr > 0)][0]
#    #U = sqrt(u10**2+v10**2)
#    U[t] = nc.variables['wspd'][:][np.where(altr > 0)][0]#wind speed
#    #Udir = nc.variables['wdir'][:]#wind direction


    timr = nc.variables['time'][:]
    latr = nc.variables['lat'][:]
    lonr = nc.variables['lon'][:]
    altr= nc.variables['alt'][:]
    
    # Variables
    Tr = nc.variables['tdry'][:]
    presr = nc.variables['pres'][:]
    ur = nc.variables['u_wind'][:]
    vr = nc.variables['v_wind'][:]
    #U = sqrt(u10**2+v10**2)
    Ur = nc.variables['wspd'][:] #wind speed
    #Udir = nc.variables['wdir'][:]#wind direction
    theta= nc.variables['theta'][:] #potential temperature
    
#    print(U[0:100], alt[0:100])
    
#    tim= timr[~Ur.mask]
#    lat= latr[~Ur.mask]
#    lon= lonr[~Ur.mask]
#    alt= altr[~Ur.mask]
#    pres= presr[~Ur.mask]
#    U= Ur[~Ur.mask]
    
    """only values where Ur is measured"""
    altr= altr[~Ur.mask]
    presr= presr[~Ur.mask]
    tim[t]= timr[~Ur.mask][np.where(altr < 30)][0]
    lat[t]= latr[~Ur.mask][np.where(altr < 30)][0]
    lon[t]= lonr[~Ur.mask][np.where(altr < 30)][0]
    pres[t]= presr[~Ur.mask][np.where(altr < 30)][0]
    theta[t]= thetar[~Ur.mask][np.where(altr < 30)][0]
    U[t]= Ur[~Ur.mask][np.where(altr < 30)][0]

    alt[t]= altr[np.where(altr < 30)][0]


"""exclude values without latitude"""
alt=alt[np.isnan(lat)== False]
tim=tim[np.isnan(lat)== False][np.isnan(alt)== False]
pres=pres[np.isnan(lat)== False][np.isnan(alt)== False]
U=U[np.isnan(lat)== False][np.isnan(alt)== False]
lon=lon[np.isnan(lat)== False][np.isnan(alt)== False]
lat=lat[np.isnan(lat)== False][np.isnan(alt)== False]
alt= alt[np.isnan(alt)== False]

UTC= [str(time.localtime(t)[3])+':'+str(time.localtime(t)[4]).zfill(2) for t in tim]

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


xpt, ypt= m(lon, lat)
m.plot(xpt,ypt,'bo')

for i in range(len(UTC)):
    plt.text(xpt[i]+100, ypt[i]+100, UTC[i])


    
from scipy.interpolate import interp2d

plt.figure(1)

loni= np.linspace(min(lon), max(lon), 30)
lati= np.linspace(min(lat), max(lat), 30)

long, latg= np.meshgrid(loni, lati)

from matplotlib.mlab import griddata
presi= griddata(lon, lat, pres, long, latg, interp)#, interp='linear')

lonm, latm= m(long, latg)

clevs = np.arange(960,1080,1)
cmslp = m.contour(lonm,latm, presi, clevs,colors='black',linewidths=1.)
plt.clabel(cmslp, fontsize=9, inline=1, fmt = '%1.0f') 

m.pcolormesh(lonm,latm, presi)