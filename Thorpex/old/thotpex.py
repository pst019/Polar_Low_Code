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
#flight='a'
#time='122309' #'102731'

flight='b'
time=['153245','154424','155225','161223','162442','163020','165104','170842','171817','172626','174236']

file=direct+date+flight+r'/D20'+date+'_'+time[1]+'QC.nc'
nc = NetCDFFile(file)


#print(nc.variables.keys()) #print(\"Variables in data set: {}\".format(\", \".join(nc.variables.keys())))
#ob= thorpex(nc) #this way all the relevant information can be imported

tim = nc.variables['time'][:]
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
alt= nc.variables['alt'][:]

# Variables
T = nc.variables['tdry'][:]
pres = nc.variables['pres'][:]
u = nc.variables['u_wind'][:]
v = nc.variables['v_wind'][:]
#U = sqrt(u10**2+v10**2)
U = nc.variables['wspd'][:] #wind speed
#Udir = nc.variables['wdir'][:]#wind direction

print(U[0:100], alt[0:100])

tim= tim[~U.mask]
lat= lat[~U.mask]
lon= lon[~U.mask]
alt= alt[~U.mask]
pres= pres[~U.mask]
U= U[~U.mask]


#tim= tim[~alt.mask]
#lat= lat[~alt.mask]
#lon= lon[~alt.mask]
#U= U[~alt.mask]
#pres= pres[~alt.mask]
#alt= alt[~alt.mask]

