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

fignr= 1
# Read in netcdf file
EXP="40h11_AA"
DATE="2015121300"
file=r"/home/'+user+'/Data/pl/"+EXP+"_"+DATE+"_fp_extract.nc"
nc = NetCDFFile(file)

#print(nc.variables.keys()) #print(\"Variables in data set: {}\".format(\", \".join(nc.variables.keys())))
#ctr= Arome(nc)

t= 0
print('time', str(time.localtime(ctr.tim[t])[3]))


plt.figure(fignr)
fignr+=1
plt.clf()

m = make_map(ncfile=nc)
x,y = m(ctr.lon,ctr.lat)

#plot surface pressure
clevs = np.arange(960,1080,3)
cmslp = m.contour(x,y,ctr.mslp[t,0, :,:],clevs,colors='black',linewidths=1.)
plt.clabel(cmslp, fontsize=9, inline=1, fmt = '%1.0f') 

#plot wind
Umax= np.max(ctr.U[:,-1])
cwind = m.pcolormesh(x,y,ctr.U[t,-1,:,:],cmap=plt.cm.YlGnBu, vmax= Umax)
cb = m.colorbar(cwind,"right", size="5%", pad="2%")
cb.set_label('Wind speed 10 m (m/s)')