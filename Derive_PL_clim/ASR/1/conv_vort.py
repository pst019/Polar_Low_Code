#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 12:00:39 2017

@author: pst019
Calculate vorticity from UU and VV (not 100% sure if the map factor m is included right)

"""

fignr= 1

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.insert(0, '/home/'+user+'/codeAROME/')
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)

month= 1
day= 1

datadir= "/media/'+user+'/1692A00D929FEF8B/ASR/2000/subasr15km.anl.3D.2000"+str(month).zfill(2)+str(day).zfill(2)+".nc"
#datadir= "/media/'+user+'/1692A00D929FEF8B/ASR/asr30km.anl.3D.20000101.nc"
#datadir= "/media/'+user+'/1692A00D929FEF8B/PL/WRF/wrfout_d01_2008-03-03_00:00:00"

nc = Dataset(datadir)

print(nc.variables.keys())

u= nc.variables['UU'][0,0,:]
v= nc.variables['VV'][0,0,:]

lat= nc.variables['XLAT'][:]
lon= nc.variables['XLONG'][:]

r= 6371E3
DX= 15000

latr= np.deg2rad(lat)
lonr= np.deg2rad(lon)

dlatr= np.gradient(latr)
dlonr= np.gradient(lonr)

dlonr[0][dlonr[0] > np.pi/2] -= np.pi
dlonr[1][dlonr[1] > np.pi/2] -= np.pi
  
dlonr[0][dlonr[0] < -np.pi/2] += np.pi
dlonr[1][dlonr[1] < -np.pi/2] += np.pi
     
xn= dlonr[0] *np.cos(latr)
yn= dlatr[0]
mx= DX/(r*np.sqrt(xn**2 + yn**2))

xn= dlonr[1] *np.cos(latr)
yn= dlatr[1]
my= DX/(r*np.sqrt(xn**2 + yn**2))

#plt.figure(fignr)
#plt.clf()
#fignr+= 1
#
#cs= plt.contour(np.arange(720), np.arange(720), mx)
#cb = plt.clabel(cs ) #, "right", size="5%", pad="2%")

#plt.figure(fignr)
#plt.clf()
#fignr+= 1
#cs= plt.contour(np.arange(720), np.arange(720), my)
#cb = plt.clabel(cs ) #, "right", size="5%", pad="2%")


m= mx

u_y= np.gradient(u, DX)[1]
v_x= np.gradient(v, DX)[0]

vort=( v_x - u_y)* m

plt.figure(fignr)
plt.clf()
fignr+= 1
cs= plt.pcolormesh(np.arange(720), np.arange(720), vort) #*1E5)
cb = plt.colorbar(cs ) #, "right", size="5%", pad="2%")

#vort=( v_x - u_y)
#
#plt.figure(fignr)
#plt.clf()
#fignr+= 1
#cs= plt.pcolormesh(np.arange(720), np.arange(720), vort*1E5)
#cb = plt.colorbar(cs ) #, "right", size="5%", pad="2%")