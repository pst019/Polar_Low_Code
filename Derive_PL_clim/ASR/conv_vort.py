#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 12:00:39 2017

@author: pst019
"""

fignr= 0

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, '/home/'+user+'/polar_low_code/Functions/')
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)

Mediadir= '/media/'+user+'/1692A00D929FEF8B/'

datadir= Mediadir+ 'ASR/2000/subasr15km.anl.3D.20001231.nc'

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
mx= DX/(r*np.sqrt(xn**2 + yn**2))  #some mapping factor can be shown by: #PlotColorMap(Lon, Lat, my, map, variable='my')


u_y= np.gradient(u, DX)[1]
v_x= np.gradient(v, DX)[0]

vort=( v_x - u_y)* mx

plt.figure(fignr)
plt.clf()
fignr+= 1

#map= Polar_map(latsouth= 30)
#map= WRF_map(lat, lon, lat0= 85, lon0= 10)
map= ASR_map()
Lon, Lat= map(lon, lat)

PlotColorMap(Lon, Lat, vort*1E5, map, variable='vorticity')