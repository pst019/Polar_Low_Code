#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 15:06:41 2017

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import sys  #to import the functions from a different directory
#sys.path.insert(0, '/home/'+user+'/codeAROME/')
sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')

from f_imp_ERA2 import *


from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
import scipy.ndimage.filters as filters

year= 2010
#month= 3
#day= 3
#hour = 240#3 3 -6
month, day, hour= 1, 1, 6
t= ((day-1)*24+hour)//6 #hours since beginning of month

"""import data"""
d= data('Vort_128', year, month)

lon= d.lon128
lat= d.lat128
vort= d.vort128

vortfilter2= FourierFilter3D(d.lon128, d.lat128, d.vort128)


T_low= 40
T_up= 100


plt.figure(1)

plt.clf()
map= Polar_map(latsouth= 30)

grid= np.meshgrid(lon, lat) 
Lon, Lat= map(grid[0], grid[1])


PlotVort(Lon, Lat, vort[t], map)
plt.title('Vorticity')
plt.tight_layout()


plt.figure(2)
plt.clf()
map= Polar_map(latsouth= 30)

from SHfilter_Tuomas import *
vort_SH= SHfilter(vort[t], lmax= T_up)- SHfilter(vort[t], lmax= T_low)

PlotVort(Lon, Lat, vort_SH, map)
plt.title('Spherical harmonic filtered')


plt.tight_layout()




plt.figure(3)
plt.clf()
map= Polar_map(latsouth= 30)

d2= vort_128('Vort_128_comb_filt', year, 'A')

grid= np.meshgrid(d2.lon, d2.lat) 
Lon, Lat= map(grid[0], grid[1])

PlotVort(Lon, Lat, d2.vort[t*2], map) #*2 since 3 hourly in combined with fc
plt.title('Spherical harmonic filtered by algorithm')
plt.tight_layout()

