#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

This is from plot_comparison
"""
import os
user = os.getcwd().split('/')[2]


import matplotlib.pyplot as plt
import numpy as np
import time
from scipy import stats

# import own modules
sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')
from f_imp_ERA2 import *
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)

year= 1987
month= 2

lllon, lllat, urlon, urlat, lat0, lon0 = -18, 55, 80, 72, 60, 0

lllon, lllat, urlon, urlat= -15, 65, 50, 75
lat0, lon0= 75, 0 #(lat0, lon0) = center point


for day in range(27,28):
#    day= 23
    hour= 0
    tim= ((day-1)*24+hour)//6 
    
    d= data(['MSLP'], year, month, tstart= tim, tend= tim+1)
    d.imp_u10()
    
    
    """MSLP, Surface temperature and surface winds"""
    plt.figure(1)
    plt.clf()
    
    
    #map = AA_map()
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)   
    #map= Polar_map(latsouth= d.lat[-1])
    grid= np.meshgrid(d.lon, d.lat)
    Lon,Lat = map(grid[0], grid[1])
    U= np.sqrt(d.u10**2 + d.v10**2)
    
    # Draw contours for mslp
    PlotContours(Lon, Lat, d.MSLP[0], map, leveldist= 1)
    #
    
    PlotWind(d.lon, d.lat, d.u10[0], d.v10[0], map)
    PlotWindVelo(Lon, Lat, U[0], map, Umax= 25)#, bar= False)


#    directory='/home/'+user+'/home/Graphs_Martin/'
#    plt.tight_layout()
#    plt.savefig(directory+ 'MSLP_Wind_'+str(year)+'_'+str(month)+'_'+str(day)+'_'+str(hour)+'.png')#, dpi= 70, bbox_inches='tight')
