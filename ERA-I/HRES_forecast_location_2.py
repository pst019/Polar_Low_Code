#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 13:25:51 2017

@author: pst019
from ASR-ERA_Inv_random_3 of folder ASR
"""


import sys
import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'

import matplotlib.pyplot as plt

sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')

from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_meteo import *
from netCDF4 import Dataset, num2date
import scipy.ndimage.filters as filters
import numpy as np

"""global variables"""

fignr= 5


plt.figure(fignr, figsize= (4.5, 4.5))
fignr += 1
plt.clf()

#map = AA_map_half()
map = Lambert_map(lllon= -5, lllat= 68, urlon=30, urlat= 70, lat0= 77.5, lon0= -25, res='i', fill=False, coastline=True, latdist=5, londist=10)


"""time"""
year=2008
month=3
filedaylist= [2, 3, 3, 4, 4]
filehourlist= [12, 0, 12, 0, 12]


point_day, point_hour= 4, 6
point_day2, point_hour2= 3, 12


ilevel= 1 #0 -950, 1- 850, 2-700, 3-500

propdistance= 100E3 *3
distance= int(500/10)  #180  100   #make distance reasonable larger than propdistance to avoid "doublematches" twice is on the safe side

latbound= [60, 71] #60, 69
lonbound= [-5, 12] #-5, 9
threshold= 8 #10
#timebound=  [8, 12]

save= True
#save=False
savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/DifferentTracks/'


for fi, fileday in enumerate(filedaylist):
    filehour= filehourlist[fi]
    
    
    timebound= list(np.array([(3- fileday)* 24 + (20- filehour), (3- fileday)* 24 + (36- filehour)])//3 ) # the time in which the PL has to propage through the box
    tpoint= ((point_day- fileday)* 24 + (point_hour - filehour) )//3
    tpoint2= ((point_day2- fileday)* 24 + (point_hour2 - filehour) ) //3

    
    """import the data"""
#    file= Mediadir + 'ECMWF/HRES/oper_SFC_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+'_'+str(filehour).zfill(2)+'.nc'
#    nc= Dataset(file)
#    SLP= nc.variables['msl'][tpoint]/100    
    
    
    file= Mediadir + 'ECMWF/HRES/oper_PL_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+'_'+str(filehour).zfill(2)+'.nc'
    
    nc= Dataset(file)
    print(nc.variables.keys())
    varlist= list(nc.variables.keys())
    
                        
    tim= nc.variables['time'][:]
    dtime= num2date(tim, nc.variables['time'].units)
    
    
    if 'latitude' in varlist:
        lat= nc.variables['latitude'][:]
        lon= nc.variables['longitude'][:]
        level= int(nc.variables['level'][ilevel])
    
    #    if file== '/media/pst019/1692A00D929FEF8B/ECMWF/HRES/oper_PL_20080302_12.nc' or file== '/media/pst019/1692A00D929FEF8B/ECMWF/HRES/oper_PL_20080304_12.nc':
    else:
        lat= nc.variables['lat'][:]
        lon= nc.variables['lon'][:]
        level= int(nc.variables['plev'][ilevel]//100)
       
                                 
    u= nc.variables['u'][:, ilevel]
    v= nc.variables['v'][:, ilevel]
    
    
    dtimemax= datetime.datetime(2008, 3, 5)
    if np.max(dtime) > dtimemax:
        tmax= np.where(dtime == dtimemax)[0][0]
        dtime= dtime[:tmax]
        u= u[:tmax]
        v= v[:tmax]
        
        
    #print(dtime)
    
    
     
    dy= (lat[0]- lat[1])*110E3 #distance along axis 1
    dx= dy* np.cos(np.deg2rad(lat)) #distance along axis 0
    
    #axis -1: latitude = y-axis  -   axis -2: longitude = x-axis
    vort= np.gradient(v, axis= -1)/dx[:, np.newaxis] - np.gradient(u, -dy, axis= -2) #-dy since the latitudes start from the north  
    #    gausfilterdist= 100
    vort= filters.gaussian_filter1d(filters.gaussian_filter1d(vort, sigma= 20, mode='nearest', truncate= 4, axis= -1), sigma= 10, mode='nearest', truncate= 4, axis= -2) #smooth more along the longitude axis since it is on .1*.1 degree grid.
    
    
    lon_gr, lat_gr= np.meshgrid(lon, lat)
          
    
    """do the tracking"""
    color = next(plt.gca()._get_lines.prop_cycler)['color']
    if fi == 4: color = next(plt.gca()._get_lines.prop_cycler)['color'] # to jump over HRES-04-06 which does not exist (for same colors as AROME)
    #    data_max = filters.maximum_filter1d(data, distance, axis= -1)
    #    data_max = filters.maximum_filter1d(data_max, distance, axis= -2)
    
    vortlist, lonlist, latlist, labellist = EasyTrack(vort*1E5, lon_gr, lat_gr, distance,
                                threshold, propdistance)
    vortlist, lonlist, latlist, labellist= Tracks_in_lat_lon(vortlist, lonlist, latlist, labellist, latbound, lonbound, list(timebound))
    
    
#    t= np.max([tpoint, 0])
#    t= np.max([timebound[0], 0])
#    cyclnratt= labellist[t]
    
    cyclnratt= remove_dublicate2D(list(labellist))

    print(cyclnratt)

    for ci, cyclnr in enumerate(cyclnratt):
        tcycl, vortcycl, loncycl, latcycl= Data_One_Track(vortlist, lonlist, latlist, labellist, cyclnr)
        xpt, ypt= map(loncycl, latcycl)
        map.plot(xpt,ypt, color=color)
    
    
    
        if ci== 0: map.plot(xpt[tcycl==tpoint],ypt[tcycl==tpoint],'s', color=color, label= 'HRES-'+str(fileday).zfill(2)+'-'+str(filehour).zfill(2), markersize=7, zorder= 5)
        else: map.plot(xpt[tcycl==tpoint],ypt[tcycl==tpoint],'s', color=color, markersize=7, zorder= 5)

        if tpoint2 >= 0: map.plot(xpt[tcycl==tpoint2],ypt[tcycl==tpoint2],'o', color=color, markersize=8, zorder= 5)

        map.plot(xpt[tcycl==0],ypt[tcycl==0],'x', color='k', zorder= 5)


plt.legend(loc= 1)

plt.tight_layout()

if save==True:
    print(savedir+'HRES_tracks_analysis_mark')
    plt.savefig(savedir+'HRES_tracks_analysis_mark')

#"""plot for t"""
##grid=  #pcolormesh needs corner points for projection
#Lon, Lat= map(lon_gr, lat_gr)
#
#
#PlotContours(Lon, Lat, SLP, map, leveldist= 2, alpha= .7, numbers= False)
#
#PlotColorMap4(Lon, Lat, vort[tpoint]*1E5, map, color='RdBu', symetric=True)
#
#PlotLocalMax(vort[tpoint]*1E5, threshold=threshold, distance=distance, map= map, lon=lon, lat=lat, typ='max',
#                     color='orange', dot=True, roundorder= 0, latbound=[67, 80])
