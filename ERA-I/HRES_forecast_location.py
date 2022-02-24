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

fignr= 4


plt.figure(fignr, figsize= (6, 4.7))
fignr += 1
plt.clf()

map = AA_map_half()


"""time"""
year=2008
month=3
fileday= 2
filehour= 12


#point_day, point_hour= 4, 6
point_day, point_hour= 3, 30


ilevel= 1 #0 -950, 1- 850, 2-700, 3-500

propdistance= 100E3 *3
distance= int(500/10) #800 #180  100   #make distance reasonable larger than propdistance to avoid "doublematches" twice is on the safe side

latbound= [60, 71] #60, 69
lonbound= [-5, 12] #-5, 9
threshold= 8 #10
#timebound=  [8, 12]

timebound= list(np.array([(3- fileday)* 24 + (20- filehour), (3- fileday)* 24 + (36- filehour)])//3 ) # the time in which the PL has to propage through the box
filterrad_latcells= 10 #*.1 (spacing of latitude)* 110km , so 10 is approximately a filter of 100km

#save= True
save=False
savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/DifferentTracks/'


tpoint= ((point_day- fileday)* 24 + (point_hour - filehour) )//3


"""import the data"""
file= Mediadir + 'ECMWF/HRES/oper_SFC_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+'_'+str(filehour).zfill(2)+'.nc'
nc= Dataset(file)
SLP= nc.variables['msl'][tpoint]/100    


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
vort= filters.gaussian_filter1d(filters.gaussian_filter1d(vort, sigma= 2* filterrad_latcells, mode='nearest', truncate= 4, axis= -1), sigma= filterrad_latcells, mode='nearest', truncate= 4, axis= -2) #smooth more along the longitude axis since it is on .1*.1 degree grid.


lon_gr, lat_gr= np.meshgrid(lon, lat)


"""plot for t"""
#grid=  #pcolormesh needs corner points for projection
Lon, Lat= map(lon_gr, lat_gr)


PlotContours(Lon, Lat, SLP, map, leveldist= 2, alpha= .7, numbers= False)

PlotColorMap4(Lon, Lat, vort[tpoint]*1E5, map, color='RdBu', symetric=True)

PlotLocalMax(vort[tpoint]*1E5, threshold=threshold, distance=distance, map= map, lon=lon, lat=lat, typ='max',
                     color='orange', dot=True, roundorder= 0, latbound=[67, 80])


"""do the tracking"""
color = next(plt.gca()._get_lines.prop_cycler)['color']
#    data_max = filters.maximum_filter1d(data, distance, axis= -1)
#    data_max = filters.maximum_filter1d(data_max, distance, axis= -2)

vortlist, lonlist, latlist, labellist = EasyTrack(vort*1E5, lon_gr, lat_gr, distance,
                            threshold, propdistance)
vortlist, lonlist, latlist, labellist= Tracks_in_lat_lon(vortlist, lonlist, latlist, labellist, latbound, lonbound, timebound)

#print(labellist)
#t= np.max([timebound[0], 0])
#cyclnratt= labellist[t]
cyclnratt= remove_dublicate2D(list(labellist))

for cyclnr in cyclnratt:
    tcycl, vortcycl, loncycl, latcycl= Data_One_Track(vortlist, lonlist, latlist, labellist, cyclnr)
    xpt, ypt= map(loncycl, latcycl)
    map.plot(xpt,ypt, color=color)



    map.plot(xpt[tcycl==tpoint],ypt[tcycl==tpoint],'s', color=color, label= 'HRES-'+str(fileday).zfill(2)+'-'+str(filehour).zfill(2), zorder= 5)

if save==True:
    plt.savefig(savedir+'HRES_signle_track')