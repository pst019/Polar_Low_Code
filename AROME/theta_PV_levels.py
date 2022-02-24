#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

This is from plot_comparison
"""

import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'



from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})


import numpy as np
import time
#from datetime import date, datetime, timedelta
from scipy import stats
import scipy.ndimage.filters as filters

# import own modules
import sys  #to import the functions from a different directory
sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')

from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_imp_AROME import *  # Read in netcdf file
import f_imp_thorpex as thorpex
from f_meteo import *
from f_useful import *


"""global variables"""
fignr= 4

#maptype='AA'
maptype='AA_half'
#maptype='Lambert'

if maptype=='AA': plt.figure(fignr, figsize= (6, 6))
elif maptype=='AA_half': plt.figure(fignr, figsize= (6, 4.7))
else: plt.figure(fignr, figsize= (6, 4.5))
plt.clf()

fignr+=1
 
"""time specification: time of the plot"""
year, month = 2008, 3
day, hour= 3, 12
    
"""name of the experiment"""

exp_name, fileday, filehour='DA_080303_CTR', 3, 0
exp_name, fileday, filehour='DA_080301_cycling', 4, 12

"""time of the file"""
#fileday, filehour= 4, 0        
#t= (day- fileday)*24 + (hour- filehour)  # -1 
t= 0

#year, month = 1987, 2
#day, hour= 26, 12    
#
#fileday, filehour= 25, 12         
#t= (day- fileday)*24 + (hour- filehour)  # -1 

tlist= [t]
#tlist= [19] #np.arange(0, 3)#42, 3)                                            


"""c) plot characteristics"""

var='tropopause_advanced'


""" presure level (if the variable ends with '_var' """
pn = 4
   

save=False
#save= True

AAres=2 #2- every second datapoint is taken


title_extra=''

"""Lambert coordinates"""


#to match with Thorpex
if maptype=='Lambert':
    if year == 2008 and month== 3 and day== 3 and hour<15:
        lllon, lllat, urlon, urlat= -4, 70, 19, 75
        lat0, lon0= 75, 0 #(lat0, lon0) = center point
    
    elif year == 2008 and month== 3 and day== 3 and hour>=15:
        lllon, lllat, urlon, urlat= -6, 69.5, 13, 74.5
        lat0, lon0= 75, 0 #(lat0, lon0) = center point
    
    elif year == 2008 and month== 3 and day== 4: 
        if var=='CTT':
            lllon, lllat, urlon, urlat= -13, 69, 12, 73
            lat0, lon0= 70, -20 #(lat0, lon0) = center point 
            title_extra +='_zoom'
        elif 'DA' in exp_name: #for original AA domain
            
            lllon, lllat, urlon, urlat= -0.5, 67.5, 12.5, 68.8
            lat0, lon0= 68, -25 #(lat0, lon0) = center point 
           
        else:  #for AA domain moved south
            lllon, lllat, urlon, urlat= -3, 63.5, 15, 69.5
            lat0, lon0= 70, 0 #(lat0, lon0) = center point 


if 'noFLXarea' in exp_name:
    boxlon, boxlat= [0, 10, 10, 0, 0], [67.9, 67.9, 74, 74, 68] #to display the area without fluxes

elif year == 2008 and month== 3 and day== 3 and hour<15:
    boxlon, boxlat= [-4, 13, 18.5, -4, -4], [70, 69.3, 75, 76, 70] #to display the first lambert map
elif year == 2008 and month== 3 and day== 3 and hour>=15:
    boxlon, boxlat= [-5, 8, 10, -6.5, -5], [69.5, 69.5, 74.5, 74.5, 69.5]   # lllon, lllat, urlon, urlat= -5, 69.5, 10, 74.5 to display the second lamber map
elif year == 2008 and month== 3 and day== 4: 
    if 'DA' in exp_name: #for original AA domain
        boxlon, boxlat= [-0.8, 8, 13, 3, -.8], [67.8, 65.9, 68.8, 70.9, 67.8] #lllon, lllat, urlon, urlat= -3, 63.5, 15, 69.5 #to display the third lambert map
    else:
        boxlon, boxlat= [-3, 12, 15, -3, -3], [63.5, 63.5, 70, 70, 63.5] #lllon, lllat, urlon, urlat= -3, 63.5, 15, 69.5 #to display the third lambert map










if maptype== 'AA': map = AA_map()
elif maptype== 'AA_half': map = AA_map_half()
else: map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0, latdist=2, londist= 5)    


AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
                                                 
d= data(filename= AAfilename, res= AAres)

Lon,Lat = map(d.lon,d.lat)
    
d.imp_standard_t(tn= t)   


plt.figure(fignr)
plt.clf()

locx, locy= 15, 15


plt.plot(d.PV[:, locx, locy]*1E6, d.pres)
plt.ylim([1000 , 100])






from scipy.interpolate import interp1d        

inc= 20
pres_int= np.arange(100, 700, inc)


PV_int= np.moveaxis( np.array([[ interp1d(d.pres, d.PV[:, x, y])(pres_int) for y in range(d.PV.shape[2]) ] for x in range(d.PV.shape[1]) ]  ), -1, 0)


plt.plot(PV_int[:, locx, locy]*1E6, pres_int)

pres_int_total= pres_int[:, np.newaxis, np.newaxis]* np.ones(np.shape(PV_int) )
pres_int_total[PV_int < 2*1E-6]= np.nan
tropheight= np.nanmax(pres_int_total, axis= 0) +inc/2

plt.plot(2, tropheight[locx, locy], 'o')


plt.figure(fignr-1)

gausfilterdist= 20
ex= 0
tropheight_filt= filters.gaussian_filter(tropheight[ex:], sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)


bounds= np.arange(200, 561, 20)
PlotColorMap4(Lon[ex:], Lat[ex:], tropheight_filt, map, bounds= bounds, color='viridis_r', label='Tropopause height [hPa]')


"""the tropopause wind"""
#U = np.sqrt(d.u**2+ d.v**2)
u_int= np.moveaxis( np.array([[ interp1d(d.pres, d.u[:, x, y])(pres_int) for y in range(d.PV.shape[2]) ] for x in range(d.PV.shape[1]) ]  ), -1, 0)[ex:]
v_int= np.moveaxis( np.array([[ interp1d(d.pres, d.v[:, x, y])(pres_int) for y in range(d.PV.shape[2]) ] for x in range(d.PV.shape[1]) ]  ), -1, 0)[ex:]

ex= 1
u_trop= u_int[:, ex:-ex, ex:-ex][tropheight[ex: -ex, ex: -ex] - inc/2 == pres_int[:, np.newaxis, np.newaxis]].reshape(np.shape(tropheight[ex: -ex, ex: -ex]).T)
v_trop= v_int[:, ex:-ex, ex:-ex][tropheight[ex: -ex, ex: -ex] - inc/2 == pres_int[:, np.newaxis, np.newaxis]].reshape(np.shape(tropheight[ex: -ex, ex: -ex]).T)

PlotWind(d.lon[ex: -ex, ex: -ex], d.lat[ex: -ex, ex: -ex], u_trop, v_trop, map, rot=False)

print('something is wrong, maybe with the transformation')

#PlotContours(Lon, Lat, tropheight, map, leveldist= 50)


#d.imp_surf(tn= t)
#
#if var in ['FLX', 'LH', 'LH_thor', 'TH', 'TH_thor', 'PBH', 'PBH_advanced', 'CAPE', 'CIN', 'LCL', 'LFC', 'OLW','Precip','Precip_acc', 'LH_release'] or 'Cloud' in var:
#    d.imp_atm(tn= t)
#
#if '_lev' in var:
#    d.imp_level(pn= pn, tn= t)
#
#if 'tropopause' in var:
#    d.imp_grib(tn= t, calcheight= True)
#
#if var in ['CTT', 'WVTC', 'WVT', 'CWR']:
#    if exp_name in ['08030312_cycling']: t//=3 #for some experiments the pseudo sat is only retrieved every 3 hour
#                   
#    d.imp_pseudo(tn= t)
#    if exp_name in ['08030312_cycling']: t*=3
##var= 'PressWind'
#
#
#
#if 'tropopause' in var:
#    Lapserate= 1E3* (d.T[1:] - d.T[:-1]) /d.dz     
#    zintermediate= (d.z[1:] + d.z[:-1])/2
#            
#    """works but extremly unefficient"""
##        tropheight, troppres= np.zeros(np.shape(Lapserate)[1:]), np.zeros(np.shape(Lapserate)[1:])
##        for xi in range(Lapserate.shape[1]):
##            print(xi)
##            for yi in range(Lapserate.shape[2]):
##                
##                tropheight[xi, yi]= next(zintermediate[i,xi, yi]
##                    for i,l in enumerate(Lapserate[:,xi, yi]) if  (l >-2 
##                        and np.min(Lapserate[i: np.argmin(np.abs(zintermediate[:, xi, yi] - zintermediate[i, xi, yi] - 2000)), xi, yi]) > -2 ) )
## 
#    """this is precises"""
#    zdist= 200
#    z_int= np.arange(3000, 13001, zdist) #this defines the range in which the tropopause height can be
#    from scipy.interpolate import interp1d        
#
#    Lapserate_int= np.moveaxis( np.array([[ interp1d(zintermediate[:, x, y], Lapserate[:, x, y])(z_int) for y in range(Lapserate.shape[2]) ] for x in range(Lapserate.shape[1]) ]  ), -1, 0)
#    p_int= np.moveaxis( np.array([[ interp1d(zintermediate[:, x, y], (d.pres[1:,x,y]+ d.pres[:-1,x,y])/2)(z_int) for y in range(Lapserate.shape[2]) ] for x in range(Lapserate.shape[1]) ]  ), -1, 0)
#
#    
#    levels2km= int(2000/zdist)
#    Lapserate_2kmabove= np.array([np.min(Lapserate_int[i: i+levels2km], axis= 0) for i in range(len(z_int) -levels2km)])
#    
#    z_int_total= z_int[:, np.newaxis, np.newaxis]* np.ones(np.shape(Lapserate_int) )
#    z_int_total[Lapserate_int < -2]= np.nan
#    z_int_total[:-levels2km][Lapserate_2kmabove < -2]= np.nan
#    tropheight= np.nanmin(z_int_total, axis= 0)
#
#
#if var == 'tropopause':
#    if display=='pres':
#        p_tropheight= np.array([[p_int[:, x,y][z_int== tropheight[x,y]][0] for y in range(Lapserate.shape[2]) ] for x in range(Lapserate.shape[1]) ]  )
#        
#        gausfilterdist= 100
#        p_tropheight= filters.gaussian_filter(p_tropheight, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)            
#        
#        bounds= np.arange(260, 361, 10)
#        PlotColorMap4(Lon, Lat, p_tropheight, map, bounds= bounds, color='viridis_r', label='Tropopause height [hPa]')
#        PlotContours(Lon, Lat, p_tropheight, map, leveldist= 10)
#
#    else: #display='z':
#        gausfilterdist= 100
#        tropheight= filters.gaussian_filter(tropheight, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#
#        bounds= np.arange(7000, 9001, 250)
#        PlotColorMap4(Lon, Lat, tropheight, map, bounds= bounds, label='Tropopause height [m]')
#        PlotContours(Lon, Lat, tropheight, map, leveldist= 250)
#
#
#if var== 'tropopause_advanced':
#    gausfilterdist= 100
#
#    tropheight= filters.gaussian_filter(tropheight, sigma= gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
#
#    bounds= np.arange(7000, 9251, 250)
#    PlotColorMap4(Lon, Lat, tropheight, map, bounds= bounds, label='Tropopause height [m]')#        PlotContours(Lon, Lat, p_tropheight, map, leveldist= 10)
#
#
#    d2= data(filename= AAfilename, res= AAres)
#    d2.imp_level(pn= 10, tn= t)
#
#    PlotContours(Lon, Lat, d2.geop, map, leveldist= 20, alpha= 1, numbers= True, color='black')
##    PlotWindVelo(Lon, Lat, np.sqrt(d.u**2+ d.v**2), map)
#
#    PlotWind(d2.lon, d2.lat, d2.u, d2.v, map, rot=False) #put this after PlotColorMap, otherwise it disappears
#
#    
##        plt.figure(fignr+1)
##        plt.clf()
##        plt.plot(Lapserate_int[:, 0, 0], p_int[:, 0, 0])
#
#
#
