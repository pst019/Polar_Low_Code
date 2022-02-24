#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

This is from plot_comparison
"""

import os
user = os.getcwd().split('/')[2]

Mediadir= '/media/'+user+'/1692A00D929FEF8B/'


import sys  #to import the functions from a different directory

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import numpy as np
import os
import time
#from datetime import date, datetime, timedelta
from scipy import stats



sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')
from f_meteo import *
# import own modules
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_imp_AROME_exp import *  # Read in netcdf file

"""global variables"""
#fignr= 13
#"""Lambert coordinates"""
#lllon, lllat, urlon, urlat= -14, 64, 50, 75
#lat0, lon0= 75, 0 #(lat0, lon0) = center point

                   
t= 11 


d= data(filename= Mediadir+'PL/ec/080303_test_3_2008030300_fp_extract.nc')
#d.imp_standard()



"""MSLP, Surface temperature and surface winds"""
plt.figure(fignr)
fignr+=1
plt.clf()

#
#map = AA_map()

map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)

Lon,Lat = map(d.lon,d.lat)

d.imp_level(pn= 0, tn= t)
geop1000 = d.geop
Temp1000 = d.Temp


"""calculate vorticity"""
#dx= 2500
#vort= (np.gradient(d.v, dx, axis= 1)- np.gradient(d.u, dx, axis= 0))*1E5
##probably one has to rotate u and v first to proper x and y coordinates
#
#PlotVort(Lon, Lat, vort, m)

# Draw contours for mslp
PlotContours(Lon, Lat, d.mslp, map)
#
#PlotSurfTemp(Lon,Lat, d.Temp1000 - 273.15, m)

#PlotWind(d.lon, d.lat, d.u, d.v, map)
PlotWindVelo(Lon, Lat, np.sqrt(d.u**2+ d.v**2), map, Umax= 25)

plt.title('Arome '+str(time.localtime(d.tim[t])[:4]))




"""for cross section"""
#x_start, y_start= 10, 30
#x_end, y_end= 110, 20
#
#x, y= mapcross_section(d.lon, d.lat, x_start, y_start, x_end, y_end, m)



#"""PV"""
##plt.figure(fignr)
##fignr+=1
##plt.clf()
##
##plt.title('PV')
###
##m = AA_map()
##d.imp_level(pn= 10, tn= t)
##PlotContours(Lon, Lat, d.mslp, m)
##PlotPV(Lon, Lat, d.PV, m)
#
#
#
#
#
"""potential temperature at 700 hPa and stability"""
plt.figure(fignr)
fignr+=1
plt.clf()

map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)


PlotContours(Lon, Lat, d.mslp, map, nrlevels= 20)

pn= 6 #6= 700hPa
d.imp_level(pn= pn, tn= t)

theta= PotTemp(d.Temp, d.press[pn])
PlotColorMap(Lon, Lat, theta, map,  variable='theta 700')

plt.title('Arome '+str(time.localtime(d.tim[t])[:4]) +'theta 700 and static stability')


#"""stability"""
pn= 8 # 8 = 500hPa
d.imp_level(pn= pn, tn= t)

theta500= PotTemp(d.Temp, d.press[pn])
thetaSST= Temp1000
PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')


#PlotStaticStability(Lon, Lat, Temp0 - d.Temp, m)

#plt.title('1000-500hPa thickness and static stability')
"""equivalent temperature and geopotential height difference"""

#pn= 4
#d.imp_level(pn= pn, tn= t)
#
#theta= d.Temp*(1000/d.press[pn])**(2/7)
#spechum= d.relhum/(.263 *d.press[pn])* np.exp(17.67*(d.Temp-273.15) / (d.Temp-29.65))
#theta_e= theta* np.exp(2501 * spechum/(1.006* d.Temp))
#
#PlotColorMap(Lon, Lat, theta_e, m,  variable='850 hPA theta_e')
#
#
## 1000 -500 hPa thickness
#pn= 8
#d.imp_level(pn= pn, tn= t)
##
#geop_diff= d.geop - geop0
#PlotContours(Lon, Lat, geop_diff, m)
#
#""" SST- T500 """
#PlotStaticStability(Lon, Lat, Temp0 - d.Temp, m)
#
#""" 700 hPa"""
#plt.figure(fignr)
#plt.clf()
#fignr += 1
#
#m = AA_map()
#d.imp_level(pn= 6, tn= t)
#PlotContours(Lon, Lat, d.geop, m)
#
#PlotPrecip(Lon, Lat, d.prec, m)
#PlotWind(d.lon, d.lat, d.u, d.v, m)
#
#
#
##"""cross section temperature"""
##plt.figure(fignr)
##fignr+=1
##plt.clf()
##
##plt.title('temperature')
##
##d.imp_cross_sec(xn= x, yn= y, tn= t)
##PlotCross_sec_T(d.press[:-5], d.lon[y,x], d.Temp[:-5])
##
##
##"""cross section potential temperature"""
##plt.figure(fignr)
##plt.clf()
##fignr += 1
##
##plt.title('potential temperature')
##
##pressm= np.tile(d.press, (np.shape(d.Temp)[1],1)).T #make a matrix with the pressure in every row
##theta= d.Temp*(1000/pressm)**(2/7)
##
##PlotCross_sec_T(d.press[:-5], d.lon[y,x], theta[:-5])
##
###fignr=4
##""" cross section eqivalent potential temperature"""
##plt.figure(fignr)
##plt.clf()
##fignr += 1
##
##plt.title('equivalent potential temperature')
##
###rel humidity to specific humidity: http://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
###spechum= relhum[0to1]*100/(26.3 * press[hPa]) * exp(17.67*(T[k]-273.15)/(T[K]-29.65))
##spechum= d.relhum/(.263 *pressm)* np.exp(17.67*(d.Temp-273.15) / (d.Temp-29.65))
##
###equivalent pot temp theta_e = theta* exp(Lc * spechum/(Cp * T))
##theta_e= theta* np.exp(2501 * spechum/(1.006* d.Temp))
##
##PlotCross_sec_T(d.press[:-5], d.lon[y,x], theta_e[:-5])
##
##PlotCross_sec_hum(d.press[:-5], d.lon[y,x], d.relhum[:-5])
##
##fignr= 5
##"""vertical velocity"""
###plt.figure(fignr)
###plt.clf()
###fignr += 1
###
###PlotCross_sec_vervel(d.press, d.lon[y,x], d.w)
##
##
##fignr= 6
##"""horizontal velocity (north- south)"""
###plt.figure(fignr)
###plt.clf()
###fignr += 1
###
###PCross_sec_horvel(d.press, d.lon[y,x], d.u, d.v)
###
###
###plt.figure(fignr)
###plt.clf()
###fignr += 1
###
###PCross_sec_horvel2(d.press, d.lon[y,x], d.u, d.v)
##
###"""heat flux"""
###fignr= 7
###plt.figure(fignr)
###plt.clf()
###fignr += 1
###
###m = AA_map()
###d.imp_level(pn= 0, tn= t)
###
###PlotHeatFlux(Lon, Lat, d.LH, m, label= 'Latent')
###
###"""heat flux"""
###fignr= 8
###plt.figure(fignr)
###plt.clf()
###fignr += 1
###
###m = AA_map()
###
###PlotHeatFlux(Lon, Lat, d.SH, m, label= 'Sensible')