#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

This is from plot_comparison
"""

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import numpy as np
import os
import time
#from datetime import date, datetime, timedelta
from scipy import stats

# import own modules
from f_imp_WRF import *  # Read in netcdf 

import sys  #to import the functions from a different directory
sys.path.insert(0, '/home/'+user+'/codeAROME/')
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)


"""global variables"""
#fignr= 2
#
#"""Lambert coordinates"""
#lllon, lllat, urlon, urlat= -15, 63, 60, 75
#lat0, lon0= 75, 0 #(lat0, lon0) = center point


domain= 1
year, month, day, hour= 2008, 3, 3, 0
t= 6 #which timestep of the data is regarded
l= 0 #level

#filename= (Mediadir+'PL/WRF/wrfout_d0'+str(domain)+'_'+str(year)
#+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+str(hour).zfill(2)+':00:00_INTRP' )
#
#nc= Dataset(filename)
#print(nc.variables.keys())
#u10= nc.variables['U10'][:]
#v10= nc.variables['V10'][:]
#geo= nc.variables['GHT'][:, l, :]
#T= nc.variables['TT'][:, l]

filenam= Mediadir+'PL/WRF/wrfout_d0'+str(domain)+'_'+str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+str(hour).zfill(2)+':00:00_INTRP'
d= WRFdata(filename= Mediadir+'PL/WRF/wrfout_d0'+str(domain)+'_'+str(year)
+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+str(hour).zfill(2)+':00:00_INTRP')
d.imp_standard()



plt.figure(fignr)
fignr+= 1
plt.clf()
#map= WRF_map(d.lat, d.lon, d.lat0, d.lon0)


map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)

Lon,Lat = map(d.lon,d.lat)


"""MSLP, Surface temperature and surface winds"""
#PlotSurfTemp(Lon,Lat, d.T[t, l]-273, map)
##"""calculate vorticity"""
#dx= 22500
#vort= (np.gradient(d.v[t,l], dx, axis= 1)- np.gradient(d.u[t,l], dx, axis= 0))* 1E5
#PlotVort(Lon, Lat, vort, map)

#PlotWind(d.lon, d.lat, d.u10[t], d.v10[t], map, rot= False)
PlotWindVelo(Lon, Lat, np.sqrt(d.u10[t]**2+ d.v10[t]**2), map, Umax= 25)

PlotContours(Lon, Lat, mslp(d.geo[t,l],d.T[t,l]), map, nrlevels= 20)


plt.title('WRF ('+str(year)+', '+str(month)+', '+str(int(day+(hour+d.tim[t])//24))+', '+str(int((hour+d.tim[t])%24))+')')


plt.figure(fignr)
fignr+= 1
plt.clf()

map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)

pn= 1 #1 - 800hPa
theta= d.T[t, pn]*(1000/d.lev[pn])**(2/7)
#spechum= d.relhum/(.263 *d.press[pn])* np.exp(17.67*(d.Temp-273.15) / (d.Temp-29.65))
#theta_e= theta* np.exp(2501 * spechum/(1.006* d.Temp))

PlotColorMap(Lon, Lat, theta, map,  variable='theta 800')
PlotContours(Lon, Lat, mslp(d.geo[t,l],d.T[t,l]), map, nrlevels= 20)


theta500= PotTemp(d.T[t, 2], d.lev[2])
thetaSST= PotTemp(d.SST[t], d.psfc[t])
PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')
plt.title('WRF ('+str(year)+', '+str(month)+', '+str(int(day+(hour+d.tim[t])//24))+', '+str(int((hour+d.tim[t])%24))+')')


# 1000 -500 hPa thickness
#pn= 8
#d.imp_level(pn= pn, tn= t)
##
#geop_diff= d.geop - geop0
#PlotContours(Lon, Lat, geop_diff, m)
#
#""" SST- T500 """
#PlotStaticStability(Lon, Lat, Temp0 - d.Temp, m)


#"""domain 2"""
#t= t*3
#domain= 2
#
#d= WRFdata(filename= Mediadir+'PL/WRF/wrfout_d0'+str(domain)+'_'+str(year)
#+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+str(hour).zfill(2)+':00:00_INTRP')
#d.imp_standard()
#
##"""calculate vorticity"""
#dx= 7500
#vort= (np.gradient(d.v[t,l], dx, axis= 1)- np.gradient(d.u[t,l], dx, axis= 0)) *1E5
##
#plt.figure(3)
#plt.clf()
#map= WRF_map(d.lat, d.lon, d.lat0, d.lon0)
#Lon,Lat = map(d.lon,d.lat)
#
#
#"""MSLP, Surface temperature and surface winds"""
##PlotSurfTemp(Lon,Lat, d.T[t, l]-273, map)
#PlotVort(Lon, Lat, vort, map)
#
#PlotWind(d.lon, d.lat, d.u10[t], d.v10[t], map, rot= False)
#
#PlotContours(Lon, Lat, mslp(d.geo[t,l],d.T[t,l]), map)


#"""domain 3"""
#domain= 3
#
#d= WRFdata(filename= Mediadir+'PL/WRF/wrfout_d0'+str(domain)+'_'+str(year)
#+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+str(hour).zfill(2)+':00:00_INTRP')
#d.imp_standard()
##
#dx= 2500
#vort= (np.gradient(d.v[t,l], dx, axis= 1)- np.gradient(d.u[t,l], dx, axis= 0)) *1E5
#
#plt.figure(4)
#plt.clf()
#map= WRF_map(d.lat, d.lon, d.lat0, d.lon0)
#Lon,Lat = map(d.lon,d.lat)
#
#
#"""MSLP, Surface temperature and surface winds"""
##PlotSurfTemp(Lon,Lat, d.T[t, l]-273, map)
#PlotVort(Lon, Lat, vort, map)
#
#PlotWind(d.lon, d.lat, d.u10[t], d.v10[t], map, rot= False)
#
#PlotContours(Lon, Lat, mslp(d.geo[t,l],d.T[t,l]), map)

#
#
#x_start, y_start= 10, 30
#x_end, y_end= 110, 20
##
#x, y= mapcross_section(d.lon, d.lat, x_start, y_start, x_end, y_end, m)
##
###d.imp_level(pn= 10, tn= t)
###PlotContours(Lon, Lat, d.mslp, m)
###PlotPV(Lon, Lat, d.PV, m)
#
#
#
#"""equivalent temperature and geopotential height difference"""
##plt.figure(2)
##plt.clf()
##
##m = AA_map()
##
##plt.title('1000-500hPa thickness and static stability')
##
##pn= 4
##d.imp_level(pn= pn, tn= t)
##
##theta= d.Temp*(1000/d.press[pn])**(2/7)
##spechum= d.relhum/(.263 *d.press[pn])* np.exp(17.67*(d.Temp-273.15) / (d.Temp-29.65))
##theta_e= theta* np.exp(2501 * spechum/(1.006* d.Temp))
##
##PlotColorMap(Lon, Lat, theta_e, m,  variable='850 hPA theta_e')
##
##
### 1000 -500 hPa thickness
##pn= 8
##d.imp_level(pn= pn, tn= t)
###
##geop_diff= d.geop - geop0
##PlotContours(Lon, Lat, geop_diff, m)
#
#""" SST- T500 """
##PlotStaticStability(Lon, Lat, Temp0 - d.Temp, m)
#
#""" 700 hPa"""
##fignr= 3
##plt.figure(fignr)
##plt.clf()
##fignr += 1
##
##m = AA_map()
##d.imp_level(pn= 6, tn= t)
##PlotContours(Lon, Lat, d.geop, m)
##
##PlotPrecip(Lon, Lat, d.prec, m)
##PlotWind(d.lon, d.lat, d.u, d.v, m)
#
#
#
#"""cross section temperature"""
##plt.figure(fignr)
##plt.clf()
##fignr += 1
##
##plt.title('temperature')
##
#d.imp_cross_sec(xn= x, yn= y, tn= t)
##PlotCross_sec_T(d.press[:-5], d.lon[y,x], d.Temp[:-5])
#
#
#"""cross section potential temperature"""
##plt.figure(fignr)
##plt.clf()
##fignr += 1
##
##plt.title('potential temperature')
##
#pressm= np.tile(d.press, (np.shape(d.Temp)[1],1)).T #make a matrix with the pressure in every row
#theta= d.Temp*(1000/pressm)**(2/7)
##
##PlotCross_sec_T(d.press[:-5], d.lon[y,x], theta[:-5])
##
#fignr=4
#""" cross section eqivalent potential temperature"""
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
#
#fignr= 5
#"""vertical velocity"""
##plt.figure(fignr)
##plt.clf()
##fignr += 1
##
##PlotCross_sec_vervel(d.press, d.lon[y,x], d.w)
#
#
#fignr= 6
#"""horizontal velocity (north- south)"""
##plt.figure(fignr)
##plt.clf()
##fignr += 1
##
##PCross_sec_horvel(d.press, d.lon[y,x], d.u, d.v)
##
##
##plt.figure(fignr)
##plt.clf()
##fignr += 1
##
##PCross_sec_horvel2(d.press, d.lon[y,x], d.u, d.v)
#
#"""heat flux"""
#fignr= 7
#plt.figure(fignr)
#plt.clf()
#fignr += 1
#
#m = AA_map()
#d.imp_level(pn= 0, tn= t)
#
#PlotHeatFlux(Lon, Lat, d.LH, m, label= 'Latent')
#
#"""heat flux"""
#fignr= 8
#plt.figure(fignr)
#plt.clf()
#fignr += 1
#
#m = AA_map()
#
#PlotHeatFlux(Lon, Lat, d.SH, m, label= 'Sensible')