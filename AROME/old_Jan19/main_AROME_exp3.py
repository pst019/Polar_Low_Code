#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

This is from plot_comparison
"""

import sys  #to import the functions from a different directory
sys.path.insert(0, '/home/'+user+'/codeAROME/')



from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import numpy as np
import os
import time
#from datetime import date, datetime, timedelta
from scipy import stats

# import own modules
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_imp_AROME_exp import *  # Read in netcdf file

t= 17 


d= data(filename= Mediadir+'PL/ec/080303_test_3_2008030300_fp_extract.nc')
#d.imp_standard()



"""MSLP, Surface temperature and surface winds"""
plt.figure(13)
plt.clf()

plt.title(str(time.localtime(d.tim[t])[:4]))
##plt.subplot(2, 2, 1)
#
m = AA_map()
Lon,Lat = m(d.lon,d.lat)

d.imp_level(pn= 0, tn= t)
geop0 = d.geop
Temp0 = d.Temp


"""calculate vorticity"""
#dx= 2500
#vort= (np.gradient(d.v, dx, axis= 1)- np.gradient(d.u, dx, axis= 0))*1E5
##probably one has to rotate u and v first to proper x and y coordinates
#
#PlotVort(Lon, Lat, vort, m)
#
## Draw contours for mslp
#PlotContours(Lon, Lat, d.mslp, m)
##
##PlotSurfTemp(Lon,Lat, d.Temp - 273.15, m)
##
##
###Umax= np.max([d.U10[t,1], d2.U10[t,1]])
###PlotWindVelo(Lon, Lat, d.U[t,-1,:,:], m)
#PlotWind(d.lon, d.lat, d.u, d.v, m)
#
x_start, y_start= 10, 30
x_end, y_end= 110, 20
##
x, y= mapcross_section(d.lon, d.lat, x_start, y_start, x_end, y_end, m)
##
###d.imp_level(pn= 10, tn= t)
###PlotContours(Lon, Lat, d.mslp, m)
###PlotPV(Lon, Lat, d.PV, m)





"""equivalent temperature and geopotential height difference"""
#plt.figure(2)
#plt.clf()
#
#m = AA_map()
#
#plt.title('1000-500hPa thickness and static stability')
#
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

""" SST- T500 """
#PlotStaticStability(Lon, Lat, Temp0 - d.Temp, m)

""" 700 hPa"""
fignr= 3
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



"""cross section temperature"""
plt.figure(fignr)
plt.clf()
fignr += 1

plt.title('temperature')

d.imp_cross_sec(xn= x, yn= y, tn= t)
PlotCross_sec_T(d.press[:-5], d.lon[y,x], d.Temp[:-5])


"""cross section potential temperature"""
plt.figure(fignr)
plt.clf()
fignr += 1

plt.title('potential temperature')

pressm= np.tile(d.press, (np.shape(d.Temp)[1],1)).T #make a matrix with the pressure in every row
theta= d.Temp*(1000/pressm)**(2/7)

PlotCross_sec_T(d.press[:-5], d.lon[y,x], theta[:-5])

#fignr=4
""" cross section eqivalent potential temperature"""
plt.figure(fignr)
plt.clf()
fignr += 1

plt.title('equivalent potential temperature')

#rel humidity to specific humidity: http://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
#spechum= relhum[0to1]*100/(26.3 * press[hPa]) * exp(17.67*(T[k]-273.15)/(T[K]-29.65))
spechum= d.relhum/(.263 *pressm)* np.exp(17.67*(d.Temp-273.15) / (d.Temp-29.65))

#equivalent pot temp theta_e = theta* exp(Lc * spechum/(Cp * T))
theta_e= theta* np.exp(2501 * spechum/(1.006* d.Temp))

PlotCross_sec_T(d.press[:-5], d.lon[y,x], theta_e[:-5])

PlotCross_sec_hum(d.press[:-5], d.lon[y,x], d.relhum[:-5])

fignr= 5
"""vertical velocity"""
#plt.figure(fignr)
#plt.clf()
#fignr += 1
#
#PlotCross_sec_vervel(d.press, d.lon[y,x], d.w)


fignr= 6
"""horizontal velocity (north- south)"""
#plt.figure(fignr)
#plt.clf()
#fignr += 1
#
#PCross_sec_horvel(d.press, d.lon[y,x], d.u, d.v)
#
#
#plt.figure(fignr)
#plt.clf()
#fignr += 1
#
#PCross_sec_horvel2(d.press, d.lon[y,x], d.u, d.v)

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