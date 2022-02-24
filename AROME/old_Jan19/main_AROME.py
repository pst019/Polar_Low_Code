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

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+'/Polar_Low/Code2/AROME/')

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import numpy as np
import time
#from datetime import date, datetime, timedelta
from scipy import stats

# import own modules
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_imp_AROME_exp import *  # Read in netcdf file

sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')

from f_meteo import *


"""global variables"""
fignr= 5

"""a) Lambert coordinates"""
lllon, lllat, urlon, urlat= -15, 65, 50, 75
lat0, lon0= 75, 0 #(lat0, lon0) = center point

lllon, lllat, urlon, urlat= -5, 67, 28, 76
lat0, lon0= 70, 0 #(lat0, lon0) = center point

                   
"""b) times"""
year, month = 2008, 3
day, hour= 3, 12                    
                    
"""c) plot characteristics"""
PlotPressWind= True
#number_plevels= 20
pleveldist= 1                   

PlotTheta700= True
pot_temp_bounds= np.arange(264, 290, 1)

PlotTheta_e850= True

PlotVorticity=False

save=True
savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/movies/'

"""end global variables"""


"""start module"""
#fileday, filehour= 3, 12        
#t= (day- fileday)*24 + (hour- filehour) -1    
#d= data(filename= Mediadir+'PL/AA/ec/'+str(year)[-2:]+str(month).zfill(2)+str(fileday).zfill(2)+'_test_3_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc')
#d.imp_standard()

exp_name= '080303_warmctr'
#exp_name= '080303_warmsens_noTH'
fileday, filehour= 3, 0        
t= (day- fileday)*24 + (hour- filehour)  # -1 
AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
d= data(filename= AAfilename)


"""MSLP, Surface temperature and surface winds"""
if PlotPressWind== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
#    map = AA_map()
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
    Lon,Lat = map(d.lon,d.lat)
    
    d.imp_level(pn= 0, tn= t)
    geop1000 = d.geop
    Temp1000 = d.Temp
    
    # Draw contours for mslp
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    #
    #PlotSurfTemp(Lon,Lat, d.Temp1000 - 273.15, m)
    #PlotWind(d.lon, d.lat, d.u, d.v, map)
    PlotWindVelo(Lon, Lat, np.sqrt(d.u**2+ d.v**2), map, Umax= 25)
    
    plt.title('Arome '+str(time.localtime(d.tim[t]-1)[:4]))
    plt.tight_layout()
    
    if save==True:
        plt.savefig(savedir+'Arome_'+str(t) )

    """for cross section"""
    x_start= thor.lon[startdrop,-20]
    x_end= thor.lon[enddrop, -20]
    y_start=thor.lat[startdrop,-20]
    y_end= thor.lat[enddrop, -20]
    
    
#    x_start, y_start= 70, 30
#    x_end, y_end= 90, 60
    
    x, y= mapcross_section(d.lon, d.lat, x_start, y_start, x_end, y_end, map, coordinates='latlon')


"""potential temperature at 700 hPa and stability"""
#if PlotTheta700== True:
#    plt.figure(fignr)
#    fignr+=1
#    plt.clf()
#    
#    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
#    
#    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
#    
#    pn= 6 #6= 700hPa
#    d.imp_level(pn= pn, tn= t)
#    
#    theta= PotTemp(d.Temp, d.press[pn])
#    PlotColorMap3(Lon, Lat, theta, map, bounds= pot_temp_bounds, label=r"theta 700")
#    
#    plt.title('Arome '+str(time.localtime(d.tim[t])[:4]) +'theta 700 and static stability')
#    
#    
#    #"""stability"""
#    pn= 8 # 8 = 500hPa
#    d.imp_level(pn= pn, tn= t)
#    
#    theta500= PotTemp(d.Temp, d.press[pn])
#    thetaSST= Temp1000
#    PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')


""" equivalent potential temperature at 850hPa"""
if PlotTheta_e850 == True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
    
    pn= 4 #850hPa
    d.imp_level(pn= pn, tn= t)
    
#    theta= d.Temp*(1000/d.press[pn])**(2/7)
    spechum= d.relhum/(.263 *d.press[pn])* np.exp(17.67*(d.Temp-273.15) / (d.Temp-29.65))
#    theta_e= theta* np.exp(2501 * spechum/(1.006* d.Temp))
    theta_e= EquiPotTemp(d.Temp, spechum, d.press[pn])
    
    PlotColorMap3(Lon, Lat, theta_e, map, bounds= pot_temp_bounds, label=r"$\theta_{e,850}$ [K]")
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    plt.title('Arome '+str(time.localtime(d.tim[t])[:4]))

"""calculate vorticity"""
#dx= 2500
#vort= (np.gradient(d.v, dx, axis= 1)- np.gradient(d.u, dx, axis= 0))*1E5
##probably one has to rotate u and v first to proper x and y coordinates
#
#PlotVort(Lon, Lat, vort, m)


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
d.imp_cross_sec(xn= x, yn= y, tn= t)


"""cross section temperature"""
#plt.figure(fignr)
#fignr+=1
#plt.clf()
#
#plt.title('temperature')
#PlotCross_sec_T(d.press[:-5], d.lon[y,x], d.Temp[:-5])


#"""cross section potential temperature"""
plt.figure(fignr)
plt.clf()
fignr += 1

plt.title('Potential temperature')

pressm= np.tile(d.press, (np.shape(d.Temp)[1],1)).T #make a matrix with the pressure in every row
theta= d.Temp*(1000/pressm)**(2/7)

PlotCross_sec_T(d.press, d.lon[x, y], theta, levels= np.arange(260, 290))#, color='default')
plt.ylim([1000, 450])


#
##fignr=4
#""" cross section eqivalent potential temperature"""
#plt.figure(fignr)
#plt.clf()
#fignr += 1
#
#plt.title('equivalent potential temperature')
#
##rel humidity to specific humidity: http://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
##spechum= relhum[0to1]*100/(26.3 * press[hPa]) * exp(17.67*(T[k]-273.15)/(T[K]-29.65))
#spechum= d.relhum/(.263 *pressm)* np.exp(17.67*(d.Temp-273.15) / (d.Temp-29.65))
#
##equivalent pot temp theta_e = theta* exp(Lc * spechum/(Cp * T))
#theta_e= theta* np.exp(2501 * spechum/(1.006* d.Temp))
#
#PlotCross_sec_T(d.press[:-5], d.lon[y,x], theta_e[:-5])
#
#PlotCross_sec_hum(d.press[:-5], d.lon[y,x], d.relhum[:-5])

#fignr= 5
"""vertical velocity"""
#plt.figure(fignr)
#plt.clf()
#fignr += 1
#
#PlotCross_sec_vervel(d.press, d.lon[y,x], d.w)


#fignr= 6
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