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
from f_imp_AROME import *  # Read in netcdf file

sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')

from f_meteo import *


"""global variables"""
fignr= 1

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
PlotTheta_e850= False
PlotVorticity=False

save=False
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
    
    d.imp_surf(tn= t)
    d.imp_level(pn= 0, tn= t) #this is not really needed here.
    geop1000 = d.geop
    Temp1000 = d.T
    
    # Draw contours for mslp
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    #
    #PlotSurfTemp(Lon,Lat, d.T1000 - 273.15, m)
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
#PlotCross_sec_T(d.press[:-5], d.lon[y,x], d.T[:-5])


#"""cross section potential temperature"""
plt.figure(fignr)
plt.clf()
fignr += 1

plt.title('Potential temperature')

pressm= np.tile(d.press, (np.shape(d.T)[1],1)).T #make a matrix with the pressure in every row
theta= d.T*(1000/pressm)**(2/7)

PlotCross_sec_T(d.press, d.lon[x, y], theta, levels= np.arange(260, 290))#, color='default')
plt.ylim([1000, 450])


"""make a profile"""
plt.figure(fignr)
plt.clf()
fignr += 1

plt.plot(theta[:, 1], d.press)
plt.gca().invert_yaxis()

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
#spechum= d.RH/(.263 *pressm)* np.exp(17.67*(d.T-273.15) / (d.T-29.65))
#
##equivalent pot temp theta_e = theta* exp(Lc * spechum/(Cp * T))
#theta_e= theta* np.exp(2501 * spechum/(1.006* d.T))
#
#PlotCross_sec_T(d.press[:-5], d.lon[y,x], theta_e[:-5])
#
#PlotCross_sec_hum(d.press[:-5], d.lon[y,x], d.RH[:-5])

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