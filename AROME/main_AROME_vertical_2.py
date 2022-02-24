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

import numpy as np
import time
#from datetime import date, datetime, timedelta
from scipy import stats

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
fignr= 3

plt.figure(fignr)
fignr+=1
plt.clf()
                  
                   
"""b) times"""
year, month = 2008, 3
day, hour= 3, 12                   


"""Lambert coordinates"""
#lllon, lllat, urlon, urlat= -15, 65, 50, 75
#lat0, lon0= 75, 0 #(lat0, lon0) = center point

#lllon, lllat, urlon, urlat= 0, 67, 60, 75
#lat0, lon0= 70, 0 #(lat0, lon0) = center point

#to match with Thorpex                    
if year == 2008 and month== 3 and day== 3 and hour<15:
    lllon, lllat, urlon, urlat= -4, 70, 19, 75
    lat0, lon0= 75, 0 #(lat0, lon0) = center point

elif year == 2008 and month== 3 and day== 3 and hour>=15:
    lllon, lllat, urlon, urlat= -5, 69.5, 10, 74.5
    lat0, lon0= 75, 0 #(lat0, lon0) = center point

elif year == 2008 and month== 3 and day== 4:
    lllon, lllat, urlon, urlat= -1, 64, 12, 70
    lat0, lon0= 70, 0 #(lat0, lon0) = center point    
    

map = AA_map()
#map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)


fileday, filehour= 3, 0        
t= (day- fileday)*24 + (hour- filehour)  # -1 

#year, month = 1987, 2
#day, hour= 27, 12    

#fileday, filehour= 25, 12        
#t= (day- fileday)*24 + (hour- filehour)  # -1 

 
    
"""name of the experiment"""
exp_name= '080303_warmctr'
#exp_name='080303_cold_sice'
#exp_name= '080303_warmsens_noQH'
#exp_name='080303_warmsens_2FLX'
#exp_name='080303_warmsens_nocondens'
#exp_name='870226_cold_ERA_can'
exp_name='080303_cold_pseudo2'
                   
"""c) plot characteristics"""
pleveldist= 1                   

var= 'PressWind'
#var= 'PressWind_advanced'


""" presure level (if the variable ends with '_var' """
pn = 6
#4- 850hPa, 6- 700hPa, 8 - 500hPa

save=False
#savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/DropComp/'

title_extra=''


Thorpex= False

"""end global variables"""



"""start module"""
AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
AAres=1 #2- every second datapoint is taken
d= data(filename= AAfilename, res= AAres)

Lon,Lat = map(d.lon,d.lat)
    
    
d.imp_surf(tn= t)

if var in ['FLX', 'LH', 'TH', 'BL', 'CAPE', 'CIN', 'LCL', 'LFC'] or 'Cloud' in var:
    d.imp_atm(tn= t)

if '_lev' in var:
    d.imp_level(pn= pn, tn= t)



#else:
#    d.imp_pseudo(tn= t)




"""MSLP, Surface temperature and surface winds"""
if var== 'PressWind':
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    #PlotSurfTemp(Lon,Lat, d.T1000 - 273.15, m)
    PlotWind(d.lon, d.lat, d.u10m, d.v10m, map) #???
    PlotWindVelo(Lon, Lat, np.sqrt(d.u10m**2+ d.v10m**2), map, Umax= 25)
    

"""MSLP, Surface temperature and surface winds -advanced"""
if var== 'PressWind_advanced':
    PlotContours(Lon, Lat, d.mslp, map, leveldist= 1, numbers=False)

    PlotWindVelo(Lon, Lat, np.sqrt(d.u**2+ d.v**2), map, Umax= 25)    
    #find local pressure min and plot them in the map   
    U= np.sqrt(d.u10m**2+d.v10m**2)
    PlotLocalMax(d.mslp, threshold=1010, distance=100/AAres, map= map, lon=d.lon, lat=d.lat,
                 data2=U, threshold2=16, distance2=80/AAres)    
  
    PlotLocalMax(U, threshold=20, distance=100/AAres, map= map, lon=d.lon, lat=d.lat, typ='max', color='orange', dot=False)    
   


"""for cross section"""
x_start= thor.lon[startdrop,-20]
x_end= thor.lon[enddrop, -20]
y_start=thor.lat[startdrop,-20]
y_end= thor.lat[enddrop, -20]


#    x_start, y_start= 70, 30
#    x_end, y_end= 90, 60

x, y= mapcross_section(d.lon, d.lat, x_start, y_start, x_end, y_end, map, coordinates='latlon')


d.imp_cross_sec(xn= x, yn= y, tn= t)





"""cross section temperature"""
#plt.figure(fignr)
#fignr+=1
#plt.clf()
#
#plt.title('temperature')
#PlotCross_sec_T(d.pres[:-5], d.lon[y,x], d.T[:-5])


#"""cross section potential temperature"""
plt.figure(fignr)
plt.clf()
fignr += 1

plt.title('Potential temperature')

pressm= np.tile(d.pres, (np.shape(d.T)[1],1)).T #make a matrix with the pressure in every row
theta= d.T*(1000/pressm)**(2/7)

PlotCross_sec_T(d.pres, d.lon[x, y], theta, levels= np.arange(260, 290))#, color='default')
plt.ylim([1000, 450])


"""make a profile"""
plt.figure(fignr)
plt.clf()
fignr += 1

plt.plot(theta[:, 1], d.pres)
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
#PlotCross_sec_T(d.pres[:-5], d.lon[y,x], theta_e[:-5])
#
#PlotCross_sec_hum(d.pres[:-5], d.lon[y,x], d.RH[:-5])

#fignr= 5
"""vertical velocity"""
#plt.figure(fignr)
#plt.clf()
#fignr += 1
#
#PlotCross_sec_vervel(d.pres, d.lon[y,x], d.w)


#fignr= 6
"""horizontal velocity (north- south)"""
#plt.figure(fignr)
#plt.clf()
#fignr += 1
#
#PCross_sec_horvel(d.pres, d.lon[y,x], d.u, d.v)
#
#
#plt.figure(fignr)
#plt.clf()
#fignr += 1
#
#PCross_sec_horvel2(d.pres, d.lon[y,x], d.u, d.v)



"""Plot Thorpex points"""
if Thorpex== True:
    thor= thorpex.data(year, month, day, hour, plevels= [d.chosen_plevel])
#    thor2= thorpex.data(year, month, day, hour, level= d.chosen_plevel)
    
    xpt, ypt= map(thor.lon, thor.lat)
    #map.plot(xpt,ypt,'bo')
    
    for i in np.arange(len(xpt)):
#        plt.text(xpt[i]+100, ypt[i]+100, str(thor.dropnr[i])+'\n'+ str(thor.datetime[i, 0])[11:16], color='r')
                
        
        if var== 'Theta_lev':
            value= thor.theta[i]
                
            colors= [(plt.cm.RdBu_r(h)) for h in range(0, 256)]
            colornow= colors[ int((value-bounds[0])/ (bounds[-1] - bounds[0]) * 256)]

        elif var== 'RH_lev':
            value= thor.RH[i]
            
            a= value - bounds
            ind= np.where(a== min(i for i in a if i >= 0))[0][0]+1
            colornow= colors[ind]

        elif var== 'Geop_wind_lev':
            dist= np.sqrt( (thor.lat[i,0]- d.lat)**2+ (np.cos(np.deg2rad(thor.lat[i,0]))* (thor.lon[i,0]- d.lon))**2)
            xpos, ypos= np.where(dist== np.min(dist))
            plt.text(xpt[i], ypt[i]+10000, str(int(thor.alt[i,0]- d.geop[xpos[0], ypos[0]])), color='r')  #depature from the geopotential height

            value= thor.U[i]
            
            a= value - bounds
            ind= np.where(a== min(i for i in a if i >= 0))[0][0]+1
            colornow= colors[ind]
#            print(thor.dropnr[i], thor.u[i], thor.v[i])

            
        else:
            break

#        plt.text(xpt[i]+100, ypt[i]+100, str(np.round(value[0],2) ) )          
        map.scatter(xpt[i],ypt[i], color= colornow, s= 150, edgecolors= 'r')


        if var== 'Geop_wind_lev':
            u,v, Lon, Lat = map.rotate_vector(thor.u[i], thor.v[i] , thor.lon[i], thor.lat[i] ,returnxy=True)
            plt.quiver(Lon, Lat, u, v, color= 'r', scale = 500, width= .004)


if save== False:
    plt.title('Arome '+str(d.datetime[t])[:-3]+ title_extra)

plt.tight_layout()

if save==True:
    if '_lev' in var:
        var += str(d.chosen_plevel)
    plt.savefig(savedir+'Arome_'+exp_name+'_'+str(t).zfill(2)+'_'+ var +title_extra)








    


    

"""for cross section definition"""
#x_start, y_start= 10, 30
#x_end, y_end= 110, 20
#
#x, y= mapcross_section(d.lon, d.lat, x_start, y_start, x_end, y_end, m)
