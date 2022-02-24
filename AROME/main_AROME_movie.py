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
from f_meteo import *

"""global variables"""
fignr= 4

"""a) Lambert coordinates"""
lllon, lllat, urlon, urlat= -15, 65, 50, 75
lat0, lon0= 75, 0 #(lat0, lon0) = center point

#lllon, lllat, urlon, urlat= -5, 67, 28, 76
#lat0, lon0= 70, 0 #(lat0, lon0) = center point

                   
"""b) times"""
#for one specific point in time:
year, month = 2008, 3
day, hour= 4, 12                    
#
fileday, filehour= 3, 0        

#year, month = 1987, 2
#day, hour= 27, 11    

#fileday, filehour= 25, 12        

#figure for the whole range, to make a movie:
movie=True
                    
"""c) plot characteristics"""
PlotPressWind= True
#number_plevels= 20
pleveldist= 1                   

PlotTheta700= False
pot_temp_bounds= np.arange(264, 290, 1)

PlotTheta_e850= False

PlotVorticity=False

save=True
savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/movies2/' #this replaces /movies/ which could be deleted

"""end global variables"""


"""start module"""
#fileday, filehour= 3, 12        
#t= (day- fileday)*24 + (hour- filehour) -1    
#d= data(filename= Mediadir+'PL/AA/ec/'+str(year)[-2:]+str(month).zfill(2)+str(fileday).zfill(2)+'_test_3_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc')
#d.imp_standard()

#exp_name= '080303_warmctr'
#exp_name= '080303_warmsens_noTH'
exp_name= '080303_warmsens_noFLX'
#exp_name= '080303_coldsens_noFLX_AREA'

#exp_name= '080303_cold'
#exp_name='870226_cold_ERA_can'


AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
d= data(filename= AAfilename, res=1)


if movie==False:
    t= np.array((day- fileday)*24 + (hour- filehour))  # -1 
else:
    t= np.arange(len(d.tim))

"""MSLP, Surface temperature and surface winds"""
if PlotPressWind== True:

    for ti in t:
        print(ti)
        plt.figure(fignr)
        plt.clf()
        
    #    map = AA_map()
        map= AA_map_half()
#        map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
        Lon,Lat = map(d.lon,d.lat)

        d.imp_surf(tn=ti)         
#        d.imp_level(pn= 0, tn= ti)
#        geop1000 = d.geop
#        Temp1000 = d.T
        
        # Draw contours for mslp
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        #
        #PlotSurfTemp(Lon,Lat, d.Temp1000 - 273.15, m)
        #PlotWind(d.lon, d.lat, d.u, d.v, map)
        PlotWindVelo(Lon, Lat, np.sqrt(d.u10m**2+ d.v10m**2), map, Umax= 25)
         
        plt.title('Arome '+str(time.localtime(d.tim[ti]-1)[:4]))
        plt.tight_layout()
        if save==True:
            plt.savefig(savedir+'Arome_'+exp_name+'_'+str(ti).zfill(2), bbox_inches='tight')

"""for cross section"""
#x_start, y_start= 10, 30
#x_end, y_end= 110, 20
#
#x, y= mapcross_section(d.lon, d.lat, x_start, y_start, x_end, y_end, m)


"""potential temperature at 700 hPa and stability"""
#if PlotTheta700== True:
#    fignr+=1
#    plt.figure(fignr)
#    plt.clf()
#    
#    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
#    
#    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
#    
#    pn= 6 #6= 700hPa
#    d.imp_level(pn= pn, tn= t)
#    
#    theta= PotTemp(d.Temp, d.pres[pn])
#    PlotColorMap3(Lon, Lat, theta, map, bounds= pot_temp_bounds, label=r"theta 700")
#    
#    plt.title('Arome '+str(time.localtime(d.tim[t])[:4]) +'theta 700 and static stability')
#    
#    
#    #"""stability"""
#    pn= 8 # 8 = 500hPa
#    d.imp_level(pn= pn, tn= t)
#    
#    theta500= PotTemp(d.Temp, d.pres[pn])
#    thetaSST= Temp1000
#    PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')


""" equivalent potential temperature at 850hPa"""
if PlotTheta_e850 == True:
    fignr+=1
    plt.figure(fignr)
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)

    d.imp_surf(tn=t)    
    pn= 4 #850hPa
    d.imp_level(pn= pn, tn= t)
    
#    theta= d.Temp*(1000/d.pres[pn])**(2/7)
    spechum= d.RH/(.263 *d.pres[pn])* np.exp(17.67*(d.Temp-273.15) / (d.Temp-29.65))
#    theta_e= theta* np.exp(2501 * spechum/(1.006* d.Temp))
    theta_e= EquiPotTemp(d.Temp, spechum, d.pres[pn])
    
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
#theta= d.Temp*(1000/d.pres[pn])**(2/7)
#spechum= d.RH/(.263 *d.pres[pn])* np.exp(17.67*(d.Temp-273.15) / (d.Temp-29.65))
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
##PlotCross_sec_T(d.pres[:-5], d.lon[y,x], d.Temp[:-5])
##
##
##"""cross section potential temperature"""
##plt.figure(fignr)
##plt.clf()
##fignr += 1
##
##plt.title('potential temperature')
##
##presm= np.tile(d.pres, (np.shape(d.Temp)[1],1)).T #make a matrix with the pressure in every row
##theta= d.Temp*(1000/presm)**(2/7)
##
##PlotCross_sec_T(d.pres[:-5], d.lon[y,x], theta[:-5])
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
###spechum= relhum[0to1]*100/(26.3 * pres[hPa]) * exp(17.67*(T[k]-273.15)/(T[K]-29.65))
##spechum= d.RH/(.263 *presm)* np.exp(17.67*(d.Temp-273.15) / (d.Temp-29.65))
##
###equivalent pot temp theta_e = theta* exp(Lc * spechum/(Cp * T))
##theta_e= theta* np.exp(2501 * spechum/(1.006* d.Temp))
##
##PlotCross_sec_T(d.pres[:-5], d.lon[y,x], theta_e[:-5])
##
##PlotCross_sec_hum(d.pres[:-5], d.lon[y,x], d.RH[:-5])
##
##fignr= 5
##"""vertical velocity"""
###plt.figure(fignr)
###plt.clf()
###fignr += 1
###
###PlotCross_sec_vervel(d.pres, d.lon[y,x], d.w)
##
##
##fignr= 6
##"""horizontal velocity (north- south)"""
###plt.figure(fignr)
###plt.clf()
###fignr += 1
###
###PCross_sec_horvel(d.pres, d.lon[y,x], d.u, d.v)
###
###
###plt.figure(fignr)
###plt.clf()
###fignr += 1
###
###PCross_sec_horvel2(d.pres, d.lon[y,x], d.u, d.v)
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