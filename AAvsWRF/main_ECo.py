#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 13:25:51 2017

@author: pst019
from ASR-ERA_Inv_random_3 of folder ASR
"""

#global user
#user='patricks'
##user='pst019'

import pickle
import sys
import os
user = os.getcwd().split('/')[2]
#print('user: ', user)

homedir= '/home/'+user+'/home/'

sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')

from f_meteo import *
from f_imp_ECo import *
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)

from random import randint


"""global variables"""
PlotPressWind= True
PlotTheta700= False
PlotVorticity= False
PlotTheta_e850= False

fignr= 14

"""Lambert coordinates"""
#lllon, lllat, urlon, urlat= -15, 63, 60, 75
#lat0, lon0= 75, 0 #(lat0, lon0) = center point

"""time"""
year=2008
month=3
day= 3
hour= 12

"""c) plot characteristics"""
#number_plevels= 20         
#pot_temp_bounds= np.arange(264, 290, 1)

"""end global variables """

teco= (day-2)*4 + hour//6 -1 #the ECo file starts at 6

do= dataECo(tstart= teco, tend=teco+1, timeperiod='2008.03.02-04')



"""ERA wind and surface pressure"""
if PlotPressWind== True:
    plt.figure(fignr)
    fignr +=1

    plt.clf()
    #plt.subplot(2, 3, 2)
    
   
    U1000= np.sqrt(do.uPL[0, 0]**2+do.vPL[0, 0]**2)
    
    #map= Polar_map(latsouth= de.lat[-1])
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
    grid= np.meshgrid(do.lon, do.lat) #pcolormesh needs corner points for projection
    Lon, Lat= map(grid[0], grid[1])
    #
    PlotWindVelo(Lon, Lat, U1000, map, Umax= 25)
    
    #a first approximation transformation from 1000hPa geopotential height to mslp (https://www.sensorsone.com/altitude-pressure-units-conversion/)
    SLP= 1000+ do.Geop[0,0]/9.81 *0.1214
    
    PlotWind(do.lon, do.lat, do.uPL[0,0], do.vPL[0,0], map)
    PlotContours(Lon, Lat, SLP, map, leveldist= pleveldist)
    
    PlotContours(Lon, Lat, do.Geop[0,3]/9.81 , map, nrlevels= 10, color='r')
    #
    #plt.title('ERA surface wind and pressure and geopotential height') #+str(t//8)+'_'+str(t%8*3))
    plt.title('ECo day'+ str(2+ do.tim[0]//24)+'-hour-'+str(do.tim[0]%24))
#    PlotTRACK_PL(TPLera, t//2, map, nowPL=False, track=True)

#
"""potential temperature at 700 hPa  + static stability"""
#if PlotTheta700== True:
#    plt.figure(fignr)
#    fignr +=1
#    plt.clf()
#    
#    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
#    
#    PlotContours(Lon, Lat, SLP, map, leveldist= pleveldist)
#    
#    Theta= PotTemp(do.T[0, 2], plev= 700)
#    PlotColorMap3(Lon, Lat, Theta, map, symetric=False, color='RdBu', bounds= pot_temp_bounds, label= r"Potential temperature at 700hPa [K]")
#    plt.title('ECo day'+ str(1+ do.tim[0]//24)+'-hour-'+str(do.tim[0]%24))#    
##    theta500= de.T[0,2]*(1000/500)**RCp
##    thetaSST= de.SST[0]*(1000/de.MSLP[0])**RCp
##    PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')
#    
#    PlotTRACK_PL(TPLera, t//2, map, track=False)
#
#
"""equivalent potential temperature at 850 hPa"""
if PlotTheta_e850== True:
    plt.figure(fignr)
    fignr +=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
    
    Theta_e= EquiPotTemp(do.T[0,1], do.sHum[0,1], plev= 850)
    
    PlotContours(Lon, Lat, SLP, map, leveldist= pleveldist)

    PlotColorMap3(Lon, Lat, Theta_e, map, symetric=False, color='RdBu', bounds= np.arange(264, 292, 1), label= r"$\theta_{e, 850}$ [K]")
    plt.title('ECo day'+ str(2+ do.tim[0]//24)+'-hour-'+str(do.tim[0]%24))    
    

