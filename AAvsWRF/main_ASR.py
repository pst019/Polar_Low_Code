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
Mediadir= '/media/'+user+'/PatsOrange/'



sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')
from f_imp_ASR import *
from f_impLists import *
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_meteo import *

#from random import randint


"""global variables"""
#var = 'PressWind'
var= 'PressWind_advanced'
#var = 'Theta700'

fignr= 3

maptype='AA_half'
#maptype='ASR'

if maptype=='AA': plt.figure(fignr, figsize= (6, 6))
elif maptype=='AA_half': plt.figure(fignr, figsize= (6, 4.7))
else: plt.figure(fignr, figsize= (6, 4.5))
fignr += 1
plt.clf()




"""Lambert coordinates"""
#lllon, lllat, urlon, urlat= -15, 63, 60, 75
#lat0, lon0= 75, 0 #(lat0, lon0) = center point

if maptype== 'AA': map = AA_map()
elif maptype == 'ASR': map= ASR_map(res='c')
elif maptype== 'AA_half': map = AA_map_half()
else: map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)


"""time"""
year=2008
month=3
day= 3
hour= 6

"""c) plot characteristics"""
#number_plevels= 20         

Track= False

"""end global variables """

save= True #this has to be further developped
savedir=   homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/Fields/' 
title_extra=''


"""module starts"""
h3 = hour// 3
t= (day-1)*8 + h3
   


    
#TPL= TRACK_PLlist_2(year, month, name='7', model='ASR')
    
print(year, month, day, hour, 't since start of month: ', t)


d= dataASR(year, month, day, sh3= h3, eh3= h3)          

Lon, Lat= map(d.lon, d.lat)



"""ASR - wind and surface pressure"""
if var== 'PressWind':
    d.impvar('SLP', level='surf')
    d.impvar('U', level='surf')
    d.impvar('V', level='surf')
    
    U10= np.sqrt(d.u10[0]**2+d.v10[0]**2)
    PlotWindVelo(Lon, Lat, U10, map, Umax= 25)
    
    PlotWind(d.lon, d.lat, d.u10[0], d.v10[0], map, nx= 30, ny= 30, rot= False, option='ASR')
    PlotContours(Lon, Lat, d.SLP[0], map, leveldist= 1)
    

if var== 'PressWind_advanced':
    d.impvar('SLP', level='surf')
    d.impvar('U', level='surf')
    d.impvar('V', level='surf')

    PlotContours(Lon, Lat, d.SLP[0], map, leveldist= 1, numbers=False)
        
    U10= np.sqrt(d.u10[0]**2+d.v10[0]**2)
    PlotWindVelo(Lon, Lat, U10, map, Umax= 25, color='YlBu')
    
    PlotLocalMax(d.SLP[0], threshold=1010, distance=30, map= map, lon=d.lon, lat=d.lat,
                     data2=U10, threshold2=10, distance2= 10) 

    PlotLocalMax(U10, threshold=15, distance=40, map= map, lon=d.lon, lat=d.lat, typ='max',
                 color='orange', dot=False, roundorder= 0, latbound= [65, 80])  
    
"""ASR- theta 700 and static stability"""
if var == 'Theta700':
    d.impvar('SLP', level='surf')   
    d.impvar('T', level=700)
    d.impvar('SST', level='surf')
    d.SST[d.SST<272]= np.nan #this sets the area with ice to land
    d.SST= np.ma.array(d.SST, mask= np.isnan(d.SST))
    d.impvar('T', level=500)
    
    Theta= PotTemp(d.T700, plev= 700)
    
    PlotContours(Lon, Lat, d.SLP[0], map, leveldist= pleveldist)
    bounds= np.arange(264, 290, 1)
    PlotColorMap3(Lon, Lat, Theta[0], map, symetric=False, color='RdBu', bounds= bounds, label= r"$\theta_{700}$ [K]")
    
    RCp= 2/7 #R/Cp
    theta500= d.T500[0]*(1000/500)**RCp
    thetaSST= d.SST[0]*(1000/d.SLP[0])**RCp
    PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')
    
    

if Track== True:
    PlotTRACK_PL(TPL, t, map, track=True, nowPL=False)



if save== False:
    plt.title('ASR '+str(year)+'_'+str(month).zfill(2)+'_'+str(day).zfill(2)+'_'+str(hour).zfill(2)+ title_extra)
    #datetime is shifted by one day

plt.tight_layout()

if save== True:
    savename= savedir + 'ASR_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+str(hour).zfill(2)+'_'+var+title_extra
    plt.savefig(savename, pad_inches=0)
    print(savename)