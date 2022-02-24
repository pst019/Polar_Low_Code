#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 15:06:41 2017

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/1692A00D929FEF8B/'


from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import sys  #to import the functions from a different directory

sys.path.insert(0, homedir +'Polar_Low/Code2/Functions/')
from f_imp_ERA2 import *
from f_imp_ASR import *
from f_impLists import *
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)

import scipy.ndimage.filters as filters

fignr= 1

syear= 2000
eyear= 2013  #this is not included


#latpos= []
#lonpos= []

#hemisp= 'NH'
model='ASR'

importdata= True

intense=True
pdiffintense= 5.071 #hPa, values are for mean(SLP)-SLP
stabtheta500intense= -3.978

"""import the data"""
if importdata== True:
    d= dataASR(syear, month=1, sday=1, level='surf')
    
    countdens= np.zeros(((eyear-syear), 12, 720, 720))
    
    
    for y, year in enumerate(range(syear, eyear)):
        for month in range(1, 13):
            TPL=TRACK_PLlist_2(year, month, model='ASR', name='7')
            
            #just take the PL points
            if intense== False:
                PLlist= TPL.PLlist[:, TPL.PLlist[-1] == 1]
            elif intense== True:
                PLlist= TPL.PLlist[:, np.logical_and.reduce((TPL.PLlist[-1] == 1, TPL.PLlist[-4] > pdiffintense, TPL.PLlist[-3] > stabtheta500intense))] 
                
            for i in range(len(PLlist[0])):
                minval= np.min((PLlist[3][i] - d.lat)**2 + (PLlist[2][i] - d.lon)**2)
                pts2D= np.where((PLlist[3][i] - d.lat)**2 + (PLlist[2][i] - d.lon)**2 == minval)
    
                countdens[y, month-1, pts2D[0], pts2D[1]] +=1
            
    
      
    print('avg annual counts: ', 1/(eyear-syear) * np.sum(countdens) )
    
    countdens_spatial= 1/(eyear-syear) *np.sum(countdens, axis= (0,1))
  
"""unsmoothed count density"""
plt.figure(fignr)
fignr+=1
plt.clf()
map= ASR_map()
Lon, Lat= map(d.lon, d.lat)

PlotColor_from0(Lon, Lat, countdens_spatial, map, label= 'Counts/year in grid cell')



"""smoothed count density new categories"""
if importdata== True:
    print('avg annual counts before filtering: ', np.sum(countdens_spatial))

    filterdist= 220 / 15
    
    import scipy.ndimage.filters as filters
    #box filter
#    countd_filter= filters.uniform_filter(countdens_spatial, size= (filterdist, filterdist), mode='constant')
#    countd_filter *= filterdist**2 *np.pi**2/4  #to get same scale as for circles

    #gaussian filter- used in first draft
#    countd_filter= filters.gaussian_filter(countdens_spatial, sigma= filterdist, mode='constant', truncate= 1.)
    #to transfer it to PLs per area
#    countd_filter *= filterdist**2 *np.pi
    
#    print('avg annual counts after filtering: ', np.sum(countd_filter)/(np.pi*filterdist**2))


  
    #used later
    countd_filter= np.zeros((d.lat.shape))
    for xi in range(d.lat.shape[0]):
        for yi in range(d.lat.shape[1]):
            mask= radMask2D(d.lat.shape, (xi, yi), radiusx= filterdist, radiusy= filterdist, continuous_y=False)
            
            countd_filter[xi, yi]= np.sum(countdens_spatial[mask])
    print('avg annual counts after filtering: ', np.sum(countd_filter)/(3.14*filterdist**2))


plt.figure(fignr)
fignr+=1
plt.clf()
map= ASR_map()
Lon, Lat= map(d.lon, d.lat)

plt.tight_layout()

maxlevel= np.max(np.abs(countd_filter))
#bounds= [0, 1, 2, 5, 10, 20, np.round(maxlevel + 0.5)]
bounds= [0, 6, 12, 24, 48, 120, 216]
if intense== True: bounds= [0,  1, 3, 6, 12, 20, 30]

#elif hemisp== 'SH': bounds= [0, 0.03, 0.1, 0.3, 1, np.round(maxlevel + 0.5)]
PlotColor_setbounds(Lon, Lat, countd_filter*3, map, bounds=bounds, label='Annual polar low hours in 220 km radius')
 
    
#Ice=np.mean(d.SST, axis= 0)
#PlotIce(Lon, Lat, Ice, map)


def makeBox(boxcoord, fignr, map):
    """plot the box of boxcoord= [west, east, north, south]"""
    west, east, north, south= boxcoord
    
    if east <0: east += 360
    if west <0: west += 360      
    
    n= int(east- west)
    x,y = map(np.concatenate(([west], np.linspace(west, east, n), [east],  np.linspace(east, west,n))),
              np.array([south]+ [north]*n+ [south]+ [south]*n))
    plt.figure(fignr)        
    map.plot(x, y, color='b', linewidth= 1.5)

boxNordic= [0, 50, 80, 55]
makeBox(boxcoord= boxNordic, fignr = fignr-1, map= map)
Sxpt, Sypt= map(39, 59)
plt.text(Sxpt, Sypt, '1', fontsize=15, ha='center',va='top')
    
boxGreen= [300, 360, 75, 50]
makeBox(boxcoord= boxGreen, fignr = fignr-1, map= map)
Sxpt, Sypt= map(-45, 73)
plt.text(Sxpt, Sypt, '2', fontsize=15, ha='center',va='top')

boxPaz= [145, 200, 65, 40]
makeBox(boxcoord= boxPaz, fignr = fignr-1, map= map)
Sxpt, Sypt= map(157, 61)
plt.text(Sxpt, Sypt, '3', fontsize=15, ha='center',va='top')

boxAla= [200, 235, 62, 48]
makeBox(boxcoord= boxAla, fignr = fignr-1, map= map)
Sxpt, Sypt= map(-144, 44)
plt.text(Sxpt, Sypt, '4', fontsize=15, ha='center',va='top')
    
boxJap= [130, 142, 50, 35]
makeBox(boxcoord= boxJap, fignr = fignr-1, map= map)
Sxpt, Sypt= map(128, 44)
plt.text(Sxpt, Sypt, '5', fontsize=15, ha='center',va='top')

#"""last minus first years"""
#nyears= 15
#countd_filter_diff= 1/nyears *(np.sum(countdens_filter[-nyears:], axis= (0,1)) - np.sum(countdens_filter[:nyears], axis= (0,1)) )
#
#plt.figure(fignr)
#fignr+=1
#plt.clf()
#map= Polar_map(latsouth= latbound, hemisp=hemisp)
#grid= np.meshgrid(d.lon, d.lat) #pcolormesh needs corner points for projection
#Lon, Lat= map(grid[0], grid[1])
#
#
##PlotColorMap3(Lon, Lat, countd_filter_diff, map, symetric= True, maxlevel= maxlevel, boxnr= 9, label='Counts change / year')
#
#bounds= [-8, -6, -4, -2, 2, 4, 6, 8]
#bounds= [-4, -3, -2, -1, 1, 2, 3, 4]
#PlotColorMap3(Lon, Lat, countd_filter_diff, map, bounds=bounds, label='Counts change / year in 220 km radius')
#
#
"""seasonal distribution of PLs"""
seasonaldistr= np.sum(countdens, axis= (0,2,3))/(eyear-syear)
seasonnorm=np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
seasonaldistr *= 30/seasonnorm

plt.figure(fignr)
fignr+=1
plt.clf()

plt.xlabel('Month', size= 14)
plt.ylabel('Average monthly polar low duration [days]', size= 14)

bins= np.arange(1, 13)
plt.bar(bins, seasonaldistr/8, label='Total') #, width=1)
plt.xticks(bins)

plt.xlim([.5,12.5])
plt.tick_params(labelsize= 14)

maskNordic= np.where(np.logical_and.reduce((d.lon%360 > boxNordic[0], d.lon%360 < boxNordic[1], d.lat > boxNordic[3], d.lat < boxNordic[2])))
sdistrNordic= np.sum(countdens[:,:, maskNordic[0], maskNordic[1]], axis= (0,2))/(eyear-syear)
sdistrNordic *= 30/seasonnorm
plt.bar(bins, sdistrNordic/8, label='Nordic Seas')

maskGreen= np.where(np.logical_and.reduce((d.lon%360 > boxGreen[0], d.lon%360 < boxGreen[1], d.lat > boxGreen[3], d.lat < boxGreen[2])))
sdistrGreen= np.sum(countdens[:,:, maskGreen[0], maskGreen[1]], axis= (0,2))/(eyear-syear)
sdistrGreen *= 30/seasonnorm
plt.bar(bins, sdistrGreen/8, bottom= sdistrNordic/8, label='Denmark Strait') #, width=1)

maskPaz= np.where(np.logical_and.reduce((d.lon%360 > boxPaz[0], d.lon%360 < boxPaz[1], d.lat > boxPaz[3], d.lat < boxPaz[2])))
sdistrPaz= np.sum(countdens[:,:, maskPaz[0], maskPaz[1]], axis= (0,2))/(eyear-syear)
sdistrPaz *= 30/seasonnorm
plt.bar(bins, sdistrPaz/8, bottom= (sdistrNordic+ sdistrGreen)/8, label='Bering Sea') #, width=1)

maskAla= np.where(np.logical_and.reduce((d.lon%360 > boxAla[0], d.lon%360 < boxAla[1], d.lat > boxAla[3], d.lat < boxAla[2])))
sdistrAla= np.sum(countdens[:,:, maskAla[0], maskAla[1]], axis= (0,2))/(eyear-syear)
sdistrAla *= 30/seasonnorm
plt.bar(bins, sdistrAla/8, bottom= (sdistrNordic+ sdistrGreen + sdistrPaz)/8, label='Gulf of Alaska') #, width=1)
    
maskJap= np.where(np.logical_and.reduce((d.lon%360 > boxJap[0], d.lon%360 < boxJap[1], d.lat > boxJap[3], d.lat < boxJap[2])))
sdistrJap= np.sum(countdens[:,:, maskJap[0], maskJap[1]], axis= (0,2))/(eyear-syear)
sdistrJap *= 30/seasonnorm
plt.bar(bins, sdistrJap/8, bottom= (sdistrNordic+ sdistrGreen + sdistrPaz + sdistrAla)/8, label='Sea of Japan') #, width=1)
    
plt.legend(loc='upper center', fontsize=13)
plt.tight_layout()


"""PL annual distribution"""
annualdistr= np.sum(countdens, axis= (1,2,3)) #/(eyear-syear)

plt.figure(fignr)
fignr+=1
plt.clf()

bins= np.arange(syear, eyear)
plt.bar(bins, annualdistr/8, label='Total')
plt.xlim([syear-.5,eyear-.5])

plt.xlabel('Year', size= 14)
plt.ylabel('Annual polar low duration [days]', size=14)

adistrNordic= np.sum(countdens[:,:, maskNordic[0], maskNordic[1]], axis= (1,2))
plt.bar(bins, adistrNordic/8, label='Nordic Seas')

adistrGreen= np.sum(countdens[:,:, maskGreen[0], maskGreen[1]], axis= (1,2))
plt.bar(bins, adistrGreen/8, bottom= adistrNordic/8, label='Denmark Strait')

adistrPaz= np.sum(countdens[:,:, maskPaz[0], maskPaz[1]], axis= (1,2))
plt.bar(bins, adistrPaz/8, bottom= (adistrNordic+ adistrGreen)/8, label='North West Pacific') #, width=1)

adistrAla= np.sum(countdens[:,:, maskAla[0], maskAla[1]], axis= (1,2))
plt.bar(bins, adistrAla/8, bottom= (adistrNordic+ adistrGreen + adistrPaz)/8, label='Gulf of Alaska') #, width=1)
    
adistrJap= np.sum(countdens[:,:, maskJap[0], maskJap[1]], axis= (1,2))
plt.bar(bins, adistrJap/8, bottom= (adistrNordic+ adistrGreen + adistrPaz + adistrAla)/8, label='Sea of Japan') #, width=1)


if intense==False: plt.legend(fontsize= 13 , loc='lower left')
else: plt.legend(fontsize= 13 , loc='upper left')

plt.tick_params(labelsize= 14)
plt.tight_layout()
#"""count density for different months"""
##for month in [1, 2, 3, 4, 5, 9, 10, 11, 12]:
##    d= data(['SST'], year=syear, month= month)
##    
##    countdens_month= 1/(eyear-syear) * areafact[:, None] *np.sum(countdens[:, month-1], axis= 0)   
##
##    import scipy.ndimage.filters as filters
##    fact= 1000 #just a factor to make the filter work
##    countd_filter= filters.uniform_filter(countdens_month*fact, size= (4, 8), mode='wrap')/fact
##    print(np.sum(countd_filter))
##    
##    plt.figure(fignr)
##    fignr+=1
##    plt.clf()
##    map= Polar_map(latsouth= d.lat[-1])
##    grid= np.meshgrid(d.lon, d.lat) #pcolormesh needs corner points for projection
##    Lon, Lat= map(grid[0], grid[1])
##    
##    plt.tight_layout()
##    
##    PlotColor_from0(Lon, Lat, countd_filter,  map, maxlevel=0.5, label= 'Counts / year')
##    #PlotColor_from0(Lon, Lat, countd_filter, map, label='Counts / year')
##    plt.title(str(month))
##    
##    Ice=np.mean(d.SST, axis= 0)
##    PlotIce(Lon, Lat, Ice, map)
##    



"""filtered count density"""
#countdens_filter= np.zeros((np.shape(countdens)))
#
#latdist= 220 //12.5
#for ilat in range(len(d.lat)):  #loop through every PL point      
#    londist= int(latdist//np.cos(np.deg2rad(d.lat[ilat])))
##    print(d.lat[ilat], londist)
#    
#    for ilon in range(len(d.lon)):       
#        mask= radMask2D((len(d.lat), len(d.lon)), (ilat, ilon), radiusx= latdist, radiusy= londist) #lat lon mask
#        if ilon == 0:
#            masksize= len(np.where(mask ==True)[0])
##            print(d.lat[ilat], londist, masksize)
#            
#        countdens_filter[:,:, ilat, ilon] = np.sum(countdens[:,:, mask], axis= -1) #* 1/masksize
##
##
##print(np.sum(countdens_filter))
#
#countd_filter_annual= 1/(eyear-syear) *np.sum(countdens_filter, axis= (0,1))
#
#plt.figure(fignr)
#fignr+=1
#plt.clf()
#map= ASR_map()
#
#
#maxlevel= np.max(np.abs(countd_filter_annual))
#if hemisp== 'NH': bounds= [0, 0.2, 1, 3, 10, 20, np.round(maxlevel + 0.5)]
#elif hemisp== 'SH': bounds= [0, 0.2, 1, 3, np.round(maxlevel + 0.5)]
#bounds= [0, 0.2, 1, 3, 10, 20, 33]
#PlotColor_setbounds(Lon, Lat, countd_filter_annual, map, bounds=bounds, label='Counts / year in 220km radius')
 
    
##Ice=np.mean(d.SST, axis= 0)
##PlotIce(Lon, Lat, Ice, map)
#
#