#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 15:06:41 2017

@author: pst019
"""
import os
user = os.getcwd().split('/')[2]

from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import sys  #to import the functions from a different directory

sys.path.insert(0, '/home/'+user+'/polar_low_code/Functions/')
from f_imp_ERA2 import *
from f_impLists import *
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
import scipy.ndimage.filters as filters


fignr= 1

syear= 1979
eyear= 2017 #this is not included

importdata= True #True #import and filter

hemisp= 'NH'
#hemisp = 'SH' #have to change it manually in some places

print('set intense = True')
#intense= True
index= 9 #9= low depth
#index= 6 #6= vortfitered
var='large'

"""import the data"""
if importdata == True:
    PLlist= np.zeros((14, 0))
    PLlistSH= np.zeros((11, 0))

    #shear= []
    #latpos= []
    #lonpos= []
    
    if hemisp == 'NH':
        d= data(['SST'], year=syear, month= 12)
        latbound= d.lat[-1]
        Tname='7'
        hemi=''
    elif hemisp == 'SH': 
        d= data(['SST'], year=syear, month= 6, hemisp='SH')
        latbound= d.lat[0]
        Tname='7'
        hemi='_SH_'
    
    
    for y, year in enumerate(range(syear, eyear)):
        TPL=TRACK_PLlist_2(year, month='all', hemi= '', name='7')
        TPLSH=TRACK_PLlist_2(year, month='all', hemi= '_SH_', name='7')

        
        #just take the PL points
        if intense== False:
            PLlist= np.hstack((PLlist, TPL.PLlist[:, TPL.PLlist[-1] == 1]))
            PLlistSH= np.hstack((PLlistSH, TPLSH.PLlist[:, TPLSH.PLlist[-1] == 1]))

            

print('shape of the PLlist: ', PLlist.shape )

unicnumber= PLlist[0]*1E8 + PLlist[1]*1E6+ PLlist[2]
            
data_hist = []
for n in remove_dublicate(unicnumber):
    if var=='small': data_hist += [np.min(PLlist[index, unicnumber == n])]
    else:data_hist+= [np.max(PLlist[index, unicnumber == n])]


unicnumberSH= PLlistSH[0]*1E8 + PLlistSH[1]*1E6+ PLlistSH[2]

data_histSH = []
for n in remove_dublicate(unicnumberSH):
    if var=='small': data_histSH += [np.min(PLlistSH[index, unicnumberSH == n])]
    else:data_histSH+= [np.max(PLlistSH[index, unicnumberSH == n])]  
print('number of PLs', len(data_histSH), len(remove_dublicate(unicnumberSH)))

data_hist= np.array(data_hist)
data_histSH= np.array(data_histSH)

if index== 9:
    data_hist *= -100
    data_histSH *= -100

"""plot the histograms"""
plt.figure(fignr)
plt.clf()

ax= plt.figure(fignr).gca()

bindist= (np.max(data_hist)-np.min(data_hist))//15 +1

if index == 9: bindist = 5
elif index == 6: bindist= 1   
mindata= np.round(np.min([np.min(data_hist), np.min(data_hist)]) -0.5)
offset= mindata%bindist #this is done such that the first bin would always start at 0
mindata -= offset

bins= np.arange(mindata, np.round(np.max(data_hist)+0.5), bindist)
bins= np.arange(mindata, np.round(np.max(data_hist)+0.5), bindist)

#weights = np.ones_like(data_track)/(len(data_track)*bindist)
#plt.hist(np.array(data_track), bins=bins_track, alpha= 0.7, weights=weights*100, label='Cyclones')

weights= np.ones_like(data_hist)/(len(data_hist)*bindist)
weightsSH= np.ones_like(data_histSH)/(len(data_histSH)*bindist)

plt.hist(np.array(data_hist), bins=bins, alpha= 0.7, weights=weights*100,  label='Polar lows NH')
plt.hist(np.array(data_histSH), bins=bins, alpha= 0.7, weights=weightsSH*100,  label='Polar lows SH')

plt.legend(fontsize= 16)

plt.ylabel('Frequency')
if index== 9:
    plt.xlim([40, 180])
    plt.xlabel(r'$\overline{SLP}$ - SLP [Pa]')
#  
#"""plot unsmoothed count density"""
##plt.figure(fignr)
##fignr+=1
##plt.clf()
##map= Polar_map(latsouth= latbound, hemisp=hemisp)
##
##grid= np.meshgrid(d.lon, d.lat) #pcolormesh needs corner points for projection
##Lon, Lat= map(grid[0], grid[1])
##
##PlotColor_from0(Lon, Lat, 1/(eyear-syear) *np.sum(countdens, axis= (0,1)), map, label= 'Counts/year in grid cell')
#
#
#
#"""filtering of the count density"""
#if importdata == True:
#    countdens_filter= np.zeros((np.shape(countdens)))
#    
#    latdist= 4
#    for ilat in range(len(d.lat)):  #loop through every PL point      
#        londist= int(latdist//np.cos(np.deg2rad(d.lat[ilat])))
#    #    print(d.lat[ilat], londist)
#        
#        for ilon in range(len(d.lon)):       
#            mask= radMask2D((len(d.lat), len(d.lon)), (ilat, ilon), radiusx= latdist, radiusy= londist) #lat lon mask
#            if ilon == 0:
#                masksize= len(np.where(mask ==True)[0])
#    #            print(d.lat[ilat], londist, masksize)
#                
#            countdens_filter[:,:, ilat, ilon] = np.sum(countdens[:,:, mask], axis= -1) #* 1/masksize
#    
#    
#    print(np.sum(countdens_filter))
#    
#    
#    countd_filter_annual= 1/(eyear-syear) *np.sum(countdens_filter, axis= (0,1))
#
#
#"""plot filtered count dens"""
#
#plt.figure(fignr)
#fignr+=1
#plt.clf()
#map= Polar_map(latsouth= latbound, hemisp=hemisp)
#grid= np.meshgrid(d.lon, d.lat) #pcolormesh needs corner points for projection
#Lon, Lat= map(grid[0], grid[1])
#
#
#maxlevel= np.max(np.abs(countd_filter_annual))
##if hemisp== 'NH': bounds= [0, 0.2, 1, 3, 10, 20, np.round(maxlevel + 0.5)]
##elif hemisp== 'SH': bounds= [0, 0.2, 1, 3, np.round(maxlevel + 0.5)]
##bounds= [0, 0.3, 1, 3, 10, 20, 33]
#
#bounds= [0, 6, 12, 24, 48, 96, 168]
#if intense== True: bounds= [0,  1, 3, 6, 12, 20, 30]
#
#PlotColor_setbounds(Lon, Lat, countd_filter_annual*6, map, bounds=bounds, label='Annual polar low hours in 220 km radius')
# 
#plt.tight_layout()
##Ice=np.mean(d.SST, axis= 0)
##PlotIce(Lon, Lat, Ice, map)
#
#
#def makeBox(boxcoord, fignr, map):
#    """plot the box of boxcoord= [west, east, north, south] - should not go over the "end" of the dataset"""
#    west, east, north, south= boxcoord
#    if east <0: east += 360
#    if west <0: west += 360      
#    
#    n= int(east- west)
#    x,y = map(np.concatenate(([west], np.linspace(west, east, n), [east],  np.linspace(east, west,n))),
#              np.array([south]+ [north]*n+ [south]+ [south]*n))
#    plt.figure(fignr)        
#    map.plot(x, y, color='b', linewidth= 1.5)
#
#
#if hemisp== 'NH':
#    boxGreen= [300, 360, 75, 50]
#    makeBox(boxcoord= boxGreen, fignr = fignr-1, map= map)
#    Sxpt, Sypt= map(-30, 50)
#    plt.text(Sxpt, Sypt, '2', fontsize=15, ha='center',va='top')
#    
#    boxNordic= [0, 50, 80, 55]
#    makeBox(boxcoord= boxNordic, fignr = fignr-1, map= map)
#    Sxpt, Sypt= map(39, 58)
#    plt.text(Sxpt, Sypt, '1', fontsize=15, ha='center',va='top')
#    
#    boxPaz= [145, 200, 65, 40]
#    makeBox(boxcoord= boxPaz, fignr = fignr-1, map= map)
#    Sxpt, Sypt= map(175, 35)
#    plt.text(Sxpt, Sypt, '3', fontsize=15, ha='center',va='top')
#    
#    boxAla= [200, 235, 62, 48]
#    makeBox(boxcoord= boxAla, fignr = fignr-1, map= map)
#    Sxpt, Sypt= map(-145, 41)
#    plt.text(Sxpt, Sypt, '4', fontsize=15, ha='center',va='top')
#    
#    boxJap= [130, 142, 50, 35]
#    makeBox(boxcoord= boxJap, fignr = fignr-1, map= map)
#    Sxpt, Sypt= map(145, 30)
#    plt.text(Sxpt, Sypt, '5', fontsize=15, ha='center',va='top')    
#    
#
#elif hemisp == 'SH':
#    boxAmu= [230, 300, -45, -70]
#    makeBox(boxcoord= boxAmu, fignr = fignr-1, map= map) 
#    Sxpt, Sypt= map(-90, -42)
#    plt.text(Sxpt, Sypt, '1', fontsize=15, ha='center',va='top')
#    
#    boxNZ= [160, 200, -45, -70]
#    makeBox(boxcoord= boxNZ, fignr = fignr-1, map= map)
#    Sxpt, Sypt= map(-175, -42)
#    plt.text(Sxpt, Sypt, '2', fontsize=15, ha='center',va='top')
#    
#    boxA= [80, 140, -45, -65]
#    makeBox(boxcoord= boxA, fignr = fignr-1, map= map)
#    Sxpt, Sypt= map(110, -42)
#    plt.text(Sxpt, Sypt, '3', fontsize=15, ha='center',va='top')
#
#
#
#
#
#    
#"""last minus first years"""
#nyears= 15
#countd_filter_diff= 1/nyears *(np.sum(countdens_filter[-nyears:], axis= (0,1)) - np.sum(countdens_filter[:nyears], axis= (0,1)) )
#
#
#"""significance calculations"""
#import scipy.stats
#statistic, pvalue= scipy.stats.ttest_ind(np.sum(countdens_filter[-nyears:], axis= (1)), np.sum(countdens_filter[:nyears], axis= (1)), axis= (0))
#
#
#plt.figure(fignr)
#fignr+=1
#plt.clf()
#map= Polar_map(latsouth= latbound, hemisp=hemisp)
#grid= np.meshgrid(d.lon, d.lat) #pcolormesh needs corner points for projection
#Lon, Lat= map(grid[0], grid[1])
#
#pvalue[np.isnan(pvalue)]= 1
##PlotContours(Lon, Lat, pvalue, map, levels=[0, 0.1])
#
##PlotColorMap3(Lon, Lat, countd_filter_diff, map, symetric= True, maxlevel= maxlevel, boxnr= 9, label='Counts change / year')
#
##plot only significant values
#countd_filter_diff[pvalue > 0.05]= 0
#
#if hemisp == 'NH': bounds= [-48, -36, -24, -12, 12, 24, 36, 48]
#else: bounds= [-24, -18, -12, -6, 6, 12, 18, 24]
#if intense== True: bounds= [-12, -8, -4, -2, 2, 4, 8, 12]
#
#PlotColorMap3(Lon, Lat, countd_filter_diff*6, map, bounds=bounds, label='Average annual difference in polar low hours')
#
#plt.tight_layout()
#
#
#"""seasonal distribution of PLs"""
#seasonaldistr= np.sum(countdens, axis= (0,2,3))/(eyear-syear)
#seasonnorm=np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
#seasonaldistr *= 30/(seasonnorm*4)
#
#plt.figure(fignr)
#fignr+=1
#plt.clf()
#ax= plt.figure(fignr-1).gca()
#
#from matplotlib.ticker import MaxNLocator
#ax.yaxis.set_major_locator(MaxNLocator(integer=True))
#
#bins= np.arange(1, 13)
#plt.bar(bins, seasonaldistr, label='Total') #, width=1)
#plt.xticks(bins)
#
#plt.xlim([.5,12.5])
#plt.xlabel('Month', size= 14)
#plt.ylabel('Average monthly polar low duration [days]', size= 14)
#plt.tick_params(labelsize= 14)
#
#"""seasonal distr of boxes"""
#if hemisp== 'NH':
#    seasonaldistrNordic= np.sum(countdens[:,:, int(np.where(d.lat== boxNordic[2])[0]): int(np.where(d.lat== boxNordic[3])[0]), int(np.where(d.lon== boxNordic[0])[0]): int(np.where(d.lon== boxNordic[1])[0])], axis= (0,2,3))/(eyear-syear)
#    seasonaldistrNordic *= 30/(seasonnorm*4)
#    plt.bar(bins, seasonaldistrNordic, label= 'Nordic Seas') #, width=1)
#    
#    seasonaldistrGreen= np.sum(countdens[:,:, int(np.where(d.lat== boxGreen[2])[0]): int(np.where(d.lat== boxGreen[3])[0]), int(np.where(d.lon== boxGreen[0]-360)[0]): int(np.where(d.lon== boxGreen[1]-360)[0])], axis= (0,2,3))/(eyear-syear)
#    seasonaldistrGreen *= 30/(seasonnorm*4)
#    plt.bar(bins, seasonaldistrGreen, bottom= seasonaldistrNordic, label='North Central Atlantic') #, width=1)
#
#    seasonaldistrPaz= np.sum(countdens[:,:, int(np.where(d.lat== boxPaz[2])[0]): int(np.where(d.lat== boxPaz[3])[0]), list(np.arange(int(np.where(d.lon== boxPaz[0])[0]), len(d.lon))) + list(np.arange(0, int(np.where(d.lon%360== boxPaz[1])[0])))], axis= (0,2,3))/(eyear-syear)
#    seasonaldistrPaz *= 30/(seasonnorm*4)
#    plt.bar(bins, seasonaldistrPaz, label='North West Pacific', bottom= seasonaldistrNordic+seasonaldistrGreen) #, width=1)
#    
#    seasonaldistrAla= np.sum(countdens[:,:, int(np.where(d.lat== boxAla[2])[0]): int(np.where(d.lat== boxAla[3])[0]), int(np.where(d.lon== boxAla[0]-360)[0]): int(np.where(d.lon== boxAla[1]-360)[0])], axis= (0,2,3))/(eyear-syear)
#    seasonaldistrAla *= 30/(seasonnorm*4)
#    plt.bar(bins, seasonaldistrAla, label='Gulf of Alaska', bottom= seasonaldistrNordic+seasonaldistrGreen+seasonaldistrPaz) #, width=1)
#        
#    seasonaldistrJap= np.sum(countdens[:,:, int(np.where(d.lat== boxJap[2])[0]): int(np.where(d.lat== boxJap[3])[0]), int(np.where(d.lon== boxJap[0])[0]): int(np.where(d.lon== boxJap[1])[0])], axis= (0,2,3))/(eyear-syear)
#    seasonaldistrJap *= 30/(seasonnorm*4)
#    plt.bar(bins, seasonaldistrJap, label='Sea of Japan', bottom= seasonaldistrNordic+seasonaldistrGreen+seasonaldistrPaz+seasonaldistrAla) #, width=1)
#
#    plt.legend(fontsize= 13)
#    
#elif hemisp == 'SH':
#    seasonaldistrAmu= np.sum(countdens[:,:, int(np.where(d.lat== boxAmu[2])[0]): int(np.where(d.lat== boxAmu[3])[0]), int(np.where(d.lon== boxAmu[0]-360)[0]): int(np.where(d.lon== boxAmu[1]-360)[0])], axis= (0,2,3))/(eyear-syear)
#    seasonaldistrAmu *= 30/(seasonnorm*4)
#    plt.bar(bins, seasonaldistrAmu, label= 'SE Pacific Ocean') #, width=1)   
#
#    seasonaldistrNZ= np.sum(countdens[:,:, int(np.where(d.lat== boxNZ[2])[0]): int(np.where(d.lat== boxNZ[3])[0]), list(np.arange(int(np.where(d.lon== boxNZ[0])[0]), len(d.lon))) + list(np.arange(0, int(np.where(d.lon%360== boxNZ[1])[0]))) ], axis= (0,2,3))/(eyear-syear)
#    seasonaldistrNZ *= 30/(seasonnorm*4)
#    plt.bar(bins, seasonaldistrNZ, bottom= seasonaldistrAmu, label= 'SW Pactic Ocean') #, width=1)  
#
#    seasonaldistrA= np.sum(countdens[:,:, int(np.where(d.lat== boxA[2])[0]): int(np.where(d.lat== boxA[3])[0]), int(np.where(d.lon== boxA[0])[0]): int(np.where(d.lon== boxA[1])[0])], axis= (0,2,3))/(eyear-syear)
#    seasonaldistrA *= 30/(seasonnorm*4)
#    plt.bar(bins, seasonaldistrA, bottom= seasonaldistrAmu+seasonaldistrNZ, label= 'SE Indian Ocean') #, width=1)  
#        
#    plt.legend(fontsize= 13, loc= 'upper left')
#    
#plt.tight_layout()
#
#
#"""annual distribution of PLs"""
#annualdistr= np.sum(countdens, axis= (1,2,3))
#
#plt.figure(fignr)
#fignr+=1
#plt.clf()
#bins= np.arange(syear, eyear)
#plt.bar(bins, annualdistr/4, label='Total') #, width=1)
#
##plt.plot(np.arange(syear, eyear), annualdistr)
#plt.xlabel('Year', size= 14)
#plt.ylabel('Annual polar low duration [days]', size= 14)
#plt.xlim([syear-.5,eyear-.5])
#plt.tick_params(labelsize= 14)
#
#if hemisp == 'NH':
#    adistrNordic= np.sum(countdens[:,:, int(np.where(d.lat== boxNordic[2])[0]): int(np.where(d.lat== boxNordic[3])[0]), int(np.where(d.lon== boxNordic[0])[0]): int(np.where(d.lon== boxNordic[1])[0])], axis= (1,2,3))
#    plt.bar(bins, adistrNordic/4, label= 'Nordic Seas') #, width=1)
#    
#    adistrGreen= np.sum(countdens[:,:, int(np.where(d.lat== boxGreen[2])[0]): int(np.where(d.lat== boxGreen[3])[0]), int(np.where(d.lon== boxGreen[0]-360)[0]): int(np.where(d.lon== boxGreen[1]-360)[0])], axis= (1,2,3))
#    plt.bar(bins, adistrGreen/4, bottom= adistrNordic/4, label='North Central Atlantic') #, width=1)
#    
#    #adistrAla= np.sum(countdens[:,:, int(np.where(d.lat== boxAla[2])[0]): int(np.where(d.lat== boxAla[3])[0]), int(np.where(d.lon== boxAla[0]-360)[0]): int(np.where(d.lon== boxAla[1]-360)[0])], axis= (1,2,3))
#    #plt.bar(bins, adistrAla/4, label='Gulv of Alaska', bottom= (adistrNordic+adistrGreen)/4) #, width=1)
#    #
#    #adistrJap= np.sum(countdens[:,:, int(np.where(d.lat== boxJap[2])[0]): int(np.where(d.lat== boxJap[3])[0]), int(np.where(d.lon== boxJap[0])[0]): int(np.where(d.lon== boxJap[1])[0])], axis= (1,2,3))
#    #plt.bar(bins, adistrJap/4, label='Sea of Japan', bottom= (adistrNordic+adistrGreen+adistrAla)/4) #, width=1)
#
#elif hemisp == 'SH':
#    seasonaldistrAmu= np.sum(countdens[:,:, int(np.where(d.lat== boxAmu[2])[0]): int(np.where(d.lat== boxAmu[3])[0]), int(np.where(d.lon== boxAmu[0]-360)[0]): int(np.where(d.lon== boxAmu[1]-360)[0])], axis= (1,2,3))
#    plt.bar(bins, seasonaldistrAmu/4, label= 'SE Pacific Ocean') #, width=1)   
#
#    seasonaldistrNZ= np.sum(countdens[:,:, int(np.where(d.lat== boxNZ[2])[0]): int(np.where(d.lat== boxNZ[3])[0]), list(np.arange(int(np.where(d.lon== boxNZ[0])[0]), len(d.lon))) + list(np.arange(0, int(np.where(d.lon== boxNZ[1]-360)[0]))) ], axis= (1,2,3))
#    plt.bar(bins, seasonaldistrNZ/4, bottom= seasonaldistrAmu/4, label= 'SW Pacific Ocean') #, width=1)  
#
#    seasonaldistrA= np.sum(countdens[:,:, int(np.where(d.lat== boxA[2])[0]): int(np.where(d.lat== boxA[3])[0]), int(np.where(d.lon== boxA[0])[0]): int(np.where(d.lon== boxA[1])[0])], axis= (1,2,3))
#    plt.bar(bins, seasonaldistrA/4, bottom= (seasonaldistrAmu+seasonaldistrNZ)/4, label= 'SE Indian Ocean') #, width=1)  
#        
#
#plt.tight_layout()
#
#"""statistic"""
#from scipy import stats
#
#slope, intercept, r_value, p_value, std_err = stats.linregress(bins, annualdistr/4)
#plt.plot(bins, intercept+ slope*bins, color='b', label= 'Trend')
#print('trend: ', slope, 'pvalue of trend:', p_value)
#
#plt.legend(fontsize= 14)#, loc='lower left')
#
#
#
##plt.ylim([0, np.max(seasonaldistr)])
#
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
##    plt.savefig('/home/'+user+'/Dropbox/Polar_Low/Documents/ERA-PLs/Graphs/Countd/'+str(month)+'.png', dpi= 50)
#
##"""shear average for every pixel"""
##shear_avg= np.zeros((110, 720))
##for i in range(720):
##    shearlist= np.array(shear)[np.array(lonpos)== i]
##    latlist= np.array(latpos)[np.array(lonpos)== i]
##
##    for j in range(110):
##        sl= shearlist[latlist== j]
##        if len(sl) > 0:
###            print(sl)
##            shear_avg[j,i]= np.mean(shearlist[latlist== j])
# 
#
#"""shear average"""
##shear_avg_filter= np.zeros((110, 720))
##latdist= 4
##for ilat in range(0, 100):
##    londist= int(latdist//np.cos(np.deg2rad(d.lat[ilat])))
##    latlistmax=  np.where(np.array(latpos)< ilat+latdist)[0]
##    latlistmin=  np.where(np.array(latpos)> ilat-latdist)[0]
##    if len(latlistmin) >0 and len(latlistmax)> 0:
##        for ilon in range(0, 720):
##            lonlistmax=  np.where(np.array(lonpos)< ilon+londist)[0]
##            lonlistmin=  np.where(np.array(lonpos)> ilon-londist)[0]
##            
##            PLlist= [PLn for PLn in np.arange(len(latpos)) if PLn in latlistmin if PLn in latlistmax if PLn in lonlistmin if PLn in lonlistmax]
##            shear_avg_filter[ilat, ilon]= np.mean(np.array(shear)[PLlist])
##
##
##plt.figure(fignr)
##fignr+=1
##plt.clf()
##map= Polar_map(latsouth= d.lat[-1])
##grid= np.meshgrid(d.lon, d.lat) #pcolormesh needs corner points for projection
##Lon, Lat= map(grid[0], grid[1])
##
##plt.tight_layout()
##
##PlotColor_from0(Lon, Lat, shear_avg_filter, map, maxlevel= 180, label='avg shear of pixel')
