#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 13:25:51 2017

@author: pst019
from Inv_time_TandS_4
"""

import os
user = os.getcwd().split('/')[2]


homedir= '/home/'+user+'/home/'

import pickle
import sys
sys.path.insert(0, code/Functions')
from f_imp_ERA2 import *
from f_impLists import *
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)

from random import randint

fignr= 4

year= randint(1979,2016)

year=2000

TPL=TRACK_PLlist_Year(year, name='3')


minlon= -90 #north atlantic
maxlon= 90 
#minlat= 30 #54
#maxlat= 70 #81

#just take the PL points
PLlist= TPL.PLlist[:, TPL.PLlist[-1] == 1]

PLlist= PLlist[:, np.logical_and(PLlist[4] <maxlon, PLlist[4] > minlon)]

#chose a random PL in this year
PLindex= randint(0, PLlist.shape[1])
print('PLindex', PLindex)

month= int(PLlist[1, PLindex])
t= int(PLlist[3, PLindex])
PLnr= int(PLlist[2, PLindex])


#reduce the PLlist, such that the automatic plotting routines are working as if only the monthly list was imported
TPL.PLlist= TPL.PLlist[2:, TPL.PLlist[1]== month]

print(year, month, (t//4)+1,  t%4*6)
print('PLnr', PLnr)
print('vort', PLlist[-5, PLindex], 'stab', PLlist[-4, PLindex], 'windspeed', PLlist[-3, PLindex])


d= data('SST', year, month, tstart= t, tend= t+1)          


"""vorticity truncated data"""
plt.figure(fignr)
fignr += 1
plt.clf()

#
#map= Lat_Lon_map(lllon= minlon, lllat= minlat, urlon= maxlon, urlat= maxlat)
#
map= Polar_map(latsouth= d.lat[-1])
grid= np.meshgrid(d.lon, d.lat) #pcolormesh needs corner points for projection
Lon, Lat= map(grid[0], grid[1])


ntrunc, ntrunc2 = 100, 40
vorttrunc_dir=Mediadir+'ERA/pickle/Vort_trunc/Vort_trunc'
vortfilter_comb= pickle.load(open(vorttrunc_dir+str(year)+'_'+str(month).zfill(2)+'_T'+str(ntrunc)+'-'+str(ntrunc2), 'rb'), encoding='latin1').astype('float')

PlotVort(Lon, Lat, vortfilter_comb[t*2], map, maxlevel= 12)
plt.title('truncated vorticity')

PlotTRACK_PL(TPL, t, map)



"""untruncated vorticity data"""
plt.figure(fignr)
fignr += 1
plt.clf()

d= data(['Vort'], year, month, tstart= t, tend= t+1)

map= Polar_map(latsouth= d.lat[-1])
PlotVort(Lon, Lat, d.vort[0], map, boxnr=15)#, maxlevel= 12)
plt.title('vorticity')

PlotTRACK_PL(TPL, t, map)



"""medium cloud"""
plt.figure(fignr)
fignr += 1
plt.clf()

d= data(['MedCloud'], year, month, tstart= t, tend= t+1)

map= Polar_map(latsouth= d.lat[-1])
PlotColorMap3(Lon, Lat, d.MedCloud[0], map, symetric=False, color='blue') #, boxnr=15)#, maxlevel= 12)
plt.title('Medium Cloud')

PlotTRACK_PL(TPL, t, map)



"""wind and surface pressure"""
plt.figure(fignr)
fignr +=1
plt.clf()
#plt.subplot(2, 3, 2)

d= data(['MSLP', 'SST', 'T'], year, month, tstart= t, tend= t+1)

d.imp_u10()
U10= np.sqrt(d.u10[0]**2+d.v10[0]**2)

map= Polar_map(latsouth= d.lat[-1])
PlotWindVelo(Lon, Lat, U10, map, Umax= 25)
PlotWind(d.lon, d.lat, d.u10[0], d.v10[0], map)
PlotContours(Lon, Lat, d.MSLP[0], map)

plt.title('surface wind and pressure')
PlotTRACK_PL(TPL, t, map)



"""static stability (SST-T500)"""    
plt.figure(fignr)
fignr +=1
plt.clf()
#plt.subplot(2, 3, 5)

map= Polar_map(latsouth= d.lat[-1])
#d= data(['SST', 'T'], year, month, tstart= t, tend= t+1)

RCp= 2/7
theta500= d.T[0,2]*(1000/500)**RCp
thetaSST= d.SST[0]*(1000/d.MSLP[0])**RCp


PlotColorMap3(Lon, Lat, d.T[0,2], map, symetric=False, bounds= np.arange(220, 285, 5), label= '$T_{500}$')
PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')


d.impvar('uPL')
d.impvar('vPL')

uPLavg= np.mean(d.uPL[0], axis=0)
vPLavg= np.mean(d.vPL[0], axis=0)

#PlotWindVelo(Lon, Lat, np.sqrt(uPLavg**2+ vPLavg**2), map, Umax= 40)
PlotWind(d.lon, d.lat, uPLavg, vPLavg, map)


plt.title('static stability and mean troposphere wind ')
PlotTRACK_PL(TPL, t, map)




"""tropopause wind"""
plt.figure(fignr)
fignr +=1
plt.clf()
#plt.subplot(2, 3, 3)

d= data(['uPV', 'vPV'], year, month, tstart= t, tend= t+1)

PlotWindVelo(Lon, Lat, np.sqrt(d.uPV[0]**2+ d.vPV[0]**2), map, Umax= 80)
PlotWind(d.lon, d.lat, d.uPV[0], d.vPV[0], map)

map= Polar_map(latsouth= d.lat[-1])


plt.title('tropopause wind ')
PlotTRACK_PL(TPL, t, map)


"""total column water"""
#plt.figure(fignr)
#fignr += 1
#plt.clf()
#
#map= Polar_map(latsouth= d.lat[-1])
#
#d= data(['Water'], year, month, tstart= t, tend= t+1)
#
#
#map= Polar_map(latsouth= d.lat[-1])
#PlotColorMap3(Lon, Lat, d.Water[0], map, symetric=False, color='blue') #, boxnr=15)#, maxlevel= 12)
#plt.title('Water')
#
#PlotTRACK_PL(TPL, t, map)

"""total cloud"""
#plt.figure(fignr)
#fignr += 1
#plt.clf()
#
#map= Polar_map(latsouth= d.lat[-1])
#
#d= data(['TotCloud'], year, month, tstart= t, tend= t+1)
#
#
#map= Polar_map(latsouth= d.lat[-1])
#PlotColorMap3(Lon, Lat, d.TotCloud[0], map, symetric=False, color='blue') #, boxnr=15)#, maxlevel= 12)
#plt.title('Total Cloud')
#
#PlotTRACK_PL(TPL, t, map)
#
#
"""High cloud"""
#plt.figure(fignr)
#fignr += 1
#plt.clf()
#
#map= Polar_map(latsouth= d.lat[-1])
#
#d= data(['HighCloud'], year, month, tstart= t, tend= t+1)
#
#
#map= Polar_map(latsouth= d.lat[-1])
#PlotColorMap3(Lon, Lat, d.HighCloud[0], map, symetric=False, color='blue') #, boxnr=15)#, maxlevel= 12)
#plt.title('High Cloud')
#
#PlotTRACK_PL(TPL, t, map)


"""Low cloud"""
#plt.figure(fignr)
#fignr += 1
#plt.clf()
#
#d= data(['LowCloud'], year, month, tstart= t, tend= t+1)
#
#map= Polar_map(latsouth= d.lat[-1])
#PlotColorMap3(Lon, Lat, d.LowCloud[0], map, symetric=False, color='blue') #, boxnr=15)#, maxlevel= 12)
#plt.title('Low Cloud')
#
#PlotTRACK_PL(TPL, t, map)


"""CAPE"""
#plt.figure(fignr)
#fignr += 1
#plt.clf()
#
#d= data(['CAPE'], year, month, tstart= t, tend= t+1)
#
#map= Polar_map(latsouth= d.lat[-1])
#PlotColorMap3(Lon, Lat, d.CAPE[0], map, symetric=False, color='blue') #, boxnr=15)#, maxlevel= 12)
#plt.title('CAPE')
#
#PlotTRACK_PL(TPL, t, map)