#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 13:25:51 2017

@author: pst019
from Inv_time_TandS_4
"""

user='patricks'
#user='pst019'

Dropboxdir= '/home/'+user+'/Dropbox/'

import pickle
import sys
sys.path.insert(0, code/Derive_PL_clim/ERA/')
from f_imp_ERA2 import *
from f_imp_ASR import *
from f_impLists import *
sys.path.insert(0, code/AROME/')
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)

from random import randint

fignr= 4

"""randomly chose a PL"""
#year= randint(2000, 2002)
year=2000
TPL=TRACK_PLlist_Year(year, name='', model='ASR')


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
day= 1+ t//8
h3= t%8

PLnr= int(PLlist[2, PLindex])

#reduce the PLlist, such that the automatic plotting routines are working as if only the monthly list was imported
TPL.PLlist= TPL.PLlist[2:, TPL.PLlist[1]== month]



print(year, month, day, h3*3)
print('PLnr', PLnr)
print('vort', PLlist[-5, PLindex], 'stab', PLlist[-3, PLindex], 'windspeed', PLlist[-4, PLindex], 'U500north', PLlist[-2, PLindex])


d= dataASR(year, month, day, sh3= h3, eh3= h3, level='surf')          


"""wind and surface pressure"""
plt.figure(fignr)
fignr += 1
plt.clf()
map= ASR_map(res='c')
Lon, Lat= map(d.lon, d.lat)


d.impvar('SLP', level='surf')
d.impvar('U', level='surf')
d.impvar('V', level='surf')


U10= np.sqrt(d.u10[0]**2+d.v10[0]**2)

#PlotWindVelo(Lon, Lat, U10, map, Umax= 25)
PlotColor_setbounds(Lon, Lat, U10, map, bounds= [0, 10, 16.7, 25], col='Blues')

PlotWind(d.lon, d.lat, d.u10[0], d.v10[0], map, nx= 30, ny= 30, rot= False, option='ASR')
PlotContours(Lon, Lat, d.SLP[0], map)
#
plt.title('surface wind and pressure')
PlotTRACK_PL(TPL, t, map)



"""static stability (SST-T500)"""    
plt.figure(fignr)
fignr +=1
plt.clf()
#plt.subplot(2, 3, 5)

map= ASR_map(res='c')

d.impvar('SST', level='surf')
d.SST[d.SST<272]= np.nan

d.impvar('T', level=500)

RCp= 2/7 #R/Cp
theta500= d.T500[0]*(1000/500)**RCp
thetaSST= d.SST[0]*(1000/d.SLP[0])**RCp


#PlotColorMap3(Lon, Lat, d.T500[0], map, symetric=False, bounds= np.arange(220, 285, 5), label= '$T_{500}$')
#PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')

PlotColor_setbounds(Lon,Lat, (thetaSST- theta500), map, bounds=[-50, -12.5, -9, -6, 0], col='Blues')

#-12.5

d.impvar('U', level=850)
d.impvar('V', level=850)

PlotWind(d.lon, d.lat, d.u850[0], d.v850[0], map, nx= 30, ny= 30, rot= False, option='ASR')

plt.title('static stability and 850 hPa wind ')
PlotTRACK_PL(TPL, t, map)




"""tropopause wind"""
plt.figure(fignr)
fignr +=1
plt.clf()
#plt.subplot(2, 3, 3)

map= ASR_map(res='c')

d.impvar('U', level=500)
d.impvar('V', level=500)

#PlotWindVelo(Lon, Lat, np.sqrt(d.u500[0]**2+ d.v500[0]**2), map, Umax= 80)
PlotColor_setbounds(Lon, Lat, np.sqrt(d.u500[0]**2+ d.v500[0]**2), map, bounds= [0, 29.7, 50, 80], col='Blues')

PlotWind(d.lon, d.lat, d.u500[0], d.v500[0], map, nx= 30, ny= 30, rot= False, option='ASR')

plt.title('tropopause wind ')
PlotTRACK_PL(TPL, t, map)


"""check tropopause wind north"""
masknorth= np.logical_and.reduce((d.lat >= PLlist[5,PLindex], d.lon > PLlist[4, PLindex] -1, d.lon < PLlist[4,PLindex] +1))

if len(np.where(masknorth == True)[0])== 0: #if there is no point north of the PL
    U500north = 0
else:
    U500north = np.max(np.sqrt(d.u500[0][masknorth]**2 + d.v500[0][masknorth]**2))

print('U500north', U500north)

"""check wind"""
latdist2= 5
londist2= int(latdist2//np.cos(np.deg2rad(PLlist[5, PLindex])))
minval= np.min((PLlist[5][PLindex] - d.lat)**2 + (PLlist[4][PLindex] - d.lon)**2)
pts2D= np.where((PLlist[5][PLindex] - d.lat)**2 + (PLlist[4][PLindex] - d.lon)**2 == minval)  

mask2= radMask2D(d.u10.shape[1:], (int(pts2D[0]), int(pts2D[1])), radiusx= latdist2, radiusy= londist2) #lat lon mask
#maskc2= [t, mask2] #tim, lat, lon mask (mask complete)
Ulist= np.max(U10[mask2])
print('U10max', Ulist)


"""compare to ERA"""
TPLera=TRACK_PLlist_2(year, month)


PLeralist= TPLera.PLlist[:, TPLera.PLlist[-1] == 1]

#PLlist= PLlist[:, np.logical_and(PLeralist[4] <maxlon, PLeralist[4] > minlon)]


"""ERA-vorticity truncated data"""
#plt.figure(fignr)
#fignr += 1
#plt.clf()
#
#ntrunc, ntrunc2 = 100, 40
#vorttrunc_dir=Mediadir+'ERA/pickle/Vort_trunc/Vort_trunc'
#vortfilter_comb= pickle.load(open(vorttrunc_dir+str(year)+'_'+str(month).zfill(2)+'_T'+str(ntrunc)+'-'+str(ntrunc2), 'rb'), encoding='latin1').astype('float')
#
#map= Polar_map(latsouth= d.lat[-1])
#PlotVort(Lon, Lat, vortfilter_comb[t], map, maxlevel= 12)
#plt.title('truncated vorticity')
#
#PlotTRACK_PL(TPLera, t//2, map)


"""ERA wind and surface pressure"""
plt.figure(fignr)
fignr +=1
plt.clf()
#plt.subplot(2, 3, 2)

d= data(['MSLP', 'SST', 'T'], year, month, tstart= t//2, tend= t//2+1)

d.imp_u10()
U10= np.sqrt(d.u10[0]**2+d.v10[0]**2)

map= Polar_map(latsouth= d.lat[-1])
grid= np.meshgrid(d.lon, d.lat) #pcolormesh needs corner points for projection
Lon, Lat= map(grid[0], grid[1])

#PlotWindVelo(Lon, Lat, U10, map, Umax= 25)
PlotColor_setbounds(Lon, Lat, U10, map, bounds= [0, 10, 14, 25], col='Blues')

PlotWind(d.lon, d.lat, d.u10[0], d.v10[0], map)
PlotContours(Lon, Lat, d.MSLP[0], map)

plt.title('surface wind and pressure_'+str(t//8+1)+'_'+str(t%8*3))
PlotTRACK_PL(TPLera, t//2, map)


"""ERA static stability (SST-T500)"""    
plt.figure(fignr)
fignr +=1
plt.clf()
#plt.subplot(2, 3, 5)

map= Polar_map(latsouth= d.lat[-1])
#d= data(['SST', 'T'], year, month, tstart= t, tend= t+1)

RCp= 2/7
theta500= d.T[0,2]*(1000/500)**RCp
thetaSST= d.SST[0]*(1000/d.MSLP[0])**RCp


#PlotColorMap3(Lon, Lat, d.T[0,2], map, symetric=False, bounds= np.arange(220, 285, 5), label= '$T_{500}$')
#PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')

PlotColor_setbounds(Lon,Lat, (thetaSST- theta500), map, bounds=[-50, -11.9, -9, -6, 0], col='Blues')


d.impvar('uPL')
d.impvar('vPL')

uPLavg= np.mean(d.uPL[0], axis=0)
vPLavg= np.mean(d.vPL[0], axis=0)

#PlotWindVelo(Lon, Lat, np.sqrt(uPLavg**2+ vPLavg**2), map, Umax= 40)
PlotWind(d.lon, d.lat, uPLavg, vPLavg, map)


plt.title('static stability and mean troposphere wind ')
PlotTRACK_PL(TPLera, t//2, map)


"""era tropopause wind"""
plt.figure(fignr)
fignr +=1
plt.clf()
#plt.subplot(2, 3, 3)

#d= data(['uPV', 'vPV'], year, month, tstart= t, tend= t+1)
d.impvar('uPV')
d.impvar('vPV')

#PlotWindVelo(Lon, Lat, np.sqrt(d.uPV[0]**2+ d.vPV[0]**2), map, Umax= 80)
PlotColor_setbounds(Lon, Lat, np.sqrt(d.uPV[0]**2+ d.vPV[0]**2), map, bounds= [0, 30.5, 50, 80], col='Blues')

PlotWind(d.lon, d.lat, d.uPV[0], d.vPV[0], map)

map= Polar_map(latsouth= d.lat[-1])


plt.title('tropopause wind ')
PlotTRACK_PL(TPLera, t//2, map)
