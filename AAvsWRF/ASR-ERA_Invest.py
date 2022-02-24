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

sys.path.insert(0, code/')
from f_meteo import *

sys.path.insert(0, code/Derive_PL_clim/ERA/')
from f_imp_ERA2 import *
sys.path.insert(0, code/Derive_PL_clim/ASR/')
from f_imp_ASR import *
from f_impLists import *
sys.path.insert(0, code/AROME/')
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)

from random import randint

fignr= 1

frontalevolution =True

"""randomly chose a PL from ASR"""
#year= randint(2000, 2008)
year=2008
month=3
day= 3
h3= 4
t= (day-1)*8 + h3
#t=151
#day= 1+ t//8
#h3= t%8  
    
TPL=TRACK_PLlist_Year(year, name='7', model='ASR')

#just take the PL points
PLlist= TPL.PLlist[:, TPL.PLlist[-1] == 1]
#reduce the PLlist, such that the automatic plotting routines are working as if only the monthly list was imported
TPL.PLlist= TPL.PLlist[2:, TPL.PLlist[1]== month]


print(year, month, day, h3*3)


d= dataASR(year, month, day, sh3= h3, eh3= h3, level='surf')          


"""Lambert coordinates"""
lllon, lllat, urlon, urlat= -15, 63, 60, 75
lat0, lon0= 75, 0 #(lat0, lon0) = center point

"""ASR - wind and surface pressure"""
plt.figure(fignr)
fignr += 1
plt.clf()
#map= ASR_map(res='c')
map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
Lon, Lat= map(d.lon, d.lat)


d.impvar('SLP', level='surf')
d.impvar('U', level='surf')
d.impvar('V', level='surf')


U10= np.sqrt(d.u10[0]**2+d.v10[0]**2)

PlotWindVelo(Lon, Lat, U10, map, Umax= 25)
#PlotColor_setbounds(Lon, Lat, U10, map, bounds= [0, 10, 16.7, 25], col='Blues')

PlotWind(d.lon, d.lat, d.u10[0], d.v10[0], map, nx= 30, ny= 30, rot= False, option='ASR')
PlotContours(Lon, Lat, d.SLP[0], map, nrlevels= 20)
#
#plt.title('ASR surface wind and pressure')
plt.title('ASR')
PlotTRACK_PL(TPL, t, map, track=True, nowPL=False)
#PlotTRACK_PL(TPL, t, map, track=True, nowPL=False)

"""ASR- theta 700 and static stability"""
plt.figure(fignr)
fignr +=1
plt.clf()

d.impvar('T', level=700)
d.impvar('SST', level='surf')
d.SST[d.SST<272]= np.nan #this sets the area with ice to land
d.SST= ma.array(d.SST, mask= isnan(d.SST))
d.impvar('T', level=500)

Theta= PotTemp(d.T700, plev= 700)

map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
PlotContours(Lon, Lat, d.SLP[0], map, nrlevels= 20)
PlotColorMap3(Lon, Lat, Theta[0], map, symetric=False, color='RdBu', bounds= np.arange(264, 290, 1), label= r"Potential temperature at 700hPa [K]")

RCp= 2/7 #R/Cp
theta500= d.T500[0]*(1000/500)**RCp
thetaSST= d.SST[0]*(1000/d.SLP[0])**RCp
PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')


plt.title(r"ASR $\theta_{700}$")
PlotTRACK_PL(TPL, t, map, track=False)


"""ASR- static stability (SST-T500)"""    
#plt.figure(fignr)
#fignr +=1
#plt.clf()
#plt.subplot(2, 3, 5)

#map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
#
#
#d.impvar('SST', level='surf')
#d.SST[d.SST<272]= np.nan #this sets the area with ice to land
#d.SST= ma.array(d.SST, mask= isnan(d.SST))
#
#d.impvar('T', level=500)
#
#RCp= 2/7 #R/Cp
#theta500= d.T500[0]*(1000/500)**RCp
#thetaSST= d.SST[0]*(1000/d.SLP[0])**RCp
#
#
#PlotColorMap3(Lon, Lat, d.T500[0], map, symetric=False, bounds= np.arange(220, 285, 5), label= '$T_{500}$')
#PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')
##PlotColor_setbounds(Lon,Lat, (thetaSST- theta500), map, bounds=[-50, -12.5, -9, -6, 0], col='Blues')
#
#
#d.impvar('U', level=850)
#d.impvar('V', level=850)
#
#PlotWind(d.lon, d.lat, d.u850[0], d.v850[0], map, nx= 30, ny= 30, rot= False, option='ASR')
#
#plt.title('ASR static stability and 850 hPa wind ')
#PlotTRACK_PL(TPL5, t, map, track=False)





"""ASR-SST"""
#plt.figure(fignr)
#fignr +=1
#plt.clf()
#
#
#PlotColorMap3(Lon, Lat, d.SST[0], map, symetric=False)#, bounds= np.arange(220, 285, 5), label= '$T_{500}$')
#PlotTRACK_PL(TPL, t, map)


"""print the SST in area"""
#TPLlistnow= TPL.PLlist[:, np.logical_and(TPL.PLlist[1]== t, TPL.PLlist[-1]==1)]
#
#ni= np.shape(TPLlistnow)[1]
##SSTlist= np.zeros(ni)
#landlist= np.zeros(ni)
#
#for i in range(ni):
#    minval= np.min((TPLlistnow[3][i] - d.lat)**2 + (TPLlistnow[2][i] - d.lon)**2)
#    pts2D= np.where((TPLlistnow[3][i] - d.lat)**2 + (TPLlistnow[2][i] - d.lon)**2 == minval)  
#    
##    t = int(TPLlistnow[1][i]%8)
#
#    latdist2= 4*55/15
#    mask2= radMask2D(d.SST.shape[1:], (int(pts2D[0]), int(pts2D[1])), radiusx= latdist2, radiusy= latdist2) #lat lon mask
#    maskc2= [0, mask2] #tim, lat, lon mask (mask complete)
##    SSTlist[i]= np.mean(d.SST[maskc2])
#    landlist[i]= len(np.where(d.SST[maskc2].mask == True)[0]) /len(d.SST[maskc2])
#
#print(TPLlistnow[0], landlist)




"""ASR- tropopause wind"""
#plt.figure(fignr)
#fignr +=1
#plt.clf()
##plt.subplot(2, 3, 3)
#
#map= ASR_map(res='c')
#
#d.impvar('U', level=500)
#d.impvar('V', level=500)
#
#PlotWindVelo(Lon, Lat, np.sqrt(d.u500[0]**2+ d.v500[0]**2), map, Umax= 80)
##PlotColor_setbounds(Lon, Lat, np.sqrt(d.u500[0]**2+ d.v500[0]**2), map, bounds= [0, 29.7, 50, 80], col='Blues')
#
#PlotWind(d.lon, d.lat, d.u500[0], d.v500[0], map, nx= 30, ny= 30, rot= False, option='ASR')
#
#plt.title('ASR tropopause wind')
#PlotTRACK_PL(TPL, t, map, track=False)



"""compare to ERA"""
TPLera=TRACK_PLlist_2(year, month, name='7')



#"""ERA- untruncated vorticity data"""
#plt.figure(fignr)
#fignr += 1
#plt.clf()
#
#
#de= data(['Vort'], year, month, tstart= t//2, tend= t//2+1)
#
#map= Polar_map(latsouth= de.lat[-1])
#grid= np.meshgrid(de.lon, de.lat) #pcolormesh needs corner points for projection
#Lon, Lat= map(grid[0], grid[1])
#
#PlotVort(Lon, Lat, de.vort[0], map, boxnr=15)#, maxlevel= 12)
#plt.title('ERA vorticity')
#
#PlotTRACK_PL(TPLera, t//2, map, track=True, nowPL=False)

#
"""ERA-vorticity truncated data"""
#plt.figure(fignr)
#fignr += 1
#plt.clf()

#ntrunc, ntrunc2 = 100, 40
#vorttrunc_dir=Mediadir+'ERA/pickle/Vort_trunc/Vort_trunc'
#vortfilter_comb= pickle.load(open(vorttrunc_dir+str(year)+'_'+str(month).zfill(2)+'_T'+str(ntrunc)+'-'+str(ntrunc2), 'rb'), encoding='latin1').astype('float')
#
#map= Polar_map(latsouth= de.lat[-1])
#PlotVort(Lon, Lat, vortfilter_comb[t], map, maxlevel= 12)
#plt.title('ERA truncated vorticity')
#
#PlotTRACK_PL(TPLera, t//2, map, track=False)



"""ERA wind and surface pressure"""
plt.figure(fignr)
fignr +=1
plt.clf()
#plt.subplot(2, 3, 2)

de= data(['MSLP', 'SST', 'T'], year, month, tstart= t//2, tend= t//2+1)
de.SST[de.SST<272]= np.nan #this sets the area with ice to land
de.SST= ma.array(de.SST, mask= isnan(de.SST))

de.imp_u10()
U10= np.sqrt(de.u10[0]**2+de.v10[0]**2)

#map= Polar_map(latsouth= de.lat[-1])
map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
grid= np.meshgrid(de.lon, de.lat) #pcolormesh needs corner points for projection
Lon, Lat= map(grid[0], grid[1])
#
PlotWindVelo(Lon, Lat, U10, map, Umax= 25)
#PlotColor_setbounds(Lon, Lat, U10, map, bounds= [0, 10, 14, 25], col='Blues')

PlotWind(de.lon, de.lat, de.u10[0], de.v10[0], map)
PlotContours(Lon, Lat, de.MSLP[0], map, nrlevels= 20)

#if year >= 2000 and year <= 2012:
#    de.impvar('Geop500')
#    PlotContours(Lon, Lat, de.Geop500[0]/9.81 , map, nrlevels= 10, color='r')
#
#plt.title('ERA surface wind and pressure and geopotential height') #+str(t//8)+'_'+str(t%8*3))
plt.title('ERA')
PlotTRACK_PL(TPLera, t//2, map, nowPL=False, track=True)




"""ERA static stability (SST-T500)"""    
#plt.figure(fignr)
#fignr +=1
#plt.clf()
##plt.subplot(2, 3, 5)
#
#map= Polar_map(latsouth= de.lat[-1])
##de= data(['SST', 'T'], year, month, , tstart= t//2, tend= t//2+1)
#
#RCp= 2/7
#theta500= de.T[0,2]*(1000/500)**RCp
#thetaSST= de.SST[0]*(1000/de.MSLP[0])**RCp
#
#PlotColorMap3(Lon, Lat, de.T[0,2], map, symetric=False, bounds= np.arange(220, 285, 5), label= '$T_{500}$')
#PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')
#
##PlotColor_setbounds(Lon,Lat, (thetaSST- theta500), map, bounds=[-50, -11.9, -9, -6, 0], col='Blues')
#
#
#de.impvar('uPL')
#de.impvar('vPL')
#
#uPLavg= np.mean(de.uPL[0], axis=0)
#vPLavg= np.mean(de.vPL[0], axis=0)
#
##PlotWindVelo(Lon, Lat, np.sqrt(uPLavg**2+ vPLavg**2), map, Umax= 40)
#PlotWind(de.lon, de.lat, uPLavg, vPLavg, map)
#
#
#plt.title('static stability and mean troposphere wind ')
#PlotTRACK_PL(TPLera, t//2, map)

#
"""ERA - tropopause wind"""
plt.figure(fignr)
fignr +=1
plt.clf()
#plt.subplot(2, 3, 3)

map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)


#d= data(['uPV', 'vPV'], year, month, tstart= t//2, tend= t//2+1)
de.impvar('uPV')
de.impvar('vPV')

PlotWindVelo(Lon, Lat, np.sqrt(de.uPV[0]**2+ de.vPV[0]**2), map, Umax= 80)
#PlotColor_setbounds(Lon, Lat, np.sqrt(de.uPV[0]**2+ de.vPV[0]**2), map, bounds= [0, 30.5, 50, 80], col='Blues')

PlotWind(de.lon, de.lat, de.uPV[0], de.vPV[0], map)
if year >= 2000 and year <= 2012:
    de.impvar('Geop500')
    PlotContours(Lon, Lat, de.Geop500[0]/9.81 , map, nrlevels= 10, color='r')

plt.title('ERA tropopause wind ')
PlotTRACK_PL(TPLera, t//2, map, track=False)


"""frontal zones"""
if year >= 2000 and year <= 2012:

    plt.figure(fignr)
    fignr +=1
    plt.clf()
    
    de.impvar('T850')
    de.impvar('sHum850')
    Theta_e= EquiPotTemp(de.T850[0], de.sHum850[0], plev= 850)
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)

    #PlotColorMap3(Lon, Lat, Theta_e, map, symetric=False) #, color='blue') #, boxnr=15)#, maxlevel= 12)
    PlotContours(Lon, Lat, Theta_e, map, nrlevels= 20, color='r')
    PlotContours(Lon, Lat, de.MSLP[0], map, nrlevels= 20)
#    PlotContours(Lon, Lat, de.Geop500[0]/9.81 , map, nrlevels= 10, color='r')
    PlotTRACK_PL(TPLera, t//2, map, track=False)
    plt.title('ERA frontal zones')
    
    
    latdist= 55.
    londist= latdist*np.cos(np.deg2rad(de.lat))
    londist= np.tile(londist, (720,1)).T
    g= np.gradient(Theta_e, latdist, londist) 
    absg= np.sqrt(g[0]**2+ g[1]**2) *100 #change in Theta_e /100km
#    
    PlotColorMap3(Lon, Lat, absg, map, symetric=False, color='blue', bounds= np.linspace(0, 15, 6), label= r"$|\nabla T|$ [K/100km]")
#    PlotColorMap3(Lon, Lat, Theta_e, map, symetric=False, color='RdBu', bounds= np.arange(264, 292, 1), label= r"Equivalent potential temperature [K]")

    PlotTRACK_PL(TPLera, t//2, map, track=False)
    

    """equivalent potential temperature at 850 hPa"""
    plt.figure(fignr)
    fignr +=1
    plt.clf()
    
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)

    #PlotColorMap3(Lon, Lat, Theta_e, map, symetric=False) #, color='blue') #, boxnr=15)#, maxlevel= 12)
#    PlotContours(Lon, Lat, Theta_e, map, nrlevels= 20, color='r')

    PlotContours(Lon, Lat, de.MSLP[0], map, nrlevels= 20)

    PlotTRACK_PL(TPLera, t//2, map, track=False)
    plt.title(r"ERA $\theta_{e,850}$")
    
    PlotColorMap3(Lon, Lat, Theta_e, map, symetric=False, color='RdBu', bounds= np.arange(264, 292, 1), label= r"Equivalent potential temperature at 850hPa [K]")

    PlotTRACK_PL(TPLera, t//2, map, track=False)

    
    """potential temperature at 700 hPa"""
    plt.figure(fignr)
    fignr +=1
    plt.clf()
    
#    de.impvar('T')
#    de.impvar('sHum')
    
#    Theta_e= EquiPotTemp(de.T[1], de.sHum[1], plev= 700)
    Theta= PotTemp(de.T[0, 1], plev= 700)
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)

    #PlotColorMap3(Lon, Lat, Theta_e, map, symetric=False) #, color='blue') #, boxnr=15)#, maxlevel= 12)
#    PlotContours(Lon, Lat, Theta_e, map, nrlevels= 20, color='r')

    PlotContours(Lon, Lat, de.MSLP[0], map, nrlevels= 20)

    PlotTRACK_PL(TPLera, t//2, map, track=False)
    plt.title(r"ERA $\theta_{700}$")
    
    PlotColorMap3(Lon, Lat, Theta, map, symetric=False, color='RdBu', bounds= np.arange(264, 290, 1), label= r"Potential temperature at 700hPa [K]")

    theta500= de.T[0,2]*(1000/500)**RCp
    thetaSST= de.SST[0]*(1000/de.MSLP[0])**RCp
    #
    #PlotColorMap3(Lon, Lat, de.T[0,2], map, symetric=False, bounds= np.arange(220, 285, 5), label= '$T_{500}$')
    PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')

    PlotTRACK_PL(TPLera, t//2, map, track=False)
    
#    if frontalevolution ==True:
#        """frontal zones earlier"""
#        if t>1:
#            tnew = t-2
#            de= data(['MSLP', 'T850', 'sHum850'], year, month, tstart= tnew//2, tend= tnew//2+1)
#        
#            plt.figure(fignr)
#            fignr +=1
#            plt.clf()
#            
#            Theta_e= EquiPotTemp(de.T850[0], de.sHum850[0], plev= 850)
#            
#            map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
#
#            PlotContours(Lon, Lat, Theta_e, map, nrlevels= 20, color='r')
#            PlotContours(Lon, Lat, de.MSLP[0], map, nrlevels= 20)
#    
#            g= np.gradient(Theta_e, latdist, londist) 
#            absg= np.sqrt(g[0]**2+ g[1]**2) *100 #change in Theta_e /100km
#            PlotColorMap3(Lon, Lat, absg, map, symetric=False, color='blue', bounds= np.linspace(0, 15, 6), label= r"$|\nabla T|$ [K/100km]")
#        
#            PlotTRACK_PL(TPLera, t//2, map, othert=tnew//2)#, track=False)
#            plt.title('ERA frontal zones t-6hours')    
#        
#        else: print('first hour of month')
        
#        if t>3:
#            tnew = t-4
#            de= data(['MSLP', 'T850', 'sHum850'], year, month, tstart= tnew//2, tend= tnew//2+1)
#        
#            plt.figure(fignr)
#            fignr +=1
#            plt.clf()
#    
#            
#            Theta_e= EquiPotTemp(de.T850[0], de.sHum850[0], plev= 850)
#            
#            map= Polar_map(latsouth= de.lat[-1])
#            PlotContours(Lon, Lat, Theta_e, map, nrlevels= 20, color='r')
#            PlotContours(Lon, Lat, de.MSLP[0], map, nrlevels= 20)
#    
#            g= np.gradient(Theta_e, latdist, londist) 
#            absg= np.sqrt(g[0]**2+ g[1]**2) *100 #change in Theta_e /100km
#            PlotColorMap3(Lon, Lat, absg, map, symetric=False, color='blue', bounds= np.linspace(0, 15, 6), label= r"$|\nabla T|$ [K/100km]")
#        
#            PlotTRACK_PL(TPLera, t//2, map, othert=tnew//2)#, track=False)
#            plt.title('ERA frontal zones t-12hours')   
#        
#        if t>7:
#            tnew = t-8
#            de= data(['MSLP', 'T850', 'sHum850'], year, month, tstart= tnew//2, tend= tnew//2+1)
#        
#            plt.figure(fignr)
#            fignr +=1
#            plt.clf()
#    
#            
#            Theta_e= EquiPotTemp(de.T850[0], de.sHum850[0], plev= 850)
#            
#            map= Polar_map(latsouth= de.lat[-1])
#            PlotContours(Lon, Lat, Theta_e, map, nrlevels= 20, color='r')
#            PlotContours(Lon, Lat, de.MSLP[0], map, nrlevels= 20)
#    
#            g= np.gradient(Theta_e, latdist, londist) 
#            absg= np.sqrt(g[0]**2+ g[1]**2) *100 #change in Theta_e /100km
#            PlotColorMap3(Lon, Lat, absg, map, symetric=False, color='blue', bounds= np.linspace(0, 15, 6), label= r"$|\nabla T|$ [K/100km]")
#    
#        
#            PlotTRACK_PL(TPLera, t//2, map, othert=tnew//2)#, track=False)
#            plt.title('ERA frontal zones t-24hours')   
#        
#        if t>11:
#            tnew = t-12
#            de= data(['MSLP', 'T850', 'sHum850'], year, month, tstart= tnew//2, tend= tnew//2+1)
#        
#            plt.figure(fignr)
#            fignr +=1
#            plt.clf()
#            
#            de.impvar('uPV')
#            de.impvar('vPV')
#    
#            
#            Theta_e= EquiPotTemp(de.T850[0], de.sHum850[0], plev= 850)
#            
#            map= Polar_map(latsouth= de.lat[-1])
#            PlotContours(Lon, Lat, Theta_e, map, nrlevels= 20, color='r')
#            PlotContours(Lon, Lat, de.MSLP[0], map, nrlevels= 20)
#    
#            g= np.gradient(Theta_e, latdist, londist) 
#            absg= np.sqrt(g[0]**2+ g[1]**2) *100 #change in Theta_e /100km
#            PlotColorMap3(Lon, Lat, absg, map, symetric=False, color='blue', bounds= np.linspace(0, 15, 6), label= r"$|\nabla T|$ [K/100km]")
#    
#            PlotWind(de.lon, de.lat, de.uPV[0], de.vPV[0], map)
#    
#            PlotTRACK_PL(TPLera, t//2, map, othert=tnew//2)#, track=False)
#            plt.title('ERA frontal zones t-36hours')
#        
        
    """ Planetary Boundary Height"""
#    plt.figure(fignr)
#    fignr +=1
#    plt.clf()
#    
#    de.impvar('PBH')
#    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
#    PlotColorMap3(Lon, Lat, de.PBH[0] , map, symetric=False) #, color='blue') #, boxnr=15)#, maxlevel= 12)
#    PlotTRACK_PL(TPLera, t//2, map, track=False)
#    plt.title('ERA planetary boundary height')



"""ERA SST"""
#plt.figure(fignr)
#fignr +=1
#plt.clf()
##plt.subplot(2, 3, 2)
#
#map= Polar_map(latsouth= de.lat[-1])
#
#PlotColorMap3(Lon, Lat, de.SST[0] , map, symetric=False) #, color='blue') #, boxnr=15)#, maxlevel= 12)
##
#plt.title('ERA SST_'+str(t//8+1)+'_'+str(t%8*3))
#PlotTRACK_PL(TPLera, t//2, map, track=False)



"""ERA -medium cloud"""
#if year == 2000:
#    plt.figure(fignr)
#    fignr += 1
#    plt.clf()
#    
#    #d= data(['MedCloud'], year, month, tstart= t//2, tend= t//2+1)
#    de.impvar('MedCloud')
#    
#    map= Polar_map(latsouth= de.lat[-1])
#    PlotColorMap3(Lon, Lat, de.MedCloud[0], map, symetric=False, color='blue') #, boxnr=15)#, maxlevel= 12)
#    plt.title('ERA Medium Cloud')
#    
#    PlotTRACK_PL(TPLera, t//2, map, track=False)
