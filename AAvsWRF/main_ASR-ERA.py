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


"""global variables"""
#PlotPressWind= False
#PlotTheta700= False
#PlotVorticity= True
#PlotTheta_e850= False

#fignr= 1

"""Lambert coordinates"""
#lllon, lllat, urlon, urlat= -15, 63, 60, 75
#lat0, lon0= 75, 0 #(lat0, lon0) = center point

"""time"""
#year=2008
#month=3
#day= 4
#hour= 12

"""c) plot characteristics"""
#number_plevels= 20         
#pot_temp_bounds= np.arange(264, 290, 1)

"""end global variables """


"""module starts"""
h3 = hour// 3
t= (day-1)*8 + h3
   

    
#TPL= TRACK_PLlist_2(year, month, name='7', model='ASR')
    
print(year, month, day, h3*3)


d= dataASR(year, month, day, sh3= h3, eh3= h3)          



"""ASR - wind and surface pressure"""
if PlotPressWind== True:
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
    PlotContours(Lon, Lat, d.SLP[0], map, leveldist= pleveldist)
    #
    #plt.title('ASR surface wind and pressure')
    plt.title('ASR day'+ str(day)+'-hour-'+str(h3*3))
#    PlotTRACK_PL(TPL, t, map, track=True, nowPL=False)
    
    plt.tight_layout()
    
    if savefigs== True:
        print('save')
        plt.savefig(homedir + savedir + 'ASR_wind_test', pad_inches=0)
    
"""ASR- theta 700 and static stability"""
if PlotTheta700== True:
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
    PlotContours(Lon, Lat, d.SLP[0], map, leveldist= pleveldist)
    PlotColorMap3(Lon, Lat, Theta[0], map, symetric=False, color='RdBu', bounds= pot_temp_bounds, label= r"$\theta_{700}$ [K]")
    
    RCp= 2/7 #R/Cp
    theta500= d.T500[0]*(1000/500)**RCp
    thetaSST= d.SST[0]*(1000/d.SLP[0])**RCp
    PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')
    
    
    plt.title(r"ASR $\theta_{700}$")
    plt.tight_layout()

#    PlotTRACK_PL(TPL, t, map, track=False)




"""compare to ERA"""
TPLera=TRACK_PLlist_2(year, month, name='7')
S= STARSList(year, smonth= month, filename='STARS_TRACKS.csv', durationlim= 3, timestep= 'threehours')


"""ERA-vorticity truncated data"""
if PlotVorticity== True:
    plt.figure(fignr)
    fignr += 1
    plt.clf()
    
    ntrunc, ntrunc2 = 100, 40
    vorttrunc_dir=Mediadir+'ERA/pickle/Vort_trunc/Vort_trunc'
    vortfilter_comb= pickle.load(open(vorttrunc_dir+str(year)+'_'+str(month).zfill(2)+'_T'+str(ntrunc)+'-'+str(ntrunc2), 'rb'), encoding='latin1').astype('float')

    de= data(['MSLP'], year, month, tstart= t//2, tend= t//2+1) #just to have de.lat

#    map= Polar_map(latsouth= de.lat[-1])
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
    grid= np.meshgrid(de.lon, de.lat) #pcolormesh needs corner points for projection
    Lon, Lat= map(grid[0], grid[1])
    PlotVort(Lon, Lat, vortfilter_comb[t], map, maxlevel= 12)
    plt.title('ERA truncated vorticity')
    plt.tight_layout()    
#    PlotTRACK_PL(TPLera, t//2, map, track=False)


    plt.figure(fignr)
    fignr += 1
    plt.clf()
    
    de= data(['Vort'], year, month, tstart= t//2, tend= t//2+1) #just to have de.lat

#    map= Polar_map(latsouth= de.lat[-1])
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
    grid= np.meshgrid(de.lon, de.lat) #pcolormesh needs corner points for projection
    Lon, Lat= map(grid[0], grid[1])
    PlotVort(Lon, Lat, de.vort[0], map, maxlevel= 12)
    plt.title('ERA vorticity')
    
#    PlotTRACK_PL(TPLera, t//2, map, track=False, number=False)
    PlotSTARS_PL(S, t//2, map)#, track=False, color= 'b')
    plt.tight_layout()

"""ERA wind and surface pressure"""
if PlotPressWind== True:
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
    PlotContours(Lon, Lat, de.MSLP[0], map, leveldist= pleveldist)
    
    de.impvar('Geop500')
    PlotContours(Lon, Lat, de.Geop500[0]/9.81 , map, nrlevels= 10, color='r')
    #
    #plt.title('ERA surface wind and pressure and geopotential height') #+str(t//8)+'_'+str(t%8*3))
    plt.title('ERA day'+ str(day)+'-hour-'+str(h3*3))
#    PlotTRACK_PL(TPLera, t//2, map, nowPL=False, track=True)
    plt.tight_layout()


"""potential temperature at 700 hPa  + static stability"""
#if PlotTheta700== True:
#    plt.figure(fignr)
#    fignr +=1
#    plt.clf()
#    
#    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
#    
#    PlotContours(Lon, Lat, de.MSLP[0], map, leveldist= pleveldist)
#    
#    Theta= PotTemp(de.T[0, 1], plev= 700)
#    PlotColorMap3(Lon, Lat, Theta, map, symetric=False, color='RdBu', bounds= pot_temp_bounds, label= r"Potential temperature at 700hPa [K]")
#    plt.title(r"ERA $\theta_{700}$")
#    
#    theta500= de.T[0,2]*(1000/500)**RCp
#    thetaSST= de.SST[0]*(1000/de.MSLP[0])**RCp
#    PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')
#    
#    PlotTRACK_PL(TPLera, t//2, map, track=False)
#    plt.tight_layout()

"""equivalent potential temperature at 850 hPa"""
if PlotTheta_e850== True:
    plt.figure(fignr)
    fignr +=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
    
    de.impvar('T850')
    de.impvar('sHum850')
    Theta_e= EquiPotTemp(de.T850[0], de.sHum850[0], plev= 850)
    
    PlotContours(Lon, Lat, de.MSLP[0], map, leveldist= pleveldist)
    
    PlotTRACK_PL(TPLera, t//2, map, track=False)
    plt.title(r"ERA $\theta_{e,850}$")
    
    PlotColorMap3(Lon, Lat, Theta_e, map, symetric=False, color='RdBu', bounds= np.arange(264, 292, 1), label= r"$\theta_{e, 850}$ [K]")
    
#    PlotTRACK_PL(TPLera, t//2, map, track=False)
    plt.tight_layout()
    
#"""ERA - tropopause wind"""
#plt.figure(fignr)
#fignr +=1
#plt.clf()
##plt.subplot(2, 3, 3)
#
#map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
#
#de.impvar('uPV')
#de.impvar('vPV')
#
#PlotWindVelo(Lon, Lat, np.sqrt(de.uPV[0]**2+ de.vPV[0]**2), map, Umax= 80)
##PlotColor_setbounds(Lon, Lat, np.sqrt(de.uPV[0]**2+ de.vPV[0]**2), map, bounds= [0, 30.5, 50, 80], col='Blues')
#
#PlotWind(de.lon, de.lat, de.uPV[0], de.vPV[0], map)
##de.impvar('Geop500')
#PlotContours(Lon, Lat, de.Geop500[0]/9.81 , map, nrlevels= 10, color='r')
#
#plt.title('ERA tropopause wind ')
#PlotTRACK_PL(TPLera, t//2, map, track=False)
#
#
#"""frontal zones"""
#plt.figure(fignr)
#fignr +=1
#plt.clf()
#
#de.impvar('T850')
#de.impvar('sHum850')
#Theta_e= EquiPotTemp(de.T850[0], de.sHum850[0], plev= 850)
#
#map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
#
##PlotColorMap3(Lon, Lat, Theta_e, map, symetric=False) #, color='blue') #, boxnr=15)#, maxlevel= 12)
#PlotContours(Lon, Lat, Theta_e, map, nrlevels= 20, color='r')
#PlotContours(Lon, Lat, de.MSLP[0], map, nrlevels= 20)
##    PlotContours(Lon, Lat, de.Geop500[0]/9.81 , map, nrlevels= 10, color='r')
#PlotTRACK_PL(TPLera, t//2, map, track=False)
#plt.title('ERA frontal zones')
#
#
#latdist= 55.
#londist= latdist*np.cos(np.deg2rad(de.lat))
#londist= np.tile(londist, (720,1)).T
#g= np.gradient(Theta_e, latdist, londist) 
#absg= np.sqrt(g[0]**2+ g[1]**2) *100 #change in Theta_e /100km
##    
#PlotColorMap3(Lon, Lat, absg, map, symetric=False, color='blue', bounds= np.linspace(0, 15, 6), label= r"$|\nabla T|$ [K/100km]")
##    PlotColorMap3(Lon, Lat, Theta_e, map, symetric=False, color='RdBu', bounds= np.arange(264, 292, 1), label= r"Equivalent potential temperature [K]")
#
#PlotTRACK_PL(TPLera, t//2, map, track=False)
#
#




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
