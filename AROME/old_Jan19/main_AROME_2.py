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
Mediadir= '/media/'+user+'/1692A00D929FEF8B/'



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
from f_imp_AROME_exp import *  # Read in netcdf file
from f_meteo import *
from f_useful import *


"""global variables"""
fignr= 3

                  
                   
"""b) times"""
year, month = 2008, 3
day, hour= 4, 12                   


"""Lambert coordinates"""
#lllon, lllat, urlon, urlat= -15, 65, 50, 75
#lat0, lon0= 75, 0 #(lat0, lon0) = center point

#lllon, lllat, urlon, urlat= 0, 67, 60, 75
#lat0, lon0= 70, 0 #(lat0, lon0) = center point

#to match with Thorpex                    
if year == 2008 and month== 3 and day== 3 and hour<15:
    lllon, lllat, urlon, urlat= -3, 70, 18, 75
    lat0, lon0= 75, 0 #(lat0, lon0) = center point

elif year == 2008 and month== 3 and day== 3 and hour>=15:
    lllon, lllat, urlon, urlat= -5, 69.5, 10, 74.5
    lat0, lon0= 75, 0 #(lat0, lon0) = center point

elif year == 2008 and month== 3 and day== 4:
    lllon, lllat, urlon, urlat= -1, 64, 12, 70
    lat0, lon0= 70, 0 #(lat0, lon0) = center point    
    

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
                   
"""c) plot characteristics"""
PlotPressWind= True
PlotPressWind_advanced= False

PlotSLP_laplacian_fitered= False
PlotSLP_fitered= False
PlotMeanSLPminSLP= False

#number_plevels= 20
pleveldist= 1                   

PlotTheta700= False
pot_temp700_bounds= np.arange(270, 282, 1)

PlotTheta_e850= True
pot_tempe850_bounds= np.arange(264, 286, 1)

PlotVorticity= False
PlotThickness= False
PlotGeop500=False
PlotPrecip=False
PlotHumidity=False
PlotTH=False
PlotLH = False
PlotFLX = False
PlotBL=False
PlotW=False
PlotPV=False
PlotTsurf=False
PlotT2m=False
PlotCloud_high=False
PlotCloud_med=False    
PlotCloud_low=False
PlotCAPE=False
PlotCIN=False #not very helpful
PlotLCL=False  #not very helpful


save=False
savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/Oslo/'

"""end global variables"""


"""start module"""
#d= data(filename= Mediadir+'PL/AA/ec/'+str(year)[-2:]+str(month).zfill(2)+str(fileday).zfill(2)+'_test_3_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc')
#d.imp_standard()

AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
AAres=1 #2- every second datapoint is taken
d= data(filename= AAfilename, res= AAres)


"""MSLP, Surface temperature and surface winds"""
if PlotPressWind== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
#    map = AA_map()
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
    Lon,Lat = map(d.lon,d.lat)
    
    d.imp_level(pn= 0, tn= t)
    geop1000 = d.geop
    Temp1000 = d.T
    
    # Draw contours for mslp
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    #
    #PlotSurfTemp(Lon,Lat, d.T1000 - 273.15, m)
    PlotWind(d.lon, d.lat, d.u, d.v, map) #???
    PlotWindVelo(Lon, Lat, np.sqrt(d.u**2+ d.v**2), map, Umax= 25)
    
    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()
    
    if save==True:
        plt.savefig(savedir+'Arome_'+exp_name+'_'+str(t).zfill(2) )

"""MSLP, Surface temperature and surface winds -advanced"""
if PlotPressWind_advanced== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map = AA_map()
#    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
    Lon,Lat = map(d.lon,d.lat)
#    
    d.imp_level(pn= 0, tn= t)
    geop1000 = d.geop
    Temp1000 = d.T
    
    # Draw contours for mslp
    PlotContours(Lon, Lat, d.mslp, map, leveldist= 1, numbers=False)

    PlotWindVelo(Lon, Lat, np.sqrt(d.u**2+ d.v**2), map, Umax= 25)    
    #find local pressure min and plot them in the map   
    U= np.sqrt(d.u**2+d.v**2)
    PlotLocalMax(d.mslp, threshold=1010, distance=100/AAres, map= map, lon=d.lon, lat=d.lat,
                 data2=U, threshold2=16, distance2=80/AAres)    
  
    PlotLocalMax(U, threshold=20, distance=100/AAres, map= map, lon=d.lon, lat=d.lat, typ='max', color='orange', dot=False)    
    
    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()
    
    if save==True:
        plt.savefig(savedir+'Arome_'+exp_name+'_'+str(t).zfill(2) )



"""mean(SLP) - SLP"""
if PlotMeanSLPminSLP== True:
    
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map = AA_map()
    AAres2= 10 #otherwise it becomes really uneffecient
    SLP_radmean= CircularFilter_samedist(d.mslp, radius= 110/(2.5* AAres2))
    
    PlotColorMap4(Lon, Lat, d.mslp - SLP_radmean, map, label=r'$\overline{SLP} - SLP$ in 110km radius')
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist, color= 'red')


"""filtered MSLP"""
if PlotSLP_fitered== True:
    
    MSLPfourierfilter= FourierFilter2d_equaldist(d.mslp, dist= 2.5* AAres, T_low= 40, T_up= 200)

    plt.figure(fignr)
    fignr+=1
    plt.clf()

    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)

    PlotContours(Lon, Lat, MSLPfourierfilter, map, leveldist= pleveldist)
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist, color= 'red')
    plt.tight_layout()
    

if PlotSLP_laplacian_fitered== True:
    SLPlap= np.gradient(np.gradient(d.mslp, dx, axis= 0, edge_order=2), dx, axis= 0, edge_order=2) + np.gradient(np.gradient(d.mslp, dx, axis= 1, edge_order=2), dx, axis= 1, edge_order=2)
    
    import scipy.ndimage.filters as filters
    SLPlapfilter= filters.gaussian_filter(SLPlap, sigma= 16/AAres, mode='constant', truncate= 1.)
    plt.figure(fignr)
    fignr+=1
    plt.clf()

    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)

    PlotColorMap4(Lon, Lat, SLPlapfilter, map)
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist, color= 'red')
    plt.tight_layout()
    

"""for cross section definition"""
#x_start, y_start= 10, 30
#x_end, y_end= 110, 20
#
#x, y= mapcross_section(d.lon, d.lat, x_start, y_start, x_end, y_end, m)


"""potential temperature at 700 hPa and stability"""
if PlotTheta700== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
    
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist*2)
    
    pn= 6 #6= 700hPa
    d.imp_level(pn= pn, tn= t)
    
    theta= PotTemp(d.T, d.press[pn])
    PlotColorMap3(Lon, Lat, theta, map, bounds= pot_temp700_bounds, label=r"$\theta_{700}$ [K]")
    
    plt.title('Arome '+str(d.datetime[t])[:-3] +' static stability')
    
    
    #"""stability"""
    pn= 8 # 8 = 500hPa
    d.imp_level(pn= pn, tn= t)
    
    theta500= PotTemp(d.T, d.press[pn])
    thetaSST= Temp1000  #should change this
    PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')
    plt.tight_layout()

""" equivalent potential temperature at 850hPa"""
if PlotTheta_e850 == True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
    Lon,Lat = map(d.lon,d.lat)

    pn= 4 #850hPa
    d.imp_level(pn= pn, tn= t)
    
#    theta= d.T*(1000/d.press[pn])**(2/7)
    spechum= d.RH/(.263 *d.press[pn])* np.exp(17.67*(d.T-273.15) / (d.T-29.65))
#    theta_e= theta* np.exp(2501 * spechum/(1.006* d.T))
    theta_e= EquiPotTemp(d.T, spechum, d.press[pn])
    
    PlotColorMap3(Lon, Lat, theta_e, map, bounds= pot_tempe850_bounds, label=r"$\theta_{e,850}$ [K]")
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()


"""calculate vorticity"""
if PlotVorticity == True:
       
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    pn= 4 #850hPa
    d.imp_level(pn= pn, tn= t)
    
#    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
    map = AA_map()
    Lon,Lat = map(d.lon,d.lat)

    dx= 2500* AAres
    vort= (np.gradient(d.v, dx, axis= 1)- np.gradient(d.u, dx, axis= 0))*1E5
    #probably one has to rotate u and v first to proper x and y coordinates
    import scipy.ndimage.filters as filters
    
    gausfilterdist= 60
    vortfilter= filters.gaussian_filter(vort, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 1.)
    PlotVort(Lon, Lat, vortfilter, map)
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    plt.title('Arome '+str(d.datetime[t])[:-3] + ' vort with Gauss filter of '+str(gausfilterdist)+'km')

    plt.tight_layout()


PlotVorticityTfiltered= False
    
"""the fourier filtered vorticity"""
if PlotVorticityTfiltered== True:
    vortfourierfilter= FourierFilter2d_equaldist(vort, dist= 2.5* AAres, T_low= 40, T_up= 100)

    plt.figure(fignr)
    fignr+=1
    plt.clf()

    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)

    PlotVort(Lon, Lat, vortfourierfilter, map)
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    plt.title('Arome '+str(d.datetime[t])[:-3] + ' vorticity with filtered T40 - T100')

    plt.tight_layout()

"""PV"""
if PlotPV==True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    plt.title('Arome '+str(d.datetime[t])[:-3])
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
    d.imp_level(pn= 10, tn= t)
    PlotContours(Lon, Lat, d.mslp, map)
    PlotPV(Lon, Lat, d.PV, map)
    plt.title('Arome '+str(d.datetime[t])[:-3] + ' potential vorticity')

    plt.tight_layout()


"""geopotential 1000 -500 hPa thickness"""

if PlotThickness==True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    plt.title('Arome '+str(d.datetime[t])[:-3] +' 1000-500 hPa thickness')
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
    Lon,Lat = map(d.lon,d.lat)

    # 1000 -500 hPa thickness
    pn= 8
    d.imp_level(pn= pn, tn= t)
    geop500= d.geop
    #
    geop_diff= d.geop - geop1000
    PlotContours(Lon, Lat, geop_diff, map)
    PlotColorMap3(Lon, Lat, geop_diff, map, symetric=False, label='1000-500 hPa thickness')
    plt.tight_layout()

""" different geopotential levels and mslp"""
if PlotGeop500==True:
    plt.figure(fignr)
    plt.clf()
    fignr += 1
    
    plt.title('Arome '+str(d.datetime[t])[:-3] +' SLP and 500hPa geopotential height')

    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
    Lon,Lat = map(d.lon,d.lat)
 
    pn= 8 #500hPa
    d.imp_level(pn= pn, tn= t)    
    
  
    
    PlotContours(Lon, Lat, d.mslp, map, leveldist= 2)
    PlotContours(Lon, Lat, d.geop, map, color='r', leveldist= 10)
#    PlotContours(Lonm, Latm, geop500, m, color=[(plt.cm.RdBu_r(h)) for h in range(0,257,32)], leveldist= 20)
    
    #d.imp_level(pn= 6, tn= t)
    #PlotContours(Lonm, Latm, d.geop, m, color='b')
    
    plt.tight_layout()


""" Precip """
if PlotPrecip== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)  
#    d.imp_atm(tn= t)
#    Lon,Lat = map(d.lon,d.lat)
    
    precipbounds= [0, 0.1, 0.2, 0.5, 1, 2, 5, 10]
    # Draw contours for mslp
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    PlotColorMap4(Lon, Lat, d.prec, map, bounds= precipbounds, color='blue', label='Precipitation [mm/h]')

#    PlotPrecip(Lon, Lat, d.prec, map)
   
    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()


""" Sensible Heat Flux"""
if PlotTH== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)  
    d.imp_level(pn= 0, tn= t)
    Lon,Lat = map(d.lon,d.lat)

    # Draw contours for mslp
    fluxbounds= [-800,  -500, -300, -200 , -100, 0, 100, 200, 300, 500, 800]
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    PlotColorMap4(Lon, Lat, d.SH, map, bounds= fluxbounds, label=r"Sensible heat flux [W/m$^2$]")
   
    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()


if PlotLH== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)  
#    d.imp_atm(tn= t)
#    Lon,Lat = map(d.lon,d.lat)

    # Draw contours for mslp
    fluxbounds= [-800,  -500, -300, -200 , -100, 0, 100, 200, 300, 500, 800]
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    PlotColorMap4(Lon, Lat, d.LH, map, bounds= fluxbounds, label=r"Latent heat flux [W/m$^2$]")
   
    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()

if PlotFLX== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)  
#    d.imp_atm(tn= t)
#    Lon,Lat = map(d.lon,d.lat)

    # Draw contours for mslp
    fluxbounds= [0,  100, 200, 300, 500, 800, 1200]
    
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    PlotColorMap4(Lon, Lat, d.SH+ d.LH, map, bounds= fluxbounds, color= 'brown', label=r"Turbulent heat flux [W/m$^2$]")
   
    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()


if PlotBL== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)  
    d.imp_atm(tn= t)
    Lon,Lat = map(d.lon,d.lat)

    # Draw contours for mslp
    BLbounds= [0,  500, 1000, 2000, 3000, 4000, 5000]
    
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    PlotColorMap4(Lon, Lat, d.BL, map, bounds=BLbounds,  color= 'blue', label=r"Boundary layer height [m]")
   
    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()



if PlotW== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)  
    d.imp_level(pn= 6, tn= t)
    Lon,Lat = map(d.lon,d.lat)

    # Draw contours for mslp
    wbounds= [-5, -1, -0.5, -0.2, -0.1, 0,  0.1, 0.2, 0.5, 1, 5]
    
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    PlotColorMap4(Lon, Lat, d.w, map, bounds=wbounds, color= 'RdBu', label=r"Vertical velocity at 700hPa [m/s]")
   
    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()



if PlotPV== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)  
    d.imp_level(pn= 8, tn= t)
    Lon,Lat = map(d.lon,d.lat)

    # Draw contours for mslp
    PVbounds= np.array([-5, -3, -2, -1, -0.5, 0,  0.5, 1, 2, 3, 5])
    
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    PlotColorMap4(Lon, Lat, d.PV *1E6, map, bounds=PVbounds, color= 'RdBu', label=r"Potential vorticity [$10^{-6}$K m$^2$ kg$^{-1}$ s$^{-1}$]")
    
#    print('wrong units')
    
    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()



if PlotTsurf== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)  
    d.imp_surf(tn= t)
    Lon,Lat = map(d.lon,d.lat)

    # Draw contours for mslp
    Tsurfbounds= np.array([-40, -20, -10, -5, -2, 0,  2,  5, 10, 20, 40])
    
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    PlotColorMap4(Lon, Lat, d.T0m -273.15 , map, bounds= Tsurfbounds, color= 'RdBu', label=r"Surface temperature [K]")
        
    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()



if PlotT2m== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)  
    d.imp_surf(tn= t)
    Lon,Lat = map(d.lon,d.lat)

    # Draw contours for mslp
    Tsurfbounds= np.array([-40, -20, -10, -5, -2, 0,  2,  5, 10, 20, 40])
    
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    PlotColorMap4(Lon, Lat, d.T2m -273.15 , map, bounds= Tsurfbounds, color= 'RdBu', label=r"2m temperature [K]")
        
    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()

if PlotHumidity== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)  
    d.imp_level(pn= 6, tn= t)

    # Draw contours for mslp
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    PlotColorMap3(Lon, Lat, d.RH, map, symetric=False, label='Relative humidity', color='blue', bounds= np.array([0, 0.2, 0.4, 0.6, 0.8, 1]) )
   
    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()


#PlotCloud_conv=False
#
#if PlotCloud_conv== True:
#    plt.figure(fignr)
#    fignr+=1
#    plt.clf()
#    
#    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)  
#    d.imp_atm(tn= t)
#    Lon,Lat = map(d.lon,d.lat)
#
#    Cloud_convbounds= np.array([0,  100,  200, 300, 500, 800])
#    
#    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
#    PlotColorMap4(Lon, Lat, d.cloud_conv , map, bounds=Cloud_convbounds, color= 'brown', label=r"CAPE [J/kg]")
#
#    plt.title('Arome '+str(d.datetime[t])[:-3])
#    plt.tight_layout()
#


if PlotCloud_high== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)  
    d.imp_atm(tn= t)
    Lon,Lat = map(d.lon,d.lat)

    Cloud_bounds= np.array([0, 0.2, 0.4, 0.6, 0.8, 1])
    
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    PlotColorMap4(Lon, Lat, d.cloud_high , map, bounds=Cloud_bounds, color= 'grey', label=r"High cloud cover")

    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()



if PlotCloud_med== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)  
    d.imp_atm(tn= t)
    Lon,Lat = map(d.lon,d.lat)

    Cloud_bounds= np.array([0, 0.2, 0.4, 0.6, 0.8, 1])
    
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    PlotColorMap4(Lon, Lat, d.cloud_med , map, bounds=Cloud_bounds, color= 'grey', label=r"Medium cloud cover")

    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()
    


if PlotCloud_low== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)  
    d.imp_atm(tn= t)
    Lon,Lat = map(d.lon,d.lat)

    Cloud_bounds= np.array([0, 0.2, 0.4, 0.6, 0.8, 1])
    
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    PlotColorMap4(Lon, Lat, d.cloud_low , map, bounds=Cloud_bounds, color= 'grey', label=r"Low cloud cover")

    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()



if PlotCAPE== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)  
    d.imp_atm(tn= t)
    Lon,Lat = map(d.lon,d.lat)

    CAPEbounds= np.array([0,  100,  200, 300, 500, 800])
    
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    PlotColorMap4(Lon, Lat, d.CAPE , map, bounds=CAPEbounds, color= 'brown', label=r"CAPE [J/kg]")

        
    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()



if PlotCIN== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)  
    d.imp_atm(tn= t)
    Lon,Lat = map(d.lon,d.lat)

    CINbounds= np.array([-100, -10, -1, -0.1, 0])
    
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    PlotColorMap4(Lon, Lat, d.CIN , map, bounds=CINbounds, color= 'brown', label=r"CIN [J/kg]")

        
    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()
 
    
PlotLFC=False #not very helpful

if PlotLFC== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)  
    d.imp_atm(tn= t)
    Lon,Lat = map(d.lon,d.lat)

    LFCbounds= np.array([0, 10, 100, 200, 500, 1000, 10000, 30000])
    
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    PlotColorMap4(Lon, Lat, d.LFC , map, bounds=LFCbounds, color= 'brown', label=r"LFC [m]")

        
    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()        



if PlotLCL== True:
    plt.figure(fignr)
    fignr+=1
    plt.clf()
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)  
    d.imp_atm(tn= t)
    Lon,Lat = map(d.lon,d.lat)

    LCLbounds= np.array([0, 10, 100, 1000, 3000, 5000, 10000, 30000])
    
    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
    PlotColorMap4(Lon, Lat, d.LCL , map, bounds=LCLbounds, color= 'brown', label=r"LCL [m]")

        
    plt.title('Arome '+str(d.datetime[t])[:-3])
    plt.tight_layout()    
    
    
##"""cross section temperature"""
##plt.figure(fignr)
##fignr+=1
##plt.clf()
##
##plt.title('temperature')
##
##d.imp_cross_sec(xn= x, yn= y, tn= t)
##PlotCross_sec_T(d.press[:-5], d.lon[y,x], d.T[:-5])
##
##
##"""cross section potential temperature"""
##plt.figure(fignr)
##plt.clf()
##fignr += 1
##
##plt.title('potential temperature')
##
##pressm= np.tile(d.press, (np.shape(d.T)[1],1)).T #make a matrix with the pressure in every row
##theta= d.T*(1000/pressm)**(2/7)
##
##PlotCross_sec_T(d.press[:-5], d.lon[y,x], theta[:-5])
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
###spechum= RH[0to1]*100/(26.3 * press[hPa]) * exp(17.67*(T[k]-273.15)/(T[K]-29.65))
##spechum= d.RH/(.263 *pressm)* np.exp(17.67*(d.T-273.15) / (d.T-29.65))
##
###equivalent pot temp theta_e = theta* exp(Lc * spechum/(Cp * T))
##theta_e= theta* np.exp(2501 * spechum/(1.006* d.T))
##
##PlotCross_sec_T(d.press[:-5], d.lon[y,x], theta_e[:-5])
##
##PlotCross_sec_hum(d.press[:-5], d.lon[y,x], d.RH[:-5])
##
##fignr= 5
##"""vertical velocity"""
###plt.figure(fignr)
###plt.clf()
###fignr += 1
###
###PlotCross_sec_vervel(d.press, d.lon[y,x], d.w)
##
##
##fignr= 6
##"""horizontal velocity (north- south)"""
###plt.figure(fignr)
###plt.clf()
###fignr += 1
###
###PCross_sec_horvel(d.press, d.lon[y,x], d.u, d.v)
###
###
###plt.figure(fignr)
###plt.clf()
###fignr += 1
###
###PCross_sec_horvel2(d.press, d.lon[y,x], d.u, d.v)
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