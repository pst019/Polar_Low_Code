#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 13:25:51 2017

@author: pst019
"""

#the data is downloaded from: https://rda.ucar.edu/#!lfd?nb=y&b=plat&v=REANALYSIS%20MODELS

import pickle
import sys
import os
user = os.getcwd().split('/')[2]
#print('user: ', user)

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'



sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')
#from f_imp_ASR import *
from f_impLists import *
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_meteo import *
import xarray as xr

#from random import randint


"""global variables"""
#var = 'PressWind'
#var= 'PressWind_advanced'
#var = 'Theta700'
#var= 'RH_lev'
var= 'SH_lev'


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
hour= 18
pressure= 850 #hPa


"""prepare dataset"""
filedir= Mediadir+'/ASR/asr15km.anl.3D.'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'.nc'

ASR= xr.open_dataset(filedir)
ASR= ASR.sel(Time= str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'T'+str(hour).zfill(2))

index_pres= np.where(ASR['PRES'] == pressure * 100)[0][0]
ASR= ASR.sel(num_metgrid_levels= index_pres)


filedir= Mediadir+'/ASR/asr15km.anl.2D.'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'.nc'

ASRsurf= xr.open_dataset(filedir)
ASRsurf= ASRsurf.sel(Time= str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'T'+str(hour).zfill(2))


"""c) plot characteristics"""

Track= False

"""end global variables """

save= True
savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/Fields/' 
title_extra='jet' #''


Lon, Lat= map(ASR['XLONG'].values, ASR['XLAT'].values)
#
#

if var== 'RH_lev':
    PlotContours(Lon, Lat, ASRsurf['PMSL']/100, map, leveldist= 1, numbers=False)
    new_map, norm= PlotColorMap4(Lon, Lat, ASR['RH']/100, map, symetric=False, label='Relative humidity at ' +str(pressure)+'hPa', color='blue', bounds= np.array([0, 0.4, 0.7, 0.8, 0.9, 0.95, 1]), return_colormap= True )


if var== 'SH_lev':
    PlotContours(Lon, Lat, ASRsurf['PMSL']/100, map, leveldist= 1, numbers=False)


    import metpy.calc as mpcalc
    from metpy.units import units 
    dSH= mpcalc.specific_humidity_from_mixing_ratio(
            mpcalc.mixing_ratio_from_relative_humidity(ASR['RH']/100, ASR['TT'].values* units.K, pressure* units.hPa))


    new_map, norm= PlotColorMap4(Lon, Lat, dSH*1E3 , map, label= 'Specific humidity [g/kg] at '+str(pressure)+' hPa',
                  color='jet_r', bounds= np.arange(0, 2.3, .2))    
    PlotContours(Lon, Lat, ASR['RH']/100, map, levels= [0.9], numbers= False, color='white')

#
#"""ASR - wind and surface pressure"""
#if var== 'PressWind':
#    d.impvar('SLP', level='surf')
#    d.impvar('U', level='surf')
#    d.impvar('V', level='surf')
#    
#    U10= np.sqrt(d.u10[0]**2+d.v10[0]**2)
#    PlotWindVelo(Lon, Lat, U10, map, Umax= 25)
#    
#    PlotWind(d.lon, d.lat, d.u10[0], d.v10[0], map, nx= 30, ny= 30, rot= False, option='ASR')
#    PlotContours(Lon, Lat, d.SLP[0], map, leveldist= 1)
#    
#
#if var== 'PressWind_advanced':
#    d.impvar('SLP', level='surf')
#    d.impvar('U', level='surf')
#    d.impvar('V', level='surf')
#
#    PlotContours(Lon, Lat, d.SLP[0], map, leveldist= 1, numbers=False)
#        
#    U10= np.sqrt(d.u10[0]**2+d.v10[0]**2)
#    PlotWindVelo(Lon, Lat, U10, map, Umax= 25, color='YlBu')
#    
#    PlotLocalMax(d.SLP[0], threshold=1010, distance=30, map= map, lon=d.lon, lat=d.lat,
#                     data2=U10, threshold2=10, distance2= 10) 
#
#    PlotLocalMax(U10, threshold=15, distance=40, map= map, lon=d.lon, lat=d.lat, typ='max',
#                 color='orange', dot=False, roundorder= 0, latbound= [65, 80])  
#    
"""ASR- theta 700 and static stability"""
#if var == 'Theta':
#    d.impvar('SLP', level='surf')   
#    d.impvar('T', level=700)
#    d.impvar('SST', level='surf')
#    d.SST[d.SST<272]= np.nan #this sets the area with ice to land
#    d.SST= np.ma.array(d.SST, mask= np.isnan(d.SST))
#    d.impvar('T', level=500)
#    
#    Theta= PotTemp(d.T700, plev= 700)
#    
#    PlotContours(Lon, Lat, d.SLP[0], map, leveldist= pleveldist)
#    bounds= np.arange(264, 290, 1)
#    PlotColorMap3(Lon, Lat, Theta[0], map, symetric=False, color='RdBu', bounds= bounds, label= r"$\theta_{700}$ [K]")
#    
#    RCp= 2/7 #R/Cp
#    theta500= d.T500[0]*(1000/500)**RCp
#    thetaSST= d.SST[0]*(1000/d.SLP[0])**RCp
#    PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')
#    
#    
#
#if Track== True:
#    PlotTRACK_PL(TPL, t, map, track=True, nowPL=False)
#
#
#
#if save== False:
#    plt.title('ASR '+str(year)+'_'+str(month).zfill(2)+'_'+str(day).zfill(2)+'_'+str(hour).zfill(2)+ title_extra)
#    #datetime is shifted by one day
#
#plt.tight_layout()
#

if save== True:
    if '_lev' in var:
        savevar = var+ str(pressure)
    else: savevar= var    
    savename= savedir + 'ASR_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+str(hour).zfill(2)+'_'+savevar+title_extra
    plt.savefig(savename, pad_inches=0)
    print(savename)