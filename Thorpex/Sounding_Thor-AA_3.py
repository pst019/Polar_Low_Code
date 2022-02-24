#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  5 14:52:50 2018

@author: pst019

conda install -c conda-forge metpy

https://github.com/Unidata/MetPy/tree/master/staticdata
"""

import os
user = os.getcwd().split('/')[2]

import sys
sys.path.insert(0, '../Functions')

from f_plot_fields import * #plot fields (onto basemap or cross section)
#from f_mapplot import * #make basemap plots
from f_imp_thorpex import data as Tdata #import the thorpex data
from f_imp_AROME import data as Adata   
from scipy.interpolate import interp1d

from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from metpy.cbook import get_test_data
import metpy.calc as mpcalc
from metpy.plots import SkewT, Hodograph
from metpy.units import units, concatenate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from metpy.calc.thermo import *

Mediadir= '/media/'+user+'/PatsOrange/'

sound = 3
year, month = 2008, 3

if sound == 1:
    day, hour= 3, 12  
    drops= [x for x in range(1,20) if x not in [3, 5]]  #for the first flight
if sound == 2:
    day, hour= 3, 17  
    drops= [x for x in range(1, 10) if x not in [1, 9]]
#    excl=[1,9,11] #for the second flight
if sound == 3:
    day, hour= 4, 12 
    drops= [x for x in range(1, 20) if x not in [1,5, 13]]
#    excl= [1, 5, 13]# for the third flight

savedir= '/home/pst019/home/Polar_Low/AromeArctic-vs-WRF/2008_03_03/Soundings4/Sound'+str(sound)+'/'

drops= [10]

for dropid in drops:  
    fignr= 3

    fig = plt.figure(fignr, figsize=(6, 5))
    plt.rcParams.update({'font.size': 15})

    fignr += 1
    plt.clf()
    
#    gs = plt.GridSpec(1, 2, width_ratios=[3, 1]) 
#    plt.subplot(gs[0])

    skew = SkewT(fig= fig, rotation=45)
    plt.title('Sounding Dropsonde '+str(dropid))
    
    
    
    thor= Tdata(year, month, day, hour, level='one_profile', dropid= dropid)
    p= thor.pres* units.hPa
    T= thor.T* units.degC
    Td= thor.Td * units.degC
    
    skew.plot(p, T, 'r', label='Thor T')
    skew.plot(p, Td, 'g', label= 'Thor Td')


#    skew.plot_barbs(plevels* units.hPa, interpu* units.m/units.s, interpv *units.m/units.s)
    skew.ax.set_ylim(1000, 400)
    skew.ax.set_xlim(-35, 2)
    
    # Calculate LCL height and plot as black dot
    l = mpcalc.lcl(p[0], T[0], Td[0]) #the lifting condensation level
    #lcl_temp = mpcalc.dry_lapse(concatenate((p[0], l[0])), T[0])[-1].to('degC') #this is l[1]
    #skew.plot(l[0], l[1], 'ko', markerfacecolor='black')
    # Calculate full parcel profile and add to plot as black line
    #prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')
    #skew.plot(p, prof, 'k', linewidth=2)
    
    # Example of coloring area between profiles
    #skew.ax.fill_betweenx(p, T, prof, where=T>=prof, facecolor='blue', alpha=0.4) #CIN
    #skew.ax.fill_betweenx(p, T, prof, where=T<prof, facecolor='red', alpha=0.4) #CAPE
    # An example of a slanted line at constant T -- in this case the 0
    # isotherm
    #l = skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)
    # Add the relevant special lines
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()
    
        
    """get the AA data and also plot it"""
    hour= thor.datetime[0].hour #the hour where the dropsonde reached the ground
    
    if sound in [1,2]:
        exp_name= '080303_cold_pseudo2'
        fileday, filehour= 3, 0 
    if sound== 3:
        exp_name= '080304_cold_pseudo'
        fileday, filehour= 4, 0        

    #exp_name= '080303_warmsens_noTH'
    
    t= (day- fileday)*24 + (hour- filehour)  # -1 
    AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
    AA= Adata(filename= AAfilename, res=1)
    
    dist= (AA.lat - thor.lat[0])**2 + (np.cos(np.deg2rad(thor.lat[0])) * (AA.lon - thor.lon[0])**2)
    x,y= np.where(dist == np.min(dist))
    
#    AA.imp_cross_sec_full(xn= x[0], yn= y[0], tn= t)
    AA.imp_cross_sec_grib(xn= x[0], yn= y[0], tn= t)
   
    
    skew.plot(AA.pres*units.hPa, (AA.T*units.K).to('degC'), 'orange', label= 'Arome T')
    
    AARH= relative_humidity_from_specific_humidity(AA.SH, AA.T* units.K, AA.pres*units.hPa)
    AATd= dewpoint_rh((AA.T*units.K).to('degC'), np.array(AARH) )
    skew.plot(AA.pres*units.hPa, AATd, 'b', label= 'Arome Td')
    
    plt.legend()
    plt.ylabel('Altitude [hPa]')
    plt.xlabel(r"Temperature [$^{\circ}$C]")
    plt.tight_layout()

    
#    plt.savefig(savedir+ 'Sounding_Dropsonde_'+str(dropid)+'_d'+str(day)+'_h'+str(hour))
#    
#    
    """profiles"""
    plt.figure(fignr, figsize=(6, 5))
    fignr +=1
    plt.clf()
    plt.rcParams.update({'font.size': 14})

    ax1= plt.subplot(131)
    
    #plot moist adiabats
    plevels= np.arange(1010, 430, -1)
    moistadiabat= mpcalc.potential_temperature(plevels*units.hPa, mpcalc.moist_lapse(plevels*units.hPa, np.arange(260,289,3)*units.K))
    for i in range(moistadiabat.shape[0]): plt.plot(moistadiabat[i], plevels, color='grey', linewidth=1, linestyle='dashed')

    
    plt.plot(thor.theta, thor.pres, 'r', label='Drop')
    
    AAtheta= potential_temperature(AA.pres*units.hPa, AA.T*units.K)   
    plt.plot(AAtheta, AA.pres, 'b', label='AA')  
    
    plt.xlim([261,285])
    plt.ylim([AA.pres[0], 450])

    plt.ylabel('Altitude [hPa]')
    plt.xlabel(r"Pot. temp. [K]")
#    plt.title('Pot. temperature')
#    plt.legend(loc='upper left')
    
    
    ax2= plt.subplot(132, sharey= ax1)
    plt.setp(ax2.get_yticklabels(), visible=False)

    plt.plot(thor.RH, thor.pres, 'r', label='Drop')   
    plt.plot(AARH, AA.pres, 'b', label='AA')
    plt.xlabel('Rel. hum.')
    
    plt.xlim([-0.05,1.05])
    plt.legend(loc='upper left', fontsize= 12)

    
    ax4= plt.subplot(133, sharey= ax1)
    plt.setp(ax4.get_yticklabels(), visible=False)
    
    plt.plot(thor.U, thor.pres, 'r', label='Drop', lw=0.5)   
    plt.plot(np.sqrt(AA.u**2 + AA.v**2), AA.pres, 'b', label='AA', lw= 0.5)
    plt.xlabel('Hor. vel. [m/s]')


    """the wind barbs"""    
#    ax4= plt.subplot(143, sharey= ax1)
#    plt.setp(ax4.get_yticklabels(), visible=False)

    plevels0= np.arange(np.min([np.max(thor.pres), np.max(AA.pres)] ) , np.min(thor.pres), -50)
    thorinterpu = interp1d(thor.pres, thor.u)(plevels0)
    thorinterpv = interp1d(thor.pres, thor.v)(plevels0)
    thorinterpU= np.sqrt(thorinterpu**2+ thorinterpv**2)
    
    plt.quiver(thorinterpU, plevels0, thorinterpu/thorinterpU, thorinterpv/thorinterpU, scale= 7, width= 0.02, color='r', pivot='mid')
#    plt.xticks([0,1], ['Thor', 'Arome']) 
    
#    plt.xlim(-.5, 1.5)
    
    """wind for AA"""
    del_lon= (AA.lon[x[0]+1, y[0]]- AA.lon[x[0], y[0]])* np.cos(np.deg2rad(AA.lat[x[0], y[0]]))
    del_lat= AA.lat[x[0]+1, y[0]]- AA.lat[x[0], y[0]]
    
    rot= np.rad2deg(np.arctan2(del_lon, del_lat)) #how much to rotate the gridcell anticlockwise to be oriented N-S
#    PlotColorMap4(Lon[1:], Lat[1:], rot, map)
    
    dir_raw= np.rad2deg(np.arctan2(AA.u, AA.v)) #the wind direction in respect to the grid (0 being upwards, 90 rightwards, -90 leftw)
#    PlotColorMap4(Lon[1:], Lat[1:], dir_raw, map)

    dir_rot= (dir_raw + rot)%360 #the wind direction with respect to the North pole (0 Northwards, 90 E, 270 W)
#    PlotColorMap4(Lon[1:], Lat[1:], dir_rot, map)

    wind_vel= np.sqrt(AA.u**2 +AA.v**2)
    AAurot= wind_vel* np.sin(np.deg2rad(dir_rot))
    AAvrot= wind_vel* np.cos(np.deg2rad(dir_rot))
       
    AAinterpu = interp1d(AA.pres, AAurot)(plevels0)
    AAinterpv = interp1d(AA.pres, AAvrot)(plevels0)
    AAinterpU= np.sqrt(AAinterpu**2+ AAinterpv**2)
    
    Q= plt.quiver(AAinterpU, #np.ones(len(plevels0))*1,
                  plevels0, AAinterpu/AAinterpU, AAinterpv/AAinterpU, scale= 7, width= 0.02, color= 'b', pivot='mid')
#    qk = plt.quiverkey(Q, 0.89, 0.93, 20, '20 m/s', coordinates='figure', labelpos='W')

#    plt.ylim([1000, 400])
#    plt.xlabel('Hor. velocity [m/s]')
    plt.xlim([2,30])


    plt.tight_layout()
    plt.savefig(savedir+ 'Profile_Dropsonde_'+str(dropid)+'_d'+str(day)+'_h'+str(hour))
  
    
    
    """difference profiles"""
#    plt.figure(fignr, figsize=(7, 5))
#    fignr +=1
#    plt.clf()
#    
##    plt.suptitle('AROME - Thorpex')
#
#    ax1= plt.subplot(141)
#    plt.axvline(0, color='grey')
#    
#    plevels= np.arange(int(np.min([np.max(thor.pres), np.max(AA.pres) ]) ), int(np.min(thor.pres)), -1) #-2 to avoid extrapolation
#    interp =interp1d(AA.pres, AAtheta)(plevels) - interp1d(thor.pres, thor.theta * units.K)(plevels)
#    plt.plot(interp, plevels)
#
#    #plt.xlim([260,300])
#    plt.ylabel('Altitude [hPa]')
#    plt.xlabel('Pot. temperature [K]')
#    plt.legend()
#
#    
#    ax2= plt.subplot(142, sharey= ax1)
#    plt.axvline(0, color='grey')    
#    plt.setp(ax2.get_yticklabels(), visible=False)
#    
#    interp =interp1d(AA.pres, AARH)(plevels) - interp1d(thor.pres, thor.RH)(plevels)
#    plt.plot(interp, plevels)
#   
#    plt.xlabel('Rel. humidity')
#    #plt.xlim([0,1])
#    plt.legend()
#
#    
#    ax3= plt.subplot(143, sharey= ax1)
#    plt.axvline(0, color='grey')    
#    plt.setp(ax3.get_yticklabels(), visible=False)
#
#    interp =interp1d(AA.pres, np.sqrt(AA.u**2 + AA.v**2))(plevels) - interp1d(thor.pres, np.sqrt(thor.u**2+ thor.v**2))(plevels)
#    plt.plot(interp, plevels)   
#
#    plt.xlabel('Hor. velocity [m/s]')
#    plt.legend()
#
#
#    ax4= plt.subplot(144, sharey= ax1)
#    plt.setp(ax4.get_yticklabels(), visible=False)
#
#
#    Q= plt.quiver(np.ones(len(plevels0))*0, plevels0, AAinterpu -thorinterpu, AAinterpv -thorinterpv,  width= 0.02,  scale= 50)    
#    qk = plt.quiverkey(Q, 0.89, 0.93, 10, '10 m/s', coordinates='figure', labelpos='W')
#
#
#    plt.xlabel('Hor. vel. vectors')
#    plt.ylim([1000, 400])
#    plt.xticks([0], ['Arome - Thor']) 
#
#
#
#    plt.tight_layout()
#    
#    plt.savefig(savedir+ 'Profile_Diff_Dropsonde_'+str(dropid)+'_d'+str(day)+'_h'+str(hour))
    

    
    
    

