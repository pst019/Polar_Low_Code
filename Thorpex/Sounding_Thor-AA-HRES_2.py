#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  5 14:52:50 2018

@author: pst019

from Sounding_Thor-AA_3.py
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
from matplotlib import gridspec

plt.rcParams.update({'font.size': 14})

import numpy as np
import xarray as xr
from metpy.cbook import get_test_data
import metpy.calc as mpcalc
from metpy.plots import SkewT, Hodograph
from metpy.units import units, concatenate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from metpy.calc.thermo import *

Mediadir= '/media/'+user+'/PatsOrange/'

sound = 2
year, month = 2008, 3

save=False
savedir= '/home/pst019/home/Polar_Low/AromeArctic-vs-WRF/2008_03_03/Soundings8/Sound'+str(sound)+'/'

ulevel= 350
ulevel= 450


if sound == 1:
    day, hour= 3, 12  
    drops= [x for x in range(1,20) if x not in [3, 5]]  #for the first flight
if sound == 2:
    day, hour= 3, 18  
    drops= [x for x in range(1, 15) if x not in [2, 12]]
#    drops= [x for x in range(1, 11) if x not in [1, 9]]
#    excl=[1,9,11] #for the second flight
if sound == 3:
    day, hour= 4, 12 
#    drops= [x for x in range(1, 20) if x not in [1,5, 13]] #moved south AA domain
    drops= [x for x in range(1, 20) if x not in [1, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]] # for the third flight - for original AA domain




if sound ==1:      
#    exp_name= 'DA_080303_CTR'
#    fileday, filehour= 3, 0 
#    HRESfileday, HRESfilehour=3, 0 

    exp_name= 'DA_080301_cycling'
    fileday, filehour= 3, 12 
    HRESfileday, HRESfilehour=3, 12

if sound ==2:      
#    exp_name= 'DA_080303_CTR'
#    fileday, filehour= 3, 0 
#    HRESfileday, HRESfilehour=3, 0 
    
    exp_name= 'DA_080301_cycling'
    fileday, filehour= 3, 18 
    HRESfileday, HRESfilehour=3, 18
    
if sound== 3:
#    exp_name= 'DA_080303_CTR'
#    fileday, filehour= 3, 0
#    HRESfileday, HRESfilehour= 3, 0         

    exp_name= 'DA_080301_cycling'
    fileday, filehour= 4, 12 
    HRESfileday, HRESfilehour=4, 12

lacktime= day*24 + hour  - (HRESfileday*24 + HRESfilehour)

"""HRES"""
HRES_name="oper_ML_200803"+str(HRESfileday).zfill(2)+'_'+str(HRESfilehour).zfill(2)+'_'+str(lacktime).zfill(2)+"_vertical.nc"
HRES= xr.open_dataset(Mediadir+'ECMWF/HRES/'+HRES_name, decode_times=False)


drops= [5]


for dropid in drops:  
    fignr= 3

    fig = plt.figure(fignr, figsize=(6, 5))

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
    
    skew.plot(p, T, 'black', label='Thor T')
    skew.plot(p, Td, 'g', label= 'Thor Td')


#    skew.plot_barbs(plevels* units.hPa, interpu* units.m/units.s, interpv *units.m/units.s)
    skew.ax.set_ylim(1000, ulevel)
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
    AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
    AA= Adata(filename= AAfilename, res=1)
    
    dist= (AA.lat - thor.lat[0])**2 + (np.cos(np.deg2rad(thor.lat[0])) * (AA.lon - thor.lon[0])**2)
    x,y= np.where(dist == np.min(dist))
    
    hour= thor.datetime[0].hour #the hour where the dropsonde reached the ground   
#    t= (day- fileday)*24 + (hour- filehour)  # -1 
    t= np.argmin(np.abs(thor.datetime[0] - AA.datetime))
    
#    AA.imp_cross_sec_full(xn= x[0], yn= y[0], tn= t)
    if t<0 : t= 0 #otherwise if t is negative, the last-t timepoint is taken
    AA.imp_cross_sec_grib(xn= x[0], yn= y[0], tn= t)
    print('AA t ', t)
    
    skew.plot(AA.pres*units.hPa, (AA.T*units.K).to('degC'), 'orange', label= 'Arome T')
    
    AARH= relative_humidity_from_specific_humidity(AA.SH, AA.T* units.K, AA.pres*units.hPa)
    AATd= dewpoint_rh((AA.T*units.K).to('degC'), np.array(AARH) )
    skew.plot(AA.pres*units.hPa, AATd, 'r', label= 'Arome Td')


    """get the HRES data"""
    HRES['plev']= HRES['hyam'] + HRES['hybm']*np.exp(HRES['lnsp'][0,0])
    HRES['plev'] /=100
    
    HRES_x= HRES.sel(time=HRES.time.values[0], lat= thor.lat[0], lon= thor.lon[0], method='nearest')

    #skew.plot(AA.pres*units.hPa, (HRES_x.T*units.K).to('degC'), 'orange', label= 'Arome T')

    
    plt.legend()
    plt.ylabel('Altitude [hPa]')
    plt.xlabel(r"Temperature [$^{\circ}$C]")
    plt.tight_layout()

    
#    plt.savefig(savedir+ 'Sounding_Dropsonde_'+str(dropid)+'_d'+str(day)+'_h'+str(hour))
#    
#    
    """profiles"""
    fig= plt.figure(fignr, figsize=(6, 5))
    gs = gridspec.GridSpec(1, 3, width_ratios=[3, 3, 2]) 
    fignr +=1
    plt.clf()

    ax1= plt.subplot(gs[0])
    
    #plot moist adiabats
    plevels= np.arange(1010, ulevel, -1)
    moistadiabat= mpcalc.potential_temperature(plevels*units.hPa, mpcalc.moist_lapse(plevels*units.hPa, np.arange(260,289,3)*units.K))
    for i in range(moistadiabat.shape[0]): plt.plot(moistadiabat[i], plevels, color='grey', linewidth=1, linestyle='dashed')

    
    plt.plot(thor.theta, thor.pres, 'black', label='Drop', lw= 2)
    
    AAtheta= potential_temperature(AA.pres*units.hPa, AA.T*units.K)   
    plt.plot(AAtheta, AA.pres, 'r', label='AA', lw= 2)  
    
    HREStheta= potential_temperature(HRES_x.plev.values* units.hPa, HRES_x.t.values*units.K)
    plt.plot(HREStheta, HRES_x.plev, 'g', label='HRES', lw= 2)
    
    plt.xlim([261,285])
    if ulevel < 380:
        plt.xlim([261,300])   

    plt.ylim([AA.pres[0], ulevel])

    plt.ylabel('Altitude [hPa]')
    plt.xlabel(r"Pot. temp. [K]")
#    plt.title('Pot. temperature')
    plt.legend(loc='upper left')
    
    
    ax2= plt.subplot(gs[1], sharey= ax1)
    plt.setp(ax2.get_yticklabels(), visible=False)

    plt.plot(thor.RH, thor.pres, 'black', label='Drop', lw= 2)   
    plt.plot(AARH, AA.pres, 'r', label='AA', lw= 2)
    
    HRESRH= relative_humidity_from_specific_humidity(HRES_x.q.values, HRES_x.t.values* units.K, HRES_x.plev.values*units.hPa)
    plt.plot(HRESRH, HRES_x.plev, 'g', label='HRES', lw= 2)
   
    plt.xlabel('Rel. humidity')
    
    plt.xlim([-0.05,1.05])

    
#    ax3= plt.subplot(gs[2], sharey= ax1)
#    plt.setp(ax3.get_yticklabels(), visible=False)
#    
#    plt.plot(thor.U, thor.pres, 'black', label='Drop', lw=2)   
#    plt.plot(np.sqrt(AA.u**2 + AA.v**2), AA.pres, 'r', label='AA', lw= 2)
#    HRESU= np.sqrt(HRES_x.u**2 + HRES_x.v**2)
#    plt.plot(HRESU, HRES_x.plev, 'g', label='HRES', lw= 2)  
#
#    plt.xlabel('Velocity [m/s]')
#    plt.xlim([0,30])




    """the wind barbs"""  
    ax4= plt.subplot(gs[2], sharey= ax1)
    plt.setp(ax4.get_yticklabels(), visible=False)
#    plt.quiver(HRESU, HRES_x.plev, HRES_x.u/HRESU, HRES_x.v/HRESU, scale= 7, width= 0.02, color='green', pivot='mid')
    
    
#    ax4= plt.subplot(143, sharey= ax1)
#    plt.setp(ax4.get_yticklabels(), visible=False)

    plevels0= np.arange(np.min([np.max(thor.pres), np.max(AA.pres)] ) , np.min(thor.pres), -50)
    thorinterpu = interp1d(thor.pres, thor.u)(plevels0)
    thorinterpv = interp1d(thor.pres, thor.v)(plevels0)
    thorinterpU= np.sqrt(thorinterpu**2+ thorinterpv**2)
    
#    plt.quiver([0]*len(plevels0), plevels0, thorinterpu/thorinterpU, thorinterpv/thorinterpU, scale= 3, width= 0.025, color='black', pivot='mid')
    plt.barbs([0]*len(plevels0), plevels0, thorinterpu, thorinterpv, barbcolor= 'k', length= 6, pivot='middle')

        
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
    
#    Q= plt.quiver([1]*len(plevels0), plevels0, AAinterpu/AAinterpU, AAinterpv/AAinterpU, scale= 3, width= 0.025, color= 'r', pivot='mid')
    plt.barbs([1]*len(plevels0), plevels0, AAinterpu, AAinterpv, barbcolor= 'r', length= 6, pivot='middle')


    HRESinterpu = interp1d(HRES_x.plev, HRES_x.u)(plevels0)
    HRESinterpv = interp1d(HRES_x.plev, HRES_x.v)(plevels0)
    HRESinterpU= np.sqrt(HRESinterpu**2+ HRESinterpv**2)

#    plt.quiver([2]*len(plevels0), plevels0, HRESinterpu/HRESinterpU, HRESinterpv/HRESinterpU, scale= 3, width= 0.025, color= 'g', pivot='mid')
    plt.barbs([2]*len(plevels0), plevels0, HRESinterpu, HRESinterpv, barbcolor= 'g', length= 6, pivot='middle')


    plt.xlim([-0.5,2.5])
    plt.xticks([0,1,2], ['Dr', 'AA', 'HR'])
    plt.xlabel('Wind')

#    qk = plt.quiverkey(Q, 0.89, 0.93, 20, '20 m/s', coordinates='figure', labelpos='W')

#    plt.ylim([1000, ulevel])
#    plt.xlabel('Hor. velocity [m/s]')


    plt.tight_layout()
    if save: plt.savefig(savedir+ 'Profile_Dropsonde_'+str(dropid)+'_d'+str(day)+'_h'+str(hour)+'HRES'+'+'+str(lacktime).zfill(2))
  
    
    
#    """difference profiles"""
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
#    plt.xlabel(r'$\theta$ [K]')
##    plt.legend()
#
#    
#    ax2= plt.subplot(142, sharey= ax1)
#    plt.axvline(0, color='grey')    
#    plt.setp(ax2.get_yticklabels(), visible=False)
#    
#    interp =interp1d(AA.pres, AARH)(plevels) - interp1d(thor.pres, thor.RH)(plevels)
#    plt.plot(interp, plevels)
#   
#    plt.xlabel('RH')
#    plt.xticks([-.3, 0, .3])
##    plt.legend()
#
#    
#    ax3= plt.subplot(143, sharey= ax1)
#    plt.axvline(0, color='grey')    
#    plt.setp(ax3.get_yticklabels(), visible=False)
#
#    interp =interp1d(AA.pres, np.sqrt(AA.u**2 + AA.v**2))(plevels) - interp1d(thor.pres, np.sqrt(thor.u**2+ thor.v**2))(plevels)
#    plt.plot(interp, plevels)   
#
#    plt.xlabel('U [m/s]')
##    plt.legend()
#
#
#    ax4= plt.subplot(144, sharey= ax1)
#    plt.setp(ax4.get_yticklabels(), visible=False)
#
#
#    Q= plt.quiver(np.ones(len(plevels0))*0, plevels0, AAinterpu -thorinterpu, AAinterpv -thorinterpv,  width= 0.02,  scale= 20)    
#    qk = plt.quiverkey(Q, 0.89, 0.93, 3, '3 m/s', coordinates='figure', labelpos='W')
#
#
#    plt.xlabel('U vectors')
#    plt.ylim([1000, ulevel])
#    plt.xticks([0], ['Arome - Thor']) 
#
#
#
#    plt.tight_layout()
#    plt.savefig(savedir+ 'Profile_Diff_Dropsonde_'+str(dropid)+'_d'+str(day)+'_h'+str(hour))
    

    
    
    

