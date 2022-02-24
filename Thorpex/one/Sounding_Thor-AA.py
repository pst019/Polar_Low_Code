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
sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_mapplot import * #make basemap plots
from f_imp_thorpex import data as Tdata #import the thorpex data

from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from metpy.cbook import get_test_data
import metpy.calc as mpcalc
from metpy.plots import SkewT, Hodograph
from metpy.units import units, concatenate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from metpy.calc.thermo import *

Mediadir= '/media/'+user+'/1692A00D929FEF8B/'
savedir= '/home/pst019/home/Polar_Low/Graphs/2008_03_03/Soundings/'

year, month = 2008, 3
day, hour= 3, 12  


dropid=1
#exclude= 3, 5, 6,7( above 500hPa)
for dropid in [x for x in range(8, 9) if x not in [3, 5]]:  #for the first flight
#for dropid in [x for x in range(1, 10) if x not in [1, 9]]: #for the second flight
#for dropid in [x for x in range(19, 20) if x not in [5, 13]]: #for the third flight

    
    fignr= 6

    fig = plt.figure(fignr, figsize=(6, 5))
    fignr += 1
    plt.clf()
    skew = SkewT(fig, rotation=45)
    plt.title('Sounding Dropsonde '+str(dropid))
    
    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dictated by the typical meteorological plot
    
    
    
    
    thor= Tdata(year, month, day, hour, level='one_profile', dropid= dropid)
    p= thor.pres* units.hPa
    T= thor.T* units.degC
    Td= thor.Td * units.degC
    
    skew.plot(p, T, 'r')
    skew.plot(p, Td, 'g')
    
    #skew.plot_barbs(p, u, v)
    skew.ax.set_ylim(1000, 200)
    skew.ax.set_xlim(-40, 10)
    
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
    
    
    thorV= Tdata(year, month, day, hour, level='vertical', plevels= np.array([1000, 950, 925, 900, 850, 800, 700, 600, 500, 400]))
    idrop= int(np.where(thorV.dropnr == dropid)[0])
    
    p= thorV.pres[idrop]* units.hPa
    T= thorV.T[idrop]* units.degC
    Td= dewpoint_rh(T, thorV.RH[idrop])
    
    
    skew.plot(p, T, 'ro', label='Thor T')
    skew.plot(p, Td, 'go', label= 'Thor Td')
    
    
    """get the AA data and also plot it"""
    hour= thor.datetime[0].hour #the hour where the dropsonde reached the ground
    
    from f_imp_AROME import data as Adata
#    exp_name= '080303_warmctr'
    exp_name= '080303_cold_pseudo2'
#    exp_name= '080303_warmsens_noTH'
    fileday, filehour= 3, 0        
    t= (day- fileday)*24 + (hour- filehour)  # -1
    t= 0
    print('change back!!!')
    AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract_12.nc'
    AA= Adata(filename= AAfilename, res=1)
    
    dist= (AA.lat - thor.lat[0])**2 + (np.cos(np.deg2rad(thor.lat[0])) * (AA.lon - thor.lon[0])**2)
    x,y= np.where(dist == np.min(dist))
    
    AA.imp_cross_sec_full(xn= x[0], yn= y[0], tn= t)
    
    #offset= -5
    #AA2=  Adata(filename= AAfilename, res=5)
    #AA2.imp_cross_sec(xn= x[0]+offset, yn= y[0]+offset, tn= t)
    
    
    skew.plot(AA.pres*units.hPa, (AA.T*units.K).to('degC'), 'orange', marker='o', label= 'AA T')
    
    Td= dewpoint_rh((AA.T*units.K).to('degC'), AA.RH)
    skew.plot(AA.pres*units.hPa, Td, 'b', marker='o', label= 'AA Td')
    
    plt.legend()
#    plt.savefig(savedir+ 'Sounding_Dropsonde_'+str(dropid)+'_d'+str(day)+'_h'+str(hour))
    
    
    """profiles"""
    plt.figure(fignr)
    fignr +=1
    plt.clf()
    
    ax1= plt.subplot(131)
    plt.plot(thor.theta, thor.pres, 'r', label='Thor')
    plt.plot(thorV.theta[idrop], thorV.pres[idrop], 'ro')
    
    AAtheta= potential_temperature(AA.pres*units.hPa, AA.T*units.K)
    
    plt.plot(AAtheta, AA.pres, 'orange', marker= 'o', ls= 'None', label='AA')
    
    #theta= potential_temperature(AA2.pres*units.hPa, AA2.T*units.K)
    #plt.plot(theta, AA2.pres, 'green', marker= 'o', ls= 'None', label='AA_2')
    
    #plt.ylim([1000, 300])
    plt.xlim([260,300])
    plt.ylabel('Altitude [hPa]')
    plt.title('Potential temperature')
    plt.legend()
    
    ax2= plt.subplot(132, sharey= ax1)
    plt.plot(thor.RH, thor.pres, 'g', label='Thor')
    plt.plot(thorV.RH[idrop], thorV.pres[idrop], 'go')
    
    plt.setp(ax2.get_yticklabels(), visible=False)
    #theta= potential_temperature(AA.pres*units.hPa, AA.T*units.K)
    
    plt.plot(AA.RH, AA.pres, 'bo', label='AA')
    plt.title('Relative humidity')
    
    plt.ylim([1000, 300])
    plt.xlim([0,1])
    plt.legend()
    
    ax3= plt.subplot(133, sharey= ax1)
    plt.plot(thor.U, thor.pres, 'g', label='Thor')
    plt.plot(thorV.U[idrop], thorV.pres[idrop], 'go')
    
    plt.setp(ax3.get_yticklabels(), visible=False)
    #theta= potential_temperature(AA.pres*units.hPa, AA.T*units.K)
    
    plt.plot(np.sqrt(AA.u**2 + AA.v**2), AA.pres, 'bo', label='AA')
    plt.title('Horizontal velocity')
    
    plt.ylim([1000, 300])
    #plt.xlim([0,1])
    plt.legend()
    plt.tight_layout()
#    plt.savefig(savedir+ 'Profile_Dropsonde_'+str(dropid)+'_d'+str(day)+'_h'+str(hour))
    
    
    
    """difference profiles"""
    plt.figure(fignr)
    fignr +=1
    plt.clf()
    
    
    iend = len(thorV.pres[idrop])
    
    ax1= plt.subplot(131)
    #plt.plot(AAtheta[:idend] - thor.theta, thor.pres, 'r', label='hor')
    plt.plot(AAtheta[:iend] -thorV.theta[idrop]* units.K, thorV.pres[idrop], 'ro')
    
    #plt.plot(theta, AA.pres, 'orange', marker= 'o', ls= 'None', label='AA')
    
    #theta= potential_temperature(AA2.pres*units.hPa, AA2.T*units.K)
    #plt.plot(theta, AA2.pres, 'green', marker= 'o', ls= 'None', label='AA_2')
    
    #plt.ylim([1000, 300])
    #plt.xlim([260,300])
    plt.ylabel('Altitude [hPa]')
    plt.title('Potential temperature')
    plt.legend()
    
    ax2= plt.subplot(132, sharey= ax1)
    #plt.plot(thor.RH, thor.pres, 'g', label='Thor')
    plt.plot(AA.RH[:iend] - thorV.RH[idrop], thorV.pres[idrop], 'go')
    
    plt.setp(ax2.get_yticklabels(), visible=False)
    #theta= potential_temperature(AA.pres*units.hPa, AA.T*units.K)
    
    #plt.plot(AA.RH, AA.pres, 'bo', label='AA')
    plt.title('Relative humidity')
    
    plt.ylim([1000, 400])
    #plt.xlim([0,1])
    plt.legend()
    
    ax3= plt.subplot(133, sharey= ax1)
    #plt.plot(thor.U, thor.pres, 'g', label='Thor')
    plt.plot(np.sqrt(AA.u**2 + AA.v**2)[:iend] - thorV.U[idrop], thorV.pres[idrop], 'go')
    
    plt.setp(ax3.get_yticklabels(), visible=False)
    #theta= potential_temperature(AA.pres*units.hPa, AA.T*units.K)
    
    #plt.plot(np.sqrt(AA.u**2 + AA.v**2), AA.pres, 'bo', label='AA')
    plt.title('Horizontal velocity')
    
    plt.ylim([1000, 300])
    #plt.xlim([0,1])
    plt.legend()
    
    plt.suptitle('AROME - Thorpex')
#    plt.savefig(savedir+ 'Profile_Diff_Dropsonde_'+str(dropid)+'_d'+str(day)+'_h'+str(hour))
    
    #plt.tight_layout()
    
    
    
    
    
    print('all fluxes upward')
#    print('Sensible heat flux: '+str(int(AA.SH)) +' W/m**2 \n Latent heat flux: ' +str(int(AA.LH)) +
#          ' W/m**2 \n Precipitation: '+str(round(AA.prec, 3)) + ' mm/h')
    
    #print('LW surf: ' + str(int(AA.LWsurf))+' W/m**2')
    #print('SW surf: ' + str(int(AA.SWsurf))+' W/m**2')
    #print('LW TOA: ' + str(int(AA.LWTOA))+' W/m**2')
    #print('SW TOA: ' + str(int(AA.SWTOA))+' W/m**2')
    
#    print('LW net atm: ' + str(int(AA.LWsurf - AA.LWTOA))+' W/m**2')
#    print('SW net atm: ' + str(int(AA.SWsurf - AA.SWTOA))+' W/m**2')
#    print('net energy atm: '+str(int(AA.SH+ AA.LH+ AA.LWsurf - AA.LWTOA + AA.SWsurf - AA.SWTOA))+' W/m**2')



#plt.figure(fignr)
#fignr +=1
#plt.clf()
#
#plt.plot(AA.PV*1E5, AA.pres)
#plt.ylim([1000, 50])


"""make the hodograph"""
#ax_hod = inset_axes(skew.ax, '40%', '40%', loc=1)
#h = Hodograph(ax_hod, component_range=80.)
#h.add_grid(increment=20)
#h.plot_colormapped(u, v, np.hypot(u, v))

#plt.show()

