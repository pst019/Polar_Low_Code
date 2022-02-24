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
from f_imp_AROME import data as Adata #import the thorpex data
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

Mediadir= '/media/'+user+'/1692A00D929FEF8B/'
savedir= '/home/pst019/home/Polar_Low/AromeArctic-vs-WRF/2008_03_03/Soundings2/Sound3/'

year, month = 2008, 3
day, hour= 4, 12  

fignr= 9

#exclude= 3, 5, 6,7( above 500hPa)
#excl= [3, 5] #for the first flight
#excl=[1,9,11] #for the second flight
excl= [1, 5, 13]# for the third flight

fig = plt.figure(fignr, figsize=(6, 5))
fignr += 1
plt.clf()
skew = SkewT(fig, rotation=45)
#plt.title('Sounding Dropsonde average')


plevels= np.arange(980, 500, -1) #np.array([1000, 950, 925, 900, 850, 800, 700, 600, 500, 400])
thorV= Tdata(year, month, day, hour, level='vertical', plevels= plevels)
#    thor= Tdata(year, month, day, hour, level='one_profile', dropid= dropid)


idrop= tuple([i for i in range(len(thorV.dropnr)) if thorV.dropnr[i] not in excl])

thorpres= thorV.pres[idrop, :]
thorT= thorV.T[idrop, :]
thorRH= thorV.RH[idrop,:]     
thortheta= thorV.theta[idrop,:]
thorU= thorV.U[idrop, :]

p= np.average(thorpres, axis= 0)* units.hPa
T= np.average(thorT, axis= 0)* units.degC       
Td= dewpoint_rh(T, np.average(thorRH, axis= 0) )

skew.plot(p, T, 'r', label= 'Thor T')
skew.plot(p, Td, 'g', label= 'Thor Td')

#skew.plot_barbs(p, u, v)
skew.ax.set_ylim(950, 500)
skew.ax.set_xlim(-30, 0)


skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
skew.plot_mixing_lines()





"""get the AA data and also plot it"""

#exp_name= '080303_warmctr'
#exp_name= '080303_cold_pseudo2'
exp_name= '080304_cold_pseudo'
#exp_name= '080303_warmsens_noTH'
fileday, filehour= 4, 0        
AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
AA= Adata(filename= AAfilename, res=1)


#AAT= np.zeros(np.shape(thorT))
#AARH= np.zeros(np.shape(thorT))
#AARH_fSH= np.zeros(np.shape(thorT))
#
#AAU= np.zeros(np.shape(thorT))
#
#for di, dropind in enumerate(idrop):
#    print(di)
#    thour= thorV.datetime[dropind, 2].hour #the hour where the dropsonde reached the ground    
##    hour= int(thorV.UTC10[dropind][:2])
#
#    t= (day- fileday)*24 + (thour- filehour)  # -1 
#    
#    dist= (AA.lat - thorV.lat[dropind, 2])**2 + (np.cos(np.deg2rad(thorV.lat[dropind, 2])) * (AA.lon - thorV.lon[dropind, 2])**2)
#    x,y= np.where(dist == np.min(dist))
#    
##    AA.imp_cross_sec_reduced(xn= x[0], yn= y[0], tn= t)
#    AA.imp_cross_sec_grib(xn= x[0], yn= y[0], tn= t)
#
#    
#    AAT[di]= interp1d(AA.pres, AA.T)(plevels)
#    AASH= interp1d(AA.pres, AA.SH)(plevels)
#    AARH[di]= relative_humidity_from_specific_humidity(AASH, AAT[di]* units.K, plevels*units.hPa)
#    
#    AAU[di]= interp1d(AA.pres, np.sqrt(AA.u**2 + AA.v**2))(plevels)




T= (np.average(AAT, axis= 0)* units.K).to('degC')
RH= np.average(AARH, axis= 0)

skew.plot(p, T, 'orange', label='Arome T')

Td= dewpoint_rh(T, RH)
skew.plot(p, Td, 'b', label= 'Arome Td')

plt.ylabel('Altitude [hPa]')
plt.xlabel(r"Temperature [$^{\circ}$C]")
plt.legend()
plt.tight_layout()


plt.savefig(savedir+ 'Sounding_Dropsonde_all_d'+str(day)+'_h'+str(hour))


"""profiles"""
plt.figure(fignr)
fignr +=1
plt.clf()

ax1= plt.subplot(131)
#plt.plot(thor.theta, thor.pres, 'r', label='Thor')
plt.plot(np.average(thortheta, axis=0), np.average(thorpres, axis= 0), 'r', label='Thor')

AAtheta= potential_temperature(plevels*units.hPa, AAT* units.K)

plt.plot(np.average(AAtheta, axis= 0), plevels, color='orange', label='Arome')

#theta= potential_temperature(AA2.pres*units.hPa, AA2.T*units.K)
#plt.plot(theta, AA2.pres, 'green', marker= 'o', ls= 'None', label='Arome_2')

plt.xlim([272,283])
#plt.xlim([269, 281])
plt.ylabel('Altitude [hPa]')
plt.xlabel('Potential temperature [K]')
plt.legend()

ax2= plt.subplot(132, sharey= ax1)
plt.plot(np.average(thorRH, axis= 0), np.average(thorpres, axis= 0), 'g', label='Thor')

plt.setp(ax2.get_yticklabels(), visible=False)
#theta= potential_temperature(AA.pres*units.hPa, AA.T*units.K)

plt.plot(np.average(AARH, axis=0), plevels, 'b', label='Arome')
plt.xlabel('Relative humidity')

plt.xlim([0,1])
plt.legend()

ax3= plt.subplot(133, sharey= ax1)
plt.plot(np.average(thorU, axis= 0), plevels, 'g', label='Thor')

plt.setp(ax3.get_yticklabels(), visible=False)
#theta= potential_temperature(AA.pres*units.hPa, AA.T*units.K)

plt.plot(np.average(AAU, axis= 0), plevels, 'b', label='Arome')
plt.xlabel('Horizontal velocity [m/s]')

plt.ylim([1000, 500])
#plt.xlim([0,1])
plt.legend()
plt.tight_layout()

plt.savefig(savedir+ 'Profile_Dropsonde_all_d'+str(day)+'_h'+str(hour))


"""difference profiles"""
plt.figure(fignr)
fignr +=1
plt.clf()

#plt.title('AROME - Thorpex')

ax1= plt.subplot(131)

plt.plot(np.average(AAtheta, axis= 0)- np.average(thortheta, axis=0), plevels, 'r', label='BIAS')
#plt.plot(np.sqrt(np.average((AAtheta - thortheta*units.K)**2, axis= 0)), plevels, 'b', label= 'RMSE')
plt.plot(np.average(np.abs(AAtheta - thortheta*units.K), axis= 0), plevels, 'g', label= 'MAE')
#plt.plot(-np.sqrt(np.average((AAtheta - thortheta*units.K)**2, axis= 0)), plevels, 'b', label= '-RMSE')
plt.plot(-np.average(np.abs(AAtheta - thortheta*units.K), axis= 0), plevels, 'g', label= '-MAE')

ax1.axvline(x=0, color='k')
plt.ylabel('Altitude [hPa]')
plt.xlabel('Potential temperature [K]')
#plt.legend()

print('RMSE T: '+ str(np.sqrt(np.nanmean((AAtheta - thortheta*units.K)**2))))

ax2= plt.subplot(132, sharey= ax1)

plt.setp(ax2.get_yticklabels(), visible=False)

plt.plot(np.average(AARH, axis=0)- np.average(thorRH, axis= 0), plevels, 'r', label='BIAS')

#plt.plot(np.sqrt(np.average((AARH - thorRH)**2, axis= 0)), plevels, 'b', label= 'RMSE')
plt.plot(np.average(np.abs(AARH - thorRH), axis= 0), plevels, 'g', label= 'MAE')
plt.plot(np.average(-np.abs(AARH - thorRH), axis= 0), plevels, 'g')


print('RMSE RH: '+ str(np.sqrt(np.nanmean((AARH - thorRH)**2))))


ax2.axvline(x=0, color='k')

plt.xlabel('Relative humidity')

#plt.xlim([0,1])

ax3= plt.subplot(133, sharey= ax1)

plt.setp(ax3.get_yticklabels(), visible=False)
#theta= potential_temperature(AA.pres*units.hPa, AA.T*units.K)

plt.plot(np.average(AAU, axis= 0)- np.average(thorU, axis= 0), plevels, 'r', label='BIAS')
#plt.plot(-np.sqrt(np.average((AAU - thorU)**2, axis= 0)), plevels, 'b', label= 'RMSE')
plt.plot(np.average(np.abs(AAU - thorU), axis= 0), plevels, 'g', label= 'MAE')
plt.plot(-np.average(np.abs(AAU - thorU), axis= 0), plevels, 'g')

ax3.axvline(x=0, color='k')
plt.legend(loc= 1)


print('RMSE U: '+ str(np.sqrt(np.nanmean((AAU - thorU)**2))))


plt.xlabel('Horizontal velocity [m/s]')

plt.ylim([1000, 500])


plt.tight_layout()

plt.savefig(savedir+ 'Profile_Diff_Dropsonde_all_d'+str(day)+'_h'+str(hour))


#plt.figure(fignr)
#fignr +=1
#plt.clf()
#
#
#iend = len(thorV.pres[idrop])
#
#ax1= plt.subplot(131)
##plt.plot(AAtheta[:idend] - thor.theta, thor.pres, 'r', label='hor')
#plt.plot(AAtheta[:iend] -thorV.theta[idrop]* units.K, thorV.pres[idrop], 'r')
#
##plt.plot(theta, AA.pres, 'orange', marker= 'o', ls= 'None', label='Arome')
#
##theta= potential_temperature(AA2.pres*units.hPa, AA2.T*units.K)
##plt.plot(theta, AA2.pres, 'green', marker= 'o', ls= 'None', label='Arome_2')
#
##plt.ylim([1000, 300])
##plt.xlim([260,300])
#plt.ylabel('Altitude [hPa]')
#plt.title('Potential temperature')
#plt.legend()
#
#ax2= plt.subplot(132, sharey= ax1)
##plt.plot(thor.RH, thor.pres, 'g', label='Thor')
#plt.plot(AA.RH[:iend] - thorV.RH[idrop], thorV.pres[idrop], 'go')
#
#plt.setp(ax2.get_yticklabels(), visible=False)
##theta= potential_temperature(AA.pres*units.hPa, AA.T*units.K)
#
##plt.plot(AA.RH, AA.pres, 'bo', label='Arome')
#plt.title('Relative humidity')
#
#plt.ylim([1000, 400])
##plt.xlim([0,1])
#plt.legend()
#
#ax3= plt.subplot(133, sharey= ax1)
##plt.plot(thor.U, thor.pres, 'g', label='Thor')
#plt.plot(np.sqrt(AA.u**2 + AA.v**2)[:iend] - thorV.U[idrop], thorV.pres[idrop], 'go')
#
#plt.setp(ax3.get_yticklabels(), visible=False)
##theta= potential_temperature(AA.pres*units.hPa, AA.T*units.K)
#
##plt.plot(np.sqrt(AA.u**2 + AA.v**2), AA.pres, 'bo', label='Arome')
#plt.title('Horizontal velocity')
#
#plt.ylim([1000, 300])
##plt.xlim([0,1])
#plt.legend()
#

#plt.tight_layout()
