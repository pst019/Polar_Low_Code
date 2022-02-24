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
fignr= 9

year, month = 2008, 3

sound= 1

if sound == 1:
    day, hour= 3, 12  
    excl= [3, 5, 12] #for the first flight
    plevel0= 985
if sound == 2:
    day, hour= 3, 17  
    excl=[1,9,11] #for the second flight
    plevel0= 985 #lowest pressure level
if sound == 3:
    day, hour= 4, 12  
    excl= [1, 5, 13]# for the third flight
    plevel0= 975

#exclude= 3, 5, 6,7( above 500hPa)

#exp_name= '080303_warmctr'
exp_name= '080303_cold_pseudo2'
#exp_name= '08030312_cycling'
#exp_name= '080304_cold_pseudo'
#exp_name= '080303_warmsens_noTH'
fileday, filehour= 3, 0

exp2= False
exp_name2= '08030312_cycling'
#exp_name2= '080303_cold_pseudo2'

fileday2, filehour2= 3, 12

impAA= True



plevels= np.arange(plevel0, 500, -1) #np.array([1000, 950, 925, 900, 850, 800, 700, 600, 500, 400])

savedir= '/home/pst019/home/Polar_Low/AromeArctic-vs-WRF/2008_03_03/Soundings3/Sound'+str(sound)+'/'


"""get Thorpex data"""
thorV= Tdata(year, month, day, hour, level='vertical', plevels= plevels)
#    thor= Tdata(year, month, day, hour, level='one_profile', dropid= dropid)


idrop= tuple([i for i in range(len(thorV.dropnr)) if thorV.dropnr[i] not in excl])
thorpres= thorV.pres[idrop, :]
thorT= thorV.T[idrop, :]
thorRH= thorV.RH[idrop,:]     
thortheta= thorV.theta[idrop,:]
thorU= thorV.U[idrop, :]

Thorpres= np.average(thorpres, axis= 0)* units.hPa
ThorTavg= np.average(thorT, axis= 0)* units.degC       
ThorTdavg= dewpoint_rh(ThorTavg, np.average(thorRH, axis= 0) )
Thorthetaavg= np.average(thortheta, axis=0)

"""get the AA data and also plot it"""        
AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
AA= Adata(filename= AAfilename, res=1)


if impAA== True:
    AAT= np.zeros(np.shape(thorT))
    AARH= np.zeros(np.shape(thorT))
    AARH_fSH= np.zeros(np.shape(thorT))
    
    AAU= np.zeros(np.shape(thorT))
    
    for di, dropind in enumerate(idrop):
        print('Dropindex ', dropind)
        thour= thorV.datetime[dropind, 2].hour #the hour where the dropsonde reached the ground    
    
        t= (day- fileday)*24 + (thour- filehour)  # -1 
        
        dist= (AA.lat - thorV.lat[dropind, 2])**2 + (np.cos(np.deg2rad(thorV.lat[dropind, 2])) * (AA.lon - thorV.lon[dropind, 2])**2)
        x,y= np.where(dist == np.min(dist))
        
        if exp_name in ['08030312_cycling']: t= int(np.round(t/3))
        AA.imp_cross_sec_grib(xn= x[0], yn= y[0], tn= t)
    
        
        AAT[di]= interp1d(AA.pres, AA.T)(plevels)
        AASH= interp1d(AA.pres, AA.SH)(plevels)
        AARH[di]= relative_humidity_from_specific_humidity(AASH, AAT[di]* units.K, plevels*units.hPa)
        AAU[di]= interp1d(AA.pres, np.sqrt(AA.u**2 + AA.v**2))(plevels)


    AATavg= (np.average(AAT, axis= 0)* units.K).to('degC')
    AARHavg= np.average(AARH, axis= 0)
    AATdavg= dewpoint_rh(AATavg, AARHavg)
    AAtheta= potential_temperature(plevels*units.hPa, AAT* units.K)
    AAthetaavg= np.average(AAtheta, axis= 0)


    if exp2== True:
        AAfilename= Mediadir+'PL/AA/ec/'+exp_name2+'_'+str(year)+str(month).zfill(2)+str(fileday2).zfill(2)+str(filehour2).zfill(2)+'_fp_extract.nc'
        AA= Adata(filename= AAfilename, res=1)
        
        AAT2= np.zeros(np.shape(thorT))
        AARH2= np.zeros(np.shape(thorT))
        AARH_fSH2= np.zeros(np.shape(thorT))
        
        AAU2= np.zeros(np.shape(thorT))
        
        for di, dropind in enumerate(idrop):
            print('Dropindex ', dropind)
            thour= thorV.datetime[dropind, 2].hour #the hour where the dropsonde reached the ground    
        
            t= (day- fileday2)*24 + (thour- filehour2)  # -1 
            
            dist= (AA.lat - thorV.lat[dropind, 2])**2 + (np.cos(np.deg2rad(thorV.lat[dropind, 2])) * (AA.lon - thorV.lon[dropind, 2])**2)
            x,y= np.where(dist == np.min(dist))
            
            if exp_name2 in ['08030312_cycling']: t= int(np.round(t/3))
            AA.imp_cross_sec_grib(xn= x[0], yn= y[0], tn= t)
        
            
            AAT2[di]= interp1d(AA.pres, AA.T)(plevels)
            AASH2= interp1d(AA.pres, AA.SH)(plevels)
            AARH2[di]= relative_humidity_from_specific_humidity(AASH2, AAT2[di]* units.K, plevels*units.hPa)
            AAU2[di]= interp1d(AA.pres, np.sqrt(AA.u**2 + AA.v**2))(plevels)
    
    
        AATavg2= (np.average(AAT2, axis= 0)* units.K).to('degC')
        AARHavg2= np.average(AARH2, axis= 0)
        AATdavg2= dewpoint_rh(AATavg2, AARHavg2)
        AAtheta2= potential_temperature(plevels*units.hPa, AAT2* units.K)
        AAthetaavg2= np.average(AAtheta2, axis= 0)




"""plot average meteogram"""
fig = plt.figure(fignr, figsize=(6, 5))
fignr += 1
plt.clf()
skew = SkewT(fig, rotation=45)
skew.plot(Thorpres, ThorTavg, 'r', label= 'Thor T')
skew.plot(Thorpres, ThorTdavg, 'g', label= 'Thor Td')

#skew.plot_barbs(p, u, v)
skew.ax.set_ylim(950, 500)
skew.ax.set_xlim(-30, 0)
skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
skew.plot_mixing_lines()

skew.plot(Thorpres, AATavg, 'orange', label='Arome T')
skew.plot(Thorpres, AATdavg, 'b', label= 'Arome Td')

plt.ylabel('Altitude [hPa]')
plt.xlabel(r"Temperature [$^{\circ}$C]")
plt.legend()
plt.tight_layout()


plt.savefig(savedir+ 'AA'+exp_name+'_Sounding_all_d'+str(day)+'_h'+str(hour))


"""plot average profiles"""
plt.figure(fignr)
fignr +=1
plt.clf()

ax1= plt.subplot(131)

#plot moist adiabats
moistadiabat= mpcalc.potential_temperature(plevels*units.hPa, mpcalc.moist_lapse(plevels*units.hPa, np.arange(260,286,2)*units.K))
for i in range(moistadiabat.shape[0]): plt.plot(moistadiabat[i], plevels, color='grey', linewidth=1, linestyle='dashed')

#plt.plot(thor.theta, thor.pres, 'r', label='Thor')
plt.plot(Thorthetaavg, plevels, 'r', label='Thorpex')

if exp2== True: plt.plot(AAthetaavg2, plevels, color='g', label='Arome +'+str(hour-filehour2))
plt.plot(AAthetaavg, plevels, color='b', label='Arome +'+str(hour-filehour))


#theta= potential_temperature(AA2.pres*units.hPa, AA2.T*units.K)
#plt.plot(theta, AA2.pres, 'green', marker= 'o', ls= 'None', label='Arome_2')

plt.xlim([np.min([Thorthetaavg, AAthetaavg])-1, np.max([Thorthetaavg, AAthetaavg])+.3 ])
plt.ylabel('Altitude [hPa]')
plt.xlabel('Potential temperature [K]')
plt.legend(loc=2)

ax2= plt.subplot(132, sharey= ax1)
plt.plot(np.average(thorRH, axis= 0), plevels, 'r', label='Thorpex')

plt.setp(ax2.get_yticklabels(), visible=False)
#theta= potential_temperature(AA.pres*units.hPa, AA.T*units.K)

plt.plot(AARHavg, plevels, 'b', label='Arome')
if exp2== True: plt.plot(AARHavg2, plevels, color='g', label='Arome +'+str(hour-filehour2))

plt.xlabel('Relative humidity')

plt.xlim([0,1])
#plt.legend()

ax3= plt.subplot(133, sharey= ax1)
plt.plot(np.average(thorU, axis= 0), plevels, 'r', label='Thor')

plt.setp(ax3.get_yticklabels(), visible=False)
#theta= potential_temperature(AA.pres*units.hPa, AA.T*units.K)

plt.plot(np.average(AAU, axis= 0), plevels, 'b', label='Arome')
if exp2== True: plt.plot(np.average(AAU2, axis= 0), plevels, color='g', label='Arome +'+str(hour-filehour2))

plt.xlabel('Horizontal velocity [m/s]')

plt.ylim([plevels[0], plevels[-1]])
#plt.legend()
plt.tight_layout()

if exp2== True: plt.savefig(savedir+ 'AA'+exp_name+'+'+exp_name2+'_Profile_all_d'+str(day)+'_h'+str(hour))
else: plt.savefig(savedir+ 'AA'+exp_name+ '_Profile_all_d'+str(day)+'_h'+str(hour))



"""difference profiles"""
plt.figure(fignr)
fignr +=1
plt.clf()

#plt.title('AROME - Thorpex')
ax1= plt.subplot(131)

plt.plot(AAthetaavg- Thorthetaavg, plevels, 'r', label='BIAS')
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

#plt.ylim([1000, 500])
plt.ylim([plevels[0], plevels[-1]])

plt.tight_layout()

plt.savefig(savedir+ 'AA'+exp_name+ '_Difference_all_d'+str(day)+'_h'+str(hour))




"""plot all profiles"""
#for dropi in idrop:
    
#plt.figure(fignr)
#fignr +=1
#plt.clf()
