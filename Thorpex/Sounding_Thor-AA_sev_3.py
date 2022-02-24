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

Mediadir= '/media/'+user+'/PatsOrange/'
fignr= 9

year, month = 2008, 3

sound= 3

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
#exp_name= '080303_cold_pseudo2'
#exp_name= '08030312_cycling'
#exp_name= '080304_cold_pseudo'
#exp_name= '080303_warmsens_noTH'
#fileday, filehour= 4, 0

if sound in [1,2]:
    exp_name= '080303_cold_pseudo2'
    fileday, filehour= 3, 0 
if sound== 3:
    exp_name= '080304_cold_pseudo'
    fileday, filehour= 4, 0       

exp2= False
exp_name2= '08030312_cycling'
#exp_name2= '080303_cold_pseudo2'
fileday2, filehour2= 3, 12

impAA= True



plevels= np.arange(plevel0, 500, -1) #np.array([1000, 950, 925, 900, 850, 800, 700, 600, 500, 400])

savedir= '/home/pst019/home/Polar_Low/AromeArctic-vs-WRF/2008_03_03/Soundings4/Sound'+str(sound)+'/'


"""get Thorpex data"""
thorV= Tdata(year, month, day, hour, level='vertical', plevels= plevels)
#    thor= Tdata(year, month, day, hour, level='one_profile', dropid= dropid)


idrop= tuple([i for i in range(len(thorV.dropnr)) if thorV.dropnr[i] not in excl])
thorpres= thorV.pres[idrop, :]
thorT= thorV.T[idrop, :]
thorRH= thorV.RH[idrop,:]     
thortheta= thorV.theta[idrop,:]
thoru= thorV.u[idrop, :]
thorv= thorV.v[idrop, :]
thorU = np.sqrt(thoru**2 + thorv**2)

Thorpres= np.average(thorpres, axis= 0)* units.hPa
ThorTavg= np.average(thorT, axis= 0)* units.degC       
ThorTdavg= dewpoint_rh(ThorTavg, np.average(thorRH, axis= 0) )
Thorthetaavg= np.average(thortheta, axis=0)
Thoruavg, Thorvavg, ThorUavg = np.average(thoru, axis= 0), np.average(thorv, axis= 0), np.average(thorU, axis= 0)

"""get the AA data and also plot it"""        
AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
AA= Adata(filename= AAfilename, res=1)


if impAA== True:
    AAT= np.zeros(np.shape(thorT))
    AARH= np.zeros(np.shape(thorT))
    AARH_fSH= np.zeros(np.shape(thorT))
    
    AAu, AAv= np.zeros(np.shape(thorT)), np.zeros(np.shape(thorT))
    
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
#        AAu[di]= interp1d(AA.pres, AA.u)(plevels)
#        AAv[di]= interp1d(AA.pres, AA.v)(plevels)
    
        del_lon= (AA.lon[x[0]+1, y[0]]- AA.lon[x[0], y[0]])* np.cos(np.deg2rad(AA.lat[x[0], y[0]]))
        del_lat= AA.lat[x[0]+1, y[0]]- AA.lat[x[0], y[0]]
        
        rot= np.rad2deg(np.arctan2(del_lon, del_lat)) #how much to rotate the gridcell anticlockwise to be oriented N-S       
        dir_raw= np.rad2deg(np.arctan2(AA.u, AA.v)) #the wind direction in respect to the grid (0 being upwards, 90 rightwards, -90 leftw)   
        dir_rot= (dir_raw + rot)%360 #the wind direction with respect to the North pole (0 Northwards, 90 E, 270 W)
    
        wind_vel= np.sqrt(AA.u**2 +AA.v**2)
        AAurot= wind_vel* np.sin(np.deg2rad(dir_rot))
        AAvrot= wind_vel* np.cos(np.deg2rad(dir_rot))
           
        AAu[di] = interp1d(AA.pres, AAurot)(plevels)
        AAv[di] = interp1d(AA.pres, AAvrot)(plevels)        



    AATavg= (np.average(AAT, axis= 0)* units.K).to('degC')
    AARHavg= np.average(AARH, axis= 0)
    AATdavg= dewpoint_rh(AATavg, AARHavg)
    AAtheta= potential_temperature(plevels*units.hPa, AAT* units.K)
    AAthetaavg= np.average(AAtheta, axis= 0)
    AAuavg, AAvavg = np.average(AAu, axis= 0), np.average(AAv, axis= 0)

    if exp2== True:
        AAfilename= Mediadir+'PL/AA/ec/'+exp_name2+'_'+str(year)+str(month).zfill(2)+str(fileday2).zfill(2)+str(filehour2).zfill(2)+'_fp_extract.nc'
        AA= Adata(filename= AAfilename, res=1)
        
        AAT2= np.zeros(np.shape(thorT))
        AARH2= np.zeros(np.shape(thorT))
        AARH_fSH2= np.zeros(np.shape(thorT))
        
        AAu2, AAv2= np.zeros(np.shape(thorT)), np.zeros(np.shape(thorT))
        
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
            AAu2[di]= interp1d(AA.pres, AA.u)(plevels)
            AAv2[di]= interp1d(AA.pres, AA.v)(plevels)
    
    
        AATavg2= (np.average(AAT2, axis= 0)* units.K).to('degC')
        AARHavg2= np.average(AARH2, axis= 0)
        AATdavg2= dewpoint_rh(AATavg2, AARHavg2)
        AAtheta2= potential_temperature(plevels*units.hPa, AAT2* units.K)
        AAthetaavg2= np.average(AAtheta2, axis= 0)
        AAuavg2, AAvavg2 = np.average(AAu2, axis= 0), np.average(AAv2, axis= 0)

AAU= np.sqrt(AAu**2 + AAv**2)
AAUavg = np.average(AAU, axis= 0)

"""plot average meteogram"""
fig = plt.figure(fignr, figsize=(6, 5))
fignr += 1
plt.clf()
plt.rcParams.update({'font.size': 14})

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
fig= plt.figure(fignr, figsize=(6, 5))
fignr +=1
plt.clf()
plt.rcParams.update({'font.size': 14})
#fig.subplots_adjust(wspace=0, hspace=0)


ax1= plt.subplot(131)

#plot moist adiabats
moistadiabat= mpcalc.potential_temperature(plevels*units.hPa, mpcalc.moist_lapse(plevels*units.hPa, np.arange(260,286,2)*units.K))
for i in range(moistadiabat.shape[0]): plt.plot(moistadiabat[i], plevels, color='grey', linewidth=1, linestyle='dashed')

#plt.plot(thor.theta, thor.pres, 'r', label='Thor')
plt.plot(Thorthetaavg, plevels, 'r', label='Drop')

if exp2== True: plt.plot(AAthetaavg2, plevels, color='g', label='Arome +'+str(hour-filehour2))
plt.plot(AAthetaavg, plevels, color='b', label='AA +'+str(hour-filehour))


#theta= potential_temperature(AA2.pres*units.hPa, AA2.T*units.K)
#plt.plot(theta, AA2.pres, 'green', marker= 'o', ls= 'None', label='Arome_2')

plt.xlim([np.min([Thorthetaavg, AAthetaavg])-1, np.max([Thorthetaavg, AAthetaavg])+.3 ])
plt.ylabel('Altitude [hPa]')
plt.xlabel('Pot. temp. [K]')

ax2= plt.subplot(132, sharey= ax1)
plt.plot(np.average(thorRH, axis= 0), plevels, 'r', label='Drop')

plt.setp(ax2.get_yticklabels(), visible=False)
#theta= potential_temperature(AA.pres*units.hPa, AA.T*units.K)

plt.plot(AARHavg, plevels, 'b', label='AA')
if exp2== True: plt.plot(AARHavg2, plevels, color='g', label='AA +'+str(hour-filehour2))

plt.xlabel('Rel. hum.')

plt.xlim([-0.05,1.05])
plt.legend(loc=2, fontsize= 12)


ax3= plt.subplot(133, sharey= ax1)

plt.plot(ThorUavg, plevels, 'r', label='Thor', lw=0.5)
plt.quiver(ThorUavg[::40], plevels[::40], Thoruavg[::40]/ThorUavg[::40], Thorvavg[::40]/ThorUavg[::40], scale= 7, width= 0.02, color='r', pivot='mid')



plt.plot(AAUavg, plevels, 'b', label='Arome', lw=0.5)
plt.quiver(AAUavg[::40], plevels[::40], AAuavg[::40]/AAUavg[::40], AAvavg[::40]/AAUavg[::40], scale= 7, width= 0.02, color='b', pivot='mid')

if exp2== True:
    AAUavg2 = np.sqrt(AAuavg2**2 + AAvavg2**2)
    plt.plot(AAUavg2, plevels, color='g', label='Arome +'+str(hour-filehour2))
    plt.quiver(AAUavg2[::40], plevels[::40], AAuavg2[::40]/AAUavg2[::40], AAvavg2[::40]/AAUavg2[::40], scale= 7, width= 0.02, color='g', pivot='mid')

plt.xlabel('Hor. vel. [m/s]')

plt.ylim([plevels[0], plevels[-1]])
plt.setp(ax3.get_yticklabels(), visible=False)

plt.xlim([5.5,21])
plt.tight_layout()

if exp2== True: plt.savefig(savedir+ 'AA'+exp_name+'+'+exp_name2+'_Profile_all_d'+str(day)+'_h'+str(hour))
else: plt.savefig(savedir+ 'AA'+exp_name+ '_Profile_all_d'+str(day)+'_h'+str(hour))



"""difference profiles"""
fig= plt.figure(fignr)
fignr +=1
plt.clf()
plt.rcParams.update({'font.size': 12})

#plt.title('AROME - Thorpex')
ax1= plt.subplot(131)

plt.plot(AAthetaavg- Thorthetaavg, plevels, 'r', label='BIAS')
#plt.plot(np.sqrt(np.average((AAtheta - thortheta*units.K)**2, axis= 0)), plevels, 'b', label= 'RMSE')
plt.plot(np.average(np.abs(AAtheta - thortheta*units.K), axis= 0), plevels, 'g', label= 'MAE')
#plt.plot(-np.sqrt(np.average((AAtheta - thortheta*units.K)**2, axis= 0)), plevels, 'b', label= '-RMSE')
plt.plot(-np.average(np.abs(AAtheta - thortheta*units.K), axis= 0), plevels, 'g', label= '-MAE')

ax1.axvline(x=0, color='k')
plt.ylabel('Altitude [hPa]')
plt.xlabel('Pot. temp. [K]')
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

plt.xlabel('Rel. hum.')

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


plt.xlabel('Hor. vel. [m/s]')

#plt.ylim([1000, 500])
plt.ylim([plevels[0], plevels[-1]])

plt.tight_layout()

plt.savefig(savedir+ 'AA'+exp_name+ '_Difference_all_d'+str(day)+'_h'+str(hour))




"""plot all profiles"""
#for dropi in idrop:
    
#plt.figure(fignr)
#fignr +=1
#plt.clf()
