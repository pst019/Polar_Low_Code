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
from f_useful import *
from scipy.interpolate import interp1d

from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import gridspec

plt.rcParams.update({'font.size': 14})

import numpy as np
import xarray as xr
from metpy.cbook import get_test_data
import metpy.calc as mpcalc
from metpy.units import units, concatenate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from metpy.calc.thermo import *
import scipy.ndimage.filters as filters

import seaborn as sns
from pandas import DataFrame


fignr= 1

Mediadir= '/media/'+user+'/PatsOrange/'

sound = 2
year, month = 2008, 3

save=True
savedir= '/home/pst019/home/Polar_Low/AromeArctic-vs-WRF/2008_03_03/Fuzzy/'#Sound'+str(sound)+'/'

#ulevel= 350
#ulevel= 450


plevel= [950, 850, 700, 500]

var= 'RH'
threshold= [.4, .6, .8, .9]

#var= 'SH'
#threshold= np.array([.1, .2, .5, 1, 2])*1E-3


var= 'U'
#threshold= [1, 2, 5, 10, 20]
threshold= [5, 10, 15, 20, 25]

#scale= [7.5, 15, 30, 60, 120, 240]
scale= [10, 20, 40, 80, 160, 320][::-1]

comptime='ana'
#comptime= '+12'
#comtime=24

skill_score='HK' #true skill score: hit rate - false alarm rate
#skill_score='POD' #hit rate (not very interesting)
#skill_score='TS' #critical success index
#skill_score='ETS' #TS but includes random 

method= 'MECT' #multi-event contigency table

smoothAA=False
#smoothAA=True

#imp=False
imp=True

if imp:

    N_obs= 0
    hit_AA, miss_AA, false_AA, correct_AA= np.zeros((len(scale), len(threshold))), np.zeros((len(scale), len(threshold))), np.zeros((len(scale), len(threshold))), np.zeros((len(scale), len(threshold)))
    hit_HRES, miss_HRES, false_HRES, correct_HRES= np.zeros((len(scale), len(threshold))), np.zeros((len(scale), len(threshold))), np.zeros((len(scale), len(threshold))), np.zeros((len(scale), len(threshold)))
    yes_AA, yes_Obs, yes_HRES= np.zeros((len(scale), len(threshold))), np.zeros((len(scale), len(threshold))), np.zeros((len(scale), len(threshold)))
    
    
    for sound in [1,2,3]:
    
        if sound == 1:
            day, hour= 3, 12  
            drops= [x for x in range(1,20) if x not in [3, 5, 12]]  #for the first flight, 12- goes only down to 890hPa
        if sound == 2:
            day, hour= 3, 18  
            drops= [x for x in range(1, 15) if x not in [2, 11, 12, 14]] #14 - no data below 750hPa, 11 - only on 514 and 336hPa
        #    drops= [x for x in range(1, 11) if x not in [1, 9]]
        #    excl=[1,9,11] #for the second flight
        if sound == 3:
            day, hour= 4, 12 
        #    drops= [x for x in range(1, 20) if x not in [1,5, 13]] #moved south AA domain
            drops= [x for x in range(1, 20) if x not in [1, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]] # for the third flight - for original AA domain
    
    
    
        
        if sound ==1: 
            if comptime=='+12':
                exp_name= 'DA_080303_CTR'
                fileday, filehour= 3, 0 
                HRESfileday, HRESfilehour=3, 0 
        
            if comptime=='ana':
                exp_name= 'DA_080301_cycling'
                fileday, filehour= 3, 12 
                HRESfileday, HRESfilehour=3, 12
#            if comptime=='+24':
#                exp_name= 'DA_08030212'
#                fileday, filehour= 2, 12 
#                HRESfileday, HRESfilehour=3, 12
            
        if sound ==2:
            if comptime=='+12':
                exp_name= 'DA_080303_CTR'
                fileday, filehour= 3, 0 
                HRESfileday, HRESfilehour=3, 0 
            
            if comptime=='ana':            
                exp_name= 'DA_080301_cycling'
                fileday, filehour= 3, 18 
                HRESfileday, HRESfilehour=3, 18
            
        if sound== 3:
            if comptime=='+12':            
                exp_name= 'DA_080303_CTR'
                fileday, filehour= 3, 0
                HRESfileday, HRESfilehour= 3, 0         
        
            if comptime=='ana':        
                exp_name= 'DA_080301_cycling'
                fileday, filehour= 4, 12 
                HRESfileday, HRESfilehour=4, 12
        
        lacktime= day*24 + hour  - (HRESfileday*24 + HRESfilehour)
        
        """HRES"""
        #HRES_name="oper_ML_200803"+str(HRESfileday).zfill(2)+'_'+str(HRESfilehour).zfill(2)+'_'+str(lacktime).zfill(2)+"_vertical.nc"
        HRES_name="oper_PL_200803"+str(HRESfileday).zfill(2)+'_'+str(HRESfilehour).zfill(2)+".nc"
        
        HRES= xr.open_dataset(Mediadir+'ECMWF/HRES/'+HRES_name, decode_times=False)
        
        if 'plev' in list(HRES.coords):
            HRES['plev']/= 100
            HRES= HRES.sel(plev= plevel)
            
        elif 'level' in list(HRES.coords):
            HRES= HRES.rename({'level': 'plev', 'latitude': 'lat', 'longitude': 'lon'})
        
        HRES= HRES.rename({'q': 'SH', 't': 'T'})
        
        tHRES= ((day - HRESfileday)*24 + (hour -HRESfilehour) ) //3
        HRES= HRES.isel(time= tHRES)
        
        
        
        #drops= [5]
        
        
        for dropid in drops:
            N_obs+=1 * len(plevel)
            """get the thorpex data"""
            thor= Tdata(year, month, day, hour, level='one_profile', dropid= dropid)
        
            if var== 'RH': thorvar= interp1d(thor.pres, thor.RH)(plevel) 
            if var== 'SH': thorvar= specific_humidity_from_mixing_ratio(mixing_ratio_from_relative_humidity(
                     interp1d(thor.pres, thor.RH)(plevel),  interp1d(thor.pres, thor.T)(plevel)*units.degC,  plevel*units.hPa))
            if var== 'T': thorvar= interp1d(thor.pres, thor.T)(plevel) +273.15
            if var== 'U': thorvar= interp1d(thor.pres, thor.U)(plevel) 
    
            
        
                
            """get the AA data """
            AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
           
            AA= xr.open_dataset(AAfilename)
            AAt= np.argmin(np.abs(np.datetime64(thor.datetime[0]) - AA.time))  
        
            AA= AA.isel(time= AAt)
            AA= AA.sel(pressure0= plevel)
            AA= AA.rename({'specific_humidity_pl': 'SH', 'x_wind_pl':'u', 'y_wind_pl': 'v', 'air_temperature_pl':'T'})
            print('AA t ', AA.time.values)
            
            if var== 'RH':
                AA['RH']= (('pressure0', 'y', 'x'), relative_humidity_from_specific_humidity(AA.SH.values , AA.T.values* units.K, AA.pressure0.values[:, np.newaxis, np.newaxis]*units.hPa) )
                HRES['RH']= (('plev', 'lat', 'lon'), relative_humidity_from_specific_humidity(HRES.SH.values , HRES.T.values* units.K, HRES.plev.values[:, np.newaxis, np.newaxis]*units.hPa) )
        
            if var=='U':
                AA['U']= (('pressure0', 'y', 'x'), np.sqrt(AA.u**2+ AA.v**2) )
                HRES['U']= (('plev', 'lat', 'lon'), np.sqrt(HRES.u**2+ HRES.v**2) )
        
            AA['dist']= (('y', 'x'), distance((AA.latitude, AA.longitude), (thor.lat[0], thor.lon[0])) )    
            HRES['dist']= (('lat', 'lon'), distance((HRES.lat, HRES.lon), (thor.lat[0], thor.lon[0])) )
        

            if smoothAA:
                #gaussfilter
                AA[var +'_smooth']= (('pressure0', 'y', 'x'),  filters.gaussian_filter1d(filters.gaussian_filter1d(AA[var], sigma= 12.5/2.5, mode='nearest', truncate= 1., axis= -1) , sigma= 12.5/2.5, mode='nearest', truncate= 1., axis= -2 ) )
                #uniform filter = local average
#                AA[var +'_smooth']= (('pressure0', 'y', 'x'),  filters.uniform_filter1d(filters.uniform_filter1d(AA[var], size= int(12.5/2.5), mode='nearest', axis= -1) , size= int(12.5/2.5), mode='nearest', axis= -2 ) )

                varAA= var +'_smooth'
                
            else: varAA= var
        
            for si, sca in enumerate(scale):
                AAy,AAx= np.where(AA.dist <= sca)
                AAx= AA.isel(x=xr.DataArray(AAx, dims='z'), y =xr.DataArray(AAy, dims='z'))
            
                
                HRx,HRy= np.where(HRES.dist <= sca)
                HRESx= HRES.isel(lat=xr.DataArray(HRx, dims='z'), lon =xr.DataArray(HRy, dims='z'))    
            
            
                for pi, plev in enumerate(plevel):
                    AAxx= AAx.sel(pressure0= plev)
                    HRESxx= HRESx.sel(plev= plev)
                    
                    for ti, thresh in enumerate(threshold):
                        I_AA, I_HRES, I_Obs= 0, 0, 0
                        
                        if method=='MECT':
                            if np.max(AAxx[varAA]) >= thresh: I_AA= 1
                            if np.max(HRESxx[var]) >= thresh: I_HRES= 1
                            if thorvar[pi] >= thresh: I_Obs= 1
                    
                        if I_AA==1 and I_Obs==1: hit_AA[si, ti]+=1
                        if I_AA==1 and I_Obs==0: false_AA[si, ti]+=1
                        if I_AA==0 and I_Obs==1: miss_AA[si, ti]+=1
                        if I_AA==0 and I_Obs==0: correct_AA[si, ti]+=1
                        
                        if I_HRES==1 and I_Obs==1: hit_HRES[si, ti]+=1
                        if I_HRES==1 and I_Obs==0: false_HRES[si, ti]+=1
                        if I_HRES==0 and I_Obs==1: miss_HRES[si, ti]+=1
                        if I_HRES==0 and I_Obs==0: correct_HRES[si, ti]+=1
        
                        if I_Obs==1: yes_Obs[si, ti] +=1
                        if I_AA==1: yes_AA[si, ti] +=1
                        if I_HRES==1: yes_HRES[si, ti] +=1                

print('threshold: ', threshold)
print('scale: ', scale)

print('hit AA \n', hit_AA)
print('miss AA \n', miss_AA)
print('false alarm AA \n', false_AA)
print('correct reject AA \n', correct_AA )

print('hit HRES \n', hit_HRES)
print('miss HRES \n', miss_HRES)
print('false alarm HRES \n', false_HRES)
print('correct reject HRES \n', correct_HRES )

POD_AA= hit_AA/(hit_AA+ miss_AA)
F_AA= false_AA/(correct_AA+ false_AA)
HK_AA= POD_AA- F_AA
TS_AA=hit_AA/(hit_AA+ miss_AA+ false_AA)

hit_rand_AA= 1/N_obs * (yes_Obs * yes_AA)
ETS_AA= (hit_AA- hit_rand_AA)/(hit_AA+ miss_AA+ false_AA- hit_rand_AA)

if skill_score=='POD': skill_AA= np.round(POD_AA, 2)
if skill_score=='HK': skill_AA= np.round(HK_AA, 2)
if skill_score=='TS': skill_AA= np.round(TS_AA, 2)
if skill_score=='ETS': skill_AA= np.round(ETS_AA, 2)

print('skill score AA: \n', skill_AA)

POD_HRES= hit_HRES/(hit_HRES+ miss_HRES)
F_HRES= false_HRES/(correct_HRES+ false_HRES)
HK_HRES= POD_HRES- F_HRES
TS_HRES=hit_HRES/(hit_HRES+ miss_HRES+ false_HRES)

hit_rand_HRES= 1/N_obs * (yes_Obs * yes_HRES)
ETS_HRES= (hit_HRES- hit_rand_HRES)/(hit_HRES+ miss_HRES+ false_HRES- hit_rand_HRES)

if skill_score=='POD': skill_HRES= np.round(POD_HRES, 2)
if skill_score=='HK': skill_HRES= np.round(HK_HRES, 2)
if skill_score=='TS': skill_HRES= np.round(TS_HRES, 2)
if skill_score=='ETS': skill_HRES= np.round(ETS_HRES, 2)

print('skill score HRES: \n', np.round(skill_HRES, 2))


"""plot AA"""
plt.close('all')

df= DataFrame(data= skill_AA, index= scale, columns= threshold)

sns.heatmap(df, annot=True, cmap='RdBu', vmin= -1, vmax= 1, cbar=False)#'viridis_r') YlGnBu'
plt.xlabel('Threshold')
plt.ylabel('Scale [km]')
plt.tight_layout()
if save: plt.savefig(savedir+ 'AA_'+comptime+'_'+method+'_'+skill_score+'_'+varAA)


"""plot HRES"""
fignr+= 1
plt.figure(fignr)
df= DataFrame(data= skill_HRES, index= scale, columns= threshold)

sns.heatmap(df, annot=True, cmap='RdBu', vmin= -1, vmax= 1, cbar=False)#'viridis_r') YlGnBu'
plt.xlabel('Threshold')
plt.ylabel('Scale [km]')
plt.tight_layout()
if save: plt.savefig(savedir+ 'HRES_'+comptime+'_'+method+'_'+skill_score+'_'+var)

"""plot AA- HRES"""
fignr+= 1
plt.figure(fignr)
df= DataFrame(data= skill_AA - skill_HRES, index= scale, columns= threshold)

sns.heatmap(df, annot=True, cmap='RdBu', vmin= -1, vmax= 1, cbar=False)#'viridis_r') YlGnBu'
plt.xlabel('Threshold')
plt.ylabel('Scale [km]')
plt.tight_layout()
if save: plt.savefig(savedir+ 'AA-HRES_'+comptime+'_'+method+'_'+skill_score+'_'+varAA)
