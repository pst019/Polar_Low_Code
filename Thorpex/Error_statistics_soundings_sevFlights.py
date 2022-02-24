#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  5 14:52:50 2018

@author: pst019

conda install -c conda-forge metpy

https://github.com/Unidata/MetPy/tree/master/staticdata

from Sounding_Thor-AA_sev
"""

import os
user = os.getcwd().split('/')[2]

import sys
sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')

from f_imp_thorpex import data as Tdata #import the thorpex data
from f_imp_AROME import data as ARdata #import the thorpex data
from f_imp_ERA2 import data as EIdata
#from f_imp_ASR import dataASR as ASRdata
import xarray as xr #for ASR
from netCDF4 import Dataset, num2date
from f_meteo import *

from datetime import datetime, timedelta
import numpy as np


Mediadir= '/media/'+user+'/PatsOrange/'




year, month = 2008, 3

plevel= 950 #no ERAI
#plevel= 850 #no ERAI wind
#plevel=700 #everything + geop
#plevel= 500 #everything
#plevel= 925 #only wind + geop for AA and ERAI (avail for ERA5)





thor_T= np.array([])
thor_U= np.array([])
thor_SH= np.array([])
thor_Geo= np.array([])

AA_T, AA_U, AA_SH, AA_Geo = np.array([]), np.array([]), np.array([]), np.array([])
AAavg_T, AAavg_U, AAavg_SH, AAavg_Geo= np.array([]), np.array([]), np.array([]), np.array([])

EI_T, EI_U, EI_SH, EI_Geo = np.array([]), np.array([]), np.array([]), np.array([])
E5_T, E5_U, E5_SH = np.array([]), np.array([]), np.array([])
HRES_T, HRES_U, HRES_SH = np.array([]), np.array([]), np.array([])
ASR_T, ASR_U, ASR_RH = np.array([]), np.array([]), np.array([])

excludecount= 0

for flight in [1,2,3]:

    """get thorpex data for this level"""
    exclude_dropsonde=[]
    if flight == 1:
        day, hour= 3, 12
        if plevel in [925, 950]: exclude_dropsonde= [12]
    #    excl= [3, 5] #for the first flight
        fileday, filehour= 3, 12
        
    if flight == 2:
        day, hour= 3, 18  
    #    excl=[1,9,11] #for the second flight
        exclude_dropsonde= [11]
        fileday, filehour= 3, 12
    
    if flight == 3:
        day, hour= 4, 12
    #    if plevel==950: exclude_dropsonde= [1]
        exclude_dropsonde= [1, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14] #the ones outside of the AROME domain
    #    excl= [5, 13]# for the third flight - exclude dropsonde 1 for 950hPa
        fileday, filehour= 4, 12
    
    
    thor= Tdata(year, month, day, hour, level='vertical', plevels= [plevel], exclude=exclude_dropsonde)
    
    thor_T= np.append(thor_T, thor.T[:,0] +273.15)
    thor_U= np.append(thor_U, thor.U[:,0])
    thor_SH= np.append(thor_SH, thor.SH[:,0])
    thor_Geo= np.append(thor_Geo, thor.alt[:,0])
    
    """get AROME data"""
    #exp_name= '080303_warmctr'
    #exp_name= '080303_cold_pseudo2'
    #exp_name= '08030312_cycling'
    #exp_name= '080304_cold_pseudo'
    #exp_name= '080303_warmsens_noTH'
    
    exp_name='DA_080301_cycling'
    
    
    filedatetime= datetime(year, month, fileday, filehour, 0)
    #tn= (thor.datetime[:,0] -filedatetime )/ timedelta(hours=1)
    #tn= [int(np.round(tni)) for tni in tn]
    #tn= remove_dublicate(tn)
    
    AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
    AA= ARdata(filename= AAfilename, res=1)
    
    if plevel== 850: pn= 4
    elif plevel== 700: pn= 6
    elif plevel== 500: pn= 8
    elif plevel== 950: pn=1
    elif plevel== 925: pn=2
    else: print('specifiy pn')
    
    
    
    for di in range(len(thor.dropnr)):
        if type(thor.datetime[di,0])== float: continue #this escludes nan values
    
        tn= int(np.round((thor.datetime[di,0] - filedatetime )/ timedelta(hours=1)))
        if exp_name in ['08030312_cycling']: tn= int(np.round(tn/3))
    
        if tn<0 : tn= 0 #otherwise if t is negative, the last-t timepoint is taken
    
        print(tn)
        AA.imp_level(pn= pn, tn= tn)
    
           
        dist= 110* np.sqrt((AA.lat - thor.lat[di, 0])**2 + (np.cos(np.deg2rad(thor.lat[di, 0])) * (AA.lon - thor.lon[di, 0]) )**2 )
        x,y= np.where(dist == np.min(dist))
            
        AA_T= np.append(AA_T, AA.T[x, y])
        AA_U= np.append(AA_U, np.sqrt(AA.u[x,y]**2 + AA.v[x,y]**2) )
        AA_SH= np.append(AA_SH, AA.SH[x,y] )
        AA_Geo= np.append(AA_Geo, AA.geop[x,y])
        
        x,y= np.where(dist <= 12.5)
        print('number cells', len(x))
        
        AAavg_T= np.append(AAavg_T, np.average(AA.T[x, y]) )
        AAavg_U= np.append(AAavg_U, np.average(np.sqrt(AA.u[x,y]**2 + AA.v[x,y]**2)) )
        AAavg_SH= np.append(AAavg_SH, np.average(AA.SH[x,y]) )
        AAavg_Geo= np.append(AAavg_Geo, np.average(AA.geop[x,y]) )
           
    
    
    """get ERA-I"""
    t= (day-1)*4 + hour//6
         
    if plevel== 850:
        EI= EIdata(['T850','sHum850'], year, month, tstart= t, tend= t+1) #have to get u
        EIT, EISH, EIU, EIGeo= EI.T850[0], EI.sHum850[0], np.zeros(np.shape(EI.T850[0])) , np.zeros(np.shape(EI.T850[0]))        
                       
    elif plevel== 700:
        EI= EIdata(['T','sHum','uPL', 'vPL', 'Geop'], year, month, tstart= t, tend= t+1) #have to get u
        EIT, EISH, EIU, EIGeo= EI.T[0,1], EI.sHum[0,1], np.sqrt(EI.uPL[0,1]**2+ EI.vPL[0,1]**2), EI.Geop[0,1]/9.81              
    
    elif plevel== 500:
        EI= EIdata(['T','sHum','uPL500', 'vPL500', 'Geop500'], year, month, tstart= t, tend= t+1) #have to get u
        EIT, EISH, EIU, EIGeo= EI.T[0,2], EI.sHum[0,2], np.sqrt(EI.uPL500[0]**2+ EI.vPL500[0]**2), EI.Geop500[0]/9.81
    
    elif plevel== 925:
        EI= EIdata(['uPL', 'vPL', 'Geop'], year, month, tstart= t, tend= t+1) #have to get u
        EIT, EISH, EIU, EIGeo= np.zeros(np.shape(EI.uPL[0,0])), np.zeros(np.shape(EI.uPL[0,0])), np.sqrt(EI.uPL[0,0]**2+ EI.vPL[0,0]**2), EI.Geop[0,0]/9.81              
    
    print(EI.datetime)
    EIlon, EIlat = np.meshgrid(EI.lon, EI.lat)
    for di in range(len(thor.dropnr)):
        if type(thor.datetime[di,0])== float: continue
        
        dist= (EIlat - thor.lat[di, 0])**2 + (np.cos(np.deg2rad(thor.lat[di, 0])) * (EIlon - thor.lon[di, 0])**2)
        x,y= np.where(dist == np.min(dist))
    
        EI_T= np.append(EI_T, EIT[x, y])
        EI_U= np.append(EI_U, EIU[x,y])
        EI_SH= np.append(EI_SH, EISH[x,y])
        EI_Geo= np.append(EI_Geo, EIGeo[x,y])
    
    """get ERA5"""
    varimp= ['T', 'SH', 'U']
    
    if plevel== 850: pn= 16
    elif plevel== 700: pn= 11
    elif plevel== 500: pn= 7
    elif plevel== 950: pn= 20
    
    for var in varimp:
        file= Mediadir + 'ECMWF/ERA5/'+str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+var+'.nc'
        nc= Dataset(file)
    
    
        if var== 'T':
            E5T= nc.variables['t'][:,pn]
    #        tim= nc.variables['time'][t]
    #        E5datetime= num2date(tim, nc.variables['time'].units)
            lat= nc.variables['latitude'][:]
            lon= nc.variables['longitude'][:]        
            E5level= nc.variables['level'][:]
            E5lon, E5lat = np.meshgrid(lon, lat)
            
        elif var== 'SH': E5SH= nc.variables['q'][:,pn]
        elif var== 'U' : E5U= np.sqrt(nc.variables['u'][:,pn]**2+ nc.variables['v'][:, pn]**2)
    
    
    for di in range(len(thor.dropnr)):
        if type(thor.datetime[di,0])== float: continue

        tn= int(np.round((thor.datetime[di,0] - datetime(year, month, day, 0, 0))/ timedelta(hours=1)))
    
        print('E5 tn', tn)
        dist= (E5lat - thor.lat[di, 0])**2 + (np.cos(np.deg2rad(thor.lat[di, 0])) * (E5lon - thor.lon[di, 0])**2)
        x,y= np.where(dist == np.min(dist))
            
        E5_T= np.append(E5_T, E5T[tn, x, y])
        E5_U= np.append(E5_U, E5U[tn,x,y])
        E5_SH= np.append(E5_SH, E5SH[tn,x,y])
        
    
    """get HRES"""
   
    
    varimp= ['T', 'SH', 'U']
    
    if plevel== 850: pn= 1
    elif plevel== 700: pn= 2
    elif plevel== 500: pn= 3
    elif plevel== 950: pn= 0
    
    file= Mediadir + 'ECMWF/HRES/oper_PL_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+'_'+str(filehour).zfill(2)+'.nc'
    nc= Dataset(file)
    
    tim= nc.variables['time'][:]
    HRESdatetime= num2date(tim, nc.variables['time'].units)
    
    if 'lat' in nc.variables.keys(): #there are two different conventions in the HRES files
        lat= nc.variables['lat'][:]
        lon= nc.variables['lon'][:]        
    #    HRESlevel= nc.variables['plev'][:]   
    else:
        lat= nc.variables['latitude'][:]
        lon= nc.variables['longitude'][:]        
    #    HRESlevel= nc.variables['level'][:]
    
    HRESlon, HRESlat = np.meshgrid(lon, lat)
    
    HREST= nc.variables['t'][:,pn]
    HRESSH= nc.variables['q'][:,pn]
    HRESU= np.sqrt(nc.variables['u'][:,pn]**2+ nc.variables['v'][:, pn]**2)
    
    
    for di in range(len(thor.dropnr)):
        if type(thor.datetime[di,0])== float: continue
        
        tn= np.argmin(np.abs(thor.datetime[di,0] - HRESdatetime)) 
    
        print('HRES tn', tn, 'HRESdatetime', HRESdatetime[tn])
        dist= (HRESlat - thor.lat[di, 0])**2 + (np.cos(np.deg2rad(thor.lat[di, 0])) * (HRESlon - thor.lon[di, 0])**2)
        x,y= np.where(dist == np.min(dist))
            
        HRES_T= np.append(HRES_T, HREST[tn, x, y])
        HRES_U= np.append(HRES_U, HRESU[tn,x,y])
        HRES_SH= np.append(HRES_SH, HRESSH[tn,x,y])



    """get ASR"""
    file= Mediadir+'/ASR/asr15km.anl.3D.'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'.nc'
    ASR= xr.open_dataset(file)
#    ASR= ASR.sel(Time= str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'T'+str(hour).zfill(2))
    
    index_pres= np.where(ASR['PRES'] == plevel * 100)[0][0]
    ASR= ASR.sel(num_metgrid_levels= index_pres)
    
    for di in range(len(thor.dropnr)):
        if type(thor.datetime[di,0])== float:
            excludecount += 1
            continue
        
        ASRhour= int(np.round(thor.datetime[di,0].hour/3)*3)
        print('ASR hour', ASRhour)
        ASRnow= ASR.sel(Time= str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'T'+str(ASRhour).zfill(2))
#        ASR= ASR.sel(Time= ASR.Time.values[tn], lat= xr.DataArray(thorlat, dims='location'), lon= xr.DataArray(thorlon, dims='location'), method='nearest')

        dist= (ASRnow.XLAT.values - thor.lat[di, 0])**2 + (np.cos(np.deg2rad(thor.lat[di, 0])) * (ASRnow.XLONG.values - thor.lon[di, 0])**2)
        x,y= np.where(dist == np.min(dist))
        ASRnow= ASRnow.isel(south_north= x, west_east= y)
        
        ASR_T= np.append(ASR_T, ASRnow.TT[0,0].values)
        ASR_U= np.append(ASR_U, np.sqrt(ASRnow.VE**2 + ASRnow.UE**2)[0,0].values)
        ASR_RH= np.append(ASR_RH, ASRnow.RH[0,0].values)        
    

"""nan exclusion"""
lenbefore= len(thor_T)
thor_SH= thor_SH[~np.isnan(thor_T)]
thor_U= thor_U[~np.isnan(thor_T)]
thor_T= thor_T[~np.isnan(thor_T)]


print('excluded nan-dropsondes in models: ', excludecount)
print('excluded in Thorpex: ', lenbefore -len(thor_T))

    
"""print statistics"""  
roundorder= 2

print('model & T_{BIAS} & T_{MAE} & RH_{BIAS} & RH_{MAE} & U_{BIAS} & U_{MAE}', r'\\')

print('AROME &', np.round( np.average(AA_T)- np.average(thor_T) , roundorder) ,
'&', np.round( np.average(np.abs(AA_T- thor_T)) , roundorder) ,
#'&', np.round( 1E3* (np.average(AA_SH)- np.average(thor_SH)) , roundorder) ,
#'&', np.round(  1E3*np.average(np.abs(AA_SH- thor_SH)) , roundorder) , r'\\' )
'&', np.round(  (np.average(SH2RH(AA_SH*1E3, plevel, AA_T))- np.average(SH2RH(thor_SH*1E3, plevel, thor_T))) , 1) ,
'&', np.round(  np.average(np.abs( SH2RH(AA_SH*1E3, plevel, AA_T) - SH2RH(thor_SH*1E3, plevel, thor_T))) , 1) , 
'&', np.round( np.average(AA_U)- np.average(thor_U) , roundorder) ,
'&', np.round( np.average(np.abs(AA_U- thor_U)) , roundorder) , r'\\' )

print('AROME avg &', np.round( np.average(AAavg_T)- np.average(thor_T) , roundorder) ,
'&', np.round( np.average(np.abs(AAavg_T- thor_T)) , roundorder) ,
#'&', np.round( 1E3* (np.average(AAavg_SH)- np.average(thor_SH)) , roundorder) ,
#'&', np.round(  1E3*np.average(np.abs(AAavg_SH- thor_SH)) , roundorder) , r'\\' )
'&', np.round( (np.average( SH2RH(AAavg_SH*1E3, plevel, AAavg_T) )- np.average( SH2RH(thor_SH*1E3, plevel, thor_T))) , 1) ,
'&', np.round( np.average(np.abs( SH2RH(AAavg_SH*1E3, plevel, AAavg_T) - SH2RH(thor_SH*1E3, plevel, thor_T))) , 1 ) ,
'&', np.round( np.average(AAavg_U)- np.average(thor_U) , roundorder) ,
'&', np.round( np.average(np.abs(AAavg_U- thor_U)) , roundorder) , r'\\')


print('HRES &', np.round( np.average(HRES_T)- np.average(thor_T) , roundorder) ,
'&', np.round( np.average(np.abs(HRES_T- thor_T)) , roundorder) ,
#'&', np.round( 1E3* (np.average(HRES_SH)- np.average(thor_SH)) , roundorder) ,
#'&', np.round(  1E3*np.average(np.abs(HRES_SH- thor_SH)) , roundorder) , r'\\')
'&', np.round(  (np.average(SH2RH(HRES_SH*1E3, plevel, HRES_T))- np.average(SH2RH(thor_SH*1E3, plevel, thor_T))) , 1) ,
'&', np.round(  np.average(np.abs(SH2RH(HRES_SH*1E3, plevel, HRES_T)- SH2RH(thor_SH*1E3, plevel, thor_T))) , 1) ,
'&', np.round( np.average(HRES_U)- np.average(thor_U) , roundorder) ,
'&', np.round( np.average(np.abs(HRES_U- thor_U)) , roundorder) , r'\\')

print('ERA-5 &', np.round( np.average(E5_T)- np.average(thor_T) , roundorder) ,
'&', np.round( np.average(np.abs(E5_T- thor_T)) , roundorder) ,
#'&', np.round( 1E3* (np.average(E5_SH)- np.average(thor_SH)) , roundorder) ,
#'&', np.round(  1E3*np.average(np.abs(E5_SH- thor_SH)) , roundorder) , r'\\')
'&', np.round(  (np.average(SH2RH(E5_SH*1E3, plevel, E5_T))- np.average(SH2RH(thor_SH*1E3, plevel, thor_T))) , 1) ,
'&', np.round(  np.average(np.abs(SH2RH(E5_SH*1E3, plevel, E5_T)- SH2RH(thor_SH*1E3, plevel, thor_T))) , 1) ,
'&', np.round( np.average(E5_U)- np.average(thor_U) , roundorder) ,
'&', np.round( np.average(np.abs(E5_U- thor_U)) , roundorder) , r'\\')

print('ERA-I &', np.round( np.average(EI_T)- np.average(thor_T) , roundorder) ,
'&', np.round( np.average(np.abs(EI_T- thor_T)) , roundorder) ,
#'&', np.round( 1E3* (np.average(EI_SH)- np.average(thor_SH)) , roundorder) ,
#'&', np.round(  1E3*np.average(np.abs(EI_SH- thor_SH)) , roundorder) , r'\\')
'&', np.round(  (np.average(SH2RH(EI_SH*1E3, plevel, EI_T))- np.average(SH2RH(thor_SH*1E3, plevel, thor_T))) , 1) ,
'&', np.round(  np.average(np.abs(SH2RH(EI_SH*1E3, plevel, EI_T)- SH2RH(thor_SH*1E3, plevel, thor_T))) , 1) ,
'&', np.round( np.average(EI_U)- np.average(thor_U) , roundorder) ,
'&', np.round( np.average(np.abs(EI_U- thor_U)) , roundorder) , r'\\')


print('ASR &', np.round( np.average(ASR_T)- np.average(thor_T) , roundorder) ,
'&', np.round( np.average(np.abs(ASR_T- thor_T)) , roundorder) ,
#'&', np.round( 1E3* (np.average(ASR_SH)- np.average(thor_SH)) , roundorder) ,
#'&', np.round(  1E3*np.average(np.abs(ASR_SH- thor_SH)) , roundorder) , r'\\')
'&', np.round(  (np.average(ASR_RH)- np.average(SH2RH(thor_SH*1E3, plevel, thor_T))) , 1) ,
'&', np.round(  np.average(np.abs(ASR_RH - SH2RH(thor_SH*1E3, plevel, thor_T))) , 1) ,
'&', np.round( np.average(ASR_U)- np.average(thor_U) , roundorder) ,
'&', np.round( np.average(np.abs(ASR_U- thor_U)) , roundorder) , r'\\')


"""print statistics"""
#print('AA_BIAS_T', np.round( np.average(AA_T)- np.average(thor_T) , 4) )
#print('AA_MAE_T', np.round( np.average(np.abs(AA_T- thor_T)) , 4) )
#print('AA_BIAS_U', np.round( np.average(AA_U)- np.average(thor_U) , 4) )
#print('AA_MAE_U', np.round( np.average(np.abs(AA_U- thor_U)) , 4) )
#print('AA_BIAS_SH', np.round( 1E3* (np.average(AA_SH)- np.average(thor_SH)) , 4) )
#print('AA_MAE_SH', np.round(  1E3*np.average(np.abs(AA_SH- thor_SH)) , 4) )
#print('AA_BIAS_Geo', np.round( np.average(AA_Geo)- np.average(thor_Geo) , 4) )
#print('AA_MAE_Geo', np.round(  np.average(np.abs(AA_Geo- thor_Geo)) , 4) )
#
#print('EI_BIAS_T', np.round( np.average(EI_T)- np.average(thor_T) , 4) )
#print('EI_MAE_T', np.round( np.average(np.abs(EI_T- thor_T)) , 4) )
#print('EI_BIAS_U', np.round( np.average(EI_U)- np.average(thor_U) , 4) )
#print('EI_MAE_U', np.round( np.average(np.abs(EI_U- thor_U)) , 4) )
#print('EI_BIAS_SH', np.round( 1E3* (np.average(EI_SH)- np.average(thor_SH)) , 4) )
#print('EI_MAE_SH', np.round(  1E3*np.average(np.abs(EI_SH- thor_SH)) , 4) )
#print('EI_BIAS_Geo', np.round( np.average(EI_Geo)- np.average(thor_Geo) , 4) )
#print('EI_MAE_Geo', np.round(  np.average(np.abs(EI_Geo- thor_Geo)) , 4) )
#
#print('E5_BIAS_T', np.round( np.average(E5_T)- np.average(thor_T) , 4) )
#print('E5_MAE_T', np.round( np.average(np.abs(E5_T- thor_T)) , 4) )
#print('E5_BIAS_U', np.round( np.average(E5_U)- np.average(thor_U) , 4) )
#print('E5_MAE_U', np.round( np.average(np.abs(E5_U- thor_U)) , 4) )
#print('E5_BIAS_SH', np.round( 1E3* (np.average(E5_SH)- np.average(thor_SH)) , 4) )
#print('E5_MAE_SH', np.round(  1E3*np.average(np.abs(E5_SH- thor_SH)) , 4) )
#
#print('HRES_BIAS_T', np.round( np.average(HRES_T)- np.average(thor_T) , 4) )
#print('HRES_MAE_T', np.round( np.average(np.abs(HRES_T- thor_T)) , 4) )
#print('HRES_BIAS_U', np.round( np.average(HRES_U)- np.average(thor_U) , 4) )
#print('HRES_MAE_U', np.round( np.average(np.abs(HRES_U- thor_U)) , 4) )
#print('HRES_BIAS_SH', np.round( 1E3* (np.average(HRES_SH)- np.average(thor_SH)) , 4) )
#print('HRES_MAE_SH', np.round(  1E3*np.average(np.abs(HRES_SH- thor_SH)) , 4) )


