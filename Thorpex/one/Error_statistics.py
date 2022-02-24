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
from f_imp_ASR import dataASR as ASRdata


from datetime import datetime, timedelta
import numpy as np


Mediadir= '/media/'+user+'/1692A00D929FEF8B/'




year, month = 2008, 3
flight= 3
#plevel= 950 #no ERAI
#plevel= 850 #no ERAI wind
plevel=700 #everything + geop
#plevel= 500 #everything
#plevel= 925 #only wind + geop for AA and ERAI (avail for ERA5)


"""get thorpex data for this level"""
exclude_dropsonde=[]
if flight == 1:
    day, hour= 3, 12
    if plevel in [925, 950]: exclude_dropsonde= [12]
#    excl= [3, 5] #for the first flight
if flight == 2:
    day, hour= 3, 18  
#    excl=[1,9,11] #for the second flight
if flight == 3:
    day, hour= 4, 12
    if plevel==950: exclude_dropsonde= [1]
#    excl= [5, 13]# for the third flight - exclude dropsonde 1 for 950hPa


thor= Tdata(year, month, day, hour, level='vertical', plevels= [plevel], exclude=exclude_dropsonde)

thor_T= thor.T[:,0]
thor_T+= 273.15

thor_U= thor.U[:,0]
thor_SH= thor.SH[:,0]
thor_Geo= thor.alt[:,0]

"""get AROME data"""
exp_name= '080303_warmctr'
#exp_name= '080303_cold_pseudo2'
#exp_name= '08030312_cycling'
exp_name= '080304_cold_pseudo'
#exp_name= '080303_warmsens_noTH'
fileday, filehour= 4, 0

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

AA_T, AA_U, AA_SH, AA_Geo = np.zeros(np.shape(thor_T)), np.zeros(np.shape(thor_T)), np.zeros(np.shape(thor_T)), np.zeros(np.shape(thor_T))
AAavg_T, AAavg_U, AAavg_SH= np.zeros(np.shape(thor_T)), np.zeros(np.shape(thor_T)), np.zeros(np.shape(thor_T))


for di in range(len(thor.dropnr)):
    
    tn= int(np.round((thor.datetime[di,0] -filedatetime )/ timedelta(hours=1)))
    if exp_name in ['08030312_cycling']: tn= int(np.round(tn/3))

    print(tn)
    AA.imp_level(pn= pn, tn= tn)

       
    dist= 110* np.sqrt((AA.lat - thor.lat[di, 0])**2 + (np.cos(np.deg2rad(thor.lat[di, 0])) * (AA.lon - thor.lon[di, 0]) )**2 )
    x,y= np.where(dist == np.min(dist))
        
    AA_T[di]= AA.T[x, y]
    AA_U[di]= np.sqrt(AA.u[x,y]**2 + AA.v[x,y]**2)
    AA_SH[di]= AA.SH[x,y]
    AA_Geo[di]= AA.geop[x,y]
    
    x,y= np.where(dist <= 10)
    print('number cells', len(x))
    
    AAavg_T[di]= np.average(AA.T[x, y])
    AAavg_U[di]= np.average(np.sqrt(AA.u[x,y]**2 + AA.v[x,y]**2))
    AAavg_SH[di]= np.average(AA.SH[x,y])
       


"""get ERA-I"""
t= (day-1)*4 + hour//6

EI_T, EI_U, EI_SH, EI_Geo = np.zeros(np.shape(thor_T)), np.zeros(np.shape(thor_T)), np.zeros(np.shape(thor_T)), np.zeros(np.shape(thor_T))
 
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

EIlon, EIlat = np.meshgrid(EI.lon, EI.lat)
for di in range(len(thor.dropnr)):
    dist= (EIlat - thor.lat[di, 0])**2 + (np.cos(np.deg2rad(thor.lat[di, 0])) * (EIlon - thor.lon[di, 0])**2)
    x,y= np.where(dist == np.min(dist))

    EI_T[di]= EIT[x, y]
    EI_U[di]= EIU[x,y]
    EI_SH[di]= EISH[x,y]
    EI_Geo[di]= EIGeo[x,y]

"""get ERA5"""
E5_T, E5_U, E5_SH = np.zeros(np.shape(thor_T)), np.zeros(np.shape(thor_T)), np.zeros(np.shape(thor_T))


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
    
    tn= int(np.round((thor.datetime[di,0] - datetime(year, month, day, 0, 0))/ timedelta(hours=1)))

    print(tn)
    dist= (E5lat - thor.lat[di, 0])**2 + (np.cos(np.deg2rad(thor.lat[di, 0])) * (E5lon - thor.lon[di, 0])**2)
    x,y= np.where(dist == np.min(dist))
        
    E5_T[di]= E5T[tn, x, y]
    E5_U[di]= E5U[tn,x,y]
    E5_SH[di]= E5SH[tn,x,y]
    

"""get HRES"""
HRES_T, HRES_U, HRES_SH = np.zeros(np.shape(thor_T)), np.zeros(np.shape(thor_T)), np.zeros(np.shape(thor_T))


varimp= ['T', 'SH', 'U']

if plevel== 850: pn= 1
elif plevel== 700: pn= 2
elif plevel== 500: pn= 3
elif plevel== 950: pn= 0

file= Mediadir + 'ECMWF/HRES/oper_PL_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+'_'+str(filehour).zfill(2)+'.nc'
nc= Dataset(file)

tim= nc.variables['time'][:]
HRESdatetime= num2date(tim, nc.variables['time'].units)
lat= nc.variables['latitude'][:]
lon= nc.variables['longitude'][:]        
HRESlevel= nc.variables['level'][:]

HRESlon, HRESlat = np.meshgrid(lon, lat)

HREST= nc.variables['t'][:,pn]
HRESSH= nc.variables['q'][:,pn]
HRESU= np.sqrt(nc.variables['u'][:,pn]**2+ nc.variables['v'][:, pn]**2)


for di in range(len(thor.dropnr)):
    
    tn= int(np.round((thor.datetime[di,0] - datetime(year, month, day, 0, 0))/ timedelta(hours=3)))

    print(tn)
    dist= (HRESlat - thor.lat[di, 0])**2 + (np.cos(np.deg2rad(thor.lat[di, 0])) * (HRESlon - thor.lon[di, 0])**2)
    x,y= np.where(dist == np.min(dist))
        
    HRES_T[di]= HREST[tn, x, y]
    HRES_U[di]= HRESU[tn,x,y]
    HRES_SH[di]= HRESSH[tn,x,y]

    
    
"""print statistics"""  
print('AROME &', np.round( np.average(AA_T)- np.average(thor_T) , 4) ,
'&', np.round( np.average(np.abs(AA_T- thor_T)) , 4) ,
'&', np.round( np.average(AA_U)- np.average(thor_U) , 4) ,
'&', np.round( np.average(np.abs(AA_U- thor_U)) , 4) ,
'&', np.round( 1E3* (np.average(AA_SH)- np.average(thor_SH)) , 4) ,
'&', np.round(  1E3*np.average(np.abs(AA_SH- thor_SH)) , 4) , r'\\' )

print('AROME avg &', np.round( np.average(AAavg_T)- np.average(thor_T) , 4) ,
'&', np.round( np.average(np.abs(AAavg_T- thor_T)) , 4) ,
'&', np.round( np.average(AAavg_U)- np.average(thor_U) , 4) ,
'&', np.round( np.average(np.abs(AAavg_U- thor_U)) , 4) ,
'&', np.round( 1E3* (np.average(AAavg_SH)- np.average(thor_SH)) , 4) ,
'&', np.round(  1E3*np.average(np.abs(AAavg_SH- thor_SH)) , 4) , r'\\' )

print('ERA-I &', np.round( np.average(EI_T)- np.average(thor_T) , 4) ,
'&', np.round( np.average(np.abs(EI_T- thor_T)) , 4) ,
'&', np.round( np.average(EI_U)- np.average(thor_U) , 4) ,
'&', np.round( np.average(np.abs(EI_U- thor_U)) , 4) ,
'&', np.round( 1E3* (np.average(EI_SH)- np.average(thor_SH)) , 4) ,
'&', np.round(  1E3*np.average(np.abs(EI_SH- thor_SH)) , 4) , r'\\')

print('ERA5 &', np.round( np.average(E5_T)- np.average(thor_T) , 4) ,
'&', np.round( np.average(np.abs(E5_T- thor_T)) , 4) ,
'&', np.round( np.average(E5_U)- np.average(thor_U) , 4) ,
'&', np.round( np.average(np.abs(E5_U- thor_U)) , 4) ,
'&', np.round( 1E3* (np.average(E5_SH)- np.average(thor_SH)) , 4) ,
'&', np.round(  1E3*np.average(np.abs(E5_SH- thor_SH)) , 4) , r'\\')

print('HRES &', np.round( np.average(HRES_T)- np.average(thor_T) , 4) ,
'&', np.round( np.average(np.abs(HRES_T- thor_T)) , 4) ,
'&', np.round( np.average(HRES_U)- np.average(thor_U) , 4) ,
'&', np.round( np.average(np.abs(HRES_U- thor_U)) , 4) ,
'&', np.round( 1E3* (np.average(HRES_SH)- np.average(thor_SH)) , 4) ,
'&', np.round(  1E3*np.average(np.abs(HRES_SH- thor_SH)) , 4) , r'\\')
    
"""print statistics"""
print('AA_BIAS_T', np.round( np.average(AA_T)- np.average(thor_T) , 4) )
print('AA_MAE_T', np.round( np.average(np.abs(AA_T- thor_T)) , 4) )
print('AA_BIAS_U', np.round( np.average(AA_U)- np.average(thor_U) , 4) )
print('AA_MAE_U', np.round( np.average(np.abs(AA_U- thor_U)) , 4) )
print('AA_BIAS_SH', np.round( 1E3* (np.average(AA_SH)- np.average(thor_SH)) , 4) )
print('AA_MAE_SH', np.round(  1E3*np.average(np.abs(AA_SH- thor_SH)) , 4) )
print('AA_BIAS_Geo', np.round( np.average(AA_Geo)- np.average(thor_Geo) , 4) )
print('AA_MAE_Geo', np.round(  np.average(np.abs(AA_Geo- thor_Geo)) , 4) )

print('EI_BIAS_T', np.round( np.average(EI_T)- np.average(thor_T) , 4) )
print('EI_MAE_T', np.round( np.average(np.abs(EI_T- thor_T)) , 4) )
print('EI_BIAS_U', np.round( np.average(EI_U)- np.average(thor_U) , 4) )
print('EI_MAE_U', np.round( np.average(np.abs(EI_U- thor_U)) , 4) )
print('EI_BIAS_SH', np.round( 1E3* (np.average(EI_SH)- np.average(thor_SH)) , 4) )
print('EI_MAE_SH', np.round(  1E3*np.average(np.abs(EI_SH- thor_SH)) , 4) )
print('EI_BIAS_Geo', np.round( np.average(EI_Geo)- np.average(thor_Geo) , 4) )
print('EI_MAE_Geo', np.round(  np.average(np.abs(EI_Geo- thor_Geo)) , 4) )

print('E5_BIAS_T', np.round( np.average(E5_T)- np.average(thor_T) , 4) )
print('E5_MAE_T', np.round( np.average(np.abs(E5_T- thor_T)) , 4) )
print('E5_BIAS_U', np.round( np.average(E5_U)- np.average(thor_U) , 4) )
print('E5_MAE_U', np.round( np.average(np.abs(E5_U- thor_U)) , 4) )
print('E5_BIAS_SH', np.round( 1E3* (np.average(E5_SH)- np.average(thor_SH)) , 4) )
print('E5_MAE_SH', np.round(  1E3*np.average(np.abs(E5_SH- thor_SH)) , 4) )

print('HRES_BIAS_T', np.round( np.average(HRES_T)- np.average(thor_T) , 4) )
print('HRES_MAE_T', np.round( np.average(np.abs(HRES_T- thor_T)) , 4) )
print('HRES_BIAS_U', np.round( np.average(HRES_U)- np.average(thor_U) , 4) )
print('HRES_MAE_U', np.round( np.average(np.abs(HRES_U- thor_U)) , 4) )
print('HRES_BIAS_SH', np.round( 1E3* (np.average(HRES_SH)- np.average(thor_SH)) , 4) )
print('HRES_MAE_SH', np.round(  1E3*np.average(np.abs(HRES_SH- thor_SH)) , 4) )


