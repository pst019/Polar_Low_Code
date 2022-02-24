#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 17:15:07 2017

@author: pst019

Import the thorpex data for one pressure level or the surface
level= 'surf' or the height of the pressure level in hPa
"""
from netCDF4 import Dataset
import numpy as np
import time

level = 'surf'

direct=Mediadir+'PL/IPY-THOTPEX_DROPSONDES_3-4-MARCH_2008/QC_'


if year == 2008 and month== 3 and day== 3 and hour<15:
    
    date='080303'
    flight='a'
    times=['102731', '104912', '110458', '111154', '112742', '113414', '113922', '114338', 
    '115202', '120832', '121359', #'121815',
    '122309', '123128', '123624', '124345', '125351',
    '130301', '131613', '132830']

elif year == 2008 and month== 3 and day== 3 and hour>=15:
    date='080303'
    flight='b'
    times=['153245','154424','155225','161223','162442','163020',
    '165104','170842','171817','172626','174236']

elif year == 2008 and month== 3 and day== 4:
    date='080304'
    flight='a'
    times=['104109','104616','105431','105935','110625','111210',
    '111659','112251','113144','114148','115208', '120229', '120656',
    '121138', '121645', '122505', '123329', '124305', '124748', '125835']
    
else:
    print('no data for that time')
    exit()
    
tim=np.zeros((len(times)))
lat=np.zeros((len(times)))
lon=np.zeros((len(times)))
alt=np.zeros((len(times)))
T=np.zeros((len(times)))
pres=np.zeros((len(times)))
u=np.zeros((len(times)))
v=np.zeros((len(times)))
U=np.zeros((len(times)))
theta=np.zeros((len(times)))
theta_e=np.zeros((len(times)))
RH=np.zeros((len(times)))

for t in [0]:#range(len(times)):

    file=direct+date+flight+r'/D20'+date+'_'+times[t]+'QC.nc'
#    nc = NetCDFFile(file)
    nc= Dataset(file)

    """raw data, gets masked with the pressure"""    
    presr = nc.variables['pres'][:]

    timr = nc.variables['time'][:][~presr.mask]
    latr = nc.variables['lat'][:][~presr.mask]
    lonr = nc.variables['lon'][:][~presr.mask]
    altr= nc.variables['alt'][:][~presr.mask]
    
    # Variables
    Tr = nc.variables['tdry'][:][~presr.mask]
    ur = nc.variables['u_wind'][:][~presr.mask]
    vr = nc.variables['v_wind'][:][~presr.mask]
    Ur = nc.variables['wspd'][:][~presr.mask] #wind speed

    #Udir = nc.variables['wdir'][:]#wind direction
    thetar= nc.variables['theta'][:][~presr.mask] #potential temperature
    theta_er= nc.variables['theta_e'][:][~presr.mask] #equivalent potential temperature    
    RHr= nc.variables['rh'][:][~presr.mask] #relative humidity ???    
    
    presr= presr[~presr.mask]
    

    presr= presr[~latr.mask]
    timr= timr[~latr.mask]
    lonr= lonr[~latr.mask]
    altr= altr[~latr.mask]
    Tr= Tr[~latr.mask]
    ur= ur[~latr.mask]
    vr= vr[~latr.mask]
    Ur= Ur[~latr.mask]
    thetar= thetar[~latr.mask]
    theta_er= theta_er[~latr.mask]
    RHr= RHr[~latr.mask] 
    latr= latr[~latr.mask]
    
    if len(presr) > 1:
        """star) only put data in the array if the array is not empty"""
        
        if level=='surf':
            """condition for surface"""
            if np.min(altr) < 30: #if there is surface information for this sonde
                tim[t]= timr[np.where(altr < 30)][0]
                lat[t]= latr[np.where(altr < 30)][0]
                lon[t]= lonr[np.where(altr < 30)][0]
                pres[t]= presr[np.where(altr < 30)][0]
                U[t]= Ur[np.where(altr < 30)][0]
                alt[t]= altr[np.where(altr < 30)][0]
                theta[t]=thetar[np.where(altr < 30)][0]
                theta_e[t]=theta_er[np.where(altr < 30)][0]
                RH[t]=RHr[np.where(altr < 30)][0]

        else:
            """condition for a defined pressure level (in hPa)"""
            if np.max(presr) > level: #only if an value for the pressure level exists
                tim[t]= timr[presr > level][-1]
                lat[t]= latr[presr > level][-1]
                lon[t]= lonr[presr > level][-1]
                alt[t]= altr[presr > level][-1]
                pres[t]= presr[presr > level][-1]
                U[t]= Ur[presr > level][-1]
                theta[t]= thetar[presr > level][-1]
                theta_e[t]= theta_er[presr > level][-1]
                RH[t]= RHr[presr > level][-1]           
            
"""if the array is empty at star), then these are excluded here"""
lat=lat[tim!=0]         
lon=lon[tim!=0]
alt=alt[tim!=0]
pres=pres[tim!=0]
U=U[tim!=0]
theta=theta[tim!=0]
theta_e=theta_e[tim!=0]
RH=RH[tim!=0]
tim=tim[tim!=0]

UTC= [str(time.localtime(t)[3])+':'+str(time.localtime(t)[4]).zfill(2) for t in tim]

