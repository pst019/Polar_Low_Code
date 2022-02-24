#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 17:15:07 2017

@author: pst019

Import the thorpex data for one pressure level or the surface
level= 'surf' or the height of the pressure level in hPa
"""
from netCDF4 import Dataset, num2date
import numpy as np
import time
from scipy.interpolate import griddata as gdata
from datetime import datetime as dtime

import os
user = os.getcwd().split('/')[2]

#from f_useful import *
Mediadir= '/media/'+user+'/PatsOrange/'

class data:
    def __init__(self, year, month, day, hour, level='vertical', filedir=Mediadir+'PL/IPY-THOTPEX_DROPSONDES_3-4-MARCH_2008/QC_',
                 plevels= np.arange(1000, 250, -10), dropid= 1, exclude=[], excl_level= 170):
        """import Thorpex data
        level= ['surf', 'vertical', height [hPa], 'one_profile'] to specify what kind of data is going to be imported
        'surf'= surface, 'vertical' = a vertical cross section 
        vertical- interpolates on pressure levels (plevels) ands also includes near surface variables
        one_profile - for level='one_profile' - imports one dropsonde with the number being specified by dropid (starting with 1)
        exclude -  for level='vertical' - specifies which dropid are excluded (starting with 1)
        """
        self.year= year
        self.month= month
        self.day= day
        self.hour= hour
        self.plevels= plevels
        
        date, flight, times= self.get_files()
        filedir= filedir+date+flight+r'/D20'+date

        if level is 'vertical':
            self.get_vertical(times, filedir, exclude, excl_level= excl_level)

        elif level is 'one_profile':
            self.get_one_profile(times[dropid-1], filedir)
            
        else:
            """the old default"""
            self.get_horizontal(level, times, filedir)

        
    def get_one_profile(self, times, filedir):
        file=filedir+'_'+times+'QC.nc'
        nc= Dataset(file)
    
        """raw data, gets masked with the pressure"""    
        presr = nc.variables['pres'][:]
        latr = nc.variables['lat'][:]
        RHr= nc.variables['rh'][:]/100 #relative humidity ???    
        
        mask= np.logical_and.reduce((~presr.mask, ~RHr.mask, ~latr.mask))

        if len(np.where(mask== True)[0]) == 0:
            import sys
            sys.exit("This dropsonde is masked everywhere")
            
                 
        self.tim = nc.variables['time'][:][mask]
#        print(nc.variables['time'][:][mask])
        self.datetime= num2date(nc.variables['time'][:][mask], nc.variables['time'].units)
        self.lon = nc.variables['lon'][:][mask]
        self.alt= nc.variables['alt'][:][mask]
        
        # Variables
        self.T = nc.variables['tdry'][:][mask]
        self.Td = nc.variables['dp'][:][mask]

        self.u = nc.variables['u_wind'][:][mask]
        self.v = nc.variables['v_wind'][:][mask]
        
        self.U = nc.variables['wspd'][:][mask] #wind speed
    
        #Udir = nc.variables['wdir'][:]#wind direction
        self.theta= nc.variables['theta'][:][mask] #potential temperature
        self.theta_e= nc.variables['theta_e'][:][mask] #equivalent potential temperature    
        
        self.pres= presr[mask]
        self.RH= RHr[mask]
        self.lat= latr[mask]
        

    def get_vertical(self, times, filedir, exclude, excl_level= 170):
        
        tim=np.zeros((0, len(self.plevels)))
        lat=np.zeros((0, len(self.plevels)))
        lon=np.zeros((0, len(self.plevels)))
        alt=np.zeros((0, len(self.plevels)))
        T=np.zeros((0, len(self.plevels)))
        pres=np.zeros((0, len(self.plevels)))
        u=np.zeros((0, len(self.plevels)))
        v=np.zeros((0, len(self.plevels)))
        U=np.zeros((0, len(self.plevels)))
        theta=np.zeros((0, len(self.plevels)))
        theta_e=np.zeros((0, len(self.plevels)))
        RH=np.zeros((0, len(self.plevels)))
        SH=np.zeros((0, len(self.plevels)))
        
        lon10m, lat10m, pres0m, U10m, T2m, dropnr10, dropnr, tim10m= [], [], [], [], [], [], [], []
        dropind= 0 #number of the dropsonde
        
        for t in range(len(times)):
            dropind = t+1
            if dropind in exclude:
                print('Dropsonde', dropind, 'is excluded by specification')
                continue
            
            file=filedir+'_'+times[t]+'QC.nc'
            print(dropind, file)
            nc= Dataset(file)
#            if t== 0:
#                print(nc.variables.keys())
        
            """raw data, gets masked with the pressure"""    
            presr = nc.variables['pres'][:]
            latr = nc.variables['lat'][:]
            RHr= nc.variables['rh'][:]/100 #relative humidity ???    
            
            mask= np.logical_and.reduce((~presr.mask, ~RHr.mask, ~latr.mask))   

#            print(np.where(mask== True)[0])
#            if np.where(mask== True)[0] == 0:
#                print('This dropsonde is masked everywhere')

            timr = nc.variables['time'][:][mask]           
            lonr = nc.variables['lon'][:][mask]
            altr= nc.variables['alt'][:][mask]
            
            # Variables
            SHr= nc.variables['mr'][:][mask]*1E-3 #mixing ration [g/kg] converted to [kg/kg]
            SHr= SHr/(1+SHr) #translation to in specific humidity [kg/kg]                           
            Tr = nc.variables['tdry'][:][mask]
            ur = nc.variables['u_wind'][:][mask]
            vr = nc.variables['v_wind'][:][mask]
            Ur = nc.variables['wspd'][:][mask] #wind speed
        
            #Udir = nc.variables['wdir'][:]#wind direction
            thetar= nc.variables['theta'][:][mask] #potential temperature
            theta_er= nc.variables['theta_e'][:][mask] #equivalent potential temperature    
            
            presr= presr[mask]
            RHr= RHr[mask] 
            latr= latr[mask]
            

#            if len(presr) < 100:
#                print(list(presr))
#            Ti= griddata(presr, Tr, self.plevels)

            if len(presr) == 0:
                """some data is empty after the masking - this is excluded here"""
                print('exclusion '+ str(dropind)+' due to lacking data')  #some of them have no latitude data
                #of the first flight sonde 3 has 75 values and sonde 4 has 0 values
                continue
            
            pdiff= np.array(presr)[:-1] - np.array(presr)[1:] #the pressure difference between two measurement points of a sonde
            
            if np.max(pdiff) > 100:                        
                print(str(dropind), 'has a pressure diff', np.max(pdiff), 'between', np.array(presr)[1:][pdiff == np.max(pdiff)][0], 'and', np.array(presr)[:-1][pdiff == np.max(pdiff)][0]  ) 

            if np.max(pdiff) > excl_level: #exclude this data
                print('exclusion '+ str(dropind) +' due to large data gap')  #some of them have no latitude data         
                continue
                            
            """interpolation on pressure levels"""
            tim= np.vstack((tim, gdata(presr, timr, self.plevels)))
            lat= np.vstack((lat, gdata(presr, latr, self.plevels)))
            lon= np.vstack((lon, gdata(presr, lonr, self.plevels)))
            pres= np.vstack((pres, gdata(presr, presr, self.plevels)))
            T= np.vstack((T, gdata(presr, Tr, self.plevels)))

            u= np.vstack((u, gdata(presr, ur, self.plevels)))
            v= np.vstack((v, gdata(presr, vr, self.plevels)))
            
            U= np.vstack((U, gdata(presr, Ur, self.plevels)))
            alt= np.vstack((alt, gdata(presr, altr, self.plevels)))
            theta= np.vstack((theta, gdata(presr, thetar, self.plevels)))
            theta_e=np.vstack((theta_e, gdata(presr, theta_er, self.plevels)))
            RH= np.vstack((RH, gdata(presr, RHr, self.plevels)))
            SH= np.vstack((SH, gdata(presr, SHr, self.plevels)))

            
            dropnr.append(dropind)
            """interpolation on 10m"""
#                print(altr[0], altr[1])
            if altr[0]< 50:
                
                lon10m.append(lonr[0])
                lat10m.append(latr[0])
                tim10m.append(timr[0])
                pres0m.append(presr[0]+ (presr[1]- presr[0])/(altr[1] - altr[0])*(0 - altr[0]) )
                U10m.append(Ur[0]+ (Ur[1]- Ur[0])/(altr[1] - altr[0])*(10 - altr[0]) )
                T2m.append(Tr[0]+ (Tr[1]- Tr[0])/(altr[1] - altr[0])*(2 - altr[0]) )
                dropnr10.append(dropind)
                
  
        self.lat=lat#[tim!=0]         
        self.lon=lon#[tim!=0]
        self.alt=alt#[tim!=0]
        self.pres=pres#[tim!=0]
        self.T=T
        self.u=u
        self.v=v
        self.U=U#[tim!=0]
        self.theta=theta#[tim!=0]
        self.theta_e=theta_e#[tim!=0]
        self.RH=RH#[tim!=0]
        self.SH=SH
        self.tim=tim#[tim!=0]
        tim[np.isnan(tim)]= 0
        datetime= num2date(tim, nc.variables['time'].units)
        datetime[datetime== dtime(1970, 1, 1, 0, 0)]= np.nan
        self.datetime= datetime
    
        
        #this can likely be excluded with num2date
        print('UTC should be replaced by datetime')
    #        self.UTC= self.datetime
    #        self.UTC= np.array([str(time.localtime(np.nanmax(tim[t]))[3])+':'+str(time.localtime(np.nanmax(tim[t]))[4]).zfill(2) for t in range(len(tim[:,0]))])
        self.datetime10= num2date(tim10m, nc.variables['time'].units)
        #        self.UTC10= np.array([str(time.localtime(tim10m[t])[3])+':'+str(time.localtime(tim10m[t])[4]).zfill(2) for t in range(len(tim10m))])
    
        self.lon10, self.lat10, self.pres0, self.T2, self.U10, self.dropnr, self.dropnr10= np.array(lon10m), np.array(lat10m), np.array(pres0m), np.array(T2m), np.array(U10m), np.array(dropnr), np.array(dropnr10)
        self.UTC10 = np.array([str(i)[11:16] for i in self.datetime10])
        
            
    def get_horizontal(self, level, times, filedir):           
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
        
        for t in range(len(times)):
        
            file=filedir+'_'+times[t]+'QC.nc'
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
            RHr= nc.variables['rh'][:][~presr.mask]/100 #relative humidity ???    
            
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
        self.lat=lat[tim!=0]         
        self.lon=lon[tim!=0]
        self.alt=alt[tim!=0]
        self.pres=pres[tim!=0]
        self.U=U[tim!=0]
        self.theta=theta[tim!=0]
        self.theta_e=theta_e[tim!=0]
        self.RH=RH[tim!=0]
        self.tim=tim[tim!=0]
        
        self.UTC= [str(time.localtime(t)[3])+':'+str(time.localtime(t)[4]).zfill(2) for t in self.tim]


    def get_files(self):
        if self.year == 2008 and self.month== 3 and self.day== 3 and self.hour<15:
            
            date='080303'
            flight='a'
            times=['102731', '104912', '110458', '111154', '112742', '113414', '113922', '114338', 
            '115202', '120832', '121359', '121815',
            '122309', '123128', '123624', '124345', '125351',
            '130301', '131613', '132830']
        
        elif self.year == 2008 and self.month== 3 and self.day== 3 and self.hour>=15:
            date='080303'
            flight='b'
#            times=['153245','154424','155225','161223','162442','163020',
#            '165104','170842','171817','172626','174236']

            times=['151935', '153245','154424', '155225','161223',
                   '161923', '162442','163020', '165104','170842',
                   '170858', '171817','172626','174236']
            #nr 2 and 12 are masked everywhere
        
        elif self.year == 2008 and self.month== 3 and self.day== 4:
            date='080304'
            flight='a'
            times=['104109','104616','105431','105935','110625','111210',
            '111659','112251','113144','114148','115208', '120229', '120656',
            '121138', '121645', '122505', '123329', '124305', '124748', '125835']
               
        else:
            print('no data for that time')
            exit()
        
        return date, flight, times
