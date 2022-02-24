#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

"""

import os
user = os.getcwd().split('/')[2]

Mediadir= '/media/'+user+'/PatsOrange/'


from netCDF4 import Dataset, num2date
import numpy as np
import time


class dataASR:
    def __init__(self, year, month, sday, eday=None, sh3=0, eh3=None, level='surf', reso= 1, filedir=Mediadir +'/ASR/'):
        """import data from ARS
        level= 'surf', 850, 700, 500 or 'all'. 'all' imports all variables from all fields, otherwise it just imports the grid and the time
        sday - startday, eday - endday (if None: same as sday), data in the range(sday, eday+1) is downloaded
        sh3 - start 3-hourly, eh3 - end 3-hourly (default: whole days) - data in the range(sh3, eh3+1) is downloaded
        reso -specifies the downscaling factor of the dataset: 1- takes every grid point, 2- every second in each horizontal direction"""
        self.filedir= filedir+str(year)
        self.year= year
        self.month= month
        self.sday=sday
        if eh3 is not None: eh3 += 1
        self.sh3, self.eh3= sh3, eh3

        self.reso= reso #this takes on
        if eday== None: self.eday=sday
        else: self.eday=eday
        
        if level == 'all':
            self.imp_all()
        
        else:
            if level== 'surf': filename='STOLL246610_2Dsubasr15km.anl.2D.'
            elif level == 850: filename='subasr15km.anl.3D.'
            elif level == 700: filename='STOLL246610_700_subTTasr15km.anl.3D.'
            elif level == 500: filename='STOLL246610_500_subasr15km.anl.3D.'
            
            nc= Dataset(self.filedir +r'/'+filename+str(self.year)+str(self.month).zfill(2)+str(self.sday).zfill(2)+'.nc')
    #        print(nc.variables.keys())
            #this one is even more interesting:
    #        print(nc.variables.values())
    
            self.lat = nc.variables['XLAT'][::self.reso, ::self.reso]
            self.lon = nc.variables['XLONG'][::self.reso, ::self.reso]
            
            if self.eday!= self.sday: eh3= None
            self.tim = nc.variables['Time'][self.sh3: eh3]%24 #this is something like the hour on a day??
#            self.tim0= nc.variables['Time'][self.sh3: eh3]  #this is tim for the other datasets, but it is one day to early - maybe a "Schalttag" is missing                 
                                   
            if self.eday > self.sday:
                for day in np.arange(self.sday+1, self.eday+1):
                    if day== self.eday: eh3= self.eh3 #the end hour
                    else: eh3= None
                    self.tim= np.hstack((self.tim, nc.variables['Time'][:eh3]%24 +24*(day-sday)))
#                    self.tim0= np.hstack((self.tim0, nc.variables['Time'][:eh3]))           
#        self.datetime = num2date(self.tim0, nc.variables['Time'].units)       

    #        if type(var)== list:
    #            for v in var:
    #                self.impvar(v, level)
    #        else:
    #            self.impvar(var, level)
        
    def impvar(self, var, level):
        """import the relevant data from the files
        var= string 
        level= ['surf', 850, 700, 500]"""
        
               
        if level== 'surf': filename='STOLL246610_2Dsubasr15km.anl.2D.'
        elif level == 850: filename='subasr15km.anl.3D.'
        elif level == 700: filename='STOLL246610_700_subTTasr15km.anl.3D.'
        elif level == 500: filename='STOLL246610_500_subasr15km.anl.3D.'
        
        nc= Dataset(self.filedir +r'/'+filename+str(self.year)+str(self.month).zfill(2)+str(self.sday).zfill(2)+'.nc')
#        print(var, nc.variables.keys())
        if self.eday== self.sday: eh3= self.eh3
        else: eh3= None

        if var== 'SST':
            self.SST= nc.variables['SST'][self.sh3: eh3 , ::self.reso, ::self.reso]
        elif var== 'SLP':
            self.SLP= nc.variables['PMSL'][self.sh3: eh3 , ::self.reso, ::self.reso]/100 #in hPa
        elif var == 'U' and level == 'surf':
            self.u10= nc.variables['U10M'][self.sh3: eh3 , ::self.reso, ::self.reso]
        elif var == 'V' and level == 'surf':
            self.v10= nc.variables['V10M'][self.sh3: eh3 , ::self.reso, ::self.reso]            
        elif var== 'U' and level == 850:
            self.u850= nc.variables['UU'][self.sh3: eh3, 0 , ::self.reso, ::self.reso]
        elif var == 'V' and level == 850:
            self.v850= nc.variables['VV'][self.sh3: eh3, 0 , ::self.reso, ::self.reso]
        elif var == 'vort' and level== 850:
            u= nc.variables['UU'][self.sh3: eh3, 0 , ::self.reso, ::self.reso]
            v= nc.variables['VV'][self.sh3: eh3, 0 , ::self.reso, ::self.reso]
            r,DX= 6371E3, 15E3         
           
            latr= np.deg2rad(self.lat)
            lonr= np.deg2rad(self.lon)
            dlatr= np.gradient(latr)
            dlonr= np.gradient(lonr)
            
            dlonr[0][dlonr[0] > np.pi/2] -= np.pi
            dlonr[1][dlonr[1] > np.pi/2] -= np.pi             
            dlonr[0][dlonr[0] < -np.pi/2] += np.pi
            dlonr[1][dlonr[1] < -np.pi/2] += np.pi
                 
            mx= DX/(r*np.sqrt((dlonr[0] *np.cos(latr))**2 + dlatr[0]**2))  #some mapping factor 
            u_y= np.gradient(u, DX, axis= 2)
            v_x= np.gradient(v, DX, axis= 1)
            
            self.vort=( v_x - u_y)* mx
            
        elif var == 'T' and level == 700:
            self.T700= nc.variables['TT'][self.sh3: eh3, 0 , ::self.reso, ::self.reso]            
        elif var == 'T' and level == 500:
            self.T500= nc.variables['TT'][self.sh3: eh3, 0 , ::self.reso, ::self.reso]           
        elif var == 'U' and level == 500:
            self.u500= nc.variables['UU'][self.sh3: eh3, 0 , ::self.reso, ::self.reso]
        elif var == 'V' and level == 500:
            self.v500= nc.variables['VV'][self.sh3: eh3, 0 , ::self.reso, ::self.reso]
        else:
            print('varibles is not downloaded')
            return

        """this part is to import data for several days"""
        if self.eday > self.sday:
            for day in np.arange(self.sday+1, self.eday+1):
                if day== self.eday: eh3= self.eh3 #the end hour
                else: eh3= None
                
                if var== 'SST':
                    self.SST= np.vstack((self.SST, nc.variables['SST'][: eh3, ::self.reso, ::self.reso]))
                elif var== 'SLP':
                    self.SLP= np.vstack((self.SLP, nc.variables['PMSL'][: eh3, ::self.reso, ::self.reso]/100)) #in hPa
                elif var == 'U' and level == 'surf':
                    self.u10= np.vstack((self.u10, nc.variables['U10M'][: eh3, ::self.reso, ::self.reso]))
                elif var == 'V' and level == 'surf':
                    self.v10= np.vstack((self.v10, nc.variables['V10M'][: eh3, ::self.reso, ::self.reso]))
                elif var== 'U' and level == 850:
                    self.u850= np.vstack((self.u850, nc.variables['UU'][: eh3, 0, ::self.reso, ::self.reso]))
                elif var == 'V' and level == 850:
                    self.v850= np.vstack((self.v850, nc.variables['VV'][: eh3, 0, ::self.reso, ::self.reso]))
                elif var == 'T' and level == 700:
                    self.T700= np.vstack((self.T700, nc.variables['TT'][: eh3, 0, ::self.reso, ::self.reso]))
                elif var == 'T' and level == 500:
                    self.T500= np.vstack((self.T500, nc.variables['TT'][: eh3, 0, ::self.reso, ::self.reso]))
                elif var == 'U' and level == 500:
                    self.u500= np.vstack((self.u500, nc.variables['UU'][: eh3, 0, ::self.reso, ::self.reso]))
                elif var == 'V' and level == 500:
                    self.v500= np.vstack((self.v500, nc.variables['VV'][: eh3, 0, ::self.reso, ::self.reso]))
                
                


    def imp_all(self):
        """imports all data for that day
        not really used"""
        if self.eday== self.sday: eh3= self.eh3
        else: eh3= None
        
        filename='STOLL246610_2Dsubasr15km.anl.2D.'
        nc= Dataset(self.filedir +r'/'+filename+str(self.year)+str(self.month).zfill(2)+str(self.sday).zfill(2)+'.nc')

        self.SST= nc.variables['SST'][self.sh3: eh3 , ::self.reso, ::self.reso]
        self.SLP= nc.variables['PMSL'][self.sh3: eh3 , ::self.reso, ::self.reso]/100 #in hPa
        self.u10= nc.variables['U10M'][self.sh3: eh3 , ::self.reso, ::self.reso]
        self.v10= nc.variables['V10M'][self.sh3: eh3 , ::self.reso, ::self.reso]
        self.lat = nc.variables['XLAT'][::self.reso, ::self.reso]
        self.lon = nc.variables['XLONG'][::self.reso, ::self.reso]
        self.tim = nc.variables['Time'][self.sh3: eh3]%24
#        self.tim0 = nc.variables['Time'][self.sh3: eh3]
            
#        filename='subasr15km.anl.3D.'
#        nc= Dataset(self.filedir +r'/'+filename+str(self.year)+str(self.month).zfill(2)+str(self.sday).zfill(2)+'.nc')
#        self.u850= nc.variables['UU'][:]
#        self.v850= nc.variables['VV'][:]

        filename='STOLL246610_700_subTTasr15km.anl.3D.'
        nc= Dataset(self.filedir +r'/'+filename+str(self.year)+str(self.month).zfill(2)+str(self.sday).zfill(2)+'.nc')
        self.T700= nc.variables['TT'][self.sh3: eh3, 0 , ::self.reso, ::self.reso]
        
        filename='STOLL246610_500_subasr15km.anl.3D.'    
        nc= Dataset(self.filedir +r'/'+filename+str(self.year)+str(self.month).zfill(2)+str(self.sday).zfill(2)+'.nc')
        self.T500= nc.variables['TT'][self.sh3: eh3, 0 , ::self.reso, ::self.reso]         
        self.u500= nc.variables['UU'][self.sh3: eh3, 0 , ::self.reso, ::self.reso]
        self.v500= nc.variables['VV'][self.sh3: eh3, 0 , ::self.reso, ::self.reso]
        
        if self.eday > self.sday:
            for day in np.arange(self.sday+1, self.eday+1):
                if day== self.eday: eh3= self.eh3 #the end hour
                else: eh3= None
                filename='STOLL246610_2Dsubasr15km.anl.2D.'
                nc= Dataset(self.filedir +r'/'+filename+str(self.year)+str(self.month).zfill(2)+str(day).zfill(2)+'.nc')               
                self.SST= np.vstack((self.SST, nc.variables['SST'][: eh3, ::self.reso, ::self.reso]))
                self.SLP= np.vstack((self.SLP, nc.variables['PMSL'][: eh3, ::self.reso, ::self.reso]/100)) #in hPa
                self.u10= np.vstack((self.u10, nc.variables['U10M'][: eh3, ::self.reso, ::self.reso]))
                self.v10= np.vstack((self.v10, nc.variables['V10M'][: eh3, ::self.reso, ::self.reso]))
                self.tim= np.hstack((self.tim, nc.variables['Time'][: eh3]%24 +24*(day-self.sday)))
                    
        #        filename='subasr15km.anl.3D.'
        #        nc= Dataset(self.filedir +r'/'+filename+str(self.year)+str(self.month).zfill(2)+str(day).zfill(2)+'.nc')
        #        self.u850= nc.variables['UU'][:]
        #        self.v850= nc.variables['VV'][:]
        
                filename='STOLL246610_700_subTTasr15km.anl.3D.'
                nc= Dataset(self.filedir +r'/'+filename+str(self.year)+str(self.month).zfill(2)+str(day).zfill(2)+'.nc')
                self.T700= np.vstack((self.T700, nc.variables['TT'][: eh3, 0, ::self.reso, ::self.reso]))
                
                filename='STOLL246610_500_subasr15km.anl.3D.'    
                nc= Dataset(self.filedir +r'/'+filename+str(self.year)+str(self.month).zfill(2)+str(day).zfill(2)+'.nc')
                self.T500= np.vstack((self.T500, nc.variables['TT'][: eh3, 0, ::self.reso, ::self.reso]))
                self.u500= np.vstack((self.u500, nc.variables['UU'][: eh3, 0, ::self.reso, ::self.reso]))
                self.v500= np.vstack((self.v500, nc.variables['VV'][: eh3, 0, ::self.reso, ::self.reso]))




class dataASR_fullfile:
    def __init__(self, year, month, day, filedir=Mediadir +'/ASR/asr15km.anl.3D.'):
        self.filedir= filedir+str(year)+str(month).zfill(2)+str(day).zfill(2)
