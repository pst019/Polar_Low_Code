#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

This is from plot_comparison to import data
"""

from netCDF4 import Dataset, num2date
import numpy as np
import time

import os
user = os.getcwd().split('/')[2]

from f_useful import *

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'

class data:
    def __init__(self, var, year, month, filedir=Mediadir +'/ERA/',
                 tstart=0, tend=None, hemisp='NH'):
        """import data from ERA interim
        var= string or a list of strings, giving the names of the variables
        tstart and tend gives the start and the end point in the timsteps of the ERA data
        hemisp - 'NH' or 'SH' for northern or southern"""
        if hemisp == 'SH': filedir += 'SouthernH/'
        self.filedir= filedir
        self.year= year
        self.month= month
        self.tstart= tstart
        
        if type(var)== list:
            v = var[0]
        else: v= var
        nc= Dataset(self.filedir+v+r'/'+v+'.'+str(year)+'.'+str(month).zfill(2)+'.nc')

#        print(nc.variables.keys())
        #this one is even more interesting:
#        print(nc.variables.values())
        """should remove height above msl"""
        
        """reduce resolution"""
        if tend ==None:
            self.tim = nc.variables['time'][tstart:]

            self.tend= len(self.tim)
        else:
            self.tim = nc.variables['time'][tstart:tend]
            self.tend=tend
        
        self.datetime = num2date(self.tim, nc.variables['time'].units)       
        self.lat = nc.variables['lat'][1:]
        self.lon = nc.variables['lon'][:]
#        self.press= nc.variables['lev'][:]/100 #hPa
        
        if type(var)== list:
            for v in var:
                self.impvar(v)
        else:
            self.impvar(var)
        
    def impvar(self, var):
        """import the relevant data from the files"""
        
        nc= Dataset(self.filedir+var+r'/'+var+'.'+str(self.year)+'.'+str(self.month).zfill(2)+'.nc')
#        print(var, nc.variables.keys())
#        if 'lev' in nc.variables.keys():
#            print('At levels: '+ str(nc.variables['lev'][:]) +'(Pa or PVU)')
        # Variables
        if var== 'Vort':
            self.vort= nc.variables['var138'][self.tstart:self.tend, 0, 1:, :]*1E5
        elif var== 'Vort_fc':
            self.vort_fc= nc.variables['var138'][self.tstart:self.tend, 0, 1:, :]*1E5
            self.tim_fc = nc.variables['time'][self.tstart:self.tend]
        elif var== 'Vort_128':
            self.vort128= nc.variables['var138'][self.tstart:self.tend, 0, :, :]*1E5
            self.lat128= nc.variables['lat'][:]
            self.lon128= nc.variables['lon'][:]
        elif var=='SST':
            self.SST= nc.variables['var34'][self.tstart:self.tend, 1:, :] # (time, lat, lon)
        elif var=='T':
            l=[0,1,2]
#            print('T at plev: '+ str(nc.variables['lev'][l]/100))
            self.T = nc.variables['var130'][self.tstart:self.tend, l, 1:, :] # (time, lat, lon)
        elif var=='T850':
            self.T850 = nc.variables['var130'][self.tstart:self.tend, 0, 1:, :] 
        elif var=='MSLP':
            self.MSLP= nc.variables['var151'][self.tstart:self.tend, 1:, :]/100 #hPa
        elif var=='uPV':
            self.uPV= nc.variables['var131'][self.tstart:self.tend,0 , 1:, :]
        elif var=='vPV':
            self.vPV= nc.variables['var132'][self.tstart:self.tend,0, 1:, :]
        elif var=='Geop':
#            print(nc.variables['lev'][:])
            self.Geop= nc.variables['var129'][self.tstart:self.tend,:, 1:, :]
        elif var=='Geop500':
            self.Geop500= nc.variables['var129'][self.tstart:self.tend, 0, 1:, :]            
        elif var=='uPL':
            self.uPL= nc.variables['var131'][self.tstart:self.tend,: , 1:, :] #at 925 and 700 hPa (for shear)
#            print('uPL', nc.variables['lev'][:])
        elif var=='vPL':
            self.vPL= nc.variables['var132'][self.tstart:self.tend,:, 1:, :]
        elif var=='uPL500':
            self.uPL500= nc.variables['var131'][self.tstart:self.tend,0 , 1:, :]
        elif var=='vPL500':
            self.vPL500= nc.variables['var132'][self.tstart:self.tend,0, 1:, :]            
        elif var=='Water':
            #total column water kg/m**2
            self.Water= nc.variables['var136'][self.tstart:self.tend, 1:, :]
        elif var=='p_PV':
            self.p_PV= nc.variables['var54'][self.tstart:self.tend, 0, 1:, :]
        elif var=='TotCloud':
            self.TotCloud= nc.variables['var164'][self.tstart:self.tend, 1:, :]
        elif var=='HighCloud':
            self.HighCloud= nc.variables['var188'][self.tstart:self.tend, 1:, :]            
        elif var=='MedCloud':
            self.MedCloud= nc.variables['var187'][self.tstart:self.tend, 1:, :]            
        elif var=='LowCloud':
            self.LowCloud= nc.variables['var186'][self.tstart:self.tend, 1:, :]
        elif var=='CAPE':
            self.CAPE= nc.variables['var59'][self.tstart:self.tend, 1:, :]            
        elif var=='sHum':
            l=[0,1,2]
#            print('sHum at plev: '+ str(nc.variables['lev'][l]/100))
            self.sHum= nc.variables['var133'][self.tstart:self.tend,l, 1:, :]   #kg H2O / kg air
        elif var=='sHum850':
            self.sHum850= nc.variables['var133'][self.tstart:self.tend, 0, 1:, :]
        elif var=='PBH':
            self.PBH= nc.variables['var159'][self.tstart:self.tend, 1:, :]
        elif var=='thetaPV':
            self.thetaPV= nc.variables['var3'][self.tstart:self.tend,0, 1:, :]
   
    def impU10(self):
        """import the wind speed at 10m and interpolate it on the same grid as the other variables"""
        from scipy import interpolate
        nc= Dataset(self.filedir+'U10/U10.'+str(self.year)+'.'+str(self.month).zfill(2)+'.nc')
#        print(nc.variables.keys())
        lat75 = nc.variables['lat'][1:]
        lon75 = nc.variables['lon'][:]
        u10_75= nc.variables['var165'][self.tstart:self.tend, 1:, :] # (time, lat, lon)
        
        nc= Dataset(self.filedir+'V10/V10.'+str(self.year)+'.'+str(self.month).zfill(2)+'.nc')
#        print(nc.variables.keys())
        v10_75= nc.variables['var166'][self.tstart:self.tend, 1:, :] # (time, lat, lon)
        
        U10_75= np.sqrt(v10_75**2 + u10_75**2)
        
        U10= np.zeros((len(self.tim), len(self.lat), len(self.lon)))
        for t in range(len(self.tim)):
            f = interpolate.interp2d(lon75, lat75, U10_75[t], kind='cubic')
            U10[t]= f(self.lon, self.lat)[::-1]
        self.U10 = U10
            
           
    def imp_u10(self):
        """import the wind speed at 10m and interpolate it on the same grid as the other variables"""
        from scipy import interpolate
        nc= Dataset(self.filedir+'U10/U10.'+str(self.year)+'.'+str(self.month).zfill(2)+'.nc')
        #        print(nc.variables.keys())
        lat75 = nc.variables['lat'][1:]
        lon75 = nc.variables['lon'][:]
        u10_75= nc.variables['var165'][self.tstart:self.tend, 1:, :] # (time, lat, lon)
        
        nc= Dataset(self.filedir+'V10/V10.'+str(self.year)+'.'+str(self.month).zfill(2)+'.nc')
        #        print(nc.variables.keys())
        v10_75= nc.variables['var166'][self.tstart:self.tend, 1:, :] # (time, lat, lon)
        
        u10= np.zeros((len(self.tim), len(self.lat), len(self.lon)))
        v10= np.zeros((len(self.tim), len(self.lat), len(self.lon)))
        
        for t in range(len(self.tim)):
            f = interpolate.interp2d(lon75, lat75, u10_75[t], kind='cubic')
            u10[t]= f(self.lon, self.lat)[::-1]
            f = interpolate.interp2d(lon75, lat75, v10_75[t], kind='cubic')
            v10[t]= f(self.lon, self.lat)[::-1]
            
        self.u10= u10
        self.v10= v10           



class vort_128:
    def __init__(self, var, year, season, filedir=Mediadir +'ERA/' ):
        """import the vort_128_comb data from ERA interim
        var= string or a list of strings, giving the names of the variables"""
        print('the data is just converted for 2010 A.')
        self.filedir= filedir
        self.year= year
        
        nc= Dataset(self.filedir+var+r'/'+var+'.'+str(year)+'.'+season+'.nc')
#        print(nc.variables.keys())
        #this one is even more interesting:
#        print(nc.variables.values())
        """should remove height above msl"""
        
        """reduce resolution"""       
        self.tim = nc.variables['time'][:]
        self.lon = nc.variables['lon'][:]
        
        if var=='Vort_128_comb_filt':
            self.lat = nc.variables['lat'][::-1]
            self.vort= nc.variables['var'][:, 0, ::-1, :]*1E5
        else:
            self.lat = nc.variables['lat'][1:]
            self.vort= nc.variables['var138'][:, 0, 1:, :]*1E5                          
        

        
#    def impvar(self, var):
#        """import the relevant data from the files"""
#        
#        nc= Dataset(self.filedir+var+r'/'+var+'.'+str(self.year)+'.'+str(self.month).zfill(2)+'.nc')
#        print(var, nc.variables.keys())
#        if 'lev' in nc.variables.keys():
#            print('At levels: '+ str(nc.variables['lev'][:]) +'(Pa or PVU)')
#        # Variables
#        if var== 'Vort_128_comb':
#            self.vort= nc.variables['var138'][:, 0, 1:, :]*1E5
##        elif var== 'Vort_fc':
##            self.vort_fc= nc.variables['var138'][:, 0, 1:, :]*1E5
##            self.tim_fc = nc.variables['time'][:]



def Vort_comb(vort, vort_fc, tim, tim_fc):
    "combine the data from the analysis and the forecast to get 3 hourly data"""
    nt, nlat, nlon= np.shape(vort)
    vort_comb= np.zeros((nt*2, nlat, nlon))
    vort_comb[::2]= vort
    vort_comb[1::2]= vort_fc
    
    tim_comb= np.zeros((nt*2))
    tim_comb[::2]= tim
    tim_comb[1::2]= tim_fc
    return vort_comb, tim_comb
    

