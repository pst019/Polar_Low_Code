#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 12:00:39 2017

@author: pst019
"""

fignr= 1

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, '/home/'+user+'/codeAROME/')
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)


for year in range(2002, 2013):
    datadir= "/media/'+user+'/1692A00D929FEF8B/ASR/"+str(year)+r"/"
    
    #for month in range(4, 7):
    for month in range(1, 13):
    
        
        if month in [1, 3, 5, 7, 8, 10, 12]:
            daynrs= 31
        elif month in [4, 6, 9, 11]:
            daynrs= 30
        elif month == 2 and year not in [2000, 2004, 2008, 2012]:
            daynrs= 28
        else:
            daynrs= 29
        
        print(year, month, daynrs)
        
        for day in np.arange(1, daynrs+1):
            filename="subasr15km.anl.3D."+str(year)+str(month).zfill(2)+str(day).zfill(2)
            
#            print(filename)
                      
            nc = Dataset(datadir+filename+".nc")            
        #    print(nc.variables.keys())
            
            
            namlon = 'west_east'
            namlat = 'south_north'
            namlng = 'XLONG'
            namlatg= 'XLAT'
            namlev = 'lev'
            namtim = 'Time'
            namvar1 = 'UU'
            namvar2 = 'VV'
            
            namx = 'x'
            namy = 'y'
            
            
            var1= nc.variables['UU'][:]
            var2= nc.variables['VV'][:]
            
            latg= nc.variables['XLAT'][:]
            lng= nc.variables['XLONG'][:]
            tim= nc.variables['Time'][:]
            
            tim_units= nc.variables['Time'].units
            tim_calendar= nc.variables['Time'].calendar
            
            nc.close()
            
            FP_PI = 0.017453292519943295 #factor deg2rad
            plon = 95.0
            plat = 90.0
            
            
            pplon = plon * FP_PI
            pplat = plat * FP_PI
            splt = np.sin(pplat)
            cplt = np.cos(pplat)
            
            ntim, nlev, nlat, nlon = var1.shape
            
            #varlev1(:, :, :) = var1(:, :, :);
            #varlev2(:, :, :) = var2(:, :, :);
            
            
            X= np.zeros((nlon))
            
            for i in range(nlon):
                xln = lng[0, i]
                yln = latg[i, 0]
                if xln < 0.0: 
                    xln = 360.0 + xln
                xln = xln * FP_PI
                yln = yln * FP_PI
                dln = xln - pplon
                sy = np.sin(yln)
                cy = np.cos(yln)
                sdl = np.sin(dln)
                cdl = np.cos(dln)
                kk = 2.0 / (1.0 + sy * splt + cy * cplt * cdl)
                xx = cy * sdl * kk
                X[i] = (cplt * sy - splt * cy * cdl) * kk
            #    print(i, xln, yln, sy, sdl, xx, X[i])
            
            
            Y = np.zeros((nlat))
            
            for i in range(nlat):
                xln = lng[i, 0]
                yln = latg[i, 0]
                if xln < 0.0: 
                    xln = 360.0 + xln
                xln = xln * FP_PI
                yln = yln * FP_PI
                dln = xln - pplon
                sy = np.sin(yln)
                cy = np.cos(yln)
                sdl = np.sin(dln)
                cdl = np.cos(dln)
                kk = 2.0 / (1.0 + sy * splt + cy * cplt * cdl)
                Y[i] = -cy * sdl * kk
                yy = (cplt * sy - splt * cy * cdl) * kk
            
            
            
            
            
            
            """ copy nc file"""
            #import shutil
            #shutil.copyfile(datadir+filename+'.nc', datadir+filename+'_added.nc')
            #f= Dataset(datadir+filename+"_added.nc", 'r+')
            
            
            """ create nc file"""
            print(datadir+filename+'_xy.nc')
            f= Dataset(datadir+filename+"_xy.nc", 'w')
            
            
            
            f.createDimension(namx, nlon)
            f.createDimension(namy, nlat)
            f.createDimension(namtim, None)
            
            
            x = f.createVariable(namx, 'f', (namx,))
            x[:]= X
            y = f.createVariable(namy, 'f', (namy,))
            y[:]= Y
            ftim = f.createVariable(namtim, 'f', (namtim,))
            ftim[:]= tim
            ftim.units= tim_units
            ftim.calendar= tim_calendar
    #        print(ftim)
            
            LAT =  f.createVariable(namlatg, 'f', (namy, namx,))
            LAT[:]= latg
             
            VAR1 =  f.createVariable(namvar1, 'f', (namtim, namy, namx,))
            VAR1[:]= var1
            VAR2 =  f.createVariable(namvar2, 'f', (namtim, namy, namx,))
            VAR2[:]= var2
        
            
            #print(f.dimensions.keys())
        #    print(f.variables.keys())
            
            f.close()
  
    

    
    
    
    
