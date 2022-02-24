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

year, month = 2008, 3
day, hour= 4, 12  

fignr= 11

#exclude= 3, 5, 6,7( above 500hPa)
excl= [3, 5] #for the first flight
#excl=[1,9,11] #for the second flight
#excl= [5, 13]# for the third flight


thorV= Tdata(year, month, day, hour, level='vertical', plevels= np.array([1000, 950, 925, 900, 850, 800, 700, 600, 500, 400]))

idrop= tuple([i for i in range(len(thorV.dropnr)) if thorV.dropnr[i] not in excl])

thorpres= thorV.pres[idrop, :]
thorT= thorV.T[idrop, :]
thorRH= thorV.RH[idrop,:]     
thortheta= thorV.theta[idrop,:]
thorU= thorV.U[idrop, :]

    






"""get the AA data and also plot it"""
#xrange= np.arange(-6, 6.01, 3)
#yrange= np.arange(-6, 6.01, 3)

#xrange= np.arange(-30, 30.01, 15)
#yrange= np.arange(-30, 30.01, 15)

xrange= np.arange(-5, 0.01, 1)
yrange= np.arange(0, 0.01, 15)

RMSET= np.zeros((len(xrange),len(yrange)))
RMSERH= np.zeros((len(xrange),len(yrange)))
RMSEU= np.zeros((len(xrange),len(yrange)))

BIAST= np.zeros((len(xrange),len(yrange)))
BIASRH= np.zeros((len(xrange),len(yrange)))
BIASU= np.zeros((len(xrange),len(yrange)))

xoffset=5
yoffset=5

for ix, xoffset in enumerate(xrange):
    for iy, yoffset in enumerate(yrange):
        print(ix, iy)
        
        from f_imp_AROME_exp import data as Adata #import the thorpex data
        exp_name= '080303_warmctr'
        #exp_name= '080303_warmsens_noTH'
        fileday, filehour= 3, 0        
        AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
        AA= Adata(filename= AAfilename, res=1)
        
        
        AAT= np.zeros(np.shape(thorT))
        AApres= np.zeros(np.shape(thorpres))
        AARH= np.zeros(np.shape(thorT))
        AAU= np.zeros(np.shape(thorT))
        
        for di, dropind in enumerate(idrop):
            print(di)
        #    hour= thorV.datetime[dropind, 2].hour #the hour where the dropsonde reached the ground    
            hour= int(thorV.UTC[dropind][:2])
        
            t= (day- fileday)*24 + (hour- filehour)  # -1 
            
            dist= (AA.lat - thorV.lat[dropind, 2])**2 + (np.cos(np.deg2rad(thorV.lat[dropind, 2])) * (AA.lon - thorV.lon[dropind, 2])**2)
            x,y= np.where(dist == np.min(dist))
            
            AA.imp_cross_sec_reduced(xn= x[0]+ int(xoffset), yn= y[0]+ int(yoffset), tn= t)
        
            AAT[di]= AA.T[:len(thorV.plevels)]
#            AApres[di]= AA.press[:len(thorV.plevels)]
            AARH[di]= AA.RH[:len(thorV.plevels)]
            AAU[di]= np.sqrt(AA.u[:len(thorV.plevels)]**2 +AA.u[:len(thorV.plevels)]**2)
        
        
        
        
        """difference profiles"""
        RMSET[ix, iy]= np.sqrt(np.nanmean((AAT - thorT -273.15)**2))
        RMSERH[ix, iy]= np.sqrt(np.nanmean((AARH - thorRH)**2))
        RMSEU[ix, iy]= np.sqrt(np.nanmean((AAU - thorU)**2))
        
        BIAST[ix, iy]= np.nanmean(AAT - thorT -273.15)
        BIASRH[ix, iy]= np.nanmean(AARH - thorRH)
        BIASU[ix, iy]= np.nanmean(AAU - thorU)

        
plt.figure(fignr, figsize=(10, 8))
#fignr+=1
plt.clf()

plt.subplot(231)
#xplot, yplot= np.meshgrid(np.append(xrange, xrange[-1]+ xrange[1]- xrange[0]) - (xrange[1]- xrange[0])/2, np.append(yrange, yrange[-1]+ yrange[1]- yrange[0]) - (yrange[1]- yrange[0])/2)
xplot, yplot= np.meshgrid(np.append(xrange, xrange[-1]+ xrange[1]- xrange[0]) - (xrange[1]- xrange[0])/2, np.append(yrange, yrange[-1]+ xrange[1]- xrange[0]) - (xrange[1]- xrange[0])/2)

xplot *= AA.res*2.5
yplot *= AA.res*2.5
cs= plt.pcolor(xplot, yplot, RMSET.T, cmap= 'Greys')
plt.colorbar(cs)
plt.title('T')
plt.ylabel('RMSE')

plt.subplot(232)
cs= plt.pcolor(xplot, yplot, RMSERH.T, cmap= 'Greys')
plt.colorbar(cs)
plt.title('RH')

plt.subplot(233)
cs= plt.pcolor(xplot, yplot, RMSEU.T, cmap= 'Greys')
plt.colorbar(cs)
plt.title('U')

plt.subplot(234)
cs= plt.pcolor(xplot, yplot, BIAST.T, cmap= 'RdBu', vmin= -np.max(np.abs(BIAST)), vmax= np.max(np.abs(BIAST)))
plt.colorbar(cs)
plt.ylabel('BIAS')
#plt.title('T')

plt.subplot(235)
cs= plt.pcolor(xplot, yplot, BIASRH.T, cmap= 'Greys')
plt.colorbar(cs)
#plt.title('RH')

plt.subplot(236)
cs= plt.pcolor(xplot, yplot, BIASU.T, cmap= 'Greys_r')
plt.colorbar(cs)
#plt.title('U')