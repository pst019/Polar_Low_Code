#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 15:06:41 2017

@author: pst019
"""

from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import sys  #to import the functions from a different directory
import pickle

#user='patricks'
user='pst019'

sys.path.insert(0, '/home/'+user+'/codeDerive_PL_clim/ERA/')
from f_imp_ERA2 import *
from f_imp_ASR import *

from f_impLists import *
sys.path.insert(0, '/home/'+user+'/codeAROME/')
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)

import scipy.ndimage.filters as filters

import csv

import datetime
import time
start_time = time.time()




fignr= 10

Tdurationlim=1 #TRACK duration limit in 3hourly time steps - higher less PLs
vortlim= 4.0 #the higher the less PLs
thetastablim= np.float64(-13.0) #the higher the less PLs

windlim= np.float64(17.5) #the higher the less PLs
polarfrontlim= np.float64(29.8) #the higher the more PLs

RCp= 2/7 #R/Cp

cnr, vnr, lnr, pnr, jnr, snr, wnr= 0,0,0,0,0,0,0

for year in np.arange(2000, 2013):
    cPLlist= np.zeros((11,0)) #complete PLlist (year, month, TrackPLnr, tPL, lon, lat, vorticity, 10mwind, stability, U500north, PLpoint(yes= 1, no= 0))

    
    for month in np.arange(1, 13):
        print(year, month)

        shearlist, stablist, U500northlist= [], [], []

        T=TRACK_list_month2(year, month, duration_limit= Tdurationlim, timestep='threehours', model='ASR')

        cnr += len(remove_dublicate(T.PLlist[0]))
        print('before exclusion', cnr)

        """vorticity exclusion"""
        PLlist= T.PLlist[:, T.PLlist[-1] > vortlim]
        vnr += len(remove_dublicate(PLlist[0]))
        print('after vorticity exclusion', vnr)
        
        """land and ice exclusion"""
        landlist= np.zeros((len(PLlist[0])))
        
        for day in range(1, daynr_month(year, month)+1):
            """import data"""
            d= dataASR(year, month, sday=day, eday=day)
            d.impvar('SST', level='surf')
            d.SST[d.SST<272]= np.nan #this sets the area with ice to land
            d.SST= ma.array(d.SST, mask= isnan(d.SST))              
            
            index= np.where(np.logical_and(PLlist[1] > (day-1)*8, PLlist[1]< day*8))[0] #indexes where PLlist[1] is on day
            for i in index:
                minval= np.min((PLlist[3][i] - d.lat)**2 + (PLlist[2][i] - d.lon)**2)
                pts2D= np.where((PLlist[3][i] - d.lat)**2 + (PLlist[2][i] - d.lon)**2 == minval)  
                
                t = int(PLlist[1][i]%8)

                dist= 4*55/15
                mask2= radMask2D(d.SST.shape[1:], (int(pts2D[0]), int(pts2D[1])), radiusx= dist, radiusy= dist) #lat lon mask
                maskc2= [t, mask2] #tim, lat, lon mask (mask complete)
                landlist[i]= len(np.where(d.SST[maskc2].mask == True)[0]) /len(d.SST[maskc2])

#        PLlist= np.vstack((PLlist, Ulist))
        PLlist = PLlist[:, landlist < .25]

        lnr += len(remove_dublicate(PLlist[0]))
        print('after land exclusion', lnr)

        """wind exclusion"""
        Ulist= np.zeros((len(PLlist[0])))
        for day in range(1, daynr_month(year, month)+1):
            """import data"""
            d= dataASR(year, month, sday=day, eday=day)
            d.impvar('U', level='surf')
            d.impvar('V', level='surf')
            U10= np.sqrt(d.u10**2 + d.v10**2)
            
            index= np.where(np.logical_and(PLlist[1] > (day-1)*8, PLlist[1]< day*8))[0] #indexes where PLlist[1] is on day
            for i in index:
                minval= np.min((PLlist[3][i] - d.lat)**2 + (PLlist[2][i] - d.lon)**2)
                pts2D= np.where((PLlist[3][i] - d.lat)**2 + (PLlist[2][i] - d.lon)**2 == minval)  
                
                t = int(PLlist[1][i]%8)

                dist= 4*55/15
                mask2= radMask2D(d.u10.shape[1:], (int(pts2D[0]), int(pts2D[1])), radiusx= dist, radiusy= dist) #lat lon mask
                maskc2= [t, mask2] #tim, lat, lon mask (mask complete)
                Ulist[i]= np.max(U10[maskc2])

        PLlist= np.vstack((PLlist, Ulist))
        PLlist = PLlist[:, Ulist > windlim]

        wnr += len(remove_dublicate(PLlist[0]))
        print('after wind exclusion', wnr)

        """stab exclusion"""
        stablist= np.zeros((len(PLlist[0])))
        for day in range(1, daynr_month(year, month)+1):
            """import data"""
            d= dataASR(year, month, sday=day, eday=day)
            d.impvar('SLP', level='surf')
            d.impvar('SST', level='surf')
            d.SST[d.SST<272]= np.nan #this sets the area with ice to land
            d.SST= ma.array(d.SST, mask= isnan(d.SST)) 
            
            d.impvar('T', level=500) 
            
            theta500= d.T500*(1000/500)**RCp
            thetaSST= d.SST*(1000/d.SLP)**RCp

            index= np.where(np.logical_and(PLlist[1] > (day-1)*8, PLlist[1]< day*8))[0] #indexes where PLlist[1] is on day
            for i in index:
                minval= np.min((PLlist[3][i] - d.lat)**2 + (PLlist[2][i] - d.lon)**2)
                pts2D= np.where((PLlist[3][i] - d.lat)**2 + (PLlist[2][i] - d.lon)**2 == minval)  
                
                t = int(PLlist[1][i]%8)

                dist= 2*55/15
                mask= radMask2D(d.SST.shape[1:], (int(pts2D[0]), int(pts2D[1])), radiusx= dist, radiusy= dist) #lat lon mask
                maskc= [t, mask] #tim, lat, lon mask (mask complete)
                stablist[i]= np.nanmean(thetaSST[maskc]- theta500[maskc])

        PLlist= np.vstack((PLlist, stablist))
        PLlist = PLlist[:, stablist > thetastablim]

        snr += len(remove_dublicate(PLlist[0]))
        print('after stab exclusion', snr)


        """u500north exclusion"""
        U500list= np.zeros((len(PLlist[0])))
        
        if len(PLlist[0]) > 0:
            for day in range(1, daynr_month(year, month)+1):
                """import data"""
                d= dataASR(year, month, sday=day, eday=day)
                d.impvar('U', level=500)
                d.impvar('V', level=500)  
                
                index= np.where(np.logical_and(PLlist[1] > (day-1)*8, PLlist[1]< day*8))[0] #indexes where PLlist[1] is on day
                for i in index:
                    masknorth= np.logical_and.reduce((d.lat >= PLlist[3][i], d.lon > PLlist[2][i] -1, d.lon < PLlist[2][i] +1))
                    t = int(PLlist[1][i]%8)
    
                    if len(np.where(masknorth == True)[0])== 0: #if there is no point north of the PL
                        U500list[i] = 0
                    else:
                        U500list[i] = np.max(np.sqrt(d.u500[t][masknorth]**2 + d.v500[t][masknorth]**2))
                    

        PLlist= np.vstack((PLlist, U500list))
        PLlist = PLlist[:, U500list < polarfrontlim]

        pnr += len(remove_dublicate(PLlist[0]))
        print('after polar front exclusion', pnr)




        """make the contribution of each PL to the PLlist"""
        if len(PLlist > 0):
            for i, PLnr in enumerate(remove_dublicate(PLlist[0])): #runs through all PLs and add also times before and after being PLpoint
                intermlist= T.PLlist[:, T.PLlist[0]== PLnr]
                intermlist=np.append(intermlist, np.zeros((4,intermlist.shape[1])), axis= 0) #append the column that indicates PL point or not
                intermlist=np.append(np.ones((1, intermlist.shape[1]))*month, intermlist, axis= 0)
                intermlist=np.append(np.ones((1, intermlist.shape[1]))*year, intermlist, axis= 0)
                
                #makes the column of PL point or not
                timones= PLlist[1, PLlist[0] ==PLnr] #times where the cyclone is a PL
                for timone in timones: #include other data for these times
                    intermlist[-1, intermlist[3]==timone] = 1
                    intermlist[-4:-1, intermlist[3]==timone]= PLlist[-3:, np.logical_and(PLlist[0]==PLnr, PLlist[1]== timone)]
                
                cPLlist = np.append(cPLlist, intermlist, axis= 1)
                
        
    WritePLlistComp2S(cPLlist, year, name='4', model='ASR')


print('cyclones:',cnr, 'vortexcl:', vnr, 'landexcl:', lnr, 'polarfrontexcl:', pnr, 'stabexcl:', snr)


print("--- %s seconds ---" % (time.time() - start_time))