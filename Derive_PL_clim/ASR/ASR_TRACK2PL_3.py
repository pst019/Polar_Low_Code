#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 15:06:41 2017

@author: pst019
"""
import os
user = os.getcwd().split('/')[2]

from netCDF4 import Dataset
import numpy as np
import sys  #to import the functions from a different directory
import pickle
import scipy.ndimage.filters as filters

import csv

import datetime
import time
start_time = time.time()



sys.path.insert(0, '/home/'+user+'/polar_low_code/Functions/')
from f_imp_ERA2 import *
from f_imp_ASR import *

from f_impLists import *
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)






fignr= 10

"""version 6"""
#versionname= '6'
#Tdurationlim=1 #TRACK duration limit in 3hourly time steps - higher less PLs
#vortlim= None #the higher the less PLs
#thetastablim= np.float64(-9.054) #the higher the less PLs
#diststab= 4*55/15
#
#windlim= None #np.float64(17.5) #the higher the less PLs
#polarfrontlim= np.float64(31.46) #the higher the more PLs
#pdifflim= np.float(2.055) 
#distpdiff= 6*55/15

"""version 7"""
versionname= '7'
Tdurationlim=1 #TRACK duration limit in 3hourly time steps - higher less PLs
vortlim= None #the higher the less PLs
thetastablim= np.float64(-9.054) #the higher the less PLs
diststab= 2*55/15

windlim= None #np.float64(17.5) #the higher the less PLs
polarfrontlim= np.float64(29.58) #the higher the more PLs
pdifflim= np.float(2.382) 
distpdiff= 6*55/15

RCp= 2/7 #R/Cp

cnr, vnr, lnr, pnr, jnr, snr, wnr, pdnr= 0,0,0,0,0,0,0,0

for year in np.arange(2000, 2013):
    cPLlist= np.zeros((12,0)) #complete PLlist (year, month, TrackPLnr, tPL, lon, lat, vorticity, land, 10mwind, stability, U500north, PLpoint(yes= 1, no= 0))

    
    for month in np.arange(1, 13):
        print(year, month)

        shearlist, stablist, U500northlist= [], [], []

        T=TRACK_list_month2(year, month, duration_limit= Tdurationlim, timestep='threehours', model='ASR')

        cnr += len(remove_dublicate(T.PLlist[0]))
        print('before exclusion', cnr)

        if vortlim != None:
            """vorticity exclusion"""
            PLlist= T.PLlist[:, T.PLlist[-1] > vortlim]
            vnr += len(remove_dublicate(PLlist[0]))
            print('after vorticity exclusion', vnr)
        else: PLlist= T.PLlist



        """land and ice exclusion"""
        landlist= np.zeros((len(PLlist[0])))
        
        for day in range(1, daynr_month(year, month)+1):
            """import data"""
            d= dataASR(year, month, sday=day, eday=day)
            d.impvar('SST', level='surf')
            d.SST[d.SST<272]= np.nan #this sets the area with ice to land
            d.SST= ma.array(d.SST, mask= isnan(d.SST))              
            
            index= np.where(np.logical_and(PLlist[1] >= (day-1)*8, PLlist[1]< day*8))[0] #indexes where PLlist[1] is on day
            for i in index:
                minval= np.min((PLlist[3][i] - d.lat)**2 + (PLlist[2][i] - d.lon)**2)
                pts2D= np.where((PLlist[3][i] - d.lat)**2 + (PLlist[2][i] - d.lon)**2 == minval)  
                
                t = int(PLlist[1][i]%8)

                dist= 4*55/15
                mask2= radMask2D(d.SST.shape[1:], (int(pts2D[0]), int(pts2D[1])), radiusx= dist, radiusy= dist) #lat lon mask
                maskc2= [t, mask2] #tim, lat, lon mask (mask complete)
                landlist[i]= len(np.where(d.SST[maskc2].mask == True)[0]) /len(d.SST[maskc2])

        PLlist= np.vstack((PLlist, landlist)) #this is maybe not so interesting in the final list
        PLlist = PLlist[:, landlist < .25]

        lnr += len(remove_dublicate(PLlist[0]))
        print('after land exclusion', lnr)
        

        """pdiff exclusion"""
        pdifflist= np.zeros((len(PLlist[0])))
        for day in range(1, daynr_month(year, month)+1):
            """import data"""
            d= dataASR(year, month, sday=day, eday=day)
            d.impvar('SLP', level='surf')

            index= np.where(np.logical_and(PLlist[1] >= (day-1)*8, PLlist[1]< day*8))[0] #indexes where PLlist[1] is on day
            for i in index:
                minval= np.min((PLlist[3][i] - d.lat)**2 + (PLlist[2][i] - d.lon)**2)
                pts2D= np.where((PLlist[3][i] - d.lat)**2 + (PLlist[2][i] - d.lon)**2 == minval)                
                t = int(PLlist[1][i]%8)

                dist= distpdiff
                mask= radMask2D(d.SLP.shape[1:], (int(pts2D[0]), int(pts2D[1])), radiusx= dist, radiusy= dist) #lat lon mask
                maskc= [t, mask] #tim, lat, lon mask (mask complete)
                pdifflist[i]= np.mean(d.SLP[maskc]) - d.SLP[t, pts2D[0], pts2D[1]]
                
        PLlist= np.vstack((PLlist, pdifflist))
        PLlist = PLlist[:, pdifflist > pdifflim]

        pdnr += len(remove_dublicate(PLlist[0]))
        print('after pdiff exclusion', pdnr)
            
        


        if windlim != None:
            """wind exclusion"""
            Ulist= np.zeros((len(PLlist[0])))
            for day in range(1, daynr_month(year, month)+1):
                """import data"""
                d= dataASR(year, month, sday=day, eday=day)
                d.impvar('U', level='surf')
                d.impvar('V', level='surf')
                U10= np.sqrt(d.u10**2 + d.v10**2)
                
                index= np.where(np.logical_and(PLlist[1] >= (day-1)*8, PLlist[1]< day*8))[0] #indexes where PLlist[1] is on day
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

            index= np.where(np.logical_and(PLlist[1] >= (day-1)*8, PLlist[1]< day*8))[0] #indexes where PLlist[1] is on day
            for i in index:
                minval= np.min((PLlist[3][i] - d.lat)**2 + (PLlist[2][i] - d.lon)**2)
                pts2D= np.where((PLlist[3][i] - d.lat)**2 + (PLlist[2][i] - d.lon)**2 == minval)  
                
                t = int(PLlist[1][i]%8)

                dist= diststab
                mask= radMask2D(d.SST.shape[1:], (int(pts2D[0]), int(pts2D[1])), radiusx= dist, radiusy= dist) #lat lon mask
                maskc= [t, mask] #tim, lat, lon mask (mask complete)
                stablist[i]= np.nanmax(thetaSST[maskc]- theta500[maskc])

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
                
                index= np.where(np.logical_and(PLlist[1] >= (day-1)*8, PLlist[1]< day*8))[0] #indexes where PLlist[1] is on day
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
                intermlist=np.append(intermlist, np.zeros((5,intermlist.shape[1])), axis= 0) #append the column that indicates PL point or not
                intermlist=np.append(np.ones((1, intermlist.shape[1]))*month, intermlist, axis= 0)
                intermlist=np.append(np.ones((1, intermlist.shape[1]))*year, intermlist, axis= 0)
                
                #makes the column of PL point or not
                timones= PLlist[1, PLlist[0] ==PLnr] #times where the cyclone is a PL
                for timone in timones: #include other data for these times
                    intermlist[-1, intermlist[3]==timone] = 1
                    intermlist[-5:-1, intermlist[3]==timone]= PLlist[-4:, np.logical_and(PLlist[0]==PLnr, PLlist[1]== timone)]
                
                cPLlist = np.append(cPLlist, intermlist, axis= 1)
                
        
    WritePLlistComp2S(cPLlist, year, name=versionname, model='ASR', rowname= ['year', 'month', 'TPLnr', 'time pnt (3hourly)', 'lon', 'lat', 'vort_filt', 'landdistr', 'meanSLP - SLP' , 'thetaSST-theta500', 'U500north', 'PLpoint'])


print('cyclones:',cnr, 'vortexcl:', vnr, 'landexcl:', lnr, 'pdiffexcl:', pdnr, 'stabexcl:', snr, 'polarfrontexcl:', pnr)


print("--- %s seconds ---" % (time.time() - start_time))