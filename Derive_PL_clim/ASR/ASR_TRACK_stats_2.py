#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 16:38:42 2017

@author: pst019
"""
import numpy as np

import sys

sys.path.insert(0, '/home/'+user+'/codeDerive_PL_clim/ERA/')
from f_imp_ERA2 import *
from f_imp_ASR import *

from f_impLists import *

sys.path.insert(0, '/home/'+user+'/codeAROME/')
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
import heapq
#import datetime

import time
start_time = time.time()

reso= 1

"""creation of cyclone statistics"""
stat_track_stab500, stat_track_stab700, stat_track_U10, stat_track_vort, stat_track_theta700, stat_track_theta500, stat_track_U500=  [], [], [], [], [], [], []
stat_track_theta500_max, stat_track_theta500_2_max, stat_track_theta500_3_max = [], [], []
stat_track_stab500_max, stat_track_stab700_max, stat_track_theta700_max, stat_track_t_500_max, stat_track_t_700_max, stat_track_SST_max =  [], [], [], [], [], []

stat_track_t_500, stat_track_t_700, stat_track_SST = [], [], []
stat_track_U500north = []

stat_track_theta500_2, stat_track_theta500_3= [], []
stat_track_U10_2, stat_track_U10_3 = [], []

stat_track_pdiff1, stat_track_pdiff2, stat_track_pdiff3= [], [], []

stat_track_TPL= []

for year in np.arange(2003, 2004):
    for month in [1, 2, 3, 4, 10, 11, 12]:
        
        T=TRACK_list_month2(year, month, duration_limit= 1, timestep='threehours', model='ASR')
#        T.local(minlon= -180, maxlon= 180, minlat= 30.5, maxlat= 84.5) #has to be within ERA data
        print('number of cyclones: ', len(remove_dublicate(T.PLlist[0])))

        print('month', month)
                
        print("--- %s seconds ---" % (time.time() - start_time))

        for day in range(1, daynr_month(year, month)+1):
            
            print(day)
            """PLlist for that day"""
            dayPLlist= T.PLlist[:, np.logical_and(T.PLlist[1] >= (day-1)*8, T.PLlist[1] < day*8)]
            """import ASR data"""
            d= dataASR(year, month, day, level= 'all', reso= reso)
            d.SST[d.SST < 272]= np.nan #replace with 0 for ice
            d.SST= ma.array(d.SST, mask= isnan(d.SST))

            print("--- %s seconds ---" % (time.time() - start_time))
            
            U10= np.sqrt(d.u10**2+ d.v10**2)
            U500= np.sqrt(d.u500**2+ d.v500**2)
            print("--- %s seconds ---" % (time.time() - start_time))
   
            
            """calculate properties for every cyclone """                          
            for i in range(len(dayPLlist[0])): #loop through every cyclone point
                minval= np.min((dayPLlist[3][i] - d.lat)**2 + (dayPLlist[2][i] - d.lon)**2)
                pts2D= np.where((dayPLlist[3][i] - d.lat)**2 + (dayPLlist[2][i] - d.lon)**2 == minval)
                t = int(dayPLlist[1][i]%8)
    #            pts3D= tuple((, int(pts2D[0]), int(pts2D[1])))
    
                dist3, dist2, dist1= 6*55/15, 4*55/15, 2*55/15

                mask= radMask2D(d.SST.shape[1:], (int(pts2D[0]), int(pts2D[1])), radiusx= dist1, radiusy= dist1) #lat lon mask
                maskc= [t, mask] #tim, lat, lon mask (mask complete)
                mask2= radMask2D(d.SST.shape[1:], (int(pts2D[0]), int(pts2D[1])), radiusx= dist2, radiusy= dist2) #lat lon mask
                maskc2= [t, mask2] #tim, lat, lon mask (mask complete)
                mask3= radMask2D(d.SST.shape[1:], (int(pts2D[0]), int(pts2D[1])), radiusx= dist3, radiusy= dist3) #lat lon mask
                maskc3= [t, mask3] #tim, lat, lon mask (mask complete)


                landistr= len(np.where(d.SST[maskc2].mask == True)[0]) /len(d.SST[maskc2])

                if landistr < .25:  #true if we are over ocean and d.SST[maskc] has a value

                    stat_track_vort += [dayPLlist[-1][i]]
                    stat_track_TPL += [year*1E7+ month*1E5 + dayPLlist[0][i]] #[yearL[n]*1E7+ monthL[n]*1E5 +TPLnr[n]]
    
                
                    """calculate potential temperature and equivalent potential temperature"""
                    RCp= 2/7 #R/Cp
                    theta700= d.T700[maskc]*(1000/700)**RCp
                    theta500= d.T500[maskc]*(1000/500)**RCp
                    thetaSST= d.SST[maskc]*(1000/d.SLP[maskc])**RCp
    
                    stat_track_stab500 += [np.nanmean(d.SST[maskc]- d.T500[maskc])]
                    stat_track_stab700 += [np.nanmean(d.SST[maskc]- d.T700[maskc])]
                    stat_track_theta500 += [np.nanmean(thetaSST- theta500)]
                    stat_track_theta700 += [np.nanmean(thetaSST- theta700)]

                    stat_track_SST += [np.nanmean(d.SST[maskc])]
                    stat_track_t_500 += [np.mean(d.T500[maskc])]
                    stat_track_t_700 += [np.mean(d.T700[maskc])]

                    stat_track_stab500_max += [np.nanmax(d.SST[maskc]- d.T500[maskc])]
                    stat_track_stab700_max += [np.nanmax(d.SST[maskc]- d.T700[maskc])]
                    stat_track_theta700_max += [np.nanmax(thetaSST- theta700)]
                    stat_track_theta500_max += [np.nanmax(thetaSST- theta500)]
                    stat_track_SST_max += [np.nanmax(d.SST[maskc])]
                    stat_track_t_500_max += [np.max(d.T500[maskc])]
                    stat_track_t_700_max += [np.max(d.T700[maskc])]

                    theta500_2= d.T500[maskc2]*(1000/500)**RCp
                    thetaSST_2= d.SST[maskc2]*(1000/d.SLP[maskc2])**RCp
                    theta500_3= d.T500[maskc3]*(1000/500)**RCp
                    thetaSST_3= d.SST[maskc3]*(1000/d.SLP[maskc3])**RCp                           
                    stat_track_theta500_2 += [np.nanmean(thetaSST_2- theta500_2)]
                    stat_track_theta500_3 += [np.nanmean(thetaSST_3- theta500_3)]
                    stat_track_theta500_2_max += [np.nanmax(thetaSST_2- theta500_2)]
                    stat_track_theta500_3_max += [np.nanmax(thetaSST_3- theta500_3)]            

                    stat_track_U10 += [np.max(U10[maskc])]
                    stat_track_U10_2 += [np.max(U10[maskc2])]                                   
                    stat_track_U10_3 += [np.max(U10[maskc3])]
                    stat_track_U500 += [np.max(U500[maskc])]

                    stat_track_pdiff1 += [np.float(d.SLP[t, pts2D[0], pts2D[1] ] - np.mean(d.SLP[maskc])) *-10]
                    stat_track_pdiff2 += [np.float(d.SLP[t, pts2D[0], pts2D[1] ] - np.mean(d.SLP[maskc2]))* -10]
                    stat_track_pdiff3 += [np.float(d.SLP[t, pts2D[0], pts2D[1] ] - np.mean(d.SLP[maskc3]))*-10] 
                            
                    masknorth= np.logical_and.reduce((d.lat >= dayPLlist[3][i], d.lon > dayPLlist[2][i] -reso, d.lon < dayPLlist[2][i] +reso))
                    if len(np.where(masknorth == True)[0])== 0:
                        stat_track_U500north += [0]
                    else:
                        stat_track_U500north += [np.max(U500[t][masknorth])]



print('number of ocean cyclones: ', len(remove_dublicate(stat_track_TPL)))
