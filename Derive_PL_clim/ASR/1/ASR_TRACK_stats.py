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
stat_track_del_theta_e500, stat_track_Water, stat_track_maxtheta_e500= [],[],[]
#stat_track_theta_e700, stat_track_theta_e500= [], []
#stat_track_p_PVd0, stat_track_p_PVd1, stat_track_p_PVd2 = [], [], []
stat_track_t_500, stat_track_t_700, stat_track_SST = [], [], []
stat_track_U500north = []
#stat_track_U500northint = []
#stat_track_U500northintabslim = []
#stat_track_U500northintabs = []
#stat_track_U500northmsqr = []
#stat_track_U500northintn = []

stat_track_TPL= []

for year in np.arange(2003, 2004):
    for month in [1, 2, 3, 4, 10, 11, 12]:
        
        T=TRACK_list_month2(year, month, duration_limit= 1, timestep='threehours', model='ASR')
#        T.local(minlon= -180, maxlon= 180, minlat= 30.5, maxlat= 84.5) #has to be within ERA data
        print('number of cyclones: ', len(remove_dublicate(T.PLlist[0])))

                
        print("--- %s seconds ---" % (time.time() - start_time))

        for day in range(1, daynr_month(year, month)+1):
            
            """PLlist for that day"""
            dayPLlist= T.PLlist[:, np.logical_and(T.PLlist[1] > (day-1)*8, T.PLlist[1] < day*8)]
            """import ASR data"""
            d= dataASR(year, month, day, level= 'all', reso= reso)
            d.SST[d.SST < 272]= 0 #replace with 0 for ice
            print("--- %s seconds ---" % (time.time() - start_time))
            
            U10= np.sqrt(d.u10**2+ d.v10**2)
            U500= np.sqrt(d.u500**2+ d.v500**2)
        
            
            """calculate properties for every cyclone """                          
            for i in range(len(dayPLlist[0])): #loop through every cyclone point
                minval= np.min((dayPLlist[3][i] - d.lat)**2 + (dayPLlist[2][i] - d.lon)**2)
                pts2D= np.where((dayPLlist[3][i] - d.lat)**2 + (dayPLlist[2][i] - d.lon)**2 == minval)
                t = int(dayPLlist[1][i]%8)
    #            pts3D= tuple((, int(pts2D[0]), int(pts2D[1])))
    
                latdist2, latdist1= 5*55/15, 2*55/15
               
#                londist2= int(latdist2//np.cos(np.deg2rad(dayPLlist[3][i])))
#                londist1= int(latdist1//np.cos(np.deg2rad(dayPLlist[3][i])))
                mask= radMask2D(d.SST.shape[1:], (int(pts2D[0]), int(pts2D[1])), radiusx= latdist1, radiusy= latdist1) #lat lon mask
                maskc= [t, mask] #tim, lat, lon mask (mask complete)
    
                if np.mean(d.SST[maskc]) > 272*3/4:  #true if we are over ocean and d.SST[maskc] has a value
                    d.SST[d.SST== 0]= np.nan #replace fill value 0 with nans

                    stat_track_vort += [dayPLlist[-1][i]]
                    stat_track_TPL += [dayPLlist[0][i]] #[yearL[n]*1E7+ monthL[n]*1E5 +TPLnr[n]]
    
                
                    """calculate potential temperature and equivalent potential temperature"""
                    RCp= 2/7 #R/Cp
                    theta700= d.T700[maskc]*(1000/700)**RCp
                    theta500= d.T500[maskc]*(1000/500)**RCp
                    thetaSST= d.SST[maskc]*(1000/d.SLP[maskc])**RCp
    
                    stat_track_stab500 += [np.nanmean(d.SST[maskc]- d.T500[maskc])]
                    stat_track_stab700 += [np.nanmean(d.SST[maskc]- d.T700[maskc])]
                    stat_track_theta500 += [np.nanmean(d.SST[maskc]- theta500)]
                    stat_track_theta700 += [np.nanmean(d.SST[maskc]- theta700)]
    ##                stat_track_p_PVmd0 += [d.p_PV[pts[0][i], pts[1][i], pts[2][i]]]
    ##                stat_track_p_PVmd1 +=[np.max(d.p_PV[maskc])]

                    stat_track_SST += [np.nanmean(d.SST[maskc])]
                    stat_track_t_500 += [np.mean(d.T500[maskc])]
                    stat_track_t_700 += [np.mean(d.T700[maskc])]
    #                
                    mask2= radMask2D(d.SST.shape[1:], (int(pts2D[0]), int(pts2D[1])), radiusx= latdist2, radiusy= latdist2) #lat lon mask
                    maskc2= [t, mask2] #tim, lat, lon mask (mask complete)
    ##                
    ##                p_PVmd2 +=[np.max(d.p_PV[maskc2])]
    #           
                    stat_track_U10 += [np.max(U10[maskc2])]
                    stat_track_U500 += [np.max(U500[maskc])]
                    
                    masknorth= np.logical_and.reduce((d.lat >= dayPLlist[3][i], d.lon > dayPLlist[2][i] -reso, d.lon < dayPLlist[2][i] +reso))
                    if len(np.where(masknorth == True)[0])== 0:
                        stat_track_U500north += [0]
                    else:
                        stat_track_U500north += [np.max(U500[t][masknorth])]

##                U500northint+= [np.sum(d.U500[pts[0][i], :Tlatraw[i]+1, Tlonraw[i]])]
##                U500northintabs+= [np.sum(np.abs(d.U500[pts[0][i], :Tlatraw[i]+1, Tlonraw[i]]))]
##                U500northintabslim+= [np.sum(np.abs(d.U500[pts[0][i], :Tlatraw[i]+1, Tlonraw[i]][np.abs(d.U500[pts[0][i], :Tlatraw[i]+1, Tlonraw[i]])> lim]))]
##                U500northmsqr+= [np.mean(np.square(d.U500[pts[0][i], :Tlatraw[i]+1, Tlonraw[i]]))]
##                U500northintn+= [np.mean(heapq.nlargest(nint, np.abs(U500[pts[0][i], :Tlatraw[i]+1, Tlonraw[i]])))]
#                stat_track_Water += [np.mean(d.Water[maskc])]
#

print('number of ocean cyclones: ', len(remove_dublicate(stat_track_TPL)))

#
#print("--- %s seconds ---" % (time.time() - start_time))
#                 
#"""analysis of the output"""                 
#
##stdfactor= 1.5
#"""the criteria are derived form TRACKmatchSTARS_stats3"""
#
#
##print('nr of cyclones: ', len(max_track_vort))
#
#
##print('p_PV d0', np.mean(stat_track_p_PVd0), p_PVd0_crit, len(np.where(stat_track_p_PVd0 < p_PVd0_crit)[0]))
##print('p_PV d1', np.mean(stat_track_p_PVd1), p_PVd1_crit, len(np.where(stat_track_p_PVd1 < p_PVd1_crit)[0]))
##print('p_PV d2', np.mean(stat_track_p_PVd2), p_PVd2_crit, len(np.where(stat_track_p_PVd2 < p_PVd2_crit)[0]))
#
##print('SST- T500', np.mean(stat_track_stab500), stab500crit, len(np.where(stat_track_stab500 < stab500crit)[0]))
##print('SST- T700', np.mean(stat_track_stab700), stab700crit, len(np.where(stat_track_stab700 < stab700crit)[0]))
##
##print('theta 700', np.mean(stat_track_theta700), t700crit, len(np.where(stat_track_theta700 < t700crit)[0]))
##print('theta 500', np.mean(stat_track_theta500), t500crit, len(np.where(stat_track_theta500 < t500crit)[0]))
##print('theta eqi 700', np.mean(stat_track_theta_e700), t_e700crit, len(np.where(stat_track_theta_e700 < t_e700crit)[0]))
##print('theta eqi 500', np.mean(stat_track_theta_e500), t_e500crit, len(np.where(stat_track_theta_e500 < t_e500crit)[0]))
###
##print('U10', np.mean(stat_track_U10), U10crit, len(np.where(stat_track_U10 < U10crit)[0]))
##print('Vort', np.mean(stat_track_vort), vortcrit, len(np.where(stat_track_vort < vortcrit)[0]))
###
###
##print('U500', np.mean(stat_track_U500), U500crit, len(np.where(stat_track_U500 > U500crit)[0]))
##print('U500 north', np.mean(stat_track_U500north), U500northcrit, len(np.where(stat_track_U500north < U500northcrit)[0]))
###print('U500 north int', np.mean(stat_track_U500northint), U500northintcrit, len(np.where(stat_track_U500northint < U500northintcrit)[0]))
###print('U500 north int abs', np.mean(stat_track_U500northintabs), U500northintabscrit, len(np.where(stat_track_U500northintabs < U500northintabscrit)[0]))
###print('U500 north int abs lim', np.mean(stat_track_U500northintabslim), U500northintabslimcrit, len(np.where(stat_track_U500northintabslim < U500northintabslimcrit)[0]))
###print('U500 north mean square', np.mean(stat_track_U500northmsqr), U500northmsqrcrit, len(np.where(stat_track_U500northmsqr < U500northmsqrcrit)[0]))
###print('U500 north', np.mean(stat_track_U500northintn), U500northintncrit, len(np.where(stat_track_U500northintn < U500northintncrit)[0]))
##
##print('Water', np.mean(stat_track_Water), Watercrit, len(np.where(stat_track_Water < Watercrit)[0]))
#
#
#print("--- %s seconds ---" % (time.time() - start_time))