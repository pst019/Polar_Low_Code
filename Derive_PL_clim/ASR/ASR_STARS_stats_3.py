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
import time
start_time = time.time()


import heapq

""" import TRACKmatchSTARS
 all the matches, tPLstart gives the starting point in time of the match"""

filename= Mediadir+'PL/PLclim/MatchList/ASR_TRACKmatchSTARS.csv' #_Sdurlim3.csv'

import pandas as pd
PLlist= np.array(pd.read_table(filename, sep="\t"))


"""combine double PLs in TRACKmatchSTARS"""
doubleTRACKS= remove_dublicate([x for x in PLlist[:, 1] if list(PLlist[:, 1]).count(x) > 1])

print('double TRACKs', doubleTRACKS)

removeind= []
for n in range(len(doubleTRACKS)): #remove the double tracks
    dind=  np.where(PLlist[:, 1] == doubleTRACKS[n])[0]
    tPLstartfirst= PLlist[:, 4][dind[0]]
    tPLstartsecond= PLlist[:, 4][dind[-1]]
    tPLendfirst= PLlist[:, 5][dind[0]]
    if tPLstartsecond >= tPLstartfirst  and tPLstartsecond <=  tPLendfirst:
        print('merge, SPLs:', PLlist[:,0][dind[0]], PLlist[:,0][dind[-1]], 'at times:', 'startPL1', tPLstartfirst, 'StartPL2', tPLstartsecond, 'endPL1', tPLendfirst)
        
        PLlist[dind[0], 5] = PLlist[dind[-1], 5]  #put the end of the second polar low for the end of the first pl
        removeind += [dind[-1]]

PLlist= np.delete(PLlist, removeind, 0) #remove the second PL

print('nr of PLs:', len(PLlist[:, 0]))
                 
SPLnr= PLlist[:, 0]
TPLnr= PLlist[:, 1]
yearL= PLlist[:, 2]
monthL= PLlist[:, 3]
tPLstart= PLlist[:, 4] #+ 1)//2
tPLend= PLlist[:, 5] #//2



"""creation of the PL statistics"""
statvort = []
statU10 = []
statstab500 = []
statstab700, stattheta700, stattheta500 = [], [], []
statt_500, statt_700 = [], []
statSST= []
statU500= []
statU500north= []

stattheta500_2, stattheta500_3= [], []
stattheta500_max, stattheta500_2_max, stattheta500_3_max= [], [], []
statstab500_max, statstab700_max, stattheta700_max, statt_500_max, statt_700_max, statSST_max =  [], [], [], [], [], []

statU10_2, statU10_3 = [], []
statpdiff1, statpdiff2, statpdiff3 =  [], [], []

statSPL = []
statTPL = []
nint= 10

for n in range(len(TPLnr)): #loop through all TRACKs with a connected STARS PLs  cyclones, n is the current row of the list
    #import data of the TRACK PL
    T= TRACK_list_month2(yearL[n], monthL[n], duration_limit= 1, timestep='threehours', model='ASR')
    
    #the index is the row in the TRACK list that are of interest - including the time information
    index= np.where(np.logical_and.reduce((T.PLlist[0]==TPLnr[n], T.PLlist[1] >= tPLstart[n], T.PLlist[1] <= tPLend[n])))[0]
    
    #the start and end time of the PL
    tim0= int(T.PLlist[1][index][0])
    day0, h3_0= tim0//8+1 , tim0%8
    tim1= int(T.PLlist[1][index][-1])
    day1, h3_1= tim1//8+1 , tim1%8

    print(n, SPLnr[n], "--- %s seconds ---" % (time.time() - start_time))

    #import ASR data
    d= dataASR(yearL[n], monthL[n], sday=day0, eday=day1, sh3= h3_0, eh3= h3_1 +1, level= 'all')
    d.SST[d.SST<272]= np.nan #this sets the area with ice to land
    d.SST= ma.array(d.SST, mask= isnan(d.SST))
    
         
    for i, indexi in enumerate(index): #loop through all time steps of the given TPLnr

        #calculate the lat lon coordinatinates of the PL in ASR grid points
        minval= np.min((T.PLlist[3][indexi] - d.lat)**2 + (T.PLlist[2][indexi] - d.lon)**2)
        pts2D= np.where((T.PLlist[3][indexi] - d.lat)**2 + (T.PLlist[2][indexi] - d.lon)**2 == minval)
        
        dist3, dist2, dist1= 6*55/15, 4*55/15, 2*55/15

        mask= radMask2D(d.SST.shape[1:], (pts2D[0], pts2D[1]), radiusx= dist1, radiusy= dist1) #lat lon mask
        maskc= [i, mask] #tim, lat, lon mask (mask complete)
        mask2= radMask2D(d.SST.shape[1:], (pts2D[0], pts2D[1]), radiusx= dist2, radiusy= dist2) #lat lon mask
        maskc2= [i, mask2] #tim, lat, lon mask (mask complete)
        mask3= radMask2D(d.SST.shape[1:], (pts2D[0], pts2D[1]), radiusx= dist3, radiusy= dist3) #lat lon mask
        maskc3= [i, mask3] #tim, lat, lon mask (mask complete)

        landistr= len(np.where(d.SST[maskc2].mask == True)[0]) /len(d.SST[maskc2])
#        print(landistr)
        if landistr < .25:  #true if we are over ocean and d.SST[maskc] has a value
        
            statvort += [T.PLlist[-1][indexi]]
            statSPL += [SPLnr[n]]
            statTPL += [yearL[n]*1E7+ monthL[n]*1E5 +TPLnr[n]]
            
            """calculate potential temperature and equivalent potential temperature"""
            RCp= 2/7 #R/Cp
            theta700= d.T700[maskc]*(1000/700)**RCp
            theta500= d.T500[maskc]*(1000/500)**RCp
            thetaSST= d.SST[maskc]*(1000/d.SLP[maskc])**RCp
     
            statstab500 += [np.nanmean(d.SST[maskc]- d.T500[maskc])]
            statstab700 += [np.nanmean(d.SST[maskc]- d.T700[maskc])]
            stattheta500 += [np.nanmean(thetaSST- theta500)]
            stattheta700 += [np.nanmean(thetaSST- theta700)]

            statSST += [np.nanmean(d.SST[maskc])]
            statt_500 += [np.mean(d.T500[maskc])]
            statt_700 += [np.mean(d.T700[maskc])]

            statstab500_max += [np.nanmax(d.SST[maskc]- d.T500[maskc])]
            statstab700_max += [np.nanmax(d.SST[maskc]- d.T700[maskc])]
            stattheta700_max += [np.nanmax(thetaSST- theta700)]
            stattheta500_max += [np.nanmax(thetaSST- theta500)]
            statSST_max += [np.nanmax(d.SST[maskc])]
            statt_500_max += [np.max(d.T500[maskc])]
            statt_700_max += [np.max(d.T700[maskc])]
                    
            theta500_2= d.T500[maskc2]*(1000/500)**RCp
            thetaSST_2= d.SST[maskc2]*(1000/d.SLP[maskc2])**RCp
            theta500_3= d.T500[maskc3]*(1000/500)**RCp
            thetaSST_3= d.SST[maskc3]*(1000/d.SLP[maskc3])**RCp                           
            stattheta500_2 += [np.nanmean(thetaSST_2- theta500_2)]
            stattheta500_3 += [np.nanmean(thetaSST_3- theta500_3)]
            stattheta500_2_max += [np.nanmax(thetaSST_2- theta500_2)]
            stattheta500_3_max += [np.nanmax(thetaSST_3- theta500_3)]


            statU10 += [np.max(np.sqrt(d.u10[maskc]**2+ d.v10[maskc]**2))]
            statU10_2 += [np.max(np.sqrt(d.u10[maskc2]**2+ d.v10[maskc2]**2))]                                   
            statU10_3 += [np.max(np.sqrt(d.u10[maskc3]**2+ d.v10[maskc3]**2))]
            statU500 += [np.max(np.sqrt(d.u500[maskc]**2 + d.v500[maskc]**2))]     
            
            masknorth= np.logical_and.reduce((d.lat > T.PLlist[3][indexi], d.lon > T.PLlist[2][indexi] -1, d.lon < T.PLlist[2][indexi] +1))
            statU500north += [np.max(np.sqrt(d.u500[i][masknorth]**2 + d.v500[i][masknorth]**2))]

            statpdiff1 += [np.float(d.SLP[i, pts2D[0], pts2D[1] ] - np.mean(d.SLP[maskc])) *-10]
            statpdiff2 += [np.float(d.SLP[i, pts2D[0], pts2D[1] ] - np.mean(d.SLP[maskc2])) *-10]
            statpdiff3 += [np.float(d.SLP[i, pts2D[0], pts2D[1] ] - np.mean(d.SLP[maskc3])) *-10]            
