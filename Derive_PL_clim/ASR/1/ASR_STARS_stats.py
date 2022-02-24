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
import datetime

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
tPLstart= (PLlist[:, 4]+ 1)//2
tPLend= PLlist[:, 5]//2



"""creation of the PL statistics"""
statvort = []
statU10 = []
statstab500 = []
statstab700, stattheta700, stattheta500 = [], [], []
statt_500, statt_700 = [], []
statSST= []
statU500= []
statU500north= []

# stattheta_e700, stattheta_e500 = [], []
#stattheta_e700_2, stattheta_e500_2, stattheta_e500_3 ,stattheta_e500_4 = [], [], [], []
#statWater, statmaxtheta_e500= [], []
#statp_PVd0 = []
#statp_PVd1, statp_PVd2 = [], []
#statuPVnorthint = []
#statuPVnorthintabslim = []
#statuPVnorthintabs = []
#statuPVnorthmsqr = []
#statUPVnorthintn= []

statSPL = []
statTPL = []
nint= 10

for n in range(len(TPLnr)): #loop through all TRACKs with a connected STARS PLs  cyclones, n is the current row of the list
    #import data of the TRACK PL
    T= TRACK_list_month2(yearL[n], monthL[n], duration_limit= 1, timestep='threehours', model='ASR')
    
    #the index is the row in the TRACK list that are of interest - including the time information
    index= np.where(np.logical_and(T.PLlist[0]==TPLnr[n], T.PLlist[1] >= tPLstart[n], T.PLlist[1] <= tPLend[n]))[0]
    
    #the start and end time of the PL
    tim0= int(T.PLlist[1][index][0])
    day0, h3_0= tim0//8+1 , tim0%8
    tim1= int(T.PLlist[1][index][-1])
    day1, h3_1= tim1//8+1 , tim1%8
    
    #import ASR data
    d= dataASR(yearL[n], monthL[n], sday=day0, eday=day1, sh3= h3_0, eh3= h3_1 +1, level= 'all')
    d.SST[d.SST < 272]= 0 #replace with 0 for ice

#    U_500= np.sqrt(d.u500**2+ d.v500**2)
    
    #the lat lon coordinatinates of the PL in ASR grid points
#    minval= np.min((T.PLlist[3][index] - d.lat)**2 + (T.PLlist[2][index] - d.lon)**2)
#    Tlatraw= np.array(np.round(-(T.PLlist[3][index]-d.lat[0])*2), dtype= int)
#    Tlonraw= np.array(np.round((T.PLlist[2][index] - d.lon[0])*2), dtype= int)
#    pts2D=(list(Tlatraw), list(Tlonraw))  #these are the cyclone points in ERA data coordinates
   
         
    for i, indexi in enumerate(index): #loop through all time steps of the given TPLnr

        #calculate the lat lon coordinatinates of the PL in ASR grid points
        minval= np.min((T.PLlist[3][indexi] - d.lat)**2 + (T.PLlist[2][indexi] - d.lon)**2)
        pts2D= np.where((T.PLlist[3][indexi] - d.lat)**2 + (T.PLlist[2][indexi] - d.lon)**2 == minval)
        
        latdist2, latdist1= 5*55/15, 2*55/15
       
#        londist2= int(latdist2//np.cos(np.deg2rad(T.PLlist[3][indexi])))
#        londist1= int(latdist1//np.cos(np.deg2rad(T.PLlist[3][indexi])))
#
#
        mask= radMask2D(d.SST.shape[1:], (pts2D[0], pts2D[1]), radiusx= latdist1, radiusy= latdist1) #lat lon mask
        maskc= [i, mask] #tim, lat, lon mask (mask complete)

        if np.mean(d.SST[maskc]) > 200:  #true if we are over ocean and d.SST[maskc] has a value
        
            d.SST[d.SST== 0]= np.nan #replace fill value 0 with nans
            statvort += [T.PLlist[-1][indexi]]
            statSPL += [SPLnr[n]]
            statTPL += [yearL[n]*1E7+ monthL[n]*1E5 +TPLnr[n]]
            
            """calculate potential temperature and equivalent potential temperature"""
            RCp= 2/7 #R/Cp
            theta700= d.T700[maskc]*(1000/700)**RCp
            theta500= d.T500[maskc]*(1000/500)**RCp
            thetaSST= d.SST[maskc]*(1000/d.SLP[maskc])**RCp
#    
#            Lc, Cp= 2501, 1                
#            theta_eSST= thetaSST* np.exp(d.sHum[:,0][maskc]* Lc/(d.SST[maskc]* Cp))
#            theta_e1000= d.T[:,0][maskc]* np.exp(d.sHum[:,0][maskc]* Lc/(d.T[:,0][maskc]* Cp))            
#            theta_e700= theta700* np.exp(d.sHum[:,1][maskc]* Lc/(d.T[:,1][maskc]* Cp))
#            theta_e500= theta500* np.exp(d.sHum[:,2][maskc]* Lc/(d.T[:,2][maskc]* Cp))
#    #        
            statstab500 += [np.nanmean(d.SST[maskc]- d.T500[maskc])]
            statstab700 += [np.nanmean(d.SST[maskc]- d.T700[maskc])]
            stattheta500 += [np.nanmean(thetaSST- theta500)]
            stattheta700 += [np.nanmean(thetaSST- theta700)]
#            stattheta_e500 += [np.mean(theta_eSST- theta_e500)]
#            stattheta_e700 += [np.mean(theta_eSST- theta_e700)]
#            stattheta_e500_2 += [np.mean(thetaSST- theta_e500)]
#            stattheta_e500_3 += [np.mean(d.SST[maskc]- theta_e500)]
#            stattheta_e500_4 += [np.mean(theta_e1000- theta_e500)]
#
##            stattheta_e700_2 += [np.mean(thetaSST- theta_e700)]            
#            
#            statmaxtheta_e500 += [np.max(theta_eSST- theta_e500)]
#    #        p_PVmd0 += [d.p_PV[i, pts2D[0][i], pts2D[1][i]]]
#    #        p_PVmd1 +=[np.max(d.p_PV[maskc])]
#    #
            statt_500 += [np.mean(d.T500[maskc])]
            statt_700 += [np.mean(d.T700[maskc])]
            statSST += [np.nanmean(d.SST[maskc])]
#            
            mask2= radMask2D(d.SST.shape[1:], (pts2D[0], pts2D[1]), radiusx= latdist2, radiusy= latdist2) #lat lon mask
            maskc2= [i, mask2] #tim, lat, lon mask (mask complete)
#            
            statU10 += [np.max(np.sqrt(d.u10[maskc2]**2+ d.v10[maskc2]**2))]
            statU500 += [np.max(np.sqrt(d.u500[maskc]**2 + d.v500[maskc]**2))]     
            
            masknorth= np.logical_and.reduce((d.lat > T.PLlist[3][indexi], d.lon > T.PLlist[2][indexi] -1, d.lon < T.PLlist[2][indexi] +1))
            statU500north += [np.max(np.sqrt(d.u500[i][masknorth]**2 + d.v500[i][masknorth]**2))]

#    ##        uPVnorthint+= [np.sum(d.uPV[i, :Tlatraw[i]+1, Tlonraw[i]])]
#    ##        uPVnorthintabs+= [np.sum(np.abs(d.uPV[i, :Tlatraw[i]+1, Tlonraw[i]]))]
#    ##        uPVnorthintabslim+= [np.sum(np.abs(d.uPV[i, :Tlatraw[i]+1, Tlonraw[i]][np.abs(d.uPV[i, :Tlatraw[i]+1, Tlonraw[i]])> lim]))]
#    ##        uPVnorthmsqr+= [np.mean(np.square(d.uPV[i, :Tlatraw[i]+1, Tlonraw[i]]))]
#    ##        UPVnorthintn+= [np.mean(heapq.nlargest(nint, np.abs(UPV[i, :Tlatraw[i]+1, Tlonraw[i]])))]
#    #
#            statWater += [np.mean(d.Water[maskc])]
#    ##        p_PVmd2 +=[np.max(d.p_PV[maskc2])]
##        
##        else:
##            print(i, d.SST[maskc]) #in this case it is masked
"""analysis of the output - mean of all values during a polar low"""                 
