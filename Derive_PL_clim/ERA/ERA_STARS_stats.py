#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 16:38:42 2017

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/1692A00D929FEF8B/'

import numpy as np

import sys

sys.path.insert(0, '/home/'+user+'/polar_low_code/Functions/')
from f_imp_ERA2 import *
from f_impLists import *
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_meteo import * #meteorological functions
from f_useful import *
import datetime


""" import TRACKmatchSTARS
 all the matches, tPLstart gives the starting point in time of the match"""

filename= Mediadir+'PL/PLclim/MatchList/TRACKmatchSTARS.csv' #_Sdurlim3.csv'

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
stattheta500_2, stattheta500_3, stattheta500_max, stattheta500_max_2, stattheta500_max_3= [], [], [], [], []
statU10_2, statU10_3 = [], []

statstab700, stattheta700, stattheta500, stattheta_e700, stattheta_e500 = [], [], [], [], []
stattheta_e700_2, stattheta_e700_3, stattheta_e700_max, stattheta_e700_max_2, stattheta_e700_max_3= [], [], [], [], []

stattheta850_1, stattheta850_2, stattheta850_3, stattheta850_max_1, stattheta850_max_2= [], [], [], [], []
stattheta_e850_1, stattheta_e850_2, stattheta_e850_3, stattheta_e850_max_1, stattheta_e850_max_2= [], [], [], [], []

statWater, statmaxtheta_e500= [], []
statt_500, statt_700, statSST = [], [], []
statattheta_e700, statattheta_e850= [], [] #the equivalent potential temperature at 700 and 850 hPa
statPBH= []
statMCAO_Kolstad = []
statMCAO_Kolstad_0 = []
statMCAO_Kolstad_0_500 = []
statgradtheta_e850north= []

statUPV= []
statUPVnorth= []
statUPVnorthM2p5, statUPVnorthM5p0, statUPVnorthM7p5= [], [], []
statUPL500north = []

statp_PV_3, statp_PV_4, statp_PV_5= [], [], []  
statMSLP= []
statgradtheta_e1, statgradtheta_e2, statgradtheta_e3, statgradtheta_e5 = [], [], [], []
statpdiff1, statpdiff2, statpdiff3, statpdiff5 =  [], [], [], []

statthetaPV_1, statthetaPV_3, statthetaPV_5, statthetaPV_1min, statthetaPV_3min, statthetaPV_5min= [], [], [], [], [], []
statthetaSSTmPV_1, statthetaSSTmPV_3, statthetaSSTmPV_1_max, statthetaSSTmPV_3_max= [], [], [], []
statCAPE, statCAPE_2, statCAPE_3, statCAPE_max, statCAPE_max_2, statCAPE_max_3= [], [], [], [], [], []

statSPL = []
statTPL = []
nint= 10

for n in range(len(TPLnr)): #loop through all TRACKs with a connected STARS PLs  cyclones, n is the current row of the list
    #import data of the TRACK PL
    T= TRACK_list_month(yearL[n], monthL[n])
    
    #the index is the row in the TRACK list that are of interest - including the time information
    index= np.where(np.logical_and.reduce((T.PLnumber==TPLnr[n], T.tPL >= tPLstart[n], T.tPL <= tPLend[n])))[0]
    
    #the start and end time of the PL
    tim0= int(T.tPL[index][0])
    tim1= int(T.tPL[index][-1]+ 1)
    
    #import ERA data
    d= data(['SST', 'T', 'MSLP', 'sHum', 'uPV', 'vPV', 'thetaPV', 'Water', 'PBH', 'T850', 'sHum850', 'Geop', 'uPL500', 'vPL500', 'p_PV', 'CAPE'], yearL[n], monthL[n], tstart= tim0, tend= tim1)
    d.SST[d.SST<272]= np.nan #this sets the area with ice to land
    d.SST= np.ma.array(d.SST, mask= np.isnan(d.SST))
    d.imp_u10()
    UPV= np.sqrt(d.uPV**2+ d.vPV**2)
    UPL500= np.sqrt(d.uPL500**2+ d.vPL500**2)

    
#    #the lat lon coordinatinates of the PL for ERA
    Tlatraw= np.array(np.round(-(T.lat[index]-d.lat[0])*2), dtype= int)
    Tlonraw= np.array(np.round((T.lon[index] - d.lon[0])*2), dtype= int)
    pts2D=(list(Tlatraw), list(Tlonraw))  #these are the cyclone points in ERA data coordinates
   
         
    for i, indexi in enumerate(index): #loop through all time steps of the given TPLnr

        
        latdist3, latdist2, latdist1= 6, 4, 2
        latdist5, latdist4= 10, 8
        
        londist1= int(latdist1//np.cos(np.deg2rad(T.lat[indexi])))
        londist2= int(latdist2//np.cos(np.deg2rad(T.lat[indexi])))
        londist3= int(latdist3//np.cos(np.deg2rad(T.lat[indexi])))
        londist4= int(latdist4//np.cos(np.deg2rad(T.lat[indexi])))
        londist5= int(latdist5//np.cos(np.deg2rad(T.lat[indexi])))

        mask= radMask2D(d.uPV.shape[1:], (pts2D[0][i], pts2D[1][i]), radiusx= latdist1, radiusy= londist1) #lat lon mask
        maskc= [i, mask] #tim, lat, lon mask (mask complete)
        mask2= radMask2D(d.SST.shape[1:], (pts2D[0][i], pts2D[1][i]), radiusx= latdist2, radiusy= londist2) #lat lon mask
        maskc2= [i, mask2] #tim, lat, lon mask (mask complete)
        mask3= radMask2D(d.SST.shape[1:], (pts2D[0][i], pts2D[1][i]), radiusx= latdist3, radiusy= londist3) #lat lon mask
        maskc3= [i, mask3] #tim, lat, lon mask (mask complete)
        mask4= radMask2D(d.SST.shape[1:], (pts2D[0][i], pts2D[1][i]), radiusx= latdist4, radiusy= londist4) #lat lon mask
        maskc4= [i, mask4] #tim, lat, lon mask (mask complete)
        mask5= radMask2D(d.SST.shape[1:], (pts2D[0][i], pts2D[1][i]), radiusx= latdist5, radiusy= londist5) #lat lon mask
        maskc5= [i, mask5] #tim, lat, lon mask (mask complete)
       
        landistr= len(np.where(d.SST[maskc2].mask == True)[0]) /len(d.SST[maskc2])
#        print(T.PLnumber[indexi], i, 'landdistr:' ,landistr)
        
        if landistr < .25:  #true if we are over ocean and d.SST[maskc] has a value
            statvort += [T.vort[indexi]]
            statSPL += [SPLnr[n]]
            statTPL += [yearL[n]*1E7+ monthL[n]*1E5 +TPLnr[n]]
            
            """calculate potential temperature and equivalent potential temperature"""
            RCp= 2/7 #R/Cp
            theta850= d.T850[maskc]*(1000/850)**RCp
            theta700= d.T[:,1][maskc]*(1000/700)**RCp
            theta700_2= d.T[:,1][maskc2]*(1000/700)**RCp
            theta700_3= d.T[:,1][maskc3]*(1000/700)**RCp
                         
            theta500= d.T[:,2][maskc]*(1000/500)**RCp
            theta500_2= d.T[:,2][maskc2]*(1000/500)**RCp
            theta500_3= d.T[:,2][maskc3]*(1000/500)**RCp
                         
            thetaSST= d.SST[maskc]*(1000/d.MSLP[maskc])**RCp
            thetaSST_2= d.SST[maskc2]*(1000/d.MSLP[maskc2])**RCp
            thetaSST_3= d.SST[maskc3]*(1000/d.MSLP[maskc3])**RCp

            Lc, Cp= 2501, 1                
            theta_eSST= thetaSST* np.exp(d.sHum[:,0][maskc]* Lc/(d.SST[maskc]* Cp))
            theta_eSST_2= thetaSST_2* np.exp(d.sHum[:,0][maskc2]* Lc/(d.SST[maskc2]* Cp))
            theta_eSST_3= thetaSST_3* np.exp(d.sHum[:,0][maskc3]* Lc/(d.SST[maskc3]* Cp))
            
            theta_e850= theta850* np.exp(d.sHum850[maskc]* Lc/(d.T850[maskc]* Cp))
            theta_e700= theta700* np.exp(d.sHum[:,1][maskc]* Lc/(d.T[:,1][maskc]* Cp))
            theta_e700_2= theta700_2* np.exp(d.sHum[:,1][maskc2]* Lc/(d.T[:,1][maskc2]* Cp))
            theta_e700_3= theta700_3* np.exp(d.sHum[:,1][maskc3]* Lc/(d.T[:,1][maskc3]* Cp))
            
            theta_e500= theta500* np.exp(d.sHum[:,2][maskc]* Lc/(d.T[:,2][maskc]* Cp))
           
            Theta_e850= EquiPotTemp(d.T850[i], d.sHum850[i], plev= 850)
            latd= 55.
            lond= np.tile(latd*np.cos(np.deg2rad(d.lat)), (720,1)).T         
            g= np.gradient(Theta_e850, latd, lond) 
            absg= np.sqrt(g[0]**2+ g[1]**2) *100 #change in Theta_e /100km
            
            statgradtheta_e1 += [np.max(absg[mask])]
            statgradtheta_e2 += [np.max(absg[mask2])]
            statgradtheta_e3 += [np.max(absg[mask3])]
            statgradtheta_e5 += [np.max(absg[mask5])]
            
            statMSLP += [d.MSLP[i, pts2D[0][i], pts2D[1][i] ] ]

            statpdiff1 += [d.MSLP[i, pts2D[0][i], pts2D[1][i] ] - np.mean(d.MSLP[maskc]) ]
            statpdiff2 += [d.MSLP[i, pts2D[0][i], pts2D[1][i] ] - np.mean(d.MSLP[maskc2]) ]
            statpdiff3 += [d.MSLP[i, pts2D[0][i], pts2D[1][i] ] - np.mean(d.MSLP[maskc3]) ]            
            statpdiff5 += [d.MSLP[i, pts2D[0][i], pts2D[1][i] ] - np.mean(d.MSLP[maskc5]) ]            

            statstab500 += [np.mean(d.SST[maskc]- d.T[:,2][maskc])]
            statstab700 += [np.mean(d.SST[maskc]- d.T[:,1][maskc])]
            
            stattheta500 += [np.mean(thetaSST- theta500)]
            stattheta500_2 += [np.mean(thetaSST_2- theta500_2)]
            stattheta500_3 += [np.mean(thetaSST_3- theta500_3)]
            stattheta500_max += [np.max(thetaSST- theta500)]
            stattheta500_max_2 += [np.max(thetaSST_2- theta500_2)]
            stattheta500_max_3 += [np.max(thetaSST_3- theta500_3)]
            
            stattheta700 += [np.mean(thetaSST- theta700)]
            stattheta_e500 += [np.mean(theta_eSST- theta_e500)]
            stattheta_e700 += [np.mean(theta_eSST- theta_e700)]
            stattheta_e700_2 += [np.mean(theta_eSST_2- theta_e700_2)]
            stattheta_e700_3 += [np.mean(theta_eSST_3- theta_e700_3)]
            stattheta_e700_max += [np.max(theta_eSST- theta_e700)]
            stattheta_e700_max_2 += [np.max(theta_eSST_2- theta_e700_2)]
            stattheta_e700_max_3 += [np.max(theta_eSST_3- theta_e700_3)]

            stattheta850_1 += [np.mean(thetaSST- theta850)]
            stattheta850_2 += [np.mean(thetaSST_2- PotTemp(d.T850[maskc2], 850) )]
            stattheta850_3 += [np.mean(thetaSST_3- PotTemp(d.T850[maskc3], 850) )]
            stattheta850_max_1 += [np.max(thetaSST- theta850)]
            stattheta850_max_2 += [np.max(thetaSST_2- PotTemp(d.T850[maskc2], 850) )]
            stattheta_e850_1 += [np.mean(theta_eSST- theta_e850)]
            stattheta_e850_2 += [np.mean(theta_eSST_2- Theta_e850[mask2])]
            stattheta_e850_3 += [np.mean(theta_eSST_3- Theta_e850[mask3])]
            stattheta_e850_max_1 += [np.max(theta_eSST- theta_e850)]
            stattheta_e850_max_2 += [np.max(theta_eSST_2-  Theta_e850[mask2])]

            L= 7.5e5
            statMCAO_Kolstad+= [np.max(L/(d.Geop[:,1][maskc]/9.81)* (np.log(thetaSST)- np.log(theta700)) )]
            statMCAO_Kolstad_0+= [np.max((thetaSST- theta700)/(d.MSLP[maskc] - 700) )]
            statMCAO_Kolstad_0_500+= [np.max((thetaSST- theta500)/(d.MSLP[maskc] - 500) )]


            statt_500 += [np.mean(d.T[:,2][maskc])]
            statt_700 += [np.mean(d.T[:,1][maskc])]
            statSST += [np.mean(d.SST[maskc])]
            statattheta_e700 += [np.mean(theta_e700)]
            statattheta_e850 += [np.mean(theta_e850)]

            statU10 += [np.max(np.sqrt(d.u10[maskc]**2+ d.v10[maskc]**2))]
            statU10_2 += [np.max(np.sqrt(d.u10[maskc2]**2+ d.v10[maskc2]**2))]
            statU10_3 += [np.max(np.sqrt(d.u10[maskc3]**2+ d.v10[maskc3]**2))]
 
            statCAPE += [np.mean(d.CAPE[maskc])]
            statCAPE_2 += [np.mean(d.CAPE[maskc2])]
            statCAPE_3 += [np.mean(d.CAPE[maskc3])]
            statCAPE_max += [np.max(d.CAPE[maskc])]
            statCAPE_max_2 += [np.max(d.CAPE[maskc2])]
            statCAPE_max_3 += [np.max(d.CAPE[maskc3])]
           
            statUPV += [np.max(UPV[maskc])]     
            
            
            statUPVnorth += [np.max(UPV[i, :Tlatraw[i]+1, Tlonraw[i]])]
            statUPVnorthM2p5 += [np.max(UPV[i, :Tlatraw[i]+6, Tlonraw[i]])]
            statUPVnorthM5p0 += [np.max(UPV[i, :Tlatraw[i]+11, Tlonraw[i]])]
            statUPVnorthM7p5 += [np.max(UPV[i, :Tlatraw[i]+16, Tlonraw[i]])]

            statUPL500north += [np.max(UPL500[i, :Tlatraw[i]+1, Tlonraw[i]])]
            statgradtheta_e850north += [np.max(absg[:Tlatraw[i]+1, Tlonraw[i]])]

            statWater += [np.mean(d.Water[maskc])]
            statPBH += [np.mean(d.PBH[maskc3])]
            statp_PV_3 += [np.max(d.p_PV[maskc3])]
            statp_PV_4 += [np.max(d.p_PV[maskc4])]
            statp_PV_5 += [np.max(d.p_PV[maskc5])]     
            
            statthetaPV_1 += [np.mean(d.thetaPV[maskc])]
            statthetaPV_3 += [np.mean(d.thetaPV[maskc3])]
            statthetaPV_5 += [np.mean(d.thetaPV[maskc5])]
            statthetaPV_1min += [np.min(d.thetaPV[maskc])]
            statthetaPV_3min += [np.min(d.thetaPV[maskc3])]
            statthetaPV_5min += [np.min(d.thetaPV[maskc5])]

            statthetaSSTmPV_1 += [np.mean(thetaSST- d.thetaPV[maskc])]
            statthetaSSTmPV_3 += [np.mean(thetaSST_3- d.thetaPV[maskc3])]
            statthetaSSTmPV_1_max += [np.max(thetaSST- d.thetaPV[maskc])]
            statthetaSSTmPV_3_max += [np.max(thetaSST_3- d.thetaPV[maskc3])]


print('nr of PLs: ', len(SPLnr))
print('nr of PLs afterwards: ', len(remove_dublicate(statSPL)))