#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 16:38:42 2017

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

import numpy as np
import sys

sys.path.insert(0, '/home/'+user+'/polar_low_code/Functions/')
from f_imp_ERA2 import *
from f_impLists import *
from f_meteo import *

import heapq
#import datetime

import time
start_time = time.time()



"""creation of cyclone statistics"""
stat_track_U10, stat_track_vort= [], []
stat_track_stab500, stat_track_stab700= [], []
stat_track_theta700, stat_track_theta500, stat_track_theta_e700, stat_track_theta_e500=  [], [], [], []
stat_track_theta500_2, stat_track_theta500_3, stat_track_theta500_max, stat_track_theta500_max_2, stat_track_theta500_max_3= [], [], [], [], []
stat_track_theta_e700_2, stat_track_theta_e700_3, stat_track_theta_e700_max, stat_track_theta_e700_max_2, stat_track_theta_e700_max_3= [], [], [], [], []

stat_track_theta850_1, stat_track_theta850_2, stat_track_theta850_3, stat_track_theta850_max_1, stat_track_theta850_max_2= [], [], [], [], []
stat_track_theta_e850_1, stat_track_theta_e850_2, stat_track_theta_e850_3, stat_track_theta_e850_max_1, stat_track_theta_e850_max_2= [], [], [], [], []

stat_track_U10_2, stat_track_U10_3= [], []

stat_track_UPV=  []
stat_track_Water= []
stat_track_t_500, stat_track_t_700, stat_track_SST = [], [], []
stat_track_UPVnorth = []
stat_track_UPVnorthM2p5, stat_track_UPVnorthM5p0, stat_track_UPVnorthM7p5= [], [], []
stat_track_UPL500north = []
stat_track_gradtheta_e850north= []

stat_track_p_PV_3, stat_track_p_PV_4, stat_track_p_PV_5 = [], [], []
stat_track_PBH= []
stat_track_attheta_e700, stat_track_attheta_e850= [], [] #the equivalent potential temperature at 700 and 850 hPa
stat_track_MCAO_Kolstad = []
stat_track_MCAO_Kolstad_0 = []
stat_track_MCAO_Kolstad_0_500 = []
                                                  
stat_track_MSLP= []
stat_track_gradtheta_e1, stat_track_gradtheta_e2, stat_track_gradtheta_e3 , stat_track_gradtheta_e5= [], [], [], []
stat_track_pdiff1, stat_track_pdiff2, stat_track_pdiff3, stat_track_pdiff5 =  [], [], [], []                                                  

stat_track_thetaPV_1, stat_track_thetaPV_3, stat_track_thetaPV_5, stat_track_thetaPV_1min, stat_track_thetaPV_3min, stat_track_thetaPV_5min= [], [], [], [], [], []
stat_track_thetaSSTmPV_1, stat_track_thetaSSTmPV_3, stat_track_thetaSSTmPV_1_max, stat_track_thetaSSTmPV_3_max=  [], [], [], []
stat_track_CAPE, stat_track_CAPE_2, stat_track_CAPE_3, stat_track_CAPE_max, stat_track_CAPE_max_2, stat_track_CAPE_max_3= [], [], [], [], [], []

stat_track_TPL= []

for year in np.arange(2003, 2004):
    for month in [1, 2, 3, 4, 10, 11, 12]:

        print(year, month)
        T=TRACK_list_month2(year, month)
        T.local(minlon= -180, maxlon= 180, minlat= 30.5, maxlat= 84.5) #has to be within ERA data
,
        print("--- %s seconds ---" % (time.time() - start_time))

        """import data"""
        d= data(['SST','uPV', 'vPV', 'thetaPV', 'T',  'MSLP', 'sHum', 'Water', 'PBH', 'T850', 'sHum850', 'uPL500', 'vPL500', 'Geop', 'p_PV', 'CAPE'], year, month)
#        d= data(['SST', 'T',  'MSLP'], year, month)

        d.SST[d.SST<272]= np.nan #this sets the area with ice to land
        d.SST= np.ma.array(d.SST, mask= np.isnan(d.SST))
        d.imp_u10()

#        print("--- %s seconds ---" % (time.time() - start_time))
        
        U10= np.sqrt(d.u10**2+ d.v10**2)
        UPV= np.sqrt(d.uPV**2+ d.vPV**2)
        UPL500= np.sqrt(d.uPL500**2+ d.vPL500**2)
    
        print('number of cyclones: ', len(remove_dublicate(T.PLlist[0])))
        
        
        """remove land cyclones and get all cyclones to consider"""
        Tlatraw= np.array(np.round(-(T.PLlist[3]-d.lat[0])*2), dtype= int)
        Tlonraw= np.array(np.round((T.PLlist[2] - d.lon[0])*2), dtype= int)%720

        pts=(list(np.array(T.PLlist[1], dtype= int)), list(Tlatraw), list(Tlonraw))  #these are the cyclone points in ERA data coordinates

        latdist5, latdist4, latdist3, latdist2, latdist1= 10, 8, 6, 4, 2

        OceanPL=remove_dublicate(T.PLlist[0][d.SST[pts]> 272])

        print('number of ocean PLs; ', len(OceanPL))
        
        Tlatraw= [Tlatraw[x] for x, PLn in enumerate(T.PLlist[0]) if PLn in OceanPL]
        Tlonraw= [Tlonraw[x] for x, PLn in enumerate(T.PLlist[0]) if PLn in OceanPL]

        T.PLlist= np.array([T.PLlist[:,x] for x, PLn in enumerate(T.PLlist[0]) if PLn in OceanPL]).T
        pts= (list(np.array(T.PLlist[1], dtype= int)), Tlatraw, Tlonraw)                          


        """calculate properties for every cyclone """                          
        
        for PLnr in OceanPL:  #loop through all Ocean cyclones, PLnr is the current cyclone
            index= np.where(T.PLlist[0]==PLnr)[0] #the index of the current cyclone
#            print(index)

            for i in index:  #loop through every time step of the current cyclone


                #calculate the distance in which the values are taken in pixels
                londist5= int(latdist5//np.cos(np.deg2rad(T.PLlist[3][i])))
                londist4= int(latdist4//np.cos(np.deg2rad(T.PLlist[3][i])))
                londist3= int(latdist3//np.cos(np.deg2rad(T.PLlist[3][i])))
                londist2= int(latdist2//np.cos(np.deg2rad(T.PLlist[3][i])))
                londist1= int(latdist1//np.cos(np.deg2rad(T.PLlist[3][i])))           

                mask= radMask2D(d.SST.shape[1:], (pts[1][i], pts[2][i]), radiusx= latdist1, radiusy= londist1) #lat lon mask
                maskc= [pts[0][i], mask] #tim, lat, lon mask (mask complete)
                mask2= radMask2D(d.SST.shape[1:], (pts[1][i], pts[2][i]), radiusx= latdist2, radiusy= londist2) #lat lon mask
                maskc2= [pts[0][i], mask2] #tim, lat, lon mask (mask complete)
                mask3= radMask2D(d.SST.shape[1:], (pts[1][i], pts[2][i]), radiusx= latdist3, radiusy= londist3) #lat lon mask
                maskc3= [pts[0][i], mask3] #tim, lat, lon mask (mask complete)
                mask4= radMask2D(d.SST.shape[1:], (pts[1][i], pts[2][i]), radiusx= latdist4, radiusy= londist4) #lat lon mask
                maskc4= [pts[0][i], mask4] #tim, lat, lon mask (mask complete)
                mask5= radMask2D(d.SST.shape[1:], (pts[1][i], pts[2][i]), radiusx= latdist5, radiusy= londist5) #lat lon mask
                maskc5= [pts[0][i], mask5] #tim, lat, lon mask (mask complete)

                landistr= len(np.where(d.SST[maskc2].mask == True)[0]) /len(d.SST[maskc2])
#                print(T.PLlist[0][i], 'landdistr:' ,landistr)


                if landistr < .25:  #true if we are over ocean and d.SST[maskc] has a value
                    stat_track_TPL += [year*1E7 + month*1E5 + T.PLlist[0][i]]
                    stat_track_vort += [T.PLlist[4][i]]
                
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
                                     
                    Lc, Cp= 2501E3, 1006
                    theta_eSST=thetaSST* np.exp(d.sHum[:,0][maskc]* Lc/(d.SST[maskc]* Cp))
                    theta_eSST_2=thetaSST_2* np.exp(d.sHum[:,0][maskc2]* Lc/(d.SST[maskc2]* Cp))
                    theta_eSST_3=thetaSST_3* np.exp(d.sHum[:,0][maskc3]* Lc/(d.SST[maskc3]* Cp))
                    
                    theta_e850= theta850* np.exp(d.sHum850[maskc]* Lc/(d.T850[maskc]* Cp))                    
                    theta_e700= theta700* np.exp(d.sHum[:,1][maskc]* Lc/(d.T[:,1][maskc]* Cp))
                    theta_e700_2= theta700_2* np.exp(d.sHum[:,1][maskc2]* Lc/(d.T[:,1][maskc2]* Cp))
                    theta_e700_3= theta700_3* np.exp(d.sHum[:,1][maskc3]* Lc/(d.T[:,1][maskc3]* Cp))
                    
                    theta_e500= theta500* np.exp(d.sHum[:,2][maskc]* Lc/(d.T[:,2][maskc]* Cp))
                    
                    Theta_e850= EquiPotTemp(d.T850[pts[0][i]], d.sHum850[pts[0][i]], plev= 850) #field without mask
                    
                    latd= 55.
                    lond= np.tile(latd*np.cos(np.deg2rad(d.lat)), (720,1)).T         
                    g= np.gradient(Theta_e850, latd, lond) 
                    absg= np.sqrt(g[0]**2+ g[1]**2) *100 #change in Theta_e /100km
                    
                    stat_track_gradtheta_e1 += [np.max(absg[mask])]
                    stat_track_gradtheta_e2 += [np.max(absg[mask2])]
                    stat_track_gradtheta_e3 += [np.max(absg[mask3])]
                    stat_track_gradtheta_e5 += [np.max(absg[mask5])]
                    
                    stat_track_MSLP += [d.MSLP[pts[0][i], pts[1][i], pts[2][i] ] ]
    
                    stat_track_pdiff1 += [d.MSLP[pts[0][i], pts[1][i], pts[2][i] ] - np.mean(d.MSLP[maskc]) ]
                    stat_track_pdiff2 += [d.MSLP[pts[0][i], pts[1][i], pts[2][i] ] - np.mean(d.MSLP[maskc2]) ]
                    stat_track_pdiff3 += [d.MSLP[pts[0][i], pts[1][i], pts[2][i] ] - np.mean(d.MSLP[maskc3]) ]            
                    stat_track_pdiff5 += [d.MSLP[pts[0][i], pts[1][i], pts[2][i] ] - np.mean(d.MSLP[maskc5]) ] 
              
                    stat_track_stab500 += [np.mean(d.SST[maskc]- d.T[:,2][maskc])]
                    stat_track_stab700 += [np.mean(d.SST[maskc]- d.T[:,1][maskc])]
                    stat_track_theta500 += [np.mean(thetaSST- theta500)]
                    stat_track_theta500_2 += [np.mean(thetaSST_2- theta500_2)]
                    stat_track_theta500_3 += [np.mean(thetaSST_3- theta500_3)]
                    stat_track_theta500_max += [np.max(thetaSST- theta500)]
                    stat_track_theta500_max_2 += [np.max(thetaSST_2- theta500_2)]
                    stat_track_theta500_max_3 += [np.max(thetaSST_3- theta500_3)]
                    stat_track_theta700 += [np.mean(thetaSST- theta700)]
                    stat_track_theta_e500 += [np.mean(theta_eSST- theta_e500)]
                    stat_track_theta_e700 += [np.mean(theta_eSST- theta_e700)]
                    stat_track_theta_e700_2 += [np.mean(theta_eSST_2- theta_e700_2)]
                    stat_track_theta_e700_3 += [np.mean(theta_eSST_3- theta_e700_3)]
                    stat_track_theta_e700_max += [np.max(theta_eSST- theta_e700)]
                    stat_track_theta_e700_max_2 += [np.max(theta_eSST_2- theta_e700_2)]
                    stat_track_theta_e700_max_3 += [np.max(theta_eSST_3- theta_e700_3)]


                    stat_track_theta850_1 += [np.mean(thetaSST- theta850)]
                    stat_track_theta850_2 += [np.mean(thetaSST_2- PotTemp(d.T850[maskc2], 850) )]
                    stat_track_theta850_3 += [np.mean(thetaSST_3- PotTemp(d.T850[maskc3], 850) )]
                    stat_track_theta850_max_1 += [np.max(thetaSST- theta850)]
                    stat_track_theta850_max_2 += [np.max(thetaSST_2- PotTemp(d.T850[maskc2], 850) )]
                    stat_track_theta_e850_1 += [np.mean(theta_eSST- theta_e850)]
                    stat_track_theta_e850_2 += [np.mean(theta_eSST_2- Theta_e850[mask2])]
                    stat_track_theta_e850_3 += [np.mean(theta_eSST_3- Theta_e850[mask3])]
                    stat_track_theta_e850_max_1 += [np.max(theta_eSST- theta_e850)]
                    stat_track_theta_e850_max_2 += [np.max(theta_eSST_2-  Theta_e850[mask2])]

                    
                    L= 7.5e5
                    stat_track_MCAO_Kolstad+= [np.max(L/(d.Geop[:,1][maskc]/9.81)* (np.log(thetaSST)- np.log(theta700)) )]
                    stat_track_MCAO_Kolstad_0+= [np.max((thetaSST- theta700)/(d.MSLP[maskc] - 700) )]
                    stat_track_MCAO_Kolstad_0_500+= [np.max((thetaSST- theta500)/(d.MSLP[maskc] - 500) )]
            
                    stat_track_t_500 += [np.mean(d.T[:,2][maskc])]
                    stat_track_t_700 += [np.mean(d.T[:,1][maskc])]
                    stat_track_SST += [np.mean(d.SST[maskc])]
                    stat_track_attheta_e700 += [np.mean(theta_e700)]
                    stat_track_attheta_e850 += [np.mean(theta_e850)]
    
                    stat_track_U10 += [np.max(np.sqrt(d.u10[maskc]**2+ d.v10[maskc]**2))]
                    stat_track_U10_2 += [np.max(np.sqrt(d.u10[maskc2]**2+ d.v10[maskc2]**2))]
                    stat_track_U10_3 += [np.max(np.sqrt(d.u10[maskc3]**2+ d.v10[maskc3]**2))]

                    stat_track_CAPE += [np.mean(d.CAPE[maskc])]
                    stat_track_CAPE_2 += [np.mean(d.CAPE[maskc2])]
                    stat_track_CAPE_3 += [np.mean(d.CAPE[maskc3])]
                    stat_track_CAPE_max += [np.max(d.CAPE[maskc])]
                    stat_track_CAPE_max_2 += [np.max(d.CAPE[maskc2])]
                    stat_track_CAPE_max_3 += [np.max(d.CAPE[maskc3])]
                    
                    stat_track_UPV += [np.max(UPV[maskc])]
                    stat_track_UPVnorth += [np.max(UPV[pts[0][i], :Tlatraw[i]+1, Tlonraw[i]])]
                    
                    Tlatmax= np.min([Tlatraw[i]+6, 109])
                    stat_track_UPVnorthM2p5 += [np.max(UPV[pts[0][i], :Tlatmax, Tlonraw[i]])]
                    Tlatmax= np.min([Tlatraw[i]+11, 109])
                    stat_track_UPVnorthM5p0 += [np.max(UPV[pts[0][i], :Tlatmax, Tlonraw[i]])]
                    Tlatmax= np.min([Tlatraw[i]+16, 109])
                    stat_track_UPVnorthM7p5 += [np.max(UPV[pts[0][i], :Tlatmax, Tlonraw[i]])]                    
                    
                    stat_track_UPL500north += [np.max(UPL500[pts[0][i], :Tlatraw[i]+1, Tlonraw[i]])]
                    stat_track_gradtheta_e850north += [np.max(absg[:Tlatraw[i]+1, Tlonraw[i]])]
                    
                    stat_track_Water += [np.mean(d.Water[maskc])]
                    stat_track_PBH += [np.mean(d.PBH[maskc3])]

                    stat_track_p_PV_3 += [np.max(d.p_PV[maskc3])]
                    stat_track_p_PV_4 += [np.max(d.p_PV[maskc4])]
                    stat_track_p_PV_5 += [np.max(d.p_PV[maskc5])] 

                    stat_track_thetaPV_1 += [np.mean(d.thetaPV[maskc])]
                    stat_track_thetaPV_3 += [np.mean(d.thetaPV[maskc3])]
                    stat_track_thetaPV_5 += [np.mean(d.thetaPV[maskc5])]
                    stat_track_thetaPV_1min += [np.min(d.thetaPV[maskc])]
                    stat_track_thetaPV_3min += [np.min(d.thetaPV[maskc3])]
                    stat_track_thetaPV_5min += [np.min(d.thetaPV[maskc5])]

                    stat_track_thetaSSTmPV_1 += [np.mean(thetaSST- d.thetaPV[maskc])]
                    stat_track_thetaSSTmPV_3 += [np.mean(thetaSST_3- d.thetaPV[maskc3])]
                    stat_track_thetaSSTmPV_1_max += [np.max(thetaSST- d.thetaPV[maskc])]
                    stat_track_thetaSSTmPV_3_max += [np.max(thetaSST_3- d.thetaPV[maskc3])]


print("--- %s seconds ---" % (time.time() - start_time))