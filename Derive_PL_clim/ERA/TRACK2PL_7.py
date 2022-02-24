#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 15:06:41 2017

@author: pst019
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import sys  #to import the functions from a different directory
import pickle

sys.path.insert(0, '/home/'+user+'/polar_low_code/Functions/')
from f_imp_ERA2 import *
from f_impLists import *
from f_meteo import * #meteorological functions
from f_useful import *

import scipy.ndimage.filters as filters

import csv

import datetime
import time
start_time = time.time()


#hemi='NH'
#approved= True

hemi='SH'
approved= False

"""10percentile - 3 crit"""
#Tdurationlim=2 #TRACK duration limit in 3hourly time steps - higher less PLs
#vortlim= 5.04 #the higher the less PLs
#thetastablim= np.float64(-9.403) #now it is the max value in circle of latitude 1 #the higher the less PLs
#
##windlim= 14 #the higher the less PLs
#polarfrontlim= 31.29 #the higher the more PLs
##jetlim=40

"""10percentile - 3 crit version 7"""
versionname= '7'
Tdurationlim=2 #TRACK duration limit in 3hourly time steps - higher less PLs
vortlim= None #the higher the less PLs
thetastablim= np.float64(-9.403) #now it is the max value in circle of latitude 1 #the higher the less PLs
radstab= 1*2

pdifflim= np.float(-0.409)
radpdiff= 1*2

windlim= None #14 #the higher the less PLs
polarfrontlim= 31.29 #the higher the more PLs

theta_e700lim= None
radtheta_e700= 2*2
gradtheta_elim= None
radgradtheta_e= 5*2

"""5pcentile version 6_2""" #to make this possible more variables have to be imported for the SH
#versionname = '6_2'
#Tdurationlim=2 #TRACK duration limit in 3hourly time steps - higher less PLs
#vortlim= 4.3535 #the higher the less PLs
#
#thetastablim= np.float64(-8.540) #now it is the max value in circle of latitude 3 #the higher the less PLs
#radstab= 3*2
#
#windlim= None  #the higher the less PLs
#polarfrontlim= np.float(35.808) #the higher the more PLs
#pdifflim= np.float(-0.2902)
#radpdiff= 1*2
#
#theta_e700lim= np.float(4.891)
#radtheta_e700= 2*2
#
#gradtheta_elim= np.float(8.390)
#radgradtheta_e= 5*2

cnr, vnr, lnr, pnr, jnr, snr, pdnr, gtnr, s700nr= 0,0,0,0,0,0,0,0,0

for year in np.arange(1979, 2016):
    print(year)
    if hemi== 'NH':    cPLlist= np.zeros((14,0)) #complete PLlist (year, month, TrackPLnr, tPL, lon, lat, PLpoint(yes= 1, no= 0))
    elif hemi== 'SH':    cPLlist= np.zeros((11,0)) #complete PLlist (year, month, TrackPLnr, tPL, lon, lat, PLpoint(yes= 1, no= 0))

    
    for month in np.arange(1, 13): #[1, 2, 3, 4, 5, 9, 10, 11, 12]:
        shearlist, stablist, UPVpolewlist, pdifflist, gradtheta_elist, theta_e700list =[], [], [], [], [], []
    
        """import data"""
        if gradtheta_elim == None and theta_e700lim == None:
            d= data(['SST', 'T', 'uPV', 'vPV', 'MSLP'], year, month, hemisp= hemi)
        else:
            d= data(['SST', 'T', 'uPV', 'vPV', 'MSLP', 'T850', 'sHum850', 'sHum'], year, month, hemisp= hemi)
        d.SST[d.SST<272]= np.nan #this sets the area with ice to land
        d.SST= ma.array(d.SST, mask= isnan(d.SST))        

        UPV= np.sqrt(d.uPV**2+ d.vPV**2)
    
        T=TRACK_list_month2(year, month, Tdurationlim, hemisp= hemi, approved= approved) #maybe have to add approved for SH
        T.local(minlon= d.lon[0], maxlon= d.lon[-1], minlat= d.lat[-1], maxlat=d.lat[0])
        #PLlist= [PLnumber, tPL, lon, lat, vort]

        cnr += len(remove_dublicate(T.PLlist[0]))
#        print('number of cyclones: '+str(len(remove_dublicate(T.PLlist[0]))))
            
        
        """extract PLs from cyclones """  
        if vortlim != None:
            PLlist= T.PLlist[:, T.PLlist[4] >= vortlim]
        else:
            PLlist= T.PLlist
            
        vnr += len(remove_dublicate(PLlist[0]))
#        print('after vorticity removal: '+str(len(remove_dublicate(PLlist[0]))))
                 

        """exclude land cyclones (new version - not applied and tested yet)"""
        Tlatraw= np.array(np.round(-(PLlist[3]-d.lat[0])*2), dtype= int)
        Tlonraw= np.array(np.round((PLlist[2] - d.lon[0])*2), dtype= int)%720
        pts=(list(PLlist[1]), list(Tlatraw), list(Tlonraw)) #the raw data points in d format

        landistrm= []
        latdist2= 4
        for i, t in enumerate(PLlist[1]):  #loop through every PL point            
            londist2= int(latdist2//np.cos(np.deg2rad(PLlist[3][i])))           
            mask2= radMask2D((len(d.lat), len(d.lon)), (pts[1][i], pts[2][i]), radiusx= latdist2, radiusy= londist2) #lat lon mask

            landistrm += [len(np.where(d.SST[t][mask2].mask == True)[0]) /len(d.SST[t][mask2])]

        PLlist= PLlist[:, np.array(landistrm) < .25]

        lnr += len(remove_dublicate(PLlist[0]))
 

        """exclude polar front cyclones"""
        if polarfrontlim != None:
            Tlatraw= np.array(np.round(-(PLlist[3]-d.lat[0])*2), dtype= int)
            Tlonraw= np.array(np.round((PLlist[2] - d.lon[0])*2), dtype= int)%720
    #        pts=(list(PLlist[1]), list(Tlatraw), list(Tlonraw)) #the raw data points in d format
            
            if hemi== 'NH':
                UPVpolarm = np.array([np.max(UPV[PLlist[1][i], : int(Tlatraw[i]+1), int(Tlonraw[i])]) for i in range(len(PLlist[0]))])
            elif hemi== 'SH':
                UPVpolarm = np.array([np.max(UPV[PLlist[1][i], int(Tlatraw[i]):, int(Tlonraw[i])]) for i in range(len(PLlist[0]))])

    #        PLlist=  np.append(PLlist, UPVnorthm.reshape(1, -1), axis= 0) #append UPVnorth to the list
            PLlist= PLlist[:, UPVpolarm< polarfrontlim]
            
            pnr +=len(remove_dublicate(PLlist[0]))       
    #        print('after polarfront removal: '+str(len(remove_dublicate(PLlist[0]))))


        """exclude theta stability cyclones"""
        if thetastablim != None:
            Tlatraw= np.array(np.round(-(PLlist[3]-d.lat[0])*2), dtype= int)
            Tlonraw= np.array(np.round((PLlist[2] - d.lon[0])*2), dtype= int)%720
            pts=(list(PLlist[1]), list(Tlatraw), list(Tlonraw)) #the raw data points in d format
            
            theta_500m= []
            for i, t in enumerate(PLlist[1]):  #loop through every PL point            
                londist= int(radstab//np.cos(np.deg2rad(PLlist[3][i])))           
                mask= radMask2D((len(d.lat), len(d.lon)), (pts[1][i], pts[2][i]), radiusx= radstab, radiusy= londist) #lat lon mask
    
                """calculate potential temperature and equivalent potential temperature"""
                RCp = 2/7 #R/Cp
                theta500= d.T[t,2][mask]*(1000/500)**RCp
                thetaSST= d.SST[t][mask]*(1000/d.MSLP[t][mask])**RCp
                             
    #            theta_500m += [np.mean(thetaSST- theta500)] #the mean in radius
                theta_500m += [np.max(thetaSST- theta500)] #the max in radius
    
    
            PLlist=  np.append(PLlist, np.array(theta_500m).reshape(1, -1), axis= 0) #append the stability to the list
            PLlist= PLlist[:, np.array(theta_500m) > thetastablim]
            
            snr += len(remove_dublicate(PLlist[0]))        


        """exclude pdiff cyclones"""
        if pdifflim != None:
            Tlatraw= np.array(np.round(-(PLlist[3]-d.lat[0])*2), dtype= int)
            Tlonraw= np.array(np.round((PLlist[2] - d.lon[0])*2), dtype= int)%720
            pts=(list(PLlist[1]), list(Tlatraw), list(Tlonraw)) #the raw data points in d format
            
            pdiffm= []
            for i, t in enumerate(PLlist[1]):  #loop through every PL point            
                londistpdiff= int(radpdiff//np.cos(np.deg2rad(PLlist[3][i])))           
                mask= radMask2D((len(d.lat), len(d.lon)), (pts[1][i], pts[2][i]), radiusx= radpdiff, radiusy= londistpdiff) #lat lon mask
    
                """calculate the depth of the low"""
                pdiffm += [d.MSLP[t, pts[1][i], pts[2][i] ] - np.mean(d.MSLP[t][mask]) ] #the difference between point and mean
    
    
            PLlist=  np.append(PLlist, np.array(pdiffm).reshape(1, -1), axis= 0) #append the stability to the list
            PLlist= PLlist[:, np.array(pdiffm) < pdifflim]
            
            pdnr += len(remove_dublicate(PLlist[0]))      


        """exclude gradtheta_e """
        if gradtheta_elim != None:
            Tlatraw= np.array(np.round(-(PLlist[3]-d.lat[0])*2), dtype= int)
            Tlonraw= np.array(np.round((PLlist[2] - d.lon[0])*2), dtype= int)%720
            pts=(list(PLlist[1]), list(Tlatraw), list(Tlonraw)) #the raw data points in d format
         
            gradtheta_em= []
            for i, t in enumerate(PLlist[1]):  #loop through every PL point            
                londistgradtheta_e= int(radgradtheta_e//np.cos(np.deg2rad(PLlist[3][i])))           
                mask= radMask2D((len(d.lat), len(d.lon)), (pts[1][i], pts[2][i]), radiusx= radgradtheta_e, radiusy= londistgradtheta_e) #lat lon mask
    
                """calculate the gradient in the equi potent temp"""
                Theta_e850= EquiPotTemp(d.T850[t], d.sHum850[t], plev= 850)
                latd= 55.
                lond= np.tile(latd*np.cos(np.deg2rad(d.lat)), (720,1)).T         
                g= np.gradient(Theta_e850, latd, lond) 
                absg= np.sqrt(g[0]**2+ g[1]**2) *100 #change in Theta_e /100km
    
                gradtheta_em += [np.max(absg[mask])]
    
    
            PLlist=  np.append(PLlist, np.array(gradtheta_em).reshape(1, -1), axis= 0) #append the stability to the list
            PLlist= PLlist[:, np.array(gradtheta_em) < gradtheta_elim]
            
            gtnr += len(remove_dublicate(PLlist[0]))  


        """exclude theta e,700 stability cyclones"""
        if theta_e700lim!= None:
            Tlatraw= np.array(np.round(-(PLlist[3]-d.lat[0])*2), dtype= int)
            Tlonraw= np.array(np.round((PLlist[2] - d.lon[0])*2), dtype= int)%720
            pts=(list(PLlist[1]), list(Tlatraw), list(Tlonraw)) #the raw data points in d format
            
            theta_e700m= []
            for i, t in enumerate(PLlist[1]):  #loop through every PL point            
                londist= int(radtheta_e700//np.cos(np.deg2rad(PLlist[3][i])))           
                mask= radMask2D((len(d.lat), len(d.lon)), (pts[1][i], pts[2][i]), radiusx= radtheta_e700, radiusy= londist) #lat lon mask
    
                """calculate potential temperature and equivalent potential temperature"""
                RCp = 2/7 #R/Cp
                Lc, Cp= 2501E3, 1006
                theta700= d.T[t,1][mask]*(1000/700)**RCp
                theta_e700= theta700* np.exp(d.sHum[t,1][mask]* Lc/(d.T[t,1][mask]* Cp))
                             
                thetaSST= d.SST[t][mask]*(1000/d.MSLP[t][mask])**RCp
                theta_eSST=thetaSST* np.exp(d.sHum[t,0][mask]* Lc/(d.SST[t][mask]* Cp))
                   
                theta_e700m += [np.max(theta_eSST- theta_e700)] #the max in radius
    
    
            PLlist=  np.append(PLlist, np.array(theta_e700m).reshape(1, -1), axis= 0) #append the stability to the list
            PLlist= PLlist[:, np.array(theta_e700m) > theta_e700lim]
            
            s700nr += len(remove_dublicate(PLlist[0]))   
        
        """exclude wind cyclones"""
        if windlim != None:
            d.imp_u10()
            U10= np.sqrt(d.u10**2+ d.v10**2)
            
            Tlatraw= np.array(np.round(-(PLlist[3]-d.lat[0])*2), dtype= int)
            Tlonraw= np.array(np.round((PLlist[2] - d.lon[0])*2), dtype= int)%720
            pts=(list(PLlist[1]), list(Tlatraw), list(Tlonraw)) #the raw data points in d format
            
            windm= []
            latdist= 5
            for i, t in enumerate(PLlist[1]):  #loop through every PL point            
                londist= int(latdist//np.cos(np.deg2rad(PLlist[3][i])))           
                mask= radMask2D((len(d.lat), len(d.lon)), (pts[1][i], pts[2][i]), radiusx= latdist, radiusy= londist) #lat lon mask
            
                windm += [np.max(U10[t][mask])]
            
            PLlist= PLlist[:, np.array(windm) > windlim]
            
            wnr += len(remove_dublicate(PLlist[0]))
            print("windspeed is not written in PLlist")
#        print('after stability removal: '+str(len(remove_dublicate(PLlist[0]))))


        """calculate shear, stability and UPVnorth for all points of the PL"""
        if hemi == 'NH': #to calculate the shear first some more variables have to be imported for the SH
            d.impvar('uPL')
            d.impvar('vPL')
            d.impvar('Geop')
        
        latdist_stab= radstab
            
        if len(PLlist > 0):
            for PLnr in remove_dublicate(PLlist[0]):
                shearinterm, stabinterm, UPVpolewinterm, pdiffinterm, gradtheta_einterm, theta_e700interm= [], [], [], [], [], []
                
                PLlistPLnr = T.PLlist[:, T.PLlist[0]== PLnr]
                Tlatraw= np.array(np.round(-(PLlistPLnr[3]-d.lat[0])*2), dtype= int)
                Tlonraw= np.array(np.round((PLlistPLnr[2] - d.lon[0])*2), dtype= int)%720
                pts=(list(PLlistPLnr[1]), list(Tlatraw), list(Tlonraw)) #the raw data points in d format
                for i, t in enumerate(PLlistPLnr[1]):  #loop through every PL point         
#                    print(t)

                    londist_stab= int(radstab//np.cos(np.deg2rad(PLlistPLnr[3][i])))           
                    mask_stab= radMask2D((len(d.lat), len(d.lon)), (pts[1][i], pts[2][i]), radiusx= radstab, radiusy= londist_stab) #lat lon mask
                    londist_pdiff= int(radpdiff//np.cos(np.deg2rad(PLlistPLnr[3][i])))           
                    mask_pdiff= radMask2D((len(d.lat), len(d.lon)), (pts[1][i], pts[2][i]), radiusx= radpdiff, radiusy= londist_pdiff) #lat lon mask
                    
                    if hemi == 'NH':
                        UPVpolewinterm += [np.max(UPV[pts[0][i], :pts[1][i] +1, pts[2][i] ])] #UPVnorth
                    elif hemi == 'SH':
                        UPVpolewinterm += [np.max(UPV[pts[0][i], pts[1][i]: , pts[2][i] ])] #UPVnorth
             
                    theta500= d.T[t,2][mask_stab]*(1000/500)**RCp  #stability
                    thetaSST= d.SST[t][mask_stab]*(1000/d.MSLP[t][mask_stab])**RCp  #stability
                    
                    stabinterm += [np.mean(thetaSST- theta500)]
                    pdiffinterm += [d.MSLP[t, pts[1][i], pts[2][i] ] - np.mean(d.MSLP[t][mask_pdiff]) ]

                    if hemi == 'NH':
                        londistgradtheta_e= int(radgradtheta_e//np.cos(np.deg2rad(PLlistPLnr[3][i])))           
                        mask_gradtheta_e= radMask2D((len(d.lat), len(d.lon)), (pts[1][i], pts[2][i]), radiusx= radgradtheta_e, radiusy= londistgradtheta_e) #lat lon mask
                        radshear= 2.5*2
                        londist_shear= int(radshear//np.cos(np.deg2rad(PLlistPLnr[3][i])))           
                        mask_shear= radMask2D((len(d.lat), len(d.lon)), (pts[1][i], pts[2][i]), radiusx= radshear, radiusy= londist_shear) #lat lon mask
                        londist= int(radtheta_e700//np.cos(np.deg2rad(PLlistPLnr[3][i])))           
                        masktheta_e= radMask2D((len(d.lat), len(d.lon)), (pts[1][i], pts[2][i]), radiusx= radtheta_e700, radiusy= londist) #lat lon mask
    
                        theta700= d.T[t,1][masktheta_e]*(1000/700)**RCp
                        theta_e700= theta700* np.exp(d.sHum[t,1][masktheta_e]* Lc/(d.T[t,1][masktheta_e]* Cp))                              
                        thetaSST= d.SST[t][masktheta_e]*(1000/d.MSLP[t][masktheta_e])**RCp
                        theta_eSST=thetaSST* np.exp(d.sHum[t,0][masktheta_e]* Lc/(d.SST[t][masktheta_e]* Cp))
                           
                        theta_e700interm += [np.max(theta_eSST- theta_e700)] #the max in radius
    
                        """gradtheta_e"""
                        Theta_e850= EquiPotTemp(d.T850[t], d.sHum850[t], plev= 850)
                        latd= 55.
                        lond= np.tile(latd*np.cos(np.deg2rad(d.lat)), (720,1)).T         
                        g= np.gradient(Theta_e850, latd, lond) 
                        absg= np.sqrt(g[0]**2+ g[1]**2) *100 #change in Theta_e /100km
            
                        gradtheta_einterm += [np.max(absg[mask_gradtheta_e])]                    
                        
                        """shear"""
                        uPLavg= np.mean(d.uPL[t], axis=0)[mask_shear]
                        vPLavg= np.mean(d.vPL[t], axis=0)[mask_shear]
                        thick=d.Geop[t,1]- d.Geop[t,0] #[mask]
                        u_T, v_T= CalculateThermalWind(thick, d.lat, d.lon)
                        alpha= np.rad2deg(np.arccos((u_T[mask_shear]*uPLavg + v_T[mask_shear]*vPLavg)/(np.sqrt(u_T[mask_shear]**2+v_T[mask_shear]**2)* np.sqrt(uPLavg**2+vPLavg**2))))
                        shearinterm += [np.mean(alpha)]


                
                UPVpolewlist += [UPVpolewinterm]
                stablist += [stabinterm]
                pdifflist += [pdiffinterm]
                if hemi == 'NH':
                    shearlist += [shearinterm]               
                    gradtheta_elist += [gradtheta_einterm]
                    theta_e700list += [theta_e700interm]
                
        """make the contribution of each PL to the PLlist"""
        if len(PLlist > 0):
            for i, PLnr in enumerate(remove_dublicate(PLlist[0])):
                intermlist= T.PLlist[:, T.PLlist[0]== PLnr]
                intermlist=np.append(intermlist, np.array(stablist[i]).reshape(1,-1), axis= 0) #append stability column              
                intermlist=np.append(intermlist, np.array(UPVpolewlist[i]).reshape(1,-1), axis= 0) #append UPVnorth column
                intermlist=np.append(intermlist, np.array(pdifflist[i]).reshape(1,-1), axis= 0) #append pdiff column                              

                if hemi == 'NH':
                    intermlist=np.append(intermlist, np.array(gradtheta_elist[i]).reshape(1,-1), axis= 0) #append gradtheta_e column                              
                    intermlist=np.append(intermlist, np.array(theta_e700list[i]).reshape(1,-1), axis= 0) #append theta_e700 column                                                   
                    intermlist=np.append(intermlist, np.array(shearlist[i]).reshape(1,-1), axis= 0) #append shear column              

                intermlist=np.append(intermlist, np.zeros((1,intermlist.shape[1])), axis= 0) #append the column that indicates PL point or not
                intermlist=np.append(np.ones((1, intermlist.shape[1]))*month, intermlist, axis= 0)
                intermlist=np.append(np.ones((1, intermlist.shape[1]))*year, intermlist, axis= 0)
                
                tones= PLlist[1, PLlist[0] ==PLnr] #makes the column of PL point or not
                for tone in tones:
                    intermlist[-1, intermlist[3]==tone] = 1
                
                cPLlist = np.append(cPLlist, intermlist, axis= 1)
        

    if hemi == 'NH':
        WritePLlistComp2S(cPLlist, year, name=versionname, rowname= ['year', 'month', 'TPLnr', 'time pnt', 'lon', 'lat', 'vort_filt', 'thetaSST-theta500', 'UPV polewards', 'low depth', 'max grad theta_e', 'theta_eSST - theta_e700', 'shear', 'PLpoint'])
    elif hemi == 'SH':
        WritePLlistComp2S(cPLlist, year, name=versionname, hemi='_SH_', rowname= ['year', 'month', 'TPLnr', 'time pnt', 'lon', 'lat', 'vort_filt', 'thetaSST-theta500', 'UPV polewards', 'low depth', 'PLpoint'])


print('cyclones:',cnr, 'vortexcl:', vnr, 'landexcl:', lnr, 'polarfrontexcl:', pnr, 'stabexcl:', snr, 'pdiffexcl:', pdnr, 'theta_e excl:', gtnr)


print("--- %s seconds ---" % (time.time() - start_time))