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

sys.path.insert(0, '/home/'+user+'/codeDerive_PL_clim/ERA/')
from f_imp_ERA2 import *
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

Tdurationlim=2 #TRACK duration limit in 3hourly time steps - higher less PLs
vortlim= 5.4 #the higher the less PLs
#theta_establim= np.float64(-12.4) #the higher the less PLs
thetastablim= np.float64(-11.2) #the higher the less PLs

#stablim= 43
#windlim= 14 #the higher the less PLs
polarfrontlim= 30.7 #the higher the more PLs
#jetlim=40


cnr, vnr, lnr, pnr, jnr, snr, wnr= 0,0,0,0,0,0,0

for year in np.arange(1981, 2017):
    print(year)
    cPLlist= np.zeros((10,0)) #complete PLlist (year, month, TrackPLnr, tPL, lon, lat, PLpoint(yes= 1, no= 0))

    
    for month in np.arange(1, 13): #[1, 2, 3, 4, 5, 9, 10, 11, 12]:
        shearlist, stablist, UPVsouthlist= [], [], []
    
        """import data"""
#        d= data(['SST', 'T', 'uPV', 'vPV', 'sHum', 'MSLP'], year, month)
        d= data(['SST', 'T', 'uPV', 'vPV', 'MSLP'], year, month, hemisp='SH')
        d.SST[d.SST<272]= np.nan #this sets the area with ice to land
        d.SST= ma.array(d.SST, mask= isnan(d.SST))        
        
#        d= data(['SST', 'T', 'uPV', 'vPV'], year, month)

        UPV= np.sqrt(d.uPV**2+ d.vPV**2)
    
        T=TRACK_list_month2(year, month, Tdurationlim, hemisp= 'SH', approved= False)
        T.local(minlon= d.lon[0], maxlon= d.lon[-1], minlat= d.lat[-1], maxlat=d.lat[0])
        #PLlist= [PLnumber, tPL, lon, lat, vort]

        cnr += len(remove_dublicate(T.PLlist[0]))
#        print('number of cyclones: '+str(len(remove_dublicate(T.PLlist[0]))))
            
        
        """extract PLs from cyclones """        
        PLlist= T.PLlist[:, T.PLlist[4] >= vortlim]

        vnr += len(remove_dublicate(PLlist[0]))
#        print('after vorticity removal: '+str(len(remove_dublicate(PLlist[0]))))
                 

        """exclude land cyclones (new version)"""
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
        Tlatraw= np.array(np.round(-(PLlist[3]-d.lat[0])*2), dtype= int)
        Tlonraw= np.array(np.round((PLlist[2] - d.lon[0])*2), dtype= int)%720
#        pts=(list(PLlist[1]), list(Tlatraw), list(Tlonraw)) #the raw data points in d format

        UPVsouthm = np.array([np.max(UPV[PLlist[1][i], int(Tlatraw[i]):, int(Tlonraw[i])]) for i in range(len(PLlist[0]))])
        
#        PLlist=  np.append(PLlist, UPVnorthm.reshape(1, -1), axis= 0) #append UPVnorth to the list
        PLlist= PLlist[:, UPVsouthm< polarfrontlim]
        
        pnr +=len(remove_dublicate(PLlist[0]))       



        """exclude theta stability cyclones"""
        Tlatraw= np.array(np.round(-(PLlist[3]-d.lat[0])*2), dtype= int)
        Tlonraw= np.array(np.round((PLlist[2] - d.lon[0])*2), dtype= int)%720
        pts=(list(PLlist[1]), list(Tlatraw), list(Tlonraw)) #the raw data points in d format
        
        theta_500m= []
        latdist= 2
        for i, t in enumerate(PLlist[1]):  #loop through every PL point            
            londist= int(latdist//np.cos(np.deg2rad(PLlist[3][i])))           
            mask= radMask2D((len(d.lat), len(d.lon)), (pts[1][i], pts[2][i]), radiusx= latdist, radiusy= londist) #lat lon mask

            """calculate potential temperature and equivalent potential temperature"""
            RCp = 2/7 #R/Cp
            theta500= d.T[t,2][mask]*(1000/500)**RCp
            thetaSST= d.SST[t][mask]*(1000/d.MSLP[t][mask])**RCp
                         
            theta_500m += [np.mean(thetaSST- theta500)]


        PLlist=  np.append(PLlist, np.array(theta_500m).reshape(1, -1), axis= 0) #append the stability to the list
        PLlist= PLlist[:, np.array(theta_500m) > thetastablim]
        
        snr += len(remove_dublicate(PLlist[0]))        



        """calculate shear, stability and UPVnorth for all points of the PL"""
#        d.impvar('uPL')
#        d.impvar('vPL')
#        d.impvar('Geop')
#        
        latdist= 2        
        if len(PLlist > 0):
            for PLnr in remove_dublicate(PLlist[0]):
                shearinterm, stabinterm, UPVsouthinterm = [], [], []
                
                PLlistPLnr = T.PLlist[:, T.PLlist[0]== PLnr]
                Tlatraw= np.array(np.round(-(PLlistPLnr[3]-d.lat[0])*2), dtype= int)
                Tlonraw= np.array(np.round((PLlistPLnr[2] - d.lon[0])*2), dtype= int)%720
                pts=(list(PLlistPLnr[1]), list(Tlatraw), list(Tlonraw)) #the raw data points in d format
                for i, t in enumerate(PLlistPLnr[1]):  #loop through every PL point         
#                    print(t)

                    londist= int(latdist//np.cos(np.deg2rad(PLlistPLnr[3][i])))           
                    mask= radMask2D((len(d.lat), len(d.lon)), (pts[1][i], pts[2][i]), radiusx= latdist, radiusy= londist) #lat lon mask
                    
                                    
                    UPVsouthinterm += [np.max(UPV[pts[0][i], pts[1][i]: , pts[2][i] ])] #UPVnorth
                    
                    theta500= d.T[t,2][mask]*(1000/500)**RCp  #stability
                    thetaSST= d.SST[t][mask]*(1000/d.MSLP[t][mask])**RCp  #stability
                    
#                    theta_e500= theta500* np.exp(d.sHum[t,2][mask]* Lc/(d.T[t,2][mask]* Cp))
                    stabinterm += [np.mean(thetaSST- theta500)]

                    """shear"""
#                    uPLavg= np.mean(d.uPL[t], axis=0)[mask]
#                    vPLavg= np.mean(d.vPL[t], axis=0)[mask]
#                    thick=d.Geop[t,1]- d.Geop[t,0] #[mask]
#                    u_T, v_T= CalculateThermalWind(thick, d.lat, d.lon)
#                    alpha= np.rad2deg(np.arccos((u_T[mask]*uPLavg + v_T[mask]*vPLavg)/(np.sqrt(u_T[mask]**2+v_T[mask]**2)* np.sqrt(uPLavg**2+vPLavg**2))))
#                    shearinterm += [np.mean(alpha)]


                
#                shearlist += [shearinterm]
                UPVsouthlist += [UPVsouthinterm]
                stablist += [stabinterm]


        """make the contribution of each PL to the PLlist"""
        if len(PLlist > 0):
            for i, PLnr in enumerate(remove_dublicate(PLlist[0])):
                intermlist= T.PLlist[:, T.PLlist[0]== PLnr]
                intermlist=np.append(intermlist, np.array(stablist[i]).reshape(1,-1), axis= 0) #append stability column              
                intermlist=np.append(intermlist, np.array(UPVsouthlist[i]).reshape(1,-1), axis= 0) #append UPVnorth column                              
#                intermlist=np.append(intermlist, np.array(shearlist[i]).reshape(1,-1), axis= 0) #append shear column              
                intermlist=np.append(intermlist, np.zeros((1,intermlist.shape[1])), axis= 0) #append the column that indicates PL point or not
                intermlist=np.append(np.ones((1, intermlist.shape[1]))*month, intermlist, axis= 0)
                intermlist=np.append(np.ones((1, intermlist.shape[1]))*year, intermlist, axis= 0)
                
                tones= PLlist[1, PLlist[0] ==PLnr] #makes the column of PL point or not
                for tone in tones:
                    intermlist[-1, intermlist[3]==tone] = 1
                
                cPLlist = np.append(cPLlist, intermlist, axis= 1)
        


    WritePLlistComp2S(cPLlist, year, hemi='_SH_')

print('cyclones:',cnr, 'vortexcl:', vnr, 'landexcl:', lnr, 'polarfrontexcl:', pnr, 'stabexcl:', snr)


print("--- %s seconds ---" % (time.time() - start_time))