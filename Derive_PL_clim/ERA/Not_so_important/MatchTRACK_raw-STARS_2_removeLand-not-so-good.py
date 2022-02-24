#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 16:38:42 2017

@author: pst019
"""
import numpy as np

import sys

sys.path.insert(0, '/home/'+user+'/polar_low_code/Functions/')
from f_imp_ERA2 import *
from f_impLists import *
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
import datetime


timestep= 'threehours'
Sdurationlim= 6 #STARS duration in hours
Tdurationlim=2 #TRACK duration limit in 3hourly time steps

hitcount, Tcount, Scount, TPL, SPL, TPL_inclland= 0, 0, 0, 0, 0, 0

matchlist= []



"""compare PLlists"""
for year in [2005]:#np.arange(2002, 2012):
#    if year== 2011: monthslist= [1,2,3]
#    else: 
#        monthslist= [1, 2, 3, 4, 9, 10, 11, 12]
    monthslist= [1]
    
    for month in monthslist:
        print(year, month)

        T=TRACK_list_month(year, month, duration_limit= Tdurationlim, timestep= timestep)
        T.local(minlon= -45, maxlon= 65, minlat= 54, maxlat= 81) #for combined
 
        S= STARSList(syear= year, smonth= month, filename='STARS_TRACKS.csv', durationlim= Sdurationlim, timestep= timestep)
        S2= STARSList(syear= year, smonth= month, filename='STARS_TRACKS_south.csv', durationlim= Sdurationlim, timestep= timestep)
        S.PLnumber= np.concatenate((S.PLnumber, S2.PLnumber))
        S.lat= np.concatenate((S.lat, S2.lat))
        S.lon= np.concatenate((S.lon, S2.lon))
        S.tPL= np.concatenate((S.tPL, S2.tPL))

        d= data(['SST'], year, month)
        d.SST[d.SST<272]= np.nan #this sets the area with ice to land
        d.SST= ma.array(d.SST, mask= isnan(d.SST))

        latdist2= 4            

        TPLs = remove_dublicate(T.PLnumber) #ERA PLs
        TPL_inclland += len(TPLs) #ERA PLs

        landistr= []
        for i in range(len(T.lat)):
            if T.lat[i] < 85 and T.lat[i] > 30:
                londist2= int(latdist2//np.cos(np.deg2rad(T.lat[i])))
                mask= radMask2D(d.SST.shape[1:], (int(-(T.lat[i]-d.lat[0])*2), int((T.lon[i] - d.lon[0])*2)), radiusx= latdist2, radiusy= londist2)
                maskc2= [int(T.tPL[i]//2 ), mask]
                landistr += [len(np.where(d.SST[maskc2].mask == True)[0]) /len(d.SST[maskc2])]
            else:
                landistr += [1]
        
        landistr = np.array(landistr)
        T.lon = T.lon[landistr <.25]
        T.lat = T.lat[landistr <.25]
        T.PLnumber = T.PLnumber[landistr <.25]
        T.tPL = T.tPL[landistr <.25]


            
        TPLs = remove_dublicate(T.PLnumber) #ERA PLs
        TPL += len(TPLs) #ERA PLs
        SPLs = remove_dublicate(S.PLnumber)
        SPL += len(SPLs)
        
        
        for sPL in SPLs:
            tsPL= S.tPL[S.PLnumber== sPL]
            
            match= 0
            for n, ntsPL in enumerate(tsPL):
                matchnow= 0
                londist= T.lon[T.tPL== ntsPL] - S.lon[S.PLnumber== sPL][n]
                latdist= T.lat[T.tPL== ntsPL] - S.lat[S.PLnumber== sPL][n]
                dist= 110* np.sqrt(latdist**2 +  (np.cos(np.deg2rad(S.lat[S.PLnumber== sPL][n])) * londist)**2)
                if len(dist) != 0:
                    if np.min(dist) < 250:
                        print('match: STARS PLnumber ', sPL, 'TRACK PLnumber ', T.PLnumber[T.tPL== ntsPL][dist == np.min(dist)][0], 'dist; ', np.min(dist) )
                        match= 1
                        matchnow= 1
#                    else:
#                        print('STARS PLnumber ', sPL, 'to far away, dist: ', np.min(dist))
                    
                    matchlist += [[int(sPL), int(ntsPL), matchnow, int(np.min(dist)), int(T.PLnumber[T.tPL== ntsPL][dist == np.min(dist)][0]), year, month]]

                else:
                    matchlist += [[int(sPL), int(ntsPL), matchnow, 1E6, 1E6, year, month]] #in this case no cyclone is around
#                    print('STARS PLnumber ', sPL, 'no PL at the same time')

                #SPLnumber, StPL, hit, distance, TPLnumber of min
#                matchlist= +=[[int(sPL), int(ntsPL), match, int(np.min(dist)), int(T.PLnumber[T.tPL== ntsPL][dist == np.min(dist)][0]), year, month]]

            hitcount += match
#            print(hitcount)


print('TPL incl land: '+str(TPL_inclland), 'TPL: '+str(TPL), 'SPL: '+str(SPL), 'hitcount: '+str(hitcount))                                             

#
#"""write matchlist"""
#matchlist = np.array(matchlist, dtype=int)
#import csv
#harddisk= Mediadir+'PL/PLclim/MatchList/Matchlist2.csv'
#with open(harddisk, 'w') as csvfile:
#    writer = csv.writer(csvfile, delimiter="\t")
##    writer.writerow(['STARS PLnr', 'S tPL', 'match', 'distance', 'TRACK PLnr', 'year', 'month'])
#    writer.writerow(['SPLnr', 'tPL', 'match', 'dist', 'TPLnr', 'year', 'month'])
#    [writer.writerow(r) for r in matchlist]
#
#
#"""analyse matchlist"""
#matchlist = np.array(matchlist)
#TPLsmatch= []
#perfectcount = 0
#crapycount = 0
#moresemicount = 0
#lesssemicount = 0
#TRACKmatchSTARS= []
#
#SPLs= remove_dublicate(matchlist[:, 0]) #all STARS PLs in matchlist
#
#for sPL in SPLs:
#    PLnrmaxmatch= -1 #default for no match
#    matchlistnow= matchlist[matchlist[:,0]== sPL]
#    matches= remove_dublicate(matchlistnow[:, 2])
#    TPLsm= remove_dublicate(matchlistnow[:, 4]) #TRACK PLs that match
#
#    
#    if len(matches)== 1 and matches[0]== 1 and len(TPLsm)==1:
#        out= 'perfect match'
#        perfectcount += 1
#        PLnrmaxmatch= TPLsm[0]
#        quality= 0
#        
#    elif len(matches)== 1 and matches[0]== 0:
#        out= 'no match at all'
#        crapycount += 1
#        quality= 3
#        
#    else:
#        PLduration= len(matchlistnow[:, 4])
#        maxmatchduration= np.max(np.bincount(np.array(matchlistnow[:, 4], dtype=int))) #the time it matches longest
#        if 2*maxmatchduration > PLduration:
#            PLnrmaxmatch= np.argmax(np.bincount(np.array(matchlistnow[:, 4], dtype=int)))
#            out= 'more semimatch with '+str(PLnrmaxmatch)
#            moresemicount += 1
#            quality= 1
#        else:
#            quality= 2
#            out= 'less semimatch'
#            lesssemicount += 1
#
#        
##    else:
##        #if only some timesteps do not fit, but all the matching steps have the same cyclone
##        TPLsm= remove_dublicate(matchlistnow[matchlistnow[:,2]==1, 4])
##        
##        if len(TPLsm) ==1:
##            out= 'rework, semi match'
##            semicount += 1
#        
#    print(sPL, matches, TPLsm, quality, out)
#    if PLnrmaxmatch != -1:
#        TRACKmatchSTARS += [[sPL, PLnrmaxmatch, matchlistnow[0, -2], matchlistnow[0, -1], matchlistnow[matchlistnow[:,4]== PLnrmaxmatch][0, 1], matchlistnow[matchlistnow[:,4]== PLnrmaxmatch][-1, 1], quality]]
##    else:
##        TRACKmatchSTARS += [[sPL, PLnrmaxmatch, matchlistnow[0, -2], matchlistnow[-1, 1], matchlistnow[matchlistnow[:,4]== PLnrmaxmatch][0, 1], matchlistnow[matchlistnow[:,4]== PLnrmaxmatch][-1, 1], quality]]
#
#        
#""" write TRACKmatchSTARS
##this writes all the matches, tPLstart gives the starting point in time of the match"""
#import csv
#harddisk= Mediadir+'PL/PLclim/MatchList/TRACKmatchSTARS.csv' #_Sdurlim3.csv'
#with open(harddisk, 'w') as csvfile:
#    writer = csv.writer(csvfile, delimiter="\t")
#    writer.writerow(['SPLnr', 'TPLnr', 'year', 'month', 'tcombPLstart', 'tcombPLend',  'quality'])
#    [writer.writerow(r) for r in TRACKmatchSTARS]
#
#
#
#print('number STARS PLs', len(SPLs))
#print('perfect matches', perfectcount)
#print('more semi matches', moresemicount)
#print('less semi matches', lesssemicount)
#print('no match at all', crapycount)
##    TPLsmatch += remove_dublicate(matchlist[matchlist[:,0]== 1, 4])
#
#
