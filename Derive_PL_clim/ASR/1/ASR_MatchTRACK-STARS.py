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


timestep= 'threehours'
Sdurationlim= 3 #STARS duration in hours
Tdurationlim=1 #TRACK duration limit in 3hourly time steps

hitcount, Tcount, Scount, TPL, SPL= 0, 0, 0, 0, 0

matchlist= []



"""compare PLlists"""
#for year in np.arange(2002, 2012):
#    if year== 2011: monthslist= [1,2,3]
#    elif year== 2008: monthslist= [10, 11, 12]

    
for year in np.arange(2002, 2012):
    if year== 2011: monthslist= [1,2,3]
    else: 
        monthslist= [1, 2, 3, 4, 9, 10, 11, 12]
    
    for month in monthslist:
        print(year, month)

        T=TRACK_list_month2(year, month, duration_limit= Tdurationlim, timestep= timestep, model='ASR')
        T.local(minlon= -45, maxlon= 65, minlat= 54, maxlat= 81) #for combined
 
        S= STARSList(syear= year, smonth= month, filename='STARS_TRACKS.csv', durationlim= Sdurationlim, timestep= timestep)
        S2= STARSList(syear= year, smonth= month, filename='STARS_TRACKS_south.csv', durationlim= Sdurationlim, timestep= timestep)
        S.PLnumber= np.concatenate((S.PLnumber, S2.PLnumber))
        S.lat= np.concatenate((S.lat, S2.lat))
        S.lon= np.concatenate((S.lon, S2.lon))
        S.tPL= np.concatenate((S.tPL, S2.tPL))



        TPLs = remove_dublicate(T.PLlist[0]) #ERA PLs
        TPL += len(TPLs) #ERA PLs
        SPLs = remove_dublicate(S.PLnumber)
        SPL += len(SPLs)
        
        
        for sPL in SPLs:
            tsPL= S.tPL[S.PLnumber== sPL]
            
            match= 0
            for n, ntsPL in enumerate(tsPL):
                matchnow= 0
                londist= T.PLlist[2][T.PLlist[1]== ntsPL] - S.lon[S.PLnumber== sPL][n]
                latdist= T.PLlist[3][T.PLlist[1]== ntsPL] - S.lat[S.PLnumber== sPL][n]
                dist= 110* np.sqrt(latdist**2 +  (np.cos(np.deg2rad(S.lat[S.PLnumber== sPL][n])) * londist)**2)
                if len(dist) != 0:
                    if np.min(dist) < 250:
#                        print('match: STARS PLnumber ', sPL, 'TRACK PLnumber ', T.PLlist[0][T.PLlist[1]== ntsPL][dist == np.min(dist)][0], 'dist; ', np.min(dist) )
                        match= 1
                        matchnow= 1
#                    else:
#                        print('STARS PLnumber ', sPL, 'to far away, dist: ', np.min(dist))
                    
                    matchlist += [[int(sPL), int(ntsPL), matchnow, int(np.min(dist)), int(T.PLlist[0][T.PLlist[1]== ntsPL][dist == np.min(dist)][0]), year, month]]

                else:
                    matchlist += [[int(sPL), int(ntsPL), matchnow, -1, -1, year, month]]
#                    print('STARS PLnumber ', sPL, 'no PL at the same time')

                #SPLnumber, StPL, hit, distance, TPLnumber of min
#                matchlist= +=[[int(sPL), int(ntsPL), match, int(np.min(dist)), int(T.PLlist[0][T.PLlist[1]== ntsPL][dist == np.min(dist)][0]), year, month]]

            hitcount += match
#            print(hitcount)


print('TPL: '+str(TPL), 'SPL: '+str(SPL), 'hitcount: '+str(hitcount))                                             

#
"""write matchlist"""
import csv
harddisk= Mediadir+'PL/PLclim/MatchList/ASR_Matchlist.csv'
with open(harddisk, 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter="\t")
#    writer.writerow(['STARS PLnr', 'S tPL', 'match', 'distance', 'TRACK PLnr', 'year', 'month'])
    writer.writerow(['SPLnr', 'tPL', 'match', 'dist', 'TPLnr', 'year', 'month'])
    [writer.writerow(r) for r in matchlist]


"""analyse matchlist"""
matchlist = np.array(matchlist)
TPLsmatch= []
perfectcount = 0
crapycount = 0
moresemicount = 0
lesssemicount = 0
TRACKmatchSTARS= []

SPLs= remove_dublicate(matchlist[:, 0]) #all STARS PLs in matchlist

for sPL in SPLs:
    PLnrmaxmatch= -1 #default for no match
    matchlistnow= matchlist[matchlist[:,0]== sPL]
    matches= remove_dublicate(matchlistnow[:, 2])
    TPLsm= remove_dublicate(matchlistnow[:, 4]) #TRACK PLs that match

    
    if len(matches)== 1 and matches[0]== 1 and len(TPLsm)==1:
        out= 'perfect match'
        perfectcount += 1
        PLnrmaxmatch= TPLsm[0]
        quality= 0
        
    elif len(matches)== 1 and matches[0]== 0:
        out= 'no match at all'
        crapycount += 1
        quality= 3
        
    else:
        PLduration= len(matchlistnow[:, 4])
        maxmatchduration= np.max(np.bincount(matchlistnow[:, 4])) #the time it matches longest
        if 2*maxmatchduration > PLduration:
            PLnrmaxmatch= np.argmax(np.bincount(matchlistnow[:, 4]))
            out= 'more semimatch with '+str(PLnrmaxmatch)
            moresemicount += 1
            quality= 1
        else:
            quality= 2
            out= 'less semimatch'
            lesssemicount += 1

        
#    else:
#        #if only some timesteps do not fit, but all the matching steps have the same cyclone
#        TPLsm= remove_dublicate(matchlistnow[matchlistnow[:,2]==1, 4])
#        
#        if len(TPLsm) ==1:
#            out= 'rework, semi match'
#            semicount += 1
        
    print(sPL, matches, TPLsm, quality, out)
    if PLnrmaxmatch != -1:
        TRACKmatchSTARS += [[sPL, PLnrmaxmatch, matchlistnow[0, -2], matchlistnow[0, -1], matchlistnow[matchlistnow[:,4]== PLnrmaxmatch][0, 1], matchlistnow[matchlistnow[:,4]== PLnrmaxmatch][-1, 1], quality]]
#    else:
#        TRACKmatchSTARS += [[sPL, PLnrmaxmatch, matchlistnow[0, -2], matchlistnow[-1, 1], matchlistnow[matchlistnow[:,4]== PLnrmaxmatch][0, 1], matchlistnow[matchlistnow[:,4]== PLnrmaxmatch][-1, 1], quality]]

        
""" write TRACKmatchSTARS
#this writes all the matches, tPLstart gives the starting point in time of the match"""
import csv
harddisk= Mediadir+'PL/PLclim/MatchList/ASR_TRACKmatchSTARS.csv' #_Sdurlim3.csv'
with open(harddisk, 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(['SPLnr', 'TPLnr', 'year', 'month', 'tcombPLstart', 'tcombPLend',  'quality'])
    [writer.writerow(r) for r in TRACKmatchSTARS]



print('number STARS PLs', len(SPLs))
print('perfect matches', perfectcount)
print('more semi matches', moresemicount)
print('less semi matches', lesssemicount)
print('no match at all', crapycount)
#    TPLsmatch += remove_dublicate(matchlist[matchlist[:,0]== 1, 4])


