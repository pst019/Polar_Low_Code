#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 16:29:52 2019

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/PatsOrange/'

homedir=Mediadir+'home/'


import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import pandas as pd
import numpy as np
import glob

import datetime as datetime
from f_useful import remove_dublicate



track_dir= Mediadir+"/data/ERA5_Clim/track_lists/fram/"
write= False

version= "fram_run1"
version180= 'fram_run180E'
ending= 'SH' #'atl-pac'
ending_name= 'SH180' #for the output

durlim= 6

start_year= 1979
end_year= 1979
i= 0

for year in range(start_year, end_year +1):
    print(year)
    for month in range(1, 2):
        i += 1
            
        tracks0= pd.read_csv(track_dir + version+'_'+str(year)+'_'+str(month).zfill(2)+'_'+ending +".csv") #.set_index(['ID', 'step'])
        tracks0['time']= pd.to_datetime(tracks0['time'])
        tracks0.ID= (tracks0.ID + 3E10).astype(int)
        tracks0 = tracks0.set_index('ID')
        
        tracks0.lon %= 360 #to shift the longitude from to 0 - 360E

        tracks0['180merged']= 0 #to label the tracks that are merged around 180 and in the end exclude the ones that are not merged
    
        """merge around 180E"""
        tracks180= pd.read_csv(track_dir + version180+'_'+str(year)+'_'+str(month).zfill(2)+'_'+ending +".csv") #.set_index(['ID', 'step'])
        tracks180['time']= pd.to_datetime(tracks180['time'])
        tracks180 = tracks180.set_index('ID')

        tracks180_red= tracks180[(tracks180.lon > 175) & (tracks180.lon < 185)]

        # tracks180_red= tracks180_red[(tracks180_red['time'] >= '1979-01-07') & (tracks180_red['time'] < '1979-01-09 00:00:00')]

        for track_ID in remove_dublicate(tracks180_red.index):
            """loop through every 180E track"""
            track180= tracks180_red.loc[track_ID]

            print('track 180: \n', track180)

            if type(track180) == pd.core.series.Series:  continue #if track180 has only one time step
        
            
            matched_track_ID_list= [] #find all matches in global tracks
            for ts in range(len(track180)): 
                matched_track_ID= tracks0[(tracks0.time== track180.time.iloc[ts]) & (tracks0.lon== track180.lon.iloc[ts])
                                      & (tracks0.lat== track180.lat.iloc[ts])]                
                if len(matched_track_ID) != 0:
                    matched_track_ID_list += [matched_track_ID.index.values[0] ]
                else:
                    continue           

            if len(matched_track_ID_list)== 0:
                print('No match')
                continue

            matched_track_ID_list= remove_dublicate(matched_track_ID_list)

            matched_track = tracks0.loc[matched_track_ID_list]

            print('matched track: \n', matched_track)
            
            new_track= matched_track[(matched_track.lon <= 175) | (matched_track.lon >= 185)]
            new_track= new_track.append(track180)
            new_track['ID']= matched_track_ID_list[0]
            new_track= new_track.set_index('ID')
            
            new_track= new_track.sort_values(by= 'time')
            new_track['step']= np.arange(len(new_track))

            # print('new_track: \n', new_track)

            tracks0= tracks0.drop(index= matched_track_ID_list).append(new_track)
            tracks0.loc[matched_track_ID_list[0], '180merged']= 1
 
            print('new_track in list: \n', tracks0.loc[matched_track_ID_list[0]])
        
 
        tracks0= tracks0.reset_index()    
        tracks0.area= tracks0.area.astype(int)   
    
        if i == 1:
            alltracks= tracks0
        else:
            """merge in time"""
            mergetime= datetime.datetime(year, month, 1)
            mergetime_tracks0 = tracks0[tracks0['time'] == mergetime]
            mergetime_alltracks= alltracks[alltracks['time']== mergetime]
            
            for it in range(len(mergetime_tracks0)):
                track0_now= mergetime_tracks0.iloc[it]
                
                merge_track= mergetime_alltracks[(mergetime_alltracks['lat'] == track0_now.lat) & (mergetime_alltracks['lon'] == track0_now.lon )]
                
                # print(len(merge_track))
                if len(merge_track)== 1:
                    

                    #have to remove the first occurance, not to have it double                    
                    tracks0.drop( tracks0[(tracks0.ID == track0_now.ID) & (tracks0.step == 0)].index, inplace =True)
                  
                    tracks0.loc[(tracks0.ID == track0_now.ID), 'step'] += merge_track.step.values[0]
                    tracks0.loc[(tracks0.ID == track0_now.ID), 'ID'] = merge_track.ID.values[0]

                    # tracks0['ID'][(tracks0.ID == track0_now.ID)] = merge_track.ID
            
            
            alltracks= pd.concat([alltracks, tracks0])



#exclude the tracks that are not merged around 180, the band is smaller since some tracks that were not merged have a short extension into the band that is not covered by the 180tracks
alltracks= alltracks[(alltracks.lon < 178) | (alltracks.lon > 182) | (alltracks['180merged'] == 1)]    


"""create the duration"""
# alltracks= alltracks.reset_index()
if durlim:
    alltracks= alltracks.join( (alltracks.groupby('ID').last()['step']).rename('Duration'), on='ID')
    alltracks= alltracks[alltracks['Duration'] > durlim]

    alltracks= alltracks.drop(['Duration', 'slp', 'vortex_type', '180merged'], axis= 1)

alltracks= alltracks.set_index(['ID', 'step'])
alltracks= alltracks.sort_index()


csv_name= 'merged_tracks_'+version+'_'+ending_name+f'_{start_year}-{end_year}'

if durlim: csv_name += 'durlim'+str(durlim)

if write:
    alltracks.to_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv')




