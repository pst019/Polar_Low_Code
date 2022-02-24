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


import pandas as pd
import numpy as np
import glob

import datetime as datetime
# from dateutil.relativedelta import relativedelta



track_dir= Mediadir+"/data/ERA5_Clim/track_lists/fram/"
write= False

version= "fram_run1"
ending_name= 'SH' #'atl-pac'

durlim= 6

start_year= 1979
end_year= 1979
i= 0

for year in range(start_year, end_year +1):
    print(year)
    for month in range(1, 2):
        for ending in ['SH']: #['atl', 'pac']:
            i += 1
        
            run= version+'_'+str(year)+'_'+str(month).zfill(2)+'_'+ending
            
            tracks0= pd.read_csv(track_dir + run +".csv") #.set_index(['ID', 'step'])
            tracks0['time']= pd.to_datetime(tracks0['time'])
            
            if ending== 'atl': tracks0.ID= (tracks0.ID + 1E10).astype(int)
            elif ending== 'pac': tracks0.ID= (tracks0.ID + 2E10).astype(int)
            elif ending== 'SH': tracks0.ID= (tracks0.ID + 3E10).astype(int)

            tracks0.area= tracks0.area.astype(int)            
            
            if i == 1:
                alltracks= tracks0
            else:
                
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




            


"""create the duration"""
# alltracks= alltracks.reset_index()
if durlim:
    alltracks= alltracks.join( (alltracks.groupby('ID').last()['step']).rename('Duration'), on='ID')
    alltracks= alltracks[alltracks['Duration'] > durlim]

    alltracks= alltracks.drop(['Duration', 'slp', 'vortex_type'], axis= 1)

alltracks= alltracks.set_index(['ID', 'step'])
alltracks= alltracks.sort_index()


csv_name= 'merged_tracks_'+version+'_'+ending_name+f'_{start_year}-{end_year}'

if durlim: csv_name += 'durlim'+str(durlim)

if write:
    alltracks.to_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv')




