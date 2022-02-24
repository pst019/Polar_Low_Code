#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 16:29:52 2019

@author: pst019
only in python 3.6
"""

import os
user = os.getcwd().split('/')[2]

if user=='pst019':
#    Mediadir= '/media/'+user+'/PatsOrange/'
    Mediadir= '/media/'+user+'/Backup1/'

else:
    Mediadir= '/run/media/pst019/PatsOrange/'

homedir=Mediadir+'home/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')


import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt

from f_useful import *
from f_ERA5_clim import *



write= False


durlim= 6


matchdist= 150
time_diff= False
time_diff= 24 #does not include the matched track if the matched track starts more than time_diff earlier or ends more than time_diff later

more_info= True
post_process= True

terrain_dist, terrain_thresh= 2, 50
#terrain_thresh= False

#split, exceed_large_to_small= True, 40
split= False

# PL_list_name, minObs= 'Rojo', 4
PL_list_name, minObs = 'Noer', 4
# PL_list_name, minObs = 'Smirnova', 4
# PL_list_name, minObs = 'Golubkin', 4

PL_list_name, minObs, time_diff = 'Yanase', 1, False #it has only one location per PL
# PL_list_name, minObs = 'Verezemskaya', 4 #it has only one location per PL


if PL_list_name == 'Verezemskaya':
    version= "" #"fram_run1_"
    ending= 'SH180'
    start_year= 1979
    end_year= 2020
    
    # time_start= False #in this case the first and last index of alltracks will be used
    time_start= '2004-1-1'
    time_end= '2005-1-1'
    
else:
    # version= "fram_run3_"
    version= ""
    ending= 'atl-pac'
    start_year= 1979
    end_year= 2020

    time_start= False #in this case the first and last index of alltracks will be used
    # time_start= '2002-1-1'
    # time_end= '2003-1-1'




"""get all tracks"""
csv_name= 'merged_tracks_'+version+ending+f'_{start_year}-{end_year}'

if durlim: csv_name += 'durlim'+str(durlim)
if more_info: csv_name += '_moreinfo'
if post_process: csv_name += '_post-process'

if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
if split: csv_name += f'_split{exceed_large_to_small}'


alltracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks2/"+csv_name+'.csv')

alltracks['time']= pd.to_datetime(alltracks['time'])

alltracks_IDs= np.array(remove_dublicate(alltracks.ID) ) #to find out if there are splits among the matched PLs
matchsplit_count = 0

alltracks= alltracks.set_index(['ID', 'step'])


if time_start == False:
    time_start = alltracks.time.iloc[0]
    time_end = alltracks.time.iloc[-1]


"""get PL list"""
if PL_list_name == 'Rojo':

    Rojo = import_Rojo_STARS(Mediadir, "data/PLclim/Rojo-etal_2019_fixes.csv")

    PLlist= Rojo


if PL_list_name == 'Noer':
    STARS= import_Gunnar_STARS(Mediadir, "data/PLclim/STARS/")

    STARS= S_Obs_nr(STARS)

    PLlist= STARS


if PL_list_name == 'Smirnova':

    Smir= pd.read_excel(Mediadir +'data/PLclim/Smirnova_PL_1995-2009.xls', 'swaths')
    Smir['ID']= np.nan

    ID_counter= 1
    
    for i in range(len(Smir)): #make the ID
        if np.isnan(Smir.iloc[i].Year):
            ID_counter += 1
        else:
            Smir['ID'].iloc[i]= ID_counter

    Smir = Smir.dropna(axis= 0, how= 'all') #remove empty rows

    #make the time
    Smir =Smir.assign(time= [f"{str(Year)}-{str(Month)}-{str(Day)} {str(Time)}" for Year, Month, Day, Time in zip(Smir['Year'].astype(int), Smir['Month'].astype(int), Smir['Day'].astype(int), Smir['Time'])])
    Smir['time']= pd.to_datetime(Smir['time'])
    Smir= Smir.drop(['Year', 'Month', 'Day', 'Time', 'Max Wind, m/s'], axis= 1)

    Smir= Smir.rename(columns= {'Latitude':'lat', 'Longitude': 'lon'})
    Smir['ID']= Smir['ID'].astype(int)
    Smir= S_Obs_nr(Smir)
    
    PLlist= Smir

if PL_list_name == 'Golubkin':

    Gol= pd.read_excel(Mediadir +'data/PLclim/Golubkin_PL_2015-2017.xls')
    # Gol['ID']= np.nan

    # ID_counter= 1
    
    # for i in range(len(Smir)): #make the ID
    #     if np.isnan(Smir.iloc[i].Year):
    #         ID_counter += 1
    #     else:
    #         Smir['ID'].iloc[i]= ID_counter
    Gol= Gol.where(Gol['Num'] != 0) #set the rows where num= 0 to nan
    
    Gol = Gol.dropna(axis= 0, how= 'all') #remove empty rows

    #make the time
    Gol =Gol.assign(time= [f"{str(Year)}-{str(Month)}-{str(Day)} {str(Time)}" for Year, Month, Day, Time in zip(Gol['Year'].astype(int), Gol['Month'].astype(int), Gol['Day'].astype(int), Gol['HH:MM'])]) 
    Gol['time']= pd.to_datetime(Gol['time'])
    
    Gol= Gol.drop(['Year', 'Month', 'Day', 'HH:MM', 'Wind, m/s', 'Diam, km', 'Stage', 'Unnamed: 11', 'Sat', 'Group'], axis= 1)

    Gol= Gol.rename(columns= {'Num': 'ID', 'Lat':'lat', 'Lon': 'lon'})
    Gol['ID']= Gol['ID'].astype(int)
    
    Gol= S_Obs_nr(Gol)
    
    PLlist= Gol

if PL_list_name == 'Yanase':
    Yan= pd.read_csv(Mediadir +'data/PLclim/Yanase_PLs.csv')
    
    Yan['time']= pd.to_datetime(Yan['time'], format= "%H%M UTC %d %b %Y")
    Yan['ID']= np.arange(1, len(Yan)+1)
    Yan['Obs']= 1
    
    PLlist= Yan
    
if PL_list_name == 'Verezemskaya':

    Ver= pd.read_csv(Mediadir +'data/PLclim/Verezemskaya_Antarctic_mesocyclone_tracks_adapt.csv')  #removed two rows that did not appear to belong to the same track  
    Ver.columns=['env_type', 'cloud_type', 'ID', 'time', 'lon', 'lat', 'cloud_diameter']
    Ver= Ver.drop(columns=['env_type', 'cloud_type', 'cloud_diameter'])
    Ver['time']= pd.to_datetime(Ver['time'], format= "%Y%m%d%H")
    Ver= S_Obs_nr(Ver)
    
    PLlist= Ver
    
"""common preparation"""

PLlist['time']= PLlist['time'].dt.round('H')
PLlist= PLlist.set_index(['ID', 'Obs'])

PLlist= PLlist[(PLlist['time'] > time_start) & (PLlist['time'] < time_end)]

PLlist= PLlist.join( (PLlist.reset_index().groupby('ID').last()['Obs']).rename('nObs'), on='ID')

PLlist= PLlist[PLlist['nObs'] >= minObs]

    


PLlist_indexes= remove_dublicate(PLlist.index.get_level_values(0)) #list of all PLs in PLlist


match_count= 0
match_track= []
match_track_PLlist= []
total= len(PLlist_indexes)

start_earlier= []
end_later= []

from collections import Counter

pd.options.mode.chained_assignment = None  # default='warn'

for index in PLlist_indexes: #['175']: #PLlist_indexes:
    
    PLlist_PL= PLlist.loc[index]
    
    
    matches= []
    print('\n PLlist nr ', index)
    for Obs in range(1, len(PLlist_PL)+1):
        PLlist_Obs= PLlist_PL.loc[Obs]
    
        track_now= alltracks[alltracks['time']== PLlist_Obs['time']]
    
        # track_now['distance'] =distance( (PLlist_Obs['lat'], PLlist_Obs['lon']) , (track_now['lat'], track_now['lon']) )
        track_now['distance'] = track_now.loc[track_now.index]['distance'] =distance( (PLlist_Obs['lat'], PLlist_Obs['lon']) , (track_now['lat'], track_now['lon']) )
        closest_track= track_now['distance'].idxmin()
        
        # print(Obs, track_now.loc[closest_track].name, np.round(track_now.loc[closest_track]['distance']))
        
        matches += list(track_now[track_now['distance']< matchdist].index.get_level_values('ID'))
    
    c= Counter(matches)
    print(c)
    npc= np.array(list(c.items()) )
    
    if len(npc)== 0: print('no match')
    else:
        maxcount= np.max(npc[:,1])
        nObs= PLlist_PL['nObs'].values[0]
        print(maxcount, ' of ', nObs, 'detected assciated to track: ', npc[np.argmax(npc[:,1]), 0])
        if maxcount/nObs >=0.5:
            matched_track_index= npc[np.argmax(npc[:,1]), 0]           
            
            start_earlier_this= (PLlist_PL.time.iloc[0] - alltracks.loc[matched_track_index].time.iloc[0])/ np.timedelta64(1, 'h')
            print('Matched PL starts earlier by:', start_earlier_this)
            end_later_this= (alltracks.loc[matched_track_index].time.iloc[-1] - PLlist_PL.time.iloc[-1])/ np.timedelta64(1, 'h')

            print('Matched PL end later by:', end_later_this)
            
            start_earlier += [start_earlier_this ]            
            end_later += [ end_later_this ]            
            if time_diff == False:
                match_count +=1
                match_track += [matched_track_index ] #add the track with the most counts
                match_track_PLlist += [index]
            else:
                if start_earlier_this <= time_diff and end_later_this <= time_diff:
                    match_count +=1
                    match_track += [matched_track_index ] #add the track with the most counts
                    match_track_PLlist += [index]

            """to find if the matched track was splitted"""
            if len(alltracks_IDs[(alltracks_IDs >= matched_track_index -3) & (alltracks_IDs <= matched_track_index +3)]) > 1:
                print('Has a split')
                matchsplit_count += 1



"""plot histograms of how much earlier and later the matched tracks occur"""
print('In total ', len(start_earlier), ' of ', total, ' were detected ', len(start_earlier)/total)


fignr= 1
plt.figure(fignr)
fignr+=1
plt.clf() 

plt.hist(start_earlier, density= True, alpha= 0.6, bins= np.arange(np.min(start_earlier)//3 *3 , np.max(start_earlier)+3, 3))

plt.xlabel('Matched track start earlier than PL [h]')
plt.ylabel('Density')
plt.plot(time_diff, 0, '^')


plt.figure(fignr)
fignr+=1
plt.clf() 

plt.hist(end_later, density= True, alpha= 0.6, bins= np.arange(np.min(end_later)//3 *3, np.max(end_later)+3, 3))

plt.xlabel('Matched track ends later than PL [h]')
plt.ylabel('Density')
plt.plot(time_diff, 0, '^')


if time_diff != False:
    start_earlier= np.array(start_earlier)
    # print(len(start_earlier[start_earlier > time_diff]), 'of', len(start_earlier), 'matched tracks start more than', time_diff, 'earlier')
    
    end_later= np.array(end_later)
    # print(len(end_later[end_later > time_diff]), 'of', len(end_later), 'matched tracks finish more than', time_diff, 'later')
    
    
    #print( len(start_earlier[(start_earlier <= time_diff) & (end_later <= time_diff)]), 'of', len(start_earlier), 'remain if these are excluded (for some tracks both is the case)' )
    print('> 24h excludes:', len(start_earlier)- len(start_earlier[(start_earlier <= time_diff) & (end_later <= time_diff)]), 'of', len(start_earlier))
    
    

"""write to list"""
print('After this ', match_count, ' of ', total, ' are detected ', match_count/total)

match_track_dict= {track_PL: track for track, track_PL in zip(match_track, match_track_PLlist)}
match_track= remove_dublicate(match_track)

print(match_count - len(match_track), ' of ', match_count, 'are double, such that ', len(match_track), ' remain')

tracks_matchPLlist= alltracks.loc[match_track]


tracks_matchPLlist["PLlist_match"] = ''
tracks_matchPLlist["PLlist_match_2"] = ''
tracks_matchPLlist["PLlist_match_3"] = ''

print('\n', matchsplit_count, 'matched PLs were splitted')


if write:
    for PL in match_track_dict:
        if tracks_matchPLlist['PLlist_match'].loc[ match_track_dict[PL] ].iloc[0] == '': #in case it is empty
            tracks_matchPLlist['PLlist_match'].loc[ match_track_dict[PL] ] = PL 
        else:
            if tracks_matchPLlist['PLlist_match_2'].loc[ match_track_dict[PL] ].iloc[0] == '': #in case it is empty
                tracks_matchPLlist['PLlist_match_2'].loc[ match_track_dict[PL] ] = PL 
            else:
                if tracks_matchPLlist['PLlist_match_3'].loc[ match_track_dict[PL] ].iloc[0] == '': #in case it is empty
                    tracks_matchPLlist['PLlist_match_3'].loc[ match_track_dict[PL] ] = PL 
                else:
                    print('There is a quadruple match')
    
    outfile= Mediadir+"/data/ERA5_Clim/track_lists/matchPLlist2/matchPLlist_"+ PL_list_name +f"_dist{matchdist}"
    if time_diff != False: outfile += f'_timediff{time_diff}'
    outfile += f'_{csv_name}.csv'
    print('\n write:', outfile)
    tracks_matchPLlist.to_csv(outfile)
