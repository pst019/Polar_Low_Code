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
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/PatsOrange/'

homedir=Mediadir+'home/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/test_Denis/mc_era5-master/code') #to get obs_tracks_api


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob

from f_useful import *



test='2013_03'
time_start= '2013-3-1'
time_end= '2013-4-1'

test= 'run3_2011_01_atl'
time_start= '2011-1-1'
time_end= '2011-2-1'


def import_STARS(Mediadir= '/media/pst019/1692A00D929FEF8B/', file= "PL/STARS/Rojo-etal_2019.csv",
                 droplist= ['Optional_diameter', 'Comment', 'Comment_visual', 'Comment_uncertainties', 'Press_min', 'U_10min_knots', 'U_3sec_kntots']):

    S = pd.read_csv(Mediadir+file, sep=',', header= 27)


    new_cols= list(S.columns)
    new_cols[1] = 'Time'
    new_cols[6:8] = ['Diameter', 'Optional_diameter']
    new_cols[10:16] = ['Comment', 'Comment_visual', 'Comment_uncertainties', 'Press_min', 'U_10min_knots', 'U_3sec_kntots']
    S.columns= new_cols
    
    S.drop(droplist, axis=1, inplace=True )
    
    S['time'] = pd.to_datetime(S['Time'])
    S.drop(['Time'], axis=1, inplace=True)
    
    S= S.rename(columns={"Latitude": "lat", "Longitude": "lon"})
    
#    S['Season']= S['Season'].str.slice(0, 4).astype(int)
    S['Month']= S['time'].dt.month
    #S = S.set_index('time')
    S['ID']= S['ID'].replace({'98': '97'}) #98 is according to Rojo a continuation of 97
    S= S.replace({'comma': 'C', 'undefined': 'U'}) # some wrong morphologies
    S['Morphology']= S['Morphology'].replace({'  ': ' ', 'H - C': 'C - H', 'MGR - C': 'C - MGR', 'T - C': 'C - T', 'U - C': 'C - U', 'H - S':'S - H', 'T - S': 'S - T', 'U - S':'S - U', 'S - W': 'W - S', 'S - MGR': 'MGR - S' }, regex= True)

    S= S[~np.isnan(S.lat)] #removes two rows where lat, lon, time is nan
    return S






tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/tracks_"+test+".csv")
tracks['time']= pd.to_datetime(tracks['time'])

tracks= tracks.set_index(['ID', 'step'])




Rojo = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv")
Rojo['time']= Rojo['time'].dt.round('H')
Rojo= Rojo.set_index(['ID', 'Obs'])

Rojo= Rojo[(Rojo['time'] > time_start) & (Rojo['time'] < time_end)]

Rojo= Rojo.join( (Rojo.reset_index().groupby('ID').last()['Obs']).rename('nObs'), on='ID')

Rojo= Rojo[Rojo['nObs'] >=3]

# just get first index: Rojo.index.get_level_values(0))
Rojo_indexes= remove_dublicate(Rojo.index.get_level_values(0))



matchdist= 200
match= 0
match_track= []
match_track_Rojo= []
total= len(Rojo_indexes)

from collections import Counter

pd.options.mode.chained_assignment = None  # default='warn'

for index in Rojo_indexes: #['175']: #Rojo_indexes:
    
    Rojo_PL= Rojo.loc[index]
    
    
    matches= []
    print('\n Rojo nr ', index)
    for Obs in range(1, len(Rojo_PL)+1):
        Rojo_Obs= Rojo_PL.loc[Obs]
    
        track_now= tracks[tracks['time']== Rojo_Obs['time']]
    
        # track_now['distance'] =distance( (Rojo_Obs['lat'], Rojo_Obs['lon']) , (track_now['lat'], track_now['lon']) )
        track_now['distance'] = track_now.loc[track_now.index]['distance'] =distance( (Rojo_Obs['lat'], Rojo_Obs['lon']) , (track_now['lat'], track_now['lon']) )
        closest_track= track_now['distance'].idxmin()
        
        print(Obs, track_now.loc[closest_track].name, np.round(track_now.loc[closest_track]['distance']))
        
        matches += list(track_now[track_now['distance']< matchdist].index.get_level_values('ID'))
    
    c= Counter(matches)
    print(c)
    npc= np.array(list(c.items()) )
    
    if len(npc)== 0: print('no match')
    else:
        maxcount= np.max(npc[:,1])
        nObs= Rojo_PL['nObs'].values[0]
        print(maxcount, ' of ', nObs, 'detected assciated to track: ', npc[np.argmax(npc[:,1]), 0])
        if maxcount/nObs >=0.5:
            match +=1
            match_track += [npc[np.argmax(npc[:,1]), 0] ]
            match_track_Rojo += [index]
            

print('In total ', match, ' of ', total, ' were detected ', match/total)


match_track= remove_dublicate(match_track)
print(match - len(match_track), ' of ', match, 'are double, such that ', len(match_track), ' remain')

tracks_matchRojo= tracks.loc[match_track]


tracks_matchRojo["Rojo_match"] = ""

for i in range(len(match_track)):
    tracks_matchRojo.loc[match_track[i]]['Rojo_match']=match_track_Rojo[i]

tracks_matchRojo.to_csv(Mediadir+"/data/ERA5_Clim/track_lists/tracks_matchRojo_"+test+".csv")
