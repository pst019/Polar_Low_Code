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







#track_list= ['2013_03']
#track_list= ['run3_2011_01_atl']
#test= track_list[0]

year= 2011
month= 1

test= 'run3_'+str(year)+'_'+str(month).zfill(2)+'_atl'

track_dir= Mediadir+"/data/ERA5_Clim/tracks/"+test+'/'



all_files = glob.glob(track_dir + "/vortrack_*1.txt")

li = []



for filename in all_files:
    df = pd.read_csv(filename, sep='\s+', names=['lon', 'lat', 'vo', 'time', 'area', 'vortex_type', 'slp'])

    df['ID']= int(filename.split("/")[-1].split('_')[1])
    df['step']= np.arange(len(df))
    
    li.append(df)
    
    
tracks = pd.concat(li, axis=0, ignore_index=True)

tracks['time']= pd.to_datetime(tracks['time'], format='%Y%m%d%H%M')
tracks['ID']= (year*1E6+ month*1E4 + tracks.ID).astype('int')

tracks= tracks.set_index(['ID', 'step'])
tracks= tracks.sort_index()



# #create the duration
# tracks= tracks.join( (tracks.reset_index().groupby('ID').last()['step']).rename('Duration'), on='ID')

# tracks= tracks[tracks['Duration'] > 6]

tracks.to_csv(Mediadir+"/data/ERA5_Clim/track_lists/my_PC/tracks_"+test+".csv")




