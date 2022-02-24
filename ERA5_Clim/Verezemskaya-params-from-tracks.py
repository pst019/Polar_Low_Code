#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 14:08:13 2021

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/PatsOrange/'


import sys  #to import the functions from a different directory
sys.path.insert(0, Mediadir+ 'home/Polar_Low/polar_low_code/Functions')

from f_useful import *
import pandas as pd
import numpy as np



version= "" #"fram_run3_"
ending= 'SH180'

start_year= 1979
end_year= 2020

durlim= 6
matchdist= 150

# time_diff= False
time_diff= 24 #does not include the matched track if the matched track starts more than time_diff earlier or ends more than time_diff later

more_info= True
post_process= True
# split, exceed_large_to_small= True, 40
split= False

terrain_dist, terrain_thresh= 2, 50


PL_name= 'Verezemskaya'



"""csv name"""
csv_name= 'merged_tracks_'+version + ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
if more_info: csv_name += '_moreinfo'
if post_process: csv_name += '_post-process'
if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
# if split: csv_name += f'_split{exceed_large_to_small}'
# if params: csv_name += '_params'
print(csv_name)

"""track list"""
time_start_SH= '2004-1-1' #the vorticity and slp is needed in this time
time_end_SH= '2005-1-1'

# tracks_params= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks2/"+csv_name.replace("atl-pac", "SH180")+'_params.csv')
tracks_params= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks_tomerge/"+csv_name+f'_params_{time_start_SH}-{time_end_SH}.csv')


"""PL list"""
PLfile= Mediadir+"/data/ERA5_Clim/track_lists/matchPLlist2/matchPLlist_"+ PL_name +f"_dist{matchdist}"
if time_diff and PL_name != 'Yanase': PLfile += f'_timediff{time_diff}'
if PL_name == 'Verezemskaya': PLfile += f'_{csv_name.replace("atl-pac", "SH180")}.csv'
else: PLfile += f'_{csv_name}.csv'

PLs=pd.read_csv(PLfile)


PLs_params= PLs.merge(tracks_params).set_index(['ID', 'step'])

PLs_params= PLs_params.replace(np.nan, 'nan', regex=True) #empty space for the theta500-thetaSST
    # PLs= PLs.replace(1000, 'nan', regex=True) #for the SST-T500

outfile= PLfile.replace(".csv", "_params.csv")
print('\n write:', outfile)

PLs_params.to_csv(outfile)
