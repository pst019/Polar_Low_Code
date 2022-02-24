#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 16:29:52 2019

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

if user=='pst019':
#    homedir= '/home/'+user+'/home/'
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/1692A00D929FEF8B/'

homedir=Mediadir+'home/'


import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/test_Denis/mc_era5-master/code') #to get obs_tracks_api

#import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from octant.core import TrackRun, OctantTrack
from pathlib import Path
from obs_tracks_api import prepare_tracks

from f_carto import *
from f_useful import *
from f_STARS import *
from datetime import timedelta, datetime
import scipy.ndimage.filters as filters





"""get STARS list"""
R = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv", droplist= ['Optional_diameter', 'Comment_visual', 'Comment_uncertainties', 'U_3sec_kntots'])
R_ind= R.groupby(['ID']).first()

R_N= R_ind[R_ind.lon < 20].index.values
R_B= R_ind[R_ind.lon >= 20].index.values

Stoll_STARS_file="ERA5_STARS/STARS-PMC-merge/"+test+"_dist_"+str(dist)+"_handfix.csv"
Stoll = pd.read_csv(Mediadir+ Stoll_STARS_file)
Stoll['time']= pd.DatetimeIndex(Stoll.time)
Stoll= Stoll_Obs_nr(Stoll) #creates Stoll['Obs']

Stoll= Stoll_excl(Stoll, lifetime_threshold= 5, double= False) #, landfraction= 0.5 )


S_ind= Stoll.groupby(['Stoll nr']).first()['Rojo nr'].values

R_N_detected = list(set(list(R_N)).intersection(list(S_ind)) )
print('detection Norwegian Sea:', len(R_N_detected)/len(R_N), len(R_N_detected), len(R_N) )

R_B_detected = list(set(list(R_B)).intersection(list(S_ind)) )
print('detection Barents Sea:', len(R_B_detected)/len(R_B) , len(R_B_detected), len(R_B) )


print('total detection rate:', len( remove_dublicate(S_ind) )/len(R_ind), len( remove_dublicate(S_ind) ), len(R_ind) )