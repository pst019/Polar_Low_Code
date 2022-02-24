#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 14:56:09 2020

@author: pst019
"""

import time
start = time.time()


import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    Mediadir= '/media/'+user+'/Backup/'

else:
    Mediadir= '/run/media/pst019/Backup/'

Dropboxdir=Mediadir+'Dropbox/'


import sys  #to import the functions from a different directory
sys.path.insert(0, Dropboxdir+ 'Polar_Low/polar_low_code/Functions')

import pandas as pd 
import numpy as np
from f_useful import *
from f_STARS import *


save= False
#save= True
#save_dir= Mediadir + '/ERA5_STARS/Stoll_list_other-landexcl/'
save_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


#Stype='system'
Stype='step'

list_type='Stoll'


system_char= 'mean'
#system_char= 'max'
#system_char='min'
#system_char='med'
#system_char='initial'
#system_char='final'




  
"""Stoll systems"""
test= 'version4'
dist= 150

Stoll_STARS_file="ERA5_STARS/STARS-PMC-merge/"+test+"_dist_"+str(dist)+"_handfix.csv"
Stoll = pd.read_csv(Mediadir+ Stoll_STARS_file)
Stoll['time']= pd.DatetimeIndex(Stoll.time)

Stoll= Stoll_Obs_nr(Stoll) #creates Stoll['Obs']