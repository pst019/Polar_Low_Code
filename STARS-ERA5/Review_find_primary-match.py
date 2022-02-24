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
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/PatsOrange/'

homedir=Mediadir+'home/'


import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

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

#Stoll= Stoll_excl(Stoll)   


Stollnrs= remove_dublicate(Stoll['Stoll nr'])
print('Stoll centers:' , len(Stollnrs) ) 

Rojonrs= remove_dublicate(Stoll['Rojo nr'])
print('of different Rojo centers:', len(Rojonrs) ) 

Rojoevents= remove_dublicate([ID.split('.')[0] for ID in Stoll['Rojo nr'] ])
print('of different Rojo events:', len(Rojoevents) ) 


Rojo_primaries= [nr for nr in Rojonrs if '.' not in nr]
print('Rojo primaries:',len(Rojo_primaries))

Rojo_others= [nr for nr in Rojonrs if '.' in nr]
print('Rojo secondaries:', len(Rojo_others))


Stoll= Stoll_excl(Stoll, lifetime_threshold= 5, double= False, land= True)
#Stoll= Stoll_excl(Stoll, lifetime_threshold= 5, double=True, land= False)   
#Stoll= Stoll_excl(Stoll)#, double= False, landfraction= 0.5 )


print('Stoll centers:' , len(remove_dublicate(Stoll['Stoll nr']) ) ) 


Rojonrs= remove_dublicate(Stoll['Rojo nr'])
print('of different Rojo centers:', len(Rojonrs) ) 

Rojoevents= remove_dublicate([ID.split('.')[0] for ID in Stoll['Rojo nr'] ])
print('of different Rojo events:', len(Rojoevents) ) 

Rojo_primaries= [nr for nr in Rojonrs if '.' not in nr]
print('Rojo primaries:', len(Rojo_primaries))

S_ind= Stoll.groupby(['Stoll nr']).first()

print('Stoll centers associated to Rojo primaries:', len([nr for nr in S_ind['Rojo nr'] if '.' not in nr]) )


Rojo_others= [nr for nr in Rojonrs if '.' in nr]
print('Rojo secondaries:', len(Rojo_others))

print('Stoll centers associated to Rojo secondaries:', len([nr for nr in S_ind['Rojo nr'] if '.' in nr]) )
