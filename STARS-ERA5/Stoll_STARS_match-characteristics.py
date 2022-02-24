#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 10:01:45 2019

@author: pst019
"""


import time
start = time.time()


import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    homedir= '/home/'+user+'/home/'
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/1692A00D929FEF8B/'

homedir=Mediadir+'home/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import xarray as xr
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import stats


import numpy as np
from f_useful import *
from f_STARS import *


lifetime= 5



"""get STARS list"""
STARS = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv", droplist= ['Optional_diameter', 'Comment_visual', 'Comment_uncertainties', 'U_3sec_kntots'])
#renaming for the trackmatch

STARS_ind= STARS_individual_systems(STARS)

STARS_ind= calc_system_duration(STARS, STARS_ind, ID_name= 'ID', Obs_name='Obs')


STARS_event_list_prior_time = remove_dublicate([int(ID.split('.')[0]) for ID in STARS_ind.index])
STARS_centre_list_prior_time = remove_dublicate(STARS_ind.index)

print('nr STARS events:', len(STARS_event_list_prior_time)) #98 is merged to 97 according to Rojo
print('nr STARS centres:', len(STARS_centre_list_prior_time))


#multi_events= [int(ID.split('.')[0]) for ID in STARS_ind.index]
#import itertools
#itertools.groupby(multi_events)


"""get multiple statistics"""
#STARS_ind['PL_ID']= STARS_ind['ID'].str.replace(r'\D', '').astype(int)
#STARS_ind['multiple']= STARS_ind.groupby('PL_ID')['PL_ID'].transform('count')
#print(STARS_ind.groupby('multiple').size())


"""after lifetime excl"""
#STARS_ind= STARS_ind[STARS_ind.Duration >= lifetime]
#STARS_event_list = remove_dublicate([int(ID.split('.')[0]) for ID in STARS_ind.index])
#STARS_centre_list = remove_dublicate(STARS_ind.index)
#
#print('nr STARS events after lifetime excl:', len(STARS_event_list))
#print('nr STARS centres after lifetime excl:', len(STARS_centre_list))







"""Stoll systems"""
test= 'version4'
dist= 150

Stoll_STARS_file="ERA5_STARS/STARS-PMC-merge/"+test+"_dist_"+str(dist)+"_handfix.csv"
Stoll = pd.read_csv(Mediadir+ Stoll_STARS_file)
Stoll['time']= pd.DatetimeIndex(Stoll.time)

Stoll= Stoll_Obs_nr(Stoll) #creates Stoll['Obs']

"""Stoll PMC track bits"""
Stoll['PMC_number']= Stoll['track file']* 1E3 + Stoll['track_idx'] + [int(ID.split('_')[0])*1E9 + int(ID.split('_')[1])*1E6 for ID in Stoll['Stoll nr'] ]

#to count how many track bits are obtained by the PMC algorithm that are merged to stoll then
print('Number of PMC track parts in Stoll:', len(remove_dublicate(Stoll['PMC_number'])))

print('nr Stoll centres:', len(remove_dublicate(Stoll['Stoll nr'])) )


Stoll_STARS_centres= remove_dublicate(Stoll['Rojo nr'])
print('nr STARS centres included in Stoll:', len(Stoll_STARS_centres) )



#to count how many Rojo centres ha


#print('Stoll exclude')
#Stoll= Stoll_excl(Stoll, give_print=False)   
#Stoll = merge_ERA5(Stoll, Mediadir, "PL/STARS/Stoll_ERA5_handfix_3.csv") #the ERA5 merge has to happen with the same stoll list as was used for the production


#Stoll_event_list_prior_time = remove_dublicate([int(ID.split('_')[0]) for ID in Stoll['Stoll nr']])
#print('nr Stoll events prior time exclusion:', len(Stoll_event_list_prior_time)) 
Stoll_STARS_events_prior_time= remove_dublicate([int(ID.split('.')[0]) for ID in Stoll['Rojo nr']])
print('nr STARS events included in Stoll:', len(Stoll_STARS_events_prior_time) )



print('\n time exclusion')
Stoll= Stoll_excl_lifetime(Stoll, lifetime_threshold = 5)

#Stoll= Stoll.rename(columns={'Stoll nr': 'ID', 'Rojo nr': 'Rojo_ID', 'Rojo nr old': 'Rojo_ID_orig',
#                             'STARS Obs nr': 'Rojo_Obs', 'STARS lat': 'Rojo_lat', 'STARS lon': 'Rojo_lon',
#                             'Press_min': 'Rojo_SLP'})
#Stoll= Stoll.drop(columns=['Unnamed: 21']) # 'STARS lat', 'STARS lon', 'STARS Obs nr', 'dist'])



#Stoll_STARS_events= remove_dublicate([int(ID.split('.')[0]) for ID in Stoll['Rojo_ID']])
#print('nr STARS events included in Stoll:', len(Stoll_STARS_events) )


print('nr Stoll centres:', len(remove_dublicate(Stoll['Stoll nr'])) )



print('\n Stoll land exclude')
Stoll= Stoll_excl(Stoll, double=False) #, give_print=False)   


print('nr Stoll centres:', len(remove_dublicate(Stoll['Stoll nr'])) )


Stoll_STARS_centres= remove_dublicate(Stoll['Rojo nr'])
print('nr STARS centres included in Stoll:', len(Stoll_STARS_centres) )


Stoll_event_list = remove_dublicate([int(ID.split('_')[0]) for ID in Stoll['Stoll nr']])
print('nr Stoll events:', len(Stoll_event_list)) 




print('\n Stoll exclude')
Stoll= Stoll_excl(Stoll) #, give_print=False)   


print('nr Stoll centres:', len(remove_dublicate(Stoll['Stoll nr'])) )


Stoll_STARS_centres= remove_dublicate(Stoll['Rojo nr'])
print('nr STARS centres included in Stoll:', len(Stoll_STARS_centres) )


Stoll_event_list = remove_dublicate([int(ID.split('_')[0]) for ID in Stoll['Stoll nr']])
print('nr Stoll events:', len(Stoll_event_list)) 


"""intersection between"""
#all STARS events, only Stoll with lifetime >5
#list1= STARS_event_list_prior_time
#list2= Stoll_STARS_events

#only STARS with lifetime >5, all Stoll
#list1= STARS_event_list
#list2= Stoll_STARS_events_prior_time

#all STARS events, all Stoll events
#list1= STARS_event_list_prior_time
#list2= Stoll_STARS_events_prior_time

#only STARS with lifetime >5, only Stoll with lifetime >5
#list1= STARS_event_list
#list2= Stoll_STARS_events
#
#intersect= list(set(list1).intersection(list2))
#print('STARS:', len(list1), 'STARS and Stoll:', len(intersect), 'fraction: ', len(intersect)/len(list1))

"""collected Stoll list"""
#imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'
#
#
##Stoll= pd.read_csv(imp_dir + 'Stoll_list_noRojo.csv')
#Stoll= pd.read_csv(imp_dir + 'Stoll_list.csv')
#Stoll= Stoll.set_index(['ID', 'time'])
##
#
#S_ERA5= pd.read_csv(imp_dir + 'Stoll_ERA5_1.csv')
#S_ERA5= S_ERA5.set_index(['ID', 'time'])
#
#S_ERA5_2= pd.read_csv(imp_dir + 'Stoll_ERA5_2.csv')
#S_ERA5_2= S_ERA5_2.set_index(['ID', 'time'])
#
#Stoll= pd.concat([Stoll, S_nodes, S_ERA5, S_ERA5_2], axis=1)


#Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps



