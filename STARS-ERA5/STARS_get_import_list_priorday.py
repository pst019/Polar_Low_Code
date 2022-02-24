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
'/home/'+user+'/home/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from f_useful import *
from datetime import timedelta

fignr= 1


S = pd.read_csv(Mediadir+"PL/STARS/Rojo-etal_2019.csv", sep=',', header= 27)


new_cols= list(S.columns)
new_cols[1] = 'Time'
new_cols[6:8] = ['Diameter', 'Optional_diameter']
new_cols[10:16] = ['Comment', 'Comment_visual', 'Comment_uncertainties', 'Press_min', 'U_10min_knots', 'U_3sec_kntots']
S.columns= new_cols

S.drop(['Optional_diameter', 'Comment', 'Comment_visual', 'Comment_uncertainties', 'Press_min', 'U_10min_knots', 'U_3sec_kntots'],axis=1, inplace=True )

S['time'] = pd.to_datetime(S['Time'])
S.drop(['Time'], axis=1, inplace=True)

S['Season']= S['Season'].str.slice(0, 4).astype(int)
S['Month']= S['time'].dt.month
#S = S.set_index('time')

S= S.replace({'comma': 'C', 'undefined': 'U'}) # some wrong morphologies
S['ID']= S['ID'].replace({'98': '97'}) #98 is according to Rojo a continuation of 97


S['dt_round']= [dt.replace(hour= 0, minute= 0) for dt in S['time']]

#S_day_list= S.groupby('dt_round').nunique()
S_day_list= pd.to_datetime(S.dt_round.unique() )
S_day_list= np.sort(S_day_list)
S_day_list= pd.DatetimeIndex(S_day_list)


#print(S_day_list[10].year)
#print(S_day_list[:10].month)
#print(S_day_list[:10].day)


"""normal daylist"""
#outF = open("STARS_daylist.txt", "w")
#for si in range(len(S_day_list)):
#  # write line to output file
#  outF.write(S_day_list[si].strftime("%Y %m %d"))
#  outF.write("\n")
#outF.close()


"""stars priorday"""
S_day_list_priorday= []

for S_day in S_day_list:
    if S_day - timedelta(1) not in S_day_list_priorday:
        S_day_list_priorday.append(S_day - timedelta(1))
    if S_day not in S_day_list_priorday:
        S_day_list_priorday.append(S_day)


#outF = open("STARS_daylist_priorday.txt", "w")
#for si in range(len(S_day_list_priorday)):
#  # write line to output file
#  outF.write(S_day_list_priorday[si].strftime("%Y %m %d"))
#  outF.write("\n")
#outF.close()
        

"""stars diff priorday - daylist"""
#S_day_list_diff_priorday= []
#for S_day in S_day_list:
#    if S_day - timedelta(1) not in S_day_list:
#        if S_day - timedelta(1) not in S_day_list_diff_priorday:
#            S_day_list_diff_priorday.append(S_day - timedelta(1))


#outF = open("STARS_daylist_diff-priorday.txt", "w")
#for si in range(len(S_day_list_diff_priorday)):
#  # write line to output file
#  outF.write(S_day_list_diff_priorday[si].strftime("%Y %m %d"))
#  outF.write("\n")
#outF.close()  
        



        
"""track list from stars priorday list"""
#outF = open("track_list_PMC.txt", "w")
#
#
#tracking_id= 0
#for S_day in S_day_list_priorday:
#    if S_day - timedelta(1) not in S_day_list_priorday:
#        tracking_id += 1
#        startdate= S_day
#        outF.write(str(tracking_id).zfill(3))
#        outF.write(startdate.strftime(" %Y %m %d "))
#    if S_day + timedelta(1) not in S_day_list_priorday:
#        enddate= S_day
#        outF.write(enddate.strftime("%Y %m %d"))
#        outF.write("\n")
#
#outF.close()




"""get the exclude dictionaries start and end from 'Plot_new_STARS' """
startdate= list(excl_dict_start.values())
startdate= [sdate - timedelta(days=1) for sdate in startdate]
enddate= list(excl_dict_end.values())
enddate= [sdate + timedelta(hours=1) for sdate in enddate]


append_date= startdate+ enddate
append_date= remove_dublicate(append_date)

#outF = open("txt_files/STARS_daylist_append.txt", "w")
#for si in append_date:
#    # write line to output file
#    outF.write(si.strftime("%Y %m %d"))
#    outF.write("\n")
#outF.close()



"""make track list from priorday list and append list"""
complete_daylist= pd.DatetimeIndex(S_day_list_priorday+ append_date).sort_values()

outF = open("txt_files/track_list_PMC_appent.txt", "w")


tracking_id= 0
for S_day in complete_daylist:
    if S_day - timedelta(1) not in complete_daylist:
        tracking_id += 1
        startdate= S_day
        outF.write(str(tracking_id).zfill(3))
        outF.write(startdate.strftime(" %Y %m %d "))
    if S_day + timedelta(1) not in complete_daylist:
        enddate= S_day
        outF.write(enddate.strftime("%Y %m %d"))
        outF.write("\n")

outF.close()