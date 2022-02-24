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

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from f_useful import *

fignr= 1


S = pd.read_csv(Mediadir+"PL/STARS/Rojo-etal_2019.csv", sep=',', header= 27)



new_cols= list(S.columns)
new_cols[1] = 'Time'
new_cols[6:8] = ['Diameter', 'Optional_diameter']
new_cols[10:16] = ['Comment', 'Comment_visual', 'Comment_uncertainties', 'Press_min', 'U_10min_knots', 'U_3sec_kntots']
S.columns= new_cols

S.drop(['Optional_diameter', 'Comment', 'Comment_visual', 'Comment_uncertainties', 'Press_min', 'U_10min_knots', 'U_3sec_kntots'],axis=1, inplace=True )

S['datetime'] = pd.to_datetime(S['Time'])
S.drop(['Time'], axis=1, inplace=True)

S['Season']= S['Season'].str.slice(0, 4).astype(int)
S['Month']= S['datetime'].dt.month
#S = S.set_index('datetime')

S= S.replace({'comma': 'C', 'undefined': 'U'}) # some wrong morphologies
S['ID']= S['ID'].replace({'98': '97'}) #98 is according to Rojo a continuation of 97

#S= S.set_index('ID')
#S.index= S.index.round('H')


S_ERA5 = pd.read_csv(Mediadir+"PL/STARS/STARS_ERA5.csv", sep=',', header= 0)
pd.concat([S, S_ERA5], axis= 1)



"""morphology"""
plt.figure(fignr)
fignr+=1
plt.clf()

morphs= S.groupby('Morphology').size()
morphs= morphs[morphs.values > 10]

morphs.plot(kind= 'bar')
#S = Spiraliform, C = Comma shaped, MGR = Merry-go-round, W = Wave system, U = Undefined, T = Transition between different forms, H = Hybrid, - = Systems don't appear completely on imagery
plt.ylabel('Number of timesteps')


"""diameter"""
#plt.figure(fignr)
#fignr+=1
#plt.clf()
#
#bindist= 50
#plt.hist(S['Diameter'], bins= np.arange(0, 1000, bindist), color= 'b', alpha= .5)
#
#from scipy import stats
#kernel= stats.gaussian_kde(S['Diameter'].dropna())
#n= len(S['Diameter'].dropna())
#plt.plot(np.arange(0, 1000, 10), bindist*n*kernel(np.arange(0, 1000, 10)), color= 'b', label='all')
#
#
#
#kernel= stats.gaussian_kde(S[S['Morphology'] == 'C']['Diameter'].dropna())
#n= len(S[S['Morphology'] == 'C']['Diameter'].dropna())
#plt.plot(np.arange(0, 1000, 10), bindist*n*kernel(np.arange(0, 1000, 10)), color='r', label='comma')
#
#
#plt.hist(pd.concat([S[S['Morphology'] == 'C']['Diameter'], S[S['Morphology'] == 'S']['Diameter']]), bins= np.arange(0, 1000, bindist), color= 'g', alpha= .5)
#kernel_1= stats.gaussian_kde(S[S['Morphology'] == 'S']['Diameter'].dropna())
#n_1= len(S[S['Morphology'] == 'S']['Diameter'].dropna())
#plt.plot(np.arange(0, 1000, 10), bindist*(n_1*kernel_1(np.arange(0, 1000, 10))+ n*kernel(np.arange(0, 1000, 10))), color='g', label= 'spirelli+ comma')
#
#plt.hist(S[S['Morphology'] == 'C']['Diameter'], bins= np.arange(0, 1000, bindist), color= 'r', alpha= .5)
#
#
#plt.xlim(0,1000)
#plt.xlabel('Diameter [km]')
#plt.ylabel('Number of timesteps')
#plt.legend()
#
#
#
#"""diameter distributions of different morphologies"""
#plt.figure(fignr)
#fignr+=1
#plt.clf()
#
#kernel= stats.gaussian_kde(S['Diameter'].dropna())
#plt.plot(np.arange(0, 1000, 10), kernel(np.arange(0, 1000, 10)), label='all')
#
#kernel= stats.gaussian_kde(S[S['Morphology'] == 'C']['Diameter'].dropna())
#plt.plot(np.arange(0, 1000, 10), kernel(np.arange(0, 1000, 10)), label='comma')
#
#kernel= stats.gaussian_kde(S[S['Morphology'] == 'MGR']['Diameter'].dropna())
#plt.plot(np.arange(0, 1000, 10), kernel(np.arange(0, 1000, 10)),  label= 'MGR')
#
#kernel= stats.gaussian_kde(S[S['Morphology'] == 'S']['Diameter'].dropna())
#plt.plot(np.arange(0, 1000, 10), kernel(np.arange(0, 1000, 10)),  label= 'spirelli')
#
#kernel= stats.gaussian_kde(S[S['Morphology'] == 'W']['Diameter'].dropna())
#plt.plot(np.arange(0, 1000, 10), kernel(np.arange(0, 1000, 10)),  label= 'wave')
#
#kernel= stats.gaussian_kde(S[S['Morphology'] == 'H']['Diameter'].dropna())
#plt.plot(np.arange(0, 1000, 10), kernel(np.arange(0, 1000, 10)),  label= 'hybrid')
#
#kernel= stats.gaussian_kde(S[S['Morphology'] == 'C - W']['Diameter'].dropna())
#plt.plot(np.arange(0, 1000, 10), kernel(np.arange(0, 1000, 10)),  label= 'comma-wave')
#
#plt.xlim(0,1000)
#plt.xlabel('Diameter [km]')
#plt.ylabel('Timestep density')
#plt.legend()
#
##
##
#
#"""get individual systems"""
#S_ind= S.groupby(['ID']).mean()
#S_ind['ID']= S_ind.index



"""system morphology"""
def remove_from_list(item_list, remove_list):
    remove_list= set(remove_list)
    return [e for e in item_list if e not in remove_list]

import itertools
def remove_repeater(inputlist):
    """remove item from list if it occurs several times in a row"""
    return [k for k,_g in itertools.groupby(inputlist)]

import operator
import functools
def split_items(inputlist, splitexpression):
    """splitexpression ' - ' splits ['C', 'C - U', 'U', 'C', '-'] to ['C', 'C', 'U', 'U', 'C', '-']"""
    listoflists= [a.split(' - ') for a in inputlist]
    return functools.reduce(operator.iconcat, listoflists, [])

#def replaceMultiple(mainString, toBeReplaces, newString):
#    """Replace a set of multiple sub strings with a new string in main string."""
#    # Iterate over the strings to be replaced
#    for elem in toBeReplaces :
#        # Check if string is in the main string
#        if elem in mainString :
#            # Replace the string
#            mainString = mainString.replace(elem, newString)
#    
#    return  mainString


# the almost original
#remove_list= ['-']
#S_ind['PL_Morph']= [" ".join(remove_from_list(remove_dublicate(S[S['ID'] == ID]['Morphology']), remove_list)) for ID in S_ind['ID']]


# respecting the order
remove_list= ['-', 'T', 'U', 'H'] #maybe should not remove hyrbid

#removes repeaters - respects the order of the development
#S['Morph_red']= S['Morphology'].replace({'T - ': '', ' - T': '', '  ': ' ', 'U - ': '', ' - U': '', ' - H': '', 'H - ': ''}, regex= True)
#S_ind['PL_Morph']= [" ".join(remove_repeater(split_items(remove_repeater(remove_from_list(S[S['ID'] == ID]['Morph_red'], remove_list)), ' - '))) for ID in S_ind['ID']]

#removes dublicates - sorted is doing alphabetic ordering
S['Morph_red']= S['Morphology'].replace({'  ': ' '}, regex= True)
S_ind['PL_Morph']= [" ".join(sorted(remove_dublicate(remove_from_list(split_items(S[S['ID']== ID]['Morph_red'], ' - '), remove_list)))) for ID in S_ind['ID']]


S_ind.loc[S_ind['PL_Morph'].str.contains('MGR'), 'PL_Morph']= 'MGR+'
S_ind.loc[S_ind['PL_Morph'].str.contains('W'), 'PL_Morph']= 'W+'


plt.figure(fignr)
fignr+=1
plt.clf()

morphs= S_ind.groupby('PL_Morph').size()
#morphs= morphs[morphs.values > 50]
morphs= morphs.drop(index= '') #remove the 3 PLs that have no signature

morphs.plot(kind= 'bar')
#S = Spiraliform, C = Comma shaped, MGR = Merry-go-round, W = Wave system, U = Undefined, T = Transition between different forms, H = Hybrid, - = Systems don't appear completely on imagery
plt.ylabel('Number of polar lows')


plt.tight_layout()



"""season"""
#plt.figure(fignr)
#fignr+=1
#plt.clf()
#
#S_ind['Season']= S_ind['Season'].astype(int)
#S_ind.groupby('Season').nunique()['ID'].plot(kind='bar', color= 'b') #the longitude seem to be unique, but it would be better by ID
#plt.ylabel('Number of polar lows')
#
#plt.tight_layout()
#
#
"""month"""
#plt.figure(fignr)
#fignr+=1
#plt.clf()
#
#S_ind['Month']= S_ind['Month'].round().astype(int)
##S_ind.groupby('Month').nunique()['ID'].plot(kind= 'bar', color= 'b')
#
#
#bot=False
#for morph in pd.crosstab(S_ind.Month, S_ind.PL_Morph).columns:
#    botn,b_dn,p_dn= plt.hist(S_ind[S_ind['PL_Morph'] == morph]['Month'], bins= np.arange(.5, 12.6, 1), label=morph, bottom= bot)
#    bot += botn
#
#plt.xlim(.5, 12.5)
#plt.xlabel('Month')
#plt.ylabel('Number of polar lows')
#plt.legend()
#plt.tight_layout()




"""multiple"""
#S_ind['PL_ID']= S_ind['ID'].str.replace(r'\D', '').astype(int)
#S_ind['multiple']= S_ind.groupby('PL_ID')['PL_ID'].transform('count')
#
#S_strong= S_ind.groupby('PL_ID').first()
#
#plt.figure(fignr)
#fignr+=1
#plt.clf()
#
#plt.hist(S_strong['multiple'], bins= np.arange(11)-.5)
#plt.hist(S_ind['multiple'], bins= np.arange(11)-.5)

#plt.xlim(0.5, 10)
#plt.xlabel('Number of simultaneous systems')
#plt.ylabel('Number of polar lows')#occurances


"""multiple in categories"""
#plt.close('all')
#plt.figure(fignr)
#fig1, ax1 = plt.subplots(figsize=(10,4))
#plt.close(fignr)
#fignr+=1
#plt.clf()




#tab= pd.crosstab(S_ind.multiple, S_ind.PL_Morph)
#tab.plot.bar(stacked=True)
#
##plt.xlim(0.5, 10)
#plt.xlabel('Number of simultaneous systems')
#plt.ylabel('Number of polar lows')
#plt.tight_layout()


"""duration"""
plt.figure(fignr)
fignr+=1
plt.clf()

#S_ind['start_time']
starttime= S.loc[S['Obs'] == 1].groupby(['ID']).first()['datetime']
S_ind= pd.concat([S_ind, starttime], axis= 1, sort=True)
S_ind= S_ind.rename(columns={"datetime":"starttime"})

#endtime
endobs= S_ind['Obs']*2 -1
endobs.loc[(endobs - endobs.astype(int)) > .01]= np.nan

#endobs.index[0], endobs[0]
#S[S['ID'] == endobs.index[0]]
endtime= S[np.logical_and(S['ID'] == endobs.index[0], S['Obs'] == endobs[0])][['ID', 'datetime']]
for i in range(1, len(endobs)):
    a= S[np.logical_and(S['ID'] == endobs.index[i], S['Obs'] == endobs[i])][['ID', 'datetime']]
    endtime= pd.concat([endtime, a])

#endtime.rename(columns={"datetime":"endtime"})
endtime= endtime.set_index('ID')

S_ind= pd.concat([S_ind, endtime], axis= 1, sort=True)
S_ind= S_ind.rename(columns={"datetime":"endtime"})


S_ind['Duration']= np.round((S_ind['endtime'] - S_ind['starttime'])/np.timedelta64(1,'h'))

#plt.hist(S_ind['Duration'], bins= np.arange(112))
#plt.xlim(0, 112)
#plt.xlabel('Duration [h]')
#plt.ylabel('Number of polar lows')

tab= pd.crosstab(S_ind.Duration, S_ind.PL_Morph)
#tab.plot.bar(stacked=True)

bindist= 5

bot=False
for morph in tab.columns:
    botn,b_dn,p_dn= plt.hist(S_ind[S_ind['PL_Morph'] == morph]['Duration'], bins= np.arange(0, 116, bindist), label=morph, bottom= bot)
    bot += botn

plt.legend()
#for index in tab.columns:
#    plt.bar(tab.index, tab[index], width= 3)#, bins= np.arange(0, 112, 2), stacked=True)

plt.xlim(0, 116)
plt.xlabel('Duration [h]')
plt.ylabel('Number of polar lows')


"""travelled distance - start-end"""
plt.figure(fignr)
fignr+=1
plt.clf()

S_ind= S_ind.rename(columns={"Latitude":"avg_lat", "Longitude":"avg_lon"})


startloc= S.loc[S['Obs'] == 1].groupby(['ID']).first()[['Latitude', 'Longitude']]
S_ind= pd.concat([S_ind, startloc], axis= 1, sort=True)
S_ind= S_ind.rename(columns={"Latitude":"start_lat", "Longitude":"start_lon"})

endobs= S_ind['Obs']*2 -1
endobs.loc[(endobs - endobs.astype(int)) > .01]= np.nan

endloc= S[np.logical_and(S['ID'] == endobs.index[0], S['Obs'] == endobs[0])][['ID', 'Latitude', 'Longitude']]
for i in range(1, len(endobs)):
    a= S[np.logical_and(S['ID'] == endobs.index[i], S['Obs'] == endobs[i])][['ID', 'Latitude', 'Longitude']]
    endloc= pd.concat([endloc, a])

#endtime.rename(columns={"datetime":"endtime"})
endloc= endloc.set_index('ID')

S_ind= pd.concat([S_ind, endloc], axis= 1, sort=True)
S_ind= S_ind.rename(columns={"Latitude":"end_lat", "Longitude":"end_lon"})

S_ind['Distance']= 110* np.sqrt( (S_ind['start_lat']- S_ind['end_lat'])**2+ (np.cos(np.deg2rad(0.5*(S_ind['start_lat']+ S_ind['end_lat'])))* (S_ind['start_lon']- S_ind['end_lon']))**2)
#
#plt.hist(S_ind['Distance'], bins= np.arange(0,2100, 50))


bot=False
for morph in tab.columns:
    botn,b_dn,p_dn= plt.hist(S_ind[S_ind['PL_Morph'] == morph]['Distance'], bins= np.arange(0, 2100, 50), label=morph, bottom= bot)
    bot += botn

plt.legend()


plt.xlim(0, 2100)
plt.xlabel('Distance [km]')
plt.ylabel('Number of polar lows')


"""prop speed"""
plt.figure(fignr)
fignr+=1
plt.clf()


S_ind['Prop_speed']= S_ind['Distance']/S_ind['Duration']
#plt.hist(S_ind['Prop_speed'], bins= np.arange(0,100, 4))

bot=False
for morph in tab.columns:
    botn,b_dn,p_dn= plt.hist(S_ind[S_ind['PL_Morph'] == morph]['Prop_speed'], bins= np.arange(0, 100, 4), label=morph, bottom= bot)
    bot += botn

plt.legend()


plt.xlim(0, 100)
plt.xlabel('Propagation speed [km/h]')
plt.ylabel('Number of polar lows')
#
#
"""prop direction - wind rose"""
#plt.figure(fignr)
#fignr+=1
#plt.clf()
#
#from f_meteo import *
#from windrose import WindroseAxes #pip install git+https://github.com/python-windrose/windrose
#
#S_ind['Prop_dir']= [calculate_initial_compass_bearing((S_ind['start_lat'][i], S_ind['start_lon'][i]), (S_ind['end_lat'][i], S_ind['end_lon'][i]) ) for i in range(len(S_ind))]
#
#ax = WindroseAxes.from_ax()
#ax.bar(S_ind['Prop_dir'], S_ind['Prop_speed'], normed=True, opening=0.8, edgecolor='white')
#ax.set_legend()




