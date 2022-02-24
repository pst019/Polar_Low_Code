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
Mediadir= '/media/'+user+'/1692A00D929FEF8B/'

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
S['ID']= S['ID'].replace({'98': '97'}) #98 is according to Rojo a continuation of 97



"""get individual systems"""
S_ind= S.groupby(['ID']).mean()
S_ind['ID']= S_ind.index


"""system morphology"""
S= S.replace({'comma': 'C', 'undefined': 'U'}) # some wrong morphologies
remove_list= ['-', 'T', 'U', 'H'] #maybe should not remove hyrbid


S['Morph_red']= S['Morphology'].replace({'T - ': '', ' - T': '', '  ': ' ', 'U - ': '', ' - U': '', ' - H': '', 'H - ': ''}, regex= True)
S_ind['PL_Morph_full']= [" ".join(remove_repeater(split_items(remove_repeater(remove_from_list(S[S['ID'] == ID]['Morph_red'], remove_list)), ' - '))) for ID in S_ind['ID']]

#splits transition, remove dublicates, sort by alphabet
S['Morph_red']= S['Morphology'].replace({'  ': ' '}, regex= True)
S_ind['PL_Morph']= [" ".join(sorted(remove_dublicate(remove_from_list(split_items(S[S['ID']== ID]['Morph_red'], ' - '), remove_list)))) for ID in S_ind['ID']]

S_ind.loc[S_ind['PL_Morph'].str.contains('MGR'), 'PL_Morph']= 'MGR+'
S_ind.loc[S_ind['PL_Morph'].str.contains('W'), 'PL_Morph']= 'W+'


"""morphology plot"""
plt.figure(fignr)
fignr+=1
plt.clf()

morphs= S_ind.groupby('PL_Morph').size()
#morphs= morphs[morphs.values > 50]
#morphs= morphs.drop(index= '') #remove the 3 PLs that have no signature

morphs.plot(kind= 'bar')
#S = Spiraliform, C = Comma shaped, MGR = Merry-go-round, W = Wave system, U = Undefined, T = Transition between different forms, H = Hybrid, - = Systems don't appear completely on imagery
plt.ylabel('Number of polar lows')

plt.tight_layout()

plt.figure(fignr)
fignr+=1
plt.clf()

#S_ind.loc[S_ind['PL_Morph_full'].str.contains(' MGR '), 'PL_Morph_full']= '..MGR..'
S_ind.loc[S_ind['PL_Morph_full'].str.contains(' MGR'), 'PL_Morph_full']= '..MGR'
S_ind.loc[S_ind['PL_Morph_full'].str.contains('MGR '), 'PL_Morph_full']= 'MGR..'

#S_ind.loc[np.logical_and.reduce((S_ind['PL_Morph_full'] != 'C W C' , S_ind['PL_Morph_full'] != 'C W S' , S_ind['PL_Morph_full'].str.contains(' W '))), 'PL_Morph_full']= '..W..'
#S_ind.loc[np.logical_and(~S_ind['PL_Morph_full'].str.contains(' W ')  , S_ind['PL_Morph_full'].str.contains(' W')), 'PL_Morph_full']= '..W'

#S_ind.loc[S_ind['PL_Morph_full'].str.contains(' W '), 'PL_Morph_full']= '..W..'
#S_ind.loc[S_ind['PL_Morph_full'].str.contains('W '), 'PL_Morph_full']= 'W..'
#S_ind.loc[S_ind['PL_Morph_full'].str.contains(' W'), 'PL_Morph_full']= '..W'

#S_ind.loc[S_ind['PL_Morph_full'].str.contains('C S C'), 'PL_Morph_full']= 'CSC+'
S_ind.loc[np.logical_or(S_ind['PL_Morph_full'] == 'S C S' , S_ind['PL_Morph_full'].str.contains('S C S C')), 'PL_Morph_full']= 'S C S ..'
#S_ind.loc[S_ind['PL_Morph_full'].str.contains('S C S'), 'PL_Morph_full']= 'CSC+'


morphs= S_ind.groupby('PL_Morph_full').size()
#morphs= morphs[morphs.values > 50]
#morphs= morphs.drop(index= '') #remove the 3 PLs that have no signature

morphs.plot(kind= 'bar')
#S = Spiraliform, C = Comma shaped, MGR = Merry-go-round, W = Wave system, U = Undefined, T = Transition between different forms, H = Hybrid, - = Systems don't appear completely on imagery
plt.ylabel('Number of polar lows')

plt.tight_layout()



fig= plt.figure(fignr)
fignr+=1

plt.clf()
ax= fig.add_subplot(111)

patterns = ('-', '+', 'x', '\\', '*', 'o', 'O', '.')

for morph in pd.crosstab(S_ind['Month'], S_ind.PL_Morph).columns:
    color = next(ax._get_lines.prop_cycler)['color']
    S_morph= S_ind[S_ind['PL_Morph'] == morph]
    print('Morph:', morph)

    S_morph_count= S_morph['PL_Morph_full'].value_counts()
    
    l= len(S_morph_count)
    if l == 1:
        plt.bar(morph, len(S_morph), color= color, edgecolor= 'k', lw= 1)

    if l > 1:
        count= 0
        for li in range(l):
            print(S_morph_count.index[li])
            b= plt.bar(morph, S_morph_count[li], bottom= count,  color= color, edgecolor= 'k', lw= 1)#hatch= patterns[li],
            
                        
            if S_morph_count[li] > 4:
                plt.text(morph, count+ S_morph_count[li]/2 -.3,  S_morph_count.index[li], ha='center', va='center')

            count+= S_morph_count[li]

            
#    botn,b_dn,p_dn= plt.hist(S_ind[S_ind['PL_Morph'] == morph][var], bins= bins , label=morph, bottom= bot)
#    bot += botn


#
#
#"""-----------------------------------------------------------------------------"""
#"""set variables"""
#save=True
##save=False
#
##var, bins='Season', np.arange(1999, 2020,1)-.5
#var, bins='Month', np.arange(.5, 12.6, 1)
##var, bins= 'multiple', np.arange(10)+.5
##var, bins= 'Duration', np.arange(0, 116, 5)
##var, bins= 'Distance', np.arange(0,2101, 100)  #"""travelled distance - start-end"""
##var, bins= 'Prop_speed', np.arange(0, 100.1, 4)
#
#if var=='Season':
#    S_ind['Season']= S_ind['Season'].astype(int)
#
#if var=='Month':
#    S_ind['Month']= S_ind['Month'].round().astype(int)
#
#if var=='multiple':
#    S_ind['PL_ID']= S_ind['ID'].str.replace(r'\D', '').astype(int)
#    S_ind['multiple']= S_ind.groupby('PL_ID')['PL_ID'].transform('count')
#
#if var in ['Duration', 'Prop_speed']:
#    #S_ind['start_time']
#    starttime= S.loc[S['Obs'] == 1].groupby(['ID']).first()['datetime']
#    S_ind= pd.concat([S_ind, starttime], axis= 1, sort=True)
#    S_ind= S_ind.rename(columns={"datetime":"starttime"})
#    
#    #endtime
#    endobs= S_ind['Obs']*2 -1
#    endobs.loc[(endobs - endobs.astype(int)) > .01]= np.nan
#    
#    #endobs.index[0], endobs[0]
#    #S[S['ID'] == endobs.index[0]]
#    endtime= S[np.logical_and(S['ID'] == endobs.index[0], S['Obs'] == endobs[0])][['ID', 'datetime']]
#    for i in range(1, len(endobs)):
#        a= S[np.logical_and(S['ID'] == endobs.index[i], S['Obs'] == endobs[i])][['ID', 'datetime']]
#        endtime= pd.concat([endtime, a])
#    
#    #endtime.rename(columns={"datetime":"endtime"})
#    endtime= endtime.set_index('ID')
#    
#    S_ind= pd.concat([S_ind, endtime], axis= 1, sort=True)
#    S_ind= S_ind.rename(columns={"datetime":"endtime"})   
#    
#    S_ind['Duration']= np.round((S_ind['endtime'] - S_ind['starttime'])/np.timedelta64(1,'h'))
#
#if var in ['Distance', 'Prop_speed']:
#    S_ind= S_ind.rename(columns={"Latitude":"avg_lat", "Longitude":"avg_lon"})
#    
#    startloc= S.loc[S['Obs'] == 1].groupby(['ID']).first()[['Latitude', 'Longitude']]
#    S_ind= pd.concat([S_ind, startloc], axis= 1, sort=True)
#    S_ind= S_ind.rename(columns={"Latitude":"start_lat", "Longitude":"start_lon"})
#    
#    endobs= S_ind['Obs']*2 -1
#    endobs.loc[(endobs - endobs.astype(int)) > .01]= np.nan
#    
#    endloc= S[np.logical_and(S['ID'] == endobs.index[0], S['Obs'] == endobs[0])][['ID', 'Latitude', 'Longitude']]
#    for i in range(1, len(endobs)):
#        a= S[np.logical_and(S['ID'] == endobs.index[i], S['Obs'] == endobs[i])][['ID', 'Latitude', 'Longitude']]
#        endloc= pd.concat([endloc, a])
#    
#    #endtime.rename(columns={"datetime":"endtime"})
#    endloc= endloc.set_index('ID')
#    
#    S_ind= pd.concat([S_ind, endloc], axis= 1, sort=True)
#    S_ind= S_ind.rename(columns={"Latitude":"end_lat", "Longitude":"end_lon"})
#    
#    S_ind['Distance']= 110* np.sqrt( (S_ind['start_lat']- S_ind['end_lat'])**2+ (np.cos(np.deg2rad(0.5*(S_ind['start_lat']+ S_ind['end_lat'])))* (S_ind['start_lon']- S_ind['end_lon']))**2)
#
#
#if var=='Prop_speed': S_ind['Prop_speed']= S_ind['Distance']/S_ind['Duration']
#    
# 
#
#
#"""------------ plot the figure------------"""   
#plt.figure(fignr)
#fignr+=1
#plt.clf() 
# 
#if var not in ['Season', 'Month']: ax1= plt.subplot(211)
#  
#bot=False
#for morph in pd.crosstab(S_ind[var], S_ind.PL_Morph).columns:
#    botn,b_dn,p_dn= plt.hist(S_ind[S_ind['PL_Morph'] == morph][var], bins= bins , label=morph, bottom= bot)
#    bot += botn
##
#plt.xlim(bins[0], bins[-1])
#plt.ylabel('Number of polar lows')
#plt.legend()
#
#
#
#if var== 'Season': plt.xticks(bins[1::2]+.5)
#
#
#if var not in ['Season', 'Month']: ax1.xaxis.set_ticks_position('both')
#
#
#"""-------------- plot density distributions ------------------"""
#if var not in ['Season', 'Month']:
#    ax2= plt.subplot(212, sharex= ax1)
#    
#    from scipy import stats
#    
#    xplot= np.linspace(bins[0], bins[-1], 100)
#    
#    kernel= stats.gaussian_kde(S_ind[var].dropna())
#    plt.plot(xplot, kernel(xplot), label='all', lw= 2)
#    
#    for morph in pd.crosstab(S_ind[var], S_ind.PL_Morph).columns[1:]:
#        kernel= stats.gaussian_kde(S_ind[S_ind['PL_Morph'] == morph][var].dropna())
#        plt.plot(xplot, kernel(xplot), label=morph)
#    
#    
#    plt.legend()
#    plt.ylabel('Normalized density')
#    ax2.xaxis.set_ticks_position('both')
#    
#
#    
#plt.xlabel(var)
#if var== 'multiple': plt.xlabel('Number of simultaneous systems')
#if var=='Duration': plt.xlabel('Life time [h]')
#if var=='Distance': plt.xlabel('Distance between start and end point [km]')
#if var=='Prop_speed': plt.xlabel('Propagation speed [km/h]')
#
#
#plt.tight_layout()
#
#
#if save:
#    savedir= homedir+ 'Polar_Low/Arome-Arctic-operational/STARS-ana/'+'Morph_'+var
#    print(savedir)
#    plt.savefig(savedir)
#
#
#
#"""prop direction - wind rose"""
##plt.figure(fignr)
##fignr+=1
##plt.clf()
##
##from f_meteo import *
##from windrose import WindroseAxes #pip install git+https://github.com/python-windrose/windrose
##
##S_ind['Prop_dir']= [calculate_initial_compass_bearing((S_ind['start_lat'][i], S_ind['start_lon'][i]), (S_ind['end_lat'][i], S_ind['end_lon'][i]) ) for i in range(len(S_ind))]
##
##ax = WindroseAxes.from_ax()
##ax.bar(S_ind['Prop_dir'], S_ind['Prop_speed'], normed=True, opening=0.8, edgecolor='white')
##ax.set_legend()
#



