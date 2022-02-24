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

Mediadir= '/run/media/pst019/1692A00D929FEF8B/'
homedir= Mediadir+'home/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from f_useful import *



save= True
#save= False

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
##S = Spiraliform, C = Comma shaped, MGR = Merry-go-round, W = Wave system, U = Undefined, T = Transition between different forms, H = Hybrid, - = Systems don't appear completely on imagery

S= S.replace({'comma': 'C', 'undefined': 'U'}) # some wrong morphologies
remove_list= ['-', 'T', 'U', 'H'] #maybe should not remove hyrbid


S['Morph_red']= S['Morphology'].replace({'T - ': '', ' - T': '', '  ': ' ', 'U - ': '', ' - U': '', ' - H': '', 'H - ': ''}, regex= True)
S_ind['PL_Morph_full']= [" ".join(remove_repeater(split_items(remove_repeater(remove_from_list(S[S['ID'] == ID]['Morph_red'], remove_list)), ' - '))) for ID in S_ind['ID']]

#splits transition, remove dublicates, sort by alphabet
S['Morph_red']= S['Morphology'].replace({'  ': ' '}, regex= True)
S_ind['PL_Morph']= [" ".join(sorted(remove_dublicate(remove_from_list(split_items(S[S['ID']== ID]['Morph_red'], ' - '), remove_list)))) for ID in S_ind['ID']]

S_ind.loc[S_ind['PL_Morph'].str.contains('MGR'), 'PL_Morph']= 'MGR+'
S_ind.loc[S_ind['PL_Morph'].str.contains('W'), 'PL_Morph']= 'W+'
S_ind['PL_Morph']= S_ind['PL_Morph'].replace({'C S': 'C-S'})


"""full morphology"""

#S_ind.loc[S_ind['PL_Morph_full'].str.contains(' MGR '), 'PL_Morph_full']= '..MGR..'
S_ind.loc[S_ind['PL_Morph_full'].str.contains(' MGR'), 'PL_Morph_full']= '.. MGR'
S_ind.loc[S_ind['PL_Morph_full'].str.contains('MGR '), 'PL_Morph_full']= 'MGR ..'

#S_ind.loc[np.logical_and.reduce((S_ind['PL_Morph_full'] != 'C W C' , S_ind['PL_Morph_full'] != 'C W S' , S_ind['PL_Morph_full'].str.contains(' W '))), 'PL_Morph_full']= '..W..'
#S_ind.loc[np.logical_and(~S_ind['PL_Morph_full'].str.contains(' W ')  , S_ind['PL_Morph_full'].str.contains(' W')), 'PL_Morph_full']= '..W'

#S_ind.loc[S_ind['PL_Morph_full'].str.contains(' W '), 'PL_Morph_full']= '..W..'
#S_ind.loc[S_ind['PL_Morph_full'].str.contains('W '), 'PL_Morph_full']= 'W..'
#S_ind.loc[S_ind['PL_Morph_full'].str.contains(' W'), 'PL_Morph_full']= '..W'

#S_ind.loc[S_ind['PL_Morph_full'].str.contains('C S C'), 'PL_Morph_full']= 'CSC+'
#S_ind.loc[np.logical_or(S_ind['PL_Morph_full'] == 'S C S' , S_ind['PL_Morph_full'].str.contains('S C S C')), 'PL_Morph_full']= '. S C S ..'
#S_ind.loc[S_ind['PL_Morph_full'].str.contains('S C S'), 'PL_Morph_full']= 'CSC+'




from f_carto import *
import cartopy.crs as ccrs

fig = plt.figure(fignr, figsize=(11,6))
fignr+=1
plt.clf()


ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (2,3,1))
scale_bar(ax, 500, location= (0.06, 0.04))



"""scatter plot of individual steps"""
ax.scatter(S.Longitude, S.Latitude, transform=ccrs.PlateCarree(), s= 1)
#ax.scatter(S[S['Morphology'] == 'C'].Longitude, S[S['Morphology'] == 'C'].Latitude, transform=ccrs.PlateCarree(), s= 1)
#ax.scatter(S[S['Morphology'] == 'MGR'].Longitude, S[S['Morphology'] == 'MGR'].Latitude, transform=ccrs.PlateCarree(), s= 1)
#ax.scatter(S[S['Morphology'] == 'S'].Longitude, S[S['Morphology'] == 'S'].Latitude, transform=ccrs.PlateCarree(), s= 1)

ax.set_title('All', position=(0.1, 1.01), color= 'k', fontsize= 15)


"""tracks of systems"""
morph_groups= pd.crosstab(S_ind['Month'], S_ind.PL_Morph).columns

for mi, morph in enumerate(morph_groups[1:]):
    
    ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (2,3,mi+2))
    
    ID_list= S_ind[S_ind.PL_Morph == morph].ID
    ax.set_title(morph, position=(0.1, 1.01), color= 'k', fontsize= 15)

    for ID in ID_list:
    
        S_ID= S[S.ID == ID]
        ax.scatter(S_ID[S_ID.Obs== 1].Longitude, S_ID[S_ID.Obs== 1].Latitude, transform= ccrs.PlateCarree(), s= 4)
        ax.plot(S_ID.Longitude, S_ID.Latitude, transform= ccrs.PlateCarree(), lw= 1)


plt.tight_layout()

if save:
    savedir= homedir+ 'Polar_Low/Arome-Arctic-operational/STARS-ana/2/'+'Morph_maps'
    print(savedir)
    plt.savefig(savedir, bbox_inches='tight')


