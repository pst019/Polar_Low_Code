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
    Mediadir= '/media/'+user+'/1692A00D929FEF8B/'

else:
    Mediadir= '/run/media/pst019/1692A00D929FEF8B/'
'/home/'+user+'/home/'

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

S['time'] = pd.to_datetime(S['Time'])
S.drop(['Time'], axis=1, inplace=True)

S['Season']= S['Season'].str.slice(0, 4).astype(int)
S['Month']= S['time'].dt.month
#S = S.set_index('time')
S['ID']= S['ID'].replace({'98': '97'}) #98 is according to Rojo a continuation of 97

S_ERA5 = pd.read_csv(Mediadir+"PL/STARS/STARS_ERA5.csv", sep=',', header= 0)
S= pd.concat([S, S_ERA5], axis= 1)

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
S_ind['Morphology']= [" ".join(sorted(remove_dublicate(remove_from_list(split_items(S[S['ID']== ID]['Morph_red'], ' - '), remove_list)))) for ID in S_ind['ID']]

S_ind.loc[S_ind['Morphology'].str.contains('MGR'), 'Morphology']= 'MGR+'
S_ind.loc[S_ind['Morphology'].str.contains('W'), 'Morphology']= 'W+'
S_ind['Morphology']= S_ind['Morphology'].replace({'C S': 'C-S'})


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



"""plot morphology"""
fig= plt.figure(fignr)
ax= fig.add_subplot(111)

fignr+=1
plt.clf()


for morph in pd.crosstab(S_ind['Month'], S_ind.PL_Morph).columns:
    color = next(ax._get_lines.prop_cycler)['color']
    S_morph= S_ind[S_ind['Morphology'] == morph]
#    print('Morph:', morph)

    S_morph_count= S_morph['PL_Morph_full'].value_counts()
    
    l= len(S_morph_count)
    if l == 1:
        plt.bar(morph, len(S_morph), color= color, edgecolor= 'k', lw= 1)

    if l > 1:
        count= 0
        for li in range(l):
#            print(S_morph_count.index[li])
            b= plt.bar(morph, S_morph_count[li], bottom= count,  color= color, edgecolor= 'k', lw= 1)#hatch= patterns[li],
            
                        
            if S_morph_count[li] > 4:
                plt.text(morph, count+ S_morph_count[li]/2 -.3,  S_morph_count.index[li], ha='center', va='center')

            count+= S_morph_count[li]

if save:
    savedir= homedir+ 'Polar_Low/STARS-ana/Figs/Sub-Morph'
    print(savedir)
    plt.savefig(savedir)            



"""-----------------------------------------------------------------------------"""
"""set variables"""
#save=True
save=False

#var, bins='Season', np.arange(1999, 2020,1)-.5
#var, bins='Month', np.arange(.5, 12.6, 1)
#var, bins= 'multiple', np.arange(10)+.5
#var, bins= 'Duration', np.arange(0, 116, 5)
#var, bins= 'Distance', np.arange(0,2101, 100)  #"""travelled distance - start-end"""
#var, bins= 'Prop_speed', np.arange(0, 100.1, 4)
#var, bins= 'Diameter', np.arange(0, 751, 50)
#var, bins= 'Diameter_max', np.arange(0, 1201, 50)
#var, bins= 'Diameter_min', np.arange(0, 651, 50)

#var, bins= 'U_400', np.arange(5, 30.1, 1)
#var, bins= 'skt_med_200', np.arange(258, 286, 1)
#var, bins= 'msl_min_50', np.arange(950, 1030.1)
#var, bins= 'blh_max_250', np.arange(1000, 3201.1, 100)
#var, bins= 'cape_max_250', np.arange(0, 301, 10)
#var, bins= 'sshf_mean_300', np.arange(0, 301, 10)
#var, bins= 'slhf_mean_300', np.arange(0, 301, 10)
var, bins= 'tp_mean_300', np.arange(0, .55, .02) #total precip


#var, bins= 'lat', np.arange(51, 82.1, 1) #55, 82
#var, bins= 'start_lat', np.arange(51, 82.1, 1) #58, 82
#var, bins= 'end_lat', np.arange(51, 82.1, 1)
#var, bins= 'lon', np.arange(-30, 70.1, 5) #-30, 55.1
#var, bins= 'start_lon', np.arange(-30, 70.1, 5) #-35, 55.1
#var, bins= 'end_lon', np.arange(-30, 70.1, 5) #-30, 85.1

#var= 'Prop_dir'


if var=='Season':
    S_ind['Season']= S_ind['Season'].astype(int)

if var=='Month':
    S_ind['Month']= S_ind['Month'].round().astype(int)

if var=='multiple':
    S_ind['PL_ID']= S_ind['ID'].str.replace(r'\D', '').astype(int)
    S_ind['multiple']= S_ind.groupby('PL_ID')['PL_ID'].transform('count')

if var in ['Duration', 'Prop_speed', 'Prop_dir']:
    #S_ind['start_time']
    starttime= S.loc[S['Obs'] == 1].groupby(['ID']).first()['time']
    S_ind= pd.concat([S_ind, starttime], axis= 1, sort=True)
    S_ind= S_ind.rename(columns={"time":"starttime"})
    
    #endtime
    endobs= S_ind['Obs']*2 -1
    endobs.loc[(endobs - endobs.astype(int)) > .01]= np.nan
    
    #endobs.index[0], endobs[0]
    #S[S['ID'] == endobs.index[0]]
    endtime= S[np.logical_and(S['ID'] == endobs.index[0], S['Obs'] == endobs[0])][['ID', 'time']]
    for i in range(1, len(endobs)):
        a= S[np.logical_and(S['ID'] == endobs.index[i], S['Obs'] == endobs[i])][['ID', 'time']]
        endtime= pd.concat([endtime, a])
    
    #endtime.rename(columns={"time":"endtime"})
    endtime= endtime.set_index('ID')
    
    S_ind= pd.concat([S_ind, endtime], axis= 1, sort=True)
    S_ind= S_ind.rename(columns={"time":"endtime"})   
    
    S_ind['Duration']= np.round((S_ind['endtime'] - S_ind['starttime'])/np.timedelta64(1,'h'))

if var in ['Distance', 'Prop_speed', 'Prop_dir', 'start_lat', 'start_lon', 'end_lat', 'end_lon']:
    S_ind= S_ind.rename(columns={"lat":"avg_lat", "lon":"avg_lon"})
    
    startloc= S.loc[S['Obs'] == 1].groupby(['ID']).first()[['lat', 'lon']]
    S_ind= pd.concat([S_ind, startloc], axis= 1, sort=True)
    S_ind= S_ind.rename(columns={"lat":"start_lat", "lon":"start_lon"})
    
    endobs= S_ind['Obs']*2 -1
    endobs.loc[(endobs - endobs.astype(int)) > .01]= np.nan
    
    endloc= S[np.logical_and(S['ID'] == endobs.index[0], S['Obs'] == endobs[0])][['ID', 'lat', 'lon']]
    for i in range(1, len(endobs)):
        a= S[np.logical_and(S['ID'] == endobs.index[i], S['Obs'] == endobs[i])][['ID', 'lat', 'lon']]
        endloc= pd.concat([endloc, a])
    
    #endtime.rename(columns={"time":"endtime"})
    endloc= endloc.set_index('ID')
    
    S_ind= pd.concat([S_ind, endloc], axis= 1, sort=True)
    S_ind= S_ind.rename(columns={"lat":"end_lat", "lon":"end_lon"})
    
    S_ind['Distance']= 110* np.sqrt( (S_ind['start_lat']- S_ind['end_lat'])**2+ (np.cos(np.deg2rad(0.5*(S_ind['start_lat']+ S_ind['end_lat'])))* (S_ind['start_lon']- S_ind['end_lon']))**2)


if var in ['Prop_speed', 'Prop_dir']: S_ind['Prop_speed']= S_ind['Distance']/S_ind['Duration']
    
if var=='Diameter_max': S_ind['Diameter_max']= S.groupby(['ID']).max()['Diameter']
if var=='Diameter_min': S_ind['Diameter_min']= S[S.Stage == 'M'].groupby(['ID']).min()['Diameter']
 
if 'U_' in var: S_ind[var]= S.groupby(['ID']).max()[var]
if 'skt_' in var: S_ind[var]= S.groupby(['ID']).median()[var]
if 'msl_' in var: S_ind[var]= S.groupby(['ID']).min()[var]/100


"""------------ plot the figure------------"""   
if var != 'Prop_dir':
    plt.figure(fignr)
    fignr+=1
    plt.clf() 
     
    if var not in ['Season', 'Month']: ax1= plt.subplot(211)
      
    bot=False
    for morph in pd.crosstab(S_ind[var], S_ind.PL_Morph).columns:
        botn,b_dn,p_dn= plt.hist(S_ind[S_ind['Morphology'] == morph][var], bins= bins , label=morph, bottom= bot, edgecolor= 'k', lw= 1)
        bot += botn
    #
    plt.xlim(bins[0], bins[-1])
    plt.ylabel('Number of polar lows')
    plt.legend()
    
    
    
    if var== 'Season': plt.xticks(bins[1::2]+.5)
    
    
    if var not in ['Season', 'Month']: ax1.xaxis.set_ticks_position('both')
    
    
    """-------------- plot density distributions ------------------"""
    if var not in ['Season', 'Month']:
        ax2= plt.subplot(212, sharex= ax1)
        
        from scipy import stats
        
        xplot= np.linspace(bins[0], bins[-1], 100)
        
        kernel= stats.gaussian_kde(S_ind[var].dropna())
        plt.plot(xplot, kernel(xplot), label='all', lw= 2)
        
        for morph in pd.crosstab(S_ind[var], S_ind.PL_Morph).columns[1:]:
            kernel= stats.gaussian_kde(S_ind[S_ind['Morphology'] == morph][var].dropna())
            plt.plot(xplot, kernel(xplot), label=morph)
        
        
        plt.legend()
        plt.ylabel('Normalized density')
        ax2.xaxis.set_ticks_position('both')
        
    
        
    plt.xlabel(var)
    if var== 'multiple': plt.xlabel('Number of simultaneous systems')
    if var=='Duration': plt.xlabel('Life time [h]')
    if var=='Distance': plt.xlabel('Distance between start and end point [km]')
    if var=='Prop_speed': plt.xlabel('Propagation speed [km/h]')
    if var=='Diameter': plt.xlabel('Diameter of the polar low [km]')
    if var=='Diameter_max': plt.xlabel('Maximum diameter of the polar low [km]')
    if var=='Diameter_min': plt.xlabel('Minimum diameter in mature stage of the polar low [km]')
    if var=='U_400': plt.xlabel('Maximum wind speed within 400 km radius [m/s]')
    if 'skt' in var: plt.xlabel('Median skin temperature within 200 km [K]')
    if 'msl' in var: plt.xlabel('Minimum sea level pressure within 50 km [hPa]')
    
    
    plt.tight_layout()


"""prop direction - wind rose"""
if var== 'Prop_dir':
    fignr+=1

    fig= plt.figure(fignr, figsize=(9.5, 6))
    fignr+=1
    plt.clf()
    
    from f_meteo import *
    from windrose import WindroseAxes #pip install git+https://github.com/python-windrose/windrose
    import matplotlib.gridspec as gridspec
    import matplotlib.ticker as tkr
    
    gs = gridspec.GridSpec(2, 3)     # (nblines, nbcol)
    # return lists of bottom and top position of rows, left and right positions of columns.

    bottom, top, left, right = gs.get_grid_positions(fig)  # [bottom, top, left, right]

    S_ind['Prop_dir']= [calculate_initial_compass_bearing((S_ind['start_lat'][i], S_ind['start_lon'][i]), (S_ind['end_lat'][i], S_ind['end_lon'][i]) ) for i in range(len(S_ind))]

    col, row= 0,0
    
    rect = [left[col],  bottom[row],  right[col]-left[col],  0.9*(top[row]-bottom[row])]
    
    ax = WindroseAxes(fig, rect)
    fig.add_axes(ax)

    ax.bar(S_ind['Prop_dir'], S_ind['Prop_speed'], normed=False, bins=np.arange(0, 61, 20), opening=0.9, edgecolor='white')
    ax.set_title('All', position=(0.1, 1.01), color= 'r', fontsize= 15)
#    ax.set_yticks(np.arange(0, 15.5, 5))
    ax.yaxis.set_major_formatter(tkr.FormatStrFormatter('%2.0f'))
    
    ax.set_legend()
    ax.legend(title="Speed [km/h]", loc= (-.3,0), decimal_places=0)

    
    for mi, morph in enumerate(pd.crosstab(S_ind[var], S_ind.PL_Morph).columns[1:]):
        col, row= (mi+1)%3, (mi+1)//3
#        print(mi, morph, col, row)
        rect = [left[col],  bottom[row],  right[col]-left[col],  0.9*(top[row]-bottom[row])]

        ax = WindroseAxes(fig, rect)
        fig.add_axes(ax)
    
        ax.bar(S_ind[S_ind['Morphology'] == morph]['Prop_dir'], S_ind[S_ind['Morphology'] == morph]['Prop_speed'], normed=False, bins=np.arange(0, 61, 20),  opening=0.9, edgecolor='white')
        ax.set_title(morph, position=(0.1, 1.01), color= 'r', fontsize= 15)
#        ax.set_yticks(np.arange(0, 15.5, 5))
        ax.yaxis.set_major_formatter(tkr.FormatStrFormatter('%2.0f'))    
    



    



if save:
    savedir= homedir+ 'Polar_Low/STARS-ana/Figs/'+'Morph_'+var
    print(savedir)
    plt.savefig(savedir, bbox_inches='tight')






