#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 15:46:14 2020

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

import xarray as xr
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib as mpl
# params = {'axes.labelsize': 15,
#           # 'axes.titlesize': 16,
#           'xtick.labelsize': 15,
#           'ytick.labelsize': 15
#           #, 'legend.fontsize': 13
#           }
# plt.rcParams.update(params)
#import matplotlib.colors as colors
from scipy import ndimage


import numpy as np
from f_useful import *
# from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs
import time

start= time.perf_counter()
#
save= True
# save= False
savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/Shear_params/'


PL_timestep= True # True #only the time steps that satisfy all PL crit

fignr= 3

hem= 'NH'
weak_shear_crit = 10 #1E-4 - criteria for the weak shear class
shear_mean_dist= 250


if hem == 'NH': ending= 'atl-pac'
elif hem == 'SH':  ending= 'SH180'

version= ""
start_year= 1979
end_year= 2020
terrain_thresh= 50
terrain_dist = 2
more_info= True
post_process= True
split, exceed_large_to_small= False, 40

durlim= 6

ylim= False

# system_based= True - does not make soo much sense for this variable
# var, system_char, varname, ylim = 'shear_strength_925-500_mean500', 'max', r'Shear strength [$10^{-3}$s$^{-1}$]', [0,9]
# var, system_char, varname, ylim= 'land_dist', 'max', 'Distance to land [km]', [0, 4000]

# var, system_char, varname, ylim= 'theta_trop_mean250', 'min', r'$\theta_{trop}$ [K]', [260, 302]
# var, system_char, varname, ylim= 'SST-T_500_mean250', 'max', r'SST - T$_{500hPa}$ [K]', [36, 56]
# var, system_char, bins, varname, ylim='theta_diff_500-sst_mean250', 'min', np.arange(-8, 42.1, 1), r'$\theta_{500hPa} - \theta_{SST}$ [K]' , [-8, 11.2] #'max'

# var, system_char, varname, ylim= 'vo', 'max', 'Relative vorticity [10$^{-5}$s$^{-1}$]', [18, 90]
# var, system_char, varname, ylim= 'diameter', 'max', 'Diameter [km]', [50, 435]

# var, system_char, varname, ylim= 'duration', False, 'Lifetime [h]', [0, 160]

# 
system_based= False
# var, varname, ylim= 'diameter', 'Diameter [km]', [0, 435]
# var, varname, ylim= 'SST-T_500_mean250', r'SST - T$_{500hPa}$ [K]', [36, 56]
var, varname, ylim='theta_diff_500-sst_mean250', r'$\theta_{500hPa} - \theta_{SST}$ [K]', [-8, 11.2]

# var, varname, ylim= 'vo', 'Relative vorticity [10$^{-5}$s$^{-1}$]', [18, 90]
# var, varname, ylim = f'shear_strength_925-500_mean{shear_mean_dist}', r'Shear strength [$10^{-3}$s$^{-1}$]', [0,9]


# var, varname, ylim= 'theta_trop_mean250',  r'$\theta_{trop}$ [K]', [260, 302]
# var, varname, ylim= 'land_dist', 'Distance to land [km]', [0, 4000]
# var, varname = 'Lifetime_fraction', 'Lifetime fraction'


time_start= False #in this case the first and last index of alltracks will be used
# time_start= '2008-1-1' #the vorticity and slp is needed in this time
# time_end= '2009-1-1'

# var_char_ind = [['SST-T_500_mean250', 'max', 41.5, 'larger'],
#             ['theta_trop_mean250', 'min', 298, 'smaller'],
#             # ['land_dist', 'max', 140, 'larger'],
#             ['vo', 'max', 20, 'larger'],
#             ['diameter', 'max', 430, 'smaller'],
#             ]

var_char_ind = [['theta_diff_500-sst_mean250', 'min', 11.0, 'smaller'],
            ['theta_trop_mean250', 'mean', 300.8, 'smaller'],
            ['vo', 'max', 20, 'larger'],
            ['diameter', 'max', 430, 'smaller']
            ]    
 



"""import PL list"""
csv_name= 'merged_tracks_'+version + ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
if more_info: csv_name += '_moreinfo'
if post_process: csv_name += '_post-process'
if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
# if split: csv_name += f'_split{exceed_large_to_small}'

PLfile= Mediadir+"/data/ERA5_Clim/track_lists/derived_PLs/PLs-from_"+csv_name
for vc in var_char_ind:
    importvar, importchar, excl_thresh, excl_type = vc
    PLfile += '_'+ importvar+ '-' + str(int(excl_thresh))
if time_start != False: PLfile += f'_{time_start}-{time_end}'
PLfile += '.csv'
print('Import:', PLfile)
tracks= pd.read_csv(PLfile)    

tracks['vo'] *= 100
tracks[f'shear_strength_925-500_mean{shear_mean_dist}']*= 1E3



"""only steps that satisfy PL-IC"""
if PL_timestep== True:
    for vc in var_char_ind:
        importvar, importchar, excl_thresh, excl_type = vc
        if importvar == 'diameter': continue
    
        if excl_type == 'larger':
            tracks= tracks[tracks[importvar] > excl_thresh] 
        elif excl_type== 'smaller':
            tracks= tracks[tracks[importvar] < excl_thresh] 
        else: print('wrong excl_type')


# tracks= tracks.dropna(axis= 'columns') #drops some columns where parameters are only calculated for some years

"""make the shear categories"""
dir_var= f'shear_angle_925-500_mean{shear_mean_dist}'
speed_var= f'shear_strength_925-500_mean{shear_mean_dist}'

tracks['shear_cat']= 0
# tracks.loc[np.logical_and(tracks[dir_var] < 3*np.pi/4, tracks[dir_var] >= 1*np.pi/4), 'shear_cat'] ='forward'
# tracks.loc[np.logical_and(tracks[dir_var] < 5*np.pi/4, tracks[dir_var] >= 3*np.pi/4), 'shear_cat'] ='right'
# tracks.loc[np.logical_and(tracks[dir_var] < 7*np.pi/4, tracks[dir_var] >= 5*np.pi/4), 'shear_cat'] ='reverse'
# tracks.loc[np.logical_or(tracks[dir_var] < 1*np.pi/4, tracks[dir_var] >= 7*np.pi/4), 'shear_cat'] ='left'
tracks.loc[np.logical_or(tracks[dir_var] < 45, tracks[dir_var] >= 315), 'shear_cat'] ='forward'
tracks.loc[np.logical_and(tracks[dir_var] < 135, tracks[dir_var] >= 45), 'shear_cat'] ='right'
tracks.loc[np.logical_and(tracks[dir_var] < 225, tracks[dir_var] >= 135), 'shear_cat'] ='reverse'
tracks.loc[np.logical_and(tracks[dir_var] < 315, tracks[dir_var] >= 225), 'shear_cat'] ='left'

tracks.loc[tracks[speed_var] <= weak_shear_crit  *1E-1 , 'shear_cat'] = 'weak'


shear_list=['forward', 'right', 'reverse', 'left', 'weak']




if var == 'diameter': tracks['diameter']= np.sqrt(tracks['area']) *4/np.pi

elif var == 'Lifetime_fraction':
    tracks['Lifetime_fraction']= -1* np.ones(len(tracks))
    PLnr= remove_dublicate(tracks.ID)
    for ID in PLnr:
        track_now= tracks[tracks.ID == ID]
        times= pd.to_datetime(track_now.index)
        dt= (times- times[0])/np.timedelta64(1,'h')
        Lifetime_frac= dt/dt[-1]    
        tracks.loc[tracks.ID == ID, 'Lifetime_fraction']= Lifetime_frac

elif var != 'duration': tracks.loc[tracks[var] == -1000, var] = np.nan




fig= plt.figure(fignr, figsize= (5,6) )
# if hem == 'SH': fig= plt.figure(fignr-1, figsize= (5,6) )

fignr+=1
plt.clf()
ax= fig.add_subplot()

#to get the boxes in the correct order
# tracks['shear_numbered']= 10

# for i, shearname in enumerate(shear_list):
#     # print(ib, boxname)
#     tracks.loc[tracks.shear_cat== shearname, 'shear_numbered']= i

# tracks= tracks[[var, 'shear_cat']].dropna(axis= 'rows')

tracks.boxplot(column = var, by = 'shear_cat', ax= ax, grid=False, rot= 20, fontsize= 14,
       color = {'whiskers' : 'black',
                         'caps' : 'black',
                         'medians' : 'black',
                         'boxes' : 'black'},
        whiskerprops = {'linewidth' : 2},
        capprops = {'linewidth' : 2},
        flierprops = {'markersize': 1, 'marker':'o', 'markeredgecolor':'grey'},
        medianprops = {'linewidth' : 2, 'color': 'k'},
        boxprops = {'linewidth' : 2})# , patch_artist = True)

# ax.set_xticklabels(shear_list)
fig.suptitle('')
plt.title('')

if ylim: plt.ylim(ylim)
plt.xlabel('')
plt.ylabel(varname, fontsize= 14)

plt.tight_layout()
   
if save:
    savefile= savedir+ f'{hem}_{var}'
    if weak_shear_crit != 10: savefile += f'_weak-shear-crit{weak_shear_crit}'
    if shear_mean_dist != 500: savefile += f'_shear-mean{shear_mean_dist}'

    if PL_timestep== True: savefile += '_PL-ts'

    for vc in var_char_ind:
        importvar, importchar, excl_thresh, excl_type = vc
        savefile += '_'+ importvar + '-'+ str(int(excl_thresh))     
    savefile += '.png'
    print('save: \n ', savefile)
    fig.savefig(savefile , bbox_inches='tight', dpi= 150)   




print(f'Time passed: {round(time.perf_counter() - start,2)} seconds')