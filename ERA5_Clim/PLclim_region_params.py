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
params = {'axes.labelsize': 16,
          # 'axes.titlesize': 16,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16
          #, 'legend.fontsize': 13
          }
plt.rcParams.update(params)
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
savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/PL_region_params/'


PL_timestep= False # True #only the time steps that satisfy all PL crit





fignr= 3

hem= 'NH'

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

system_based= True
# var, system_char, varname, ylim = 'shear_strength_925-500_mean500', 'max', r'Shear strength [$10^{-3}$s$^{-1}$]', [0,9]
# var, system_char, varname, ylim= 'land_dist', 'max', 'Distance to land [km]', [0, 4000]

# var, system_char, varname, ylim= 'theta_trop_mean250', 'mean', r'$\theta_{trop}$ [K]', [260, 302]
# var, system_char, varname, ylim= 'SST-T_500_mean250', 'max', r'SST - T$_{500hPa}$ [K]', [36, 56]
var, system_char, bins, varname, ylim='theta_diff_500-sst_mean250', 'min', np.arange(-8, 42.1, 1), r'$\theta_{500hPa} - \theta_{SST}$ [K]' , [-8, 11.2] #'max'

# var, system_char, varname, ylim= 'vo', 'max', 'Relative vorticity [10$^{-5}$s$^{-1}$]', [18, 90]
# var, system_char, varname, ylim= 'diameter', 'max', 'Diameter [km]', [50, 435]

# var, system_char, varname, ylim= 'duration', False, 'Lifetime [h]', [0, 160]


# system_based= False
# var, varname= 'diameter', 'Diameter [km]'
# var, varname= 'SST-T_500_mean250', r'SST - T$_{500hPa}$ [K]'
# var, varname= 'land_dist_correct', 'Distance to land [km]'



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
 
if hem== 'NH':
    boxlist= [
          ['Nordic Seas', [-15, 70, 80, 55] ] #W, E, N, S
        , ['Irminger Sea', [-45, -15, 75, 50] ]
        , ['Labrador Sea', [-70, -45, 75, 50] ]
        , ['G. of Alaska', [200, 235, 65, 40] ]       
        , ['Bering Sea', [160, 200, 65, 40] ]       
        , ['S. of Okhotsk', [142, 160, 65, 40] ]       
        , ['S. of Japan',  [127, 142, 50, 35] ]
        ]
elif hem == 'SH':
    # boxlist= [ ['Amundsen Sea', [-140, -50, -45, -70] ] #W, E, N, S
    # , ['New Zealand', [150, 210, -45, -70] ]
    # , ['Indian Ocean', [40, 120, -45, -70] ]
    # ]

    boxlist= [
    ['SE Pacific', [-140, -60, -38, -70] ] #Amundsen Sea #W, E, N, S
    , ['SW Pacific', [150, 220, -38, -70] ]
    , ['Indian Oc.', [20, 150, -38, -70] ]
    , ['Atlantic', [-60, 20, -38, -70] ]
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
tracks['shear_strength_925-500_mean500']*= 1E3

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







if var == 'diameter': tracks['diameter']= np.sqrt(tracks['area']) *4/np.pi

if var != 'duration': tracks.loc[tracks[var] == -1000, var] = np.nan

"""system based"""
if system_based:
    tracks_ind= tracks.groupby(['ID']).mean()
   
    if var == 'duration': tracks_ind[var]= tracks.groupby(['ID']).max()['step']   

    if system_char == 'max':
        tracks_ind[var]= tracks.groupby(['ID']).max()[var]
    if system_char == 'min':
        tracks_ind[var]= tracks.groupby(['ID']).min()[var]
    tracks= tracks_ind


"""make the regions"""
tracks['region']= 'Rest'

for box_name, box_loc in boxlist:
    tracks.loc[((tracks.lon <= box_loc[1]) & (tracks.lon > box_loc[0]) &
                          (tracks.lat <= box_loc[2]) & (tracks.lat > box_loc[3]) ), 'region'] = box_name





if hem == 'NH': fig= plt.figure(fignr, figsize= (9,6) )
if hem == 'SH': fig= plt.figure(fignr-1, figsize= (5,6) )

fignr+=1
plt.clf()
ax= fig.add_subplot()

#to get the boxes in the correct order
boxnamelist= [box[0] for box in boxlist] +['Rest']
tracks['region_numbered']= 10

for ib, boxname in enumerate(boxnamelist):
    # print(ib, boxname)
    tracks.loc[tracks.region== boxname, 'region_numbered']= ib


tracks.boxplot(column = var, by = 'region_numbered', ax= ax, grid=False, rot= 23, fontsize= 16,
       color = {'whiskers' : 'black',
                         'caps' : 'black',
                         'medians' : 'black',
                         'boxes' : 'black'},
        whiskerprops = {'linewidth' : 2},
        capprops = {'linewidth' : 2},
        flierprops = {'markersize': 1, 'marker':'o', 'markeredgecolor':'grey'},
        medianprops = {'linewidth' : 2, 'color': 'k'},
        boxprops = {'linewidth' : 2})# , patch_artist = True)

ax.set_xticklabels(boxnamelist)
fig.suptitle('')
plt.title('')

if ylim: plt.ylim(ylim)
plt.xlabel('')
plt.ylabel(varname, fontsize= 14)

   
if save:
    savefile= savedir+ f'{hem}_{var}'
    if system_based:
        savefile += '_system'
        if system_char: savefile += '-'+system_char
    # if PL_timestep== True: savefile += '_PL-ts'



    for vc in var_char_ind:
        importvar, importchar, excl_thresh, excl_type = vc
        savefile += '_'+ importvar + '-'+ str(int(excl_thresh))     
    savefile += '.png'
    print(savefile)
    fig.savefig(savefile , bbox_inches='tight', dpi= 150)   




print(f'Time passed: {round(time.perf_counter() - start,2)} seconds')