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
save= False
savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/PL_shear_distribution/'

write= False

import_PLs= False #if True the matched PL list is imported
Plot_take_time= True

if import_PLs: write= False

per_dist = 250 #km

# version='excl_prop-dist'

fignr= 1

#track_dir= Mediadir+"/data/ERA5_Clim/track_lists/fram/"

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



time_start= False #in this case the first and last index of alltracks will be used
# time_start= '2008-1-1' #the vorticity and slp is needed in this time
# time_end= '2009-1-1'

var_char_ind = [['SST-T_500_mean250', 'max', 41.5, 'larger'],
            ['theta_trop_mean250', 'min', 298, 'smaller'],
            ['land_dist', 'max', 140, 'larger'],
            ['vo', 'max', 20, 'larger'],
            ['diameter', 'max', 430, 'smaller'],
            ]

    
if hem== 'NH':
    boxlist= [ ['Nordic Seas', [-15, 70, 80, 55] ] #W, E, N, S
        , ['Irminger Sea', [-45, -15, 75, 50] ]
        , ['Labrador Sea', [-70, -45, 75, 50] ]
        , ['Gulf of Alaska', [200, 235, 65, 40] ]       
        , ['Bering Sea', [160, 200, 65, 40] ]       
        , ['Sea of Okhotsk', [140, 160, 65, 40] ]       
        , ['Sea of Japan',  [130, 140, 50, 35] ]
        ]

elif hem == 'SH':
    boxlist= [ ['Amundsen Sea', [-140, -50, -45, -70] ] #W, E, N, S
    , ['New Zealand', [150, 210, -45, -70] ]
    , ['Indian Ocean', [40, 120, -45, -70] ]
    ]

# elif hem == 'SH':
#     boxlist= [ ['Amundsen Seas', [-110, -60, -45, -70] ] #W, E, N, S
#     , ['New Zealand', [160, 180, -45, -70] ]
#     , ['Australia', [80, 140, -45, -70] ]
#     ]

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

tracks= tracks.dropna(axis= 'columns')


"""make the shear categories"""
dir_var= 'shear_angle_925-500_mean500'
speed_var= 'shear_strength_925-500_mean500'

tracks['shear_cat']= 0
# tracks.loc[np.logical_and(tracks[dir_var] < 3*np.pi/4, tracks[dir_var] >= 1*np.pi/4), 'shear_cat'] ='forward'
# tracks.loc[np.logical_and(tracks[dir_var] < 5*np.pi/4, tracks[dir_var] >= 3*np.pi/4), 'shear_cat'] ='right'
# tracks.loc[np.logical_and(tracks[dir_var] < 7*np.pi/4, tracks[dir_var] >= 5*np.pi/4), 'shear_cat'] ='reverse'
# tracks.loc[np.logical_or(tracks[dir_var] < 1*np.pi/4, tracks[dir_var] >= 7*np.pi/4), 'shear_cat'] ='left'
tracks.loc[np.logical_or(tracks[dir_var] < 45, tracks[dir_var] >= 315), 'shear_cat'] ='forward'
tracks.loc[np.logical_and(tracks[dir_var] < 135, tracks[dir_var] >= 45), 'shear_cat'] ='right'
tracks.loc[np.logical_and(tracks[dir_var] < 225, tracks[dir_var] >= 135), 'shear_cat'] ='reverse'
tracks.loc[np.logical_and(tracks[dir_var] < 315, tracks[dir_var] >= 225), 'shear_cat'] ='left'

tracks.loc[tracks[speed_var] <= 0.0015 , 'shear_cat'] = 'weak'




# count_shear_cats= tracks.groupby('shear_cat').size()

# count_shear_cats.plot.bar()

    # plt.bar(nr_box.index.values, nr_box.values, bottom= nr_box_prev, label= box_name) #, width=1)

fig= plt.figure(fignr, figsize= (9,6) )
fignr+=1
plt.clf()
ax= fig.add_subplot()


nr_box_prev= 0
for box_name, box_loc in boxlist:
    count_shear_box = tracks.where((tracks.lon <= box_loc[1]) & (tracks.lon > box_loc[0]) &
                          (tracks.lat <= box_loc[2]) & (tracks.lat > box_loc[3]) ).dropna().groupby('shear_cat').size()
    #                       ).dropna().groupby(tracks.time.dt.year).size().reindex(np.arange(start_year, end_year+1), fill_value=0)/24

    # count_shear_box.plot.bar()

    plt.bar(count_shear_box.index.values, count_shear_box.values, bottom= nr_box_prev, label= box_name) #, width=1)

    nr_box_prev += count_shear_box


#the rest - still includes zeros right now
count_shear_cats= tracks.groupby('shear_cat').size()
plt.bar(count_shear_cats.index.values, count_shear_cats.values- nr_box_prev, bottom= nr_box_prev, label= 'Rest') #, width=1)


handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1])
   
if save:
    savefile= savedir+ f'{hem}_shear_distr_box'
    for vc in var_char_ind:
        importvar, importchar, excl_thresh, excl_type = vc
        savefile += '_'+ importvar + '-'+ str(int(excl_thresh))     
    savefile += '.png'
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight', dpi= 150)   




print(f'Time passed: {round(time.perf_counter() - start,2)} seconds')