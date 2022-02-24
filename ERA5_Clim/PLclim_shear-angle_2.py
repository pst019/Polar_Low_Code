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
params = {'axes.labelsize': 15,
          # 'axes.titlesize': 16,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15
          #, 'legend.fontsize': 13
          }
plt.rcParams.update(params)
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
# save= False
savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/PL_shear_distribution/'
fignr= 11


display_fraction= True #calculate the fraction of the total activity that each category is responsible for
weak_shear_crit = 10 #1E-4 - criteria for the weak shear class

shear_mean_dist= 250

hem= 'NH'
PL_timestep= True  # True #only the time steps that satisfy all PL crit

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
    boxlist= [ ['Nordic Seas', [-15, 70, 80, 55] ] #W, E, N, S
        , ['Irminger Sea', [-45, -15, 75, 50] ]
        , ['Labrador Sea', [-70, -45, 75, 50] ]
        , ['Gulf of Alaska', [200, 235, 65, 40] ]       
        , ['Bering Sea', [160, 200, 65, 40] ]       
        , ['Sea of Okhotsk', [142, 160, 65, 40] ]       
        , ['Sea of Japan',  [127, 142, 50, 35] ]
        ]

elif hem == 'SH':
    boxlist= [
    ['SE Pacific', [-140, -60, -38, -70] ] #Amundsen Sea #W, E, N, S
    , ['SW Pacific', [150, 220, -38, -70] ]
    , ['Indian Ocean', [20, 150, -38, -70] ]
    , ['Atlantic', [-60, 20, -38, -70] ]
    ]    
    # boxlist= [ ['Amundsen Sea', [-140, -50, -45, -70] ] #W, E, N, S
    # , ['New Zealand', [150, 210, -45, -70] ]
    # , ['Indian Ocean', [40, 120, -45, -70] ]
    # ]

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


"""only steps that satisfy PL-IC"""
if PL_timestep== True:
    for vc in var_char_ind:
        importvar, importchar, excl_thresh, excl_type = vc
        if importvar == 'diameter': continue
        if importvar == 'vo':  tracks['vo'] *= 100
    
        if excl_type == 'larger':
            tracks= tracks[tracks[importvar] > excl_thresh] 
        elif excl_type== 'smaller':
            tracks= tracks[tracks[importvar] < excl_thresh] 
        else: print('wrong excl_type')


tracks= tracks.dropna(axis= 'columns') #drops some columns where parameters are only calculated for some years




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

tracks.loc[tracks[speed_var] <= weak_shear_crit  *1E-4 , 'shear_cat'] = 'weak'



# count_shear_cats= tracks.groupby('shear_cat').size()

# count_shear_cats.plot.bar()

    # plt.bar(nr_box.index.values, nr_box.values, bottom= nr_box_prev, label= box_name) #, width=1)

fig= plt.figure(fignr, figsize= (9,6) )
fignr+=1
plt.clf()
ax= fig.add_subplot()


#calculate the rest - still includes zeros right now
count_shear_cats= tracks.groupby('shear_cat').size()/(end_year - start_year +1)/24
# if display_fraction: count_shear_cats /= (count_shear_cats.sum()/100)

nr_box_prev= 0
for box_name, box_loc in boxlist:
    count_shear_box = tracks.where((tracks.lon <= box_loc[1]) & (tracks.lon > box_loc[0]) &
                          (tracks.lat <= box_loc[2]) & (tracks.lat > box_loc[3]) ).dropna().groupby('shear_cat').size()/(end_year - start_year +1)/24
    #                       ).dropna().groupby(tracks.time.dt.year).size().reindex(np.arange(start_year, end_year+1), fill_value=0)/24

    # count_shear_box.plot.bar()
    

    if display_fraction:
        ax.bar(count_shear_box.index.values, count_shear_box.values/(count_shear_cats.sum()/100), bottom= nr_box_prev/(count_shear_cats.sum()/100), label= box_name) #, width=1)
    else:
        ax.bar(count_shear_box.index.values, count_shear_box.values, bottom= nr_box_prev, label= box_name) #, width=1)

    nr_box_prev += count_shear_box


    """the individual plot for each basin"""
    sfig= plt.figure(fignr, figsize= (5,3.5) )  
    fignr+=1
    plt.clf()
    sax= sfig.add_subplot()


    sax.bar(count_shear_box.index.values, count_shear_box.values/(count_shear_box.sum()/100)) #, width=1)

    # handles, labels = sax.get_legend_handles_labels()
    # sax.legend(handles[::-1], labels[::-1])
    if display_fraction:
        sax.set_ylabel('Fraction of activity [%]')
        sax.set_ylim(0,42)        
    
    else: sax.set_ylabel('Time of activity [d]')


    if save:
        savefile= savedir+ 'shear_distr_'+ box_name.replace(' ', '_')
        if weak_shear_crit != 15: savefile += f'_weak-shear-crit{weak_shear_crit}'

        if PL_timestep== True: savefile += '_PL-ts'
        if display_fraction: savefile += '_disp-frac'
            
        for vc in var_char_ind:
            importvar, importchar, excl_thresh, excl_type = vc
            savefile += '_'+ importvar + '-'+ str(int(excl_thresh))
                
        savefile += '.png'
        print(savefile)
        plt.savefig(savefile , bbox_inches='tight', dpi= 150)  


# plot the rest which is calculated befor the loop
if display_fraction:
    ax.bar(count_shear_cats.index.values, (count_shear_cats.values- nr_box_prev)/(count_shear_cats.sum()/100), bottom= nr_box_prev/(count_shear_cats.sum()/100), label= 'Rest') #, width=1)
else:
    ax.bar(count_shear_cats.index.values, count_shear_cats.values- nr_box_prev, bottom= nr_box_prev, label= 'Rest') #, width=1)


handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1])
if display_fraction: ax.set_ylabel('Fraction of polar-low activity [%]')
else: ax.set_ylabel('Time of polar-low activity [d]')

   
if save:
    savefile= savedir+ f'{hem}_shear_distr_box'
    if weak_shear_crit != 15: savefile += f'_weak-shear-crit{weak_shear_crit}'
    if shear_mean_dist != 500: savefile += f'_shear-mean{shear_mean_dist}'
    if PL_timestep== True: savefile += '_PL-ts'
    if display_fraction: savefile += '_disp-frac'

    for vc in var_char_ind:
        importvar, importchar, excl_thresh, excl_type = vc
        savefile += '_'+ importvar + '-'+ str(int(excl_thresh))     
    savefile += '.png'
    print(savefile)
    fig.savefig(savefile , bbox_inches='tight', dpi= 150)   




print(f'Time passed: {round(time.perf_counter() - start,2)} seconds')