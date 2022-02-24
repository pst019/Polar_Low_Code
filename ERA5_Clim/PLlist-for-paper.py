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


import numpy as np
from f_useful import *
# from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs
import time

start= time.perf_counter()
#


weak_shear_crit = 10 #1E-4 - criteria for the weak shear class

shear_mean_dist= 250

hem= 'SH'
PL_timestep= True  #add PL timestep as parameter to list

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


"""add whether a time step is a PL time step"""
if PL_timestep== True:
    tracks['PL_time_step']= 0

    for vc in var_char_ind:
        importvar, importchar, excl_thresh, excl_type = vc
        if importvar == 'diameter': continue
    
        if excl_type == 'larger':
            tracks.loc[tracks[importvar] > excl_thresh , 'PL_time_step'] += .1
        elif excl_type== 'smaller':
            tracks.loc[tracks[importvar] < excl_thresh , 'PL_time_step'] += .1
        else: print('wrong excl_type')

tracks.loc[tracks['PL_time_step'] > 0.25, 'PL_time_step'] = 1
tracks.loc[tracks['PL_time_step'] <= 0.25, 'PL_time_step'] = 0

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

tracks.loc[tracks[speed_var] <= weak_shear_crit  *1E-4 , 'shear_cat'] = 'weak'





tracks= tracks.set_index(['ID', 'step'])

tracks['vortex_diameter']= np.sqrt(tracks['area']) *4/np.pi
tracks= tracks.rename(columns={'vo': 'rel_vort_850_smth'})


tracks= tracks.drop(columns=['Unnamed: 0', 'area', 'land_dist_old', 'vortex_type', 'shear_angle_925-500_mean500',
       'shear_strength_925-500_mean500', 'near_land', 
       'theta_diff_500-925_mean250' ])

if hem == 'SH':
    tracks= tracks.drop(columns=['U_trop_polew_False250', 'U_surface_max250'])

tracks= tracks.replace(np.nan, 'nan', regex=True) #empty space for the theta500-thetaSST
tracks= tracks.replace(-1000, 'nan', regex=True) #for the SST-T500

tracks= tracks[['lon', 'lat', 'time', 'PL_time_step', 'shear_cat',
                  'rel_vort_850_smth' , 'vortex_diameter', 'theta_trop_mean250', 'theta_diff_500-sst_mean250',
                   'shear_angle_925-500_mean250', 'shear_strength_925-500_mean250',
                   'slp', 'SST-T_500_mean250', 'land_dist']]
    

print('write:', Mediadir+f"/data/ERA5_Clim/track_lists/PL_list_for_paper/{hem}")
tracks.to_csv(Mediadir+f"/data/ERA5_Clim/track_lists/PL_list_for_paper/{hem}.csv")



print(f'Time passed: {round(time.perf_counter() - start,2)} seconds')