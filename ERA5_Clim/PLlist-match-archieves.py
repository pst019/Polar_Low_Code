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
import time

start= time.perf_counter()
#


weak_shear_crit = 10 #1E-4 - criteria for the weak shear class

shear_mean_dist= 250

PL_timestep= True  #add PL timestep as parameter to list


version= ""
start_year= 1979
end_year= 2020
terrain_thresh= 50
terrain_dist = 2
more_info= True
post_process= True
split, exceed_large_to_small= False, 40

durlim= 6


ending= 'atl-pac'
SH_ending= 'SH180' #merging tracks along 180E

durlim= 6
matchdist= 150

# time_diff= False
time_diff= 24 #does not include the matched track if the matched track starts more than time_diff earlier or ends more than time_diff later



time_start= False #in this case the first and last index of alltracks will be used
# time_start= '2008-1-1' #the vorticity and slp is needed in this time
# time_end= '2009-1-1'

# PL_name_list= ['Noer', 'Rojo', 'Smirnova', 'Yanase', 'Golubkin', 'Verezemskaya']
PL_name_list= ['Verezemskaya']


var_char_ind = [['theta_diff_500-sst_mean250', 'min', 11.0, 'smaller'],
            ['theta_trop_mean250', 'mean', 300.8, 'smaller'],
            ['vo', 'max', 20, 'larger'],
            ['diameter', 'max', 430, 'smaller']
            ]    
 
csv_name= 'merged_tracks_'+version + ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
if more_info: csv_name += '_moreinfo'
if post_process: csv_name += '_post-process'
if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
# if split: csv_name += f'_split{exceed_large_to_small}'
print(csv_name)


for PL_name in PL_name_list:
    """import PL list"""
    # for PL_name in PL_name_list:
    PLfile= Mediadir+"/data/ERA5_Clim/track_lists/matchPLlist2/matchPLlist_"+ PL_name +f"_dist{matchdist}"
    if time_diff and PL_name != 'Yanase': PLfile += f'_timediff{time_diff}'
    if PL_name == 'Verezemskaya': PLfile += f'_{csv_name.replace("atl-pac", SH_ending)}'
    else: PLfile += '_'+ csv_name
    
    print(PLfile)
    PLs=pd.read_csv(PLfile +'.csv').set_index(['ID', 'step'])
    PLs.time= pd.to_datetime(PLs.time)
    
    if PL_name == 'Verezemskaya': PLs_params= pd.read_csv(PLfile+ '_params_corrected_2.csv').set_index(['ID', 'step'])
    else: PLs_params= pd.read_csv(PLfile+ '_params.csv').set_index(['ID', 'step'])
    
    PLs= PLs.merge(PLs_params, left_index=True, right_index=True)
    
    
    PLs['vo'] *= 100
    
    
    # PLs= PLs.reset_index()
    
    
    
    
    # """make the shear categories"""
    # dir_var= f'shear_angle_925-500_mean{shear_mean_dist}'
    # speed_var= f'shear_strength_925-500_mean{shear_mean_dist}'
    
    # PLs['shear_cat']= 0
    # # PLs.loc[np.logical_and(PLs[dir_var] < 3*np.pi/4, PLs[dir_var] >= 1*np.pi/4), 'shear_cat'] ='forward'
    # # PLs.loc[np.logical_and(PLs[dir_var] < 5*np.pi/4, PLs[dir_var] >= 3*np.pi/4), 'shear_cat'] ='right'
    # # PLs.loc[np.logical_and(PLs[dir_var] < 7*np.pi/4, PLs[dir_var] >= 5*np.pi/4), 'shear_cat'] ='reverse'
    # # PLs.loc[np.logical_or(PLs[dir_var] < 1*np.pi/4, PLs[dir_var] >= 7*np.pi/4), 'shear_cat'] ='left'
    # PLs.loc[np.logical_or(PLs[dir_var] < 45, PLs[dir_var] >= 315), 'shear_cat'] ='forward'
    # PLs.loc[np.logical_and(PLs[dir_var] < 135, PLs[dir_var] >= 45), 'shear_cat'] ='right'
    # PLs.loc[np.logical_and(PLs[dir_var] < 225, PLs[dir_var] >= 135), 'shear_cat'] ='reverse'
    # PLs.loc[np.logical_and(PLs[dir_var] < 315, PLs[dir_var] >= 225), 'shear_cat'] ='left'
    
    # PLs.loc[PLs[speed_var] <= weak_shear_crit  *1E-4 , 'shear_cat'] = 'weak'
    
    
    
    """ add wheather it is a PL track"""
    PLs_ind= PLs.groupby(['ID']).max()
    PLs_ind= PLs_ind.rename(columns={'vo': 'vo_max', 'area': 'area_max'})
    
    PLs_ind['diameter_max']= np.sqrt(PLs_ind['area_max']) *4/np.pi
    PLs_ind['vo_max'] *= 100
    # PLs_ind['prop_dist'] = distance( (PLs.groupby(['ID']).first()['lat'], PLs.groupby(['ID']).first()['lon']),
    #           (PLs.groupby(['ID']).last()['lat'], PLs.groupby(['ID']).last()['lon']) )
    
    
    PLs= PLs.replace(-1000, 1000, regex=True) #for the SST-T500
    
    
    for vc in var_char_ind:
        importvar, importchar, excl_thresh, excl_type = vc
        varname= importvar+'_'+importchar
        
        if varname not in PLs_ind.columns:
            
            if importchar== 'max':
                PLs_ind_param=  PLs[importvar].groupby(['ID']).max().rename( importvar+'_'+importchar)
            if importchar== 'min':
                PLs_ind_param=  PLs[importvar].groupby(['ID']).min().rename( importvar+'_'+importchar)
            if importchar== 'mean':
                PLs_ind_param=  PLs[importvar].groupby(['ID']).mean().rename( importvar+'_'+importchar)
            
            PLs_ind= PLs_ind.merge(PLs_ind_param, left_index=True, right_index=True)
    
    
    PLs_ind_crit = PLs_ind
    for vc in var_char_ind:
        importvar, importchar, excl_thresh, excl_type = vc
        varname= importvar+'_'+importchar
    
        if excl_type == 'larger':
            PLs_ind_crit= PLs_ind_crit[PLs_ind_crit[varname] > excl_thresh] 
        elif excl_type== 'smaller':
            PLs_ind_crit= PLs_ind_crit[PLs_ind_crit[varname] < excl_thresh] 
        else: print('wrong excl_type')
            
    
    sel= PLs_ind_crit.index.values
    
    PLs['PL_track']= 0
    
    for s in sel:
        PLs.loc[s, 'PL_track']= 1
    
    
    """add whether a time step is a PL time step"""
    if PL_timestep== True:
        PLs['PL_time_step']= 0
    
        for vc in var_char_ind:
            importvar, importchar, excl_thresh, excl_type = vc
            if importvar == 'diameter': continue
        
            if excl_type == 'larger':
                PLs.loc[PLs[importvar] > excl_thresh , 'PL_time_step'] += .1
            elif excl_type== 'smaller':
                PLs.loc[PLs[importvar] < excl_thresh , 'PL_time_step'] += .1
            else: print('wrong excl_type')
    
    PLs.loc[PLs['PL_time_step'] > 0.25, 'PL_time_step'] = 1
    PLs.loc[PLs['PL_time_step'] <= 0.25, 'PL_time_step'] = 0
    
    
    
    
    PLs['vortex_diameter']= np.sqrt(PLs['area']) *4/np.pi
    PLs= PLs.rename(columns={'vo': 'rel_vort_850_smth'})
    if PL_name != 'Yanase':
        PLs= PLs.rename(columns={ 'U_trop_polew_False250':  'U_trop_polew'})
    
    PLs= PLs.drop(columns=['Unnamed: 0', 'area', 'vortex_type', 'shear_angle_925-500_mean500',
            'shear_strength_925-500_mean500', 'near_land' ])
    
    # if hem == 'SH':
    #     PLs= PLs.drop(columns=['U_trop_polew_False250', 'U_surface_max250'])
    
    PLs= PLs.replace(np.nan, 'nan', regex=True) #empty space for the theta500-thetaSST
    PLs= PLs.replace(1000, 'nan', regex=True) #for the SST-T500
    
    # if PL_name not in ['Yanase', 'Verezemskaya']:
    PLs= PLs[['lon', 'lat', 'time', 'PL_track', 'PL_time_step', 'PLlist_match', 'PLlist_match_2', 'PLlist_match_3', 
                  'rel_vort_850_smth' , 'vortex_diameter', 'theta_trop_mean250', 'theta_diff_500-sst_mean250',
                  'slp', 'SST-T_500_mean250', 'U_surface_max250', 'U_trop_polew', 'land_dist', 'theta_diff_500-925_mean250']]
    
    print('write:', Mediadir+f"/data/ERA5_Clim/track_lists/Detected_systems_from_PL_archives/{PL_name}")
    PLs.to_csv(Mediadir+f"/data/ERA5_Clim/track_lists/Detected_systems_from_PL_archives/{PL_name}.csv")
    
    
    
    # print(f'Time passed: {round(time.perf_counter() - start,2)} seconds')