#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 16:20:33 2021

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/PatsOrange/'

homedir=Mediadir+'home/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

from f_useful import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import scipy.signal as signal
import matplotlib as mpl
from cycler import cycler

# save=True
save= False
savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/PLvsTRACKS_distr_2/'

fignr = 4

version= "" #"fram_run3_"
ending= 'atl-pac'

start_year= 1979
end_year= 2020

durlim= 6
matchdist= 150

# time_diff= False
time_diff= 24 #does not include the matched track if the matched track starts more than time_diff earlier or ends more than time_diff later

more_info= True
post_process= True
# split, exceed_large_to_small= True, 40
split= False # is included in post_process

terrain_dist, terrain_thresh= 2, 50
# terrain_thresh= False



percentile= 10
PLs_larger= False

# PL_name= 'Rojo'
# PL_name = 'Noer'
# PL_name = 'Smirnova'
# PL_name = 'Yanase' #it has only one location per PL
# PL_name_list= [PL_name]

PL_name_list= ['Noer', 'Rojo', 'Smirnova', 'Yanase', 'Golubkin', 'Verezemskaya']
# PL_name_list= ['Smirnova', 'Yanase']
# PL_name_list= ['Noer', 'Rojo']
# PL_name_list= ['Noer', 'Yanase', 'Golubkin']


var_char_ind = [
                # ['theta_diff_500-sst_mean250', 'min'],
                ['SST-T_500_mean250', 'max'],

                ['theta_trop_mean250', 'min'],
                # ['land_dist', 'max']
                ]


#for the climatology
excl_crit_list= [
                ['vo_max', 20, 'larger'],
                # ['vo_tendency', 0.5, 'larger'],                
                ['diameter_max', 430, 'smaller'], #429
                # ['theta_diff_500-sst_mean250_min', 10.5, 'smaller'],
                # ['theta_diff_500-sst_mean250_min', 13, 'smaller'],

                ['SST-T_500_mean250_max', 41.5, 'larger'],
                # ['SST-T_500_mean250_max', 38.4, 'larger'],

                ['theta_trop_mean250_min', 298.4, 'smaller'],
                # ['land_dist_max', 140, 'larger']
                ]


#for strong cases
# excl_crit_list= [
#                 ['vo_max', 20, 'larger'],
#                 ['diameter_max', 430, 'smaller'], #429
#                 ['SST-T_500_mean250_max', 41.5, 'larger'],
#                 ['theta_trop_mean250_min', 298.4, 'smaller'],
#                 ]





bins, plotname= [], ''


SH_ending= 'SH180' #merging tracks along 180E
# SH_ending= 'SH' #the old one



# plotvar, plotchar, bins, plotname='U_trop_polew', 'mean', np.arange(0, 100.1, 1), r'U$_{trop, polew.}$ [m s$^{-1}$]'
# plotvar, plotchar, bins, plotname='theta_trop_mean250', 'min' , np.arange(265, 380.1, 1), r'$\theta_{trop}$ [K]'
# plotvar, plotchar, bins, plotname='theta_trop_mean250', 'max' , np.arange(265, 380.1, 1), r'$\theta_{trop}$ [K]'

# plotvar, plotchar, bins, plotname='theta_trop_mean100', 'min' , np.arange(265, 380.1, 1), r'$\theta_{trop}$ [K]'
# plotvar, plotchar, bins, plotname='theta_trop_mean500', 'min' , np.arange(265, 380.1, 1), r'$\theta_{trop}$ [K]'
# plotvar, plotchar, bins, plotname='theta_trop_mean250', 'min' , np.arange(265, 380.1, .2), r'$\theta_{trop}$ [K]'

# plotvar, plotchar='theta_trop_mean500', 'mean' #'max'

# plotvar, plotchar, bins, plotname='theta_diff_500-sst_mean250', 'min', np.arange(-8, 42.1, 1), r'$\theta_{500hPa} - \theta_{SST}$ [K]'  #'max'
## plotvar, plotchar, bins, plotname='theta_diff_500-sst_mean250', 'mean', np.arange(-8, 42.1, 1), r'$\theta_{500hPa} - \theta_{SST}$ [K]'  #'max'
## plotvar, plotchar, bins, plotname='theta_diff_500-sst_mean100', 'min', np.arange(-5, 40.1, 0.5), r'$\theta_{500hPa} - \theta_{SST}$ [K]'  #'max'
## plotvar, plotchar, bins, plotname='theta_diff_500-sst_min250', 'min', np.arange(-10, 40.1, 0.5), r'$\theta_{500hPa} - \theta_{SST}$ [K]'  #'max'

# plotvar, plotchar, bins, plotname, PLs_larger='SST-T_500_mean250', 'max', np.arange(15, 55.1, 1), r'SST - T$_{500hPa}$ [K]' , True


# plotvar, plotchar='theta_diff_500-skt_mean250', 'min' #'max'
# plotvar, plotchar, bins, plotname='theta_diff_500-925_mean250', 'min',  np.arange(0, 35.1, 0.5), r'$\theta_{500hPa} - \theta_{925hPa}$ [K]'

# plotvar, plotchar, PLs_larger, bins, plotname='U_surface_max250', 'max', True, np.arange(0, 30.1, 0.5), 'U$_{10m}$ [m s$^{-1}$]'
# plotvar, plotchar, PLs_larger, bins, plotname='vo_max', False, True, np.arange(0, 90.1, 0.5), 'Relative vorticity [10$^{-5}$s$^{-1}$]'
# plotvar, plotchar, PLs_larger, plotname= 'vo_tendency', False, True , r'Vorticity tend [10^{-5}$s$^{-1}$/h]'

# plotvar, plotchar, bins, plotname='slp', False, np.arange(945, 1030.1, 1), 'Sea-level pressure [hPa]'
# plotvar, plotchar='vortex_type', False
# plotvar, plotchar, bins, plotname='diameter_max', False, np.arange(0, 901, 25), 'Diameter [km]'
plotvar, plotchar, bins, PLs_larger, plotname= 'duration', False, np.arange(0, 90.1, 1), True, 'Life time [h]'

# plotvar, plotchar, PLs_larger, bins= 'near_land', False, False, np.arange(0, 1.01, 0.025)
# plotvar, plotchar, bins='month', False, np.arange(0.5, 13, 1)

# plotvar, plotchar, plotname = 'prop_dist', False, 'Propagation distance [km]'
# plotvar, plotchar = 'prop_speed', False
# plotvar, plotchar, plotname, bins, PLs_larger= 'land_dist', 'max', 'Distance to land [km]', np.arange(0, 701, 10), True #only accurate for short distances


# plotvar, plotchar, PLs_larger= 'shear_strength_925-500_mean500', 'max', True



"""the plot variable"""
if plotchar == False: plotvarchar= plotvar
else:
    plotvarchar= plotvar+'_'+plotchar
    if [plotvar, plotchar] not in var_char_ind: var_char_ind += [[plotvar, plotchar]]


#excludes the plotted criteria from the list of exclusion criteria (if included)
# if plotvarchar in [e[0] for e in excl_crit_list]:
excl_crit_list = [e for e in excl_crit_list if e[0] != plotvarchar]

if plotname == '': plotname= plotvarchar

"""import lists"""
print('two different track lists')
csv_name= 'merged_tracks_'+version + ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
if more_info: csv_name += '_moreinfo'
if post_process: csv_name += '_post-process'
if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
# if split: csv_name += f'_split{exceed_large_to_small}'
print(csv_name)

tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks2/"+csv_name+'.csv')

tracks_SH= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks2/"+csv_name.replace("atl-pac", SH_ending)+'.csv')

time_start= False #in this case the first and last index of alltracks will be used
time_start= '2008-1-1' #the vorticity and slp is needed in this time
time_end= '2009-1-1'

time_start_SH= '2004-1-1' #the vorticity and slp is needed in this time
time_end_SH= '2005-1-1'

if time_start != False:
    tracks.time= pd.to_datetime(tracks.time)
    tracks= tracks[(tracks['time'] >= time_start) & (tracks['time'] < time_end)]
    tracks_SH.time= pd.to_datetime(tracks_SH.time)
    tracks_SH= tracks_SH[(tracks_SH['time'] >= time_start_SH) & (tracks_SH['time'] < time_end_SH)]


# print('now it is only compared against atlantic tracks')

systems= [tracks, tracks_SH] #list that gets tracks and PLs
systems_name_list= ['tracks NH', 'tracks SH'] + PL_name_list
# systems_name_list= ['tracks NH'] + PL_name_list

systems_ind= []

for PL_name in PL_name_list:
    PLfile= Mediadir+"/data/ERA5_Clim/track_lists/matchPLlist2/matchPLlist_"+ PL_name +f"_dist{matchdist}"
    if time_diff and PL_name != 'Yanase': PLfile += f'_timediff{time_diff}'
    if PL_name == 'Verezemskaya': PLfile += f'_{csv_name.replace("atl-pac", SH_ending)}.csv'
    else: PLfile += f'_{csv_name}.csv'
    
    PLs=pd.read_csv(PLfile)
    PLs.time= pd.to_datetime(PLs.time)
    systems += [PLs]


# csv_name= 'merged_tracks_'+version+'_'+ending+f'_{start_year}-{end_year}'
# if durlim: csv_name += 'durlim'+str(durlim)
# if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
# # if split: csv_name += f'_split{exceed_large_to_small}'
# print(csv_name)

"""import the parameters"""
for system, system_name in zip(systems, systems_name_list):
    # print(system_name)
    system_ind= system.groupby(['ID']).max()
    system_ind= system_ind.rename(columns={'step':'duration', 'vo': 'vo_max', 'area': 'area_max'})
    

    """get individual params for systems (tracks + PLs) - calculate the system_ind variable """
    if time_start != False and system_name == 'tracks NH':
        system_params= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks_tomerge/"+csv_name+f'_params_{time_start}-{time_end}.csv').set_index(['ID', 'step'])
        #have to exclude the ones that are not in the tracks
    elif time_start_SH != False and system_name == 'tracks SH':
        system_params= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks_tomerge/"+csv_name.replace("atl-pac", SH_ending)+f'_params_{time_start_SH}-{time_end_SH}.csv').set_index(['ID', 'step'])
    elif system_name== 'tracks NH':
        system_params= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks2/"+csv_name+'_params.csv').set_index(['ID', 'step'])
    
    else: #for the PLs
        PLfile_params = Mediadir+"/data/ERA5_Clim/track_lists/matchPLlist2/matchPLlist_"+system_name+f'_dist{matchdist}'
        if time_diff and system_name != 'Yanase': PLfile_params += f'_timediff{time_diff}'
        if system_name== 'Verezemskaya':
            system_params= pd.read_csv(PLfile_params+ f'_{csv_name.replace("atl-pac", SH_ending)}_params.csv').set_index(['ID', 'step'])
        else:            
            system_params= pd.read_csv(PLfile_params+ f'_{csv_name}_params.csv').set_index(['ID', 'step'])
    
    if 'U_trop_polew_False250' in system_params.columns:
        system_params= system_params.rename(columns={'U_trop_polew_False250':'U_trop_polew'})
    
    system_params[system_params['theta_diff_500-sst_mean250'] == -1000] = np.nan

    for vc in var_char_ind:
        importvar, importchar= vc
        if importchar== 'max':
            system_ind_param=  system_params[importvar].groupby(['ID']).max().rename( importvar+'_'+importchar)
        if importchar== 'min':
            system_ind_param=  system_params[importvar].groupby(['ID']).min().rename( importvar+'_'+importchar)
        if importchar== 'mean':
            system_ind_param=  system_params[importvar].groupby(['ID']).mean().rename( importvar+'_'+importchar)
        if importchar== 'median':
            system_ind_param=  system_params[importvar].groupby(['ID']).median().rename( importvar+'_'+importchar)
        if importchar== 'initial':
            system_ind_param=  system_params[importvar].groupby(['ID']).first().rename( importvar+'_'+importchar)
        if importchar== 'final':
            system_ind_param=  system_params[importvar].groupby(['ID']).last().rename( importvar+'_'+importchar)  
            
        system_ind= system_ind.merge(system_ind_param, left_index=True, right_index=True)


    """create additional variables"""
    # if plotvarchar == 'diameter_max':
    system_ind['diameter_max']= np.sqrt(system_ind['area_max']) *4/np.pi
    system_ind['vo_max'] *= 100

    if plotvar== 'month':
        system[plotvar]= pd.to_datetime(system['time']).dt.month
        system_ind[plotvar]= system.groupby(['ID']).median()[plotvar]

    if plotvar== 'near_land':
        system_ind[plotvar]= system.groupby(['ID']).mean()[plotvar]

    if plotvar== 'slp':
        system_ind[plotvar]= system.groupby(['ID']).mean()[plotvar]

    if plotvar == 'prop_dist':
        system_ind[plotvar] = distance( (system.groupby(['ID']).first()['lat'], system.groupby(['ID']).first()['lon']),
                  (system.groupby(['ID']).last()['lat'], system.groupby(['ID']).last()['lon']) )

        # system_ind[plotvar] = (110* np.sqrt( (system.groupby(['ID']).first()['lat'] - system.groupby(['ID']).last()['lat'])**2 +
        #                     (np.cos(np.deg2rad(system.groupby(['ID']).first()['lat'])) *(system.groupby(['ID']).first()['lon']- system.groupby(['ID']).last()['lon']) )**2 ) )



    if plotvar == 'prop_speed':
        print(system)
        system_ind[plotvar] = distance( (system.groupby(['ID']).first()['lat'], system.groupby(['ID']).first()['lon']),
                 (system.groupby(['ID']).last()['lat'], system.groupby(['ID']).last()['lon']) )  /(
                 (system.groupby(['ID']).last()['time'] - system.groupby(['ID']).first()['time'])/ np.timedelta64(1,'h') ) #in m/s


    if plotvar == 'vo_tendency' or 'vo_tendency' in [e[0] for e in excl_crit_list]:
        system['vo_tendency']= np.nan
        for index in system_ind.index:
            track_now=system[system.ID == index]
            
            #just when tracks are time reduced:
            if len(track_now) >= 6:
                if len(track_now)>= 11:  window= 11
                else: window= len(track_now) -(len(track_now)+1)%2 #to get the largest uneven number smaller than the length of vorticity
            
                # track_now= track_now.assign(vo_tendency= signal.savgol_filter(track_now['vo'].values, window_length=window, polyorder=2, deriv= 1) )
                vo_tend= signal.savgol_filter(track_now['vo'].values*100, window_length=window, polyorder=2, deriv= 1)
    
                system.loc[system.ID == index, 'vo_tendency']= vo_tend
            
        system_ind['vo_tendency']=  system.set_index(['ID', 'step'])['vo_tendency'].groupby(['ID']).max()


    systems_ind += [system_ind]







"""tracks after crit"""
# tracks_ind_crit = systems_ind[0]
# tracks_ind_SH_crit = systems_ind[1]



# for excl_crit_dict in excl_crit_list:
#     excl_crit, excl_thresh, excl_type= excl_crit_dict

#     if excl_type == 'smaller':
#         tracks_ind_crit= tracks_ind_crit[tracks_ind_crit[excl_crit] < excl_thresh]
#         tracks_ind_SH_crit= tracks_ind_crit[tracks_ind_SH_crit[excl_crit] < excl_thresh]

#     elif excl_type == 'larger':
#         tracks_ind_crit= tracks_ind_crit[tracks_ind_crit[excl_crit] > excl_thresh]
#         tracks_ind_SH_crit= tracks_ind_crit[tracks_ind_SH_crit[excl_crit] > excl_thresh]
        
#     else: print('excl_type is wrong')


# tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['theta_diff_500-925_mean250_min'] < 20] #conservative Yanase
# tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['theta_diff_500-sst_mean250_min'] < 12.5] #conservative Yanase

# # tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['theta_diff_500-sst_min250_min'] < 8.] #all are stricter beside Smironova

# # tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['theta_trop_mean250_min'] < 300] #conservative - Yan /Smir/Noer

# tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['vo_max'] > 22] #in between
# tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['diameter_max'] < 450] #approx middle


# systems_ind = systems_ind[:2]+ [tracks_ind_crit, tracks_ind_SH_crit]+ systems_ind[2:]
# systems_name_list =systems_name_list[:2] + ['tracks NH after crit', 'tracks SH after crit'] + systems_name_list[2:]
# print('Tracks remaining with criteria:', np.round(100*len(tracks_ind_crit)/len(systems_ind[0]), 1), '%' , len(tracks_ind_crit), 'of', len(systems_ind[0]) ,'\n')



# for system, system_name in zip(systems_ind[4:], systems_name_list[4:]): 
#     sys_crit= system[#(system['theta_diff_500-sst_mean250_min'] < 12.5) &
#                      # (system['theta_diff_500-sst_min250_min'] < 8.) &
#                      # (system['theta_trop_mean250_min'] < 300) &
#                       (system['vo_max'] > 22) &
#                      (system['diameter_max'] < 450)]
#     print(system_name, np.round(100*len(sys_crit)/len(system)),'% (', len(sys_crit), 'of', len(system),') remain with criteria')

# print('\n')


for system, system_name in zip(systems_ind, systems_name_list): 
    system_crit= system
    for excl_crit_dict in excl_crit_list:
        excl_crit, excl_thresh, excl_type= excl_crit_dict
        if excl_type == 'smaller':
            system_crit= system_crit[system_crit[excl_crit] < excl_thresh]   
        elif excl_type == 'larger':
            system_crit= system_crit[system_crit[excl_crit] > excl_thresh]
    
    print(system_name, np.round(100*len(system_crit)/len(system), 1),'% (', len(system_crit), 'of', len(system),') remain with criteria')

    if system_name =='tracks NH': tracks_ind_crit= system_crit
    elif system_name =='tracks SH': tracks_ind_SH_crit= system_crit
        
print('\n')

systems_ind = systems_ind[:2]+ [tracks_ind_crit, tracks_ind_SH_crit]+ systems_ind[2:]
systems_name_list =systems_name_list[:2] + ['tracks NH after crit', 'tracks SH after crit'] + systems_name_list[2:]
# systems_ind = systems_ind[:1]+ [tracks_ind_crit]+ systems_ind[1:]
# systems_name_list =systems_name_list[:1] + ['tracks NH after crit'] + systems_name_list[1:]




system_var_list= [system[plotvarchar] for system in systems_ind]
   
if len(bins)== 0:
    bins= np.linspace(system_var_list[0].min(), system_var_list[0].max(), 20)


"""plot histogram"""

# mpl.rcParams['axes.prop_cycle']= cycler(color= ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])
# plt.figure(fignr)
# fignr+=1
# plt.clf() 

# for system_var, system_name in zip(system_var_list, systems_name_list):
#     plt.hist(system_var, density= True, bins= bins, alpha= 0.6, label= system_name)


# plt.legend()
# plt.xlabel(plotname)
# plt.ylabel('Density')


"""plot density distribution"""
# mpl.rcParams['axes.prop_cycle'] = cycler(color='kkbgrcmy')
mpl.rcParams['axes.prop_cycle']= cycler(color= ['k', 'k', 'k', 'k', "#800000FF", "#767676FF", "#FFA319FF" ,"#8A9045FF", "#155F83FF", "#C16622FF"])
# mpl.rcParams['axes.prop_cycle']= cycler(color= ['k', 'k',"#800000FF", "#767676FF", "#FFA319FF" ,"#8A9045FF", "#155F83FF", "#C16622FF"])


plt.figure(fignr, figsize= (8, 6))
fignr+=1
plt.clf() 
# fig,ax = plt.subplots(1)

xplot= np.linspace(bins[0], bins[-1], 100)
 

crit_list, crit_list_name = [], []

for system_var, system_name in zip(system_var_list, systems_name_list):
    system_var= system_var.dropna()
    # if plotvar== 'vo_max': system_var *= 100
    
    kernel= stats.gaussian_kde(system_var)
    if system_name== 'tracks NH after crit':    p= plt.plot(xplot, kernel(xplot), '--', label= system_name, linewidth= 2 )
    elif system_name== 'tracks SH':    p= plt.plot(xplot, kernel(xplot), '-.', label= system_name, linewidth= 2 )
    elif system_name== 'tracks SH after crit':    p= plt.plot(xplot, kernel(xplot), ':', label= system_name, linewidth= 2 )

    else:    p= plt.plot(xplot, kernel(xplot), label= system_name, linewidth= 2 )

    if system_name in systems_name_list[4:]:
    # if system_name in systems_name_list[2:]:
        
        if PLs_larger:
            crit= np.percentile( system_var, percentile, interpolation='lower')
            crit_list += [crit]
            crit_list_name += [system_name]
            
            plt.scatter(crit, 0, marker='^', c= p[0].get_color(), s= 400 )
        
            #if PL value is larger than track value
            # print(PL_name, np.round(crit, 3), 'excludes: ', np.round(len(np.where(tracks_var < crit)[0])/len(tracks_var), 2))
            print(system_name, 'crit:', np.round(crit, 3), 'excludes:' , 100*np.round(len(np.where(system_var_list[0] < crit)[0])/len(system_var_list[0]), 3), '% (', len(np.where(system_var_list[0] < crit)[0]), 'of', len(system_var_list[0]), ')' )
            print('Of remaining tracks excludes:' , np.round(100*(len(np.where(system_var_list[2] < crit)[0])/len(system_var_list[2]) ), 1), '% (', len(np.where(system_var_list[2] < crit)[0]), 'of', len(system_var_list[2]), ')\n' )
    
    
        else:
            crit= np.percentile( system_var, 100- percentile, interpolation='higher')
            crit_list += [crit]
            crit_list_name += [system_name]
            
            plt.scatter(crit, 0, marker='^', c= p[0].get_color(), s= 400 )
        
            #if PL value is larger than track value
            print(system_name, 'crit:', np.round(crit, 3), 'excludes:' , 100*np.round(len(np.where(system_var_list[0] > crit)[0])/len(system_var_list[0]), 3), '% (', len(np.where(system_var_list[0] > crit)[0]), 'of', len(system_var_list[0]), ')' )
            print('Of remaining tracks excludes:' , np.round(100* (len(np.where(system_var_list[2] > crit)[0])/len(system_var_list[2])), 1), '% (', len(np.where(system_var_list[2] > crit)[0]), 'of', len(system_var_list[2]), ')\n' )


str_for_latex= ''

if PLs_larger:
    crit_weak= np.min(crit_list[:-1])   
    crit_weak_name= crit_list_name[np.argmin(crit_list[:-1])]
    # crit_weak= 44

    # crit_weak= np.min(crit_list[:3] + [crit_list[4]])
    # crit_weak_name= crit_list_name[np.argmin(crit_list[:3] + [crit_list[4]])]
    
    print('Weakest criteria:', crit_weak, crit_weak_name)
    print('Exclusion of the PL lists:')    
    for system_var,  system_name in zip(system_var_list[4:], systems_name_list[4:]):
        print(system_name, np.round(100* (len(np.where(system_var < crit_weak)[0])/len(system_var)), 1), '% (', len(np.where(system_var < crit_weak)[0]), 'of', len(system_var), ')' )
        str_for_latex += str(int(np.round(100* (len(np.where(system_var < crit_weak)[0])/len(system_var)))))+'/'
    str_for_latex= str_for_latex[:-1]+' & '
    
    print('Excludes NH:' , np.round(100* (len(np.where(system_var_list[0] < crit_weak)[0])/len(system_var_list[0])), 1), '% (', len(np.where(system_var_list[0] < crit_weak)[0]), 'of', len(system_var_list[0]), ')' )
    str_for_latex += str(int(100*np.round(len(np.where(system_var_list[0] < crit_weak)[0])/len(system_var_list[0]), 2))) + '/ '  

    print('Excludes SH:' , np.round(100*(len(np.where(system_var_list[1] < crit_weak)[0])/len(system_var_list[1])), 1), '% (', len(np.where(system_var_list[1] < crit_weak)[0]), 'of', len(system_var_list[1]), ')' )
    str_for_latex += str(int(100*np.round(len(np.where(system_var_list[1] < crit_weak)[0])/len(system_var_list[1]), 2))) + ' & '  

    print('Remaining tracks NH excludes:' , np.round(100*(len(np.where(system_var_list[2] < crit_weak)[0])/len(system_var_list[2])), 1), '% (', len(np.where(system_var_list[2] < crit_weak)[0]), 'of', len(system_var_list[2]), ')' )
    str_for_latex += str(int(100*np.round(len(np.where(system_var_list[2] < crit_weak)[0])/len(system_var_list[2]), 2))) + '/ '  

    print('Remaining tracks SH excludes:' , np.round(100*(len(np.where(system_var_list[3] < crit_weak)[0])/len(system_var_list[3])), 1), '% (', len(np.where(system_var_list[3] < crit_weak)[0]), 'of', len(system_var_list[3]), ')\n' )
    str_for_latex += str(int(100*np.round(len(np.where(system_var_list[3] < crit_weak)[0])/len(system_var_list[3]), 2))) + ' \\ '  
        

else:
    crit_weak= np.max(crit_list[:-1])
    # crit_weak= 300
    crit_weak_name= crit_list_name[np.argmax(crit_list[:-1])]
    # crit_weak= np.max(crit_list[:3] + [crit_list[4]]) #excludes Yanase and Verezemskaya
    # crit_weak_name= crit_list_name[np.argmax(crit_list[:3] + [crit_list[4]])]    
    
    print('Weakest criteria:', crit_weak, crit_weak_name)
    print('Exclusion of the PL lists:')
    for system_var,  system_name in zip(system_var_list[4:], systems_name_list[4:]):
        print(system_name, np.round(100* (len(np.where(system_var > crit_weak)[0])/len(system_var)), 1), '% (', len(np.where(system_var > crit_weak)[0]), 'of', len(system_var), ')' )
        str_for_latex += str(int(np.round(100* (len(np.where(system_var > crit_weak)[0])/len(system_var)))))+'/'
    str_for_latex= str_for_latex[:-1]+' & '
    
    print('Excludes NH:' , 100*np.round(len(np.where(system_var_list[0] > crit_weak)[0])/len(system_var_list[0]), 3), '% (', len(np.where(system_var_list[0] > crit_weak)[0]), 'of', len(system_var_list[0]), ')' )
    str_for_latex += str(int(100*np.round(len(np.where(system_var_list[0] > crit_weak)[0])/len(system_var_list[0]), 2))) + '/ '  
    print('Excludes SH:' , 100*np.round(len(np.where(system_var_list[1] > crit_weak)[0])/len(system_var_list[1]), 3), '% (', len(np.where(system_var_list[1] > crit_weak)[0]), 'of', len(system_var_list[1]), ')' )
    str_for_latex += str(int(100*np.round(len(np.where(system_var_list[1] > crit_weak)[0])/len(system_var_list[1]), 2))) + ' & '  

    print('Remaining tracks NH excludes:' , 100*np.round(len(np.where(system_var_list[2] > crit_weak)[0])/len(system_var_list[2]), 3), '% (', len(np.where(system_var_list[2] > crit_weak)[0]), 'of', len(system_var_list[2]), ')' )
    str_for_latex += str(int(100*np.round(len(np.where(system_var_list[2] > crit_weak)[0])/len(system_var_list[2]), 2))) + '/ '  
    print('Remaining tracks SH excludes:' , 100*np.round(len(np.where(system_var_list[3] > crit_weak)[0])/len(system_var_list[3]), 3), '% (', len(np.where(system_var_list[3] > crit_weak)[0]), 'of', len(system_var_list[3]), ')\n' )
    str_for_latex += str(int(100*np.round(len(np.where(system_var_list[3] > crit_weak)[0])/len(system_var_list[3]), 2))) + ' \\ '  
    
print('str for latex:', str_for_latex)

plt.ylim(ymin=0)
plt.legend()
plt.xlabel(plotname)
if len(bins)>1: plt.xlim([bins[0], bins[-1]])
plt.ylabel('Normalized frequency')
# plt.gca().axes.get_yaxis().set_visible(False)
# ax.set_yticklabels([])
plt.gca().axes.get_yaxis().set_ticklabels([])

print(plotvarchar)

if save:
    savefile= savedir+ plotvarchar
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight', dpi= 100)        
