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

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import scipy.signal as signal
import matplotlib as mpl
from cycler import cycler

save=True
# save= False
savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/PLvsTRACKS_distr/'

fignr = 2

version= "" #"fram_run3"
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
split= False

terrain_dist, terrain_thresh= 2, 50
# terrain_excl= False



percentile= 5
PLs_larger= False

# PL_name= 'Rojo'
PL_name = 'Noer'
# PL_name = 'Smirnova'
# PL_name = 'Yanase' #it has only one location per PL
PL_name_list= [PL_name]

PL_name_list= ['Noer', 'Rojo', 'Smirnova', 'Yanase', 'Golubkin']
# PL_name_list= ['Smirnova', 'Yanase']
# PL_name_list= ['Noer', 'Rojo']
# PL_name_list= ['Noer', 'Yanase', 'Golubkin']


var_char_ind = [#['theta_diff_500-925_mean250', 'min'],
                ['theta_diff_500-sst_mean250', 'min'],
                ['SST-T_500_mean250', 'max'],
            ['theta_trop_mean250', 'mean'] ]


# Stype = 'system'
# plotchar= 'mean'
# plotchar= 'max'
# plotchar='min'
#plotchar='med'
#plotchar='initial'
#plotchar='final'



bins, plotname, plotname2= [], '', ''






# plotvar, plotchar='U_trop_polew_mean250', 'min' #'max'
plotvar, plotchar, plotname='theta_trop_mean250', 'mean',  r'$\theta_{trop}$ [K]' #'max'
# plotvar, plotchar='theta_trop_mean500', 'min' #'max'

# 
# plotvar, plotchar='theta_diff_500-sst_mean250', 'mean' #'max'
# plotvar, plotchar='theta_diff_500-skt_mean250', 'min' #'max'
# plotvar, plotchar='theta_diff_500-925_mean250', 'min'

# plotvar, plotchar, PLs_larger='U_surface_max250', 'max', True
# plotvar, plotchar, PLs_larger, plotname='vo_max', False, True, 'Rel. vorticity, lifetime max [10$^{-3}$s${-1}$]'
# plotvar, plotchar, PLs_larger= 'vo_tendency', False, True 

# plotvar, plotchar='diameter_max', False
# plotvar, plotchar, bins, plotname='diameter_max', False, np.arange(0, 901, 25), 'Diameter [km]'

# plotvar, plotchar, bins, PLs_larger= 'duration', False, np.arange(6, 200, 6), True

# plotvar, plotchar, PLs_larger, bins= 'near_land', False, False, np.arange(0, 1.01, 0.025)
# plotvar, plotchar, bins='month', False, np.arange(0.5, 13, 1)



# plotvar2, plotchar2='theta_trop_mean250', 'min' #'max'
# plotvar2, plotchar2, plotname2='SST-T_500_mean250', 'max', r'SST - T$_{500hPa}$ [K]'

plotvar2, plotchar2, plotname2='theta_diff_500-sst_mean250', 'min',  r'$\theta_{500hPa} - \theta_{SST}$ [K]' #'max'
# plotvar2, plotchar2, bins2, plotname2='vo_max', False, np.arange(0, 90.1, 0.5), 'Relative vorticity [10$^{-5}$s$^{-1}$]'


"""the plot variable"""
if plotchar == False: plotvarchar= plotvar
else:
    plotvarchar= plotvar+'_'+plotchar
    if [plotvar, plotchar] not in var_char_ind: var_char_ind += [[plotvar, plotchar]]

if plotchar2 == False: plotvarchar2= plotvar2
else:
    plotvarchar2= plotvar2+'_'+plotchar2
    if [plotvar2, plotchar2] not in var_char_ind: var_char_ind += [[plotvar2, plotchar2]]

if plotname == '': plotname= plotvarchar
if plotname2 == '': plotname2= plotvarchar2


"""import lists"""
csv_name= 'merged_tracks_'+version+ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
if more_info: csv_name += '_moreinfo'
if post_process: csv_name += '_post-process'
if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'

print(csv_name)

tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks2/"+csv_name+'.csv')

time_start= False #in this case the first and last index of alltracks will be used
# time_start= '2008-1-1' #the vorticity and slp is needed in this time
# time_end= '2009-1-1'

if time_start != False:
    tracks.time= pd.to_datetime(tracks.time)
    tracks= tracks[(tracks['time'] >= time_start) & (tracks['time'] < time_end)]



print('now it is only compared against NH tracks')

systems= [tracks] #list that gets tracks and PLs
systems_name_list= ['tracks'] + PL_name_list
systems_ind= []

for PL_name in PL_name_list:
    PLfile= Mediadir+"/data/ERA5_Clim/track_lists/matchPLlist2/matchPLlist_"+ PL_name +f"_dist{matchdist}"
    if time_diff and PL_name != 'Yanase': PLfile += f'_timediff{time_diff}'
    PLfile += f'_{csv_name}.csv'
    PLs=pd.read_csv(PLfile)
    systems += [PLs]


for system, system_name in zip(systems, systems_name_list):
    print(system_name)
    system_ind= system.groupby(['ID']).max()
    system_ind= system_ind.rename(columns={'step':'duration', 'vo': 'vo_max', 'area': 'area_max'})
    

    """get individual params for systems (tracks + PLs) - calculate the system_ind variable """
    if time_start != False and system_name == 'tracks':
        system_params= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks2/"+csv_name+f'_params_{time_start}-{time_end}.csv').set_index(['ID', 'step'])
    elif system_name== 'tracks':
        system_params= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks2/"+csv_name+'_params.csv').set_index(['ID', 'step'])
    else: #for the PLs
        PLfile_params = Mediadir+"/data/ERA5_Clim/track_lists/matchPLlist2/matchPLlist_"+system_name+f'_dist{matchdist}'
        if time_diff and system_name != 'Yanase': PLfile_params += f'_timediff{time_diff}'
        system_params= pd.read_csv(PLfile_params+ f'_{csv_name}_params.csv').set_index(['ID', 'step'])
    
    system_params[system_params['theta_diff_500-sst_mean250'] == -1000] = np.nan
    system_params[system_params['SST-T_500_mean250'] == -1000] = np.nan
    
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

    if plotvar == 'vo_tendency':
        system[plotvar]= np.nan
        for index in system_ind.index:
            track_now=system[system.ID == index]
            
            #just when tracks are time reduced:
            if len(track_now) >= 6:
                if len(track_now)>= 11:  window= 11
                else: window= len(track_now) -(len(track_now)+1)%2 #to get the largest uneven number smaller than the length of vorticity
            
                # track_now= track_now.assign(vo_tendency= signal.savgol_filter(track_now['vo'].values, window_length=window, polyorder=2, deriv= 1) )
                vo_tend= signal.savgol_filter(track_now['vo'].values, window_length=window, polyorder=2, deriv= 1)
    
                system.loc[system.ID == index, plotvar]= vo_tend
            
        system_ind[plotvar]=  system.set_index(['ID', 'step'])[plotvar].groupby(['ID']).max()


    systems_ind += [system_ind]







"""tracks after crit"""
tracks_ind_crit = systems_ind[0]
# tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['theta_diff_500-925_mean250_min'] < 20] #conservative Yanase
tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['theta_diff_500-sst_mean250_min'] < 11] #conservative Yanase
# tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['SST-T_500_mean250_max'] > 41.5] #conservative Yanase

tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['theta_trop_mean250_mean'] < 300.8] #conservative - Yan /Smir/Noer

tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['vo_max'] > 0.2] #in between
tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['diameter_max'] < 430] #approx middle


# systems_ind += [tracks_ind_crit]
systems_ind = systems_ind[:1]+ [tracks_ind_crit]+ systems_ind[1:]

# systems_name_list += ['tracks after crit'] 
systems_name_list =systems_name_list[:1] + ['tracks after crit'] + systems_name_list[1:]


print('Tracks remaining with criteria:', np.round(100*len(tracks_ind_crit)/len(systems_ind[0]), 1), '%' , len(tracks_ind_crit), 'of', len(systems_ind[0]) ,'\n')













# system_var_list= [system[plotvarchar] for system in systems_ind]
# system_var_list2= [system[plotvarchar2] for system in systems_ind]

system_var_list= [system[[plotvarchar, plotvarchar2]] for system in systems_ind]

   
if len(bins)== 0:
    bins= np.linspace(system_var_list[0].min(), system_var_list[0].max(), 20)





"""plot density distribution"""
mpl.rcParams['axes.prop_cycle'] = cycler(color='kbgrcmy')

plt.figure(fignr,  figsize= (8, 8))
fignr+=1
plt.clf() 


xplot= np.linspace(bins[0], bins[-1], 100)
 

for system_var, system_name in zip(system_var_list, systems_name_list):
# for system_var, system_var2, system_name in zip(system_var_list, system_var_list2, systems_name_list):
    # kernel= stats.gaussian_kde(system_var.dropna())
    plt.scatter(system_var[plotvarchar], system_var[plotvarchar2], s= 4, label= system_name)
    # plt.scatter(system_var.dropna(), Stoll[var2])
    
    # print(system_name, 'corr: ', np.round(np.corrcoef(system_var[plotvarchar], system_var[plotvarchar2])[0,1], 2))
    print(system_name, 'corr: ', np.round(system_var[[plotvarchar,plotvarchar2]].corr().iloc[0,1], 2))

    # if system_name== 'tracks after crit':    p= plt.plot(xplot, kernel(xplot), '--', label= system_name )
    # else:    p= plt.plot(xplot, kernel(xplot), label= system_name )

    # if system_name in systems_name_list[2:]:
    #     if PLs_larger:
    #         crit= np.percentile( system_var, percentile)
    #         plt.scatter(crit, 0, marker='^', c= p[0].get_color(), s= 200 )
        
    #         #if PL value is larger than track value
    #         # print(PL_name, np.round(crit, 3), 'excludes: ', np.round(len(np.where(tracks_var < crit)[0])/len(tracks_var), 2))
    #         print(system_name, 'crit:', np.round(crit, 3), 'excludes:' , 100*np.round(len(np.where(system_var_list[0] < crit)[0])/len(system_var_list[0]), 3), '% (', len(np.where(system_var_list[0] < crit)[0]), 'of', len(system_var_list[0]), ')' )
    #         print('Of remaining tracks excludes:' , 100*np.round(len(np.where(system_var_list[1] < crit)[0])/len(system_var_list[1]), 3), '% (', len(np.where(system_var_list[1] < crit)[0]), 'of', len(system_var_list[1]), ')\n' )
    
    
    #     else:
    #         crit= np.percentile( system_var, 100- percentile)
    #         plt.scatter(crit, 0, marker='^', c= p[0].get_color(), s= 200 )
        
    #         #if PL value is larger than track value
    #         print(system_name, 'crit:', np.round(crit, 3), 'excludes:' , 100*np.round(len(np.where(system_var_list[0] > crit)[0])/len(system_var_list[0]), 3), '% (', len(np.where(system_var_list[0] > crit)[0]), 'of', len(system_var_list[0]), ')' )
    #         print('Of remaining tracks excludes:' , 100*np.round(len(np.where(system_var_list[1] > crit)[0])/len(system_var_list[1]), 3), '% (', len(np.where(system_var_list[1] > crit)[0]), 'of', len(system_var_list[1]), ')\n' )
    

# plt.ylim(ymin=0)
plt.legend()
plt.xlabel(plotname)
plt.ylabel(plotname2)   

# plt.ylabel('Normalized density')

print(plotvarchar)

if save:
    savefile= savedir+ plotvarchar+'_'+ plotvarchar2
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight')        
