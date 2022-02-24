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
from scipy import stats


start= time.perf_counter()
#
save= True
# save= False
savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/PL_time_evol/'



fignr= 1

#track_dir= Mediadir+"/data/ERA5_Clim/track_lists/fram/"

hem= 'SH'

if hem == 'NH':
    ending= 'atl-pac'

elif hem == 'SH':
    ending= 'SH180'

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


var_char_ind = [
    #['theta_diff_500-sst_mean250', 'min', 11.0, 'smaller'],
    ['SST-T_500_mean250', 'max', 41.5, 'larger'],
            ['theta_trop_mean250', 'min', 298, 'smaller'],
            # ['land_dist', 'max', 140, 'larger'],
            ['vo', 'max', 20, 'larger'],
            ['diameter', 'max', 430, 'smaller']
            ]
    


"""import track list"""

csv_name= 'merged_tracks_'+version + ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
if more_info: csv_name += '_moreinfo'
if post_process: csv_name += '_post-process'
if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
# if split: csv_name += f'_split{exceed_large_to_small}'


""" import_PLs"""
PLfile= Mediadir+"/data/ERA5_Clim/track_lists/derived_PLs/PLs-from_"+csv_name
for vc in var_char_ind:
    importvar, importchar, excl_thresh, excl_type = vc
    PLfile += '_'+ importvar+ '-' + str(int(excl_thresh))
if time_start != False: PLfile += f'_{time_start}-{time_end}'
PLfile += '.csv'
print('Import:', PLfile)
tracks= pd.read_csv(PLfile)    

tracks.time= pd.to_datetime(tracks.time)


tracks= tracks.dropna(axis= 'columns')



if hem== 'NH':
    boxlist= [ ['Nordic Seas', [-15, 70, 80, 55] ] #W, E, N, S
        , ['Denmark Strait', [-45, -15, 75, 50] ]
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
#     boxlist= [ ['Amundsen Sea', [-130, -60, -45, -70] ] #W, E, N, S
#     , ['New Zealand', [160, 200, -45, -70] ]
#     , ['Australia', [80, 140, -45, -70] ]
#     ]





"""annual"""
fig= plt.figure(fignr, figsize= (9,6) )
fignr+=1
plt.clf()
ax= fig.add_subplot()




nr_box_prev= 0
for box_name, box_loc in boxlist:
    if box_name== 'New Zealand': tracks.lon%=360
        
    nr_box = tracks.where((tracks.lon <= box_loc[1]) & (tracks.lon > box_loc[0]) &
                          (tracks.lat <= box_loc[2]) & (tracks.lat > box_loc[3]) 
                          ).dropna().groupby(tracks.time.dt.year).size().reindex(np.arange(start_year, end_year+1), fill_value=0)/24
    plt.bar(nr_box.index.values, nr_box.values, bottom= nr_box_prev, label= box_name) #, width=1)

    slope, intercept, r_value, p_value, std_err = stats.linregress(nr_box.index.values, nr_box.values)
    # ax.plot(tracks_annual.index, intercept+ slope*tracks_annual.index, color='k', label= 'Trend')
    print(box_name, ' trend: ', np.round(slope, 2), 'pvalue of trend:', np.round(p_value,4))

    nr_box_prev += nr_box

    if box_name== 'New Zealand': tracks.lon[tracks.lon > 180] -= 360 #maybe this needs to be done differenly


tracks_annual= tracks.groupby(tracks.time.dt.year).size().reindex(np.arange(start_year, end_year+1), fill_value=0)/24
plt.bar(tracks_annual.index.values, tracks_annual.values- nr_box_prev, bottom= nr_box_prev, label= 'Rest') #, width=1)

handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles[::-1], labels[::-1]) #, title='Line', loc='upper left')


plt.xlim([1979-.5, 2020+.5])
plt.xlabel('Year')
plt.ylabel('Time of polar-low activity [d]')



"""statistic"""
print('mean', tracks_annual.mean(), ' std', tracks_annual.std())

slope, intercept, r_value, p_value, std_err = stats.linregress(tracks_annual.index, tracks_annual.values)
plt.plot(tracks_annual.index, intercept+ slope*tracks_annual.index, color='k', label= 'Trend')
print('Total trend: ', np.round(slope, 2), 'pvalue of trend:', np.round(p_value,4))


plt.legend(handles[::-1], labels[::-1]) #fontsize= 14)#, loc='lower left')

if save:
    savefile= savedir+ f'{hem}_time_evol_year'
    for vc in var_char_ind:
        importvar, importchar, excl_thresh, excl_type = vc
        savefile += '_'+ importvar + '-'+ str(int(excl_thresh))    
    savefile += '.png'
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight', dpi= 150)   


    
"""seasonal"""
fig= plt.figure(fignr, figsize= (9,6) )
fignr+=1
plt.clf()
ax= fig.add_subplot()

season_norm= np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]) /30

nr_box_prev= 0
for box_name, box_loc in boxlist:
    if box_name== 'New Zealand': tracks.lon%=360

    nr_box = tracks.where((tracks.lon <= box_loc[1]) & (tracks.lon > box_loc[0]) &
                          (tracks.lat <= box_loc[2]) & (tracks.lat > box_loc[3]) 
                          ).dropna().groupby(tracks.time.dt.month).size().reindex(np.arange(1,13)
                        , fill_value=0)/(end_year - start_year +1)/24/season_norm
    plt.bar(nr_box.index.values, nr_box.values, bottom= nr_box_prev, label= box_name) #, width=1)

    #have to add 0 in table
    nr_box_prev += nr_box
    
    if box_name== 'New Zealand': tracks.lon[tracks.lon > 180] -= 360 #maybe this needs to be done differenly



tracks_seasonal= tracks.groupby(tracks.time.dt.month).size().reindex(np.arange(1,13), fill_value=0)/(end_year - start_year +1)/24/season_norm
plt.bar(tracks_seasonal.index.values, tracks_seasonal.values- nr_box_prev, bottom= nr_box_prev, label= 'Rest') #, width=1)


handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1])

print('Should normalize by length of month')
    
# nordic_seasonal= tracks.where((tracks.lon <= boxNordic[1]) & (tracks.lat >= boxNordic[0]) & (tracks.lat <= boxNordic[2]) & (tracks.lat > boxNordic[3]) ).dropna().groupby(tracks.time.dt.month).size()/(end_year - start_year +1)/24
# plt.bar(tracks_seasonal.index.values, nordic_seasonal.values, label= 'Nordic Seas') #, width=1)

# plt.bar(bins, seasonaldistrNordic, label= 'Nordic Seas') #, width=1)


# adistrGreen= np.sum(countdens[:,:, int(np.where(d.lat== boxGreen[2])[0]): int(np.where(d.lat== boxGreen[3])[0]), int(np.where(d.lon== boxGreen[0]-360)[0]): int(np.where(d.lon== boxGreen[1]-360)[0])], axis= (1,2,3))
# plt.bar(bins, adistrGreen/4, bottom= adistrNordic/4, label='North Central Atlantic') #, width=1)
   
# plt.legend()
plt.xticks(np.arange(1,13))
plt.xlim([0.5, 12.5])
plt.xlabel('Month')
plt.ylabel('Time of polar-low activity [d]')

if save:
    savefile= savedir+ f'{hem}_time_seasonal'
    for vc in var_char_ind:
        importvar, importchar, excl_thresh, excl_type = vc
        savefile += '_'+ importvar + '-'+ str(int(excl_thresh))
            
    savefile += '.png'
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight', dpi= 150)   


    
print(f'Time passed: {round(time.perf_counter() - start,2)} seconds')