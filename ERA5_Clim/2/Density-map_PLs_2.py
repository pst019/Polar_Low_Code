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
#import matplotlib.colors as colors
from scipy import ndimage


import numpy as np
from f_useful import *
# from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs

#
save= True
# save= False
savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/Investigate_situation/'

write= False

fignr= 6

#track_dir= Mediadir+"/data/ERA5_Clim/track_lists/fram/"

hem= 'NH'

if hem == 'NH':
    version= "" #"fram_run3_"
    ending= 'atl-pac'
    start_year= 1979
    end_year= 2020
    terrain_thresh= 50
    terrain_dist = 2
    more_info= True
    post_process= True
    split, exceed_large_to_small= True, 40

if hem == 'SH':
    version= "fram_run1"
    ending= 'SH'
    start_year= 1988
    end_year= 2020
    terrain_thresh= 50
    terrain_dist=4
    more_info= True
    split, exceed_large_to_small= False, 40



durlim= 6




#matchdist= 150


"""import track list"""
# csv_name= 'merged_tracks_'+version+'_'+ending+f'_{start_year}-{end_year}'
# if durlim: csv_name += 'durlim'+str(durlim)
# if more_info: csv_name += '_moreinfo'

# if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'

# tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv')



csv_name= 'merged_tracks_'+version + ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
if more_info: csv_name += '_moreinfo'
if post_process: csv_name += '_post-process'
if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
# if split: csv_name += f'_split{exceed_large_to_small}'
print(csv_name)

tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv').set_index(['ID', 'step'])

time_start= False #in this case the first and last index of alltracks will be used
# time_start= '2008-1-1' #the vorticity and slp is needed in this time
# time_end= '2009-1-1'

if time_start != False:
    tracks.time= pd.to_datetime(tracks.time)
    tracks= tracks[(tracks['time'] >= time_start) & (tracks['time'] < time_end)]

print('nr tracks:', len(tracks))


start_year= 1980
end_year= 2019
version= "fram_run3_"
terrain_dist= 2

csv_name= 'merged_tracks_'+version + ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
if split: csv_name += f'_split{exceed_large_to_small}'
print(csv_name)

if time_start != False:
    tracks_params= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+f'_params_{time_start}-{time_end}.csv').set_index(['ID', 'step'])
else:
    tracks_params= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'_params.csv').set_index(['ID', 'step'])

print('nr tracks with param:', len(tracks_params))


print('Params list is longer by:', len(set(tracks_params.index.values)- set(tracks.index.values)) )  #find out if one has more indexes


tracks= tracks.merge(tracks_params, left_index=True, right_index=True)

# tracks= tracks.reset_index()

"""filter out some tracks"""
tracks_ind= tracks.groupby(['ID']).max()
tracks_ind= tracks_ind.rename(columns={'step':'duration', 'vo': 'vo_max', 'area': 'area_max'})

tracks_ind['diameter_max']= np.sqrt(tracks_ind['area_max']) *4/np.pi

var_char_ind = [['theta_diff_500-sst_mean250', 'min'],
            ['theta_trop_mean250', 'min'] ]

for vc in var_char_ind:
    importvar, importchar= vc
    if importchar== 'max':
        tracks_ind_param=  tracks[importvar].groupby(['ID']).max().rename( importvar+'_'+importchar)
    if importchar== 'min':
        tracks_ind_param=  tracks[importvar].groupby(['ID']).min().rename( importvar+'_'+importchar)
    if importchar== 'mean':
        tracks_ind_param=  tracks[importvar].groupby(['ID']).mean().rename( importvar+'_'+importchar)
    
    tracks_ind= tracks_ind.merge(tracks_ind_param, left_index=True, right_index=True)

tracks_ind[tracks_ind['theta_diff_500-sst_mean250'] == -1000] = np.nan



tracks_ind_crit = tracks_ind
# tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['theta_diff_500-925_mean250_min'] < 20] #conservative Yanase
tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['theta_diff_500-sst_mean250_min'] < 12.5] #conservative Yanase

tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['theta_trop_mean250_min'] < 300] #conservative - Yan /Smir/Noer

tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['vo_max'] > 0.22] #in between
tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['diameter_max'] < 450] #approx middle

sel= tracks_ind_crit.index.values

tracks= tracks.reset_index().set_index('ID').loc[sel].reset_index()

print('Nr PLs:', len(sel))
print('Nr PLpoints:', len(tracks))

if hem == 'SH':
    """180E exclusion"""
    tracks['180E']= 0
    tracks.loc[np.abs(tracks.lon +.1) > 179.8, '180E']= 1
    
    tracks_ind= tracks[['ID', '180E']].groupby(['ID']).mean() 
    
    thresh180E= 30
    tracks= tracks.set_index('ID').loc[tracks_ind[tracks_ind['180E'] <= thresh180E/100].index].reset_index()
    tracks = tracks.drop(columns= '180E')
    
    outfile= Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+f'_180E-{thresh180E}.csv'
    print('write:', outfile)
    tracks.set_index('ID').to_csv(outfile)    


if hem == 'NH':
    """80N exclusion"""
    tracks['80N']= 0
    tracks.loc[np.abs(tracks.lat ) == 80, '80N']= 1
    
    tracks_ind= tracks[['ID', '80N']].groupby(['ID']).mean() 
    
    thresh80N= 30
    tracks= tracks.set_index('ID').loc[tracks_ind[tracks_ind['80N'] <= thresh80N/100].index].reset_index()
    tracks = tracks.drop(columns= '80N')
    
    # outfile= Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+f'_80N-{thresh80N}.csv'
    # print('write:', outfile)
    # tracks.set_index('ID').to_csv(outfile)    


"""make frequency distribution"""
loc_freq= tracks.groupby(['lat', 'lon']).size()
loc_freq= loc_freq.to_xarray().fillna(0)

loc_freq /= (end_year -start_year +1) #per year


"""simple density map per pixel"""
fig = plt.figure(fignr, figsize= (13,10) )
fignr+=1
plt.clf()

if hem == 'NH':
    ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, 30, 90], scalebar=False, subplot= (1,1,1), hem='NH', circular= True)
elif hem == 'SH':
    ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, -90, -30], scalebar=False, subplot= (1,1,1), hem='SH', circular= True)


cf= ax.pcolormesh(loc_freq.lon, loc_freq.lat, loc_freq.values, transform= ccrs.PlateCarree(), cmap='Reds' ) #, levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)
cb= fig.colorbar(cf, ax= ax, shrink=0.7)#, orientation= 'horizontal')
cb.set_label('Track points per grid cell per year')

plt.tight_layout()



"""simple density map per box, noramized"""
per_dist = 500 #km
normalize= True


if per_dist:
    loc_freq.values= ndimage.uniform_filter(loc_freq, size=(1+int(per_dist/(110/4)), 1+ int(1/np.cos(np.deg2rad(50)) *per_dist/(110/4))) )
else:
    loc_freq.values= ndimage.uniform_filter(loc_freq, size=(6, 10))

if normalize: loc_freq /=np.cos(np.deg2rad(loc_freq.lat))
if per_dist: loc_freq *= (per_dist/ (0.25*110))**2


fig = plt.figure(fignr, figsize= (13,10) )
fignr+=1
plt.clf()

if hem == 'NH':
    ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, 30, 90], scalebar=False, subplot= (1,1,1), hem='NH', circular= True)   
elif hem == 'SH':
    ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, -90, -30], scalebar=False, subplot= (1,1,1), hem='SH', circular= True)


# cf= ax.pcolormesh(ds.lon - 0.25, ds.lat + 0.125, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
# logarithmic
# cf= ax.pcolormesh(loc_freq.lon, loc_freq.lat, loc_freq.values, transform= ccrs.PlateCarree(), norm=colors.LogNorm(vmin=loc_freq.min(), vmax=loc_freq.max()) ) #, levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)

cf= ax.pcolormesh(loc_freq.lon, loc_freq.lat, loc_freq.values, transform= ccrs.PlateCarree(), cmap='Reds' ) #, levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)
cb= fig.colorbar(cf, ax= ax, shrink=0.7)#, orientation= 'horizontal')

if per_dist:
    cb.set_label(f'Annual track points within {per_dist}km distance')
else:
    cb.set_label('Track points per grid cell')

plt.tight_layout()


if save:
    savefile= savedir+ f'{hem}_track_density.png'
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight', dpi= 100)   
    
    
    
if write:
    print('write')
    outfile= Mediadir+"/data/ERA5_Clim/track_lists/derived_PLs/PLs-from_"+csv_name
    if time_start != False: outfile += f'_{time_start}-{time_end}'
    outfile += '.csv'
    print(outfile)
    tracks.set_index('ID').to_csv(outfile)    
