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

fignr= 2

#track_dir= Mediadir+"/data/ERA5_Clim/track_lists/fram/"

hem= 'NH'

if hem == 'NH':
    version= "fram_run3"
    ending= 'atl-pac'
    start_year= 1980
    end_year= 2019
    terrain_thresh= 50
    terrain_dist = 4
    more_info= True

if hem == 'SH':
    version= "fram_run1"
    ending= 'SH'
    start_year= 1988
    end_year= 2020
    terrain_thresh= 50
    terrain_dist=4
    more_info= True



durlim= 6




#matchdist= 150


"""import track list"""
csv_name= 'merged_tracks_'+version+'_'+ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
if more_info: csv_name += '_moreinfo'

if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'

tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv')





"""simple density map"""
fig = plt.figure(fignr, figsize= (13,10) )
fignr+=1
plt.clf()

if hem == 'NH':
    # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-50, 80, 40, 80], subplot= (1,1,1))
    # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, 30, 80], subplot= (1,1,1))
    ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, 30, 90], scalebar=False, subplot= (1,1,1), hem='NH', circular= True)
    
if hem == 'SH':
    # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, -80, -30], subplot= (1,1,1))   
    ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, -90, -30], scalebar=False, subplot= (1,1,1), hem='SH', circular= True)


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
    
    outfile= Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+f'_80N-{thresh80N}.csv'
    print('write:', outfile)
    tracks.set_index('ID').to_csv(outfile)    


loc_freq= tracks.groupby(['lat', 'lon']).size()
loc_freq= loc_freq.to_xarray().fillna(0)


# loc_freq.values= gaussian_filter(loc_freq, sigma=(3, 5))
# loc_freq.values= ndimage.gaussian_filter(loc_freq, sigma=(6, 10))

per_dist = 200 #km





loc_freq /= (end_year -start_year +1) #per year

if per_dist:
    loc_freq.values= ndimage.uniform_filter(loc_freq, size=(1+int(per_dist/(110/4)), 1+ int(1/np.cos(np.deg2rad(50)) *per_dist/(110/4))) )
else:
    loc_freq.values= ndimage.uniform_filter(loc_freq, size=(6, 10))

normalize= True
if normalize: loc_freq /=np.cos(np.deg2rad(loc_freq.lat))

if per_dist: loc_freq *= (per_dist/ (0.25*110))**2

    # cf= ax.pcolormesh(ds.lon - 0.25, ds.lat + 0.125, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
# logarithmic
# cf= ax.pcolormesh(loc_freq.lon, loc_freq.lat, loc_freq.values, transform= ccrs.PlateCarree(), norm=colors.LogNorm(vmin=loc_freq.min(), vmax=loc_freq.max()) ) #, levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)

cf= ax.pcolormesh(loc_freq.lon, loc_freq.lat, loc_freq.values, transform= ccrs.PlateCarree(), cmap='Reds' ) #, levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)

#it appears like I have to exclude some stationary systems
#
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
