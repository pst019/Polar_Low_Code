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
savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/PL_density_map/'

write= True
write= False

import_PLs= False #if True the matched PL list is imported
Plot_take_time= False # True

if import_PLs: write= False

per_dist = 250 #km

# version='excl_prop-dist'

fignr= 4

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

# extra= '' 
extra= '_land-dist-140'


time_start= False #in this case the first and last index of alltracks will be used
# time_start= '2008-1-1' #the vorticity and slp is needed in this time
# time_end= '2009-1-1'

#matchdist= 150

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
    boxlist= [ ['Amundsen Sea', [-130, -60, -45, -70] ] #W, E, N, S
    , ['New Zealand', [160, 200, -45, -70] ]
    , ['Australia', [80, 140, -45, -70] ]
    ]

"""import track list"""

csv_name= 'merged_tracks_'+version + ending+f'_{start_year}-{end_year}'
if durlim: csv_name += 'durlim'+str(durlim)
if more_info: csv_name += '_moreinfo'
if post_process: csv_name += '_post-process'
if terrain_thresh: csv_name += f'_terrain{terrain_dist}-excl{terrain_thresh}'
# if split: csv_name += f'_split{exceed_large_to_small}'


if import_PLs:
    PLfile= Mediadir+"/data/ERA5_Clim/track_lists/derived_PLs/PLs-from_"+csv_name+ extra
    if time_start != False: PLfile += f'_{time_start}-{time_end}'
    PLfile += '.csv'
    print('Import:', PLfile)
    tracks= pd.read_csv(PLfile)    


else:
    """make the PL list"""
    print(csv_name)

    tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks2/"+csv_name+'.csv').set_index(['ID', 'step'])
    
    
    if time_start != False:
        tracks.time= pd.to_datetime(tracks.time)
        tracks= tracks[(tracks['time'] >= time_start) & (tracks['time'] < time_end)]
    
    print('nr track points:', len(tracks))
     

    if time_start != False:
        tracks_params= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks2/"+csv_name+f'_params_{time_start}-{time_end}.csv').set_index(['ID', 'step'])
    else:
        tracks_params= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks2/"+csv_name+'_params.csv').set_index(['ID', 'step'])
    
    print('nr track points with param:', len(tracks_params))
    
    
    print('Params list is longer by:', len(set(tracks_params.index.values)- set(tracks.index.values)) )  #find out if one has more indexes
    
    
    tracks= tracks.merge(tracks_params, left_index=True, right_index=True)
    

    
    # tracks= tracks.reset_index()
    
    """filter out some tracks"""
    tracks_ind= tracks.groupby(['ID']).max()
    tracks_ind= tracks_ind.rename(columns={'step':'duration', 'vo': 'vo_max', 'area': 'area_max'})
    
    tracks_ind['diameter_max']= np.sqrt(tracks_ind['area_max']) *4/np.pi

    tracks_ind['prop_dist'] = distance( (tracks.groupby(['ID']).first()['lat'], tracks.groupby(['ID']).first()['lon']),
              (tracks.groupby(['ID']).last()['lat'], tracks.groupby(['ID']).last()['lon']) )

    
    # var_char_ind = [['theta_diff_500-sst_mean250', 'min'],
                # ['theta_trop_mean250', 'min'],
                # ['land_dist', 'max']] #, 140, 'larger']]

    var_char_ind = [['SST-T_500_mean250', 'max'],
                ['theta_trop_mean250', 'min'],
                ['land_dist', 'max']] #, 140, 'larger']]
    
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
    tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['SST-T_500_mean250_max'] > 41.5] 
    tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['theta_trop_mean250_min'] < 298] 
    tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['vo_max'] > 0.20] 
    tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['diameter_max'] < 430] 

    if extra== '_land-dist-140': tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['land_dist_max'] > 140] 



    # tracks_ind_crit= tracks_ind_crit[tracks_ind_crit['prop_dist'] > 300] 
    
    sel= tracks_ind_crit.index.values
    
    tracks= tracks.reset_index().set_index('ID').loc[sel].reset_index()
    
    print('Nr PLs:', len(sel))
    print('Nr PLpoints:', len(tracks))


    # if hem == 'SH':
    #     """180E exclusion"""
    #     tracks['180E']= 0
    #     tracks.loc[np.abs(tracks.lon +.1) > 179.8, '180E']= 1
        
    #     tracks_ind= tracks[['ID', '180E']].groupby(['ID']).mean() 
        
    #     thresh180E= 30
    #     tracks= tracks.set_index('ID').loc[tracks_ind[tracks_ind['180E'] <= thresh180E/100].index].reset_index()
    #     tracks = tracks.drop(columns= '180E')
        
    #     outfile= Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+f'_180E-{thresh180E}.csv'
    #     print('write:', outfile)
    #     tracks.set_index('ID').to_csv(outfile)    
    
    
    # if hem == 'NH':
    #     """80N exclusion"""
    #     tracks['80N']= 0
    #     tracks.loc[np.abs(tracks.lat ) == 80, '80N']= 1
        
    #     tracks_ind= tracks[['ID', '80N']].groupby(['ID']).mean() 
        
    #     thresh80N= 30
    #     tracks= tracks.set_index('ID').loc[tracks_ind[tracks_ind['80N'] <= thresh80N/100].index].reset_index()
    #     tracks = tracks.drop(columns= '80N')
        
    #     # outfile= Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+f'_80N-{thresh80N}.csv'
    #     # print('write:', outfile)
    #     # tracks.set_index('ID').to_csv(outfile)    




"""make frequency distribution"""
loc_freq= tracks.groupby(['lat', 'lon']).size()
loc_freq= loc_freq.to_xarray().fillna(0)

loc_freq /= (end_year -start_year +1) #per year
loc_freq0= loc_freq.copy()


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

# plt.tight_layout()


if save:
    savefile= savedir+ f'{hem}_track_density{extra}_pertrackpoint.png'
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight', dpi= 100)   
    

"""simple density map per box, normalized"""
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


if save:
    savefile= savedir+ f'{hem}_track_density{extra}_uniform-filter-per-{per_dist}km.png'
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight', dpi= 100)   
    


"""complex density map per distance"""
if Plot_take_time:
    
    fig = plt.figure(fignr, figsize= (13,10) )
    fignr+=1
    plt.clf()
    
    if hem == 'NH':
        ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, 30, 90], scalebar=False, subplot= (1,1,1), hem='NH', circular= True)   
    elif hem == 'SH':
        ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, -90, -30], scalebar=False, subplot= (1,1,1), hem='SH', circular= True)
    
    loc_freq_dist= xr.DataArray(np.zeros(loc_freq.shape), coords=[loc_freq.lat, loc_freq.lon], dims=["lat", "lon"])
    
    
    print(f'Time passed: {round(time.perf_counter() - start,2)} seconds')
    
    if hem == 'NH':
        # the old way, does not work for the SH at 180E
        for lat in loc_freq0.lat.values:
            print(lat)
            for lon in loc_freq0.lon.values:
                dist= (110* np.sqrt( (loc_freq0.lat - lat)**2 + (np.cos(np.deg2rad(lat)) *(loc_freq0.lon- lon) )**2 ) )
            
                loc_freq_dist.loc[lat, lon]= loc_freq0.where(dist < per_dist).sum()


    elif hem == 'SH':
        nlon= len(loc_freq0.lon.values)
        for lon in loc_freq0.lon.values[int(nlon/4): int(3/4*nlon)]:
            print(lon)
            for lat in loc_freq0.lat.values:
                dist= (110* np.sqrt( (loc_freq0.lat - lat)**2 + (np.cos(np.deg2rad(lat)) *(loc_freq0.lon- lon) )**2 ) )
            
                loc_freq_dist.loc[lat, lon]= loc_freq0.where(dist < per_dist).sum()
        
        # loc_freq_dist.lon %= 360
        for lon in list(loc_freq0.lon.values[:int(nlon/4)]) + list(loc_freq0.lon.values[int(3/4*nlon):]):
            print(lon)
            for lat in loc_freq0.lat.values:
                lon_shift= lon%360
                dist= (110* np.sqrt( (loc_freq0.lat - lat)**2 + (np.cos(np.deg2rad(lat)) *((loc_freq0.lon)%360- lon_shift) )**2 ) )
            
                loc_freq_dist.loc[lat, lon]= loc_freq0.where(dist < per_dist).sum()
    

    
    loc_freq_dist /= 24
    
    plot_log=True
    
    if plot_log== True:
        # bounds= [0, 6, 12, 24, 48, 96, 144, 192, ]
        bounds= [0, .25, .5, 1, 2, 4, 7, 10]
        
        colors= ['white']+[(plt.cm.YlOrBr(h)) for h in range(256)]
        
        new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=len(bounds))
        norm= mpl.colors.BoundaryNorm(bounds, new_map.N)
        
        cf= ax.pcolormesh(loc_freq_dist.lon, loc_freq_dist.lat, loc_freq_dist.values, 
                          transform= ccrs.PlateCarree(), norm= norm, cmap=new_map ) #, levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)
    
    
    else: cf= ax.pcolormesh(loc_freq_dist.lon, loc_freq_dist.lat, loc_freq_dist.values, transform= ccrs.PlateCarree(), cmap='YlOrBr' ) #, levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)
    
    
    cb= fig.colorbar(cf, ax= ax, shrink=0.9)#, orientation= 'horizontal')
    
    cb.set_label(f'Annual time with polar-low activity within {per_dist}km distance [d]', size=14)
    cb.ax.tick_params(labelsize=14)
    
    
    if save:
        savefile= savedir+ f'{hem}_track_density{extra}_per{per_dist}km'
        if plot_log: savefile += '_log'
        savefile += '.png'
        print(savefile)
        plt.savefig(savefile , bbox_inches='tight', dpi= 150)   
    

    for box_name, loc in boxlist:
        ax.plot(list(np.linspace(0,1, 50)* [loc[1]-.3  - (loc[0]+.3)] + loc[0]+.3) + list(np.linspace(1,0, 50)* [loc[1]-.3 - (loc[0]+.3)] + loc[0]+.3) + [loc[0]+.3],
                50* [loc[3] ] + 50*[ loc[2] ] + [loc[3]], transform=ccrs.PlateCarree(), lw=2, label= box_name)
        
    plt.legend(loc= 'upper left')     
    
    if save:
        savefile= savedir+ f'{hem}_track_density{extra}_per{per_dist}km_box'
        if plot_log: savefile += '_log'
        savefile += '.png'
        print(savefile)
        plt.savefig(savefile , bbox_inches='tight', dpi= 150)   

    

    
if write:
    print('write')
    outfile= Mediadir+"/data/ERA5_Clim/track_lists/derived_PLs/PLs-from_"+csv_name+ extra
    if time_start != False: outfile += f'_{time_start}-{time_end}'
    outfile += '.csv'
    print(outfile)
    # tracks = tracks.drop(columns= ['Unnamed: 0', 'vortex_type','near_land'] )
    tracks.set_index(['ID', 'step']).to_csv(outfile)    

    
print(f'Time passed: {round(time.perf_counter() - start,2)} seconds')