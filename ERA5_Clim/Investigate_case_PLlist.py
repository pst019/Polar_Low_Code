#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 16:29:52 2019

@author: pst019
only in python 3.6
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
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
from f_carto import *
from f_useful import *
from f_ERA5_clim import *
import scipy
from scipy.ndimage import uniform_filter
import random as rd

fignr= 1

save= True

# 
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


select_year= 2013
number_cases= 20 #number of randomly selected cases

savedir= homedir + f'/Polar_Low/ERA5_PLclim/Figs/2-PLcase-investigation_{select_year}_{hem}_2/'
if not os.path.exists(savedir): os.makedirs(savedir)

time_sel= 'orig_middle' #middle time step of the original PL for the PL list
# time_sel= 'orig_start'
# time_sel= 'orig_end'


if hem== 'NH':
    boxlist= [ ['Nordic Seas', [-15, 70, 80, 55] ] #W, E, N, S
        , ['Denmark Strait', [-45, -15, 75, 60] ]
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
PLfile += '.csv'
print('Import:', PLfile)
tracks= pd.read_csv(PLfile) 

tracks.time= pd.to_datetime(tracks.time)

# tracks= tracks.dropna(axis= 'columns') #theta diff is already set to nan
tracks[tracks['SST-T_500_mean250'] == -1000] = np.nan


tracks= tracks[tracks.time.dt.year == select_year]

"""just a box"""
# box_name, box_loc = boxlist[3]
# print(box_name)
# tracks= tracks.where((tracks.lon <= box_loc[1]) & (tracks.lon > box_loc[0]) &
#                           (tracks.lat <= box_loc[2]) & (tracks.lat > box_loc[3]) 
#                           ).dropna()


random_sample= rd.sample( remove_dublicate(tracks.ID), number_cases )

"""landdist < 140"""
# tracks_landdistmax= tracks.groupby(['ID']).max()['land_dist']
# tracks_landdistmax[tracks_landdistmax < 140]

# random_sample= rd.sample( list(tracks_landdistmax[tracks_landdistmax < 140].index.values), number_cases )


for index in random_sample:
    print(index)
    track= tracks[tracks.ID== index]
    
    
    """the track evolution"""
    fig = plt.figure(fignr, figsize= (10,10))
    # fignr += 1
    plt.clf()
    
    axs=fig.subplots(8, 1, sharex= 'col')
    
    axs[0].plot(track.time, track.vo*100)
    axs[0].set_ylabel('Vorticity')
    axs[0].plot([track.time.iloc[0], track.time.iloc[-1]], [20, 20], c= 'k', linewidth= 0.5, linestyle= '--')
    
    axs[1].plot(track.time, np.sqrt(track['area']) *4/np.pi)
    axs[1].set_ylabel('Diameter')
    axs[1].plot([track.time.iloc[0], track.time.iloc[-1]], [430, 430], c= 'k', linewidth= 0.5, linestyle= '--')
    
    axs[2].plot(track.time, track['SST-T_500_mean250'])
    axs[2].set_ylabel(r'SST-T$_{500}$')
    axs[2].plot([track.time.iloc[0], track.time.iloc[-1]], [41, 41], c= 'k', linewidth= 0.5, linestyle= '--')

    axs[3].plot(track.time, track['theta_diff_500-sst_mean250'])
    axs[3].set_ylabel(r'$\theta_{500} - \theta{SST}$')
    axs[3].plot([track.time.iloc[0], track.time.iloc[-1]], [11, 11], c= 'k', linewidth= 0.5, linestyle= '--')

    
    axs[4].plot(track.time, track['theta_trop_mean250'])
    axs[4].set_ylabel(r'$\theta_{trop}$')
    axs[4].plot([track.time.iloc[0], track.time.iloc[-1]], [301, 301], c= 'k', linewidth= 0.5, linestyle= '--')
    
    # axs[4].plot(track.time, track['slp'])
    # axs[4].set_ylabel('SLP')
        
    axs[5].plot(track.time, track['shear_strength_925-500_mean500'])
    axs[5].plot(track.time, track['shear_strength_925-500_mean250'], 'r')
    axs[5].set_ylabel('Shear strength')
        
    axs[6].plot(track.time, track[ 'shear_angle_925-500_mean500'])
    axs[6].plot(track.time, track[ 'shear_angle_925-500_mean250'], 'r')
    
    axs[6].set_ylabel('Shear angle')
    axs[6].set_ylim([0,360])
 
    axs[7].plot(track.time, track['land_dist'])
    axs[7].set_ylabel('Land dist')
    axs[7].set_ylim([0,360])   
 
    axs[7].set_xlabel('Time')
    
    if save:
        savefile= savedir+ str(index)+'evolution'
        savefile+= '.png'      
        print(savefile)
        plt.savefig(savefile , bbox_inches='tight', dpi= 100)   
        
    
    
    """map of a track"""
    for time_sel in [time_sel]: # ['orig_start']: #[ 'orig_middle' ]: #, 'orig_start', 'orig_end', 'matched_start', 'matched_end', 'matched_middle', 'matched_2_3']:
        # for index in PL_indexes:
            # print(index)
        if time_sel== 'orig_start':
            track_step= track.iloc[0]
        elif time_sel== 'orig_end':
            track_step= track.iloc[-1]
        elif time_sel== 'orig_middle':
            track_step= track.iloc[len(track)//2]        
        
        track_time= track_step.time
        
        """map"""
        fig = plt.figure(fignr +1, figsize= (10,8))
        # fignr += 1
        plt.clf()
        if hem == 'NH':
            ax= Plot_PlateCarree(fig, central_longitude= track_step.lon, extent= 
                                 [track_step.lon-25, track_step.lon+25, track_step.lat-15, track_step.lat+15], subplot= (1,1,1))
            # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, 30, 80], subplot= (1,1,1))
            # ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, 30, 90], scalebar=False, subplot= (1,1,1), hem='NH', circular= True)
        
        if hem == 'SH':
            ax= Plot_PlateCarree(fig, central_longitude= track_step.lon, extent= 
                                 [track_step.lon-35, track_step.lon+35, track_step.lat-25, track_step.lat+25], subplot= (1,1,1))

            # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, -80, -30], subplot= (1,1,1))   
            # ax= Plot_Polar_Stereo(fig, central_longitude= 0, extent= [-180, 180, -90, -30], scalebar=False, subplot= (1,1,1), hem='SH', circular= True)
    
      
    
        
        """plot map"""
        # if Plot_fields:
        # var='var138'#'vo'
        # ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/vorticity_era5_{track_time.year}_{track_time.month:02}.nc")
        if hem== 'NH' and select_year==2013:
            ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/vorticity_era5_{track_time.year}_{track_time.month:02}.nc")
        elif hem== 'NH':
            ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/vorticity_shift_era5_{track_time.year}_{track_time.month:02}.nc")
        elif hem== 'SH':
            ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/vorticity_SHm1_era5_{track_time.year}_{track_time.month:02}.nc")
        var= 'var138'
    
        ds[var]*= 1E5
        
        ds= ds.sel(time= track_time)
        ds= ds.isel(plev= 0)
        
        # vextr= np.max([np.max(ds[var]), -np.min(ds[var])])*.7
        #    # cf= ax.pcolormesh(ds.lon - 0.25, ds.lat + 0.125, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
        # cf= ax.contourf(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)
        #
        ax.contour(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), colors= 'red', levels=np.array([15, 20]) )
        # cb= fig.colorbar(cf, ax= ax, shrink=0.7)
        # # cb= fig.colorbar(cf, ax= ax, shrink=0.5, orientation='horizontal') #, extend='both')
        # varlabel= r'Relative vorticity [10$^{-5}$ s$^{-1}$]'
        #if plevel != None: varlabel = var+#'_'+str(plevel)
        # cb.set_label(varlabel , size=11)    
        # cb.set_label(ds[var].long_name + ' ['+ ds[var].units + ']', size=11)    
        
        # var= 'var151'  #'msl'
        # ds2= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/slp_era5_{track_time.year}_{track_time.month:02}.nc")
        # var='msl'
        # ds2= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/slp_shift_era5_{track_time.year}_{track_time.month:02}.nc")
        if hem== 'NH' and select_year==2013:
            ds2= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/slp_era5_{track_time.year}_{track_time.month:02}.nc")
        elif hem== 'NH':
            ds2= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/slp_shift_era5_{track_time.year}_{track_time.month:02}.nc")
        if hem== 'SH':
            ds2= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/slp_SH_era5_{track_time.year}_{track_time.month:02}.nc")
        var= 'var151'  #'msl'
        ds2= ds2.sel(time= track_time)
        
        cs= ax.contour(ds2.lon, ds2.lat, ds2[var]/100, np.arange(900, 1100, 5), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
        #if numbers_cont:
        plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')
        
    
        if hem== 'NH':
            ds3= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/levels_500_era5_{track_time.year}_{track_time.month:02}.nc")
            if select_year==2013: var='t'
            else: var='var130'
        if hem== 'SH':
            ds3= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/levels_500_SH_era5_{track_time.year}_{track_time.month:02}.nc")
            var= 'var130'  #'msl'
        ds3= ds3.sel(time= track_time)
        ds3= ds3.sel(plev= 92500)
        
        cf= ax.contourf(ds3.lon, ds3.lat, ds3[var]-273.15, transform= ccrs.PlateCarree(), levels= 41, cmap= 'RdBu_r', vmin= -20, vmax= 10)
    
        cb= fig.colorbar(cf, ax= ax, shrink=0.7)
        # cb= fig.colorbar(cf, ax= ax, shrink=0.5, orientation='horizontal') #, extend='both')
        varlabel= r'Temperature at 925hPa [K]' 
        cb.set_label(varlabel , size=11)    
       
        # cs= ax.contour(ds2.lon, ds2.lat, ds2[var]/100, np.arange(900, 1100, 5), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
        #if numbers_cont:
        # plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')
    
    
        """plot track"""
        ax.plot(track.lon.values, track.lat.values, c= 'g',  transform= ccrs.PlateCarree())
        ax.plot(track[track.time== track_time].lon.values, track[track.time== track_time].lat.values, marker='s', c= 'g',  transform= ccrs.PlateCarree())
        ax.plot(track.lon.values[0], track.lat.values[0], c= 'g', marker= 'o',  transform= ccrs.PlateCarree())
        
    
        
        if save:
            savefile= savedir+ str(index)+time_sel+'_map'
            savefile+= '.png'      
            print(savefile)
            plt.savefig(savefile , bbox_inches='tight', dpi= 100)  
    

