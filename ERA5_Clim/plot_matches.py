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



fignr= 3
savedir= homedir + '/Polar_Low/ERA5_PLclim/Figs/Match_compare/'
save= True

test= '2013_03'
# 

version= "fram_run3"
ending= 'atl'

durlim= 6

start_year= 2004
end_year= 2019


matchdist= 150


# PL_list_name, minObs= 'Rojo', 4
# PL_list_name, minObs = 'Gunnar', 4
# PL_list_name, minObs = 'Smirnova', 2
PL_list_name, minObs, ending = 'Yanase', 1, 'pac' #it has only one location per PL


time_start= False #in this case the first and last index of alltracks will be used
time_start= '2008-1-1'
time_end= '2009-1-1'




"""get all tracks"""
csv_name= 'merged_tracks_'+version+'_'+ending+f'_{start_year}-{end_year}'

if durlim: csv_name += 'durlim'+str(durlim)
tracks= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/mergedtracks/"+csv_name+'.csv')
tracks['time']= pd.to_datetime(tracks['time'])

tracks= tracks.set_index(['ID', 'step'])


"""alltracks - for a test"""
# tracks= pd.read_csv(Mediadir+"data/ERA5_Clim/track_lists/my_PC/tracks_"+test+".csv").set_index(['ID', 'step'])
# tracks['time']= pd.to_datetime(tracks['time'])






"""get PL list"""
if PL_list_name == 'Rojo':

    Rojo = import_Rojo_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv")
    PLlist= Rojo


if PL_list_name == 'Gunnar':
    STARS= import_Gunnar_STARS(Mediadir= Mediadir)

    STARS= S_Obs_nr(STARS)
    PLlist= STARS



if PL_list_name == 'Smirnova':

    Smir= pd.read_excel(Mediadir +'PL/PLclim/Smirnova_PL_1995-2009.xls', 'swaths')
    Smir['ID']= np.nan

    ID_counter= 1
    
    for i in range(len(Smir)): #make the ID
        if np.isnan(Smir.loc[i].Year):
            ID_counter += 1
        else:
            Smir.loc[i, 'ID']= ID_counter

    Smir = Smir.dropna(axis= 0, how= 'all') #remove empty rows

    #make the time
    Smir =Smir.assign(time= [f"{str(Year)}-{str(Month)}-{str(Day)} {str(Time)}" for Year, Month, Day, Time in zip(Smir['Year'].astype(int), Smir['Month'].astype(int), Smir['Day'].astype(int), Smir['Time'])])
    Smir['time']= pd.to_datetime(Smir['time'])
    Smir= Smir.drop(['Year', 'Month', 'Day', 'Time', 'Max Wind, m/s'], axis= 1)

    Smir= Smir.rename(columns= {'Latitude':'lat', 'Longitude': 'lon'})
    Smir['ID']= Smir['ID'].astype(int)
    Smir= S_Obs_nr(Smir)
    
    PLlist= Smir

if PL_list_name == 'Yanase':
    Yan= pd.read_csv(Mediadir +'PL/PLclim/Yanase_PLs.csv')
    
    Yan['time']= pd.to_datetime(Yan['time'], format= "%H%M UTC %d %b %Y") #does not work with pands version at university
    # Yan['time']= pd.to_datetime(Yan['time'], format= "%H%M UTC %d %b %Y")

    Yan['ID']= np.arange(1, len(Yan)+1)
    Yan['Obs']= 1
    
    PLlist= Yan



"""common preparation"""

PLlist['time']= PLlist['time'].dt.round('H')
PLlist= PLlist.set_index(['ID', 'Obs'])

PLlist= PLlist[(PLlist['time'] > time_start) & (PLlist['time'] < time_end)]

PLlist= PLlist.join( (PLlist.reset_index().groupby('ID').last()['Obs']).rename('nObs'), on='ID')

PLlist= PLlist[PLlist['nObs'] >= minObs]

print('remove this when September is done')
PLlist= PLlist[PLlist.time.dt.month != 9]


PL_indexes= remove_dublicate(PLlist.index.get_level_values(0))


"""match list"""
PL_matches= pd.read_csv(Mediadir+"/data/ERA5_Clim/track_lists/matchPLlist_"+ PL_list_name +f"_dist{matchdist}.csv")
PL_matches= PL_matches.set_index(['ID', 'step'])
PL_matches['time']= pd.to_datetime(PL_matches['time'])






"""map"""
fig = plt.figure(fignr, figsize= (11,8))
plt.clf()
ax= Plot_Polar_Stereo(fig, central_longitude= 20, extent= [-10, 40, 60, 77], subplot= (1,1,1), scalebar=False)
# ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-45, 45, 50, 80], subplot= (1,1,1))
# ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, 30, 80], subplot= (1,1,1))
#ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-50, 80, 40, 80], subplot= (1,1,1))



for index in PL_indexes:
    print(index)
    
    """map"""
    fig = plt.figure(fignr, figsize= (10,8))
    plt.clf()
    ax= Plot_Polar_Stereo(fig, central_longitude= 20, extent= [-10, 40, 60, 77], subplot= (1,1,1), scalebar=False)
    # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-45, 45, 50, 80], subplot= (1,1,1))
    # ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, 30, 80], subplot= (1,1,1))
    #ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-50, 80, 40, 80], subplot= (1,1,1))


    
    "PL"
    PL= PLlist.loc[index]
    ax.plot(PL.lon.values, PL.lat.values, c= 'g',  lw= 2, transform= ccrs.PlateCarree())
    ax.plot(PL.lon.values, PL.lat.values, 'x', c= 'g',  transform= ccrs.PlateCarree())
    ax.plot(PL.lon.values[0], PL.lat.values[0], marker='o', c= 'g',  transform= ccrs.PlateCarree() , label= PL_list_name+ ' '+str(index) )
    
    PL_middletime= PL['time'].iloc[PL['nObs'].iloc[0]//2]
    ax.plot(PL.lon.values[PL['nObs'].iloc[0]//2], PL.lat.values[PL['nObs'].iloc[0]//2], marker='s', c= 'g',  transform= ccrs.PlateCarree())



    """plot map"""
    # var='var138'#'vo'
    # ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/vorticity_era5_{PL_middletime.year}_{PL_middletime.month:02}.nc")
    var= 'vo'
    ds= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/vorticity_shift_era5_{PL_middletime.year}_{PL_middletime.month:02}.nc")

    ds= ds.sel(time= PL_middletime)
    ds= ds.isel(plev= 0)

    vextr= np.max([np.max(ds[var]), -np.min(ds[var])])*.7
#    # cf= ax.pcolormesh(ds.lon - 0.25, ds.lat + 0.125, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
    cf= ax.contourf(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), levels= 15, cmap= 'RdBu_r', vmin= -vextr, vmax= vextr)
#
    ax.contour(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), colors= 'red', levels=np.array([1.5, 2.0])*1E-4 )
    cb= fig.colorbar(cf, ax= ax, shrink=0.7)
    # cb= fig.colorbar(cf, ax= ax, shrink=0.5, orientation='horizontal') #, extend='both')
    varlabel= r'Relative vorticity [s$^{-1}$]'
    #if plevel != None: varlabel = var+#'_'+str(plevel)
    cb.set_label(varlabel , size=11)    
    # cb.set_label(ds[var].long_name + ' ['+ ds[var].units + ']', size=11)    

    # var= 'var151'  #'msl'
    # ds2= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/slp_era5_{PL_middletime.year}_{PL_middletime.month:02}.nc")
    var='msl'
    ds2= xr.open_dataset(Mediadir + f"data/ERA5_Clim/ERA5_data/slp_shift_era5_{PL_middletime.year}_{PL_middletime.month:02}.nc")
    ds2= ds2.sel(time= PL_middletime)
    cs= ax.contour(ds2.lon, ds2.lat, ds2[var]/100, np.arange(900, 1100, 5), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
    #if numbers_cont:
    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')
    
    
    """now tracks"""
    tracks_now= tracks[tracks['time']== PL_middletime]
    
    ax.scatter(tracks_now.lon.values, tracks_now.lat.values, marker= 's', c= 'purple',  transform= ccrs.PlateCarree(), zorder= 2)
    
    for ti, track_ID in enumerate(tracks_now.index.get_level_values('ID')):
        ax.plot(tracks.loc[track_ID].lon.values, tracks.loc[track_ID].lat.values, c='purple',  transform= ccrs.PlateCarree())
        ax.plot(tracks.loc[track_ID].lon.values, tracks.loc[track_ID].lat.values, 'x',  c='purple',  transform= ccrs.PlateCarree())
        if ti== 0: ax.plot(tracks.loc[track_ID].lon.values[0], tracks.loc[track_ID].lat.values[0], marker= 'o', c='purple',  transform= ccrs.PlateCarree(), label= 'all tracks')
        else: ax.plot(tracks.loc[track_ID].lon.values[0], tracks.loc[track_ID].lat.values[0], marker= 'o', c='purple',  transform= ccrs.PlateCarree())

    """matched PL"""
    PLmatched= PL_matches[PL_matches['PLlist_match'] == index]
    if len(PLmatched)== 0:
        PLmatched= PL_matches[PL_matches['PLlist_match_2'] == index]
        if len(PLmatched)== 0:
            if 'PLlist_match_3' in PL_matches:
                PLmatched= PL_matches[PL_matches['PLlist_match_3'] == index]
            if len(PLmatched)== 0:
                print('No match')

    if len(PLmatched)!= 0:
        ax.plot(PLmatched.lon.values, PLmatched.lat.values, c='b',  transform= ccrs.PlateCarree())
        ax.plot(PLmatched.lon.values, PLmatched.lat.values, 'x',  c='b',  transform= ccrs.PlateCarree())
        ax.plot(PLmatched.lon.values[0], PLmatched.lat.values[0], marker= 'o', c='b',  transform= ccrs.PlateCarree(), label= 'match')
        ax.plot(PLmatched[PLmatched['time'] == PL_middletime].lon.values,
                   PLmatched[PLmatched['time'] == PL_middletime].lat.values, marker= 's', c= 'b',  transform= ccrs.PlateCarree())
        
        
    
    plt.legend()
    
    if save:
        savefile= savedir+ PL_list_name +'_'+str(index)+ '.png'
        print(savefile)
        plt.savefig(savefile , bbox_inches='tight')
        
