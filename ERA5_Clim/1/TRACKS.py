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
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/test_Denis/mc_era5-master/code') #to get obs_tracks_api


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cartopy.crs as ccrs

# from octant.core import TrackRun, OctantTrack
# from pathlib import Path

from f_carto import *
from f_useful import *
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/STARS-ERA5')
from f_STARS import *

from datetime import datetime, timedelta

# from obs_tracks_api import prepare_tracks

from octant.core import TrackRun, OctantTrack #conda install -c dennissergeev octant

from pathlib import Path
#from obs_tracks_api import prepare_tracks


fignr= 12

# save=False
# savedir= homedir+ 'Polar_Low/STARS-ana/Figs/Track_eval/'

#track_list=['test_2_1']
#now= np.datetime64("2010-12-05 00")
#ds_filetime="2010_12"

##track_list= [ 'test_050_7', 'test_050_9', 'test_025_9'] #'test_050', 'test_050_2', 'test_025_2']#, 'test_025', 'test_0525']
##now= np.datetime64("2010-12-09 13")
#
#track_list= ['test_050_11', 'test_025_11', 'test_041210_2']#, 'test_025', 'test_0525'] , 'test_050_10' 'test_050_3', 'test_050_4', 
#now= np.datetime64("2010-12-05 23")
##
###track_list= [ 'test_050_6', 'test_025_6']#, 'test_025', 'test_0525']'test_050_5',
###now= np.datetime64("2010-12-21 12")

# track_list= ['test_041210_5', 'test_041210_2']#, 'test_025', 'test_0525'] , 'test_050_10' 'test_050_3', 'test_050_4', 
# now= np.datetime64("2010-12-09 18")

# ds_filetime='2010_12'


#track_list= [ 'test_100311']#, 'test_025', 'test_0525']'test_050_5',
#now= np.datetime64("2011-03-11 06")

#track_list= [ 'test_120311']#, 'test_025', 'test_0525']'test_050_5',
#now= np.datetime64("2011-03-13 12")

# track_list= ['test_200311',  'test_200311_1', 'test_200311_2' ]#, 'test_025', 'test_0525']'test_050_5',
# now= np.datetime64("2011-03-20 07")
# 
# ds_filetime='2011_03'

#track_list= ['test_140313',  'test_140313_1', 'test_140313_2', 'test_180313_3']#, 'test_200311_2' ]#, 'test_025', 'test_0525']'test_050_5',
#now= np.datetime64("2013-03-19 12")

#track_list= [ 'test_140313_2', 'test_0313_4', 'test_0313_5']#, 'test_200311_2' ]#, 'test_025', 'test_0525']'test_050_5',
#now= np.datetime64("2013-03-15 12")
#
#track_list= [ 'test_110313_2', 'test_0313_4', 'test_0313_5']#, 'test_200311_2' ]#, 'test_025', 'test_0525']'test_050_5',
#now= np.datetime64("2013-03-12 06")
#
##track_list= ['test_260313_3', 'test_0313_4', 'test_0313_5']#, 'test_0313_6']
##now= np.datetime64("2013-03-27 15")
#
## track_list= ['test_030313_2']
## now= np.datetime64("2013-03-08 10")
#
#track_list= ['test_030313_2', 'test_0313_4', 'test_0313_5', 'test_0313_6']
#now= np.datetime64("2013-03-05 00")
#
#ds_filetime='2013_03'


#track_list= ['test_global', 'test_atlantic', 'test_pacific_2', '2013_03']#, 'test_025', 'test_0525'] , 'test_050_10' 'test_050_3', 'test_050_4', 
#now= np.datetime64("2013-03-01 02")

#ds_filetime='2013_03'


track_list= ['test6']#, 'test_025', 'test_0525'] , 'test_050_10' 'test_050_3', 'test_050_4', 
now= np.datetime64("2013-01-01 02")

ds_filetime='2013_01'


track_list= ['run3_2011_01_atl']
now= np.datetime64("2011-01-02 05")


#track_list= ['test2_total', 'test2_atlantic', 'test2_pacific']#, 'test_025', 'test_0525'] , 'test_050_10' 'test_050_3', 'test_050_4', 
#now= np.datetime64("2014-01-01 01")
#
#ds_filetime='2014_01'

vo_filter= False #False
Plot_movement= False

# """map"""
fig = plt.figure(fignr)
plt.clf()
# ax= Plot_Polar_Stereo(fig, central_longitude= 20, extent= [-10, 40, 60, 77], subplot= (1,1,1), scalebar=False)
#ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-45, 45, 50, 80], subplot= (1,1,1))
ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-180, 180, 30, 80], subplot= (1,1,1))
#ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-50, 80, 40, 80], subplot= (1,1,1))
    




##
##
#var='var138'#'vo'
#ds= xr.open_dataset(Mediadir + "data/ERA5_Clim/ERA5_data/vorticity_era5_"+ ds_filetime + '.nc')
#ds= ds.sel(time= now)
#ds= ds.isel(plev= 0)
#
#cmap= 'RdBu_r'
#
#if vo_filter:
#    import scipy.ndimage.filters as filters
#    variable= filters.uniform_filter(ds[var], size= (5, 5))
#
#    vextr= np.max([np.max(variable), -np.min(variable)])*.8
#    cf= ax.pcolormesh(ds.lon - 0.25, ds.lat + 0.125, variable, transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
#    ax.contour(ds.lon, ds.lat, variable, transform= ccrs.PlateCarree(), colors= 'red', levels=np.array([1.2, 1.5])*1E-4 )
#
#
#else:
#    vextr= np.max([np.max(ds[var]), -np.min(ds[var])])*.7
#    # cf= ax.pcolormesh(ds.lon - 0.25, ds.lat + 0.125, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
#    cf= ax.contourf(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), levels= 15, cmap= cmap, vmin= -vextr, vmax= vextr)
#
#    ax.contour(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), colors= 'red', levels=np.array([1.5, 2.0])*1E-4 )
#
##clevels= np.array([0, 0.3, .6, .9, 1.2, 1.5, 2])*1E-4
##cmap='Reds'
##norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=False)
##cf= ax.pcolormesh(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap)#, vmin= clevels[0], vmax=clevels[-1])#, extend='both')
#
#
#
#
#
#
#
#cb= fig.colorbar(cf, ax= ax, shrink=0.5, orientation='horizontal') #, extend='both')
#varlabel= var
##if plevel != None: varlabel = var+#'_'+str(plevel)
##    cb.set_label(varlabel + ' ['+ ds[var].units + ']', size=11)    
##cb.set_label(ds[var].long_name + ' ['+ ds[var].units + ']', size=11)    
#
##
#
##var='msl'
#ds2= xr.open_dataset(Mediadir + "data/ERA5_Clim/ERA5_data/slp_era5_"+ ds_filetime + '.nc')
#ds2= ds2.sel(time= now)
#var= 'var151'  #'msl'
#cs= ax.contour(ds2.lon, ds2.lat, ds2[var]/100, np.arange(900, 1100, 5), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
##if numbers_cont: plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')
#
#if Plot_movement:
#    ds3= xr.open_dataset(Mediadir + "data/ERA5_Clim/ERA5_data/wind_era5_"+ ds_filetime + '.nc')
#    ds3= ds3.sel(time= now)
#    ds3= ds3.isel(plev= 0)


track_color_list= ['b', 'g', 'purple', 'cyan']

for it, test in enumerate(track_list):
    track_dir= Mediadir+"/data/ERA5_Clim/tracks/"+test
    
    t_dir= Path(".") / track_dir
    tr = TrackRun(t_dir, columns=['lon', 'lat', 'vo', 'time', 'area', 'vortex_type', 'slp'])
    
    """conditions"""
    conditions = [ ("n", [lambda ot: ot.lifetime_h >= 1]),]
    tr.classify(conditions)
    
    """list of octant tracks"""
    tr_tracks= []
    
    for (_, track) in tr.gb:
        tr_tracks.append(track)
    
    ilegend= 0
    for track in tr_tracks:
        #all tracks that intersect in time with the STARS track
        #they do not end before STARS PL start and do not start after the STARS PL has ended
    #    if (not (track.time.values[-1] < now) and (not (track.time.values[0] > now)) ):
        if now in track.time.values:
            ilegend+=1
            track.lat+=0.1*it
            if ilegend== 1:  track.plot_track(ax=ax, label= test, color= track_color_list[it])
            else: track.plot_track(ax=ax, label='', color= track_color_list[it])
            
            tr_i_now = np.argmin(np.abs(track.time - now).values)
    #            print(tr_i_now)
            tracklon, tracklat= track.lonlat[tr_i_now]
    #        if ilegend== 1: ax.scatter(tracklon, tracklat, transform=ccrs.PlateCarree(), s= 50, marker="s", color= 'b', zorder= 2, label= 'Location now')
            ax.scatter(tracklon, tracklat, transform=ccrs.PlateCarree(), s= 50, marker="s", color= track_color_list[it], zorder= 2)
#           the vortex type
#            ax.text(tracklon, tracklat+ 0.1, track.vortex_type[tr_i_now], transform=ccrs.PlateCarree(), color=track_color_list[it], fontsize=13, fontweight= 'bold')
#            ax.text(tracklon, tracklat- 0.1, track.index[0][0], transform=ccrs.PlateCarree(), color=track_color_list[it], fontsize=13, fontweight= 'bold')
    
    
            ax.scatter(track.lonlat[:,0], track.lonlat[:,1], transform=ccrs.PlateCarree(), s= 50, marker="x", color= 'k', zorder= 2)
    
            if Plot_movement:
                ds3track= ds3.sel(lat= tracklat, lon= tracklon)
                tracklat_d = ds3track['v'].values*3.6/110
                tracklon_d = ds3track['u'].values*3.6/(np.cos(np.deg2rad(tracklat)) *110)
                
                ax.plot([tracklon, tracklon+tracklon_d], [tracklat, tracklat+tracklat_d], color='yellow')
    

#PlotWind(ds3.lon, ds3.lat, ds3['u'].values, ds3['v'].values, ax, nx= 15, ny=30, arrowtype= 'quiver', alen=15, scale=False)




Rojo = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv")
tdelta= np.abs((Rojo['time'] - now) / np.timedelta64(1, 'h'))
Rojo_now= Rojo[tdelta<= 4]
Rojo_nr_list= remove_dublicate(Rojo_now.ID)



#    Rojo_nr_list= Rojo[Rojo['time']== ts_datetime_now].ID.values
#    Rojo_nr_list= remove_dublicate(Rojo_nr_list)
#    print('Rojo nr list: ', Rojo_nr_list)
#    
#    if len(Rojo_nr_list)== 0: Rojo_nr_list= [ID.split('_')[0]]

for Rojo_nr in Rojo_nr_list:
    Rojo_now= Rojo[Rojo['ID'] == Rojo_nr] 
   
    R_tracks= prepare_tracks(Rojo_now, system_id_label="ID")
    
    for track in R_tracks:
        track.plot_track(ax=ax, color='r', label='Rojo track')   
        
        i_R_now= (np.abs(track.time - now)).idxmin()

        ax.scatter(track.lon.values[i_R_now], track.lat.values[i_R_now], transform=ccrs.PlateCarree(), s= 70, marker= "s", color= 'r', zorder= 2)
        ax.text(track.lon.values[i_R_now], track.lat.values[i_R_now]+ 0.1, track.ID[0], transform=ccrs.PlateCarree(), color='r', fontsize=13, fontweight= 'bold')
    
        time_plot_iRojo= track.time.values[i_R_now]
        print(track.ID.values[0], 'Rojo time: ', time_plot_iRojo)



plt.legend()
plt.title( str(now))

#plt.tight_layout()
# tr_D_now= tr_tracks_D[Stars_nr_D-1]
# print('Denis match time: ', tr_D_now.time[0], tr_D_now.time[-1])    
# #now is the mean time of the start and the end time of denis track
# datetime_now= tr_D_now.time[0] + (tr_D_now.time[-1]- tr_D_now.time[0])/2
# datetime_now= datetime_now.round('h')

# datetime_now= datetime_now + timedelta(hours= hourdiff)
# print('Datetime now: ', datetime_now)


# tr_D_now.plot_track(ax=ax, zorder= 1, color='g', label='Denis match')
# ind_now= np.where((tr_D_now.time == datetime_now).values)[0]
# if len(ind_now)== 1:
# #    track.plot_track(ax=ax, zorder= 1)   
#     tracklon, tracklat= tr_D_now.lonlat[ind_now[0]]
#     ax.scatter(tracklon, tracklat, transform=ccrs.PlateCarree(), s= 40, marker="s", color= 'g', zorder= 2)




# """match to STARS"""
# S = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv")
# #S_ind= STARS_individual_systems(S)

# S= S.rename(columns={"ID": "N"}) #"Obs": "row_idx",


# #all the PL IDs of now
# ID_now= remove_dublicate(S[np.abs((S.time - datetime_now)) < timedelta(hours= 12)].N)

# for ID in ID_now:
#     S_now= S[S['N']== ID]

#     ax.plot(S_now.lon, S_now.lat, transform=ccrs.PlateCarree(), lw= 1.5, color= 'r', zorder= 2, label='Rojo track')
#     ax.scatter(S_now.lon.values[0], S_now.lat.values[0], transform=ccrs.PlateCarree(), s= 40, color= 'r', zorder= 2)

#     i_now= np.argmin(np.abs(S_now.time -datetime_now))
    
#     ax.scatter(S_now.lon[i_now], S_now.lat[i_now], transform=ccrs.PlateCarree(), s= 40, marker= "s", color= 'r', zorder= 2)
#     ax.text(S_now.lon[i_now], S_now.lat[i_now]- 0.3, S_now.N[i_now], transform=ccrs.PlateCarree(), color='r', fontsize=13, fontweight= 'bold')
    
#     print('Rojo time now: ', S_now.time[i_now])


# """match to STARS- raw"""
# STARSdir= Mediadir+ "/PL/PLclim/STARS/"
# S_org= read_Gunnar_stars(STARSdir)

# ID_now= remove_dublicate(S_org[np.abs((S_org.time - datetime_now)) < timedelta(hours= 12)].N)

# for ID in ID_now:
#     S_now= S_org[S_org['N']== ID]

#     ax.plot(S_now.lon, S_now.lat, transform=ccrs.PlateCarree(), lw= 1.5, color= 'y', zorder= 2, label='Gunnar track')
#     ax.scatter(S_now.lon.values[0], S_now.lat.values[0], transform=ccrs.PlateCarree(), s= 40, color= 'y', zorder= 2)

#     i_now= np.argmin(np.abs(S_now.time -datetime_now))
#     ax.scatter(S_now.lon[i_now], S_now.lat[i_now], transform=ccrs.PlateCarree(), s= 40, marker= "s", color= 'y', zorder= 2)
    
#     print('Gunnar time now: ', S_now.time[i_now])

# #S_tracks= prepare_tracks(S_org)


# """my tracks"""
# if Plot_all_my_tracks:
#     test, fignr= 'test', 4
#     #test, fignr= 'test_1', 3
#     #test, fignr= 'Denis', 5
#     #test, fignr= 'version2', 7
#     test, fignr= 'version4', 6
    
#     #get the correct tracktime_id for which datetime_now is between startdate and enddate
    
#     file= "txt_files/track_list_PMC_appent.txt"
#     tracktime_id, startyear, startmonth, startday, endyear, endmonth, endday= np.loadtxt(file).T.astype(int)
#     startdate= np.array([datetime(startyear[i], startmonth[i], startday[i]) for i in range(len(startyear))])
#     enddate= np.array([datetime(endyear[i], endmonth[i], endday[i], 23) for i in range(len(startyear))])
    
    
    
#     tracktime_id_now= np.where(np.logical_and(datetime_now > startdate, datetime_now < enddate))[0]
#     if len(tracktime_id_now) > 0:
#         i = tracktime_id[tracktime_id_now ][0]
#         print('PMC timegroup: ',str(i))
#         track_dir= Mediadir+"/ERA5_STARS/tracks/"+test+"/tracks_"+str(i).zfill(3)
    
#         t_dir= Path(".") / track_dir
#         tr = TrackRun(t_dir, columns=['lon', 'lat', 'vo', 'time', 'area', 'vortex_type', 'slp'])
    
#         ilegend= 0
#         for (_, track) in tr.gb:
#             datetime_now_round= datetime_now.replace(hour= int(np.round(datetime_now.hour/3)*3))
#             ind_now= np.where((track.time == datetime_now_round).values)[0]
#             if len(ind_now)== 1:
#                 ilegend+=1
#                 if ilegend== 1: track.plot_track(ax=ax, zorder= 1, label='my tracks')
#                 else: track.plot_track(ax=ax, zorder= 1, label='')
                
#                 tracklon, tracklat= track.lonlat[ind_now[0]]
#                 ax.scatter(tracklon, tracklat, transform=ccrs.PlateCarree(), s= 50, marker="s", color= 'b', zorder= 2)
    
    
    
#     else: print('Not included in '+file)
    
    

# """trackmatch"""

# #
# #    
# #
# #S_tracks= prepare_tracks(S)


# plt.legend()
# plt.title('STARS nr: '+str(Stars_nr_D) +", Date: {}".format(datetime_now) )



# """plot ERA5 data"""
# filetype='tracking'
# plevel= 850
# var='vo'

# sym= True
# cmap= 'RdBu_r'

# cont_var='msl'


# year, month, day, hour = datetime_now.year, datetime_now.month, datetime_now.day, datetime_now.hour
# filetime= str(year)+'_'+str(month).zfill(2)+'_'+str(day).zfill(2)
# d0= xr.open_dataset(Mediadir + "ERA5_STARS/data/surface_era5_"+ filetime + '.nc')
# d0= d0.isel(time= hour)
# ds= xr.open_dataset(Mediadir + "ERA5_STARS/data/"+filetype+ "_era5_"+ filetime + '.nc')
# ds= ds.isel(time= hour)
# ds= xr.merge([d0, ds])


# if filetype == 'tracking':
#     ds['plev']/= 100
#     ds= ds.sel(plev= plevel)
    
# if var=='vo':
#     ds[var]*= 1E4
#     ds[var].attrs['units']= '10$^{-4}$ 1/s'   
    
# if sym== True:
# #    vextr= np.max([np.max(ds[var]), -np.min(ds[var])])
# #    cf= ax.pcolormesh(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
#     vextr= 6
#     cf= ax.contourf(ds.lon, ds.lat, ds[var], levels= [-4, -3, -2, -1, 1, 2, 3, 4, 6], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr, extend='max')


# cb= fig.colorbar(cf, ax= ax, shrink=0.8)
# varlabel= var
# if plevel != None: varlabel = var+'_'+str(plevel)
# cb.set_label(varlabel + ' ['+ ds[var].units + ']', size=14) 

# if 'msl' in list(ds.keys()):
#     ds['msl']/= 100
#     ds['msl'].attrs['units']= 'hPa'

# if cont_var== 'msl' :
#     cs= ax.contour(ds.lon, ds.lat, ds['msl'], np.arange(950, 1030, 2), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
#     plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')



# if save:
#     savefile= savedir+ 'Denis_nr_'+str(Stars_nr_D)+'_timediff'+str(hourdiff)
#     print(savefile)
#     plt.savefig(savefile , bbox_inches='tight') 
# #plt.tight_layout()