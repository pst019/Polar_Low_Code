#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 16:29:52 2019

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    homedir= '/home/'+user+'/home/'
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/1692A00D929FEF8B/'
'/home/'+user+'/home/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cartopy.crs as ccrs

from octant.core import TrackRun, OctantTrack
from pathlib import Path

from f_carto import *
from f_useful import *
from f_STARS import *

from datetime import datetime, timedelta

from obs_tracks_api import prepare_tracks



fignr= 1

save=False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/Track_eval/'

Stars_nr_D= 2
#6 nice case to illustrate difference

###it plots the average time of the denis track
hourdiff= 0 # how much before or after the average time the plot is done

Plot_all_my_tracks=False


"""map"""
fig = plt.figure(fignr)
plt.clf()
ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (1,1,1), scalebar=False)
#ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-45, 45, 50, 80], subplot= (1,1,1))
#ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-40, 85, 50, 85], subplot= (1,1,1))


"""file from Denis"""
filename = "txt_files/era5_run000_2000_2018__matched_to_stars_6h__bs2000_beta100.h5"
track_D = TrackRun.from_archive(filename)
#print(track_D.data)

tr_tracks_D= []
for (_, track) in track_D.gb:
    tr_tracks_D.append(track)
    



tr_D_now= tr_tracks_D[Stars_nr_D-1]
print('Denis match time: ', tr_D_now.time[0], tr_D_now.time[-1])    
#now is the mean time of the start and the end time of denis track
datetime_now= tr_D_now.time[0] + (tr_D_now.time[-1]- tr_D_now.time[0])/2
datetime_now= datetime_now.round('h')

datetime_now= datetime_now + timedelta(hours= hourdiff)
print('Datetime now: ', datetime_now)


tr_D_now.plot_track(ax=ax, zorder= 1, color='g', label='Denis match')
ind_now= np.where((tr_D_now.time == datetime_now).values)[0]
if len(ind_now)== 1:
#    track.plot_track(ax=ax, zorder= 1)   
    tracklon, tracklat= tr_D_now.lonlat[ind_now[0]]
    ax.scatter(tracklon, tracklat, transform=ccrs.PlateCarree(), s= 40, marker="s", color= 'g', zorder= 2)




"""match to STARS"""
S = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv")
#S_ind= STARS_individual_systems(S)

S= S.rename(columns={"ID": "N"}) #"Obs": "row_idx",


#all the PL IDs of now
ID_now= remove_dublicate(S[np.abs((S.time - datetime_now)) < timedelta(hours= 12)].N)

for ID in ID_now:
    S_now= S[S['N']== ID]

    ax.plot(S_now.lon, S_now.lat, transform=ccrs.PlateCarree(), lw= 1.5, color= 'r', zorder= 2, label='Rojo track')
    ax.scatter(S_now.lon.values[0], S_now.lat.values[0], transform=ccrs.PlateCarree(), s= 40, color= 'r', zorder= 2)

    i_now= np.argmin(np.abs(S_now.time -datetime_now))
    
    ax.scatter(S_now.lon[i_now], S_now.lat[i_now], transform=ccrs.PlateCarree(), s= 40, marker= "s", color= 'r', zorder= 2)
    ax.text(S_now.lon[i_now], S_now.lat[i_now]- 0.3, S_now.N[i_now], transform=ccrs.PlateCarree(), color='r', fontsize=13, fontweight= 'bold')
    
    print('Rojo time now: ', S_now.time[i_now])


"""match to STARS- raw"""
STARSdir= Mediadir+ "/PL/PLclim/STARS/"
S_org= read_Gunnar_stars(STARSdir)

ID_now= remove_dublicate(S_org[np.abs((S_org.time - datetime_now)) < timedelta(hours= 12)].N)

for ID in ID_now:
    S_now= S_org[S_org['N']== ID]

    ax.plot(S_now.lon, S_now.lat, transform=ccrs.PlateCarree(), lw= 1.5, color= 'y', zorder= 2, label='Gunnar track')
    ax.scatter(S_now.lon.values[0], S_now.lat.values[0], transform=ccrs.PlateCarree(), s= 40, color= 'y', zorder= 2)

    i_now= np.argmin(np.abs(S_now.time -datetime_now))
    ax.scatter(S_now.lon[i_now], S_now.lat[i_now], transform=ccrs.PlateCarree(), s= 40, marker= "s", color= 'y', zorder= 2)
    
    print('Gunnar time now: ', S_now.time[i_now])

#S_tracks= prepare_tracks(S_org)


"""my tracks"""
if Plot_all_my_tracks:
    test, fignr= 'test', 4
    #test, fignr= 'test_1', 3
    #test, fignr= 'Denis', 5
    #test, fignr= 'version2', 7
    test, fignr= 'version4', 6
    
    #get the correct tracktime_id for which datetime_now is between startdate and enddate
    
    file= "txt_files/track_list_PMC_appent.txt"
    tracktime_id, startyear, startmonth, startday, endyear, endmonth, endday= np.loadtxt(file).T.astype(int)
    startdate= np.array([datetime(startyear[i], startmonth[i], startday[i]) for i in range(len(startyear))])
    enddate= np.array([datetime(endyear[i], endmonth[i], endday[i], 23) for i in range(len(startyear))])
    
    
    
    tracktime_id_now= np.where(np.logical_and(datetime_now > startdate, datetime_now < enddate))[0]
    if len(tracktime_id_now) > 0:
        i = tracktime_id[tracktime_id_now ][0]
        print('PMC timegroup: ',str(i))
        track_dir= Mediadir+"/ERA5_STARS/tracks/"+test+"/tracks_"+str(i).zfill(3)
    
        t_dir= Path(".") / track_dir
        tr = TrackRun(t_dir, columns=['lon', 'lat', 'vo', 'time', 'area', 'vortex_type', 'slp'])
    
        ilegend= 0
        for (_, track) in tr.gb:
            datetime_now_round= datetime_now.replace(hour= int(np.round(datetime_now.hour/3)*3))
            ind_now= np.where((track.time == datetime_now_round).values)[0]
            if len(ind_now)== 1:
                ilegend+=1
                if ilegend== 1: track.plot_track(ax=ax, zorder= 1, label='my tracks')
                else: track.plot_track(ax=ax, zorder= 1, label='')
                
                tracklon, tracklat= track.lonlat[ind_now[0]]
                ax.scatter(tracklon, tracklat, transform=ccrs.PlateCarree(), s= 50, marker="s", color= 'b', zorder= 2)
    
    
    
    else: print('Not included in '+file)
    
    

"""trackmatch"""

#
#    
#
#S_tracks= prepare_tracks(S)


plt.legend()
plt.title('STARS nr: '+str(Stars_nr_D) +", Date: {}".format(datetime_now) )



"""plot ERA5 data"""
filetype='tracking'
plevel= 850
var='vo'

sym= True
cmap= 'RdBu_r'

cont_var='msl'


year, month, day, hour = datetime_now.year, datetime_now.month, datetime_now.day, datetime_now.hour
filetime= str(year)+'_'+str(month).zfill(2)+'_'+str(day).zfill(2)
d0= xr.open_dataset(Mediadir + "ERA5_STARS/data/surface_era5_"+ filetime + '.nc')
d0= d0.isel(time= hour)
ds= xr.open_dataset(Mediadir + "ERA5_STARS/data/"+filetype+ "_era5_"+ filetime + '.nc')
ds= ds.isel(time= hour)
ds= xr.merge([d0, ds])


if filetype == 'tracking':
    ds['plev']/= 100
    ds= ds.sel(plev= plevel)
    
if var=='vo':
    ds[var]*= 1E4
    ds[var].attrs['units']= '10$^{-4}$ 1/s'   
    
if sym== True:
#    vextr= np.max([np.max(ds[var]), -np.min(ds[var])])
#    cf= ax.pcolormesh(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
    vextr= 6
    cf= ax.contourf(ds.lon, ds.lat, ds[var], levels= [-4, -3, -2, -1, 1, 2, 3, 4, 6], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr, extend='max')


cb= fig.colorbar(cf, ax= ax, shrink=0.8)
varlabel= var
if plevel != None: varlabel = var+'_'+str(plevel)
cb.set_label(varlabel + ' ['+ ds[var].units + ']', size=14) 

if 'msl' in list(ds.keys()):
    ds['msl']/= 100
    ds['msl'].attrs['units']= 'hPa'

if cont_var== 'msl' :
    cs= ax.contour(ds.lon, ds.lat, ds['msl'], np.arange(950, 1030, 2), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')



if save:
    savefile= savedir+ 'Denis_nr_'+str(Stars_nr_D)+'_timediff'+str(hourdiff)
    print(savefile)
    plt.savefig(savefile , bbox_inches='tight') 
#plt.tight_layout()