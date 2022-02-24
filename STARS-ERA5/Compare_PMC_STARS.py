#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 16:29:52 2019

@author: pst019
"""


#import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from octant.core import TrackRun, OctantTrack
from pathlib import Path

from f_carto import *
from f_useful import *
from f_STARS import *

import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    homedir= '/home/'+user+'/home/'
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/1692A00D929FEF8B/'
'/home/'+user+'/home/'



fignr= 4

fig = plt.figure(fignr)
plt.clf()
ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (1,1,1))
#ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-45, 45, 50, 80], subplot= (1,1,1))
#ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-40, 85, 50, 85], subplot= (1,1,1))


file= "txt_files/track_list_PMC.txt"
tracktime_id, startyear, startmonth, startday, endyear, endmonth, endday= np.loadtxt(file).T.astype(int)


"""import the PMCtracks"""
i= 11 #2 - there could be one more #this is python counting
#6 is weird with high pressure - matching works
#7 ERA5 has several tracks, in sat images it looks like it not a multiple system..
#8 no PMC track even though vorticity is high
#9 - weak, not recognized in ERA5
#10 - matching detects wrong PMC. It does not merge the correct one to one PMC - interesting for see if algorithm works
#11 only works with ds2000 method for larger time


#track_dir= Mediadir+"/ERA5_STARS/tracks/tracks"+str(tracktime_id[i]).zfill(3)
track_dir= Mediadir+"/ERA5_STARS/tracks/tracks_h_"+str(tracktime_id[i]).zfill(3)
#track_dir= Mediadir+"/ERA5_STARS/tracks/tracks_smth"+str(tracktime_id[i]).zfill(3)

t_dir= Path(".") / track_dir
tr = TrackRun(t_dir, columns=['lon', 'lat', 'vo', 'time', 'area', 'vortex_type', 'slp'])


"""conditions"""
conditions = [ ("n", [lambda ot: ot.lifetime_h >= 3]),]
tr.classify(conditions)


    
"""plot only some tracks"""
year, month, day, hour= startyear[i], startmonth[i], startday[i]+1, 6

datetime_now= np.datetime64(str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'T'+str(hour).zfill(2) )
datetime_now += np.timedelta64(1, 'D')
datetime_now= pd.to_datetime(datetime_now)
#plt.title("{} {:02d} {:02d} {:02d}UTC".format(year, month, day, hour))
plt.title("{}".format(datetime_now))
 
for (_, track) in tr.gb:
    ind_now= np.where((track.time == datetime_now).values)[0]
    if len(ind_now)== 1:
        track.plot_track(ax=ax, zorder= 1)
        
        tracklon, tracklat= track.lonlat[ind_now[0]]
        ax.scatter(tracklon, tracklat, transform=ccrs.PlateCarree(), s= 50, marker="s", color= 'b', zorder= 2)


"""list of octant tracks"""
tr_tracks= []
for (_, track) in tr.gb:
    tr_tracks.append(track)





"""match to STARS"""
S = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv")
#renaming for the trackmatch
S= S.rename(columns={"ID": "N"}) #"Obs": "row_idx",



#all the PL IDs of now
ID_now= remove_dublicate(S[np.abs((S.time - datetime_now)) < timedelta(hours= 12)].N)
print('STARS ID_now:', ID_now)


for ID in ID_now:
    S_now= S[S['N']== ID]

    ax.plot(S_now.lon, S_now.lat, transform=ccrs.PlateCarree(), lw= 1.5, color= 'r', zorder= 2)
    
    i_now= np.argmin(np.abs(S_now.time -datetime_now))
    
    ax.scatter(S_now.lon[i_now], S_now.lat[i_now], transform=ccrs.PlateCarree(), s= 50, marker= "s", color= 'r', zorder= 2)


print('STARS now: ', S_now['time'])


"""trackmatch"""
#S_tracks= prepare_tracks(S)
#match= tr.match_tracks(S_tracks)

S_tracks= prepare_tracks(S_now)
match= tr.match_tracks(S_tracks, method='bs2000') #, return_dist_matrix=True)
#print(match.values()) #maybe this gives the distance matrix
#match= tr.match_tracks(S_tracks, method='intersection')


#get the track that matches
nr_matches= len(match['n'])
print('Matched PMC with STARS: ', match['n'])


if nr_matches > 0:
    for ni in range(nr_matches):
        match_track= tr_tracks[match['n'][ni][0] ]
        match_track.plot_track(ax=ax, zorder= 1, color='g')

        print('Time of matched PMCtrack: ', match_track.time)

        print('Radius of matched track: ', np.sqrt(match_track.area)* 4/np.pi)

"""plot ERA5 data"""
filetype='vorticity'
plevel= 850
var='vo'

sym= True
cmap= 'RdBu_r'

cont_var='msl'


filetime= str(datetime_now.year)+'_'+str(datetime_now.month).zfill(2)+'_'+str(datetime_now.day).zfill(2)
d0= xr.open_dataset(Mediadir + "ERA5_STARS/data/surface_era5_"+ filetime + '.nc')
d0= d0.isel(time= datetime_now.hour)
ds= xr.open_dataset(Mediadir + "ERA5_STARS/data/"+filetype+ "_era5_"+ filetime + '.nc')
ds= ds.isel(time= datetime_now.hour)
ds= xr.merge([d0, ds])


if filetype in ['plevels', 'vorticity']:
    ds['plev']/= 100
    ds= ds.sel(plev= plevel)
    
if var=='vo':
    ds[var]*= 1E5
    ds[var].attrs['units']= '10$^{-5}$ 1/s'   
    
if sym== True:
    vextr= np.max([np.max(ds[var]), -np.min(ds[var])])
    cf= ax.pcolormesh(ds.lon, ds.lat, ds[var], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
    
    ax.contour(ds.lon, ds.lat, ds[var], [15], transform= ccrs.PlateCarree(), colors= 'r', linewidths= 1)
    ax.contour(ds.lon, ds.lat, ds[var], [25], transform= ccrs.PlateCarree(), colors= 'r', linewidths= 2)



cb= fig.colorbar(cf, ax= ax, shrink=0.8)
varlabel= var
if plevel != None: varlabel = var+'_'+str(plevel)
cb.set_label(varlabel + ' ['+ ds[var].units + ']', size=14) 

if 'msl' in list(ds.keys()):
    ds['msl']/= 100
    ds['msl'].attrs['units']= 'hPa'

if cont_var== 'msl' :
    cs= ax.contour(ds.lon, ds.lat, ds['msl'], np.arange(950, 1050, 2), transform= ccrs.PlateCarree(), colors='k', linewidths= 1)
    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')
    
