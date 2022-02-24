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
import cartopy.crs as ccrs


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



fignr= 6

fig = plt.figure(fignr)
plt.clf()
#ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (1,1,1))
#ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-45, 45, 50, 80], subplot= (1,1,1))
ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-40, 85, 50, 85], subplot= (1,1,1))


"""import the tracks"""
#track_dir= Mediadir+"/ERA5_STARS/tracks"
#track_dir= Mediadir+"/ERA5_STARS/tracks/conf"
#track_dir= Mediadir+"/ERA5_STARS/tracks/test"
track_dir= Mediadir+"/ERA5_STARS/tracks/version4/tracks_006"

t_dir= Path(".") / track_dir
#t_dir= Path(".") / "/media/pst019/1692A00D929FEF8B/ERA5_STARS/pmctrack-master/results/test6"

#tr= TrackRun(t_dir)
tr = TrackRun(t_dir, columns=['lon', 'lat', 'vo', 'time', 'area', 'vortex_type', 'slp'])
print(tr)
print(tr.columns)



"""conditions"""
conditions = [
    ("n", [lambda ot: ot.lifetime_h >= 3]),        
#    ("long_lived", [lambda ot: ot.lifetime_h >= 6]),
#    ("far_travelled_and_very_long_lived", [lambda ot: ot.lifetime_h >= 36, lambda ot: ot.gen_lys_dist_km > 300.0]),
#    ("strong", [lambda x: x.max_vort > 5e-4]),
]
tr.classify(conditions)



"""make dataframe"""
#trxar= tr.data.to_xarray()
#trdf= trxar.to_dataframe()

"""plot all tracks"""
#for (_, track) in tr.gb:
#    track.plot_track(ax=ax)
    
"""plot only some tracks"""
#for (_, track) in tr['long_lived'].gb:
#    track.plot_track(ax=ax)


year, month, day, hour= 1999, 12, 18, 12
datetime_now= np.datetime64(str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'T'+str(hour).zfill(2) )

for (_, track) in tr.gb:
#for (_, track) in tr['long_lived'].gb:   
    ind_now= np.where((track.time == datetime_now).values)[0]
    print(ind_now)
    if len(ind_now)== 1:
        track.plot_track(ax=ax, zorder= 1)
        
        tracklon, tracklat= track.lonlat[ind_now[0]]
        ax.scatter(tracklon, tracklat, transform=ccrs.PlateCarree(), s= 40, marker="s", color= 'b', zorder= 2)


"""list of octant tracks"""
tr_tracks= []
for (_, track) in tr.gb:
    tr_tracks.append(track)


"""plot ERA5 data"""
filetype='tracking'
plevel= 850
var='vo'

sym= True
cmap= 'RdBu_r'

cont_var='msl'


filetime= str(year)+'_'+str(month).zfill(2)+'_'+str(day).zfill(2)
d0= xr.open_dataset(Mediadir + "ERA5_STARS/data/surface_era5_"+ filetime + '.nc')
d0= d0.isel(time= hour)
ds= xr.open_dataset(Mediadir + "ERA5_STARS/data/"+filetype+ "_era5_"+ filetime + '.nc')
ds= ds.isel(time= hour)
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


"""match to STARS"""
S = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv")
#S_ind= STARS_individual_systems(S)


#all the PL IDs of now
ID_now= remove_dublicate(S[(S.time - datetime_now) < timedelta(hours= 3)].ID)

for ID in ID_now:
    S_now= S[S['ID']== ID]

    ax.plot(S_now.lon, S_now.lat, transform=ccrs.PlateCarree(), lw= 1.5, color= 'r', zorder= 2)
    
    i_now= np.argmin(np.abs(S_now.time -datetime_now))
    
    ax.scatter(S_now.lon[i_now], S_now.lat[i_now], transform=ccrs.PlateCarree(), s= 40, marker= "s", color= 'r', zorder= 2)




"""trackmatch"""
S= S.rename(columns={"ID": "N"}) #"Obs": "row_idx",


    
#from obs_tracks_api import prepare_tracks

S_tracks= prepare_tracks(S)


match= tr.match_tracks(S_tracks)

#get the track that matches
match_track= tr_tracks[match['n'][0][0] ]

match_track.plot_track(ax=ax, zorder= 1, color='g')