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
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/test_Denis/mc_era5-master/code') #to get obs_tracks_api

#import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from octant.core import TrackRun, OctantTrack
from pathlib import Path
from obs_tracks_api import prepare_tracks

from f_carto import *
from f_useful import *
from f_STARS import *
from datetime import timedelta
import scipy.ndimage.filters as filters




fignr= 4


#apply_filter = True #for the vorticity
apply_filter = False #for the vorticity
Plot_ERA5= True
#Plot_ERA5= False
#Plot_Steering=False
Plot_Steering=True

file= "txt_files/track_list_PMC_appent.txt"
tracktime_id, startyear, startmonth, startday, endyear, endmonth, endday= np.loadtxt(file).T.astype(int)

#test, fignr= 'test', 4
#test, fignr= 'test_smth40', 3
#test, fignr= 'test_nolsm', 2
test, fignr= 'test_1', 3
#test, fignr= 'Denis', 5
#test, fignr= 'version1', 6
#test, fignr= 'version1_nosmth', 8
#test, fignr= 'version2', 7
#test, fignr= 'version2_gamma01', 8
#test, fignr= 'version2_gamma01smth30', 9
test, fignr= 'version3', 6
#test, fignr= 'version3_rlink120', 8
test, fignr= 'version4', 6



#test, fignr= 'Denis_1h', 6


#print(tr.conf)

"""import the PMCtracks"""
i= 51 #this is non- python counting
#3 - there could be one more 
#7 is weird with high pressure - matching works
#8 ERA5 has several tracks, in sat images it looks like it not a multiple system..
#9 no PMC track even though vorticity is high  - tracking with smoothing detects 2 which is quite good compared to sat image
#10 - weak, not recognized in ERA5
#11 - something is wrong - matching detects wrong PMC. It does not merge the correct one to one PMC - interesting for see if algorithm works
#12 only works with ds2000 method for larger time
#13 -works fine
#14 - hopping a lot around .. test with smoothing 60km - gets a bit more straight, but still zig-zaggy
#15 - dual , again a weird track
#16 - split up, this should be reproduced

#track_dir= Mediadir+"/ERA5_STARS/tracks/tracks_"+str(tracktime_id[i]).zfill(3)
#track_dir= Mediadir+"/ERA5_STARS/tracks/tracks_ec_incl_"+str(tracktime_id[i]).zfill(3)
#track_dir= Mediadir+"/ERA5_STARS/tracks/tracks_smth"+str(tracktime_id[i]).zfill(3)'
i -=1
year, month, day, hour= startyear[i], startmonth[i], startday[i] , 23
datetime_now= np.datetime64(str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'T'+str(hour).zfill(2) )

datetime_now += np.timedelta64(1, 'D')
datetime_now= pd.to_datetime(datetime_now)
#plt.title("{} {:02d} {:02d} {:02d}UTC".format(year, month, day, hour))

track_dir= Mediadir+"/ERA5_STARS/tracks/"+test+"/tracks_"+str(tracktime_id[i]).zfill(3)
print(track_dir)


t_dir= Path(".") / track_dir
tr = TrackRun(t_dir, columns=['lon', 'lat', 'vo', 'time', 'area', 'vortex_type', 'slp'])


fig = plt.figure(fignr)
plt.clf()
ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (1,1,1), scalebar=False)
#ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-45, 45, 50, 80], subplot= (1,1,1))
#ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-40, 85, 50, 85], subplot= (1,1,1))

plt.title(test +" {}".format(datetime_now))



"""conditions"""
conditions = [ ("n", [lambda ot: ot.lifetime_h >= 3]),]
tr.classify(conditions)


    
"""plot ERA5 data"""
if Plot_ERA5:
    filetype='tracking'
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
    
    
    if filetype in ['plevels', 'vorticity', 'tracking']:
        ds['plev']/= 100
    #    ds= ds.sel(plev= plevel)
        
    if var=='vo':
        ds[var]*= 1E4
        ds[var].attrs['units']= '10$^{-4}$ 1/s'   
    
    variable= ds[var].sel(plev= plevel)    
    
    if 'Denis' in test: apply_filter= False
    if apply_filter:
        variable= filters.uniform_filter(variable, size= (2, 3), mode='nearest')
        
    if sym== True:
        vextr= np.max([np.max(variable), -np.min(variable)])
#        cf= ax.pcolormesh(ds.lon, ds.lat, variable, transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
#        ax.contour(ds.lon, ds.lat, variable, [1.5], transform= ccrs.PlateCarree(), colors= 'orange', linewidths= 1)
#        ax.contour(ds.lon, ds.lat, variable, [2.5], transform= ccrs.PlateCarree(), colors= 'orange', linewidths= 2)

        #much faster
#        cf= ax.contourf(ds.lon, ds.lat, variable, transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)

        vextr= 8
        cf= ax.contourf(ds.lon, ds.lat, variable, levels= [-4, -3, -2, -1, 1, 2, 3, 4, 6, 8, 12], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
        
#    
#    
#    
#    
#    
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
        


"""plot only some tracks"""
ilegend= 0
for (_, track) in tr.gb:
    ind_now= np.where((track.time == datetime_now).values)[0]
    if len(ind_now)== 1:
        ilegend+=1
        if ilegend== 1: track.plot_track(ax=ax, zorder= 1, label='all tracks now')
        else: track.plot_track(ax=ax, zorder= 1, label='')
    
        tracklon, tracklat= track.lonlat[ind_now[0]]
#        if ilegend== 1: ax.scatter(tracklon, tracklat, transform=ccrs.PlateCarree(), s= 50, marker="s", color= 'b', zorder= 2, label= 'Location now')
        ax.scatter(tracklon, tracklat, transform=ccrs.PlateCarree(), s= 50, marker="s", color= 'b', zorder= 2)


        steering_time= 1 #hours
        if Plot_ERA5 and Plot_Steering:
            radius= 200
            u_mean= value_in_rad(ds['u'].sel(plev= [1000, 700]).mean(axis= 0), ds.lat, ds.lon, tracklat, tracklon, radius)
            v_mean= value_in_rad(ds['v'].sel(plev= [1000, 700]).mean(axis= 0), ds.lat, ds.lon, tracklat, tracklon, radius)
    
#            if ilegend== 1: ax.scatter(tracklon+ (u_mean* 3.6*steering_time)/(110* np.cos(np.deg2rad(tracklat))), tracklat+ (v_mean* 3.6*3)/110, transform=ccrs.PlateCarree(), s= 50, marker="v", color= 'b', zorder= 2, label= 'Loc 3h steer r='+str(radius))
            est_lon= tracklon+ (u_mean* 3.6*steering_time)/(110* np.cos(np.deg2rad(tracklat))) #estimated longitude with steering
            est_lat=  tracklat+ (v_mean* 3.6*steering_time)/110
            ax.scatter(est_lon, est_lat, transform=ccrs.PlateCarree(), s= 50, marker="v", color= 'b', zorder= 2)
            plot_circle(ax, est_lon, est_lat, 120, edgecolor= 'b')
    
        if Plot_Steering:
            ax.text(tracklon, tracklat+ 0.1, track.index.labels[0][0], transform=ccrs.PlateCarree(), color='b', fontsize=13, fontweight= 'bold')
    
            ind_now= np.where((track.time == datetime_now + timedelta(hours= steering_time)).values)[0]
            if len(ind_now)== 1:   
                tracklon, tracklat= track.lonlat[ind_now[0]]
    #            if ilegend== 1: ax.scatter(tracklon, tracklat, transform=ccrs.PlateCarree(), s= 50, marker="x", color= 'b', zorder= 2, label= 'Location +3h')
                ax.scatter(tracklon, tracklat, transform=ccrs.PlateCarree(), s= 50, marker="x", color= 'b', zorder= 2)



"""list of octant tracks"""
tr_tracks= []
for (_, track) in tr.gb:
    tr_tracks.append(track)


#tr_tracks[15].plot_track(ax)



"""match to STARS"""
S = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv")
#renaming for the trackmatch
S= S.rename(columns={"ID": "N" }) #, "Obs": "row_idx")


#all the PL IDs of now
ID_now= remove_dublicate(S[np.abs((S.time - datetime_now)) < timedelta(hours= 12)].N)
print('STARS ID_now:', ID_now)


for ID in ID_now:
    S_now= S[S['N']== ID]

    ax.plot(S_now.lon, S_now.lat, transform=ccrs.PlateCarree(), lw= 1.5, color= 'r', zorder= 2, label= 'Rojo track')
    ax.scatter(S_now.lon.values[0], S_now.lat.values[0], transform=ccrs.PlateCarree(), s= 40, color= 'r', zorder= 2)
    
    i_now= np.argmin(np.abs(S_now.time -datetime_now))
    
    ax.scatter(S_now.lon[i_now], S_now.lat[i_now], transform=ccrs.PlateCarree(), s= 50, marker= "s", color= 'r', zorder= 2)


print('STARS now: ', S_now['time'])


"""trackmatch"""
#S_tracks= prepare_tracks(S)
#match= tr.match_tracks(S_tracks)

S_tracks= prepare_tracks(S_now)
#print(match.values()) #maybe this gives the distance matrix
#match= tr.match_tracks(S_tracks, method='intersection')
match= tr.match_tracks(S_tracks) #, method='intersection')


#get the track that matches
nr_matches= len(match['n'])
print('Matched PMC with STARS (simple): ', match['n'])


if nr_matches > 0:
    for ni in range(nr_matches):
        match_track= tr_tracks[match['n'][ni][0] ]
        match_track.plot_track(ax=ax, zorder= 1, color='g', label='simple match')

        print('Time of matched PMCtrack: \n', match_track.time[0], match_track.time[-1])
#
#
"""trackmatch with BS"""
matchbs= tr.match_tracks(S_tracks, method='bs2000') #, return_dist_matrix=True)
nr_matchesBS= len(matchbs['n'])
print('Matched PMC with STARS by BS: ', matchbs['n'])


if nr_matchesBS > 0:
    for ni in range(nr_matchesBS):
        
        match_track= tr_tracks[matchbs['n'][ni][0] ]
        match_track.plot_track(ax=ax, zorder= 1, color='purple', label='BS match')

        print('Time of matched PMCtrack: \n', match_track.time[0],  match_track.time[-1])

#        print('Radius of matched track: ', np.sqrt(match_track.area)* 4/np.pi)





plt.legend()
plt.tight_layout()