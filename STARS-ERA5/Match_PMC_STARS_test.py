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
from datetime import timedelta, datetime
import scipy.ndimage.filters as filters




fignr= 1
Plot=False

Plot=True
plot_iSTARS= 1 #only for Plot= True
plot_tracktime_id= 102


test= 'version4'
dist= 150 #merge STARS to PMC
merge_distance= 150 #merge two PMC that are close enough to STARS


outfile=Mediadir+"ERA5_STARS/STARS-PMC-merge/"+test+"_dist_"+str(dist)+"_test.csv"
write_i= 0 #running index to create the outfile only at first time
S_nr_for_tr_count= 0 #the full number of the Rojo PL for the Stoll number

"""get tracktime intervals - the time intervals in which PLs occur"""
#file= "txt_files/track_list_PMC.txt"
file= "txt_files/track_list_PMC_appent.txt" #used from version 4 - non-python length: 170

tracktime_id_list, startyear, startmonth, startday, endyear, endmonth, endday= np.loadtxt(file).T.astype(int)
#startdate_list= np.array([datetime(startyear[i], startmonth[i], startday[i]) for i in range(len(startyear))])
#enddate_list= np.array([datetime(endyear[i], endmonth[i], endday[i], 23) for i in range(len(startyear))])


"""get STARS list"""
S = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv", droplist= ['Optional_diameter', 'Comment_visual', 'Comment_uncertainties', 'U_3sec_kntots'])
#renaming for the trackmatch
S= S.rename(columns={"ID": "N"})

#round the time to nearest hour
S['time']= S.time.dt.round("H")





if Plot:
    fig = plt.figure(fignr)
    plt.clf()
    ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (1,1,1), scalebar=False)
    #ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-45, 45, 50, 80], subplot= (1,1,1))
    #ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-40, 85, 50, 85], subplot= (1,1,1))
    



"""loop through the tracktime_ids"""
if Plot: tr_id_start, tr_id_end= plot_tracktime_id-1, plot_tracktime_id
else: tr_id_start, tr_id_end= 0, len(tracktime_id_list)

for tracktime_id in tracktime_id_list[tr_id_start: tr_id_end]: #for testing
    ti= tracktime_id -1 # just for testing if only one case is taken
    print('tracktime id: ', tracktime_id)
    
    startdate= datetime(startyear[ti], startmonth[ti], startday[ti]) #startdate of this tracktime interval
    enddate= datetime(endyear[ti], endmonth[ti], endday[ti], 23)


    S_now_idlist= remove_dublicate(S[np.logical_and(S.time > startdate, S.time < enddate)].N) #all PLs from STARS that occur now
#    print('S_now: ', S_now_idlist)

    #get the PMC tracks from this tracktime interval
    track_dir= Mediadir+"/ERA5_STARS/tracks/"+test+"/tracks_"+str(tracktime_id).zfill(3)

    t_dir= Path(".") / track_dir
    tr = TrackRun(t_dir, columns=['lon', 'lat', 'vo', 'time', 'area', 'vortex_type', 'slp'])

    """conditions"""
    conditions = [ ("n", [lambda ot: ot.lifetime_h >= 1]),]
    tr.classify(conditions)

        
    """list of octant tracks"""
    tr_tracks= []
    ilegend= 0

    for (_, track) in tr.gb:
        tr_tracks.append(track)

       

    """loop through the STARS PLs that occur in this tracktime interval"""
#    for Si, S_now_id in enumerate(S_now_idlist):
    tr_count= 0
    for S_now_id in S_now_idlist:
        
        """sets the tr_count for the STOLL nr to 0 if a new full Rojo number is looped"""
        S_nr_for_tr_count_old= S_nr_for_tr_count
        S_nr_for_tr_count= int(S_now_id.split('.')[0])
        if S_nr_for_tr_count != S_nr_for_tr_count_old: tr_count= 0
        
        print('S_now_id: ', S_now_id)
        
        
        S_now= S[S['N']== S_now_id]
        
        S_now_track= prepare_tracks(S_now)
       
        if Plot:
            S_now_track[0].plot_track(ax=ax, color='r', label= 'STARS')
            ax.scatter(S_now_track[0].lon.values[plot_iSTARS], S_now_track[0].lat.values[plot_iSTARS], transform=ccrs.PlateCarree(), s= 50, marker= "s", color= 'r', zorder= 2)
            ax.text(S_now_track[0].lon.values[plot_iSTARS], S_now_track[0].lat.values[plot_iSTARS]+ 0.1, S_now_id, transform=ccrs.PlateCarree(), color='r', fontsize=13, fontweight= 'bold')

            time_plot_iSTARS= S_now_track[0].time.values[plot_iSTARS]
            
            for track in tr_tracks:
                #all tracks that intersect in time with the STARS track
                #they do not end before STARS PL start and do not start after the STARS PL has ended
                if (not (track.time.values[-1] < S_now_track[0].time.values[0])) and (not (track.time.values[0] > S_now_track[0].time.values[-1])):
                    ilegend+=1
                    if ilegend== 1:  track.plot_track(ax=ax, label= 'all tracks')
                    track.plot_track(ax=ax, label='')
                    
                    tr_i_now = np.argmin(np.abs(track.time - time_plot_iSTARS).values)
        #            print(tr_i_now)
                    tracklon, tracklat= track.lonlat[tr_i_now]
            #        if ilegend== 1: ax.scatter(tracklon, tracklat, transform=ccrs.PlateCarree(), s= 50, marker="s", color= 'b', zorder= 2, label= 'Location now')
                    ax.scatter(tracklon, tracklat, transform=ccrs.PlateCarree(), s= 50, marker="s", color= 'b', zorder= 2)
                    ax.text(tracklon, tracklat+ 0.1, track.index.labels[0][0], transform=ccrs.PlateCarree(), color='b', fontsize=13, fontweight= 'bold')




        match= tr.match_tracks(S_now_track)
        print('match: ', match)
        if len(match['n']) > 0:
            match_track= tr_tracks[match['n'][0][0] ]
            if Plot: match_track.plot_track(ax=ax, zorder= 3, color='g', label='simple match')
        
        matchbs= tr.match_tracks(S_now_track, method='bs2000') #, return_dist_matrix=True)
        print('match BS: ', matchbs)
        if len(matchbs['n']) > 0:
            match_track= tr_tracks[matchbs['n'][0][0] ]
            if Plot: match_track.plot_track(ax=ax, zorder= 3, color='purple', label='BS match')

#        match_track.plot_track(ax=ax, zorder= 1, color='g', label='simple match')
        


        """find all tracks that are within distance (150km) for at least one time step"""

        tr_match= []
        
        for t_S in range(len(S_now)):
            datetime_now= S_now.time.values[t_S]
            
            for track in tr_tracks: #maybe loop only over intersecting tracks
                i_track= np.where(track.time == datetime_now)[0]
                if len(i_track) > 0: #track occurs at this time
                    ind_track= i_track[0]
                    
                    dist_now= distance((S_now.lat.values[t_S], S_now.lon.values[t_S]), (track.lat.values[ind_track], track.lon.values[ind_track]) )
                    if dist_now <= dist:
                        print('S time', datetime_now,  'Track label', track.index.labels[0][0])
                        tr_match+= [track.index.labels[0][0] ]
                        
        tr_match= remove_dublicate(tr_match) #the tracks ids that match to STARS
        
        if Plot:
            for tr_match_nr in tr_match:
                tr_tracks[tr_match_nr].plot_track(ax=ax, color= 'orange', label= 'distance match', zorder= 3)
                        
#                    S_now.lon.values[t_S], S_now.lat.values[t_S],  track.lon[ind_track], track.lat[ind_track],
         
        """create a match matrix that labels the PMC tracks that can be merged by extrapolation in time"""
        match_matrix= np.zeros((len(tr_match), len(tr_match))) #matrix that label which of the tr_match ids are matching to each other
        for tr_i, tr_match_i in enumerate(tr_match):
            for tr_j, tr_match_j in enumerate(tr_match[:tr_i]):
                #check for time overlap - only match if they are not time overlapping
                overlaptime= np.intersect1d(tr_tracks[tr_match_i].time.values, tr_tracks[tr_match_j].time.values)
                if len(overlaptime) == 0:
                    #find which tr occurs first
                    first= np.argmin([tr_tracks[tr_match_i].time.max(), tr_tracks[tr_match_j].time.max()])
                    print('first (should be 1): ', first)
                    
                    if first==1: tr_match_i, tr_match_j= tr_match_j, tr_match_i #rotate such that i occurs before j
                        
                    time_diff= (tr_tracks[tr_match_j].time[0] - tr_tracks[tr_match_i].time[-1] ).seconds//3600 #the time difference in hours
                    
                    #extrapolate the first
                    if time_diff < 6: #only if the track is intersected for a few hours
                        if len(tr_tracks[tr_match_i]) > time_diff: #the track must at least be as long as the time_diff
                            
                            lat_extr, lon_extr= 2* tr_tracks[tr_match_i].lat[-1] - tr_tracks[tr_match_i].lat[-1 - time_diff], 2* tr_tracks[tr_match_i].lon[-1] - tr_tracks[tr_match_i].lon[-1 - time_diff]
                            if Plot: ax.scatter(lon_extr, lat_extr, transform=ccrs.PlateCarree(), s= 50, marker="v", color= 'orange', zorder= 2)
        
                            #compare distance to second
                            dist_extr= distance((lat_extr, lon_extr), (tr_tracks[tr_match_j].lat[0], tr_tracks[tr_match_j].lon[0]) )
                            if dist_extr < merge_distance: match_matrix[tr_i, tr_j]= 1
                        
                        #the extrapolation of the second backward in time
                        if len(tr_tracks[tr_match_j]) > time_diff: #the track must at least be as long as the time_diff

                            lat_extr, lon_extr= 2* tr_tracks[tr_match_j].lat[0] - tr_tracks[tr_match_j].lat[time_diff], 2* tr_tracks[tr_match_j].lon[0] - tr_tracks[tr_match_j].lon[time_diff]
                            if Plot: ax.scatter(lon_extr, lat_extr, transform=ccrs.PlateCarree(), s= 50, marker="^", color= 'orange', zorder= 2)
    
                            dist_extr= distance((lat_extr, lon_extr), (tr_tracks[tr_match_i].lat[-1], tr_tracks[tr_match_i].lon[-1]) )
                            if dist_extr < merge_distance: match_matrix[tr_i, tr_j]= 1
            
        
        
        match_matrix= match_matrix.T #to make mi occur before mj
        mi, mj = np.where(match_matrix == 1)
        for mx in range(len(mi)):            
            print('Merge: ', tr_match[mi[mx]], tr_match[mj[mx]])
        
        
        
        """write the csv file"""
#        print('STARS: ', S_now_id)
        tr_match_vector= [] #vector to the second PMC that occurs in the match_matix the same tr_count as the first
        start_obs_tr_match_vec = [] #the start observation of the Stoll track. normally it is 1, but if two tracks are matched it is one larger than the last number of the first part of the match
        
        for tr_i, tr_nr in enumerate(tr_match):
            
            """to make the index of the Stoll number, also for cases of matches"""
            PMCmatch= np.where(match_matrix[:, tr_i]== 1)[0]
            if len(PMCmatch) == 0: #no PMC match, so the track count is increased
                tr_count+=1
                tr_match_vector+= [tr_count] 
                start_obs_tr_match_vec += [1]
            elif len(PMCmatch) == 1:
                tr_match_vector+= [tr_match_vector[PMCmatch[0]]] #PMC match is put into the vector
                start_obs_tr_match_vec += [start_obs_tr_match_vec[PMCmatch[0]] +tr_tracks[tr_match[PMCmatch[0]]].index.values[-1][1] +1 ] #this is the row_idx +2 of the first part of the match
            else:
                print('more than two matches: have to adapt') #this has not been tested - since it does not occur
                break
                tr_match_vector+= [tr_match_vector[PMCmatch[-1]]] #PMC match is put into the vector
                start_obs_tr_match_vec += [start_obs_tr_match_vec[PMCmatch[-1]] + tr_tracks[tr_match[PMCmatch[-1]]].index.values[-1][1] +1 ] #this should be the start Obs of the previous part of the matched track + the row_idx +2 of the previous part of the match
                
            out_track= tr_tracks[tr_nr][:]
#            out_track.insert(0, "row_idx", out_track.index.labels[1])
            out_track.insert(0, "Rojo nr", S_now_id)
            out_track.insert(0, "track file", tracktime_id)

            out_track_pd_0= OctantTrack_to_df(out_track)
            out_track_pd= out_track_pd_0.set_index('time')
            out_track_pd.insert(0, "track_row_idx", out_track_pd_0.index.values)
            out_track_pd.insert(0, "Obs", start_obs_tr_match_vec[tr_i]+ out_track_pd_0.index.values )
#            out_track_pd.insert(0, "row_idx", list(out_track.index.labels[1])

            
            """get dataframe of STARS info for the tr_nr"""
            overlap_time= np.intersect1d(tr_tracks[tr_nr].time.values, S_now.time.values)
            STARS_now= S_now.loc[S_now['time'].isin(overlap_time)]
            STARS_now= STARS_now.set_index('time')
            STARS_now= STARS_now.loc[~STARS_now.index.duplicated(keep='first')] #if a time index in STARS occurs more than once (due to time rounding)
            STARS_now= STARS_now.drop(columns=['Season', 'Month', 'Stage'])
            STARS_now= STARS_now.rename(columns={"lat": "STARS lat", "lon": "STARS lon", "Obs": "STARS Obs nr", "Diameter": "Cloud diameter", "N": "Rojo nr old"})

            merged_track= pd.concat([out_track_pd, STARS_now], axis=1)
            
            merged_track.insert(0, "Stoll nr", S_now_id.split('.')[0]+ '_' + str(tr_match_vector[tr_i]))
            
            if write_i == 0:
                merged_track.to_csv(outfile)
                write_i += 1
            else:  merged_track.to_csv(outfile, mode='a', header=False)

                
plt.legend()










"""plot ERA5 data"""
if Plot:
    filetype='tracking'
    plevel= 850
    var='vo'
    
    sym= True
    cmap= 'RdBu_r'
    
    cont_var='msl'
    
    
    datetime_now= pd.DatetimeIndex([time_plot_iSTARS])[0]
    
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
    
#    if 'Denis' in test: apply_filter= False
#    if apply_filter:
#        variable= filters.uniform_filter(variable, size= (2, 3), mode='nearest')
        
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

plt.title(test +" {}".format(datetime_now))
