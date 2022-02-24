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
import pandas as pd

from octant.core import TrackRun, OctantTrack
from pathlib import Path
from obs_tracks_api import prepare_tracks

from f_carto import *
from f_useful import *
from f_STARS import *
from datetime import timedelta, datetime
import scipy.ndimage.filters as filters


fignr= 2

save=False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/Track_eval/'

"""version of my tracking"""
test= 'version4'
dist= 150


"""either specify the Rojo number and the index of the Rojo PL"""
specify='Rojonr'
Rojo_nr='6' #and all that occur simultaneously
R_Obs= 6 #Observation nr of the Rojo system, starts with 1

#or the datetime at which everything occurs - does not exactly work
#specify='datetime'
#year, month, day, hour= 2001, 2, 4, 8
#datetime_now= np.datetime64(str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'T'+str(hour).zfill(2) )

"""ERA 5 parameters"""
Plot_ERA5= True
#Plot_ERA5= False

filetype='tracking'
plevel= 850
var='vo'
#apply_filter = True #for the vorticity
apply_filter = False #for the vorticity

sym= True
cmap= 'RdBu_r'

cont_var='msl'


"""additional tracks"""
Denis=True #specify if the Denis track for this time should be plotted - often no Denis track




"""map"""
fig = plt.figure(fignr)
fignr+=1
plt.clf()
ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 55, 80], subplot= (1,1,1), scalebar=False)
#ax= Plot_Polar_Stereo(fig, central_longitude= 10, extent= [-10, 35, 68, 80], subplot= (1,1,1), scalebar=False)

#ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-45, 45, 50, 80], subplot= (1,1,1), scalebar=False)
#ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-40, 85, 50, 85], subplot= (1,1,1))


"""import Stoll matches"""
Stoll_STARS_file="ERA5_STARS/STARS-PMC-merge/"+test+"_dist_"+str(dist)+"_handfix.csv"
#Stoll_STARS_file="ERA5_STARS/STARS-PMC-merge/"+test+"_dist_"+str(dist)+".csv"

Stoll = pd.read_csv(Mediadir+ Stoll_STARS_file)
Stoll['time']= pd.DatetimeIndex(Stoll.time)
Stoll= Stoll_Obs_nr(Stoll) #get the Obs nr



"""Rojo PL"""
Rojo = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv")
Rojo= Rojo.rename(columns={"ID": "N"})
Rojo['time']= Rojo.time.dt.round("H") #this can lead to a Rojo PL having two values for one time 


if specify == 'Rojonr': #get datetime_now from the Rojonr and id
    datetime_now = Rojo[np.logical_and(Rojo['N'] == Rojo_nr, Rojo['Obs']== R_Obs)].time.values[0]
    #pd.DatetimeIndex(Rojo[np.logical_and(Rojo['N'] == Rojo_nr, Rojo['Obs']== i_R_now)].time.values)
    print('datatime now: ', datetime_now)


datetime_now= pd.DatetimeIndex([datetime_now])[0]
plt.title(test +" {}".format(datetime_now))
    
Rojo_nr_list= Rojo[Rojo['time']== datetime_now].N.values
Rojo_nr_list= remove_dublicate(Rojo_nr_list)
print('Rojo nr list: ', Rojo_nr_list)



""" exclude the dublicate - sometimes the Stoll systems matched to different Rojo PLs are the same """
excl_list_double=[] #this is the list of the Stoll matches that are double in the Stoll list
excl_dict_double={} #the Stoll nr of excluded system and the Rojo number that was assigned to it

for i_track_file in np.arange(np.min(Stoll['track file'].values), np.max(Stoll['track file'].values)+1):
    Stoll_trnr= Stoll[Stoll['track file'] == i_track_file]    

    Stoll_nrs_in_trnr= remove_dublicate(Stoll_trnr['Stoll nr'])
    
    d= {} #dictionary in which the Stoll nr and the corresponding tr_id_nr are saved
    for Stoll_nr in Stoll_nrs_in_trnr:
        d[Stoll_nr]= remove_dublicate(Stoll_trnr[Stoll_trnr['Stoll nr'] == Stoll_nr]['track_idx'])

    for i_S, Stoll_nr in enumerate(list(d.keys())):
        #if two Stoll_nr have a common member in the tr_id_nr it is added to the excl_list
        if common_member(list(d.values())[i_S], flat_list(list(d.values())[:i_S])):
            excl_list_double += [Stoll_nr]
            excl_dict_double[Stoll_nr]= remove_dublicate(Stoll[Stoll['Stoll nr'] == Stoll_nr]['Rojo nr'])[0]
 
print('Nr PLs excluded due to double occurrance: ', len(excl_list_double))


           
"""exclude systems with short lifetime"""
lifetime_threshold= 2 #tracks have to have at least this lifetime

Stoll_tracks= prepare_tracks(Stoll, system_id_label="Stoll nr")
excl_list_lifetime= [track['Stoll nr'][0] for track in Stoll_tracks if track.lifetime_h < lifetime_threshold]

print('Nr PLs excluded by lifetime: ', len(excl_list_lifetime))

"""exclude systems occuring mainly on land"""
lsm = xr.open_dataarray(Mediadir+ "ERA5_STARS/data/lsm.nc")
lsm = lsm.squeeze()  # remove singular time dimension

lsm.values[lsm > 0]= 1
#extremely unefficient
#Stoll['land']= [int(lsm.sel(latitude= Stoll.lat.values[i], longitude=Stoll.lon.values[i]) ) for i in range(len(Stoll))]

excl_list_land=[]
for track in Stoll_tracks:
    land_start=int(lsm.sel(latitude= track.lat.values[0], longitude=track.lon.values[0]) )
    land_final=int(lsm.sel(latitude= track.lat.values[-1], longitude=track.lon.values[-1]) )
    i_middle= int(track.lifetime_h//2)
    land_middle=int(lsm.sel(latitude= track.lat.values[i_middle], longitude=track.lon.values[i_middle]) )

    if land_start+ land_middle+ land_final== 3: excl_list_land += [track['Stoll nr'][0]]

print('Nr PLs excluded by land: ', len(excl_list_land))


"""systems that are existing already at the start of the file or still at the end of the file"""
#start of the file
file= "txt_files/track_list_PMC.txt"
tracktime_id, startyear, startmonth, startday, endyear, endmonth, endday= np.loadtxt(file).T.astype(int)
startdate= np.array([datetime(startyear[i], startmonth[i], startday[i]) for i in range(len(startyear))])
enddate= np.array([datetime(endyear[i], endmonth[i], endday[i], 23) for i in range(len(startyear))])

excl_dict_start, excl_dict_end = {}, {}
for track in Stoll_tracks:
    time0= np.datetime64(track.time.values[0], 's') #units ms to s
    time0=  datetime.utcfromtimestamp(time0.astype(int)) #convert np.datetime64 to datetime.datetime
    if time0 in startdate: excl_dict_start[track['Stoll nr'][0]]= time0

    timeend= np.datetime64(track.time.values[-1], 's') #units ms to s
    timeend=  datetime.utcfromtimestamp(timeend.astype(int)) #convert np.datetime64 to datetime.datetime
    
    if timeend in enddate: excl_dict_end[track['Stoll nr'][0]]= timeend



#exclude dublicate detections
excl_list= excl_list_double + excl_list_lifetime + excl_list_land       
Stoll_excl= Stoll[[Stoll_nr in excl_list for Stoll_nr in Stoll['Stoll nr']]]
Stoll= Stoll[[Stoll_nr not in excl_list for Stoll_nr in Stoll['Stoll nr']]]



"""plot all Rojo PLs that occur now"""
for Rojo_nr in Rojo_nr_list:
    
    Rojo_now= Rojo[Rojo['N'] == Rojo_nr]
    print('len Rojo: ', len(Rojo_now))

#    Rojo_now= Rojo_now[~Rojo_now.time.duplicated(keep='first')]
    
    R_tracks= prepare_tracks(Rojo_now, system_id_label="N")
    
    for track in R_tracks:
        track.plot_track(ax=ax, color='red', label='Rojo track')   
        
#        i_R_now= np.argmin(np.abs(track.time - datetime_now))
        i_R_now= (np.abs(track.time - datetime_now)).idxmin()

        ax.scatter(track.lon.values[i_R_now], track.lat.values[i_R_now], transform=ccrs.PlateCarree(), s= 50, marker= "s", color= 'r', zorder= 2)
        ax.text(track.lon.values[i_R_now], track.lat.values[i_R_now]+ 0.1, track.N[0], transform=ccrs.PlateCarree(), color='r', fontsize=13, fontweight= 'bold')
    
        time_plot_iRojo= track.time.values[i_R_now]
    
    


    
    """Stoll matches"""
    if Rojo_nr in list(excl_dict_double.values()):
        print('this PL was excluded due to double occurrence for Rojo nr ', Rojo_nr)
        #could try to save in exclusion which systems occurred simultaneously
    
    Stoll_now= Stoll[Stoll['Rojo nr']== Rojo_nr]
    
    S_tracks= prepare_tracks(Stoll_now, system_id_label="Stoll nr")
    
    for track in S_tracks:
        track.plot_track(ax=ax, color='b', label='Stoll match')    
    
        i_S_now= (np.abs(track.time - datetime_now)).idxmin()
        
        ax.scatter(track.lon.values[i_S_now], track.lat.values[i_S_now], transform=ccrs.PlateCarree(), s= 50, marker= "s", color= 'b', zorder= 2)
        ax.text(track.lon.values[i_S_now], track.lat.values[i_S_now]- 0.3, track['Stoll nr'][0], transform=ccrs.PlateCarree(), color='b', fontsize=13, fontweight= 'bold')

    """the excluded Stoll matches"""
    Stoll_now= Stoll_excl[Stoll_excl['Rojo nr']== Rojo_nr]
    
    S_tracks= prepare_tracks(Stoll_now, system_id_label="Stoll nr")
    
    for track in S_tracks:
        track.plot_track(ax=ax, color='cornflowerblue', label='Stoll excluded')    
    
#        i_S_now= np.argmin(np.abs(track.time - datetime_now))
        i_S_now= (np.abs(track.time - datetime_now)).idxmin()       
        ax.scatter(track.lon.values[i_S_now], track.lat.values[i_S_now], transform=ccrs.PlateCarree(), s= 50, marker= "s", color='cornflowerblue', zorder= 2)
        ax.text(track.lon.values[i_S_now], track.lat.values[i_S_now]- 0.3, track['Stoll nr'][0], transform=ccrs.PlateCarree(), color='cornflowerblue', fontsize=13, fontweight= 'bold')




"""match to STARS- raw"""
Gunnardir= Mediadir+ "/PL/PLclim/STARS/"
Gunnar_org= read_Gunnar_stars(Gunnardir)

ID_now= remove_dublicate(Gunnar_org[np.abs((Gunnar_org.time - datetime_now)) < timedelta(hours= 12)].N)

for ID in ID_now:
    Gunnar_now= Gunnar_org[Gunnar_org['N']== ID]

    ax.plot(Gunnar_now.lon, Gunnar_now.lat, transform=ccrs.PlateCarree(), lw= 1.5, color= 'y', zorder= 2, label='Gunnar track')
    ax.scatter(Gunnar_now.lon.values[0], Gunnar_now.lat.values[0], transform=ccrs.PlateCarree(), s= 40, color= 'y', zorder= 2)

#    i_now= np.argmin(np.abs(Gunnar_now.time -datetime_now))
    i_now= (np.abs(Gunnar_now.time -datetime_now)).idxmin()
    
    ax.scatter(Gunnar_now.lon[i_now], Gunnar_now.lat[i_now], transform=ccrs.PlateCarree(), s= 40, marker= "s", color= 'y', zorder= 2)
    
    Gunnar_nr= Gunnar_now.N.iloc[0]
    print('Gunnar time now: ', Gunnar_now.time[i_now])







"""file from Denis"""
if Denis:
    filename = "txt_files/era5_run000_2000_2018__matched_to_stars_6h__bs2000_beta100.h5"
    track_D = TrackRun.from_archive(filename)
    #print(track_D.data)
    
    Denis_tracks= []
    Denis_array= np.zeros((0,8))
    for (_, track) in track_D.gb:
        Denis_tracks.append(track)
        Denis_array=np.vstack([Denis_array, track.values])
    
    columns= Denis_tracks[0].columns.values
    
    i_Denis_now= np.where(Denis_array[:, 3] == datetime_now)[0]
    if len(i_Denis_now)== 0: print('no Denis track at this time')
    elif len(i_Denis_now) >1: print('more than one Denis track')
    else:
        i_Denis_now= i_Denis_now[0]
        Denis_nr= Denis_array[i_Denis_now][-1]
        Denis_now= Denis_array[Denis_array[:, -1]== Denis_nr]
        Denis_now= pd.DataFrame(Denis_now, columns=columns)
        Denis_track_now= prepare_tracks(Denis_now, system_id_label="STARS_N")
        Denis_track_now= Denis_track_now[0]
        
        print('Denis match time: ', Denis_track_now.time[0], Denis_track_now.time.iloc[-1])    

        
        Denis_track_now.plot_track(ax=ax, zorder= 1, color='g', label='Denis match')
        ind_now= np.where((Denis_track_now.time == datetime_now).values)[0]
        if len(ind_now)== 1:
        #    track.plot_track(ax=ax, zorder= 1)   
            tracklon, tracklat= Denis_track_now.lonlat[ind_now[0]]
            ax.scatter(tracklon, tracklat, transform=ccrs.PlateCarree(), s= 40, marker="s", color= 'g', zorder= 2)
    





"""plot ERA5 data"""
if Plot_ERA5:    
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

        vextr= 6
#        cf= ax.contourf(ds.lon, ds.lat, variable, levels= [-4, -3, -2, -1, 1, 2, 3, 4, 6, 8, 12], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr)
        cf= ax.contourf(ds.lon, ds.lat, variable, levels= [-3, -2, -1, 1, 2, 3, 4, 6], transform= ccrs.PlateCarree(), cmap= cmap, vmin= -vextr, vmax= vextr, extend='max')
        
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
        plt.clabel(cs, fontsize=10, fmt='%1.0f')
        


plt.legend()

if save:
    print(savedir)
    plt.savefig(savedir+ specify +Rojo_nr +'_Obsnr'+str(R_Obs), bbox_inches='tight') 
#plt.tight_layout()



"""plot the time evolution of the systems"""
plt.figure(fignr, figsize=(6,8))
fignr +=1
plt.clf()


ax1= plt.subplot(511)
ax2= plt.subplot(512, sharex= ax1)
ax3= plt.subplot(513, sharex= ax1)
ax4= plt.subplot(514, sharex= ax1)
ax5= plt.subplot(515, sharex= ax1) #for each Stoll low, the cloud morph, number and distance to the Rojo PL is printed


i= -1 #index for the cloud morphology of the Rojo low connected to each Stoll low
for Rojo_nr in Rojo_nr_list:
    Stoll_now= Stoll[Stoll['Rojo nr']== Rojo_nr]
    
    list_Stoll_nr= remove_dublicate(list(Stoll_now['Stoll nr']))
    
    for Stoll_nr in list_Stoll_nr:
        color = next(ax1._get_lines.prop_cycler)['color']
        i+=1
        Stoll_track= Stoll_now[Stoll_now['Stoll nr']== Stoll_nr]
    
        ax1.plot(Stoll_track.time, Stoll_track.vo*1e5, label= Stoll_nr, color=color)
        ax2.plot(Stoll_track.time, np.sqrt(Stoll_track.area)*4/np.pi, label= Stoll_nr +' Vort', color=color)
        ax2.plot(Stoll_track.time, Stoll_track['Cloud diameter']/2, 'x', label= Rojo_nr +' Cloud', color=color)
        
        print('Rojo nr: ', Rojo_nr, remove_dublicate(Stoll_track['Comment'].dropna()))
    
        ax3.plot(Stoll_track.time, Stoll_track.slp, label= Stoll_nr, color=color)
        ax4.plot(Stoll_track.time, Stoll_track.vortex_type, label=Stoll_nr, color=color)
        
        Stoll_track_Rdata= Stoll_track[~Stoll_track.Morphology.isna()] #get only the parts of the Stoll track where Rojo data is also available
        for ti in range(len(Stoll_track_Rdata)):       
            dist= distance((Stoll_track_Rdata.lat.values[ti], Stoll_track_Rdata.lon.values[ti]),
                           (Stoll_track_Rdata['STARS lat'].values[ti], Stoll_track_Rdata['STARS lon'].values[ti]))
            ax5.text(Stoll_track_Rdata.time.values[ti], i,
                     Stoll_track_Rdata.Morphology.values[ti] +'\n'+ Stoll_track_Rdata['Rojo nr'].values[ti] +'\n'+ str(int(dist)),
                     color=color, verticalalignment='center', horizontalalignment='center')


#        ax3.plot(Stoll_track.time, Stoll_track['Press_min'], 'x', label= Rojo_nr)
#
    
ax1.set_ylabel('Vorticity [10$^{-5}$s$^{-1}$]')
ax2.set_ylabel('Radius [km]')
ax3.set_ylabel('SLP min [hPa]')
ax4.set_ylabel('Vortex type')
ax5.set_ylabel('Cloud morphology')
ax5.set_ylim([-.5, i+ 0.5]) #on the axis should be the Stoll nr

import matplotlib.dates as mdates
#xfmt = mdates.DateFormatter('%d-%m-%y %H:%M')
xfmt = mdates.DateFormatter('%d.%m %H:00')

ax1.tick_params(axis='x', labelbottom=False)
ax2.tick_params(axis='x', labelbottom=False)
ax3.tick_params(axis='x', labelbottom=False)
ax4.tick_params(axis='x', labelbottom=False)

ax5.xaxis.set_major_formatter(xfmt)
plt.setp(ax5.get_xticklabels(), rotation=15, horizontalalignment='center')

ax2.legend(ncol=2)
ax1.legend()


if save:
    print(savedir)
    plt.savefig(savedir+ specify +Rojo_nr +'_param-evol', bbox_inches='tight') 



plt.tight_layout()
