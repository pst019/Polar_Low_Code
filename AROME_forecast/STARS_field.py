#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 10:01:45 2019

@author: pst019
"""


import time
start = time.time()


import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'

Mediadir= '/run/media/pst019/1692A00D929FEF8B/'
homedir= Mediadir+'home/'


import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from f_useful import *



#save= True
save= False

fignr= 1


S = pd.read_csv(Mediadir+"PL/STARS/Rojo-etal_2019.csv", sep=',', header= 27)


new_cols= list(S.columns)
new_cols[1] = 'Time'
new_cols[6:8] = ['Diameter', 'Optional_diameter']
new_cols[10:16] = ['Comment', 'Comment_visual', 'Comment_uncertainties', 'Press_min', 'U_10min_knots', 'U_3sec_kntots']
S.columns= new_cols

S.drop(['Optional_diameter', 'Comment', 'Comment_visual', 'Comment_uncertainties', 'Press_min', 'U_10min_knots', 'U_3sec_kntots'],axis=1, inplace=True )

S['datetime'] = pd.to_datetime(S['Time'])
S.drop(['Time'], axis=1, inplace=True)

S['Season']= S['Season'].str.slice(0, 4).astype(int)
S['Month']= S['datetime'].dt.month
#S = S.set_index('datetime')
S['ID']= S['ID'].replace({'98': '97'}) #98 is according to Rojo a continuation of 97


"""get individual systems"""
S_ind= S.groupby(['ID']).mean()
S_ind['ID']= S_ind.index


"""system morphology"""
##S = Spiraliform, C = Comma shaped, MGR = Merry-go-round, W = Wave system, U = Undefined, T = Transition between different forms, H = Hybrid, - = Systems don't appear completely on imagery

S= S.replace({'comma': 'C', 'undefined': 'U'}) # some wrong morphologies
remove_list= ['-', 'T', 'U', 'H'] #maybe should not remove hyrbid


S['Morph_red']= S['Morphology'].replace({'T - ': '', ' - T': '', '  ': ' ', 'U - ': '', ' - U': '', ' - H': '', 'H - ': ''}, regex= True)
S_ind['PL_Morph_full']= [" ".join(remove_repeater(split_items(remove_repeater(remove_from_list(S[S['ID'] == ID]['Morph_red'], remove_list)), ' - '))) for ID in S_ind['ID']]

#splits transition, remove dublicates, sort by alphabet
S['Morph_red']= S['Morphology'].replace({'  ': ' '}, regex= True)
S_ind['PL_Morph']= [" ".join(sorted(remove_dublicate(remove_from_list(split_items(S[S['ID']== ID]['Morph_red'], ' - '), remove_list)))) for ID in S_ind['ID']]

S_ind.loc[S_ind['PL_Morph'].str.contains('MGR'), 'PL_Morph']= 'MGR+'
S_ind.loc[S_ind['PL_Morph'].str.contains('W'), 'PL_Morph']= 'W+'
S_ind['PL_Morph']= S_ind['PL_Morph'].replace({'C S': 'C-S'})


"""full morphology"""

#S_ind.loc[S_ind['PL_Morph_full'].str.contains(' MGR '), 'PL_Morph_full']= '..MGR..'
S_ind.loc[S_ind['PL_Morph_full'].str.contains(' MGR'), 'PL_Morph_full']= '.. MGR'
S_ind.loc[S_ind['PL_Morph_full'].str.contains('MGR '), 'PL_Morph_full']= 'MGR ..'

#S_ind.loc[np.logical_and.reduce((S_ind['PL_Morph_full'] != 'C W C' , S_ind['PL_Morph_full'] != 'C W S' , S_ind['PL_Morph_full'].str.contains(' W '))), 'PL_Morph_full']= '..W..'
#S_ind.loc[np.logical_and(~S_ind['PL_Morph_full'].str.contains(' W ')  , S_ind['PL_Morph_full'].str.contains(' W')), 'PL_Morph_full']= '..W'

#S_ind.loc[S_ind['PL_Morph_full'].str.contains(' W '), 'PL_Morph_full']= '..W..'
#S_ind.loc[S_ind['PL_Morph_full'].str.contains('W '), 'PL_Morph_full']= 'W..'
#S_ind.loc[S_ind['PL_Morph_full'].str.contains(' W'), 'PL_Morph_full']= '..W'

#S_ind.loc[S_ind['PL_Morph_full'].str.contains('C S C'), 'PL_Morph_full']= 'CSC+'
#S_ind.loc[np.logical_or(S_ind['PL_Morph_full'] == 'S C S' , S_ind['PL_Morph_full'].str.contains('S C S C')), 'PL_Morph_full']= '. S C S ..'
#S_ind.loc[S_ind['PL_Morph_full'].str.contains('S C S'), 'PL_Morph_full']= 'CSC+'




from f_carto import *
import cartopy.crs as ccrs

fig = plt.figure(fignr)
fignr+=1
plt.clf()


ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (1,1,1))
scale_bar(ax, 500, location= (0.06, 0.04))




"""get the ERA-5 data"""
import xarray as xr

year, month, day= 2019, 5, 3 
#year, month, day= 2000, 3, 9
#year, month, day= 1999, 12, 20
#year, month, day= 2002, 1, 26

hour= 12


filetype='surface'
#var= '2t'
#var= '10u'
#var= '10v'
var= 'U'





filetime= str(year)+'_'+str(month).zfill(2)+'_'+str(day).zfill(2)
ds= xr.open_dataset(Mediadir + "ERA5_STARS/"+filetype+ "_era5_"+ filetime + '.nc')
ds= ds.isel(time= hour)



if var== 'U': ds['U']= np.sqrt(ds['10u']**2 + ds['10v']**2)


print(ds[var].max())


now_datetime= np.datetime64(str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'T'+str(hour).zfill(2) )


tdelta= np.abs((S['datetime'] - now_datetime) / np.timedelta64(1, 'h'))
S_now= S[tdelta== tdelta.min()]
#S_now= S[np.abs((S['datetime'] - now_datetime) / np.timedelta64(1, 'h')) < 1]


data= ds[var].values
dist= 300

ds_lon= np.tile(ds.lon.values, (len(ds.lat), 1))
ds_lat= np.tile(ds.lat.values, (len(ds.lon), 1)).T


#ds['dist']= (dims=('lon', 'lat'), data= 110* np.sqrt( (S_now.Latitude.values[0]- ds_lat)**2+ (np.cos(np.deg2rad(S_now.Latitude.values[0]))* (S_now.Longitude.values[0]- ds_lon))**2)  )

dist_S_now= 110* np.sqrt( (S_now.Latitude.values[0]- ds_lat)**2+ (np.cos(np.deg2rad(S_now.Latitude.values[0]))* (S_now.Longitude.values[0]- ds_lon))**2)

localmax= np.max(ds[var].values[dist_S_now < dist])
localmin= np.min(ds[var].values[dist_S_now < dist])

print(localmax, localmin)

#import scipy.ndimage.filters as filters
#newmax= filters.maximum_filter(data, dist)

#def LocalMax(ds[var].values, )


"""scatter plot of systems at this time"""

##ds[var].plot.contourf(ax= ax, transform= ccrs.PlateCarree())
ds[var].plot(ax= ax, transform= ccrs.PlateCarree())
##
##ax.scatter(S_now.Longitude, S_now.Latitude, transform=ccrs.PlateCarree(), s= 3, color= 'k')
#
#
#fig = plt.figure(fignr)
#fignr+=1
#plt.clf()
#
#
#ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (1,1,1))
##newmax.plot(ax= ax, transform= ccrs.PlateCarree())
#
#cf= ax.pcolormesh(ds.lon, ds.lat, newmax, transform= ccrs.PlateCarree())
#cb= fig.colorbar(cf, ax= ax)
#
#cb.set_label(var, size=14)    
#cb.ax.tick_params(labelsize=14)


"""plot radius"""
lon, lat= S_now.Longitude.values[0], S_now.Latitude.values[0]
radius= S_now.Diameter.values[0]/2
#
plot_circle(ax, lon, lat, radius, edgecolor= 'k')














