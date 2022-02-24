#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 15:46:14 2020

@author: pst019
"""

import time
start = time.time()


import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    homedir= '/home/'+user+'/home/'
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/1692A00D929FEF8B/'

homedir=Mediadir+'home/'


import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import xarray as xr
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import numpy as np
from f_useful import *
from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs
from windrose import WindroseAxes

#
save= True
#save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM_tracks/'

fignr= 2


from scipy.stats import gaussian_kde


x= 3
y= 3

lifelim= 6 #lifetime limit


"""import Stoll list"""
Stoll_imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


Stoll= pd.read_csv(Stoll_imp_dir + 'Stoll_list_noRojo.csv')
Stoll= Stoll.set_index(['ID', 'time'])

S_nodes= pd.read_csv(Stoll_imp_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')
S_nodes= S_nodes.set_index(['ID', 'time'])

S_ERA5_2= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5.csv')
S_ERA5_2= S_ERA5_2.set_index(['ID', 'time'])

Stoll= pd.concat([Stoll, S_nodes, S_ERA5_2], axis=1)


Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps

Stoll= Stoll.reset_index()


fig = plt.figure(fignr, figsize= (13,10) )
fignr+=1
plt.clf()

n_all_tsteps= len(Stoll)


"""simple density map"""
for iSOM in range(1, x*y+1):
#    ax2= Plot_PlateCarree(fig, central_longitude= 0, extent= [-40, 85, 50, 85], subplot= (3,3,iSOM+1))
    ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (3,3,iSOM))

    Stoll_node= Stoll[Stoll.node == iSOM]
    
    
    perc_of_SOM= len(Stoll_node)/len(Stoll[Stoll.node > 0])
    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)', size = 15 )# , '+ str(int(np.round(after_tend_excl/before_tend_excl *100))) + '%)')


    Stoll_lon, Stoll_lat = Stoll_node.lon.values,Stoll_node.lat.values

    xy = np.vstack([Stoll_lon, Stoll_lat])
    z = gaussian_kde(xy)(xy)
    
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    Stoll_lon, Stoll_lat, z = Stoll_lon[idx], Stoll_lat[idx], z[idx]
    
#    fig, ax = plt.subplots()
    ax.scatter(Stoll_lon, Stoll_lat, c=z, s=50, edgecolor='', transform=ccrs.PlateCarree())


    #add all points of the SOM
#    ax.scatter(Stoll_node.lon, Stoll_node.lat, transform=ccrs.PlateCarree(), s= 8, color= 'y')


#
#
fig.tight_layout()
if save:
    save_file= savedir+ 'SOM_density_map'
    print(save_file)
    plt.savefig(save_file, bbox_inches='tight')



"""more advanced density map"""
dist= 250
lon_list= np.arange(-40, 85.1, 2) #0.5)
lat_list= np.arange(55, 85.1, 1)#0.25)

fig = plt.figure(fignr, figsize= (13,10) )
fignr+=1
plt.clf()



for iSOM in range(1, x*y+1):
    ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (3,3,iSOM))

    Stoll_node= Stoll[Stoll.node == iSOM]


    density= np.zeros((len(lon_list), len(lat_list)))
    
    
    Stoll_lon, Stoll_lat = Stoll_node.lon.values,Stoll_node.lat.values
    
    for ilat, lat in enumerate(lat_list):
        for ilon, lon in enumerate(lon_list):
            Stoll_dist= distance((lat, lon), (Stoll_lat, Stoll_lon))
            
            density[ilon, ilat]= len(Stoll_dist[Stoll_dist <= dist])
        
    
    ax.contourf(lon_list, lat_list, density.T, transform=ccrs.PlateCarree(), cmap='Blues')




fig.tight_layout()
