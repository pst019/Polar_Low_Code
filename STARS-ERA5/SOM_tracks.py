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

#
save= True
save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM_tracks/'

fignr= 1





x= 3
y= 3

lifelim= 6 #lifetime limit
add_PLnr = False
add_allPL_points = False

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
plt.clf()

n_all_tsteps= len(Stoll)

for iSOM in range(1, x*y+1):
#    ax2= Plot_PlateCarree(fig, central_longitude= 0, extent= [-40, 85, 50, 85], subplot= (3,3,iSOM+1))
    ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (3,3,iSOM))

    Stoll_node= Stoll[Stoll.node == iSOM]
    
    
    perc_of_SOM= len(Stoll_node)/len(Stoll[Stoll.node > 0])
    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)', size = 15 )# , '+ str(int(np.round(after_tend_excl/before_tend_excl *100))) + '%)')


    #add all points of the SOM
    if add_allPL_points: ax.scatter(Stoll_node.lon, Stoll_node.lat, transform=ccrs.PlateCarree(), s= 8, color= 'y')

#    Stoll_init= Stoll_node[Stoll_node['Obs']== 1]
#    ax.scatter(Stoll_init.lon, Stoll_init.lat, transform=ccrs.PlateCarree(), s= 12, color= 'b')


    for ID in remove_dublicate(Stoll_node['ID']):
        Stoll_ID= Stoll_node[Stoll_node['ID'] == ID]

        split_loc= np.where(Stoll_ID.Obs[1:].values - Stoll_ID.Obs[:-1].values != 1)[0]
#        if len(split_loc[0]) == 0: #no split of the track
#            ax.plot(Stoll_ID.lon, Stoll_ID.lat, transform=ccrs.PlateCarree(), color= 'k')
            
        if len(split_loc) > 0:
            split_loc += 1 #since the way of calculating the split_loc shifts everything
            for i_split_loc in range(len(split_loc)):            
                value_split_loc= split_loc[i_split_loc]
#        
                Stoll_ID_1= Stoll_ID[:value_split_loc]
                ax.plot(Stoll_ID_1.lon, Stoll_ID_1.lat, transform=ccrs.PlateCarree(), color= 'k', lw= 1)
                ax.scatter(Stoll_ID_1.lon.values[0], Stoll_ID_1.lat.values[0], transform=ccrs.PlateCarree(), s= 10, color= 'b')
                if add_PLnr: ax.text(Stoll_ID_1.lon.values[0], Stoll_ID_1.lat.values[0], Stoll_ID_1.ID.values[0], transform=ccrs.PlateCarree(), color= 'g')

                Stoll_ID= Stoll_ID[value_split_loc:]
                
                split_loc -= value_split_loc
        
        ax.plot(Stoll_ID.lon, Stoll_ID.lat, transform=ccrs.PlateCarree(), color= 'k', lw= 1)

        ax.scatter(Stoll_ID.lon.values[0], Stoll_ID.lat.values[0], transform=ccrs.PlateCarree(), s= 10, color= 'b')
        if add_PLnr: ax.text(Stoll_ID.lon.values[0], Stoll_ID.lat.values[0], Stoll_ID.ID.values[0], transform=ccrs.PlateCarree(), color= 'g')



    """all initial time steps"""
    Stoll_init= Stoll_node[Stoll_node['Obs']== 1]
    ax.scatter(Stoll_init.lon, Stoll_init.lat, transform=ccrs.PlateCarree(), s= 10, color= 'r')

#        
#        ax.plot(Stoll_ID.lon, Stoll_ID.lat, transform=ccrs.PlateCarree(), color= 'k')
#        ax.scatter(Stoll_ID['Obs'.lon, Stoll_node.lat, transform=ccrs.PlateCarree(), s= 1, color= 'y')


#

#    S_tracks= prepare_tracks(Stoll_node, system_id_label="ID")
#    for track in S_tracks:
#        track.plot_track(ax=ax, color='y')#, label='Stoll excluded')   
#
#
fig.tight_layout()
if save:
    save_file= savedir+ 'SOM_tracks'
    print(save_file)
    plt.savefig(save_file, bbox_inches='tight')
