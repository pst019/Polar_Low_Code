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
    Mediadir= '/run/media/pst019/PatsOrange/'

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
import matplotlib as mpl
from windrose import WindroseAxes

#
save= True
#save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM_tracks/'

fignr= 4

typ='SOM'
typ='shear_cat'
node_label= ['Forward', 'Right', 'Reverse', 'Left', 'Weak']

x= 3
y= 3


vo_intensifying= True
vo_decaying= False
lat_lon_excl= False
#lifelim= 6 #lifetime limit


"""import Stoll list"""
Stoll_imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


Stoll= pd.read_csv(Stoll_imp_dir + 'Stoll_list_noRojo.csv')
Stoll= Stoll.set_index(['ID', 'time'])

Stoll_dir_speed= pd.read_csv(Stoll_imp_dir + 'Stoll_list_dir-speed_smth1E-3.csv')
Stoll_dir_speed= Stoll_dir_speed.set_index(['ID', 'time'])


S_nodes= pd.read_csv(Stoll_imp_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')
S_nodes= S_nodes.set_index(['ID', 'time'])

S_ERA5_2= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5.csv')
S_ERA5_2= S_ERA5_2.set_index(['ID', 'time'])

S_shear_cats= pd.read_csv(Stoll_imp_dir + 'shear_categories.csv')
S_shear_cats= S_shear_cats.set_index(['ID', 'time'])


Stoll= pd.concat([Stoll, S_nodes, S_ERA5_2, Stoll_dir_speed, S_shear_cats], axis=1)


#Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps

Stoll= Stoll.reset_index()



    

speed_var= 'prop_speed'
dir_var= 'prop_dir'
bins= np.arange(0, 16, 5)
title= "Propagation speed [km$\cdot$h$^{-1}$]"


#from windrose import *


"""more advanced density map"""
dist= 250
density_levels= np.arange(100, 501, 50)
#density_levels= np.arange(25, 501, 25)

lon_list= np.arange(-40, 85.1, 2) #0.5)
lat_list= np.arange(55, 85.1, 1)#0.25)

#fig = plt.figure(fignr, figsize= (13,10) )
#fignr+=1
#plt.clf()

n_all_tsteps= len(Stoll)





"""get the ERA-5 data"""
imp_data=True

plevel=None
vmin,vmax=None,None

"""specify field"""
#filetype='surface'
#var= '2t'
#var= '10u'
#var= '10v'
#var= 'U'
#var= 'skt'
#var='msl'

filetype='tracking'
plevel= 850
#var='U'
#var= 'u'
#var='Vort'
var= 't'
if plevel== 850: color_levels=  np.arange(250, 271, 2)
if plevel== 700: color_levels=  np.arange(244, 265, 2)

"""the contour variable"""
cont_var='z'


cmap= 'RdBu_r'







if typ == 'SOM': cat_list= np.arange(1, x*y +1)
elif typ == 'shear_cat': cat_list= np.arange(1, 6)

for node in [5]: #cat_list:
#    ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (3,3,node))

    fig = plt.figure(fignr, figsize= (6,5.5) )
#    fignr+=1
    plt.clf()
    
    ax= Plot_Polar_Stereo(fig, central_longitude= 20, extent= [-10, 50, 60, 80])

    
    if typ == 'SOM':
        Stoll_node= Stoll[Stoll.node == node]
        perc_of_SOM= len(Stoll_node)/len(Stoll[Stoll.node > 0])
        plt.title('SOM '+str(node)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)', size = 15 )# , '+ str(int(np.round(after_tend_excl/before_tend_excl *100))) + '%)')

    elif typ == 'shear_cat':
        Stoll_node= Stoll[Stoll['shear_category'] == node]
        

#        Stoll_node= Stoll_node[np.logical_and(Stoll_node['prop_dir'] < 90, Stoll_node['prop_dir'] >= 0)]
#        Stoll_node= Stoll_node[np.logical_or(Stoll_node['prop_dir'] < 0, Stoll_node['prop_dir'] >= 270)]

    if lat_lon_excl== True:
#        Stoll_node= Stoll_node[Stoll_node['lon'] > 30]
#        Stoll_node= Stoll_node[Stoll_node['lat'] < 70]
        Stoll_node= Stoll_node[np.logical_and(Stoll_node['lat'] > 70, Stoll_node['lon'] < 30)]

    if vo_intensifying== True:
        Stoll_node= Stoll_node[Stoll_node['vo_tendency'] > 0]
    elif vo_decaying== True:
        Stoll_node= Stoll_node[Stoll_node['vo_tendency'] < 0]

        perc_of_SOM= len(Stoll_node)/len(Stoll)
        plt.title(node_label[node-1]+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)', size = 15 )# , '+ str(int(np.round(after_tend_excl/before_tend_excl *100))) + '%)')



    """import data for the environment"""
    Stoll_node['date']=  pd.to_datetime(Stoll_node['time']).dt.strftime('%Y_%m_%d')
    
    date_list_different= remove_dublicate(Stoll_node['date'] )
     
    if imp_data:
        for i, date in enumerate(date_list_different):
            
            ds0= xr.open_dataset(Mediadir + "ERA5_STARS/data/" +filetype+'_era5_'+ date + '.nc')
            
            ds0= ds0.sel( time= np.array(pd.to_datetime(Stoll_node[Stoll_node['date'] == date].time) ) ) #select the correct hours
            ds0= ds0[[var, cont_var]]
        
            if plevel:
                ds0= ds0.sel(plev= plevel*100)
        
            if i == 0:
                ds= ds0
                
            else:
                ds= xr.concat([ds, ds0], 'time')
                
            if i%10 == 0: print(i)



    """plot the environmental composites"""
    dsmean= ds[var].mean(dim='time')    
    dscontmean = ds[cont_var].mean(dim='time')  /9.81  
    
    cf= ax.contourf(ds.lon, ds.lat, dsmean, transform= ccrs.PlateCarree(), levels= color_levels, cmap= cmap, extend='both')
    cs= ax.contour(ds.lon, ds.lat, dscontmean, transform= ccrs.PlateCarree(), colors='k', linewidths= 1, levels= np.arange( 0, 5E3, 20))
    plt.clabel(cs, fontsize=12, inline=1, fmt='%1.0f')



    """the density maps"""
    density= np.zeros((len(lon_list), len(lat_list)))
    
    
    Stoll_lon, Stoll_lat = Stoll_node.lon.values,Stoll_node.lat.values
    
    for ilat, lat in enumerate(lat_list):
        for ilon, lon in enumerate(lon_list):
            Stoll_dist= distance((lat, lon), (Stoll_lat, Stoll_lon))
            
            density[ilon, ilat]= len(Stoll_dist[Stoll_dist <= dist])
        
    
    cf= ax.contour(lon_list, lat_list, density.T, transform=ccrs.PlateCarree(), levels= density_levels, colors='g')
    plt.clabel(cf, density_levels[::2], fontsize=12, fmt='%1.0f')


#    rect=[0.57,0.21,0.3,0.3] 
    rect=[0.61,0.2,0.3,0.3] 
    
    wa=WindroseAxes(fig, rect) #fig, rect)
    fig.add_axes(wa)
    
    wr= wa.bar(Stoll_node[dir_var], Stoll_node[speed_var], normed=True, bins=bins, opening=1 , 
                         edgecolor='k', cmap=mpl.cm.viridis_r)
#    wa.set_xticklabels(['E', '', 'N', '', 'W', '', 'S', ''])
#    wa.set_xticklabels(['', 'NE', '', 'NW', '', 'SW', '', 'SE'])
    wa.set_xticklabels(['', '', 'N', '', 'W', '', '', ''])
    
    wa.set_rgrids([]) #remove the circular axis
    

#    fig.tight_layout()


#    fig.subplots_adjust(bottom=0.11)
#    cbar_ax = fig.add_axes([0.09, 0.06, 0.4, 0.015])
#    cb= fig.colorbar(cf, cax=cbar_ax, orientation="horizontal")

    wa.set_legend()
    wa.legend(title=title, loc= (-1.5, -.55), ncol= 2, decimal_places=0)

    
    if save:
        save_file= savedir+ 'environment+density_map_'+typ+'_'+str(node)
        if vo_intensifying: save_file += '_intensifying'
        elif vo_decaying: save_file += '_decaying'

        print(save_file)
        plt.savefig(save_file, bbox_inches='tight')
        


"""colorbar figure"""
#fignr+= 1
#figbar = plt.figure(fignr, figsize= (4,3))
#plt.clf()
#
#
#cb= figbar.colorbar(cf, orientation="horizontal")#, labelsize= 13)
#cb.set_label('Number of trackpoints within '+str(dist)+'km distance')#, size= 15)
#
#
#if save:
#    save_file= savedir+ 'environment+density_map_'+typ+'_legend'
#    print(save_file)
#    plt.savefig(save_file, bbox_inches='tight')

"""Windrose label"""
#fignr+= 1
#figrose = plt.figure(fignr)
#plt.clf()
#
#figrose.legend(wa.get_legend())
#
#wr.legend(title=title, loc= (-1.1, 0), ncol= 2, decimal_places=0)
#
#figrose.figlegend(*wa.get_legend_handles_labels() )