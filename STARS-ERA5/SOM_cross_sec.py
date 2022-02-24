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

if user=='pst019':
    homedir= '/home/'+user+'/home/'
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/1692A00D929FEF8B/'
'/home/'+user+'/home/'

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

plt.rcParams.update({'font.size': 14})


#
save= True
save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM4/'
#
fignr= 7


load_ds=True
load_ds=False


"""SOM var"""
Obs= 'allObs_SOMprep' #the other versions were deleted

x=3 #3
y=3 #x+1

mirror_top_bottom=False
rotate_clock=False

if x== 3 and y == 4: mirror_top_bottom=True
if x== 3 and y == 3: rotate_clock=True

PLCG_type, smooth_param= 'track_smth', '1E-3'

lifelim= 6 #lifetime limit

vo_tend_excl= True
vo_tend_thresh= 0



var='z'
var2='w'



cross_loc= 0
cross= 'x'
#cross= 'y' 
    
    
"""import Stoll list"""
if load_ds:

    Stoll_imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'
    
    
    Stoll= pd.read_csv(Stoll_imp_dir + 'Stoll_list_noRojo.csv')
    Stoll= Stoll.set_index(['ID', 'time'])
    
    S_nodes= pd.read_csv(Stoll_imp_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')
    S_nodes= S_nodes.set_index(['ID', 'time'])
    
    Stoll= pd.concat([Stoll, S_nodes], axis=1)
    
    
    Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps
    
    Stoll= Stoll.reset_index()
    Stoll['PLnr']= [int(ID.split('_')[0])*1E3 + int(ID.split('_')[1]) for ID in Stoll.ID.values]


    print('import from all time step list')

    file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +var + '_all_levs_allObs_track-smth-'+smooth_param+'.nc'
    ds= xr.open_dataset(file)



    file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +var2 + '_all_levs_allObs_track-smth-'+smooth_param+'.nc'
    ds2= xr.open_dataset(file)
    
    ds= xr.merge([ds, ds2]) 


    
    ds2= imp_ds('adv_t', 500, 'adv_t_500', lsm_mask=False, compute=True, imp_var=['t', 'u', 'v'], data_gausfilter=4, smooth_param= '1E-3')
    ds= xr.merge([ds, ds2['adv_t_500']])

    ds2= imp_ds('adv_t', 925, 'adv_t_925', lsm_mask=False, compute=True, imp_var=['t', 'u', 'v'], data_gausfilter=4, smooth_param= '1E-3')
    ds= xr.merge([ds, ds2['adv_t_925']])


    
    """add the node to ds """
    ds_node= np.zeros(len(ds.time))
    ds_vo_tendency= np.zeros(len(ds.time))

    for s in range(len(Stoll)):
        S= Stoll.iloc[s]
        ds_index= np.argwhere(np.logical_and(ds.PLnr.values== S.PLnr, ds.time.values== np.datetime64(S.time)))[0][0]
        ds_node[ds_index]= S.node
        ds_vo_tendency[ds_index]= S.vo_tendency

    
    ds['node']= (('time'), ds_node)
    ds['vo_tendency']= (('time'), ds_vo_tendency)






iSOM = 1






for iSOM in range(1, x*y +1):
    
    print('SOM ', iSOM)
    

    ds_node= ds.where(ds.node == iSOM, drop=True)
    
    
    
    
    #    
    if vo_tend_excl:
        ds_node= ds_node.where(ds_node.vo_tendency > vo_tend_thresh, drop=True)
    
    
    
    if cross== 'y':
        ds_node = ds_node.sel(x= cross_loc)
        dsx= ds_node.y
    #ds.where(np.abs(ds.y) < 400, drop= True)
    if cross== 'x':
        ds_node = ds_node.sel(y= cross_loc)
        dsx= ds_node.x
    
    
    
    
    
    
    
    fig= plt.figure(fignr, figsize= (7,5) )
    plt.clf()
    
    
    """1000hPa"""
    ax1= plt.subplot(3, 1, 3)
    
    #plt.plot(ds_node.x, ds_node.z.sel(plev=1000).mean(dim='time'))
    ds_var = ds_node.z.sel(plev=925).mean(dim='time')
    ds_var -= ds_var.mean()
    
    
    color = 'k' #'tab:blue'
    ax1.plot(dsx, ds_var, color=color)
    ax1.set_ylabel('z$^{\prime}$ [m]', color=color)
    ax1.tick_params(axis='y', labelcolor=color)
#    ax1.set_yticks(np.arange(np.round(np.min(ds_var/2),-1)*2, np.round(np.max(ds_var/2),-1)*2+.1, 20))
    ax1.locator_params(axis='y', nbins=5)

    
    
    
    ax1_1 = ax1.twinx()
    ds_var = ds_node['adv_t_925'].mean(dim='time')
    #ds_var -= ds_var.mean()
    
    color = 'tab:red'
    ax1_1.plot(dsx, ds_var, color=color)
    ax1_1.set_ylabel(r"$-u \cdot \nabla_h T$"+'\n'+r"[K$\cdot$h$^{-1}$]", color=color)
    ax1_1.tick_params(axis='y', labelcolor=color)
#    ax1_1.set_yticks(np.arange(np.round(np.min(ds_var),1), np.round(np.max(ds_var),1)+.05, 0.1))
    ax1_1.locator_params(axis='y', nbins=5)

    
    """700hPa"""
    ax2= plt.subplot(3, 1, 2)#, sharex=ax1)
    
    ds_var = ds_node.w.sel(plev=700).mean(dim='time')
    #ds_var -= ds_var.mean()
#    ax2.plot(dsx, -ds_var, 'x')

    ax2.quiver(dsx, np.zeros(len(dsx)), np.zeros(len(dsx)), -ds_var, angles='xy', scale_units='xy', scale=1) #angle, scale and unit make the arrow the correct length
    ax2.set_ylim(-0.2,1.3)
    ax2.set_yticks([0, 0.5, 1])
    ax2.set_ylabel(r"$\omega$ [-Pa$\cdot$s$^{-1}$]")
    
    
    """500hPa"""
    ax3= plt.subplot(3, 1, 1)#, sharex=ax1)
    #plt.plot(ds_node.x, ds_node.z.sel(plev=1000).mean(dim='time'))
    ds_var = ds_node.z.sel(plev=500).mean(dim='time')
    ds_var -= ds_var.mean()
    
    
    color = 'k'#'tab:blue'
    ax3.plot(dsx, ds_var, color=color)
    ax3.set_ylabel('z$^{\prime}$ [m]', color=color)
    ax3.tick_params(axis='y', labelcolor=color)
#    ax3.set_yticks(np.arange(np.round(np.min(ds_var/2),-1)*2, np.round(np.max(ds_var/2),-1)*2+.1, 20))
    ax3.locator_params(axis='y', nbins=5)

    
    
    
    ax3_1 = ax3.twinx()
    ds_var = ds_node['adv_t_500'].mean(dim='time')
    #ds_var -= ds_var.mean()
    
    color = 'tab:red'
    ax3_1.plot(dsx, ds_var, color=color)
    ax3_1.set_ylabel(r"$-u \cdot \nabla_h T$"+'\n'+r"[K$\cdot$h$^{-1}$]", color=color)
    ax3_1.tick_params(axis='y', labelcolor=color)
#    ax3_1.set_yticks(np.arange(np.round(np.min(ds_var),1), np.round(np.max(ds_var),1)+.05, 0.1))
    ax3_1.locator_params(axis='y', nbins=5)
    
    #ds_var = ds_node.w.sel(plev=500).mean(dim='time')
    #ds_var -= ds_var.mean()
    #plt.quiver(dsx, np.zeros(len(dsx)), np.zeros(len(dsx)), -ds_var, scale= 10)
    
    
    
    
    
    """ticks and stuff"""
    xlim= 500
    
    ax1.set_xlim([-xlim, xlim])
    ax2.set_xlim([-xlim, xlim])
    ax3.set_xlim([-xlim, xlim])
    
    
    if cross== 'y': ax1.set_xlabel('Distance in propagation direction [km]')
    if cross== 'x': ax1.set_xlabel('Distance across propagation direction [km]')
    
    ax2.set_xticklabels([])
    ax3.set_xticklabels([])
    
    
    #    
    #
#    """make the arrow"""
#    if cross== 'y':
#        arrow_ax = fig.add_axes([0.7, 0.03, 0.18, 0.04])
#        import matplotlib.patches as mpatches
#        
#        arrow = mpatches.FancyArrowPatch((0.5, 0.5), (1, 0.5), mutation_scale=20)
#        arrow_ax.add_patch(arrow)
#        arrow_ax.text(0. ,0.3 , 'Prop', fontweight= 'bold')#, fontsize=13)
#        arrow_ax.set_frame_on(False)
#        arrow_ax.get_xaxis().set_visible(False)
#        arrow_ax.get_yaxis().set_visible(False)
    plt.subplots_adjust(left=0.17, right=0.83)
    
    fig.text(0.01, 0.23, r"925 hPa", va='center', rotation='vertical')
    fig.text(0.01, 0.5, r"700 hPa", va='center', rotation='vertical')
    fig.text(0.01, 0.75, r"500 hPa", va='center', rotation='vertical')

    
#    fig.tight_layout()

    #
    if save:
        save_name=savedir+ 'vertical_SOM'+str(iSOM)+'_cr-'+cross
    
        print(save_name)
        fig.savefig(save_name , bbox_inches='tight')
    
    
    #
    #

