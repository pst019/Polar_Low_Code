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
    Mediadir= '/media/'+user+'/1692A00D929FEF8B/'

else:
    Mediadir= '/run/media/pst019/1692A00D929FEF8B/'
'/home/'+user+'/home/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import xarray as xr
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from f_useful import *
from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs


save= True
#save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM/'

fignr= 13

imp_list=True
#imp_list=False

imp_rot_nc=True
#imp_rot_nc=False


"""SOM var"""
Obs= "Obsnr1"
#Obs='mature'

plevel=850
ano=True
#ano= False
var= 't'
#var='z'
#var='vo'
#var='U'
#var='q'

if ano==True: var_full= var+ '_ano'
else: var_full= var

"""evolve var"""
#Obs_evol=1
Obs_evol= 'mature'
mat_evol_step= 100 #in percent how much of time the PL has overcome towards the mature phase (0- initial stage, 100- mature stage)
mat_evol_step= 75

ano_evol=True
#ano_evol=False
var_evol='t'

if ano_evol==True: var_evol_full= var_evol+ '_ano'
else: var_evol_full= var_evol

if Obs_evol== 'mature': Obs_evol_str= Obs_evol +str(mat_evol_step)
else: Obs_evol_str= 'Obsnr'+str(Obs_evol)

x=3 #3
y=4 #4




sym= False
#cmap='jet'
#cmap= 'viridis'
#cmap= 'Greys'
#cmap= 'Greys_r'
#cmap= 'Blues'
#cmap= 'Reds'

if "ano" in var_full: sym= True
cmap= 'RdBu_r'


filedir= Mediadir + "ERA5_STARS/PL_centred_fields/SOM/"+var_full+"_"+str(plevel)+'_'+Obs+"_x"+str(x)+"_y"+str(y)
ds_init= xr.open_dataset(filedir+"_cluster.nc")


"""import Stoll list"""
if imp_list:
    """Stoll systems"""
    test= 'version4'
    dist= 150
    system_char='initial' # for the Stoll_ind
    
    Stoll_STARS_file="ERA5_STARS/STARS-PMC-merge/"+test+"_dist_"+str(dist)+"_handfix.csv"
    Stoll = pd.read_csv(Mediadir+ Stoll_STARS_file)
    Stoll['time']= pd.DatetimeIndex(Stoll.time)
    
    Stoll= Stoll_Obs_nr(Stoll) #creates Stoll['Obs']
    Stoll= Stoll_excl(Stoll)
    Stoll= Stoll.rename(columns={'Stoll nr': 'ID'})
    
    Stoll = merge_ERA5(Stoll, Mediadir, "PL/STARS/Stoll_ERA5_handfix.csv")
    Stoll= Stoll.rename(columns={
        'U_max_200': 'Wind Speed 10m (max)',
        'vo': 'Vorticity$_{850}$ (centre)',
        'slp': 'Sea-level pressure (min)',
        'blh_max_200': 'Boundary layer height (max)',
        'cape_max_200': 'CAPE (max)', 
        'skt_med_200': 'Skin temperature (med)',
       'skt-t500_max_200': 'SST- T$_{500}$ (max)', 
       'skt-t700_max_200': 'SST -T$_{700}$ (max)',
        'tp_mean_200': 'Total precip. (mean) ',
        'cp_mean_200': 'Convective precip. (mean)',
        'sf_mean_200': 'Snow fall (mean)',
        'lsp_mean_200': 'Large-scale precip. (mean)', 
        'sshf_mean_200': 'Sensible heat flux (mean)',
        'slhf_mean_200': 'Latent heat flux (mean)',
       'grad_t850_max_200': 'Grad T$_{850}$ (max)',
       'baroclinic_gr_filter4_925-700_max_200': 'Baroclinic growth (max)' ,
       'barotropic_gr_filter4_850_max_200': 'Barotropic growth (max)',
       'vert_shear_angle925-700_mean_200': 'Vertical shear angle (mean)'
       })
    
    
    Stoll= Stoll.drop(columns=['Comment', 'track file', 'Rojo nr', 'Rojo nr old', 'Press_min', 'U_10min_knots', 'row_idx', 'track_idx'])
       
    
    """Stoll_individual_systems"""
    Stoll_ind= Stoll_individual_systems(Stoll, ID_name='ID', system_char=system_char)
    Stoll_ind= calc_system_duration(Stoll, Stoll_ind, ID_name= 'ID', Obs_name='Obs')
    
    Stoll_ind= Stoll_ind.drop(columns=[ 'STARS lat', 'STARS lon', 'STARS Obs nr'])
    Stoll= Stoll.drop(columns=[ 'STARS lat', 'STARS lon', 'STARS Obs nr', 'dist'])
    Stoll_ind= calc_Obs_mature(Stoll, Stoll_ind, intensity_var='Vorticity$_{850}$ (centre)')

    df_nodes= pd.read_csv(filedir+"_node_nr.txt", sep=" ")
    df_dist= pd.read_csv(filedir+"_distances.txt", sep=" ")
    
    df_nodes= pd.concat([df_nodes, df_dist], axis= 1)
    df_nodes['ID']= [str(int(df_nodes.PLnr[n]))[:-3]+'_'+str(int(df_nodes.PLnr[n])%1000) for n in range(len(df_nodes.PLnr))]
    df_nodes= df_nodes.set_index('ID')
    df_nodes= df_nodes.drop(columns='as.character(dates)')
    
    Stoll_ind=Stoll_ind.join(df_nodes['node'], how='outer')
    Stoll_ind=Stoll_ind.join(df_nodes['K_SOM.distances'], how='outer')





"""make the figure on which the original SOMs are displayed"""
fig= plt.figure(fignr, figsize= (2.7*x,2.7*y) )
fignr+=1
plt.clf()

if sym: vextr= np.max([np.max(ds_init.field), -np.min(ds_init.field)])
else: vmax, vmin= float(np.max(ds_init.field)), float(np.min(ds_init.field))
    
for iSOM in range(1, x*y +1):
    ax1= plt.subplot(y, x, iSOM)
    if sym:
        cf= plt.contourf(ds_init.x, ds_init.y, ds_init.field.sel(SOM=iSOM), cmap= cmap, vmin= -vextr, vmax= vextr)
        cs= plt.contour(ds_init.x, ds_init.y, ds_init.field.sel(SOM=iSOM), colors='k', vmin= -vextr, vmax= vextr, linewidth= 1)
    else:
        cf= plt.contourf(ds_init.x, ds_init.y, ds_init.field.sel(SOM=iSOM), cmap= cmap, vmin= vmin, vmax= vmax)
        cs= plt.contour(ds_init.x, ds_init.y, ds_init.field.sel(SOM=iSOM), colors='k', vmin= vmin, vmax= vmax, linewidth= 1)

    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.1f')

    perc_of_SOM= len(np.where(df_nodes.node == iSOM)[0])/len(df_nodes.node)
    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')

#fig.subplots_adjust(right=0.87)
#cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
#cb= fig.colorbar(cf, cax=cbar_ax)

plt.tight_layout()

fig.subplots_adjust(bottom=0.1)
cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.015])
cb= fig.colorbar(cf, cax=cbar_ax, orientation="horizontal")
#cb.set_clim(vmin, vmax)

if var_full== 'z_ano': labelvar= 'Geopotential height anomaly [m]'
elif var_full== 'z': labelvar= 'Geopotential height [m]'

elif var_full== 't_ano': labelvar= 'Temperature anomaly [K]'
elif var_full== 't': labelvar= 'Temperature [K]'

elif var_full== 'U_ano': labelvar= 'Wind speed anomaly [m/s]'
elif var_full== 'U': labelvar= 'Wind speed [m/s]'

elif var_full== 'q_ano': labelvar= 'Specific humidity anomaly [g(kg)]'
elif var_full== 'q': labelvar= 'Specific humidity [g/kg]'
else: labelvar= var_full

cb.set_label(labelvar, size=14) 

if save:
    save_name=savedir+ 'SOM_'+var_full+"_"+str(plevel)+"_"+Obs+"_x"+str(x)+"_y"+str(y)
    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')



"""evolution"""

fig= plt.figure(fignr, figsize= (2.7*x,2.7*y) )
fignr+=1
plt.clf()




if imp_rot_nc:
    print('import from all time step list')
    ds= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields/" +var_evol + '_'+ str(plevel) + '_allObs.nc')
    ds['ID']= ('time', [str(int(ds.PLnr[n]))[:-3]+'_'+str(int(ds.PLnr[n])%1000) for n in range(len(ds.PLnr))] )

    if type(Obs_evol)== int: ds_now= ds.where(ds.Obs== Obs_evol, drop=True)
    else:
        mature_index= [np.where(np.logical_and(ds.ID== Stoll_ind.ID[n], ds.Obs== Stoll_ind.Obs_mature[n]))[0][0] for n in range(len(Stoll_ind.ID)) ]
        ds_now= ds.isel(time= mature_index)

    ds_now= ds_now.dropna(dim='time') #to exclude nan values
    print('excluded PLs due to nan values:', set(Stoll_ind.ID.values) - set(ds_now.ID.values) )
    var_interp= ds_now[var_evol].values
    PLID= ds_now['ID'].values


if type(Obs_evol)== int: ds_now= ds.where(ds.Obs== Obs_evol, drop=True)
else:
#    mature_index= [np.where(np.logical_and(ds.ID== Stoll_ind.ID[n], ds.Obs== Stoll_ind.Obs_mature[n]))[0][0] for n in range(len(Stoll_ind.ID)) ]
    mature_index= [np.where(np.logical_and(ds.ID== Stoll_ind.ID[n], ds.Obs== int(np.round(1+ mat_evol_step/100 *(Stoll_ind.Obs_mature[n] -1))) ))[0][0] for n in range(len(Stoll_ind.ID)) ]

    ds_now= ds.isel(time= mature_index)

ds_Obs_evol= ds_now.dropna(dim='time')

if ano_evol == True:
    ds_Obs_evol[var_evol_full]= (('time', 'x', 'y'), ds_Obs_evol[var_evol].values- np.mean(ds_Obs_evol[var_evol].values, axis= (1,2))[:, np.newaxis, np.newaxis])
    

#if sym: vextr= np.max([np.max(ds_Obs_evol[var_evol_full]), -np.min(ds_Obs_evol[var_evol_full])])
#else:
#    vmax, vmin= float(np.max(ds_Obs_evol[var_evol_full])), float(np.min(ds_Obs_evol[var_evol_full]))
#

for iSOM in range(1, x*y +1):
    ax1= plt.subplot(y, x, iSOM)

    node_list= Stoll_ind[Stoll_ind.node== iSOM].index.values
    weight_list= 1/Stoll_ind[Stoll_ind.node== iSOM]['K_SOM.distances'].values
    ds_now= ds_Obs_evol.where(ds_Obs_evol.ID.isin(node_list), drop=True)
    

    
    if sym:
#        cf= plt.contourf(ds_now.x, ds_now.y, np.average(ds_now[var_evol_full], axis= 0, weights=weight_list), cmap= cmap, vmin= -vextr, vmax= vextr)
        cf= plt.contourf(ds_now.x, ds_now.y, np.average(ds_now[var_evol_full], axis= 0), cmap= cmap, vmin= -vextr, vmax= vextr)
        
        cs= plt.contour(ds_now.x, ds_now.y, np.mean(ds_now[var_evol_full], axis= 0), colors='k', vmin= -vextr, vmax= vextr)
    else:
        cf= plt.contourf(ds_now.x, ds_now.y, np.mean(ds_now[var_evol_full], axis= 0), cmap= cmap, vmin= vmin, vmax= vmax)
        cs= plt.contour(ds_now.x, ds_now.y, np.mean(ds_now[var_evol_full], axis= 0), colors='k', vmin= vmin, vmax= vmax)

    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.1f')

    perc_of_SOM= len(np.where(df_nodes.node == iSOM)[0])/len(df_nodes.node)
    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')



plt.tight_layout()

fig.subplots_adjust(bottom=0.1)
cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.015])
cb= fig.colorbar(cf, cax=cbar_ax, orientation="horizontal")



if var_evol_full== 'z_ano': labelvar= 'Geopotential height anomaly [m]'
elif var_evol_full== 'z': labelvar= 'Geopotential height [m]'

elif var_evol_full== 't_ano': labelvar= 'Temperature anomaly [K]'
elif var_evol_full== 't': labelvar= 'Temperature [K]'

elif var_evol_full== 'U_ano': labelvar= 'Wind speed anomaly [m/s]'
elif var_evol_full== 'U': labelvar= 'Wind speed [m/s]'

elif var_evol_full== 'q_ano': labelvar= 'Specific humidity anomaly [g(kg)]'
elif var_evol_full== 'q': labelvar= 'Specific humidity [g/kg]'
else: labelvar= var_evol_full

cb.set_label(labelvar, size=14) 

if save:
    save_name=savedir+ 'SOM_'+var_full+"_"+str(plevel)+"_"+Obs+'--mean_'+var_evol_full+"_"+str(plevel)+"_"+Obs_evol_str+"_x"+str(x)+"_y"+str(y)
    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')




"""standard deviation"""
fig= plt.figure(fignr, figsize= (2.7*x,2.7*y) )
fignr+=1
plt.clf()


cmap= 'viridis'


for iSOM in range(1, x*y +1):
    ax1= plt.subplot(y, x, iSOM)

    node_list= Stoll_ind[Stoll_ind.node== iSOM].index.values
    ds_now= ds_Obs_evol.where(ds_Obs_evol.ID.isin(node_list), drop=True)
    

    

    cf= plt.contourf(ds_now.x, ds_now.y, np.std(ds_now[var_evol_full], axis= 0), cmap= cmap, vmin= 0)#, vmax= vmax)
    cs= plt.contour(ds_now.x, ds_now.y, np.std(ds_now[var_evol_full], axis= 0), colors='k')#, vmin= vmin, vmax= vmax)

    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.1f')

    perc_of_SOM= len(np.where(df_nodes.node == iSOM)[0])/len(df_nodes.node)
    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')



plt.tight_layout()

fig.subplots_adjust(bottom=0.1)
cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.015])
cb= fig.colorbar(cf, cax=cbar_ax, orientation="horizontal")
cb.set_label('Standard deviation', size=14) 

if save:
    save_name=savedir+ 'SOM_'+var_full+"_"+str(plevel)+"_"+Obs+'--std_'+var_evol_full+"_"+str(plevel)+"_"+Obs_evol_str+"_x"+str(x)+"_y"+str(y)
    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')