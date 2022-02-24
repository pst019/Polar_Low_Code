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
    Mediadir= '/media/'+user+'/PatsOrange/'

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

#
save= False
#save= True
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM2/'
#
fignr= 14

imp_list=True
#imp_list=False

imp_rot_nc=True
#imp_rot_nc=False


"""SOM var"""
#Obs= "Obsnr1"
#Obs= "Obsnr1_prep"
#Obs='mature'
Obs= 'allObs_SOMprep'

x=3 #3
y=3 #x+1

mirror_top_bottom=False
rotate_clock=False

if x== 3 and y == 4: mirror_top_bottom=True
if x== 3 and y == 3: rotate_clock=True

#PLCG_type, proplim= 'stearing_flow', 3 #propagation speed limit
PLCG_type, smooth_param= 'track_smth', '1E-3'

lifelim= 6 #lifetime limit

Splevel=850
Sano=True
#ano= False
Svar= 't'
#Svar='z'
#Svar='vo'
#Svar='U'
#Svar='q'

if Sano==True: Svar_full= Svar+ '_ano'
else: Svar_full= Svar


S_sym= False
if "ano" in Svar_full: S_sym= True
S_cmap= 'RdBu_r'





"""evolve var"""
#Obs_evol=1
#Obs_evol= 'mature'
#mat_evol_step= 100 #in percent how much of time the PL has overcome towards the mature phase (0- initial stage, 100- mature stage)
#mat_evol_step= 25
Obs_evol='allObs'

if Obs_evol== 'mature': Obs_evol_str= Obs_evol +str(mat_evol_step)
elif Obs_evol =='allObs': Obs_evol_str= Obs_evol
else: Obs_evol_str= 'Obsnr'+str(Obs_evol)

#ano_var=True
ano_var=False
#var='t'
#var='z'
#var='vo'
#var='d'
#var='w'
#var='pv'
#plevel_var=850


lsm_mask=False

"""variables without plevel"""
plevel_var=''

#var='lcc'
#var='mcc'
#var='hcc'
#cmap= 'Greys_r'

var='cp'
#var='lsp'
cmap= 'Blues'

#var='slhf'
#var='sshf'

#var='pres_pv'
#var='z_pv'
#var='cin'
#var='cape'
#var='skt'
#cmap= 'Reds'

#var, land='ci', 1 #land= [0,1] specifies if land is set to 0 or to 1
#cmap='Blues'

#variables that needs to be computed
var= 'N'
#plevel= [925, 500]

var_full= var
if ano_var==True:  var_full= var_full+ '_ano'
if type(plevel_var) == int: var_full= var_full+ '_'+str(plevel_var)




"""specify the levels"""
levels=[]


if var_full== 't_1000': levels= np.arange(258, 274.1, 2)
elif var_full== 't_500': levels= np.arange(230, 240.1, 2)
elif var_full== 'z_500': levels= np.arange(5000, 5160.1, 20)
elif var_full== 'z_1000': levels= np.arange(-60, 100.1, 20)
elif var_full== 'z_pv': levels= np.arange(6200, 8400.1, 200)
elif var_full== 'pres_pv': levels= np.arange(300, 421, 20)
elif var_full== 'cin': levels= np.arange(0, 2.1, .25)
elif var_full== 'cape': levels= np.arange(0, 45.1, 5)
elif var_full== 'ci': levels= np.arange(0, 1.01, .1)
elif var_full== 'slhf': levels= np.arange(0, 200.1, 10)
elif var_full== 'sshf': levels= np.arange(0, 200.1, 10)
elif var_full== 'pv_500':  levels = np.arange(0, 20.1, 1)*10E-8
elif var_full=='skt': levels= np.arange(-8, 8.1, 1)

elif 't_ano' in var_full:  levels = np.linspace(-7, 7, 15)
elif 'z_ano' in var_full:  levels = np.arange(-200, 201, 20)
elif 'vo' in var_full:  levels = np.arange(-50, 50.1, 5)
elif 'pv' in var_full:  levels = np.arange(0, 7.1, 1)*10E-8

elif 'cp' in var_full or 'lsp' in var_full: levels= np.arange(0, .51, .1) # levels= np.linspace(0, 1, 11)
elif 'tp' in var_full : levels= np.linspace(0, 1, 11)

elif 'cc' in var_full: levels= np.arange(0.2, 1.01, .1)
elif 'd' in var_full:  levels = [-4, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 4]
elif 'w' in var_full: levels= [-1.2, -1, -.8, -.4, -.2, .2, .4, .6, .8, 1, 1.2]



scale= .7 #scale factor for the color bar. it takes scale* [min, max] of the all time step values, only relevant if levels = []
sym= False # only relevant if levels=[]
#cmap='jet'
#cmap= 'viridis'
#cmap= 'Greys_r'
#cmap= 'RdBu_r'

if ano_var: sym= True
if sym: cmap= 'RdBu_r'



"""import Stoll list"""
if imp_list:
    if PLCG_type== 'stearing_flow':
        SOM_filedir= Mediadir + "ERA5_STARS/PL_centred_fields/SOM/"+Svar_full+"_"+str(Splevel)+'_'+Obs+'_vel'+str(proplim)+'_dur'+str(lifelim)+"_x"+ str(x)+"_y"+str(y)
    elif PLCG_type== 'track_smth':
        SOM_filedir= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/SOM/"+Svar_full+"_"+str(Splevel)+'_'+Obs+ '_track-smth-'+smooth_param+'_dur'+str(lifelim)+"_x"+ str(x)+"_y"+str(y)

    ds_init= xr.open_dataset(SOM_filedir+"_cluster.nc")

    Stoll, Stoll_ind= imp_standard_Stoll() #SOM_filedir)#, Obs_evol)
    df_nodes= pd.read_csv(SOM_filedir+"_node_nr.txt", sep=" ")

    
    if Obs_evol!='allObs':
        df_nodes['ID']= [str(int(df_nodes.PLnr[n]))[:-3]+'_'+str(int(df_nodes.PLnr[n])%1000) for n in range(len(df_nodes.PLnr))]
        df_nodes= df_nodes.set_index('ID')
#        Stoll_ind=Stoll_ind.join(df_nodes['node'], how='outer')

  
    if Obs_evol=='allObs':
        df_nodes['time']= pd.to_datetime(df_nodes.date.apply(str), format='%Y%m%d%H')
        df_nodes['ID']= [str(int(df_nodes.PLnr[n]))[:-3]+'_'+str(int(df_nodes.PLnr[n])%1000) for n in range(len(df_nodes.PLnr))]
#        df_nodes= df_nodes.drop(columns=['date', 'PLnr'])
        Stoll= Stoll.set_index(['ID', 'time'])
        df_nodes= df_nodes.set_index(['ID', 'time'])
        
        """to rotate/mirror the nodes"""
        df_nodes_orig= df_nodes.copy()
        if mirror_top_bottom:
            df_nodes['node']= (y-1 -(df_nodes.node.values-1)//x)*x + (df_nodes.node.values -1)%x +1
        if rotate_clock:
            x_n= (df_nodes['node'].values-1)%x +1
            y_n= (df_nodes['node'].values-1)//x +1
            y_n= y+1 - y_n #reverse the y_nodes
            
            df_nodes['node']= (x_n -1)*x + y_n
            
#        Stoll= pd.concat([Stoll, df_nodes], axis= 1)



"""make the figure on which the original SOMs are displayed"""
fig= plt.figure(fignr, figsize= (2.5*x,2.5*y +0) )
fignr+=1
plt.clf()

if S_sym: vextr= np.max([np.max(ds_init.field), -np.min(ds_init.field)])
else: vmax, vmin= float(np.max(ds_init.field)), float(np.min(ds_init.field))
    
for iSOM in range(1, x*y +1):
    ax1= plt.subplot(y, x, iSOM)
    if S_sym:
        cf= plt.contourf(ds_init.x, ds_init.y, ds_init.field.sel(SOM=iSOM), cmap= S_cmap, vmin= -vextr, vmax= vextr)
        cs= plt.contour(ds_init.x, ds_init.y, ds_init.field.sel(SOM=iSOM), colors='k', vmin= -vextr, vmax= vextr, linewidth= 1)
    else:
        cf= plt.contourf(ds_init.x, ds_init.y, ds_init.field.sel(SOM=iSOM), cmap= S_cmap, vmin= vmin, vmax= vmax)
        cs= plt.contour(ds_init.x, ds_init.y, ds_init.field.sel(SOM=iSOM), colors='k', vmin= vmin, vmax= vmax, linewidth= 1)

    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.1f')

    perc_of_SOM= len(np.where(df_nodes_orig.node == iSOM)[0])/len(df_nodes_orig.node)
    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')


plt.tight_layout()

fig.subplots_adjust(bottom=0.1)
cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.015])
cb= fig.colorbar(cf, cax=cbar_ax, orientation="horizontal")

if Svar_full== 'z_ano': Slabelvar= 'Geopotential height anomaly [m]'
elif Svar_full== 'z': Slabelvar= 'Geopotential height [m]'

elif Svar_full== 't_ano': Slabelvar= 'Temperature anomaly [K]'
elif Svar_full== 't': Slabelvar= 'Temperature [K]'

elif Svar_full== 'U_ano': Slabelvar= 'Wind speed anomaly [m/s]'
elif Svar_full== 'U': Slabelvar= 'Wind speed [m/s]'

elif Svar_full== 'q_ano': Slabelvar= 'Specific humidity anomaly [g(kg)]'
elif Svar_full== 'q': Slabelvar= 'Specific humidity [g/kg]'
else: Slabelvar= Svar_full

cb.set_label(Slabelvar, size=14) 

if save:
    if PLCG_type== 'stearing_flow': save_name=savedir+ 'SOM_'+Svar_full+"_"+str(Splevel)+"_"+Obs+"_x"+str(x)+"_y"+str(y)
    elif PLCG_type== 'track_smth': save_name=savedir+ 'SOM_track-smth'+smooth_param+'_'+Svar_full+"_"+str(Splevel)+"_"+Obs+"_x"+str(x)+"_y"+str(y)
    
    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')



"""evolution"""
#fig= plt.figure(fignr, figsize= (2.7*x,2.7*y) )
fig= plt.figure(fignr, figsize= (2.5*x,2.5*y +0) )

fignr+=1
plt.clf()

if imp_rot_nc:
    print('import from all time step list')
    
    if plevel_var != '': plev_str= '_'+ str(plevel_var)
    else: plev_str= plevel_var
    
    if PLCG_type== 'stearing_flow': file= Mediadir + "ERA5_STARS/PL_centred_fields/" +var + plev_str + '_allObs.nc'
    elif PLCG_type== 'track_smth': file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +var + plev_str + '_allObs_track-smth-'+smooth_param+'.nc'


    if os.path.isfile(file):
        ds= xr.open_dataset(file)
        if 'tp' in list(ds.keys()): #for the variables without plevels
            ds['tp']*= 1E3
            ds['tp'].attrs['units']= 'mm'
        if 'cp' in list(ds.keys()): #for the variables without plevels
            ds['cp']*= 1E3
            ds['cp'].attrs['units']= 'mm'
        if 'lsp' in list(ds.keys()): #for the variables without plevels
            ds['lsp']*= 1E3
            ds['lsp'].attrs['units']= 'mm'            
        if 'pres_pv' in list(ds.keys()):
            ds['pres_pv']*= 1E-2
            ds['pres_pv'].attrs['units'] = 'h'+ ds['pres_pv'].attrs['units']            
        if 'z_pv' in list(ds.keys()):
            ds['z_pv']/= 9.81
            ds['z_pv'].attrs['units'] = 'm' 
        if 'slhf' in list(ds.keys()):
            ds['slhf']/= -60**2
            ds['slhf'].attrs['units'] = 'W/m**2'
        if 'sshf' in list(ds.keys()):
            ds['sshf']/= -60**2
            ds['sshf'].attrs['units'] = 'W/m**2'              
        if 'skt' in list(ds.keys()):
            ds['skt'] -= 273.15
            ds['skt'].attrs['units'] = 'C'                 
        if 'cin' in list(ds.keys()):  ds['cin']= ds['cin'].fillna(value= 0)
        if 'ci' in list(ds.keys()):
            if land== 1: ds['ci']= ds['ci'].fillna(value= 1) #replace nan values for the land mask with 1
            if land== 0: ds['ci']= ds['ci'].fillna(value= 0)
            
        
    else:
        ds= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields/" +var + '_all_levs_allObs.nc')
        ds= ds.sel(plev= plevel_var)
        if 'd' in list(ds.keys()): #for the variables with plevels
            ds['d']*= 1E5
            ds['d'].attrs['units'] = '1e-5 '+ ds['d'].attrs['units']
        
    ds= ds.rename({var: var_full}) #this should not include the anomaly

    """do land-sea ice mask""" 
    if lsm_mask:
        ds2= xr.open_dataset( Mediadir + "ERA5_STARS/PL_centred_fields/ci_allObs.nc")
        ds2['ci']= ds2['ci'].fillna(value= 1)
           
        ds[var_full]= ds[var_full].where(ds2['ci'] < 0.3)
      

    """select the time steps from allObs"""
    if Obs_evol != 'allObs':
        ds['ID']= ('time', [str(int(ds.PLnr[n]))[:-3]+'_'+str(int(ds.PLnr[n])%1000) for n in range(len(ds.PLnr))] )
        
        if type(Obs_evol)== int: ds_now= ds.where(ds.Obs== Obs_evol, drop=True)
        elif Obs_evol == 'mature':
            mature_index= [np.where(np.logical_and(ds.ID== Stoll_ind.ID[n], ds.Obs== Stoll_ind.Obs_mature[n]))[0][0] for n in range(len(Stoll_ind.ID)) ]
            ds_now= ds.isel(time= mature_index)
        #    elif Obs_evol == 'allObs':
        
        ds_now= ds_now.dropna(dim='time') #to exclude nan values
        print('excluded PLs due to nan values:', set(Stoll_ind.ID.values) - set(ds_now.ID.values) )
        var_interp= ds_now[var_full].values
        PLID= ds_now['ID'].values
        
        
        if type(Obs_evol)== int: ds_now= ds.where(ds.Obs== Obs_evol, drop=True)
        else:
        #    mature_index= [np.where(np.logical_and(ds.ID== Stoll_ind.ID[n], ds.Obs== Stoll_ind.Obs_mature[n]))[0][0] for n in range(len(Stoll_ind.ID)) ]
            mature_index= [np.where(np.logical_and(ds.ID== Stoll_ind.ID[n], ds.Obs== int(np.round(1+ mat_evol_step/100 *(Stoll_ind.Obs_mature[n] -1))) ))[0][0] for n in range(len(Stoll_ind.ID)) ]
        
            ds_now= ds.isel(time= mature_index)
        
        ds_Obs_evol= ds_now.dropna(dim='time')
    
    
    """add the node to ds - only for allObs"""
    if Obs_evol == 'allObs':
        ds_node_var= np.zeros(len(ds.time))
        
        for s in range(len(df_nodes)):
            df_nodes_s= df_nodes.iloc[s]
            ds_index= np.argwhere(np.logical_and(ds.PLnr.values== df_nodes_s.PLnr, ds.time.values== np.datetime64(df_nodes_s.name[1])))[0][0]
            ds_node_var[ds_index]= df_nodes_s.node
        
        ds['node']= (('time'), ds_node_var)
        ds_Obs_evol= ds
    

    if ano_var == True:
        ds_Obs_evol[var_full]= (('time', 'x', 'y'), ds_Obs_evol[var_full].values- np.nanmean(ds_Obs_evol[var_full].values, axis= (1,2))[:, np.newaxis, np.newaxis])       

plt.clf()


if sym: vextr= scale* np.max([np.max(ds_Obs_evol[var_full]), -np.min(ds_Obs_evol[var_full])])
else: vmax, vmin= scale*float(np.max(ds_Obs_evol[var_full])), scale*float(np.min(ds_Obs_evol[var_full]))



for iSOM in range(1, x*y +1):
    ax1= plt.subplot(y, x, iSOM)
    
    
    if Obs_evol != 'allObs':
        node_list= Stoll_ind[Stoll_ind.node== iSOM].index.values
        ds_node= ds_Obs_evol.where(ds_Obs_evol.ID.isin(node_list), drop=True)
    
    if Obs_evol == 'allObs':
        ds_node= ds_Obs_evol.where(ds_Obs_evol.node == iSOM, drop=True)
    
    if len(levels) > 0:
        cf= plt.contourf(ds_node.x, ds_node.y, np.nanmean(ds_node[var_full], axis= 0), levels= levels, cmap= cmap, extend='both')      
        cs= plt.contour(ds_node.x, ds_node.y, np.nanmean(ds_node[var_full], axis= 0), levels= levels, colors='k', linestyles='solid')

        
    elif sym:
        cf= plt.contourf(ds_node.x, ds_node.y, np.mean(ds_node[var_full], axis= 0), vmin= -vextr, vmax= vextr, cmap= cmap)#, extend='both')      
        cs= plt.contour(ds_node.x, ds_node.y, np.mean(ds_node[var_full], axis= 0), colors='k', linestyles='solid')

    else:
        cf= plt.contourf(ds_node.x, ds_node.y, np.mean(ds_node[var_full], axis= 0), cmap= cmap, vmin= vmin, vmax= vmax)
        cs= plt.contour(ds_node.x, ds_node.y, np.mean(ds_node[var_full], axis= 0), colors='k', vmin= vmin, vmax= vmax)

    if len(levels) > 0: #make the label of the contours
        if np.max(levels) > 5: plt.clabel(cs, fontsize=10, fmt='%1.0f')#, inline=1)
        else: plt.clabel(cs, fontsize=10, fmt='%1.1f')#, inline=1)
    else: plt.clabel(cs, fontsize=10, fmt='%1.1f')

    perc_of_SOM= len(np.where(df_nodes.node == iSOM)[0])/len(df_nodes.node)
    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')

    ax1.scatter(0,0, color='r', zorder= 2)

    ax1.set_yticks(np.arange(-500,501, 250))
    ax1.set_xticks(np.arange(-500,501, 250))
#    ax1.set_xticklabels('')
#|||||||    ax1.set_yticklabels('')    

plt.tight_layout()

fig.subplots_adjust(bottom=0.11)
cbar_ax = fig.add_axes([0.09, 0.06, 0.6, 0.015])
cb= fig.colorbar(cf, cax=cbar_ax, orientation="horizontal")

"""make the arrow"""
arrow_ax = fig.add_axes([0.75, 0.03, 0.2, 0.04])
import matplotlib.patches as mpatches

arrow = mpatches.FancyArrowPatch((0.2, 0.7), (.7, 0.7), mutation_scale=20)
arrow_ax.add_patch(arrow)
arrow_ax.text(0,0 , 'Propagation direction', fontweight= 'bold')#, fontsize=13)
arrow_ax.set_frame_on(False)
arrow_ax.get_xaxis().set_visible(False)
arrow_ax.get_yaxis().set_visible(False)



if 'z_ano' in var_full: labelvar= 'Geopotential height anomaly [m]'
elif 'z' in var_full: labelvar= 'Geopotential height [m]'
elif 't_ano' in var_full: labelvar= 'Temperature anomaly [K]'
elif 't_' in var_full: labelvar= 'Temperature [K]'
elif 'U_ano' in var_full: labelvar= 'Wind speed anomaly [m/s]'
elif 'U_' in var_full: labelvar= 'Wind speed [m/s]'
elif 'q_ano' in var_full: labelvar= 'Specific humidity anomaly [g(kg)]'
elif 'q_' in var_full: labelvar= 'Specific humidity [g/kg]'
elif 'pres_pv' in var_full: labelvar= 'Tropopause level [hPa]'


elif ano_var== False: labelvar= ds[var_full].long_name + ' ['+ ds[var_full].units+ ']'
else: labelvar= ''

cb.set_label(labelvar, size=12) 

if save:
    if PLCG_type== 'stearing_flow':
        save_name=savedir+ 'SOM_'+Svar_full+"_"+str(Splevel)+"_"+Obs+'--mean_'+var_full+'_'+Obs_evol_str+"_x"+str(x)+"_y"+str(y)
    elif PLCG_type== 'track_smth':
        save_name=savedir+ 'SOM_track-smth'+smooth_param+'_'+Svar_full+"_"+str(Splevel)+"_"+Obs+'--mean_'+var_full+'_'+Obs_evol_str+"_x"+str(x)+"_y"+str(y)

    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')



