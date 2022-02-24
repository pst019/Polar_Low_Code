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

#
save= True
#save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM2/'
#
fignr= 7

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


"""shading variable"""
ano_var=True
#ano_var=False
var='t'
plevel_var=925
cmap= 'RdBu_r'


var_full= var
if ano_var==True:  var_full += '_ano'
if type(plevel_var) == int: var_full += '_'+str(plevel_var)

if 't_ano' in var_full:  levels = np.linspace(-7, 7, 15)
elif 'pres_pv' in var_full: levels= np.arange(300, 421, 20)
elif 'tp' in var_full: levels= np.linspace(0, 1, 11)
elif 'cc' in var_full: levels= np.arange(0.2, 1.01, .1)
elif 'd' in var_full:  levels = [-4, -3, -2, -1, 1, 2, 3, 4]


"""contour variable"""
ano_cont=True
cont='z'
plevel_cont=925


cont_full= cont
if ano_cont==True:  cont_full += '_ano'
if type(plevel_cont) == int: cont_full += '_'+str(plevel_cont)

if 'z_ano' in cont_full: cont_levels = np.arange(-200, 201, 20)  

var_list= [var, cont]
var_full_list= [var_full, cont_full]


"""contour variable 2"""
plot_cont2=True
ano_cont2=False
cont2='d'
plevel_cont2=925

#ano_cont2=False
#cont2='cape'
#plevel_cont2=''

if plot_cont2:
    cont2_full= cont2
    if ano_cont2==True:  cont2_full += '_ano'
    if type(plevel_cont2) == int: cont2_full += '_'+str(plevel_cont2)
    
    if cont2 == 'd': cont_levels2 = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]  
    if cont2 == 'cape': cont_levels2= [30, 40, 50, 60] #here warnings are displayed since somethimes no contours are found
    
    var_list+= [cont2]
    var_full_list += [cont2_full]
#
#
#ano_cont3=False
##cont3='blh'
##cont3='cape'
##cont3='tp'
##cont3='lcc'
#cont3='hcc'
#
#
#plevel_cont3=False
##cont2_p= cont2+ '_'+str(plevel_cont2)
#
#
#if ano_cont3==True:
#    cont3_full= cont3+ '_ano'
#else:
#    cont3_full= cont3
##    cont2_p_full= cont2_p





#if cont3_full == 'blh': cont_levels3 = np.arange(0, 1501, 100) 
#if cont3_full == 'cape': cont_levels3= [30, 40, 50, 60] #here warnings are displayed since somethimes no contours are found
#if cont3_full == 'tp': cont_levels3= np.arange(0.2, 1, .2) #[30, 40, 50, 60] #here warnings are displayed since somethimes no contours are found
#if 'cc' in cont3_full: cont_levels3= np.arange(0.2, 1, .1) #[30, 40, 50, 60] #here warnings are displayed since somethimes no contours are found




 
    
    
"""import Stoll list"""
if imp_list:
#    SOM_filedir= Mediadir + "ERA5_STARS/PL_centred_fields/SOM/"+Svar_full+"_"+str(Splevel)+'_'+Obs+"_x"+str(x)+"_y"+str(y)
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
#        Stoll= pd.concat([Stoll, df_nodes], axis= 1)

        """to rotate/mirror the nodes"""
        df_nodes_orig= df_nodes.copy()
        if mirror_top_bottom:
            df_nodes['node']= (y-1 -(df_nodes.node.values-1)//x)*x + (df_nodes.node.values -1)%x +1
        if rotate_clock:
            x_n= (df_nodes['node'].values-1)%x +1
            y_n= (df_nodes['node'].values-1)//x +1
            y_n= y+1 - y_n #reverse the y_nodes
            
            df_nodes['node']= (x_n -1)*x + y_n


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

    """file for shading variable"""
    if plevel_var != '': plev_str= '_'+ str(plevel_var)
    else: plev_str= plevel_var
    
    if PLCG_type== 'stearing_flow': file= Mediadir + "ERA5_STARS/PL_centred_fields/" +var + plev_str + '_allObs.nc'
    elif PLCG_type== 'track_smth': file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +var + plev_str + '_allObs_track-smth-'+smooth_param+'.nc'


    if os.path.isfile(file): ds= xr.open_dataset(file)
    else:
        ds= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields/" +var + '_all_levs_allObs.nc')
        ds= ds.sel(plev= plevel_var)
    ds= ds.rename({var: var_full})
    

    """file for cont variable"""
    if plevel_cont != '': plev_str= '_'+ str(plevel_cont)
    else: plev_str= plevel_cont
    
    if PLCG_type== 'stearing_flow': file= Mediadir + "ERA5_STARS/PL_centred_fields/" +cont + plev_str + '_allObs.nc'
    elif PLCG_type== 'track_smth': file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +cont + plev_str + '_allObs_track-smth-'+smooth_param+'.nc'


    if os.path.isfile(file): ds2= xr.open_dataset(file)
    else:
        ds2= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields/" +cont + '_all_levs_allObs.nc')
        ds2= ds2.sel(plev= plevel_cont)
    
    ds2= ds2.rename({cont: cont_full})
    ds= xr.merge([ds, ds2])
    
    
    """file for cont variable 2"""
    if plot_cont2:
        if plevel_cont2 != '': plev_str= '_'+ str(plevel_cont2)
        else: plev_str= plevel_cont2
        
        if PLCG_type== 'stearing_flow': file= Mediadir + "ERA5_STARS/PL_centred_fields/" +cont2 + plev_str + '_allObs.nc'
        elif PLCG_type== 'track_smth': file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +cont2 + plev_str + '_allObs_track-smth-'+smooth_param+'.nc'
    
    
        if os.path.isfile(file): ds2= xr.open_dataset(file)
        else:
            ds2= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields/" +cont2 + '_all_levs_allObs.nc')
            ds2= ds2.sel(plev= plevel_cont2)
            ds2= ds2.rename({'plev': 'plev2'})
        ds2= ds2.rename({cont2: cont2_full})
        ds= xr.merge([ds, ds2])
#    
#    
#    
#    ds2= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields/" +cont3 + '_allObs.nc')
#    ds= xr.merge([ds, ds2])


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
        var_interp= ds_now[var].values
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
        ds_Obs_evol[var_full]= (('time', 'x', 'y'), ds_Obs_evol[var_full].values- np.mean(ds_Obs_evol[var_full].values, axis= (1,2))[:, np.newaxis, np.newaxis])
    if ano_cont == True:
        ds_Obs_evol[cont_full]= (('time', 'x', 'y'), ds_Obs_evol[cont_full].values- np.mean(ds_Obs_evol[cont_full].values, axis= (1,2))[:, np.newaxis, np.newaxis])
#    if ano_cont2 == True:
#        ds_Obs_evol[cont2_p_full]= (('time', 'x', 'y'), ds_Obs_evol[cont2_p].values- np.mean(ds_Obs_evol[cont2_p].values, axis= (1,2))[:, np.newaxis, np.newaxis])
#    if ano_cont3 == True:
#        ds_Obs_evol[cont3_full]= (('time', 'x', 'y'), ds_Obs_evol[cont3_full].values- np.mean(ds_Obs_evol[cont3_full].values, axis= (1,2))[:, np.newaxis, np.newaxis])
#        
    if 'tp' in list(ds_Obs_evol.keys()):
        ds_Obs_evol['tp']*= 1E3
        ds_Obs_evol['tp'].attrs['units']= 'mm'

    if 'd' in var_list:
        ind= np.where(np.array(var_list) == 'd')[0][0]
        ds_Obs_evol[var_full_list[ind]]*= 1E5
        ds_Obs_evol[var_full_list[ind]].attrs['units'] = '1e-5 '+ ds_Obs_evol[var_full_list[ind]].attrs['units']


plt.clf()


for iSOM in range(1, x*y +1):
    ax1= plt.subplot(y, x, iSOM)
    
    
    if Obs_evol != 'allObs':
        node_list= Stoll_ind[Stoll_ind.node== iSOM].index.values
        ds_node= ds_Obs_evol.where(ds_Obs_evol.ID.isin(node_list), drop=True)
    
    if Obs_evol == 'allObs':
        ds_node= ds_Obs_evol.where(ds_Obs_evol.node == iSOM, drop=True)


    cf= plt.contourf(ds_node.x, ds_node.y, np.mean(ds_node[var_full], axis= 0), levels= levels, cmap= cmap, extend='both')      
    
    cs= plt.contour(ds_node.x, ds_node.y, np.mean(ds_node[cont_full], axis= 0), levels= cont_levels, colors='k')#, linestyles='solid')

    if plot_cont2:
        cs2= plt.contour(ds_node.x, ds_node.y, np.mean(ds_node[cont2_full], axis= 0), levels= cont_levels2, colors= 'g')
      

#        cont_levels3= np.arange(0, 76, 5) #[1000,  1200, 1400]
#        cf= plt.contourf(ds_node.x, ds_node.y, np.average(ds_node[cont3_full], axis= 0), levels= cont_levels3, cmap= cmap, extend='both')
#        cs3= plt.contour(ds_node.x, ds_node.y, np.average(ds_node[cont3_full], axis= 0), levels= cont_levels3, colors= 'g')



    if np.max(cont_levels) > 3: plt.clabel(cs, cont_levels[::2], fontsize=10, fmt='%1.0f')#, inline=1)
    else: plt.clabel(cs, cont_levels[::2], fontsize=10, fmt='%1.1f')#, inline=1)

    if plot_cont2:
        if np.max(np.mean(ds_node[cont2_full], axis= 0) ) > min(cont_levels2):
            if max(cont_levels2) > 3: plt.clabel(cs2, cont_levels2[::2], fontsize=10, fmt='%1.0f')
            else: plt.clabel(cs2, cont_levels2[::2], fontsize=10, fmt='%1.1f')


    perc_of_SOM= len(np.where(df_nodes.node == iSOM)[0])/len(df_nodes.node)
    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')

    ax1.scatter(0,0, color='r')
    

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
elif 't' in var_full and var_full != 'tp': labelvar= 'Temperature [K]'
elif 'U_ano' in var_full: labelvar= 'Wind speed anomaly [m/s]'
elif 'U' in var_full: labelvar= 'Wind speed [m/s]'
elif 'q_ano' in var_full: labelvar= 'Specific humidity anomaly [g(kg)]'
elif 'q' in var_full: labelvar= 'Specific humidity [g/kg]'
elif 'pres_pv' in var_full: labelvar= 'Tropopause level [hPa]'

elif ano_var== False: labelvar= ds[var_full].long_name + ' ['+ ds[var_full].units+ ']'
else: labelvar= ''
cb.set_label(labelvar, size=12) 


if save:
    if PLCG_type== 'stearing_flow':
        save_name=savedir+ 'SOM_'+Svar_full+"_"+str(Splevel)+"_"+Obs+'--mean_'+var_full+"_"+cont_full
    elif PLCG_type== 'track_smth':
        save_name=savedir+ 'SOM_track-smth'+smooth_param+'_'+Svar_full+"_"+str(Splevel)+"_"+Obs+'--mean_'+var_full+"_"+cont_full


    if plot_cont2: save_name+= "_"+cont2_full 
    
    save_name += '_'+Obs_evol_str+"_x"+str(x)+"_y"+str(y)

    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')



