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
import matplotlib.colors as colors
import numpy as np
import cartopy.crs as ccrs
import scipy.ndimage.filters as filters

from f_useful import *
from f_meteo import *
from f_STARS import *
from f_carto import *

#
save= False
#save= True
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM2/'
#
fignr= 16

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

#PLCG_type= 'track_smth'
smooth_param= '1E-3'

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

compute=False

#ano_var=True
ano_var=False
#var='t'
#var='z'
var='q'
#var, cmap='vo', 'RdBu_r'
#var='d'
#var, cmap='u', 'RdBu'
#var, cmap='v', 'RdBu'

#var, cmap='w', 'RdBu'
#var, cmap='pv', 'Reds'
plevel_var=850


lsm_mask=False

"""variables without plevel"""
#plevel_var=''

#var='lcc'
#var='mcc'
#var='hcc'
#cmap= 'Greys_r'

#var='cp'
#var='lsp'
#cmap= 'Blues'

#var='slhf'
#var='sshf'

#var='pres_pv'
#var='z_pv'
#var='cin'
#var='cape'
#var='skt'
#var='sst'
#var='blh'
#var='msl'
#cmap= 'Reds'

#var, land='ci', 1 #land= [0,1] specifies if land is set to 0 or to 1
#cmap='Blues'

"""variables that needs to be computed"""
#compute=True

#var, imp_vars, plevel_var= 'N', ['t', 'z'], [925, 700]
#var, imp_vars, plevel_var= 't_diff', ['t'], [1000, 500] #first minus second
#var, imp_vars, plevel_var= 'sst-t', ['sst', 't'], 500
#var, imp_vars, plevel_var= 'sst-theta', ['sst', 't'], 500
#var, imp_vars, cmap, plevel_var= 'theta_sst-theta', ['sst', 'msl', 't'],'Reds', 700
#var, imp_vars, cmap, plevel_var= 'theta_e_sst-theta_e', ['sst', 'msl', 't', 'q'], 'Reds', 500


#var, imp_vars, plevel_var= 'theta_diff', ['t'], [925, 500] #first minus second
#var, imp_vars, plevel_var= 'theta_e_diff', ['t', 'q'], [925, 500] #first minus second
#cmap= 'Reds_r'

#var, imp_vars, cmap, plevel_var= 'rh', ['q', 't'], 'Blues',  925
#var, imp_vars, cmap, data_gausfilter, plevel_var= 'grad_t', ['t'], 'Reds', 4, 850
#var, imp_vars, cmap, data_gausfilter, plevel_var= 'adv_t', ['t', 'u', 'v'], 'RdBu_r', 4, 700
#var, imp_vars, cmap, data_gausfilter, plevel_var= 'adv_t_g', ['t', 'z'], 'RdBu_r', 4, 850

#var, imp_vars, cmap, data_gausfilter, plevel_var= 'adv_q', ['q', 'u', 'v'], 'RdBu', 4, 925
#var, imp_vars, cmap, data_gausfilter, plevel_var= 'adv_q_g', ['q', 'z'], 'RdBu', 4, 925

#var, imp_vars, cmap, data_gausfilter, plevel_var= 'u_g', ['z'], 'RdBu_r', '', 925
#var, imp_vars, cmap, data_gausfilter, plevel_var= 'v_g', ['z'], 'RdBu_r', '', 850

#var, imp_vars, cmap, data_gausfilter, plevel_var= 'u_r', ['u' ,'v'], 'RdBu_r', '', 925 #tangential wind
#var, imp_vars, cmap, data_gausfilter, plevel_var= 'v_r', ['u' ,'v'], 'RdBu_r', '', 500 #tangential wind

#var, imp_vars, cmap, plevel_var= 'u_a', ['u' ,'v', 'z'], 'RdBu_r',  850
#var, imp_vars, cmap, plevel_var= 'v_a', ['u' ,'v', 'z'], 'RdBu_r',  500


#var, imp_vars, cmap, plevel_var= 'U', ['u', 'v'], 'Reds', 850
#var, imp_vars, cmap, plevel_var= 'U_g', ['z'], 'Reds',  850
#var, imp_vars, cmap, plevel_var= 'U_r', ['u' ,'v'], 'Reds', 850 #it is the same as U which makes sense
#var, imp_vars, cmap, plevel_var= 'U_a', ['u' ,'v', 'z'], 'Reds',  500


#var, imp_vars, cmap, data_gausfilter, plevel_var= 'vo_g', ['z'], 'RdBu_r', '', 850


var_full= var
if ano_var==True:  var_full= var_full+ '_ano'
if type(plevel_var) == int: var_full= var_full+ '_'+str(plevel_var)
if type(plevel_var) == list:
    if len(plevel_var)== 2: var_full += str(plevel_var[0])+'-'+str(plevel_var[1])
    elif len(plevel_var)== 3: var_full += str(plevel_var[0])+'-'+str(plevel_var[1])+'-'+str(plevel_var[2])




"""specify the levels"""
levels = []
levels_different = False #normally this is false, only if the levels are differently spaced, e.g. logarithmic it should be true - in this case the first and last dispayed level are cut for the "extent" to work

if var_full== 't_1000': levels= np.arange(258, 274.1, 2)
elif var_full== 't_500': levels= np.arange(230, 240.1, 2)
elif var_full== 'z_500': levels= np.arange(5000, 5160.1, 20)
elif var_full== 'z_1000': levels= np.arange(-60, 100.1, 20)
elif var_full== 'q_850': levels= np.arange(0.6, 1.71, .1)
elif var_full== 'q_925': levels= np.arange(0.6, 2.41, .2)
elif var_full== 'q_700': levels= np.arange(0.3, 0.76, .05)
elif var_full== 'q_500': levels= np.arange(0.07, 0.21, .01)


elif var_full== 'z_pv': levels= np.arange(6200, 8400.1, 200)
elif var_full== 'pres_pv': levels= np.arange(300, 421, 20)
elif var_full== 'cin': levels= np.arange(0, 2.1, .25)
elif var_full== 'cape': levels= np.arange(0, 45.1, 5)
elif var_full== 'ci': levels= np.arange(0, 1.01, .1)
elif var_full== 'blh': levels= np.arange(600, 1301, 100)
elif var_full== 'msl': levels= np.arange(992, 1019, 2)

elif var_full== 'slhf': levels= np.arange(0, 200.1, 10)
elif var_full== 'sshf': levels= np.arange(0, 200.1, 10)
elif var_full== 'pv_500':  levels = np.arange(0, 20.1, 1)*10E-8
elif var_full=='skt': levels= np.arange(-8, 8.1, 1)
elif var_full=='sst': levels= np.arange(273, 280.1, 1)

elif var_full== 't_diff1000-500': levels= np.arange(30, 45.1, 1)
elif var_full== 'theta_diff1000-500': levels= np.arange(10, 25.1, 1)
elif var_full== 'theta_diff925-500': levels= np.arange(9, 23.1, 1)
elif var_full== 'theta_e_diff925-500': levels= np.arange(5, 19.1, 1)
elif var_full== 'sst-t_500': levels= np.arange(40, 46.1, 1)
elif var_full== 'sst-theta_500': levels= np.arange(-11, -4.9, 1)
elif var_full== 'theta_e_sst-theta_e_500': levels= np.arange(-1, 9.1, 1)
elif var_full== 'theta_sst-theta_500': levels= np.arange(-11, -4.9, 1)
elif var_full== 'theta_sst-theta_700': levels= np.arange(-4, 2.1, 1)

elif 'adv_t_' in  var_full: levels= np.arange(-1, 1.01, .1)
elif 'adv_q_' in  var_full: levels= np.arange(-0.1, 0.101, .02)
elif 'grad_t' in  var_full: levels= np.arange(1.0, 3.1, .2)

elif 'u_' in var_full or 'v_' in var_full: levels= np.arange(-20, 20.1, 2)
elif 'U_' in var_full: levels= np.arange(0, 20.1, 2)

#elif var_full== 'N925-500': levels= np.arange(0.0042, .0067, 0.0002)

elif 'N' in var_full: levels= np.arange(0.0042, .0067, 0.0002)
#elif 'adv_t_' in  var_full: levels= np.arange(-15, 15, 1)*1E-5


elif 't_ano' in var_full:  levels = np.linspace(-7, 7, 15)
elif 'z_ano' in var_full:  levels = np.arange(-200, 201, 20)
elif 'vo_' in var_full:  levels = np.arange(-50, 50.1, 5)
elif 'pv_' in var_full:  levels = np.arange(0, 7.1, 1)*10E-8

elif 'cp' in var_full or 'lsp' in var_full: levels= np.arange(0, .51, .1) # levels= np.linspace(0, 1, 11)
elif 'tp' in var_full : levels= np.linspace(0, 1, 11)

elif 'cc' in var_full: levels= np.arange(0.2, 1.01, .1)
elif 'd_' in var_full:  levels = np.array([-4, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 4])
elif 'w_' in var_full: levels_different, levels = True, np.array([-1.6, -.8, -.4, -.2, -.1, -.05, 0.05, .1, .2, .4, .8, 1.6])



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
    save_name=savedir+ 'SOM_track-smth'+smooth_param+'_'+Svar_full+"_"+str(Splevel)+"_"+Obs+"_x"+str(x)+"_y"+str(y)
    
    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')



"""evolution"""
#fig= plt.figure(fignr, figsize= (2.7*x,2.7*y) )
fig= plt.figure(fignr, figsize= (2.5*x,2.5*y +0) )

fignr+=1
plt.clf()

if imp_rot_nc:
    print('import from all time step list')
    
    if compute== False:
        if plevel_var != '': plev_str= '_'+ str(plevel_var)
        else: plev_str= plevel_var
        
        file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +var + plev_str + '_allObs_track-smth-'+smooth_param+'.nc'
    
    
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
            if 'msl' in list(ds.keys()):
                ds['msl'] /= 100
                ds['msl'].attrs['units'] = 'hPa'                 
            if 'cin' in list(ds.keys()):  ds['cin']= ds['cin'].fillna(value= 0)
            if 'ci' in list(ds.keys()):
                if land== 1: ds['ci']= ds['ci'].fillna(value= 1) #replace nan values for the land mask with 1
                if land== 0: ds['ci']= ds['ci'].fillna(value= 0)
                
            
        else:
            ds= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +var + '_all_levs_allObs_track-smth-'+smooth_param+'.nc')
            ds= ds.sel(plev= plevel_var)
            
            if 'd' in list(ds.keys()): #for the variables with plevels
                ds['d']*= 1E5
                ds['d'].attrs['units'] = '1e-5 '+ ds['d'].attrs['units']
            
        ds= ds.rename({var: var_full}) #this should not include the anomaly


    elif compute:
        for ni, i_var in enumerate(imp_vars):
            file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +i_var +'_allObs_track-smth-'+smooth_param+'.nc'
            if os.path.isfile(file):
                ds1= xr.open_dataset(file)
            else:
                ds1= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +i_var + '_all_levs_allObs_track-smth-'+smooth_param+'.nc')
                ds1= ds1.sel(plev= plevel_var)
                
            if ni == 0: ds= ds1
            else: ds= xr.merge([ds, ds1])                

        if var=='t_diff':
            ds[var_full]= (('time', 'x', 'y'), ds['t'].isel(plev=0) - ds['t'].isel(plev=1) )
            ds[var_full].attrs['units'] = 'K'

        if var=='sst-t':
            ds[var_full]= (('time', 'x', 'y'), ds['sst'] - ds['t'])
            ds[var_full].attrs['units'] = 'K'            

        if var=='sst-theta':
            ds[var_full]= (('time', 'x', 'y'), ds['sst'] - PotTemp(ds.t, plevel_var) )
            ds[var_full].attrs['units'] = 'K'     
            ds[var_full].attrs['long_name']='SST - theta'+str(plevel_var)

        if var=='theta_sst-theta':
            theta= PotTemp(ds.t, plevel_var)
            theta_sst= PotTemp(ds.sst, ds.msl/100)

            ds[var_full]= (('time', 'x', 'y'), theta_sst - theta )
            ds[var_full].attrs['units'] = 'K'     
            ds[var_full].attrs['long_name']='theta_sst-theta'+str(plevel_var)


        if var=='theta_e_sst-theta_e':
            theta_e= EquiPotTemp(ds.t, ds.q/1E3, plevel_var)
#            theta_e_sst= EquiPotTemp(ds.sst, RH2SH(100, 1000, ds.sst)/1E3, 1000)
            theta_e_sst= EquiPotTemp(ds.sst, RH2SH(100, ds.msl/100, ds.sst)/1E3, ds.msl/100)

            ds[var_full]= (('time', 'x', 'y'), theta_e_sst - theta_e )
            ds[var_full].attrs['units'] = 'K'     
            ds[var_full].attrs['long_name']='theta_e_sst-theta_e'+str(plevel_var)


        if var=='theta_diff':
            theta_1, theta_0= PotTemp(ds.t.isel(plev=1), plevel_var[1]), PotTemp(ds.t.isel(plev=0), plevel_var[0])

            ds[var_full]= (('time', 'x', 'y'), theta_1 - theta_0)
            ds[var_full].attrs['units'] = 'K'

        if var=='theta_e_diff':
            theta_e_1, theta_e_0= EquiPotTemp(ds.t.isel(plev=1), ds.q.isel(plev=1)/1E3, plevel_var[1]), EquiPotTemp(ds.t.isel(plev=0), ds.q.isel(plev=0)/1E3, plevel_var[0])

            ds[var_full]= (('time', 'x', 'y'), theta_e_1 - theta_e_0)
            ds[var_full].attrs['units'] = 'K'

        if var=='rh':
            rh= SH2RH(ds.q, plevel_var, ds.t)

            ds[var_full]= (('time', 'x', 'y'), rh)
            ds[var_full].attrs['units'] = '%'
            ds[var_full].attrs['long_name']='Relative humidity'
            
            
        if var=='N':
            h_diff=  (ds['z'].isel(plev=1) - ds['z'].isel(plev=0) )

            theta_1, theta_0= PotTemp(ds.t.isel(plev=1), plevel_var[1]), PotTemp(ds.t.isel(plev=0), plevel_var[0])
            N= np.sqrt(9.81/(theta_1+theta_0)/2 *  (theta_1- theta_0)/h_diff)        

            ds[var_full]= (('time', 'x', 'y'),  N)
            ds[var_full].attrs['units']= '1/s' 

        if var== 'grad_t':
            if data_gausfilter:
                ds['t'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.t, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
                ds['t'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.t, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )

            grad_t_vec= np.gradient(ds.t, 0.25) #two components of the temperature gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km
            ds[var_full]= (('time', 'x', 'y'), np.sqrt(grad_t_vec[0]**2 + grad_t_vec[1]**2) )
            ds[var_full].attrs['units']= 'K/100km'

        if var== 'u_g':
#            if data_gausfilter:
#                ds['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.z, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
#                ds['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.z, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )
#
            f= np.sin(np.deg2rad(Stoll.lat.values)) * 2*2*np.pi/(24*60**2)

            grad_z_vec = np.gradient(ds.z *9.81 , 25E3)
            u_g= -1/f[:, np.newaxis, np.newaxis]* grad_z_vec[1]
            
            ds[var_full]= (('time', 'x', 'y'), u_g )
            ds[var_full].attrs['units'] = 'm/s'
            ds[var_full].attrs['long_name']='Tangential geostrophic wind'
            
        if var== 'v_g':
#            if data_gausfilter:
#                ds['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.z, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
#                ds['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.z, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )

            f= np.sin(np.deg2rad(Stoll.lat.values)) * 2*2*np.pi/(24*60**2)
            grad_z_vec = np.gradient(ds.z *9.81 , 25E3)
            v_g= 1/f[:, np.newaxis, np.newaxis]* grad_z_vec[2] 
            
            ds[var_full]= (('time', 'x', 'y'), v_g )
            ds[var_full].attrs['units'] = 'm/s'
            ds[var_full].attrs['long_name']='Azimuthal geostrophic wind'            

        if var== 'vo_g':
#            if data_gausfilter:
#                ds['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.z, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
#                ds['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.z, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )

            f= np.sin(np.deg2rad(Stoll.lat.values)) * 2*2*np.pi/(24*60**2)
            grad_z_vec = np.gradient(ds.z *9.81 , 25E3)
            u_g= -1/f[:, np.newaxis, np.newaxis]* grad_z_vec[1]
            v_g= 1/f[:, np.newaxis, np.newaxis]* grad_z_vec[2] 
            
            vo_g = np.gradient(v_g, 25E3, axis= 2) - np.gradient(u_g, 25E3, axis= 1)
            
            ds[var_full]= (('time', 'x', 'y'), vo_g*1E5 )
            ds[var_full].attrs['units'] = '1/s'
            ds[var_full].attrs['long_name']='Geostrophic vorticity'  


        if var== 'u_r':
            U= np.sqrt(ds.u**2 + ds.v**2)
            wind_beering= UV2Direction(ds.u, ds.v)
            track_beering = ds.beering.values
            
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='') #maybe u_zonal and v_meridional are exchanged in WindSpeedDirection2UV, but this seem to work"""

            ds[var_full]= (('time', 'x', 'y'), u_r )
            ds[var_full].attrs['units'] = 'm/s'
            ds[var_full].attrs['long_name']='Tangential wind'

        if var== 'v_r':
            U= np.sqrt(ds.u**2 + ds.v**2)
            wind_beering= UV2Direction(ds.u, ds.v)
            track_beering = ds.beering.values
            
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='')
            ds[var_full]= (('time', 'x', 'y'), v_r )
            ds[var_full].attrs['units'] = 'm/s'
            ds[var_full].attrs['long_name']='Azimuthal wind'

        if var== 'u_a':
            U= np.sqrt(ds.u**2 + ds.v**2)
            wind_beering= UV2Direction(ds.u, ds.v)
            track_beering = ds.beering.values           
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='') #maybe u_zonal and v_meridional are exchanged in WindSpeedDirection2UV, but this seem to work"""

            f= np.sin(np.deg2rad(Stoll.lat.values)) * 2*2*np.pi/(24*60**2)
            grad_z_vec = np.gradient(ds.z *9.81 , 25E3)
            u_g= -1/f[:, np.newaxis, np.newaxis]* grad_z_vec[1]
#            v_g= 1/f[:, np.newaxis, np.newaxis]* grad_z_vec[2] 
            
            u_a= u_r - u_g
#            v_a= v_r - v_g
            
            ds[var_full]= (('time', 'x', 'y'), u_a )
            ds[var_full].attrs['units'] = 'm/s'
            ds[var_full].attrs['long_name']='Ageostrophic tangential wind' 

        if var== 'v_a':
            U= np.sqrt(ds.u**2 + ds.v**2)
            wind_beering= UV2Direction(ds.u, ds.v)
            track_beering = ds.beering.values           
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='') #maybe u_zonal and v_meridional are exchanged in WindSpeedDirection2UV, but this seem to work"""

            f= np.sin(np.deg2rad(Stoll.lat.values)) * 2*2*np.pi/(24*60**2)
            grad_z_vec = np.gradient(ds.z *9.81 , 25E3)
            v_g= 1/f[:, np.newaxis, np.newaxis]* grad_z_vec[2] 
            
            v_a= v_r - v_g
            
            ds[var_full]= (('time', 'x', 'y'), v_a )
            ds[var_full].attrs['units'] = 'm/s'
            ds[var_full].attrs['long_name']='Ageostrophic azimuthal wind' 
            
        if var== 'U':
            U= np.sqrt(ds.u**2 + ds.v**2)

            ds[var_full]= (('time', 'x', 'y'), U )
            ds[var_full].attrs['units'] = 'm/s'
            ds[var_full].attrs['long_name']= 'Wind speed'

#        if var== 'U_r':
#            U= np.sqrt(ds.u**2 + ds.v**2)
#            wind_beering= UV2Direction(ds.u, ds.v)
#            track_beering = ds.beering.values
#            
#            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
#            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='') #maybe u_zonal and v_meridional are exchanged in WindSpeedDirection2UV, but this seem to work"""
#
#            ds[var_full]= (('time', 'x', 'y'), np.sqrt(u_r**2 + v_r**2) )
#            ds[var_full].attrs['units'] = 'm/s'
#            ds[var_full].attrs['long_name']='Rotated wind speed'

        if var== 'U_g':
            f= np.sin(np.deg2rad(Stoll.lat.values)) * 2*2*np.pi/(24*60**2)
            grad_z_vec = np.gradient(ds.z *9.81 , 25E3)
            u_g= -1/f[:, np.newaxis, np.newaxis]* grad_z_vec[1]
            v_g= 1/f[:, np.newaxis, np.newaxis]* grad_z_vec[2] 
            
            ds[var_full]= (('time', 'x', 'y'), np.sqrt(u_g**2 + v_g**2) )
            ds[var_full].attrs['units'] = 'm/s'
            ds[var_full].attrs['long_name']='Geostrophic wind speed' 

        if var== 'U_a':
            U= np.sqrt(ds.u**2 + ds.v**2)
            wind_beering= UV2Direction(ds.u, ds.v)
            track_beering = ds.beering.values           
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='') #maybe u_zonal and v_meridional are exchanged in WindSpeedDirection2UV, but this seem to work"""

            f= np.sin(np.deg2rad(Stoll.lat.values)) * 2*2*np.pi/(24*60**2)
            grad_z_vec = np.gradient(ds.z *9.81 , 25E3)
            u_g= -1/f[:, np.newaxis, np.newaxis]* grad_z_vec[1]
            v_g= 1/f[:, np.newaxis, np.newaxis]* grad_z_vec[2] 
            
            u_a= u_r - u_g
            v_a= v_r - v_g
            
            ds[var_full]= (('time', 'x', 'y'), np.sqrt(u_a**2 + v_a**2) )
            ds[var_full].attrs['units'] = 'm/s'
            ds[var_full].attrs['long_name']='Ageostrophic wind speed'  

           
        if var== 'adv_t_g':
            if data_gausfilter:
                ds['t'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.t, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
                ds['t'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.t, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )
#                ds['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.z, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
#                ds['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.z, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )

            grad_t_vec= np.gradient(ds.t, 25E3) #two components of the temperature gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km
            f= np.sin(np.deg2rad(Stoll.lat.values)) * 2*2*np.pi/(24*60**2)
            grad_z_vec = np.gradient(ds.z *9.81 , 25E3)
            u_g= -1/f[:, np.newaxis, np.newaxis]* grad_z_vec[1]
            v_g= 1/f[:, np.newaxis, np.newaxis]* grad_z_vec[2] 
            
            ds[var_full]= (('time', 'x', 'y'), -1* (u_g * grad_t_vec[2] + v_g * grad_t_vec[1] ) * 3600 )
            ds[var_full].attrs['units'] = 'K/h'
            ds[var_full].attrs['long_name']='Horizontal geostrophic temperature advection' 


        if var== 'adv_t':
            if data_gausfilter:
                ds['t'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.t, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
                ds['t'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.t, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )

            grad_t_vec= np.gradient(ds.t, 25E3) #two components of the temperature gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km

            U= np.sqrt(ds.u**2 + ds.v**2)
            wind_beering= UV2Direction(ds.u, ds.v)
            track_beering = ds.beering.values
            
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='')
            
            ds[var_full]= (('time', 'x', 'y'), -1* (u_r * grad_t_vec[2] + v_r * grad_t_vec[1] ) * 3600 )
            ds[var_full].attrs['units'] = 'K/h'
            ds[var_full].attrs['long_name']='Horizontal temperature advection' 



      
        if var== 'adv_q':
            if data_gausfilter:
                ds['q'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.q, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
                ds['q'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.q, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )
#                ds['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.z, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
#                ds['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.z, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )

            grad_q_vec= np.gradient(ds.q, 25E3) #two components of the temperature gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km
            
            U= np.sqrt(ds.u**2 + ds.v**2)
            wind_beering= UV2Direction(ds.u, ds.v)
            track_beering = ds.beering.values
            
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='')
            
            ds[var_full]= (('time', 'x', 'y'), -1* (u_r * grad_q_vec[2] + v_r * grad_q_vec[1] ) * 3600 )
            ds[var_full].attrs['units'] = '(g/kg)/h'
            ds[var_full].attrs['long_name']='Horizontal humidity advection' 



        if var== 'adv_q_g':
            if data_gausfilter:
                ds['q'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.q, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
                ds['q'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.q, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )

            grad_q_vec= np.gradient(ds.q, 25E3) #two components of the temperature gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km
            f= np.sin(np.deg2rad(Stoll.lat.values)) * 2*2*np.pi/(24*60**2)
            grad_z_vec = np.gradient(ds.z *9.81 , 25E3)
            u_g= -1/f[:, np.newaxis, np.newaxis]* grad_z_vec[1]
            v_g= 1/f[:, np.newaxis, np.newaxis]* grad_z_vec[2] 
            
            ds[var_full]= (('time', 'x', 'y'), -1* (u_g * grad_q_vec[2] + v_g * grad_q_vec[1] ) * 3600 )
            ds[var_full].attrs['units'] = '(g/kg)/h'
            ds[var_full].attrs['long_name']='Horizontal humidity advection' 

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


print('plot the SOMs')

for iSOM in range(1, x*y +1):
#for iSOM in range(1, 2):    
    print('SOM ', iSOM)
    ax1= plt.subplot(y, x, iSOM)
    
    
    if Obs_evol != 'allObs':
        node_list= Stoll_ind[Stoll_ind.node== iSOM].index.values
        ds_node= ds_Obs_evol.where(ds_Obs_evol.ID.isin(node_list), drop=True)
    
    if Obs_evol == 'allObs':
        ds_node= ds_Obs_evol.where(ds_Obs_evol.node == iSOM, drop=True)
    
    if levels_different:
        norm = colors.BoundaryNorm(boundaries=levels, ncolors=256) #this is for the case when the levels are not equally distanced, e.g. for w where I also want to display small values
        cf= plt.contourf(ds_node.x, ds_node.y, np.nanmean(ds_node[var_full], axis= 0), levels= levels[1:-1], cmap= cmap, norm=norm, extend='both')      
        cs= plt.contour(ds_node.x, ds_node.y, np.nanmean(ds_node[var_full], axis= 0), levels= levels[1:-1], colors='k', linestyles='solid')

    elif len(levels) > 0:
        cf= plt.contourf(ds_node.x, ds_node.y, np.nanmean(ds_node[var_full], axis= 0), levels= levels, cmap= cmap, extend='both')      
        cs= plt.contour(ds_node.x, ds_node.y, np.nanmean(ds_node[var_full], axis= 0), levels= levels, colors='k', linestyles='solid')
       
    elif sym:
        cf= plt.contourf(ds_node.x, ds_node.y, np.mean(ds_node[var_full], axis= 0), vmin= -vextr, vmax= vextr, cmap= cmap)#, extend='both')      
        cs= plt.contour(ds_node.x, ds_node.y, np.mean(ds_node[var_full], axis= 0), colors='k', linestyles='solid')

    else:
        cf= plt.contourf(ds_node.x, ds_node.y, np.mean(ds_node[var_full], axis= 0), cmap= cmap, vmin= vmin, vmax= vmax)
        cs= plt.contour(ds_node.x, ds_node.y, np.mean(ds_node[var_full], axis= 0), colors='k', vmin= vmin, vmax= vmax)

    if len(levels) > 0: #make the label of the contours
        if np.max(abs(levels)) > 5: plt.clabel(cs, fontsize=10, fmt='%1.0f')#, inline=1)
        elif np.max(abs(levels)) > 0.5 and var != 'w': plt.clabel(cs, fontsize=10, fmt='%1.1f')#, inline=1)
        elif np.max(abs(levels)) > 0.05 or var== 'w': plt.clabel(cs, fontsize=10, fmt='%1.2f')#, inline=1)
        
        elif np.max(abs(levels)) > 0.005 : plt.clabel(cs, fontsize=10, fmt='%1.4f')#, inline=1)

    else: plt.clabel(cs, fontsize=10, fmt='%1.1f')

    perc_of_SOM= len(np.where(df_nodes.node == iSOM)[0])/len(df_nodes.node)
    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')

    ax1.scatter(0,0, color='r', zorder= 2)

    ax1.set_yticks(np.arange(-500,501, 250))
    ax1.set_xticks(np.arange(-500,501, 250))
#    ax1.set_xticklabels('')
#    ax1.set_yticklabels('')    

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
elif 't_diff' in var_full: labelvar= 'Temperature difference '+str(plevel_var[0])+'-'+str(plevel_var[1])+'hPa [K]'
elif 'theta_diff' in var_full: labelvar= 'Pot. temp. diff. '+str(plevel_var[0])+'-'+str(plevel_var[1])+'hPa [K]'
elif 'theta_e_diff' in var_full: labelvar= 'Equiv. pot. temp. diff. '+str(plevel_var[0])+'-'+str(plevel_var[1])+'hPa [K]'
elif 'sst-t_' in var_full: labelvar= 'SST - T'+str(plevel_var)+' [K]'

    
elif 'N' in var_full: labelvar= 'N '+str(plevel_var[0])+'-'+str(plevel_var[1])+'hPa [1/s]'
elif 'grad_t' in var_full: labelvar= r"$\nabla_h$ T "+str(plevel_var)+' [K/100km]'
elif 'adv_t' in var_full: labelvar= "Horizontal temperature advection "+str(plevel_var)+' [K/h]'
elif 'adv_q' in var_full: labelvar= "Horizontal humidity advection "+str(plevel_var)+' [(g/kg)/h]'

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
    save_name=savedir+ 'SOM_track-smth'+smooth_param+'_'+Svar_full+"_"+str(Splevel)+"_"+Obs+'--mean_'+var_full+'_'+Obs_evol_str+"_x"+str(x)+"_y"+str(y)

    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')



