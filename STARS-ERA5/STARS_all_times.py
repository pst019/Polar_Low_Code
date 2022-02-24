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
#    homedir= '/home/'+user+'/home/'
#    Mediadir= '/media/'+user+'/PatsOrange/'
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
#    Mediadir= '/run/media/pst019/1692A00D929FEF8B/'
    Mediadir= '/run/media/pst019/PatsOrange/'
    
homedir=Mediadir+'home/'


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
save= True
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM_composites/'
#
fignr= 24

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



#these are not important and the script could be cleaned
Splevel=850
Sano=True
#Sano= False
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
Obs_evol='allObs'



calculate='mean'
calculate='std'


compute=False

cmap= 'RdBu'
ano_var=True
#ano_var=False
var, cmap='t', 'RdBu'
#var, cmap='z', 'RdBu'
#var, cmap='q', 'Blues'
#var, cmap='vo', 'RdBu_r'
#var='d'
#var, cmap='u', 'RdBu'
#var, cmap='v', 'RdBu'

#var, cmap='w', 'RdBu'
#var, cmap='pv', 'Reds'
plevel_var=850


lsm_mask=False
size= 500 #250
contour_numbers= True #plot the numbers of the contours in the plot
only_contours=False
data_gausfilter=False

vo_tend_thresh= False
#vo_tend_thresh= 0 #1E-5

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
#var, cmap= 'tcwv', 'RdBu_r'
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
#var, imp_vars, cmap, plevel_var= 'theta_sst-theta', ['sst', 'msl', 't'],'Reds', 925
#var, imp_vars, cmap, plevel_var= 'theta_e_sst-theta_e', ['sst', 'msl', 't', 'q'], 'Reds', 700
#var, imp_vars, cmap, plevel_var= 'theta_e_2m-theta_e', ['2t', 'msl', 't', 'q'], 'Reds', 500

#var, imp_vars, cmap, plevel_var= 'div', ['u', 'v'], 'RdBu_r', 850 #the self calculated divergence - to test if I get the same as for the ERA-5 divergence


#var, imp_vars, plevel_var= 'theta_diff', ['t'], [925, 500] #first minus second
#var, imp_vars, plevel_var= 'theta_e_diff', ['t', 'q'], [925, 500] #first minus second
#cmap= 'Reds_r'

#var, imp_vars, cmap, plevel_var= 'rh', ['q', 't'], 'Blues',  925
#var, imp_vars, cmap, data_gausfilter, plevel_var= 'grad_t', ['t'], 'Reds', 4, 850
#var, imp_vars, cmap, data_gausfilter, plevel_var= 'adv_t', ['t', 'u', 'v'], 'RdBu_r', 4, 700
#var, imp_vars, cmap, data_gausfilter, plevel_var= 'adv_t_g', ['t', 'z'], 'RdBu_r', 4, 850

#var, imp_vars, cmap, data_gausfilter, plevel_var= 'adv_q', ['q', 'u', 'v'], 'RdBu', 4, 925
#var, imp_vars, cmap, data_gausfilter, plevel_var= 'adv_q_g', ['q', 'z'], 'RdBu', 4, 925


#var, imp_vars, cmap, plevel_var, contour_numbers= 'adv_rel_vo', ['vo', 'u', 'v'], 'RdBu_r', 850, False
#var, imp_vars, cmap, plevel_var, contour_numbers= 'adv_plan_vo', ['lat', 'u', 'v'], 'RdBu_r', 850, False

#var, imp_vars, cmap, plevel_var, contour_numbers= 'stretch_vo', ['vo', 'd', 'lat'], 'RdBu_r', 850, False
#var, imp_vars, cmap, plevel_var, contour_numbers= 'stretch_rel_vo', ['vo', 'd'], 'RdBu_r', 850, False

#var, imp_vars, cmap, plevel_var, contour_numbers= 'tilt_vo', ['u', 'v', 'w'], 'RdBu_r', [925, 850, 700], False
#var, imp_vars, cmap, plevel_var, contour_numbers= 'vert_adv_vo', ['vo', 'w'], 'RdBu_r', [925, 850, 700], False
#var, imp_vars, cmap, plevel_var= 'stretch_vo_self', ['vo', 'u', 'v'], 'RdBu_r', 850


#var, imp_vars, cmap, data_gausfilter, plevel_var= 'u_g', ['z'], 'RdBu_r', '', 925
#var, imp_vars, cmap, data_gausfilter, plevel_var= 'v_g', ['z'], 'RdBu_r', '', 850

#var, imp_vars, cmap, data_gausfilter, plevel_var= 'u_r', ['u' ,'v'], 'RdBu_r', '', 925 #tangential wind
#var, imp_vars, cmap, data_gausfilter, plevel_var= 'v_r', ['u' ,'v'], 'RdBu_r', '', 500 #tangential wind

#var, imp_vars, cmap, plevel_var= 'u_a', ['u' ,'v', 'z'], 'RdBu_r',  850
#var, imp_vars, cmap, plevel_var= 'v_a', ['u' ,'v', 'z'], 'RdBu_r',  500


#var, imp_vars, cmap, plevel_var= 'U', ['u', 'v'], 'Reds', 850
#var, imp_vars, cmap, plevel_var= '10U', ['10u', '10v'], 'Reds', ''

#var, imp_vars, cmap, plevel_var= 'U_g', ['z'], 'Reds',  850
#var, imp_vars, cmap, plevel_var= 'U_r', ['u' ,'v'], 'Reds', 850 #it is the same as U which makes sense
#var, imp_vars, cmap, plevel_var= 'U_a', ['u' ,'v', 'z'], 'Reds',  500


#var, imp_vars, cmap, data_gausfilter, plevel_var= 'vo_g', ['z'], 'RdBu_r', '', 850


var_full= var
if ano_var==True:  var_full= var_full+ '_ano'
if type(plevel_var) == int: var_full= var_full+ '_'+str(plevel_var)
if type(plevel_var) == list:
    var_full += ('').join(['-'+str(p) for p in plevel_var])

#    if len(plevel_var)== 2: var_full += str(plevel_var[0])+'-'+str(plevel_var[1])
#    elif len(plevel_var)== 3: var_full += str(plevel_var[0])+'-'+str(plevel_var[1])+'-'+str(plevel_var[2])




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
elif var_full== 'theta_diff1000-500': levels= np.arange(-25, -10, 1)
elif var_full== 'theta_diff925-500': levels= np.arange(-23, -9, 1)
elif var_full== 'theta_e_diff925-500': levels= np.arange(-19, 0, 1)
elif var_full== 'theta_e_diff1000-500': levels= np.arange(-19, 0, 1)

elif var_full== 'sst-t_500': levels= np.arange(40, 46.1, 1)
elif var_full== 'sst-theta_500': levels= np.arange(-11, -4.9, 1)
elif var_full== 'theta_e_sst-theta_e_500': levels= np.arange(-1, 9.1, 1)
elif var_full== 'theta_e_2m-theta_e_500': levels= np.arange(-10, 5, 2)

elif var_full== 'theta_sst-theta_500': levels= np.arange(-11, -4.9, 1)
elif var_full== 'theta_sst-theta_700': levels= np.arange(-4, 2.1, 1)

elif 'adv_t_' in  var_full: levels= np.arange(-1, 1.01, .1)
elif 'adv_q_' in  var_full: levels= np.arange(-0.1, 0.101, .02)
elif 'grad_t' in  var_full: levels= np.arange(1.0, 3.1, .2)

elif var in ['u', 'v']: levels= np.arange(-20, 20.1, 2)
elif 'U_' in var_full: levels= np.arange(0, 20.1, 2)
elif var_full == '10U': levels= np.arange(6, 16.1, 1)

#elif var_full== 'N925-500': levels= np.arange(0.0042, .0067, 0.0002)

elif 'N' in var_full: levels= np.arange(0.0042, .0067, 0.0002)
#elif 'adv_t_' in  var_full: levels= np.arange(-15, 15, 1)*1E-5
elif var == 'tcwv': levels= np.arange(2, 6.1, 0.5)

elif 't_ano' in var_full:  levels = np.linspace(-7, 7, 15)
elif 'z_ano' in var_full:  levels = np.arange(-200, 201, 20)
elif var==  'vo':  levels = np.arange(-50, 50.1, 5)
elif var in ['adv_rel_vo', 'adv_plan_vo']:
    levels = np.arange(-10, 10.1, 1)
    levels= levels[levels!= 0]
    
elif var in ['stretch_vo', 'stretch_rel_vo', 'stretch_vo_self', 'tilt_vo', 'vert_adv_vo', 'adv_vo']:
    levels = np.arange(-10, 10.1, 1)
    levels= levels[levels!= 0]

elif 'pv_' in var_full:  levels = np.arange(0, 7.1, 1)*10E-8

elif 'cp' in var_full or 'lsp' in var_full: levels= np.arange(0, .51, .1) # levels= np.linspace(0, 1, 11)
elif 'tp' in var_full : levels= np.linspace(0, 1, 11)

elif 'cc' in var_full: levels= np.arange(0.2, 1.01, .1)
elif var in ['d', 'div']:  levels = np.array([-4, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 4])
elif 'w_' in var_full: levels_different, levels = True, np.array([-1.6, -.8, -.4, -.2, -.1, -.05, 0.05, .1, .2, .4, .8, 1.6])




scale= .7 #scale factor for the color bar. it takes scale* [min, max] of the all time step values, only relevant if levels = []
sym= False # only relevant if levels=[]
#cmap='jet'
#cmap= 'viridis'
#cmap= 'Greys_r'
#cmap= 'RdBu_r'

if ano_var: sym= True
if sym: cmap= 'RdBu_r'


if calculate == 'std':
    levels= np.arange(1, 5.1, 0.2)
    cmap= 'viridis_r'


"""import Stoll list"""
Stoll_imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


Stoll= pd.read_csv(Stoll_imp_dir + 'Stoll_list_noRojo.csv')
Stoll= Stoll.set_index(['ID', 'time'])

S_nodes= pd.read_csv(Stoll_imp_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')
S_nodes= S_nodes.set_index(['ID', 'time'])

Stoll= pd.concat([Stoll, S_nodes], axis=1)


Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps

Stoll= Stoll.reset_index()
Stoll['PLnr']= [int(ID.split('_')[0])*1E3 + int(ID.split('_')[1]) for ID in Stoll.ID.values]







"""evolution"""
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

#        if size != 500:
#            ds = ds.where(np.logical_and(np.logical_and(ds.x<=size, ds.x >=-size) , np.logical_and(ds.y<=size, ds.y >=-size) ), drop= True)

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


        if var=='theta_e_2m-theta_e':
            theta_e= EquiPotTemp(ds.t, ds.q/1E3, plevel_var)
#            theta_e_sst= EquiPotTemp(ds.sst, RH2SH(100, 1000, ds.sst)/1E3, 1000)
            theta_e_2m= EquiPotTemp(ds['2t'], RH2SH(100, ds.msl/100, ds['2t'])/1E3, ds.msl/100)

            ds[var_full]= (('time', 'x', 'y'), theta_e_2m - theta_e )
            ds[var_full].attrs['units'] = 'K'     
            ds[var_full].attrs['long_name']='theta_e_2m-theta_e'+str(plevel_var)



        if var=='theta_diff':
            theta_1, theta_0= PotTemp(ds.t.isel(plev=1), plevel_var[1]), PotTemp(ds.t.isel(plev=0), plevel_var[0])

            ds[var_full]= (('time', 'x', 'y'), theta_0 - theta_1)
            ds[var_full].attrs['units'] = 'K'

        if var=='theta_e_diff':
            theta_e_1, theta_e_0= EquiPotTemp(ds.t.isel(plev=1), ds.q.isel(plev=1)/1E3, plevel_var[1]), EquiPotTemp(ds.t.isel(plev=0), ds.q.isel(plev=0)/1E3, plevel_var[0])

            ds[var_full]= (('time', 'x', 'y'), theta_e_0- theta_e_1)
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

        if var== '10U':
            U= np.sqrt(ds['10u']**2 + ds['10v']**2)

            ds[var_full]= (('time', 'x', 'y'), U )
            ds[var_full].attrs['units'] = r"m $\cdot$ s$^{-1}$"
            ds[var_full].attrs['long_name']= 'Wind speed at 10m'

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

            grad_q_vec= np.gradient(ds.q, 25E3) #two components of the temperature gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km
            
            U= np.sqrt(ds.u**2 + ds.v**2)
            wind_beering= UV2Direction(ds.u, ds.v)
            track_beering = ds.beering.values
            
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='')
            
            ds[var_full]= (('time', 'x', 'y'), -1* (u_r * grad_q_vec[2] + v_r * grad_q_vec[1] ) * 3600 )
            ds[var_full].attrs['units'] = '(g/kg)/h'
            ds[var_full].attrs['long_name']='Horizontal humidity advection' 


        if var== 'adv_rel_vo':
            if data_gausfilter:
                ds['vo'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.vo, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
                ds['vo'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds.vo, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )

            grad_vo_vec= np.gradient(ds.vo, 25E3) #two components of the gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km
            
            U= np.sqrt(ds.u**2 + ds.v**2)
            wind_beering= UV2Direction(ds.u, ds.v)
            track_beering = ds.beering.values
            
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='')
            
            ds[var_full]= (('time', 'x', 'y'), -1* (u_r * grad_vo_vec[2] + v_r * grad_vo_vec[1] ) * 3600 )
            ds[var_full].attrs['units'] = '(10$^{-5}$s$^{-1}$ h$^{-1}$'
            ds[var_full].attrs['long_name']='Horizontal advection of rel. vorticity' 


        if var== 'adv_plan_vo':
            f= np.sin(np.deg2rad(ds.lat)) * 2*2*np.pi/(24*60**2) *1E5

            grad_f_vec= np.gradient(f, 25E3) #two components of the gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km
            
            U= np.sqrt(ds.u**2 + ds.v**2)
            wind_beering= UV2Direction(ds.u, ds.v)
            track_beering = ds.beering.values
            
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='')
            
            ds[var_full]= (('time', 'x', 'y'), -1* (u_r * grad_f_vec[2] + v_r * grad_f_vec[1] ) * 3600 )
            ds[var_full].attrs['units'] = '(10$^{-5}$s$^{-1}$ h$^{-1}$'
            ds[var_full].attrs['long_name']='Horizontal advection of plan. vorticity' 



        if var== 'stretch_vo':
            f= np.sin(np.deg2rad(ds.lat)) * 2*2*np.pi/(24*60**2) *1E5
            ds[var_full]= (('time', 'x', 'y'), -1* (ds.vo+ f) * ds.d * 3600 )

            ds[var_full].attrs['units'] = '(10$^{-5}$s$^{-1}$ h$^{-1}$'
            ds[var_full].attrs['long_name']='Vorticity stretching' 

        if var== 'stretch_rel_vo':
            ds[var_full]= (('time', 'x', 'y'), -1* ds.vo* ds.d * 3600 )

            ds[var_full].attrs['units'] = '(10$^{-5}$s$^{-1}$ h$^{-1}$'
            ds[var_full].attrs['long_name']='Relative vorticity stretching' 

        if var== 'tilt_vo':           
            U= np.sqrt(ds.u**2 + ds.v**2)
            wind_beering= UV2Direction(ds.u, ds.v)
            track_beering = ds.beering.values
            
            rot_beering= (track_beering[:, np.newaxis, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='')

            dudp= (u_r[:, 2,:,:]- u_r[:,0,:,:])/(100*(ds.plev[2]-ds.plev[0]))
            dvdp= (v_r[:, 2,:,:]- v_r[:,0,:,:])/(100*(ds.plev[2]-ds.plev[0]))

            grad_w_vec= np.gradient(ds.w.isel(plev=1), 25E3) #two components of the gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km
           
            ds[var_full]= (('time', 'x', 'y'), (-1* dudp* grad_w_vec[1] + dvdp* grad_w_vec[2] ) * 3600* 1E5 )
            ds[var_full].attrs['units'] = '(10$^{-5}$s$^{-1}$ h$^{-1}$'
            ds[var_full].attrs['long_name']='Vorticity tilting' 


        if var== 'vert_adv_vo':           
            DvoDp= (ds.vo.isel(plev=2)- ds.vo.isel(plev= 0))/(100*(ds.plev[2]-ds.plev[0]))

            ds[var_full]= (('time', 'x', 'y'), (-1* ds.w.isel(plev=1) * DvoDp) * 3600 )
            ds[var_full].attrs['units'] = '(10$^{-5}$s$^{-1}$ h$^{-1}$'
            ds[var_full].attrs['long_name']='Vertical vorticity advection' 


        if var== 'stretch_vo_self':
            U= np.sqrt(ds.u**2 + ds.v**2)
            wind_beering= UV2Direction(ds.u, ds.v)
            track_beering = ds.beering.values
            
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='')

            grad_u_vec= np.gradient(u_r, 25E3) #two components of the gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km
            grad_v_vec= np.gradient(v_r, 25E3) #two components of the gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km
            
            
            ds[var_full]= (('time', 'x', 'y'), -1* ds.vo* (grad_u_vec[2] + grad_v_vec[1] ) * 3600 )
            ds[var_full].attrs['units'] = '(10$^{-5}$s$^{-1}$ h$^{-1}$'
            ds[var_full].attrs['long_name']='Vorticity stretching' 


        if var== 'div': #the self calculated divergence to test if I get the same
            U= np.sqrt(ds.u**2 + ds.v**2)
            wind_beering= UV2Direction(ds.u, ds.v)
            track_beering = ds.beering.values
            
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='')

            grad_u_vec= np.gradient(u_r, 25E3) #two components of the gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km
            grad_v_vec= np.gradient(v_r, 25E3) #two components of the gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km

            
            ds[var_full]= (('time', 'x', 'y'), grad_u_vec[2] + grad_v_vec[1] )
            ds[var_full]*= 1E5
            ds[var_full].attrs['units'] = '(10$^{-5}$s$^{-1}$'
            ds[var_full].attrs['long_name']='Div (self calc)' 


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
        ds2= xr.open_dataset( Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/ci_allObs_track-smth-"+smooth_param+'.nc')
        ds2['ci']= ds2['ci'].fillna(value= 1)
           
        ds[var_full]= ds[var_full].where(ds2['ci'] < 0.3)
      
    
    
    
    
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
    

    if ano_var == True:
        ds[var_full]= (('time', 'x', 'y'), ds[var_full].values- np.nanmean(ds[var_full].values, axis= (1,2))[:, np.newaxis, np.newaxis])       



if sym: vextr= scale* np.max([np.max(ds[var_full]), -np.min(ds[var_full])])
else: vmax, vmin= scale*float(np.max(ds[var_full])), scale*float(np.min(ds[var_full]))


print('plot the SOMs')
#fig= plt.figure(fignr, figsize= (2.7*x,2.7*y) )
fig= plt.figure(fignr, figsize= (5,5.5) )

fignr+=1
plt.clf()





if calculate =='mean':
    ds_var_calc= np.nanmean(ds[var_full], axis= 0)
if calculate =='std':
    ds_var_calc= np.nanstd(ds[var_full], axis= 0)

if levels_different:
    norm = colors.BoundaryNorm(boundaries=levels, ncolors=256) #this is for the case when the levels are not equally distanced, e.g. for w where I also want to display small values
    cf= plt.contourf(ds.x, ds.y, ds_var_calc, levels= levels[1:-1], cmap= cmap, norm=norm, extend='both')      
    cs= plt.contour(ds.x, ds.y, ds_var_calc, levels= levels[1:-1], colors='k', linestyles='solid')

elif len(levels) > 0:
    if only_contours== False:
        cf= plt.contourf(ds.x, ds.y, ds_var_calc, levels= levels, cmap= cmap, extend='both')      
        cs= plt.contour(ds.x, ds.y, ds_var_calc, levels= levels, colors='k', linestyles='solid')

    if only_contours== True:
        cs= plt.contour(ds.x, ds.y, ds_var_calc, levels= levels, colors='k')

   
elif sym:
    cf= plt.contourf(ds.x, ds.y, ds_var_calc, vmin= -vextr, vmax= vextr, cmap= cmap)#, extend='both')      
    cs= plt.contour(ds.x, ds.y, ds_var_calc, colors='k', linestyles='solid')

else:
    cf= plt.contourf(ds.x, ds.y, ds_var_calc, cmap= cmap, vmin= vmin, vmax= vmax)
    cs= plt.contour(ds.x, ds.y, ds_var_calc, colors='k', vmin= vmin, vmax= vmax)

if contour_numbers:
    if len(levels) > 0: #make the label of the contours
        if np.max(abs(levels)) > 5: plt.clabel(cs, fontsize=10, fmt='%1.0f')#, inline=1)
        elif np.max(abs(levels)) > 0.5 and var != 'w': plt.clabel(cs, fontsize=10, fmt='%1.1f')#, inline=1)
        elif np.max(abs(levels)) > 0.05 or var== 'w': plt.clabel(cs, fontsize=10, fmt='%1.2f')#, inline=1)
        
        elif np.max(abs(levels)) > 0.005 : plt.clabel(cs, fontsize=10, fmt='%1.4f')#, inline=1)

    else: plt.clabel(cs, fontsize=10, fmt='%1.1f')


plt.scatter(0,0, color='r', zorder= 2)

plt.xlim(-size, size)
plt.ylim(-size, size)

plt.yticks(np.arange(-size,size+1, size/2))
plt.xticks(np.arange(-size,size+1, size/2))
#    ax1.set_xticklabels('')
#    ax1.set_yticklabels('')    

plt.tight_layout()

if only_contours== False:
    fig.subplots_adjust(bottom=0.2)
    cbar_ax = fig.add_axes([0.09, 0.12, 0.6, 0.02])
    cb= fig.colorbar(cf, cax=cbar_ax, orientation="horizontal")

    """make the arrow"""
    arrow_ax = fig.add_axes([0.75, 0.05, 0.2, 0.1])
    import matplotlib.patches as mpatches
    
    arrow = mpatches.FancyArrowPatch((0.2, 0.7), (.7, 0.7), mutation_scale=20)
    arrow_ax.add_patch(arrow)
    arrow_ax.text(0,0 , 'Propagation\n direction', fontweight= 'bold')#, fontsize=13)
    arrow_ax.set_frame_on(False)
    arrow_ax.get_xaxis().set_visible(False)
    arrow_ax.get_yaxis().set_visible(False)



if 'z_ano' in var_full: labelvar= 'Geopotential height anomaly [m]'
elif var == 'z': labelvar= 'Geopotential height [m]'
elif var == 'z_pv': labelvar= 'Tropopause height [m]'

elif var == 'vo': labelvar= 'Relative vorticity [10$^{-5}$ s$^{-1}$]'

elif 't_diff' in var_full: labelvar= 'Temperature difference '+str(plevel_var[0])+'-'+str(plevel_var[1])+'hPa [K]'
elif 'theta_diff' in var_full: labelvar= 'Pot. temp. diff. '+str(plevel_var[0])+'-'+str(plevel_var[1])+'hPa [K]'
elif 'theta_e_diff' in var_full: labelvar= 'Equiv. pot. temp. diff. '+str(plevel_var[0])+'-'+str(plevel_var[1])+'hPa [K]'
elif 'sst-t_' in var_full: labelvar= 'SST - T'+str(plevel_var)+' [K]'

    
elif 'N' in var_full: labelvar= 'N '+str(plevel_var[0])+'-'+str(plevel_var[1])+'hPa [1/s]'
elif 'grad_t' in var_full: labelvar= r"$\nabla_h$ T "+str(plevel_var)+' [K/100km]'
elif 'adv_t' in var_full: labelvar= "Horizontal temperature advection "+str(plevel_var)+' [K/h]'
elif 'adv_q' in var_full: labelvar= "Horizontal humidity advection "+str(plevel_var)+' [(g/kg)/h]'

elif 't_ano' in var_full: labelvar= 'Temperature anomaly [K]'
elif var =='t': labelvar= 'Temperature [K]'
elif 'U_ano' in var_full: labelvar= 'Wind speed anomaly [m/s]'
elif var == 'U' : labelvar= 'Wind speed [m/s]'
elif 'q_ano' in var_full: labelvar= 'Specific humidity anomaly [g(kg)]'
elif var == 'q': labelvar= 'Specific humidity [g/kg]'
elif 'pres_pv' in var_full: labelvar= 'Tropopause level [hPa]'


elif ano_var== False: labelvar= ds[var_full].long_name + ' ['+ ds[var_full].units+ ']'
else: labelvar= ''

if only_contours== False:
    if calculate == 'std':
        labelvar= 'Standard deviation '+ labelvar
    cb.set_label(labelvar, size=12) 

if save:
    save_name=savedir+ calculate+'_'+var_full+ '_track-smth'+smooth_param +"_"+Obs
    
    if vo_tend_thresh: save_name+= '_vo-tend-'+str(vo_tend_thresh)
    
    save_name+='_size'+str(size)

    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')



