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
import scipy.ndimage.filters as filters

from f_useful import *
from f_STARS import *
from f_meteo import *

plt.rcParams.update({'font.size': 13})


#
save= True
#save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM3/'
#
fignr= 8

imp_list=True
#imp_list=Falsez 

imp_ds=True
#imp_ds=False


"""SOM var"""
Obs= 'allObs_SOMprep' #the other versions were deleted

x=3 #3
y=3 #x+1

mirror_top_bottom=False
rotate_clock=False


SOM=1


if x== 3 and y == 4: mirror_top_bottom=True
if x== 3 and y == 3: rotate_clock=True

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
#var='2t'
#plevel_var= '' #925
var='t'
plevel_var= '700'


cmap= 'RdBu_r'


var_full= var
if ano_var==True:  var_full += '_ano'
if type(plevel_var) == int: var_full += '_'+str(plevel_var)

if 't_ano' in var_full:  levels = np.linspace(-8, 8, 17)
elif 'pres_pv' in var_full: levels= np.arange(300, 421, 20)
elif 'tp' in var_full: levels= np.linspace(0, 1, 11)
elif 'cc' in var_full: levels= np.arange(0.2, 1.01, .1)
elif 'd' in var_full:  levels = [-4, -3, -2, -1, 1, 2, 3, 4]


"""contour variable"""
ano_cont=True
cont='z'
plevel_cont= plevel_var

#cont='msl'

cont_full= cont
if ano_cont==True:  cont_full += '_ano'
if type(plevel_cont) == int: cont_full += '_'+str(plevel_cont)

if 'z_ano' in cont_full: cont_levels = np.arange(-300, 301, 20)  
if cont_full== 'msl': cont_levels = np.arange(960, 1030, 1) 
    
var_list= [var, cont]
var_full_list= [var_full, cont_full]




"arrows"
plot_arrow=False
#arrow='w'
#plevel_arrow=700
#
#
##cont_= cont
##if type(plevel_cont) == int: cont_full += '_'+str(plevel_cont)
#
#var_list= [var, cont, arrow]
#var_full_list= [var_full, cont_full]


"""contour variable 2"""
plot_cont2=True
ano_cont2=False

compute_cont2=False
#cont2='mcc'
#cont2= 'sshf'
#plevel_cont2='' #925

cont2='w'
plevel_cont2=700 #925


#compute_cont2=True
#cont2, imp_cont2, data_gausfilter, plevel_cont2= 'adv_t', ['t', 'u', 'v'], 4, 300



if plot_cont2:
    cont2_full= cont2
    if ano_cont2==True:  cont2_full += '_ano'
    if type(plevel_cont2) == int: cont2_full += '_'+str(plevel_cont2)
    
    if cont2 == 'd': cont_levels2 = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]  
    if cont2 == 'sshf': cont_levels2 = np.arange(80, 200.1, 20)

    if cont2 == 'mcc': cont_levels2= [0.6, 0.7, 0.8] #here warnings are displayed since somethimes no contours are found
    if cont2 == 'adv_t': cont_levels2= np.arange(-2, 2.01, .2)
    if cont2 == 'w': cont_levels2 = np.array([-1.6, -.8, -.4, -.2, -.1, -.05, 0.05, .1, .2, .4, .8, 1.6])

    var_list+= [cont2]
    var_full_list += [cont2_full]


"""contour variable 3"""
plot_cont3=True

compute_cont3=False
cont3='mcc'
#cont3= 'sshf'
plevel_cont3='' #925


if plot_cont3:
    cont3_full= cont3
    if type(plevel_cont3) == int: cont3_full += '_'+str(plevel_cont3)
    
    if cont3 == 'mcc': cont_levels3= [0.7, 0.8, 0.9] #here warnings are displayed since somethimes no contours are found

    var_list+= [cont3]
    var_full_list += [cont3_full]



"""contour variable 3"""
plot_cont4=True

compute_cont4=True
cont4, imp_cont4 = 'tp', ['cp', 'lsp']
#cont3= 'sshf'
plevel_cont4='' #925


if plot_cont4:
    cont4_full= cont4
    if type(plevel_cont4) == int: cont4_full += '_'+str(plevel_cont4)
    
    if cont4 == 'mcc': cont_levels4= [0.7, 0.8, 0.9] #here warnings are displayed since somethimes no contours are found
    if cont4 == 'tp': cont_levels4= np.arange(.2, 2, .2)

    var_list+= [cont4]
    var_full_list += [cont4_full]

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

    SOM_filedir= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/SOM/"+Svar_full+"_"+str(Splevel)+'_'+Obs+ '_track-smth-'+smooth_param+'_dur'+str(lifelim)+"_x"+ str(x)+"_y"+str(y)

    ds_init= xr.open_dataset(SOM_filedir+"_cluster.nc")

    Stoll, Stoll_ind= imp_standard_Stoll() #SOM_filedir)#, Obs_evol)

    df_nodes= pd.read_csv(SOM_filedir+"_node_nr.txt", sep=" ") 
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


    Stoll= pd.concat([Stoll, df_nodes], axis= 1)
    Stoll= Stoll_Vort_tend(Stoll)




if imp_ds:
    print('import from all time step list')

    """file for shading variable"""
    if plevel_var != '': plev_str= '_all_levs'
    else: plev_str= ''
    
    file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +var + plev_str + '_allObs_track-smth-'+smooth_param+'.nc'
    ds= xr.open_dataset(file)
    if plevel_var != '': ds= ds.sel(plev= plevel_var)
    
    ds= ds.rename({var: var_full})
    

    """file for cont variable"""
    if plevel_cont != '': plev_str= '_all_levs'
    else: plev_str= ''
    
    file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +cont + plev_str + '_allObs_track-smth-'+smooth_param+'.nc'
    ds2= xr.open_dataset(file)
    if plevel_var != '': ds2= ds2.sel(plev= plevel_cont)
    
    ds2= ds2.rename({cont: cont_full})
    ds= xr.merge([ds, ds2])
    


    """file for arrow"""
    if plot_arrow:
        
        if plevel_arrow != '': plev_str= '_all_levs'
        else: plev_str= ''
        
        file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +arrow + plev_str + '_allObs_track-smth-'+smooth_param+'.nc'
        ds2= xr.open_dataset(file)
        if plevel_var != '':
            ds2= ds2.sel(plev= plevel_arrow)
            ds2= ds2.rename({'plev': 'plev2'})

        
#        ds2= ds2.rename({cont2: cont2_full})
        ds= xr.merge([ds, ds2])


    
    """file for cont variable 2"""
    if plot_cont2:
        
        if compute_cont2== False:
        
            if plevel_cont2 != '': plev_str= '_all_levs'
            else: plev_str= ''
            
            file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +cont2 + plev_str + '_allObs_track-smth-'+smooth_param+'.nc'
            ds2= xr.open_dataset(file)
            if plevel_cont2 != '':
                ds2= ds2.sel(plev= plevel_cont2)
                ds2= ds2.rename({'plev': 'plev2'})
    
            
            ds2= ds2.rename({cont2: cont2_full})
            ds= xr.merge([ds, ds2])
    

        elif compute_cont2:
            for ni, i_var in enumerate(imp_cont2):
                file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +i_var +'_allObs_track-smth-'+smooth_param+'.nc'
                if os.path.isfile(file):
                    ds1= xr.open_dataset(file)
                else:
                    ds1= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +i_var + '_all_levs_allObs_track-smth-'+smooth_param+'.nc')
                    ds1= ds1.sel(plev= plevel_cont2)
                    
                if ni == 0: ds2= ds1
                else: ds2= xr.merge([ds2, ds1]) 
    
    
            if cont2== 'adv_t':
                if data_gausfilter:
                    ds2['t'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.t, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
                    ds2['t'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.t, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )
    
                grad_t_vec= np.gradient(ds2.t, 25E3) #two components of the temperature gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km
    
                U= np.sqrt(ds2.u**2 + ds2.v**2)
                wind_beering= UV2Direction(ds2.u, ds2.v)
                track_beering = ds2.beering.values
                
                rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
                v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='')
                
                ds[cont2_full]= (('time', 'x', 'y'), -1* (u_r * grad_t_vec[2] + v_r * grad_t_vec[1] ) * 3600 )
                ds[cont2_full].attrs['units'] = 'K/h'
                ds[cont2_full].attrs['long_name']='Horizontal temperature advection' 



    """file for cont variable 3"""
    if plot_cont3:
        
        if compute_cont3== False:
        
            if plevel_cont3 != '': plev_str= '_all_levs'
            else: plev_str= ''
            
            file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +cont3 + plev_str + '_allObs_track-smth-'+smooth_param+'.nc'
            ds2= xr.open_dataset(file)
            if plevel_cont3 != '':
                ds2= ds2.sel(plev= plevel_cont3)
                ds2= ds2.rename({'plev': 'plev3'})
    
            
            ds2= ds2.rename({cont3: cont3_full})
            ds= xr.merge([ds, ds2])



    """file for cont variable 4"""
    if plot_cont4:
        
        if compute_cont4== False:
        
            if plevel_cont4 != '': plev_str= '_all_levs'
            else: plev_str= ''
            
            file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +cont4 + plev_str + '_allObs_track-smth-'+smooth_param+'.nc'
            ds2= xr.open_dataset(file)
            if plevel_cont4 != '':
                ds2= ds2.sel(plev= plevel_cont4)
                ds2= ds2.rename({'plev': 'plev2'})
    
            
            ds2= ds2.rename({cont4: cont4_full})
            ds= xr.merge([ds, ds2])
    

        elif compute_cont4:
            for ni, i_var in enumerate(imp_cont4):
                file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +i_var +'_allObs_track-smth-'+smooth_param+'.nc'
                if os.path.isfile(file):
                    ds1= xr.open_dataset(file)
                else:
                    ds1= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +i_var + '_all_levs_allObs_track-smth-'+smooth_param+'.nc')
                    ds1= ds1.sel(plev= plevel_cont4)
                    
                if ni == 0: ds2= ds1
                else: ds2= xr.merge([ds2, ds1]) 
    
    
            if cont4== 'tp':                  
                ds[cont4_full]= (('time', 'x', 'y'), (ds2['cp']+ ds2['lsp'] ) * 1E3 )
                ds[cont4_full].attrs['units'] = 'mm'
                ds[cont4_full].attrs['long_name']='Total precipitation' 

    
    
#    ds2= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields/" +cont3 + '_allObs.nc')
#    ds= xr.merge([ds, ds2])



    
    """add the node to ds """
    ds_node_var= np.zeros(len(ds.time))
    
    for s in range(len(df_nodes)):
        df_nodes_s= df_nodes.iloc[s]
        ds_index= np.argwhere(np.logical_and(ds.PLnr.values== df_nodes_s.PLnr, ds.time.values== np.datetime64(df_nodes_s.name[1])))[0][0]
        ds_node_var[ds_index]= df_nodes_s.node
    
    ds['node']= (('time'), ds_node_var)
    
    
    """select only the specific SOM node"""
    ds= ds.where(ds.node == SOM, drop=True) #
    
 
    
    """calculate with data"""
    if ano_var == True:
        ds[var_full]= (('time', 'x', 'y'), ds[var_full].values- np.mean(ds[var_full].values, axis= (1,2))[:, np.newaxis, np.newaxis])
    if ano_cont == True:
        ds[cont_full]= (('time', 'x', 'y'), ds[cont_full].values- np.mean(ds[cont_full].values, axis= (1,2))[:, np.newaxis, np.newaxis])
#    if ano_cont2 == True:
#        ds[cont2_p_full]= (('time', 'x', 'y'), ds[cont2_p].values- np.mean(ds[cont2_p].values, axis= (1,2))[:, np.newaxis, np.newaxis])
#    if ano_cont3 == True:
#        ds[cont3_full]= (('time', 'x', 'y'), ds[cont3_full].values- np.mean(ds[cont3_full].values, axis= (1,2))[:, np.newaxis, np.newaxis])
#        
        
        
#    if 'tp' in list(ds.keys()):
#        ds['tp']*= 1E3
#        ds['tp'].attrs['units']= 'mm'

    if 'msl' in list(ds.keys()):
        ds['msl']/= 1E2
        ds['msl'].attrs['units']= 'hPa'
        
    if 'slhf' in list(ds.keys()):
        ds['slhf']/= -60**2
        ds['slhf'].attrs['units'] = 'W/m**2'
    if 'sshf' in list(ds.keys()):
        ds['sshf']/= -60**2
        ds['sshf'].attrs['units'] = 'W/m**2'         
        
    if 'd' in var_list:
        ind= np.where(np.array(var_list) == 'd')[0][0]
        ds[var_full_list[ind]]*= 1E5
        ds[var_full_list[ind]].attrs['units'] = '1e-5 '+ ds[var_full_list[ind]].attrs['units']



"""plot"""
fig= plt.figure(fignr, figsize= (7,7) )

fignr+=1
plt.clf()
    

cf= plt.contourf(ds.x, ds.y, np.mean(ds[var_full], axis= 0), levels= levels, cmap= cmap, extend='both')      

cs= plt.contour(ds.x, ds.y, np.nanmean(ds[cont_full], axis= 0), levels= cont_levels, colors='k')#, linestyles='solid')

if plot_cont2:
    cs2= plt.contour(ds.x, ds.y, np.mean(ds[cont2_full], axis= 0), levels= cont_levels2, colors= 'r')
  

if plot_cont3:
    cs3= plt.contour(ds.x, ds.y, np.mean(ds[cont3_full], axis= 0), levels= cont_levels3, colors= 'grey')

if plot_cont4:
    cs4= plt.contour(ds.x, ds.y, np.mean(ds[cont4_full], axis= 0), levels= cont_levels4, colors= 'b')


#        cont_levels3= np.arange(0, 76, 5) #[1000,  1200, 1400]
#        cf= plt.contourf(ds.x, ds.y, np.average(ds[cont3_full], axis= 0), levels= cont_levels3, cmap= cmap, extend='both')
#        cs3= plt.contour(ds.x, ds.y, np.average(ds[cont3_full], axis= 0), levels= cont_levels3, colors= 'g')



if np.max(cont_levels) > 3: plt.clabel(cs, cont_levels[::2], fmt='%1.0f')#, inline=1)
else: plt.clabel(cs, cont_levels[::2],  fmt='%1.1f')#, inline=1)

if plot_cont2:
    if np.max(np.mean(ds[cont2_full], axis= 0) ) > min(cont_levels2):
        if max(cont_levels2) > 3: plt.clabel(cs2, cont_levels2[::2], fmt='%1.0f')
        else: plt.clabel(cs2, cont_levels2[::2], fmt='%1.1f')

if plot_cont3:
    if np.max(np.mean(ds[cont3_full], axis= 0) ) > min(cont_levels3):
        if max(cont_levels3) > 3: plt.clabel(cs3, cont_levels3[::2], fmt='%1.0f')
        else: plt.clabel(cs3, cont_levels3[::2], fmt='%1.1f')
        
if plot_cont4:
    if np.max(np.mean(ds[cont4_full], axis= 0) ) > min(cont_levels3):
        if max(cont_levels4) > 3: plt.clabel(cs4, cont_levels4[::2], fmt='%1.0f')
        else: plt.clabel(cs4, cont_levels4[::2], fmt='%1.1f')


plt.scatter(0,0, color='r')
    

plt.tight_layout()

fig.subplots_adjust(bottom=0.14)
cbar_ax = fig.add_axes([0.12, 0.08, 0.5, 0.015])
cb= fig.colorbar(cf, cax=cbar_ax, orientation="horizontal")





"""make the arrow"""
arrow_ax = fig.add_axes([0.75, 0.03, 0.2, 0.08])
import matplotlib.patches as mpatches

arrow = mpatches.FancyArrowPatch((0.2, 0.8), (.7, 0.8), mutation_scale=20)
arrow_ax.add_patch(arrow)
arrow_ax.text(0,0 , 'Propagation\n direction')#, fontweight= 'bold')#, fontsize=13)
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
cb.set_label(labelvar) 


if save:
    save_name=savedir+ 'SOM_track-smth'+smooth_param+'_'+Svar_full+"_"+str(Splevel)+"_"+Obs+'_SOM_'+ str(SOM)+'--mean_'+var_full+"_"+cont_full


    if plot_cont2: save_name+= "_"+cont2_full 
    if plot_cont3: save_name+= "_"+cont3_full 
    if plot_cont4: save_name+= "_"+cont4_full 
    
    save_name += '_'+Obs_evol_str+"_x"+str(x)+"_y"+str(y)

    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')



