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
import matplotlib.patches as mpatches

import numpy as np
from f_useful import *
from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs

#
save= True
#save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM4/'
#
fignr= 11


load_ds=True


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
vo_tend_thresh= 1E-5


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






#selection= 'surf'
#selection, s_plevel, s_cont2_lev= 't_adv', 500, 925
#selection, s_plevel, s_cont2_lev= 't_adv', 925, 500
#selection, s_plevel, s_cont2_lev= 'mid-level', 700, 500
#selection, s_plevel, s_cont2_lev= 'mid-level', 850, 500
selection= 't+tilt'
#selection= 'moist'
#selection= 'flux'
#selection= 'flux_Bowen'

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
var='t'
plevel_var=500
compute_var, imp_var, lsm_var= False, None, False
cmap= 'RdBu_r'

if selection== 'surf': ano_var, var, plevel_var = True, '2t', ''
if selection== 't_adv': ano_var, var, plevel_var = True, 't', s_plevel
if selection=='mid-level': ano_var, var, plevel_var = False, 'w', s_plevel
if selection== 't+tilt': ano_var, var, plevel_var = True, 't', 850
if selection== 'moist': ano_var, var, plevel_var, cmap = False, 'tcwv', '', 'YlGnBu'
#if selection== 'flux': ano_var, var, plevel_var, compute_var, imp_var = False, '10U', '', True, ['10u', '10v']
if selection=='flux': ano_var, var, plevel_var, lsm_var = False, 'sshf', '', True
if selection=='flux_Bowen': ano_var, var, plevel_var, lsm_var, compute_var, imp_var = False, 'flux', '', True, True, ['sshf', 'slhf']


var_full= var
if ano_var==True:  var_full += '_ano'
if type(plevel_var) == int: var_full += '_'+str(plevel_var)


levels_different = False
if 't_ano' in var_full:  levels = np.linspace(-6, 6, 13) #np.linspace(-8, 8, 17)
if var =='2t':  levels = np.linspace(-8, 8, 17)
elif var== 'w' : levels_different, levels = True, np.array([-1.6, -.8, -.4, -.2, -.1, -.05, 0.05, .1, .2, .4, .8, 1.6])
elif var_full== 'q_850': levels= np.arange(0.4, 2, .1)
elif var_full== '10U': levels, cmap= np.arange(4, 18.1), 'Reds'
elif var_full in ['sshf', 'slhf']: levels, cmap= np.arange( 0, 201, 20), 'Reds'
elif var_full == 'flux': levels, cmap= np.arange( 50, 401, 50), 'Reds'

elif var == 'tcwv': levels= np.arange(2, 6.1, 0.5)

elif 'pres_pv' in var_full: levels= np.arange(300, 421, 20)
elif 'tp' in var_full: levels= np.linspace(0, 1, 11)
elif 'cc' in var_full: levels= np.arange(0.2, 1.01, .1)
elif 'd' in var_full:  levels = [-4, -3, -2, -1, 1, 2, 3, 4]


"""contour variable"""
ano_cont=True
cont='z' 
cont_color='k'
cont_label= True
plevel_cont=plevel_var
compute_cont, imp_cont, lsm_cont= False, None, False

if selection in ['surf']: ano_cont, cont = False, 'msl'
if selection== 't_adv': ano_cont, cont = True, 'z'
if selection=='mid-level': ano_cont, cont = True, 'z'
if selection=='t+tilt': ano_cont, cont, plevel_cont = False, 'msl', ''
if selection=='moist': ano_cont, cont, plevel_cont, compute_cont, imp_cont, cont_color, cont_label = False, 'tp', '', True, ['cp', 'lsp'], 'w', False
if selection=='flux': ano_cont, cont, plevel_cont, lsm_cont, cont_color = False, 'slhf', '', True, 'w'
if selection=='flux_Bowen': cont, ano_cont= 'Bowen' , False #ano_cont, cont, plevel_cont, lsm_cont = False, 'slhf', '', True

#if selection=='flux': ano_cont, cont, plevel_cont, compute_cont, imp_cont, lsm_cont = False, 'Bowen', '', True, ['sshf', 'slhf'], True


cont_full= cont
if ano_cont==True:  cont_full += '_ano'
if type(plevel_cont) == int: cont_full += '_'+str(plevel_cont)

if cont_full== 'msl': cont_levels = np.arange(960, 1030, 2) 
elif 'z_ano' in cont_full: cont_levels = np.arange(-200, 201, 20)  
elif 'tp' in cont_full: cont_levels= np.arange(0, 1.1, .2)
elif cont_full in ['sshf', 'slhf']: cont_levels= np.arange( 0, 500, 20)
#elif cont_full== 'Bowen': cont_levels= np.arange( 0, 2, 0.1)


"""contour variable 2"""
color_cont2= 'g'
plot_cont2=True
ano_cont2=False
imp_cont2, compute_cont2, lsm_cont2, data_gausfilter2= None, False, False, ''
#cont2='d'
cont2='w'
plevel_cont2=plevel_var

#if selection== 'surf': cont2, imp_cont2, plevel_cont2, plot_cont2, ano_cont2, compute_cont2 = '10U', ['10u', '10v'], '', True, False, True
if selection== 'surf': cont2, imp_cont2, plevel_cont2, compute_cont2, ano_cont2, plot_cont2, data_gausfilter2= 'adv_t', ['t', 'u', 'v'], 500, True, False, True, 4
if selection== 't_adv': cont2, imp_cont2, plevel_cont2, compute_cont2, ano_cont2, plot_cont2, data_gausfilter2= 'adv_t', ['t', 'u', 'v'], s_cont2_lev, True, False, True, 4
if selection=='mid-level': cont2, imp_cont2, plevel_cont2, compute_cont2, ano_cont2, plot_cont2, data_gausfilter2, lsm_mask2 = 'theta_e_sst-theta_e', ['t', 'q', 'sst', 'msl'], s_cont2_lev, True, False, True, False, True
#if selection=='mid-level': cont2, imp_cont2, plevel_cont2, compute_cont2, ano_cont2, plot_cont2, data_gausfilter2, lsm_mask2 = 'theta_e_2m-theta_e', ['t', 'q', '2t', 'msl'], s_cont2_lev, True, False, True, False, True
#if selection=='mid-level': cont2, imp_cont2, plevel_cont2, compute_cont2, ano_cont2, plot_cont2, data_gausfilter2, lsm_mask2 = 'theta_e_diff', ['t', 'q'], [925, s_cont2_lev], True, False, True, False, True
if selection== 't+tilt': cont2, plevel_cont2, ano_cont2, plot_cont2= 'z', 500,  False, True
if selection== 'moist': cont2, plevel_cont2, ano_cont2, plot_cont2= 'z', 500,  False, False
if selection in ['flux', 'flux_Bowen']: plot_cont2= False #ano_cont2, cont2, plevel_cont2, plot_cont2, color_cont2 = False, 'slhf', '', True, 'b'


#ano_cont2=False
#cont2='cape'
#plevel_cont2=''


#var, imp_vars, plevel_var= 'theta_e_diff', ['t', 'q'], [925, 500]


if plot_cont2:
    cont2_full= cont2
    if ano_cont2==True:  cont2_full += '_ano'
    if type(plevel_cont2) == int: cont2_full += '_'+str(plevel_cont2)
    if type(plevel_cont2) == list: cont2_full += ('').join(['-'+str(p) for p in plevel_cont2])
    
    if cont2 == 'd': cont_levels2 = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]  
    if cont2 == 'cape': cont_levels2= [30, 40, 50, 60] #here warnings are displayed since somethimes no contours are found
    if cont2 == 'adv_t': cont_levels2= np.hstack(( np.arange(-2, 0, .2), np.arange(0.2, 2.01, .2) ))
    if cont2_full in ['adv_t_500', 'adv_t_300'] : cont_levels2= np.hstack(( np.arange(-2, 0, .1), np.arange(0.1, 2.01, .1) ))
    if cont2_full== 'theta_e_sst-theta_e_500': cont_levels2= np.arange(2, 20, 2)
    if cont2_full== 'theta_e_2m-theta_e_500': cont_levels2= np.arange(-14, 4, 2)
    elif cont2_full== 'theta_e_diff925-500': cont_levels2= np.arange(2, 19.1, 2)
    elif 'w_' in cont2_full: cont_levels2 = np.array([-1.6, -.8, -.4, -.2, -.1, -.05, 0.05, .1, .2, .4, .8, 1.6])
    elif cont2_full== 'z_500': cont_levels2= np.arange(4500, 5500, 20)
    elif cont2_full == 'slhf': cont_levels2= np.arange( 0, 500, 20)





"""contour variable 3"""
plot_cont3=True
ano_cont3=False
compute_cont3=True
#cont2='d'
cont3='w'
plevel_cont3=plevel_var

#if selection== 'surf': cont3, imp_cont3, plevel_cont3, plot_cont3, compute_cont3 = 'tp', ['cp', 'lsp'], '', True, True
if selection in ['surf', 't_adv', 'mid-level', 't+tilt', 'moist', 'flux', 'flux_Bowen']: plot_cont3=False


if plot_cont3:
    cont3_full= cont3
    if type(plevel_cont3) == int: cont3_full += '_'+str(plevel_cont3)
    
    if cont3 == 'mcc': cont_levels3= [0.7, 0.8, 0.9] #here warnings are displayed since somethimes no contours are found
    if cont3 == 'tp': cont_levels3= np.arange(.2, 2, .2)




"""wind_variable"""
plot_wind=False


if selection in ['flux', 'flux_wind']:
    plot_wind = True
    arrow_type= 'barb'
    uwind, vwind = '10u', '10v'
    plevel_wind= ''
       
    wind= 'wind_r'
    compute_wind= True
    imp_wind=[uwind, vwind]
    wind_full=[uwind+'_r', vwind+'_r']

if selection== 'moist':
    plot_wind = True
    arrow_type= 'quiver'

    uwind, vwind = 'viwve', 'viwvn'
    plevel_wind= ''
       
    wind= 'wind_r'
    compute_wind= True
    imp_wind=[uwind, vwind]
    wind_full=[uwind+'_r', vwind+'_r']
    


if selection== 't+tilt': insert_shear_vec= True
else: insert_shear_vec= False
 
    
"""import Stoll list"""
Stoll_imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


Stoll= pd.read_csv(Stoll_imp_dir + 'Stoll_list_noRojo.csv')
Stoll= Stoll.set_index(['ID', 'time'])

S_nodes= pd.read_csv(Stoll_imp_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')
S_nodes= S_nodes.set_index(['ID', 'time'])

S_ERA5= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5.csv')
S_ERA5= S_ERA5.set_index(['ID', 'time'])

S_ERA5_2= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5_2.csv')
S_ERA5_2= S_ERA5_2.set_index(['ID', 'time'])

Stoll= pd.concat([Stoll, S_nodes, S_ERA5, S_ERA5_2], axis=1)


Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps

Stoll= Stoll.reset_index()
Stoll['PLnr']= [int(ID.split('_')[0])*1E3 + int(ID.split('_')[1]) for ID in Stoll.ID.values]







"""evolution"""
if load_ds:
    print('import from all time step list')

    """file for shading variable"""   
    ds= imp_ds(var, plevel_var, var_full, compute= compute_var, imp_var= imp_var, lsm_mask= lsm_var)

    """file for cont variable"""
    if cont != 'Bowen':
        ds2= imp_ds(cont, plevel_cont, cont_full, compute= compute_cont, imp_var= imp_cont, lsm_mask= lsm_cont)
        if plevel_cont != '':
            if plevel_cont != plevel_var: ds2= ds2.rename({'plev': 'plev1'})
        
        ds= xr.merge([ds, ds2])
    

    
    """file for cont variable 2"""
    if plot_cont2:
        ds2= imp_ds(cont2, plevel_cont2, cont2_full,
                    compute= compute_cont2, imp_var= imp_cont2, data_gausfilter= data_gausfilter2, lsm_mask= lsm_cont2)

#        if plevel_cont2 != plevel_var or plevel_cont2 != plevel_cont: ds2= ds2.rename({'plev': 'plev2'})
        
        ds= xr.merge([ds, ds2])



    if plot_cont3:
        ds2= imp_ds(cont3, plevel_cont3, cont3_full,
                    compute= compute_cont3, imp_var= imp_cont3)#, data_gausfilter= data_gausfilter3)

        ds= xr.merge([ds, ds2])



    if plot_wind:
        ds2= imp_ds(wind, plevel_wind, wind_full, #var_full is irrelevant
                    compute= compute_wind, imp_var= imp_wind)#, data_gausfilter= data_gausfilter3)

        ds= xr.merge([ds, ds2])




    
    """add the node to ds """
    ds_node= np.zeros(len(ds.time))
    ds_vo_tendency= np.zeros(len(ds.time))
    ds_grad_t= np.zeros(len(ds.time))

    for s in range(len(Stoll)):
        S= Stoll.iloc[s]
        ds_index= np.argwhere(np.logical_and(ds.PLnr.values== S.PLnr, ds.time.values== np.datetime64(S.time)))[0][0]
        ds_node[ds_index]= S.node
        ds_vo_tendency[ds_index]= S.vo_tendency
        ds_grad_t[ds_index]= S['grad_t850_dgauss-10_mean-1']

    
    ds['node']= (('time'), ds_node)
    ds['vo_tendency']= (('time'), ds_vo_tendency)
    ds['grad_t850_dgauss-10']= (('time'), ds_grad_t)
    
    

    
    
    if ano_var == True:
        ds[var_full]= (('time', 'x', 'y'), ds[var_full].values- np.mean(ds[var_full].values, axis= (1,2))[:, np.newaxis, np.newaxis])
    if ano_cont == True:
        ds[cont_full]= (('time', 'x', 'y'), ds[cont_full].values- np.mean(ds[cont_full].values, axis= (1,2))[:, np.newaxis, np.newaxis])
#    if ano_cont2 == True:
#        ds[cont2_p_full]= (('time', 'x', 'y'), ds[cont2_p].values- np.mean(ds[cont2_p].values, axis= (1,2))[:, np.newaxis, np.newaxis])
#    if ano_cont3 == True:
#        ds[cont3_full]= (('time', 'x', 'y'), ds[cont3_full].values- np.mean(ds[cont3_full].values, axis= (1,2))[:, np.newaxis, np.newaxis])
#        
        

    if var_full=='w_700':
        ds['w_700'] *= 100
        levels *= 100
        

print('plot the SOMs')
#fig= plt.figure(fignr, figsize= (2.7*x,2.7*y) )
fig= plt.figure(fignr, figsize= (2.5*x,2.5*y +1) )

fignr+=1
plt.clf()

#grad_t_thresh= .4

for iSOM in range(1, x*y +1):
    
    print('SOM ', iSOM)
    
    fig= plt.figure(fignr-1)
    ax1= plt.subplot(y, x, iSOM)
    
    
    ds_node= ds.where(ds.node == iSOM, drop=True)
    
    if vo_tend_excl:
        before_tend_excl= len(ds_node.time)
        ds_node= ds_node.where(ds_node.vo_tendency > vo_tend_thresh, drop=True)

#        ds_node= ds_node.where(ds_node['grad_t850_dgauss-10'] < grad_t_thresh, drop=True)
        after_tend_excl= len(ds_node.time)
    
    
    if levels_different: #for w the color levels are uneven, but the colors should be even within each bin
        norm = colors.BoundaryNorm(boundaries=levels, ncolors=256) #this is for the case when the levels are not equally distanced, e.g. for w where I also want to display small values
        cf= plt.contourf(ds_node.x, ds_node.y, np.nanmean(ds_node[var_full], axis= 0), levels= levels[1:-1], cmap= cmap, norm=norm, extend='both')      
    else:
        cf= plt.contourf(ds_node.x, ds_node.y, np.mean(ds_node[var_full], axis= 0), levels= levels, cmap= cmap, extend='both', alpha= 0.7)      
    
    if cont == 'Bowen':
#        cont_levels= np.arange(0, 10, 0.1)
#        cont_levels= [1.25**n for n in np.arange(-4, 4.1)]# [0.25, 0.5, 1, 2, 4]
        cont_levels= [1.5**n for n in np.arange(-3, 3.1)]# [0.25, 0.5, 1, 2, 4]

        cs= plt.contour(ds_node.x, ds_node.y, np.mean(ds_node['sshf'], axis= 0)/np.mean(ds_node['slhf']), axis= 0, levels= cont_levels, colors='k')#, linestyles='solid')
        plt.clabel(cs, cont_levels, fontsize=10, fmt='%1.2f')
    else:
        cs= plt.contour(ds_node.x, ds_node.y, np.mean(ds_node[cont_full], axis= 0), levels= cont_levels, colors=cont_color)#, linestyles='solid')
        if cont_label:
            if cont_full== 'msl': plt.clabel(cs, cont_levels[cont_levels>=1000][::2], fontsize=10, fmt='%1.0f')
            elif np.max(np.abs(cont_levels)) > 3: plt.clabel(cs, cont_levels[::2], fontsize=10, fmt='%1.0f')#, inline=1)
            else: plt.clabel(cs, cont_levels[::2], fontsize=10, fmt='%1.1f')#, inline=1)


    if plot_cont2:
        cs2= plt.contour(ds_node.x, ds_node.y, np.mean(ds_node[cont2_full], axis= 0), levels= cont_levels2, colors= color_cont2)

        if np.max(np.abs(np.mean(ds_node[cont2_full], axis= 0) )) > min(np.abs(cont_levels2)):
            if max(np.abs(cont_levels2)) > 3: plt.clabel(cs2, cont_levels2[::2], fontsize=10, fmt='%1.0f')
            else: plt.clabel(cs2, cont_levels2[::2], fontsize=10, fmt='%1.1f')


    if plot_cont3:
        cs3= plt.contour(ds_node.x, ds_node.y, np.mean(ds_node[cont3_full], axis= 0), levels= cont_levels3, colors= 'b')      

#        if np.max(np.abs(np.mean(ds_node[cont3_full], axis= 0) )) > min(np.abs(cont_levels3)):
#            if max(cont_levels3) > 3: plt.clabel(cs3, cont_levels3[::2], fontsize=10, fmt='%1.0f')
#            else: plt.clabel(cs3, cont_levels3[::2], fontsize=10, fmt='%1.1f')

    if plot_wind:
        if arrow_type == 'quiver':
            nth= 4
            plt.quiver(ds_node.x[nth//2:][::nth], ds_node.y[nth//2:][::nth], np.mean(ds_node[wind_full[0]], axis= 0)[nth//2:, nth//2:][::nth,::nth], np.mean(ds_node[wind_full[1]], axis= 0)[nth//2:, nth//2:][::nth,::nth], pivot='mid', width= 0.01)
        if arrow_type == 'barb':
            nth= 5
            plt.barbs(ds_node.x[nth//2:][::nth], ds_node.y[nth//2:][::nth], np.mean(ds_node[wind_full[0]], axis= 0)[nth//2:, nth//2:][::nth,::nth], np.mean(ds_node[wind_full[1]], axis= 0)[nth//2:, nth//2:][::nth,::nth], pivot='middle', length= 5.5, barbcolor= 'k')



    perc_of_SOM= len(np.where(Stoll.node == iSOM)[0])/len(Stoll.node[Stoll.node > 0])
    if vo_tend_excl: plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%, '+ str(int(np.round(after_tend_excl/before_tend_excl *100))) + '%)')
    else: plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')


    ax1.scatter(0,0, color='r')
 
    ax1.set_yticks(np.arange(-250,501, 250))
    ax1.set_xticks(np.arange(-500,501, 250))
    plt.yticks(rotation=90, verticalalignment="center")



    
    """the single Figure"""
    fig2= plt.figure(fignr, figsize= (3,3) )
    fig2.clf()
    
    if levels_different:
        norm = colors.BoundaryNorm(boundaries=levels, ncolors=256) #this is for the case when the levels are not equally distanced, e.g. for w where I also want to display small values
        Cf= plt.contourf(ds_node.x, ds_node.y, np.nanmean(ds_node[var_full], axis= 0), levels= levels[1:-1], cmap= cmap, norm=norm, extend='both')      
    else:
        Cf= plt.contourf(ds_node.x, ds_node.y, np.mean(ds_node[var_full], axis= 0), levels= levels, cmap= cmap, extend='both', alpha= 0.7)      
        
#    Cs= plt.contour(ds_node.x, ds_node.y, np.mean(ds_node[cont_full], axis= 0), levels= cont_levels, colors='k')#, linestyles='solid')

    if cont == 'Bowen':
        Cs= plt.contour(ds_node.x, ds_node.y, np.mean(ds_node['sshf'], axis= 0)/np.mean(ds_node['slhf']), axis= 0, levels= cont_levels, colors='k')#, linestyles='solid')
        plt.clabel(Cs, cont_levels, fontsize=10, fmt='%1.2f')

    else:
        Cs= plt.contour(ds_node.x, ds_node.y, np.mean(ds_node[cont_full], axis= 0), levels= cont_levels, colors=cont_color)#, linestyles='solid')
        if cont_label:
            if cont_full== 'msl': plt.clabel(Cs, cont_levels[cont_levels>=1000][::2], fontsize=10, fmt='%1.0f')
            elif max(np.abs(cont_levels)) > 3: plt.clabel(Cs, cont_levels[::2], fontsize=10, fmt='%1.0f')#, inline=1)
            else: plt.clabel(Cs, cont_levels[::2], fontsize=10, fmt='%1.1f')#, inline=1)

    if plot_cont2:
        Cs2= plt.contour(ds_node.x, ds_node.y, np.mean(ds_node[cont2_full], axis= 0), linestyles= '--', linewidths= 2, levels= cont_levels2, colors= color_cont2)
        if np.max(np.abs(np.mean(ds_node[cont2_full], axis= 0) )) > min(np.abs(cont_levels2)):
            if max(np.abs(cont_levels2)) > 3: plt.clabel(Cs2, cont_levels2[::2], fontsize=10, fmt='%1.0f')
            else: plt.clabel(Cs2, cont_levels2[::2], fontsize=10, fmt='%1.1f')
            
    if plot_cont3:
        Cs3= plt.contour(ds_node.x, ds_node.y, np.mean(ds_node[cont3_full], axis= 0), levels= cont_levels3, colors= 'b')
      
#        if np.max(np.abs(np.mean(ds_node[cont3_full], axis= 0) )) > min(np.abs(cont_levels3)):
#            if max(cont_levels3) > 3: plt.clabel(Cs3, cont_levels3[::2], fontsize=10, fmt='%1.0f')
#            else: plt.clabel(Cs3, cont_levels3[::2], fontsize=10, fmt='%1.1f')


    if plot_wind:
        if arrow_type == 'quiver':
            nth= 4
            Q= plt.quiver(ds_node.x[nth//2:][::nth], ds_node.y[nth//2:][::nth], np.mean(ds_node[wind_full[0]], axis= 0)[nth//2:, nth//2:][::nth,::nth], np.mean(ds_node[wind_full[1]], axis= 0)[nth//2:, nth//2:][::nth,::nth], pivot='mid', width= 0.01, zorder= 2)
            
            if iSOM== 1:
                qk = plt.quiverkey(Q, 0.5, 1.05, 50, r"50 kg m$^{-1}$ s$^{-1}$", coordinates='axes', labelpos='E')
#            qk.text.set_backgroundcolor('w')    

        if arrow_type == 'barb':
            nth= 5
            plt.barbs(ds_node.x[nth//2:][::nth], ds_node.y[nth//2:][::nth], np.mean(ds_node[wind_full[0]], axis= 0)[nth//2:, nth//2:][::nth,::nth], np.mean(ds_node[wind_full[1]], axis= 0)[nth//2:, nth//2:][::nth,::nth], pivot='middle', length= 5.3, barbcolor= 'k')



    plt.scatter(0,0, color='r')
 
    plt.yticks(np.arange(-250,501, 250))
    plt.xticks(np.arange(-500,501, 250))
    plt.yticks(rotation=90, verticalalignment="center")
    fig2.tight_layout()


    if insert_shear_vec:
        S_node= Stoll[Stoll.node== iSOM]
        if vo_tend_excl: S_node= S_node[S_node.vo_tendency > vo_tend_thresh]

        l1= 925
        l2= 500
        mean= 250
        speed_var= 'vert_shear_strength_vec'+str(l1)+'-'+str(l2)+'_mean-'+str(mean)
        dir_var= 'vert_shear_angle_vec3_'+str(l1)+'-850-700-'+str(l2)+'_mean-'+str(mean)
#        dir_var= 'vert_shear_angle_vec'+str(l1)+'-'+str(l2)+'_mean-'+str(mean)        

        S_node[dir_var] = (S_node[dir_var]+90)%360 #to rotate to the right

        u_therm, v_therm= WindSpeedDirection2UV(S_node[speed_var], S_node[dir_var], orientation='')
        u_mean, v_mean= np.mean(u_therm), np.mean(v_therm)
        ax_insert = fig2.add_axes([0.75, 0.75, 0.15, 0.15]) # inset axes
        ax_insert.add_patch(mpatches.FancyArrowPatch((-u_mean/2, -v_mean/2), (u_mean/2, v_mean/2),  mutation_scale=15) )

        ax_insert.set_xlim(-.0014, .0014)
        ax_insert.set_ylim(-.0014, .0014)

        ax_insert.get_xaxis().set_visible(False)
        ax_insert.get_yaxis().set_visible(False)

    
    if save:
        save_name=savedir+ 'SOM_'+ str(iSOM) +'_mean_'+var_full+"_"+cont_full
    
        if plot_wind: save_name += "_"+wind_full[0]
        if plot_cont2: save_name+= "_"+cont2_full 
        if plot_cont3: save_name+= "_"+cont3_full 

        if vo_tend_excl: save_name += '_vo-tend-'+str(vo_tend_thresh)
        save_name += "_x"+str(x)+"_y"+str(y)
    
        print(save_name)
        fig2.savefig(save_name , bbox_inches='tight', dpi= 150)
    
    
    
    
    
    

fig.tight_layout()

fig.subplots_adjust(bottom=0.11)
cbar_ax = fig.add_axes([0.09, 0.06, 0.4, 0.015])
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
elif var == 't': labelvar= 'Temperature [K]'
elif 'U_ano' in var_full: labelvar= 'Wind speed anomaly [m/s]'
elif var == 'U': labelvar= 'Wind speed [m/s]'
elif 'q_ano' in var_full: labelvar= 'Specific humidity anomaly [g(kg)]'
elif var == 'q': labelvar= 'Specific humidity [g/kg]'
elif 'pres_pv' in var_full: labelvar= 'Tropopause level [hPa]'
elif var_full== 'w_700': labelvar= 'Vertical velocity $100^{-1} Pa$\cdot$s$^{-1}$'

elif ano_var== False:
    if '**' in ds[var_full].units:
        units= ds[var_full].units.split('**')
        ds[var_full].attrs['units']= units[0]+ r"$^{" + units[1]+ r"}$"

    labelvar= ds[var_full].long_name + ' ['+ ds[var_full].units+ ']'

else: labelvar= ''

cb.set_label(labelvar, size=12) 


if save:
    save_name=savedir+ 'mean_'+var_full+"_"+cont_full

    if plot_wind: save_name += "_"+wind_full[0]
    if plot_cont2: save_name+= "_"+cont2_full 
    if plot_cont3: save_name+= "_"+cont3_full 

    if vo_tend_excl: save_name += '_vo-tend-'+str(vo_tend_thresh)
    save_name += "_x"+str(x)+"_y"+str(y)

    print(save_name)
    fig.savefig(save_name , bbox_inches='tight')



