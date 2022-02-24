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
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/diff_SOMs/'
#
fignr= 3

imp_list=True
#imp_list=False

imp_rot_nc=True
#imp_rot_nc=False

"""SOM var"""
Obs= "Obs1_SOMprep"
Obs= "Obsmature_SOMprep"
#Obs= "Obslast_SOMprep"
#Obs= 'allObs_SOMprep'

lifelim= 6
wind_thresh = False
matchdist= False
Nmatches= False

test=''
#test = 'primary+strong'
#test='strict_matching'


if test == 'primary+strong':
    Obs= 'allObs_primary_SOMprep'
    lifelim= 12 #lifetime limit
    wind_thresh= 20 #None    

if test == 'strict_matching':
    matchdist= 75
    Nmatches= 5


PLCG_type= 'track_smth'

x=3 #3
y=3 #x+1

mirror_top_bottom=False
mirror_left_right=False
mirror_diag_topleft_bottomright=False
rotate_clock=False
rotate_anti_clock=False



#PLCG_type= 'track_smth'
smooth_param= '1E-3'





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


if Svar_full=='t_ano':
    if x== 3 and y == 4: mirror_top_bottom=True
    if x== 3 and y == 3:
        if test == '' and Obs == 'allObs_SOMprep': rotate_clock=True
        elif test == 'strict_matching':
            mirror_top_bottom=True
            mirror_left_right=True
            mirror_diag_topleft_bottomright=True
        elif Obs== "Obs1_SOMprep":
            rotate_anti_clock=True
            mirror_diag_topleft_bottomright=True
        elif Obs== "Obslast_SOMprep":
            rotate_clock=True
        elif Obs== "Obsmature_SOMprep":
            mirror_top_bottom=True
            mirror_left_right=True            
            mirror_diag_topleft_bottomright=True

#    if Splevel== 1000: should diagnonally mirror

if Svar_full=='t' and Splevel== 850:
    mirror_top_bottom= True

if Svar_full=='z_ano' and Splevel== 500:
    rotate_anti_clock= True

if Svar_full=='z_ano' and Splevel== 700:
    mirror_top_bottom= True

if Svar_full=='z_ano' and Splevel== 1000:
    mirror_left_right= True
#    mirror_top_bottom= True
    mirror_diag_topleft_bottomright=True

    
if Svar_full=='q_ano' and Splevel== 850:
    mirror_diag_topleft_bottomright=True


"""evolve var"""
Obs_evol='allObs'


var, ano_var, plevel_var= Svar, Sano, Splevel



#ano_var=True
#ano_var=False
#var, cmap='t', 'RdBu'
#var, cmap='z', 'RdBu'
#var, cmap='q', 'Blues'
#var, cmap='vo', 'RdBu_r'
#var='d'
#var, cmap='u', 'RdBu'
#var, cmap='v', 'RdBu'

#var, cmap='w', 'RdBu'
#var, cmap='pv', 'Reds'
#plevel_var=850






var_full= var
if ano_var==True:  var_full= var_full+ '_ano'
if type(plevel_var) == int: var_full= var_full+ '_'+str(plevel_var)
if type(plevel_var) == list:
    var_full += ('').join(['-'+str(p) for p in plevel_var])


"""specify the levels"""
levels = []
levels_different = False #normally this is false, only if the levels are differently spaced, e.g. logarithmic it should be true - in this case the first and last dispayed level are cut for the "extent" to work

if var_full== 't_1000': levels= np.arange(258, 274.1, 2)
elif var_full== 't_850': levels= np.arange(251, 266.1, 1)
elif var_full== 't_500': levels= np.arange(230, 240.1, 2)
elif var_full== 'z_500': levels= np.arange(5000, 5160.1, 20)
elif var_full== 'z_1000': levels= np.arange(-60, 100.1, 20)
elif var_full== 'q_850': levels= np.arange(0.6, 1.71, .1)
elif var_full== 'q_ano_850': levels= np.arange(-0.6, 0.61, .1)





elif var in ['u', 'v']: levels= np.arange(-20, 20.1, 2)
elif 'U_' in var_full: levels= np.arange(0, 20.1, 2)
elif var_full == '10U': levels= np.arange(6, 16.1, 1)

#elif var_full== 'N925-500': levels= np.arange(0.0042, .0067, 0.0002)

elif 'N' in var_full: levels= np.arange(0.0042, .0067, 0.0002)
#elif 'adv_t_' in  var_full: levels= np.arange(-15, 15, 1)*1E-5


elif 't_ano' in var_full:  levels = np.linspace(-7, 7, 15)
elif 'z_ano' in var_full:  levels = np.arange(-200, 201, 20)
elif var==  'vo':  levels = np.arange(-50, 50.1, 5)

    


scale= .7 #scale factor for the color bar. it takes scale* [min, max] of the all time step values, only relevant if levels = []
sym= False # only relevant if levels=[]
#cmap='jet'
#cmap= 'viridis'
#cmap= 'Greys_r'
#cmap= 'RdBu_r'

if ano_var: sym= True
if sym: cmap= 'RdBu_r'






"""import SOM result"""
if imp_list:
    if PLCG_type== 'stearing_flow':
        SOM_filedir= Mediadir + "ERA5_STARS/PL_centred_fields/SOM/"+Svar_full+"_"+str(Splevel)+'_'+Obs+'_vel'+str(proplim)+'_dur'+str(lifelim)+"_x"+ str(x)+"_y"+str(y)
    elif PLCG_type== 'track_smth':
        SOM_filedir= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/SOM/"+Svar_full+"_"+str(Splevel)+'_'+Obs+ '_track-smth-'+smooth_param+'_dur'+str(lifelim)
        if wind_thresh: SOM_filedir += '_Umax'+str(wind_thresh)
        if matchdist: SOM_filedir += '_matchdist'+str(matchdist)
        if Nmatches: SOM_filedir += '_Nmatches' + str(Nmatches)
        SOM_filedir +="_x"+ str(x)+"_y"+str(y)

    ds_init= xr.open_dataset(SOM_filedir+"_cluster.nc")

#    Stoll, Stoll_ind= imp_standard_Stoll() #SOM_filedir)#, Obs_evol)
    df_nodes= pd.read_csv(SOM_filedir+"_node_nr.txt", sep=" ")

    

    df_nodes['time']= pd.to_datetime(df_nodes.date.apply(str), format='%Y%m%d%H')
    df_nodes['ID']= [str(int(df_nodes.PLnr[n]))[:-3]+'_'+str(int(df_nodes.PLnr[n])%1000) for n in range(len(df_nodes.PLnr))]
#        df_nodes= df_nodes.drop(columns=['date', 'PLnr'])
#    Stoll= Stoll.set_index(['ID', 'time'])
    df_nodes= df_nodes.set_index(['ID', 'time'])
    
    """to rotate/mirror the nodes"""
    df_nodes_orig= df_nodes.copy()
    if mirror_top_bottom:
        df_nodes['node']= (y-1 -(df_nodes.node.values-1)//x)*x + (df_nodes.node.values -1)%x +1
    if mirror_left_right:
        df_nodes['node']= ((df_nodes.node.values-1)//x +1)*x - (df_nodes.node.values-1)%x

    if mirror_diag_topleft_bottomright:
        x_n= (df_nodes['node'].values-1)%x +1
        y_n= (df_nodes['node'].values-1)//x +1
#        y_n= y+1 - y_n #reverse the y_nodes

        df_nodes['node']= (x_n-1 )*(x) + y_n
    
    
    if rotate_clock:
        x_n= (df_nodes['node'].values-1)%x +1
        y_n= (df_nodes['node'].values-1)//x +1
        y_n= y+1 - y_n #reverse the y_nodes
        
        df_nodes['node']= (x_n -1)*x + y_n

    if rotate_anti_clock:
        x_n= (df_nodes['node'].values-1)%x +1
        y_n= (df_nodes['node'].values-1)//x +1
        x_n= x+1 - x_n #reverse the x_nodes
        
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

elif Svar_full== 'q_ano': Slabelvar= 'Specific humidity anomaly [g (kg)$^{-1}$]'
elif Svar_full== 'q': Slabelvar= 'Specific humidity [g/kg]'
else: Slabelvar= Svar_full

cb.set_label(Slabelvar, size=14) 

if save:
    if PLCG_type== 'stearing_flow': save_name=savedir+ 'SOM_'+Svar_full+"_"+str(Splevel)+"_"+Obs+"_x"+str(x)+"_y"+str(y)
    elif PLCG_type== 'track_smth': save_name=savedir+ 'SOM_track-smth'+smooth_param+'_'+Svar_full+"_"+str(Splevel)+"_"+Obs+"_x"+str(x)+"_y"+str(y)
    
    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')









"""evolution"""
if imp_rot_nc:
    print('import from all time step list')
    
    if plevel_var != '': plev_str= '_'+ str(plevel_var)
    else: plev_str= plevel_var
    
    file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +var + plev_str + '_allObs_track-smth-'+smooth_param+'.nc'


    if os.path.isfile(file):
        ds= xr.open_dataset(file)

    else:
        ds= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +var + '_all_levs_allObs_track-smth-'+smooth_param+'.nc')
        ds= ds.sel(plev= plevel_var)
        
        if 'd' in list(ds.keys()): #for the variables with plevels
            ds['d']*= 1E5
            ds['d'].attrs['units'] = '1e-5 '+ ds['d'].attrs['units']
        
    ds= ds.rename({var: var_full}) #this should not include the anomaly



    


    """add the node to ds - only for allObs"""
    ds_node_var= np.zeros(len(ds.time))
    
    for s in range(len(df_nodes)):
        df_nodes_s= df_nodes.iloc[s]
        ds_index= np.argwhere(np.logical_and(ds.PLnr.values== df_nodes_s.PLnr, ds.time.values== np.datetime64(df_nodes_s.name[1])))[0][0]
        ds_node_var[ds_index]= df_nodes_s.node
    
    ds['node']= (('time'), ds_node_var)

    

    if ano_var == True:
        ds[var_full]= (('time', 'x', 'y'), ds[var_full].values- np.nanmean(ds[var_full].values, axis= (1,2))[:, np.newaxis, np.newaxis])       



if sym: vextr= scale* np.max([np.max(ds[var_full]), -np.min(ds[var_full])])
else: vmax, vmin= scale*float(np.max(ds[var_full])), scale*float(np.min(ds[var_full]))


print('plot the SOMs')
#fig= plt.figure(fignr, figsize= (2.7*x,2.7*y) )
fig= plt.figure(fignr, figsize= (2.5*x,2.5*y +0) )

fignr+=1
plt.clf()


for iSOM in range(1, x*y +1):
    print('SOM ', iSOM)
    ax1= plt.subplot(y, x, iSOM)
    
    
    ds_node= ds.where(ds.node == iSOM, drop=True)
    
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



    perc_of_SOM= len(ds_node.node)/len(np.where(ds.node > 0)[0])
    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')



    ax1.scatter(0,0, color='r', zorder= 2)

    size= 500    
    ax1.set_xlim(-size, size)
    ax1.set_ylim(-size, size)

    ax1.set_yticks(np.arange(-size,size+1, size/2))
    ax1.set_xticks(np.arange(-size,size+1, size/2))
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
elif var == 'z': labelvar= 'Geopotential height [m]'
elif var == 'z_pv': labelvar= 'Tropopause height [m]'

elif var == 'vo': labelvar= 'Relative vorticity [10$^{-5}$ s$^{-1}$]'


elif 't_ano' in var_full: labelvar= 'Temperature anomaly [K]'
elif var =='t': labelvar= 'Temperature [K]'
elif 'U_ano' in var_full: labelvar= 'Wind speed anomaly [m/s]'
elif var == 'U' : labelvar= 'Wind speed [m/s]'
elif 'q_ano' in var_full: labelvar= 'Specific humidity anomaly [g (kg)$^{-1}$]'
elif var == 'q': labelvar= 'Specific humidity [g/kg]'


elif ano_var== False: labelvar= ds[var_full].long_name + ' ['+ ds[var_full].units+ ']'
else: labelvar= ''

cb.set_label(labelvar, size=12) 


if save:
    save_name=savedir+ 'SOM_track-smth'+smooth_param+'_'+Svar_full+"_"+str(Splevel)+"_"+Obs

    if matchdist: save_name += '_matchdist'+str(matchdist)
    if Nmatches: save_name += '_Nmatches' + str(Nmatches)
       
    save_name+='_size'+str(size)+'--mean_'+var_full+"_x"+str(x)+"_y"+str(y)


    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')



