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
#save= True
save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM/'
#
fignr= 10

imp_list=True
#imp_list=False


"""SOM var"""
#Obs= "Obsnr1"
#Obs= "Obsnr1_prep"
#Obs='mature_prep'
Obs= 'allObs_SOMprep'

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
#Obs_evol= 'mature'
#mat_evol_step= 100 #in percent how much of time the PL has overcome towards the mature phase (0- initial stage, 100- mature stage)
#mat_evol_step= 25
Obs_evol='allObs'

if Obs_evol== 'mature': Obs_evol_str= Obs_evol +str(mat_evol_step)
elif Obs_evol =='allObs': Obs_evol_str= Obs_evol
else: Obs_evol_str= 'Obsnr'+str(Obs_evol)


x=3 #3
y=x+1


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
        'tp_mean_200': 'Total precip. (mean)',
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
#    df_dist= pd.read_csv(filedir+"_distances.txt", sep=" ")
#    df_nodes= pd.concat([df_nodes, df_dist], axis= 1)

    if Obs_evol!='allObs':
        df_nodes['ID']= [str(int(df_nodes.PLnr[n]))[:-3]+'_'+str(int(df_nodes.PLnr[n])%1000) for n in range(len(df_nodes.PLnr))]
        df_nodes= df_nodes.set_index('ID')
#        df_nodes= df_nodes.drop(columns='as.character(dates)')

        Stoll_ind=Stoll_ind.join(df_nodes['node'], how='outer')
#        Stoll_ind=Stoll_ind.join(df_nodes['K_SOM.distances'], how='outer')

    
    if Obs_evol=='allObs':
        df_nodes['time']= pd.to_datetime(df_nodes.date.apply(str), format='%Y%m%d%H')
        df_nodes['ID']= [str(int(df_nodes.PLnr[n]))[:-3]+'_'+str(int(df_nodes.PLnr[n])%1000) for n in range(len(df_nodes.PLnr))]

#        df_nodes= df_nodes.drop(columns=['date', 'PLnr'])
        Stoll= Stoll.set_index(['ID', 'time'])
        df_nodes= df_nodes.set_index(['ID', 'time'])
        
        Stoll= pd.concat([Stoll, df_nodes], axis= 1)
    






"""make the figure on which the original SOMs are displayed"""
fig= plt.figure(fignr, figsize= (2.5*x,2.5*y +1) )
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


plt.tight_layout()

fig.subplots_adjust(bottom=0.1)
cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.015])
cb= fig.colorbar(cf, cax=cbar_ax, orientation="horizontal")

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




"""the other parameters"""
fig= plt.figure(fignr, figsize= (2.7*x,2.7*y) )
fignr+=1
plt.clf()


if type(Obs_evol)== int: Stoll_now= Stoll[Stoll.Obs== Obs_evol]
elif Obs_evol== 'mature':
    mature_index= [np.where(np.logical_and(Stoll.ID== Stoll_ind.ID[n], Stoll.Obs== Stoll_ind.Obs_mature[n]))[0][0] for n in range(len(Stoll_ind.ID)) ]
    Stoll_now=Stoll.iloc[mature_index]
#elif Obs_evol == 'allObs':
    


for iSOM in range(1, x*y +1):
    ax1= plt.subplot(y, x, iSOM)

    if Obs_evol != 'allObs':
        node_list= Stoll_ind[Stoll_ind.node== iSOM].index.values
        Stoll_node= Stoll_now[Stoll_now.ID.isin(node_list)]
    
    if Obs_evol == 'allObs':
        Stoll_node= Stoll[Stoll.node== iSOM]
        
    yplot= 0.9
    
    variable=Stoll_node['Vertical shear angle (mean)']
    color= 'k'
    if variable.median() < 45: color= 'r'
    if variable.median() > 135: color= 'b'
    ax1.text(0,yplot, 'Shear angle: '+ str(int(np.round(variable.median())))+
            ' ('+ str(int(np.round(np.percentile(variable, 10)))) + ', '+ str(int(np.round(np.percentile(variable, 90)))) +')', color=color  )

#    yplot -=.1
#    variable=Stoll_node['vert_wind_grad925-700_mean_200']*1E3
#    color= 'k'
#    if variable.median() > 2: color= 'r'
#    if variable.median() < -2: color= 'b'     
#    ax1.text(0,yplot, 'Shear speed: '+ str(np.round(variable.median(),1))+
#            ' ('+ str(np.round(np.percentile(variable, 10),1)) + ', '+ str(np.round(np.percentile(variable, 90),1)) +')' , color=color )


    yplot -=.1
    variable=Stoll_node['Baroclinic growth (max)']*60**2*24
    color= 'k'
    if variable.median() > 3.5: color= 'r'
    if variable.median() < 2.8: color= 'b'
    ax1.text(0,yplot , 'Baroclinic growth: '+ str(np.round(variable.median(),1))+
            ' ('+ str(np.round(np.percentile(variable, 10),1)) + ', '+ str(np.round(np.percentile(variable, 90),1)) +')' , color=color  )

    yplot -=.1
    variable=Stoll_node['Grad T$_{850}$ (max)']
    color= 'k'
    if variable.median() > 5: color= 'r'
    if variable.median() < 3: color= 'b'
    ax1.text(0,yplot, 'Grad T$_{850}$: '+ str(np.round(variable.median(),1))+
            ' ('+ str(np.round(np.percentile(variable, 10),1)) + ', '+ str(np.round(np.percentile(variable, 90),1)) +')', color=color   )


    yplot -=.1    
    variable=Stoll_node['Barotropic growth (max)']*60**2*24
    color= 'k'
#    if variable.median() > 2: color= 'r'
#    if variable.median() < 1.3: color= 'b'    
    ax1.text(0,yplot , 'Barotropic growth: '+ str(np.round(variable.median(),1))+
            ' ('+ str(np.round(np.percentile(variable, 10),1)) + ', '+ str(np.round(np.percentile(variable, 90),1)) +')' , color=color  )

    yplot -=.1
    variable=Stoll_node['CAPE (max)']
    color= 'k'
    if variable.median() > 130: color= 'r'
    if variable.median() < 90: color= 'b'    
    ax1.text(0,yplot, 'CAPE: '+ str(int(np.round(variable.median())))+
            ' ('+ str(int(np.round(np.percentile(variable, 10)))) + ', '+ str(int(np.round(np.percentile(variable, 90)))) +')' , color=color )

#    yplot -=.1
#    variable=100*Stoll_node['Convective precip. (mean)']/Stoll_node['Total precip. (mean)']
#    color= 'k'
#    if variable.median() > 80: color= 'r'
#    if variable.median() < 60: color= 'b'   
#    ax1.text(0,yplot, 'Conv precip rate: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.percentile(variable, 10)))) + ', '+ str(int(np.round(np.percentile(variable, 90)))) +')' , color=color )

#    yplot -=.1
#    variable=Stoll_node['Sensible heat flux (mean)']
#    ax1.text(0,yplot, 'Sensible flux: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.percentile(variable, 10)))) + ', '+ str(int(np.round(np.percentile(variable, 90)))) +')' , color=color )

#    yplot -=.1
#    variable=Stoll_node['Latent heat flux (mean)']
#    ax1.text(0,yplot, 'Latent flux: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.percentile(variable, 10)))) + ', '+ str(int(np.round(np.percentile(variable, 90)))) +')' , color=color )

    yplot -=.1
    variable=Stoll_node['Sensible heat flux (mean)'] + Stoll_node['Latent heat flux (mean)']
    color= 'k'
    if variable.median() > 220: color= 'r'
    if variable.median() < 150: color= 'b'      
    ax1.text(0,yplot, 'Surface fluxes: '+ str(int(np.round(variable.median())))+
            ' ('+ str(int(np.round(np.percentile(variable, 10)))) + ', '+ str(int(np.round(np.percentile(variable, 90)))) +')' , color=color )

    yplot -=.1
    variable=Stoll_node['SST- T$_{500}$ (max)']
    color= 'k'
    if variable.median() > 50: color= 'r'
    if variable.median() < 45: color= 'b' 
    ax1.text(0,yplot, 'SST- T$_{500}$: '+ str(int(np.round(variable.median())))+
            ' ('+ str(int(np.round(np.percentile(variable, 10)))) + ', '+ str(int(np.round(np.percentile(variable, 90)))) +')' , color=color )

#    yplot -=.1
#    variable=Stoll_node['SST -T$_{700}$ (max)']
#    ax1.text(0,yplot, 'SST- T$_{700}$: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.percentile(variable, 10)))) + ', '+ str(int(np.round(np.percentile(variable, 90)))) +')' , color=color )
#
#    yplot -=.1
#    variable=Stoll_node['Boundary layer height (max)']
#    ax1.text(0,yplot, 'BLH: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.percentile(variable, 10)))) + ', '+ str(int(np.round(np.percentile(variable, 90)))) +')' , color=color )
#
#    yplot -=.1
#    variable=Stoll_node['Vorticity$_{850}$ (centre)']*1e5
#    ax1.text(0,yplot, 'Vorticity: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.percentile(variable, 10)))) + ', '+ str(int(np.round(np.percentile(variable, 90)))) +')' , color=color )


#    yplot -=.1
#    variable=Stoll_node['Wind Speed 10m (max)']
#    color= 'k'
#    if variable.median() > 19: color= 'r'
##    if variable.median() < 16: color= 'b' 
#    ax1.text(0,yplot , 'Wind speed: '+ str(np.round(variable.median(),1))+
#            ' ('+ str(np.round(np.percentile(variable, 10),1)) + ', '+ str(np.round(np.percentile(variable, 90),1)) +')'  , color=color )

    yplot -=.1
    variable=Stoll_node['Skin temperature (med)']- 273.15
    color= 'k'
    if variable.median() > 6: color= 'r'
    if variable.median() < 4: color= 'b'     
    ax1.text(0,yplot, 'Skin temp: '+ str(np.round(variable.median(),1))+
            ' ('+ str(np.round(np.percentile(variable, 10),1)) + ', '+ str(np.round(np.percentile(variable, 90),1)) +')' , color=color )

    
#    yplot -=.1
#    variable=np.sqrt(Stoll_node['area'])*4/np.pi
#    color= 'k'
##    if variable.median() > 6: color= 'r'
##    if variable.median() < 4: color= 'b'     
#    ax1.text(0,yplot, 'Diameter: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.percentile(variable, 10)))) + ', '+ str(int(np.round(np.percentile(variable, 90)))) +')' , color=color )
#    
    
#
#    if Obs_evol == 'allObs':
#        yplot -=.1
#        variable=Stoll_node['Obs']
#        color= 'k'
#        if variable.median() > 25: color= 'r'
#        if variable.median() < 15: color= 'b'     
#        ax1.text(0,yplot, 'Observ number: '+ str(np.round(variable.median(),1))+
#                ' ('+ str(np.round(np.percentile(variable, 10),1)) + ', '+ str(np.round(np.percentile(variable, 90),1)) +')' , color=color )

#    yplot -=.1
#    variable=Stoll_node['U850_mean_200']
#    color= 'k'
#    if variable.median() > 11: color= 'r'
#    if variable.median() < 8: color= 'b' 
#    ax1.text(0,yplot, 'U$_{850}$: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.percentile(variable, 10)))) + ', '+ str(int(np.round(np.percentile(variable, 90)))) +')' , color=color )


#    yplot -=.1
#    variable=Stoll_node['stearing_vel1000-700_mean_200']
#    color= 'k'
##    if variable.median() > 7: color= 'r'
##    if variable.median() < 6: color= 'b' 
#    ax1.text(0,yplot, 'Prop speed: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.percentile(variable, 10)))) + ', '+ str(int(np.round(np.percentile(variable, 90)))) +')' , color=color )
#

    


#    yplot -=.1
#    variable=Stoll_ind[Stoll_ind.node== iSOM]['Duration']
#    color= 'k'
#    if variable.median() > 35: color= 'r'
#    if variable.median() < 20: color= 'b'     
#    ax1.text(0,yplot, 'Duration: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.percentile(variable, 10)))) + ', '+ str(int(np.round(np.percentile(variable, 90)))) +')' , color=color )

#    yplot -=.1
#    variable=Stoll_ind[Stoll_ind.node== iSOM]['Obs_mature']
#    ax1.text(0,yplot, 'Time to mature: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.percentile(variable, 10)))) + ', '+ str(int(np.round(np.percentile(variable, 90)))) +')' , color=color )

    perc_of_SOM= len(np.where(df_nodes.node == iSOM)[0])/len(df_nodes.node)
    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')
    ax1.axis('off')




plt.tight_layout()


if save:
    save_name=savedir+ 'SOM_'+var_full+"_"+str(plevel)+"_"+Obs+'--other_vars_'+Obs_evol_str+"2_x"+str(x)+"_y"+str(y)
    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')




"""cloud distributions"""
fig= plt.figure(fignr, figsize= (2.7*x,2.7*y) )
fignr+=1
plt.clf()

if Obs_evol != 'allObs':
    morphlist= list(Stoll_ind.groupby('Morphology').size().sort_values(ascending=False).index)

if Obs_evol == 'allObs':
    #should maybe adapt this
    morphlist= list(Stoll.groupby('Morphology').size().sort_values(ascending=False).index)[:7]

#if type(Obs_evol)== int: Stoll_now= Stoll[Stoll.Obs== Obs_evol]
#else:
#    mature_index= [np.where(np.logical_and(Stoll.ID== Stoll_ind.ID[n], Stoll.Obs== Stoll_ind.Obs_mature[n]))[0][0] for n in range(len(Stoll_ind.ID)) ]
#    Stoll_now=Stoll.ix[mature_index]



for iSOM in range(1, x*y +1):
    ax1= plt.subplot(y, x, iSOM)
    
    if Obs_evol != 'allObs':
        Morph_node= Stoll_ind[Stoll_ind.node== iSOM]['Morphology']
    if Obs_evol == 'allObs':
        Morph_node= Stoll[Stoll.node== iSOM]['Morphology']


    for morph in morphlist:
#        plt.hist(Sx[Sx['Morphology'] == morph], bins= bins , label=morph, bottom= bot, edgecolor= 'k', lw= 1)
        ax1.bar(morph, len(Morph_node[Morph_node == morph]), edgecolor= 'k', lw= 1)

    perc_of_SOM= len(np.where(df_nodes.node == iSOM)[0])/len(df_nodes.node)
    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')
#    ax1.axis('off')
#
#
#
plt.tight_layout()

if save:
    save_name=savedir+ 'SOM_'+var_full+"_"+str(plevel)+"_"+Obs+'--cloud_morphs'+str(x)+"_y"+str(y)
    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')

