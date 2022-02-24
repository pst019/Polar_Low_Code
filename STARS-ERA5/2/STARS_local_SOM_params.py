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
save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM/'
#
fignr= 19

imp_list=True
#imp_list=False


"""SOM var"""
#Obs= "Obsnr1"
#Obs= "Obsnr1_prep"
#Obs='mature_prep'
Obs= 'allObs_SOMprep'

x=3 #3
y=3 #x+1

proplim= 3 #propagation speed limit
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


"""evolve params"""

#Obs_evol=1
#Obs_evol= 'mature'
#mat_evol_step= 100 #in percent how much of time the PL has overcome towards the mature phase (0- initial stage, 100- mature stage)
#mat_evol_step= 25
Obs_evol='allObs'

if Obs_evol== 'mature': Obs_evol_str= Obs_evol +str(mat_evol_step)
elif Obs_evol =='allObs': Obs_evol_str= Obs_evol
else: Obs_evol_str= 'Obsnr'+str(Obs_evol)



"""import Stoll list"""
if imp_list:
#    SOM_filedir= Mediadir + "ERA5_STARS/PL_centred_fields/SOM/"+Svar_full+"_"+str(Splevel)+'_'+Obs+"_x"+str(x)+"_y"+str(y)
    SOM_filedir= Mediadir + "ERA5_STARS/PL_centred_fields/SOM/"+Svar_full+"_"+str(Splevel)+'_'+Obs+'_vel'+str(proplim)+'_dur'+str(lifelim)+"_x"+ str(x)+"_y"+str(y)

    ds_init= xr.open_dataset(SOM_filedir+"_cluster.nc")

    Stoll, Stoll_ind= imp_standard_Stoll(drop=False, rename=False) #SOM_filedir)#, Obs_evol)
    df_nodes= pd.read_csv(SOM_filedir+"_node_nr.txt", sep=" ")


    Stoll= Stoll.rename(columns={
        'U_max_200': 'Wind Speed 10m (max)',
        'vo': 'Vorticity$_{850}$ (centre)',
        'slp': 'Sea-level pressure (min)',
        'blh_max_200_msk': 'Boundary layer height (max)',       
        'cape_max_200_msk': 'CAPE (max)', 
        'skt_med_200': 'Skin temperature (med)',
        'sst_mean_200': 'SST (mean)',

        'sst-t500_max_200': 'SST -T$_{500}$ (max)',
        'sst-t500_mean_200': 'SST -T$_{500}$ (mean)',
        'sst-t700_mean_200': 'SST -T$_{700}$ (mean)',
        'N1000-500_dgauss-4_msk_min-200': 'N$_{500-1000}$ (min)',
        'tp_mean_200_msk': 'Total precip. (mean)',
        'cp_mean_200_msk': 'Convective precip. (mean)',
        'sf_mean_200_msk': 'Snow fall (mean)',
        'lsp_mean_200_msk': 'Large-scale precip. (mean)', 
        'sshf_mean_200_msk': 'Sensible heat flux (mean)',
        'slhf_mean_200_msk': 'Latent heat flux (mean)',

       'grad_t850_max_200': 'Grad T$_{850}$ (max)',
       'baroclinic_dUdz_gr_filter4_925-700_max_200': 'Baroclinic dU/dz (max)',
       'baroclinic_dTdy_gr_filter4_1000-850-500_max_200_msk': 'UL Baroclinic dT/dy (max)',
       'baroclinic_dTdy_gr1000-925-700_dgauss-4_msk_max-200': 'LL Baroclinic dT/dy (max)',

       'barotropic_gr_filter4_850_max_200': 'Barotropic growth (max)',
       'vert_shear_angle925-700_mean_200': 'Vertical shear angle (mean)',
       'vert_wind_grad925-700_mean_200': 'Vertical shear strength (mean)',
       })

    
    if Obs_evol!='allObs':
        df_nodes['ID']= [str(int(df_nodes.PLnr[n]))[:-3]+'_'+str(int(df_nodes.PLnr[n])%1000) for n in range(len(df_nodes.PLnr))]
        df_nodes= df_nodes.set_index('ID')
        Stoll_ind=Stoll_ind.join(df_nodes['node'], how='outer')

  
    if Obs_evol=='allObs':
        df_nodes['time']= pd.to_datetime(df_nodes.date.apply(str), format='%Y%m%d%H')
        df_nodes['ID']= [str(int(df_nodes.PLnr[n]))[:-3]+'_'+str(int(df_nodes.PLnr[n])%1000) for n in range(len(df_nodes.PLnr))]
#        df_nodes= df_nodes.drop(columns=['date', 'PLnr'])
        Stoll= Stoll.set_index(['ID', 'time'])
        df_nodes= df_nodes.set_index(['ID', 'time'])       
        Stoll= pd.concat([Stoll, df_nodes], axis= 1)






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

    perc_of_SOM= len(np.where(df_nodes.node == iSOM)[0])/len(df_nodes.node)
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
    save_name=savedir+ 'SOM_'+Svar_full+"_"+str(Splevel)+"_"+Obs+"_x"+str(x)+"_y"+str(y)
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
            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')', color=color  )

#    yplot -=.1
#    variable=Stoll_node['vert_wind_grad925-700_mean_200']*1E3
#    color= 'k'
#    if variable.median() > 2: color= 'r'
#    if variable.median() < -2: color= 'b'     
#    ax1.text(0,yplot, 'Shear speed: '+ str(np.round(variable.median(),1))+
#            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')' , color=color )


    yplot -=.1
    variable=Stoll_node['UL Baroclinic dT/dy (max)']*24
    color= 'k'
    if variable.median() > 5.5: color= 'r'
    if variable.median() < 4.5: color= 'b'
    ax1.text(0,yplot , 'UL Baroclinic growth: '+ str(np.round(variable.median(),1))+
            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')' , color=color  )

    yplot -=.1
    variable=Stoll_node['LL Baroclinic dT/dy (max)']*24
    color= 'k'
    if variable.median() > 5.5: color= 'r'
    if variable.median() < 4.5: color= 'b'
    ax1.text(0,yplot , 'LL Baroclinic growth: '+ str(np.round(variable.median(),1))+
            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')' , color=color  )


    yplot -=.1
    variable=Stoll_node['Grad T$_{850}$ (max)']
    color= 'k'
    if variable.median() > 5: color= 'r'
    if variable.median() < 3.6: color= 'b'
    ax1.text(0,yplot, 'Grad T$_{850}$: '+ str(np.round(variable.median(),1))+
            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')', color=color   )


#    yplot -=.1
#    variable=Stoll_node['N$_{500-1000}$ (min)']
#    color= 'k'
##    if variable.median() > 5: color= 'r'
##    if variable.median() < 3.6: color= 'b'
#    ax1.text(0,yplot, 'N$_{500-1000}$: '+ str(np.round(variable.median(),3))+
#            ' ('+ str(np.round(np.nanpercentile(variable, 10),3)) + ', '+ str(np.round(np.nanpercentile(variable, 90),3)) +')', color=color   )


#    yplot -=.1    
#    variable=Stoll_node['Barotropic growth (max)']*60**2*24
#    color= 'k'
##    if variable.median() > 2: color= 'r'
##    if variable.median() < 1.3: color= 'b'    
#    ax1.text(0,yplot , 'Barotropic growth: '+ str(np.round(variable.median(),1))+
#            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')' , color=color  )

    yplot -=.1
    variable=Stoll_node['CAPE (max)']
    color= 'k'
    if variable.median() > 130: color= 'r'
    if variable.median() < 90: color= 'b'    
    ax1.text(0,yplot, 'CAPE: '+ str(int(np.round(variable.median())))+
            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')' , color=color )

#    yplot -=.1
#    variable=100*Stoll_node['Convective precip. (mean)']/Stoll_node['Total precip. (mean)']
#    color= 'k'
#    if variable.median() > 80: color= 'r'
#    if variable.median() < 60: color= 'b'   
#    ax1.text(0,yplot, 'Conv precip rate: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')' , color=color )

#    yplot -=.1
#    variable=Stoll_node['Sensible heat flux (mean)']
#    ax1.text(0,yplot, 'Sensible flux: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')' , color=color )

#    yplot -=.1
#    variable=Stoll_node['Latent heat flux (mean)']
#    ax1.text(0,yplot, 'Latent flux: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')' , color=color )

    yplot -=.1
    variable=Stoll_node['Sensible heat flux (mean)'] + Stoll_node['Latent heat flux (mean)']
    color= 'k'
    if variable.median() > 220: color= 'r'
    if variable.median() < 150: color= 'b'      
    ax1.text(0,yplot, 'Surface fluxes: '+ str(int(np.round(variable.median())))+
            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')' , color=color )

#    yplot -=.1
#    variable=Stoll_node['SST -T$_{500}$ (mean)']
#    color= 'k'
##    if variable.median() > 50: color= 'r'
##    if variable.median() < 45: color= 'b' 
#    ax1.text(0,yplot, 'SST- T$_{500}$: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')' , color=color )

#    yplot -=.1
#    variable=Stoll_node['SST -T$_{700}$ (max)']
#    ax1.text(0,yplot, 'SST- T$_{700}$: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')' , color=color )
#
#    yplot -=.1
#    variable=Stoll_node['Boundary layer height (max)']
#    ax1.text(0,yplot, 'BLH: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')' , color=color )

#    yplot -=.1
#    variable=Stoll_node['Vorticity$_{850}$ (centre)']*1e5
#    ax1.text(0,yplot, 'Vorticity: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')' , color=color )


#    yplot -=.1
#    variable=Stoll_node['Wind Speed 10m (max)']
#    color= 'k'
#    if variable.median() > 19: color= 'r'
##    if variable.median() < 16: color= 'b' 
#    ax1.text(0,yplot , 'Wind speed: '+ str(np.round(variable.median(),1))+
#            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')'  , color=color )

#    yplot -=.1
#    variable=Stoll_node['Skin temperature (med)']- 273.15
#    color= 'k'
#    if variable.median() > 6: color= 'r'
#    if variable.median() < 4: color= 'b'     
#    ax1.text(0,yplot, 'Skin temp: '+ str(np.round(variable.median(),1))+
#            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')' , color=color )

    yplot -=.1
    variable=Stoll_node['SST (mean)']- 273.15
    color= 'k'
    if variable.median() > 6: color= 'r'
    if variable.median() < 4: color= 'b'     
    ax1.text(0,yplot, 'SST: '+ str(np.round(variable.median(),1))+
            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')' , color=color )


    
#    yplot -=.1
#    variable=np.sqrt(Stoll_node['area'])*4/np.pi
#    color= 'k'
##    if variable.median() > 6: color= 'r'
##    if variable.median() < 4: color= 'b'     
#    ax1.text(0,yplot, 'Diameter: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')' , color=color )
#    
    
#
#    if Obs_evol == 'allObs':
#        yplot -=.1
#        variable=Stoll_node['Obs']
#        color= 'k'
#        if variable.median() > 25: color= 'r'
#        if variable.median() < 15: color= 'b'     
#        ax1.text(0,yplot, 'Observ number: '+ str(np.round(variable.median(),1))+
#                ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')' , color=color )

#    yplot -=.1
#    variable=Stoll_node['U850_mean_200']
#    color= 'k'
#    if variable.median() > 11: color= 'r'
#    if variable.median() < 8: color= 'b' 
#    ax1.text(0,yplot, 'U$_{850}$: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')' , color=color )


#    yplot -=.1
#    variable=Stoll_node['stearing_vel1000-700_mean_200']
#    color= 'k'
##    if variable.median() > 7: color= 'r'
##    if variable.median() < 6: color= 'b' 
#    ax1.text(0,yplot, 'Prop speed: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')' , color=color )
#

    


#    yplot -=.1
#    variable=Stoll_ind[Stoll_ind.node== iSOM]['Duration']
#    color= 'k'
#    if variable.median() > 35: color= 'r'
#    if variable.median() < 20: color= 'b'     
#    ax1.text(0,yplot, 'Duration: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')' , color=color )

#    yplot -=.1
#    variable=Stoll_ind[Stoll_ind.node== iSOM]['Obs_mature']
#    ax1.text(0,yplot, 'Time to mature: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')' , color=color )

    perc_of_SOM= len(np.where(df_nodes.node == iSOM)[0])/len(df_nodes.node)
    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')
    ax1.axis('off')




plt.tight_layout()


if save:
    save_name=savedir+ 'SOM_'+Svar_full+"_"+str(Splevel)+"_"+Obs+'--other_vars_'+Obs_evol_str+"2_x"+str(x)+"_y"+str(y)+'_msk'
    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')




"""cloud distributions"""
fig= plt.figure(fignr, figsize= (2.7*x,2.7*y) )
fignr+=1
plt.clf()

if Obs_evol != 'allObs':
    morphlist= list(Stoll_ind.groupby('Morphology').size().sort_values(ascending=False).index)

if Obs_evol == 'allObs':

    Sdist= 100
    #only the steps where the cloud morphology is observed and the distance is lower than Sdist
    Stoll= Stoll.loc[Stoll['STARS lat'] > 1]
    Stoll= Stoll.loc[distance((Stoll.lat, Stoll.lon), (Stoll['STARS lat'], Stoll['STARS lon'])) < Sdist]
#    Stoll= remove_morphs(Stoll, 'Morphology', count= 25, how='remove', excludelist=['-', 'U'])

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
#        Morph_node= Stoll[Stoll.node== iSOM]['Morphology']
        Morph_node= Stoll[Stoll.node== iSOM]['Morph_red']


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
    save_name=savedir+ 'SOM_'+Svar_full+"_"+str(Splevel)+"_"+Obs+'--cloud_morphs'+str(x)+"_y"+str(y)
    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')

