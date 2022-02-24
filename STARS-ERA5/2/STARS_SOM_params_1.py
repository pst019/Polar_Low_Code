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
fignr= 9

imp_list=True
#imp_list=False


cloud_hist= 'percentage' #[percentage, number] - determines if the cloud morphology histogram is displayed in number of time steps or in percentage within the SOM

"""SOM var"""
#Obs= "Obsnr1"
#Obs= "Obsnr1_prep"
#Obs='mature_prep'
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
    if PLCG_type== 'stearing_flow':
        SOM_filedir= Mediadir + "ERA5_STARS/PL_centred_fields/SOM/"+Svar_full+"_"+str(Splevel)+'_'+Obs+'_vel'+str(proplim)+'_dur'+str(lifelim)+"_x"+ str(x)+"_y"+str(y)
    elif PLCG_type== 'track_smth':
        SOM_filedir= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/SOM/"+Svar_full+"_"+str(Splevel)+'_'+Obs+ '_track-smth-'+smooth_param+'_dur'+str(lifelim)+"_x"+ str(x)+"_y"+str(y)

    ds_init= xr.open_dataset(SOM_filedir+"_cluster.nc")

    Stoll, Stoll_ind= imp_standard_Stoll(drop=False, ERA5_file= "PL/STARS/Stoll_ERA5_handfix_3.csv",
                                         rename=False) #SOM_filedir)#, Obs_evol) #rename - of the ERA5 parameters
    df_nodes= pd.read_csv(SOM_filedir+"_node_nr.txt", sep=" ")


    Stoll= Stoll.rename(columns={
        'U_msk_max-200': 'Wind Speed 10m (max)',
        'vo': 'Vorticity$_{850}$ (centre)',
        'slp': 'Sea-level pressure (min)',
        'blh_msk_mean-200': 'Boundary layer height (mean)*',       
        'cape_msk_mean-200': 'CAPE (mean)*', 

#        'skt_med-200': 'Skin temperature (med)',
        'sst_msk_mean-200': 'SST (mean)*',

#       'N1000-700_dgauss-4_msk_mean-200': 'N$_{700-1000}$ (mean)',
       'N1000-500_msk_mean-200': 'N$_{500-1000}$ (mean)*',
       'N925-500_dgauss-4_msk_mean-200': 'N$_{500-925}$ (mean)*',
       'N1000-850_msk_mean-200': 'N$_{850-1000}$ (mean)*',
       
       'sst-t500_msk_mean-200': 'SST -T$_{500}$ (mean)*',
       'sst-t700_msk_mean-200': 'SST -T$_{700}$ (mean)*',
       
        'tp_msk_mean-200': 'Total precip. (mean)*',
        'cp_msk_mean-200': 'Convective precip. (mean)*',
        'sf_msk_mean-200': 'Snow fall (mean)*',
        'lsp_msk_mean-200': 'Large-scale precip. (mean)*', 
        'sshf_msk_mean-200': 'Sensible heat flux (mean)*',
        'slhf_msk_mean-200': 'Latent heat flux (mean)*',

       'grad_t850_max-200': 'Grad T$_{850}$ (max)',
       'grad_t850_dgauss-4_msk_max-200': 'Grad T$_{850}$ (max)*',
#       'baroclinic_dUdz_gr_filter4_925-700_max-200': 'Baroclinic dU/dz (max)',
#       'baroclinic_dTdy_gr_filter4_1000-850-500_max-200_msk': 'UL Baroclinic dT/dy (max)',
#       'baroclinic_dTdy_gr1000-925-700_dgauss-4_msk_max-200': 'LL Baroclinic dT/dy (max)',
       'baroclinic_dTdy_gr1000-925-850_dgauss-4_msk_mean-200': 'LL baroclinic (mean)*',
#       'baroclinic_dTdy_gr1000-850-500_dgauss-4_msk_mean-200': 'Deep baroclinic (mean)*',
       'baroclinic_dTdy_gr925-850-500_dgauss-4_msk_mean-200': 'Deep baroclinic (mean)*',

#       'barotropic_gr_filter4_850_max-200': 'Barotropic growth (max)',
       'vert_shear_angle925-700_mean-400': 'Vertical shear angle (mean)',
       'vert_shear_strength925-700_mean-200': 'Vertical shear strength (mean)',
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


    perc_of_SOM= len(np.where(df_nodes.node == iSOM)[0])/len(df_nodes.node)
    ax1.set_title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')
    ax1.axis('off')
        
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
    variable=Stoll_node['Deep baroclinic (mean)*']*24
    color= 'k'
    if variable.median() > 3: color= 'r'
    if variable.median() < 2.5: color= 'b'
    ax1.text(0,yplot , 'Deep baroclinic gr: '+ str(np.round(variable.median(),1))+
            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')' , color=color  )

    yplot -=.1
    variable=Stoll_node['LL baroclinic (mean)*']*24
    color= 'k'
    if variable.median() > 4: color= 'r'
    if variable.median() < 2.5: color= 'b'
    ax1.text(0,yplot , 'LL baroclinic gr: '+ str(np.round(variable.median(),1))+
            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')' , color=color  )


    yplot -=.1
    variable=Stoll_node['Grad T$_{850}$ (max)']
    color= 'k'
    if variable.median() > 5: color= 'r'
    if variable.median() < 3.6: color= 'b'
    ax1.text(0,yplot, r"$ \nabla$ T$_{850}$: "+ str(np.round(variable.median(),1))+
            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')', color=color   )


    yplot -=.1
    variable=Stoll_node['N$_{500-925}$ (mean)*']
    color= 'k'
    if variable.median() < 0.0042: color= 'r'
    if variable.median() > 0.0047: color= 'b'
#    ax1.text(0,yplot, 'N$_{500-1000}$: '+ str(np.round(variable.median(),4))+
#            ' ('+ str(np.round(np.nanpercentile(variable, 10),4)) + ', '+ str(np.round(np.nanpercentile(variable, 90),4)) +')', color=color   )
    ax1.text(0,yplot, 'N: '+ str(np.round(variable.median(),4))+
            ' ('+ str(np.round(np.nanpercentile(variable, 10),4)) + ', '+ str(np.round(np.nanpercentile(variable, 90),4)) +')', color=color   )


#    yplot -=.1
#    variable=Stoll_node['N$_{850-1000}$ (mean)*']
#    color= 'k'
##    if variable.median() > 5: color= 'r'
##    if variable.median() < 3.6: color= 'b'
#    ax1.text(0,yplot, 'N$_{850-1000}$: '+ str(np.round(variable.median(),4))+
#            ' ('+ str(np.round(np.nanpercentile(variable, 10),4)) + ', '+ str(np.round(np.nanpercentile(variable, 90),4)) +')', color=color   )
#

#    yplot -=.1    
#    variable=Stoll_node['Barotropic growth (max)']*60**2*24
#    color= 'k'
##    if variable.median() > 2: color= 'r'
##    if variable.median() < 1.3: color= 'b'    
#    ax1.text(0,yplot , 'Barotropic growth: '+ str(np.round(variable.median(),1))+
#            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')' , color=color  )

    yplot -=.1
    variable=Stoll_node['CAPE (mean)*']
    color= 'k'
    if variable.median() > 30: color= 'r'
    if variable.median() < 21: color= 'b'    
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
    variable=Stoll_node['Sensible heat flux (mean)*'] + Stoll_node['Latent heat flux (mean)*']
    color= 'k'
    if variable.median() > 250: color= 'r'
    if variable.median() < 200: color= 'b'      
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
    variable=Stoll_node['SST (mean)*']- 273.15
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


    """find out what fraction is occuring with another PL time step simultaneously
    and with less than 500km distance"""
    time_list_node= Stoll_node.index.get_level_values(1)
    time_list_total= Stoll.index.get_level_values(1) #get the time index from the Stoll list
    time_list_total= np.array(time_list_total.tolist())
    
    simultaneous_occurence= 0   
    for i_node in range(len(Stoll_node)):
        time_now= time_list_node[i_node]
        i_now_all= np.where(time_list_total== time_now)[0]
        occurence_now= len(i_now_all)

        if occurence_now > 1:
            lat_node, lon_node= Stoll_node['lat'].iloc[i_node], Stoll_node['lon'].iloc[i_node]
            
            lat_now_all, lon_now_all= Stoll.iloc[i_now_all].lat.values, Stoll.iloc[i_now_all].lon.values
    
            dist_to_others= distance((lat_node, lon_node), (lat_now_all, lon_now_all))
            dist_to_others= dist_to_others[dist_to_others != 0] #remove the distance to itself
            
            if np.min(dist_to_others) < 500: simultaneous_occurence +=1
            
            
    yplot -=.1
    variable= simultaneous_occurence/len(time_list_node) *100

    color= 'k'
#    if variable.median() < 0.0042: color= 'r'
#    if variable.median() > 0.0047: color= 'b'
#    ax1.text(0,yplot, 'N$_{500-1000}$: '+ str(np.round(variable.median(),4))+
#            ' ('+ str(np.round(np.nanpercentile(variable, 10),4)) + ', '+ str(np.round(np.nanpercentile(variable, 90),4)) +')', color=color   )
    ax1.text(0,yplot, 'Simultaneous occurrence: '+ str(int(variable))+'%' )

    
    
    
    
    






plt.tight_layout()


if save:
    if PLCG_type== 'stearing_flow': save_name=savedir+ 'SOM_'+Svar_full+"_"+str(Splevel)+"_"+Obs+'--other_vars_'+Obs_evol_str+"_x"+str(x)+"_y"+str(y)+'_msk'
    elif PLCG_type== 'track_smth': save_name=savedir+ 'SOM_track-smth'+smooth_param+'_'+Svar_full+"_"+str(Splevel)+"_"+Obs+'--other_vars_'+Obs_evol_str+"_x"+str(x)+"_y"+str(y)+'_msk'

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
    Stoll_cloud= Stoll.loc[Stoll['STARS lat'] > 1]
    Stoll_cloud= Stoll_cloud.loc[distance((Stoll_cloud.lat, Stoll_cloud.lon), (Stoll_cloud['STARS lat'], Stoll_cloud['STARS lon'])) < Sdist]
#    Stoll_cloud= remove_morphs(Stoll_cloud, 'Morphology', count= 25, how='remove', excludelist=['-', 'U'])

    #should maybe adapt this
    morphlist= list(Stoll_cloud.groupby('Morph_red').size().sort_values(ascending=False).index)[:7]
    if 'other' in morphlist: morphlist.sort(key = 'other'.__eq__) #this moves other to the end of the list




for iSOM in range(1, x*y +1):
    ax1= plt.subplot(y, x, iSOM)
    
    if Obs_evol != 'allObs':
        Morph_node= Stoll_ind[Stoll_ind.node== iSOM]['Morphology']
    if Obs_evol == 'allObs':
#        Morph_node= Stoll_cloud[Stoll_cloud.node== iSOM]['Morphology']
        Morph_node= Stoll_cloud[Stoll_cloud.node== iSOM]['Morph_red']


    for morph in morphlist:
#        plt.hist(Sx[Sx['Morphology'] == morph], bins= bins , label=morph, bottom= bot, edgecolor= 'k', lw= 1)
        if cloud_hist== 'number':
            ax1.bar(morph, len(Morph_node[Morph_node == morph]), edgecolor= 'k', lw= 1)
        
        if cloud_hist== 'percentage':
            ax1.bar(morph, len(Morph_node[Morph_node == morph])/len(Morph_node), edgecolor= 'k', lw= 1)
            ax1.set_ylim(0, 0.65)

    perc_of_SOM= len(np.where(df_nodes.node == iSOM)[0])/len(df_nodes.node)
    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')
#    ax1.axis('off')
#
#
#
plt.tight_layout()

if save:
    if PLCG_type== 'stearing_flow': save_name=savedir+ 'SOM_'+Svar_full+"_"+str(Splevel)+"_"+Obs+'--cloud_morphs_red'+str(x)+"_y"+str(y)
    elif PLCG_type== 'track_smth': save_name=savedir+ 'SOM_track-smth'+smooth_param+'_'+Svar_full+"_"+str(Splevel)+"_"+Obs+'--cloud_morphs_red'+str(x)+"_y"+str(y)
    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')

