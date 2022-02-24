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
save= True
#save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM4/'
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


"""evolve params"""
Obs_evol='allObs'






"""import Stoll list"""
imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


Stoll= pd.read_csv(imp_dir + 'Stoll_list.csv')
Stoll= Stoll.set_index(['ID', 'time'])

S_ERA5= pd.read_csv(imp_dir + 'Stoll_ERA5.csv')
S_ERA5= S_ERA5.set_index(['ID', 'time'])

S_nodes= pd.read_csv(imp_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')
S_nodes= S_nodes.set_index(['ID', 'time'])

Stoll= pd.concat([Stoll, S_ERA5, S_nodes], axis=1)


Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps



#vo_tend_thresh= 0
#Stoll= Stoll[Stoll.vo_tendency > vo_tend_thresh]


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
   'grad_t850_dgauss-4_mean-200': 'Grad T$_{850}$ (mean)',
#       'baroclinic_dUdz_gr_filter4_925-700_max-200': 'Baroclinic dU/dz (max)',
#       'baroclinic_dTdy_gr_filter4_1000-850-500_max-200_msk': 'UL Baroclinic dT/dy (max)',
#       'baroclinic_dTdy_gr1000-925-700_dgauss-4_msk_max-200': 'LL Baroclinic dT/dy (max)',
   'baroclinic_dTdy_gr1000-925-850_dgauss-4_msk_mean-200': 'LL baroclinic (mean)*',
#       'baroclinic_dTdy_gr1000-850-500_dgauss-4_msk_mean-200': 'Deep baroclinic (mean)*',
   'baroclinic_dTdy_gr925-850-500_dgauss-4_msk_mean-200': 'Deep baroclinic (mean)*',

#       'barotropic_gr_filter4_850_max-200': 'Barotropic growth (max)',
   'vert_shear_angle925-700_mean-400': 'Vertical shear angle (mean)',
#   'vert_shear_angle850-700_mean-500': 'Vertical shear angle (mean)',

   'vert_shear_strength925-700_mean-200': 'Vertical shear strength (mean)',
   })









"""the other parameters"""
fig= plt.figure(fignr, figsize= (2.7*x,2.7*y) )
fignr+=1
plt.clf()





for iSOM in range(1, x*y +1):
    ax1= plt.subplot(y, x, iSOM)

    
    Stoll_node= Stoll[Stoll.node== iSOM]


    perc_of_SOM= len(np.where(Stoll.node == iSOM)[0])/len(Stoll.node[Stoll.node > 0])
    ax1.set_title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')
    ax1.axis('off')
        
    yplot= 0.9
    
    variable=Stoll_node['Vertical shear angle (mean)']
    varall= Stoll['Vertical shear angle (mean)']
    color= 'k'
    if variable.median() < np.nanpercentile(varall, 30): color= 'b'
    if variable.median() > np.nanpercentile(varall, 70): color= 'r'
    ax1.text(0,yplot, 'Shear angle: '+ str(int(np.round(variable.median())))+
            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')', color=color  )

#    yplot -=.1
#    variable=Stoll_node['vert_wind_grad925-700_mean_200']*1E3
#    color= 'k'
#    if variable.median() > 2: color= 'r'
#    if variable.median() < -2: color= 'b'     
#    ax1.text(0,yplot, 'Shear speed: '+ str(np.round(variable.median(),1))+
#            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')' , color=color )


#    yplot -=.1
#    variable=Stoll_node['Deep baroclinic (mean)*']*24
#    varall=Stoll['Deep baroclinic (mean)*']*24
#
#    color= 'k'
#    if variable.median() > np.nanpercentile(varall, 70): color= 'r'
#    if variable.median() < np.nanpercentile(varall, 30): color= 'b'
#    ax1.text(0,yplot , 'Baroclinic growth: '+ str(np.round(variable.median(),1))+
#            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')' , color=color  )

#    yplot -=.1
#    variable=Stoll_node['LL baroclinic (mean)*']*24
#    color= 'k'
#    if variable.median() > 4: color= 'r'
#    if variable.median() < 2.5: color= 'b'
#    ax1.text(0,yplot , 'LL baroclinic gr: '+ str(np.round(variable.median(),1))+
#            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')' , color=color  )


#    yplot -=.1
#    variable=Stoll_node['Grad T$_{850}$ (mean)']
#    varall=Stoll['Grad T$_{850}$ (mean)']
#
#    color= 'k'
#    if variable.median() > np.nanpercentile(varall, 70): color= 'r'
#    if variable.median() < np.nanpercentile(varall, 30): color= 'b'
#    ax1.text(0,yplot, r"$ \nabla$ T$_{850}$: "+ str(np.round(variable.median(),1))+
#            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')', color=color   )


    yplot -=.1
    variable=Stoll_node['grad_t850_dgauss-10_mean-1']
    varall=Stoll['grad_t850_dgauss-10_mean-1']

    color= 'k'
    if variable.median() > np.nanpercentile(varall, 70): color= 'r'
    if variable.median() < np.nanpercentile(varall, 30): color= 'b'
    ax1.text(0,yplot, r"$ \nabla$ T$_{850}$: "+ str(np.round(variable.median(),1))+
            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')', color=color   )



    yplot -=.1
    variable=Stoll_node['N$_{500-925}$ (mean)*']
    varall=Stoll['N$_{500-925}$ (mean)*']

    color= 'k'
    if variable.median() < np.nanpercentile(varall, 30): color= 'r'
    if variable.median() > np.nanpercentile(varall, 70): color= 'b'
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


#    yplot -=.1
#    variable=Stoll_node['SST -T$_{500}$ (mean)*']
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
#    variable=Stoll_node['CAPE (mean)*']
#    varall=Stoll['CAPE (mean)*']
#
#    color= 'k'
#    if variable.median() > np.nanpercentile(varall, 70): color= 'r'
#    if variable.median() < np.nanpercentile(varall, 30): color= 'b'   
#    ax1.text(0,yplot, 'CAPE: '+ str(int(np.round(variable.median())))+
#            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')' , color=color )

    yplot -=.1
    variable=Stoll_node['Total precip. (mean)*']
    varall=Stoll['Total precip. (mean)*']

    color= 'k'
    if variable.median() > np.nanpercentile(varall, 70): color= 'r'
    if variable.median() < np.nanpercentile(varall, 30): color= 'b'   
    ax1.text(0,yplot, 'Precipitation: '+ str(np.round(variable.median(),2))+
            ' ('+ str(np.round(np.nanpercentile(variable, 10),2)) + ', '+ str(np.round(np.nanpercentile(variable, 90),2)) +')', color=color   )


    yplot -=.1
    variable=100*Stoll_node['Convective precip. (mean)*']/Stoll_node['Total precip. (mean)*']
    varall=100*Stoll['Convective precip. (mean)*']/Stoll['Total precip. (mean)*']

    color= 'k'
    if variable.median() > np.nanpercentile(varall, 70): color= 'r'
    if variable.median() < np.nanpercentile(varall, 30): color= 'b'   
    ax1.text(0,yplot, 'Conv precip rate: '+ str(int(np.round(variable.median())))+
            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')' , color=color )


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
    varall=Stoll['Sensible heat flux (mean)*'] + Stoll['Latent heat flux (mean)*']

    color= 'k'
    if variable.median() > np.nanpercentile(varall, 70): color= 'r'
    if variable.median() < np.nanpercentile(varall, 30): color= 'b'     
    ax1.text(0,yplot, 'Surface fluxes: '+ str(int(np.round(variable.median())))+
            ' ('+ str(int(np.round(np.nanpercentile(variable, 10)))) + ', '+ str(int(np.round(np.nanpercentile(variable, 90)))) +')' , color=color )



    yplot -=.1
    variable=Stoll_node['Wind Speed 10m (max)']
    varall=Stoll['Wind Speed 10m (max)']

    color= 'k'
    if variable.median() > np.nanpercentile(varall, 70): color= 'r'
    if variable.median() < np.nanpercentile(varall, 30): color= 'b' 
    ax1.text(0,yplot , 'U$_{10m}$: '+ str(np.round(variable.median(),1))+
            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')'  , color=color )



#    yplot -=.1
#    variable=Stoll_node['Skin temperature (med)']- 273.15
#    color= 'k'
#    if variable.median() > 6: color= 'r'
#    if variable.median() < 4: color= 'b'     
#    ax1.text(0,yplot, 'Skin temp: '+ str(np.round(variable.median(),1))+
#            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')' , color=color )


    yplot -=.1
    variable=Stoll_node['t850_mean-200']- 273.15
    varall=Stoll['t850_mean-200']- 273.15
    
    color= 'k'
    if variable.median() > np.nanpercentile(varall, 70): color= 'r'
    if variable.median() < np.nanpercentile(varall, 30): color= 'b'
    ax1.text(0,yplot, r"T$_{850}$: "+ str(np.round(variable.median(),1))+
            ' ('+ str(np.round(np.nanpercentile(variable, 10),1)) + ', '+ str(np.round(np.nanpercentile(variable, 90),1)) +')' , color=color )




    yplot -=.1
    variable=Stoll_node['SST (mean)*']- 273.15
    varall=Stoll['SST (mean)*']- 273.15
    
    color= 'k'
    if variable.median() > np.nanpercentile(varall, 70): color= 'r'
    if variable.median() < np.nanpercentile(varall, 30): color= 'b'
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
    save_name=savedir+ 'SOM_other_parameters_x'+str(x)+"_y"+str(y)

    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')




"""cloud distributions"""
fig= plt.figure(fignr, figsize= (2.7*x,2.7*y) )
fignr+=1
plt.clf()


Sdist= 100
#only the steps where the cloud morphology is observed and the distance is lower than Sdist
Stoll_cloud= Stoll.loc[Stoll['Rojo_lat'] > 1]
Stoll_cloud= Stoll_cloud.loc[distance((Stoll_cloud.lat, Stoll_cloud.lon), (Stoll_cloud['Rojo_lat'], Stoll_cloud['Rojo_lon'])) < Sdist]
#    Stoll_cloud= remove_morphs(Stoll_cloud, 'Morphology', count= 25, how='remove', excludelist=['-', 'U'])

#should maybe adapt this
morphlist= list(Stoll_cloud.groupby('Morph_red').size().sort_values(ascending=False).index)[:7]
if 'other' in morphlist: morphlist.sort(key = 'other'.__eq__) #this moves other to the end of the list




for iSOM in range(1, x*y +1):
    ax1= plt.subplot(y, x, iSOM)
    
    Morph_node= Stoll_cloud[Stoll_cloud.node== iSOM]['Morph_red']


    for morph in morphlist:
#        plt.hist(Sx[Sx['Morphology'] == morph], bins= bins , label=morph, bottom= bot, edgecolor= 'k', lw= 1)
        if cloud_hist== 'number':
            ax1.bar(morph, len(Morph_node[Morph_node == morph]), edgecolor= 'k', lw= 1)
        
        if cloud_hist== 'percentage':
            ax1.bar(morph, len(Morph_node[Morph_node == morph])/len(Morph_node), edgecolor= 'k', lw= 1)
            ax1.set_ylim(0, 0.65)

    perc_of_SOM= len(np.where(Stoll.node == iSOM)[0])/len(Stoll.node[Stoll.node > 0])
    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')
#    ax1.axis('off')
#
#
#
plt.tight_layout()

if save:
    save_name=savedir+ 'SOM_clouds_x'+str(x)+"_y"+str(y)
    print(save_name)

    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')

