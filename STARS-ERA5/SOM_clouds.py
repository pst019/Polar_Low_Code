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

#plt.rcParams.update({'font.size': 10})

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




"""reduced cloud distribution"""
fig= plt.figure(fignr, figsize= (6,2.7) )
fignr+=1
plt.clf()

Stoll_c2= Stoll_cloud
Stoll_c2.loc[Stoll_c2.node.isin([1,2,3,7,9]), 'node']= 1
Stoll_c2.loc[Stoll_c2.node.isin([5,6,8]), 'node']= 5

SOM_label_list=['1+2+3+7+9', '5+6+8']

for iS, iSOM in enumerate([1,5]):
    print(iSOM)
    ax1= plt.subplot(1, 2, iS+1)
    
    Morph_node= Stoll_c2[Stoll_c2.node== iSOM]['Morph_red']


    for morph in morphlist:
#        plt.hist(Sx[Sx['Morphology'] == morph], bins= bins , label=morph, bottom= bot, edgecolor= 'k', lw= 1)
        if cloud_hist== 'number':
            ax1.bar(morph, len(Morph_node[Morph_node == morph]), edgecolor= 'k', lw= 1)
        
        if cloud_hist== 'percentage':
            cloud_share= len(Morph_node[Morph_node == morph])/len(Morph_node)
            print(morph, cloud_share)
            ax1.bar(morph, cloud_share, edgecolor= 'k', lw= 1)
            ax1.set_ylim(0, 0.6)

    perc_of_SOM= len(np.where(Stoll_c2.node == iSOM)[0])/len(Stoll_c2.node[Stoll_c2.node > 0])
    plt.title('SOM '+SOM_label_list[iS]+ ' ('+str(int(np.round(perc_of_SOM*100, 0)))+'%)')
#    ax1.axis('off')
#
#
#
    if iS== 0: plt.ylabel('Frequency')    
    
plt.tight_layout()    



if save:
    save_name=savedir+ 'SOM_clouds_reduced'
    print(save_name)

    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')

