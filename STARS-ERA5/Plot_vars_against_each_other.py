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
    Mediadir= '/run/media/pst019/PatsOrange/'

homedir=Mediadir+'home/'


import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import xarray as xr
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import stats

plt.rcParams.update({'font.size': 15})

import numpy as np
from f_useful import *
from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs

#
save= True
save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/new_class_3/'
#
fignr= 4

lifelim= 6

#x=3 #3
#y=3 #x+1

#
"""SOM var"""
imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


Stoll= pd.read_csv(imp_dir + 'Stoll_list_noRojo.csv')
Stoll= Stoll.set_index(['ID', 'time'])

#S_nodes= pd.read_csv(imp_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')
#Stoll_nodes_full
#S_nodes= pd.read_csv(imp_dir + 'Stoll_nodes_full_x'+str(x)+'_y'+str(y)+'.csv')
#S_nodes= S_nodes.set_index(['ID', 'time'])

S_shear_cats= pd.read_csv(imp_dir + 'shear_categories.csv')
S_shear_cats= S_shear_cats.set_index(['ID', 'time'])

S_ERA5= pd.read_csv(imp_dir + 'Stoll_ERA5.csv')
S_ERA5= S_ERA5.set_index(['ID', 'time'])

S_ERA5_2= pd.read_csv(imp_dir + 'Stoll_ERA5_2.csv')
S_ERA5_2= S_ERA5_2.set_index(['ID', 'time'])

Stoll= pd.concat([Stoll, S_shear_cats, S_ERA5, S_ERA5_2], axis=1)


Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps






fig= plt.figure(fignr, figsize= (7,3.5) )

fignr+=1
plt.clf()

#extra=''
#extra='reduce' #reduces the number of SOM nodes to display
#extra='reduce_more'




#var= 'vo'
var= 'vo_tendency'

#var= 'pres_pv_mean-250'
#var= 'grad_t850_dgauss-4_mean-200'
#var= 'grad_t850_mean-250'
#var= 'baroclinic_dTdy_gr925-850-500_dgauss-4_msk_mean-200'
#var= 'grad_t850_dgauss-10_mean-1'

#var= 'cape_msk_mean-200'
#var= 'sst-t500_msk_mean-200'
#var= 'slhf_msk_mean-200'
#var= 'sshf_msk_mean-200'
#var=  'blh_msk_mean-200'
#var= 't850_mean-250'
#var= 'sst_msk_mean-200'
#var= 'N925-500_dgauss-4_msk_mean-200'
#var= 'N925-500_msk_mean-250'
#var= 'Lifetime fraction'
#var= 'tp_msk_mean-250'
#var= 'vert_shear_angle925-700_mean-200'
#var= 'vert_shear_angle850-700_dgauss-10_mean-250'
#var2= 'vert_shear_strength925-700_mean-200'
var2= 'vert_shear_strength_vec925-500_mean-250'

#var= 'Heat_flux'

#var= 'U_msk_max-250'
#var= "Convrate"

#if 'grad_t850' in var: Stoll[var] *= 5
if var== "Convrate": Stoll[var]= Stoll['cp_msk_mean-250']/Stoll['tp_msk_mean-250']
if 'sst_' in var: Stoll[var] -= 273.15
if var== 'Heat_flux': Stoll[var]= Stoll['sshf_msk_mean-250']+ Stoll['slhf_msk_mean-250']




plt.scatter(Stoll[var], Stoll[var2])
print('Corr: ', np.corrcoef(Stoll[var], Stoll[var2])[0,1])

"""reducing Stollnodes"""


if 'grad_t850' in var: varlabel= r"$|\nabla_h T_{850}|$ [K $\cdot$ (100 km)$^{-1}$]"
elif 't850_' in var: varlabel= r"$T_{850}$ [K]"
elif 'sst_' in var: varlabel= r"SST [$^{\circ}$C]"
elif 'N925-500' in var: varlabel='N$_{500-925}$ [s$^{-1}$]'
elif 'tp_' in var: varlabel= 'Total precip. [mm $\cdot$ h$^{-1}$]'
elif 'vert_shear_angle' in var: varlabel= 'Vertical shear angle [$^{\circ}$]'
elif 'vert_shear_strength' in var: varlabel= 'Vertical shear strength [s$^{-1}$]'
elif 'U_msk' in var: varlabel= 'U$_{10m}$ [m$\cdot$s$^{-1}$]'
elif 'Convrate' in var: varlabel= 'Conv. precip. rate'
elif 'Heat_flux' in var: varlabel= 'Surface heat flux [W$\cdot$m$^{-2}$]'
else: varlabel= var
plt.xlabel(varlabel)
plt.ylabel(var2)

plt.tight_layout()



if save:
    save_name=savedir+ 'Var-'+var+'_against'+var2
    
#    if extra=='reduce': save_name += '_red-SOMs'
#    if extra=='reduce_more': save_name += '_red-SOMs-more'

    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')
