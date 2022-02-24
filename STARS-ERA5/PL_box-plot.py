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
#save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/distr/'
#
fignr= 4

lifelim= 6

x=3 #3
y=3 #x+1

#
"""SOM var"""
imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


Stoll= pd.read_csv(imp_dir + 'Stoll_list_noRojo.csv')
Stoll= Stoll.set_index(['ID', 'time'])

#S_nodes= pd.read_csv(imp_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')
#Stoll_nodes_full
S_nodes= pd.read_csv(imp_dir + 'Stoll_nodes_full_x'+str(x)+'_y'+str(y)+'.csv')

S_nodes= S_nodes.set_index(['ID', 'time'])

S_ERA5= pd.read_csv(imp_dir + 'Stoll_ERA5.csv')
S_ERA5= S_ERA5.set_index(['ID', 'time'])

S_ERA5_2= pd.read_csv(imp_dir + 'Stoll_ERA5_2.csv')
S_ERA5_2= S_ERA5_2.set_index(['ID', 'time'])

Stoll= pd.concat([Stoll, S_nodes, S_ERA5, S_ERA5_2], axis=1)


Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps






fig= plt.figure(fignr, figsize= (7,3.5) )

fignr+=1
plt.clf()

extra=''
#extra='reduce' #reduces the number of SOM nodes to display
extra='reduce_more'




#var= 'vo'
#var= 'vo_tendency'

#var= 'pres_pv_mean-250'
#var= 'grad_t850_dgauss-4_mean-200'
#var= 'grad_t850_mean-250'
#var= 'baroclinic_dTdy_gr925-850-500_dgauss-4_msk_mean-200'
#var= 'grad_t850_dgauss-10_mean-1'

#var= 't850_mean-250'
#var= 'sst_msk_mean-200'
#var= 'N925-500_dgauss-4_msk_mean-200'
#var= 'N925-500_msk_mean-250'

#var= 'tp_msk_mean-250'
#var= 'vert_shear_angle925-700_mean-200'
#var= 'vert_shear_angle850-700_dgauss-10_mean-250'
#var= 'vert_shear_strength925-700_mean-200'
#
#var= 'Heat_flux'

#var= 'U_msk_max-250'
var= "Convrate"

#if 'grad_t850' in var: Stoll[var] *= 5
if var== "Convrate": Stoll[var]= Stoll['cp_msk_mean-250']/Stoll['tp_msk_mean-250']
if 'sst_' in var: Stoll[var] -= 273.15
if var== 'Heat_flux': Stoll[var]= Stoll['sshf_msk_mean-250']+ Stoll['slhf_msk_mean-250']


x= np.linspace( Stoll[var].min(), Stoll[var].max(), 60)

#kernel= stats.gaussian_kde(Stoll[var].dropna())
#plt.plot(xplot, kernel(xplot))

#for node in range(1,10):
#    kernel= stats.gaussian_kde(Stoll[Stoll.node==node][var].dropna())
#    plt.plot(xplot, kernel(xplot), label= str(node))


nodelist= [[1,2], [3], [7], [9], [8], [5,6]]
for node in nodelist:
    kernel= stats.gaussian_kde(Stoll[Stoll.node.isin(node)][var].dropna())
    plt.plot(x, kernel(x), label= str(node))

#plt.hist(Stoll[var])
plt.legend()
plt.xlabel(var)
plt.xlim([x[0], x[-1]])
plt.ylabel('Frequency')

    
fig.tight_layout()



fig2= plt.figure(fignr, figsize= (6,3.5))
#fig2= plt.figure(fignr, figsize= (5,3.5))
plt.clf()
ax= plt.subplot()

"""reducing Stollnodes"""
if extra=='reduce':
    Stoll.loc[Stoll.node== 2, 'node']= 1
if extra=='reduce_more':
    Stoll.loc[Stoll.node== 2, 'node']= 1
    Stoll= Stoll[Stoll.node != 4]
    Stoll.loc[Stoll.node== 6, 'node']= 5
    Stoll.loc[Stoll.node== 8, 'node']= 5


Stoll= Stoll.dropna()
Stoll.boxplot(column = var, by = 'node', ax= ax, grid=False,
              color = {'whiskers' : 'black',
                            'caps' : 'black',
                            'medians' : 'black',
                            'boxes' : 'black'},
           whiskerprops = {'linewidth' : 2},
           capprops = {'linewidth' : 2},
                   flierprops = {'markersize': 1, 'marker':'o', 'markeredgecolor':'grey'},
                   medianprops = {'linewidth' : 2},
                   boxprops = {'linewidth' : 2})# , patch_artist = True)



print(Stoll.groupby(['node']).median()[var]) #print the median in each node

#plt.xticks(np.arange(1,10).astype(int))

if extra=='reduce':
    ax.set_xticklabels(['1+2', 3, 4, 5, 6, 7, 8, 9])
if extra=='reduce_more':
    ax.set_xticklabels(['1+2', 3, '5+6+8', 7, 9])

else:
    from matplotlib.ticker import FormatStrFormatter
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    
ax.tick_params(bottom= False)

ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.8)

ax.set_title('')
plt.suptitle('')

ax.set_xlabel('')
#ax.set_xlabel('SOM node')

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
ax.set_ylabel(varlabel)

fig2.tight_layout()



if save:
    save_name=savedir+ 'boxplot_'+var
    
    if extra=='reduce': save_name += '_red-SOMs'
    if extra=='reduce_more': save_name += '_red-SOMs-more'

    print(save_name)
    fig2.savefig(save_name , bbox_inches='tight')