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

plt.rcParams.update({'font.size': 13})

import numpy as np
from f_useful import *
from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs

#
save= True
save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/evol/'
#
fignr= 1

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

Stoll= pd.concat([Stoll, S_nodes, S_ERA5], axis=1)


Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps






fig= plt.figure(fignr, figsize= (7,3.5) )

fignr+=1
plt.clf()

from scipy.interpolate import interp1d

#var= 'vo'
#var= 'grad_t850_dgauss-4_mean-200'
var= 'baroclinic_dTdy_gr925-850-500_dgauss-4_msk_mean-200'
#var= 'N925-500_dgauss-4_msk_mean-200'
#var= 'tp_msk_mean-200'

sel= 'node_evol_hopp_rem_nan_rm'

t_interp= np.linspace(0,1, 21)
PLnr= remove_dublicate(Stoll.index.get_level_values(0))

ds = xr.Dataset( coords={'PLnr': (['PLnr'], PLnr), 'time_fraction': (['tf'], t_interp)})
ds[var]= (('PLnr', 'tf'),  np.zeros((len(PLnr), len(t_interp))) )
ds[sel]= (('PLnr'), Stoll[sel].groupby('ID').first() )

for ID in PLnr:

    S_now= Stoll.xs(ID)
    times= pd.to_datetime(S_now.index)
    dt= (times- times[0])/np.timedelta64(1,'h')
    Lifetime_frac= dt/dt[-1]
    
    ds[var].loc[ID]= interp1d(S_now['Lifetime fraction'], S_now[var])(t_interp)


if var== 'vo': ds[var]*= 1E5


color = 'tab:blue'

plt.plot(t_interp, ds[var].mean(dim='PLnr'), color= color) 
plt.plot(t_interp, np.nanpercentile(ds[var].values, q= 10, axis= 0), '--', color= color)
plt.plot(t_interp, np.nanpercentile(ds[var].values, q= 90, axis= 0), '--', color= color)


"""plot different variable evolutions"""
for isel in ['[1]', '[9]', '[9, 8]', '[9, 6]', '[3]', '[7]', '[7, 4]']:
    plt.plot(t_interp, np.nanmean(ds.where(ds[sel] == isel, drop=True)[var], axis= 0), label= isel) 

plt.legend()


plt.xlabel('Lifetime fraction')
plt.xlim(0,1)


if var== 'vo': labelvar= 'Relative vorticity [10$^{-5}$ s$^{-1}$]'

else: labelvar= var

#plt.ylabel(labelvar)

    
fig.tight_layout()

if save:
    save_name=savedir+ 'evol_'+var

    print(save_name)
    fig.savefig(save_name , bbox_inches='tight')
