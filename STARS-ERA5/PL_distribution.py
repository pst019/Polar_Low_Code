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
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/distr/'
#
fignr= 4

lifelim= 6

x=3 #3
y=3 #x+1

#
"""SOM var"""
imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


#Stoll= pd.read_csv(imp_dir + 'Stoll_list_noRojo.csv')
Stoll= pd.read_csv(imp_dir + 'Stoll_list.csv')
Stoll= Stoll.set_index(['ID', 'time'])

#S_nodes= pd.read_csv(imp_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')
#Stoll_nodes_full
S_nodes= pd.read_csv(imp_dir + 'Stoll_nodes_full_x'+str(x)+'_y'+str(y)+'.csv')

S_nodes= S_nodes.set_index(['ID', 'time'])

S_ERA5= pd.read_csv(imp_dir + 'Stoll_ERA5_1.csv')
S_ERA5= S_ERA5.set_index(['ID', 'time'])

S_ERA5_2= pd.read_csv(imp_dir + 'Stoll_ERA5_2.csv')
S_ERA5_2= S_ERA5_2.set_index(['ID', 'time'])

Stoll= pd.concat([Stoll, S_nodes, S_ERA5, S_ERA5_2], axis=1)


Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps



Smin, Smax, varlabel= '','',''


fig= plt.figure(fignr, figsize= (5,4) )
ax= fig.add_subplot(111)

fignr+=1
plt.clf()





""" chose variable"""
#var= 'vo'
#var= 'grad_t850_dgauss-4_mean-200'
#var= 'grad_t850_mean-250'

#var= 'baroclinic_dTdy_gr925-850-500_dgauss-4_msk_mean-200'
#var= 'grad_t850_dgauss-10_mean-1'
#if var== 'grad_t850_dgauss-10_mean-1': Stoll[var] *= 5

#var= 't850_mean-250'
#var= 'sst_msk_mean-200'
#if 'sst_' in var: Stoll[var] -= 273.15
#var= 'N925-500_dgauss-4_msk_mean-200'
var= 'N925-500_msk_mean-250'

#var= 'tp_msk_mean-250'
#var= 'vert_shear_angle925-700_mean-200'
#var= 'vert_shear_angle850-700_dgauss-10_mean-250'
#var= 'vert_shear_strength925-700_mean-200'
#
#var, varlabel, Smin, Smax= 'grad_t850_dgauss-10_mean-1', r" $\nabla$ T$_{850}$ [K $\cdot$ (500km)$^{-1}$]", 0, 4*5

#var= 'Heat_flux'

#var, varlabel= 'U_msk_max-250', 'Max U$_{10m}$ [m/s]'
#var= "Convrate"
#var= 'area'
#var= 'Cloud_diameter'
#var='Growth_rate'


"""other things"""
system_max= False
#system_max= True #calculate the maximum of the system



"""calculate variable"""
if var== 'area':
    Stoll[var]= np.sqrt(Stoll[var])/np.pi *4
    Smin, Smax, Swidth= 0, 500 , 25
    varlabel= 'Vortex diameter [km]'

if var== 'Cloud_diameter':
    Stoll[var]= Stoll['Cloud diameter']
    Smin, Smax, Swidth= 0, 800 , 50
    varlabel= 'Cloud diameter [km]'

if var== 'Growth_rate':
    Stoll[var]= 24*Stoll['vo_tendency']/Stoll['vo']
    varlabel= 'Growth rate [day$^{-1}$]'

#if var== 'grad_t850_dgauss-10_mean-1': Stoll[var]*= 5

if var== 'Heat_flux': Stoll[var]= Stoll['sshf_msk_mean-250']+ Stoll['slhf_msk_mean-250']
if var== "Convrate": Stoll[var]= Stoll['cp_msk_mean-250']/Stoll['tp_msk_mean-250']


"""make system variable"""
#for all time steps
Svar = Stoll[var].dropna()

if system_max:
    Svar= Svar.groupby('ID').max()


"""get the correct bounds"""
if Smin== '':
    Smin, Smax= Svar.min(), Svar.max()
    Swidth= (Smax - Smin)/10
    Smin= Smin -Swidth
    Smax= Smax +Swidth


x= np.linspace(Smin, Smax, 100)

#kernel= stats.gaussian_kde(Svar.dropna())
#plt.plot(xplot, kernel(xplot))

#for node in range(1,10):
#    kernel= stats.gaussian_kde(Stoll[Stoll.node==node][var].dropna())
#    plt.plot(xplot, kernel(xplot), label= str(node))


#nodelist= [[1,2], [3], [7], [9], [8], [5,6]]
#for node in nodelist:
kernel= stats.gaussian_kde(Svar.dropna())

#plt.hist(Svar, bins= np.arange(Smin, Smax + 20, 20) )
#plt.plot(x, kernel(x)* len(Svar.dropna()) *20)

#plt.hist(Svar, density=True, bins= np.arange(Smin, Smax + Swidth/2, Swidth)
#    , edgecolor='black', alpha= 0.6 )
plt.plot(x, kernel(x), lw= 3 , color= 'k' )

plt.plot(Svar.median(), np.max(kernel(x))*.02 , 'o', color='k')
print('Median', np.median(Svar) )

plt.plot(np.percentile(Svar, 10), np.max(kernel(x))*.02 , 'v', color='k')
print('10th percentile', np.percentile(Svar, 10 ) )

plt.plot(np.percentile(Svar, 90), np.max(kernel(x))*.02 , 'v', color='k')
print('90th percentile', np.percentile(Svar, 90 ) )

#if var== 'area':
#    plt.plot(Svar[Svar != 0].min(), np.max(kernel(x))*.02 , 's', color='k')
#else:    
#    plt.plot(Svar.min(), np.max(kernel(x))*.02 , 's', color='k')
#plt.plot(Svar.max(), np.max(kernel(x))*.02 , 'x', color='k')

plt.yticks([])
#import pylab
#ax.yaxis.set_major_locator(pylab.NullLocator())
#
#plt.legend()
#plt.xlabel(var)
plt.ylim(0, np.max(kernel(x)*1.05))
if var == 'grad_t850_mean-250': plt.xlim([0,5])
else: plt.xlim([x[0], x[-1]])

#
#    
##
#
#
#
#
if 'grad_t850' in var: varlabel= r"$|\nabla_h T_{850}|$ [K $\cdot$ (100 km)$^{-1}$]"
##elif 't850_' in var: varlabel= r"$T_{850}$ [K]"
##elif 'sst_' in var: varlabel= r"SST [$^{\circ}$C]"
##elif 'N925-500' in var: varlabel='N$_{500-925}$ [s$^{-1}$]'
##elif 'tp_' in var: varlabel= 'Total precip. [mm $\cdot$ h$^{-1}$]'
##elif 'vert_shear_angle' in var: varlabel= 'Vertical shear angle [$^{\circ}$]'
##elif 'vert_shear_strength' in var: varlabel= 'Vertical shear strength [s$^{-1}$]'
##elif 'U_msk' in var: varlabel= 'U$_{10m}$ [m$\cdot$s$^{-1}$]'
##elif 'Convrate' in var: varlabel= 'Conv. precip. rate [%]'
##elif 'Heat_flux' in var: varlabel= 'Surface heat flux [W$\cdot$m$^{-2}$]'

if varlabel=='': varlabel= var
plt.xlabel(varlabel)
#
plt.ylabel('Frequency')
plt.tight_layout()



if save:
    save_name=savedir+ 'distr_'+var
    if system_max: save_name += '_system_max'    

    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')
