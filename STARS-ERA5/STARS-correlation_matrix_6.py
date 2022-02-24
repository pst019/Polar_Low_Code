#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 18:05:34 2019

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

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr, pearsonr
from f_useful import *
from f_STARS import *


plt.rcParams.update({'font.size': 13})

fignr= 4

save= False
#save= True



Stype='step'

list_type='Stoll'



imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


  
"""Stoll systems"""
Stoll= pd.read_csv(imp_dir + 'Stoll_list_noRojo.csv')
Stoll= Stoll.set_index(['ID', 'time'])

S_ERA5= pd.read_csv(imp_dir + 'Stoll_ERA5.csv')
S_ERA5= S_ERA5.set_index(['ID', 'time'])

S_PCs= pd.read_csv(imp_dir + 'Stoll_PCs.csv')
S_PCs= S_PCs.set_index(['ID', 'time'])

Stoll= pd.concat([Stoll, S_ERA5, S_PCs], axis=1)


Stoll= Stoll.loc[Stoll['Duration'] >= 5] #Duration gives actually the number of time steps


   
Stoll= Stoll.rename(columns={
    'lat': 'Latitude',
    'U_msk_max-200': 'Wind Speed 10m (max)',
    'vo': 'Vorticity$_{850}$ (centre)',
    'slp': 'Sea-level pressure (min)',
    'blh_msk_mean-200': 'Boundary layer height*',       
    'cape_msk_mean-200': 'CAPE*', 
    'sst_msk_mean-200': 'SST*',
   'N1000-500_msk_mean-200': 'N$_{500-1000}$ (mean)*',
   'N1000-850_msk_mean-200': 'N$_{850-1000}$ (mean)*',
   'N925-500_dgauss-4_msk_mean-200': 'N$_{500-925}$*$^{\dagger}$',
   'sst-t500_msk_mean-200': 'SST -T$_{500}$*',
   'sst-t700_msk_mean-200': 'SST -T$_{700}$*',
   'tp_msk_mean-200': 'Total precip.*',
    'cp_msk_mean-200': 'Convective precip.*',
    'lsp_msk_mean-200': 'Large-scale precip.*', 
    'sshf_msk_mean-200': 'Sensible heat flux*',
    'slhf_msk_mean-200': 'Latent heat flux*',

   'grad_t850_max-200': r"$\nabla_h$ T$_{850}$ (max)",
   'grad_t850_dgauss-4_msk_max-200': r"$\nabla_h$ T$_{850}$ (max)*",
   'grad_t850_dgauss-4_mean-200': r"$|\nabla_h T_{850}|$ $^{\dagger}$",
#       'baroclinic_dUdz_gr_filter4_925-700_max-200': 'Baroclinic dU/dz (max)',
#       'baroclinic_dTdy_gr_filter4_1000-850-500_max-200_msk': 'UL Baroclinic dT/dy (max)',
#       'baroclinic_dTdy_gr1000-925-700_dgauss-4_msk_max-200': 'LL Baroclinic dT/dy (max)',
   'baroclinic_dTdy_gr1000-925-850_dgauss-4_msk_mean-200': 'LL baroclinic (mean)*',
#       'baroclinic_dTdy_gr1000-850-500_dgauss-4_msk_mean-200':'Deep baroclinic (mean)*',
   'baroclinic_dTdy_gr925-850-500_dgauss-4_msk_mean-200':'Baroclinic growth*$^{\dagger}$',

#       'barotropic_gr_filter4_850_max-200': 'Barotropic growth (max)',
   'vert_shear_angle925-700_mean-200': 'Vertical shear angle',
   'vert_shear_strength925-700_mean-200': 'Vertical shear strength (mean)',
   })

    """drop some columns"""
#    Stoll= Stoll.drop(columns=['Comment', 'track file', 'Rojo nr', 'Rojo nr old', 'Press_min', 'U_10min_knots', 'row_idx', 'track_idx'])
   
    
    


        
        #the "full" matrix
#        Sx= Sx[ ['Morphology', 'Latitude', #'Cloud diameter', 
#                 'Vortex area', #'vortex_type',
#                 'Vorticity$_{850}$ (centre)', 'Wind Speed 10m (max)',
#               'Sea-level pressure (min)', 'Boundary layer height (mean)*', 'CAPE (mean)*',
##                'Skin temperature (med)',
## 'SKT -T$_{500}$ (max)', 'SKT -T$_{700}$ (max)', 'SST -T$_{500}$ (max)',
#                'SST (mean)*', 'SST -T$_{500}$ (mean)*', 'SST -T$_{700}$ (mean)*',
##                'N$_{500-1000}$ (min)', 'N$_{700-1000}$ (mean)',
#                'N$_{500-1000}$ (mean)*',
#                'N$_{850-1000}$ (mean)*',
#
#               'Total precip. (mean)*', 'Convective precip. (mean)*', 'Large-scale precip. (mean)*', 'Snow fall (mean)*',
#               'Sensible heat flux (mean)*', 'Latent heat flux (mean)*', 
##               'Grad T$_{850}$ (max)*',
#               r"$|\nabla_h T_{850}|$ $^{\dagger}$", #'Baroclinic growth (max)' , #'Barotropic growth (max)',
##               'Baroclinic dU/dz (max)' ,
##               'UL Baroclinic dT/dy (max)', 
#               'Deep baroclinic (mean)*',
##               'LL Baroclinic dT/dy (max)', 
#               'LL baroclinic (mean)*',
#               'Vertical shear angle (mean)', 'Vertical shear strength (mean)',
#               't_PC1', 't_PC2', 't_PC3'] ]

#some chosen parameters"
Sx0= Stoll[ ['Lifetime fraction',  'Latitude', 
         'Vorticity$_{850}$ (centre)', 'vo_tendency', 'Wind Speed 10m (max)',
         'Sensible heat flux*', 'Latent heat flux*', 
       'Total precip.*', 'Convective precip.*', 'Large-scale precip.*', #'Snow fall (mean)*',
       'Boundary layer height*', 'CAPE*',
        'SST*',  'SST -T$_{700}$*', 'SST -T$_{500}$*',
        'N$_{500-925}$*$^{\dagger}$',
      r"$|\nabla_h T_{850}|$ $^{\dagger}$", #'Baroclinic growth (max)' , #'Barotropic growth (max)',
       'Baroclinic growth*$^{\dagger}$',
       'Vertical shear angle',
       't_PC1', 't_PC2', 't_PC3'] ]

Sx0= Sx0.rename(columns={'vo_tendency': 'Vorticity tendency', 't_PC1': 'PC 1 (warm - cold)', 't_PC2': 'PC 2 (forward - reverse)', 't_PC3': 'PC 3 (warm - cold propag.)'})  

        
Sx0= Sx0.dropna() #remove all time steps where the land+ice mask is nan everywhere within the chosen radii

        
corr= Sx0.corr(method= 'spearman')





for n in range(len(corr)):
    corr.iloc[n, :n] = 0
    corr.iloc[n, n] = 0
    
corr= corr.drop(index=corr.index[-1])
corr= corr.drop(columns=corr.columns[0])

def discrete_cmap(base_cmap=None, N=11):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)


def discrete_cmap_even(base_cmap=None, N=10):
    """Create an N-bin discrete colormap from the specified input map
    N must be even"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = np.vstack((base(np.linspace(0, .5, N/2)), base(np.linspace(0.5, 1, N/2)) ))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)


corr = pd.melt(corr.reset_index(), id_vars='index') # Unpivot the dataframe, so we can get pair of arrays for x and y
corr.columns = ['x', 'y', 'value']
#


"""make the matrix"""
x=corr['x']
y=corr['y']
value=corr['value']
sortalph= False
xlabel='bottom'
half= True

""" sortalph= [True, False] - True -the columns should be ordered alphabetically, False - the column order of x, y is taken
xlabel= ['top', 'bootom'] - the position of the xticks and labels 
half- specifies if only the elements below the diagonal are plotted """

#    if xlabel =='top':
#        plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
#        plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True

#fig, ax = plt.subplots(num= fignr)
fig =plt.figure(fignr, figsize=(9,8))
plt.clf()

ax= fig.add_subplot(111)
N= 10 #number of levels in the correlation coeff color bar



#     Mapping from column names to integer coordinates
if sortalph:
    x_labels = [v for v in sorted(x.unique())]
    y_labels = [v for v in sorted(y.unique())][::-1]
else:
    x_labels = [v for v in x.unique()]
    y_labels = [v for v in y.unique()][::-1]
    
    
    
x_to_num = {p[1]:p[0] for p in enumerate(x_labels)} 
y_to_num = {p[1]:p[0] for p in enumerate(y_labels)} 

size_scale = 350
sc= ax.scatter(
    x=x.map(x_to_num), # Use mapping for x
    y=y.map(y_to_num), # Use mapping for y
    s=value.abs() * size_scale, # Vector of square sizes, proportional to size parameter
    c= value, vmin= -1, vmax= 1, 
#        cmap= plt.cm.get_cmap('RdYlGn'),
#        cmap= plt.cm.get_cmap('RdBu_r'),
    cmap= discrete_cmap_even('RdBu_r', N= N),
    marker='s' # Use square as scatterplot marker
)

# Show column labels on the axes


ax.set_xticks([x_to_num[v] for v in x_labels])
#    ax.set_xticks([x_to_num[v].split('(')[0] for v in x_labels])

if xlabel =='top':
    ax.tick_params(bottom=False, top=True, labelbottom=False, labeltop=True)
    ax.set_xticklabels(x_labels, rotation=45, horizontalalignment='left')
if xlabel =='bottom': ax.set_xticklabels(x_labels, rotation=90, horizontalalignment='center')
if xlabel =='both': #has to be developped
    x_labels_red= [x_lab.split(' (')[0] for x_lab in x_labels] #remove the '()'
    #bottom
    ax.set_xticklabels(x_labels_red, rotation=45, horizontalalignment='right')
    #top
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(ax.get_xticks())
#        ax2.set_xticklabels(x_labels, rotation=45, horizontalalignment='left')
    ax2.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='left')
#        ax2.xaxis.set_major_formatter(ax.xaxis.get_major_formatter)

ax.set_yticks([y_to_num[v] for v in y_labels])
ax.set_yticklabels(y_labels)
#ax.get_yaxis().set_visible(False)
ax.set_frame_on(False)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')

#cb_ax = fig.add_axes([0.04, 0.2, 0.25, 0.02])
cb_ax = fig.add_axes([0.65, 0.92, 0.25, 0.02])

#plt.colorbar(sc, ticks= np.linspace(-1, 1, N+1), label= 'Correlation coeffient', shrink= 0.7)


plt.colorbar(sc, ticks= np.linspace(-1, 1, N+1)[::2], label= 'Correlation coeffient', orientation = 'horizontal',
             shrink= 0.7, cax=cb_ax)


















plt.tight_layout()
fignr+=1


if save:
#    if Stype== 'system': Stype += '_'+system_char
    
    savedir= homedir+ 'Polar_Low/STARS-ana/Figs/2/'+'Correlation-Matrix_'+Stype +'_msk2'
    print(savedir)
    plt.savefig(savedir, bbox_inches='tight')





