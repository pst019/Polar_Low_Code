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
from windrose import WindroseAxes

import numpy as np
from f_useful import *
from f_STARS import *

#
save= True
#save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM4/'
#
fignr= 3


load_ds=True


"""SOM var"""
Obs= 'allObs_SOMprep' #the other versions were deleted

x=3 #3
y=3 #x+1

mirror_top_bottom=False
rotate_clock=False

if x== 3 and y == 4: mirror_top_bottom=True
if x== 3 and y == 3: rotate_clock=True

PLCG_type, smooth_param= 'track_smth', '1E-3'

lifelim= 6 #lifetime limit

vo_tend_excl= False
vo_tend_thresh= 0


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






 
    
    
"""import Stoll list"""
Stoll_imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


Stoll= pd.read_csv(Stoll_imp_dir + 'Stoll_list_noRojo.csv')
Stoll= Stoll.set_index(['ID', 'time'])

S_nodes= pd.read_csv(Stoll_imp_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')
S_nodes= S_nodes.set_index(['ID', 'time'])

Stoll= pd.concat([Stoll, S_nodes], axis=1)


Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps

Stoll= Stoll.reset_index()
Stoll['PLnr']= [int(ID.split('_')[0])*1E3 + int(ID.split('_')[1]) for ID in Stoll.ID.values]






file= Mediadir + 'ERA5_STARS/PL_centred_fields_smooth-tracks/ci_allObs_track-smth-'+smooth_param+'.nc'
ds= xr.open_dataset(file)


"""add the node to ds """
ds_node= np.zeros(len(ds.time))
ds_vo_tendency= np.zeros(len(ds.time))

for s in range(len(Stoll)):
    S= Stoll.iloc[s]
    ds_index= np.argwhere(np.logical_and(ds.PLnr.values== S.PLnr, ds.time.values== np.datetime64(S.time)))[0][0]
    ds_node[ds_index]= S.node
    ds_vo_tendency[ds_index]= S.vo_tendency


ds['node']= (('time'), ds_node)
ds['vo_tendency']= (('time'), ds_vo_tendency)



print('plot the SOMs')
#fig= plt.figure(fignr, figsize= (2.7*x,2.7*y) )
fig= plt.figure(fignr, figsize= (2.5*x,2.5*y +0) )

fignr+=1
plt.clf()

import matplotlib.gridspec as gridspec

gs = gridspec.GridSpec(x, y)     # (nblines, nbcol)
# return lists of bottom and top position of rows, left and right positions of columns.

bottom, top, left, right = gs.get_grid_positions(fig)


for iSOM in range(1, x*y +1):
    
    print('SOM ', iSOM)
    
#    fig= plt.figure(fignr-1)
#    ax= plt.subplot(y, x, iSOM)
    
    ds_SOM= ds.where(ds.node == iSOM, drop=True) #

#    col, row= (iSOM)%3, (iSOM)//3
#        print(mi, morph, col, row)
#    rect = [left[col],  bottom[row],  right[col]-left[col],  0.9*(top[row]-bottom[row])]

    
    ax = WindroseAxes.from_ax() #fig, rect)
#    fig.add_axes(ax2)
    
    wd= ds_SOM.beering.values
    ws= np.ones((len(wd)))
    ax.bar(wd, ws, normed=False,  opening=0.9 , edgecolor='white')

    ax.set_title('SOM '+str(iSOM), position=(0.1, 1.01), color= 'r', fontsize= 15)

#    ax.set_legend()
#    ds_node= ds.where(ds.node == iSOM, drop=True)
    


    if save:
        save_name=savedir+ 'Wind_dir_SOM-'+str(iSOM)
       
        print(save_name)
        plt.savefig(save_name , bbox_inches='tight')
    

  





