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


save= True
#save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM/'

fignr= 14


"""get the ERA-5 data"""

plevel=850
var= 't_ano'
#var= 't'

#var='z'
#var='z_ano'

#var='vo'

#var='U'
#var='U_ano'

#var='q'
#var='q_ano'

x=5
y=6

#Obs= "Obsnr1"
#Obs= "Obsnr1_prep"
#Obs='mature'
Obs= 'allObs_SOMprep'

filedir= Mediadir + "ERA5_STARS/PL_centred_fields/SOM/"+var+"_"+str(plevel)+'_'+Obs+"_x"+str(x)+"_y"+str(y)
ds= xr.open_dataset(filedir+"_cluster.nc")

txt = pd.read_csv(filedir+"_node_nr.txt", sep=" ")
PLnr_vec= [str(txt.PLnr[n])[:-3]+'_'+str(txt.PLnr[n])[-1] for n in range(len(txt))]
"""the contour variable"""


sym= False
#cmap='jet'
#cmap= 'viridis'
#cmap= 'Greys'
#cmap= 'Greys_r'
#cmap= 'Blues'
#cmap= 'Reds'

if "ano" in var: sym= True
cmap= 'RdBu_r'



fig= plt.figure(fignr, figsize= (2.7*x,2.7*y) )
fignr+=1
plt.clf()

if sym: vextr= np.max([np.max(ds.field), -np.min(ds.field)])
else: vmax, vmin= float(np.max(ds.field)), float(np.min(ds.field))
    
for iSOM in range(1, x*y +1):
    ax1= plt.subplot(y, x, iSOM)
    if sym:
        cf= plt.contourf(ds.x, ds.y, ds.field.sel(SOM=iSOM), cmap= cmap, vmin= -vextr, vmax= vextr)
        cs= plt.contour(ds.x, ds.y, ds.field.sel(SOM=iSOM), colors='k', vmin= -vextr, vmax= vextr, linewidth= 1)
    else:
        cf= plt.contourf(ds.x, ds.y, ds.field.sel(SOM=iSOM), cmap= cmap, vmin= vmin, vmax= vmax)
        cs= plt.contour(ds.x, ds.y, ds.field.sel(SOM=iSOM), colors='k', vmin= vmin, vmax= vmax, linewidth= 1)

    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.1f')

    perc_of_SOM= len(np.where(txt.node == iSOM)[0])/len(txt.node)
    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')

#fig.subplots_adjust(right=0.87)
#cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
#cb= fig.colorbar(cf, cax=cbar_ax)

plt.tight_layout()

fig.subplots_adjust(bottom=0.1)
cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.015])
cb= fig.colorbar(cf, cax=cbar_ax, orientation="horizontal")
#cb.set_clim(vmin, vmax)

if var== 'z_ano': labelvar= 'Geopotential height anomaly [m]'
elif var== 'z': labelvar= 'Geopotential height [m]'

elif var== 't_ano': labelvar= 'Temperature anomaly [K]'
elif var== 't': labelvar= 'Temperature [K]'

elif var== 'U_ano': labelvar= 'Wind speed anomaly [m/s]'
elif var== 'U': labelvar= 'Wind speed [m/s]'

elif var== 'q_ano': labelvar= 'Specific humidity anomaly [g(kg)]'
elif var== 'q': labelvar= 'Specific humidity [g/kg]'
else: labelvar= var

cb.set_label(labelvar, size=14) 




if save:
    save_name=savedir+ 'SOM_'+var+"_"+str(plevel)+"_"+Obs+"_x"+str(x)+"_y"+str(y)
    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')
