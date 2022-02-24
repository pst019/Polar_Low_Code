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
#    homedir= '/home/'+user+'/home/'
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/PatsOrange/'

homedir=Mediadir+'home/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from f_useful import *
from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs
import xarray as xr


save= True
#save= False
savedir= homedir + 'Polar_Low/STARS-ana/Figs/SOM_geographic_composites/'

fignr= 9

imp_data=True
#imp_data=False

#maptype= 'PlateCarree'
maptype= 'Polar_Stereo'


typ='SOM'
typ='shear_cat'

x, y= 3, 3
lifelim= 6


node= 1

"""import Stoll list"""
Stoll_imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


Stoll= pd.read_csv(Stoll_imp_dir + 'Stoll_list_noRojo.csv')
Stoll= Stoll.set_index(['ID', 'time'])

S_nodes= pd.read_csv(Stoll_imp_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')
S_nodes= S_nodes.set_index(['ID', 'time'])

S_ERA5= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5.csv')
S_ERA5= S_ERA5.set_index(['ID', 'time'])

S_ERA5_2= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5_2.csv')
S_ERA5_2= S_ERA5_2.set_index(['ID', 'time'])

S_shear_cats= pd.read_csv(Stoll_imp_dir + 'shear_categories.csv')
S_shear_cats= S_shear_cats.set_index(['ID', 'time'])

Stoll= pd.concat([Stoll, S_nodes, S_ERA5, S_ERA5_2, S_shear_cats], axis=1)


Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps

Stoll= Stoll.reset_index()
Stoll['PLnr']= [int(ID.split('_')[0])*1E3 + int(ID.split('_')[1]) for ID in Stoll.ID.values]








"""get the ERA-5 data"""

plevel=None
vmin,vmax=None,None

"""specify field"""
#filetype='surface'
#var= '2t'
#var= '10u'
#var= '10v'
#var= 'U'
#var= 'skt'
#var='msl'

filetype='tracking'
plevel= 850
#var='U'
#var= 'u'
#var='Vort'
var= 't'
if plevel== 850: levels=  np.arange(250, 271, 2)
if plevel== 700: levels=  np.arange(244, 265, 2)

#var='grad_t'
#vmin,vmax= 0, 10

#var= 'uvec_mean

#var='skt-t'
#var='sst-t'
#vmin,vmax= 35, 52
#
#
#filetype='vorticity'
#plevel= 850
#var='vo'
#var='z'

#filetype='boundary'
#var= 'cape'
#var= 'blh'
#var= 'tcc'
#var= 'hcc'

#filetype='forecast'
##var='sshf'
##var='slhf'
#var= 'tp'
##var= 'sf'
##var= 'lsp'
##var= 'cp'
##var= 'lsp'
#var='sf'
#var= 'ttr' #toa_outgoing_longwave_flux


"""the contour variable"""
cont_var='z'
#cont_plev= 850


sym= False
#cmap='jet'
#cmap= 'viridis'
#cmap= 'Greys'
#cmap= 'Greys_r'
#cmap= 'Blues'
#cmap= 'Reds'

#sym= True
cmap= 'RdBu_r'


if typ == 'SOM': cat_list= np.arange(1, x*y +1)
elif typ == 'shear_cat': cat_list= np.arange(1, 6)

for node in cat_list:
    print('node', node)
    if typ == 'SOM': Stoll_node= Stoll[Stoll['node'] == node]
    elif typ == 'shear_cat': Stoll_node= Stoll[Stoll['shear_category'] == node]
    
    Stoll_node['date']=  pd.to_datetime(Stoll_node['time']).dt.strftime('%Y_%m_%d')
    
    date_list_different= remove_dublicate(Stoll_node['date'] )
    
    
    
    if imp_data:
        for i, date in enumerate(date_list_different):
            
            ds0= xr.open_dataset(Mediadir + "ERA5_STARS/data/" +filetype+'_era5_'+ date + '.nc')
            
            ds0= ds0.sel( time= np.array(pd.to_datetime(Stoll_node[Stoll_node['date'] == date].time) ) ) #select the correct hours
            ds0= ds0[[var, cont_var]]
        
            if plevel:
                ds0= ds0.sel(plev= plevel*100)
        
            if i == 0:
                ds= ds0
                
            else:
                ds= xr.concat([ds, ds0], 'time')
                
            if i%10 == 0: print(i)
    
    
    
    
    """plot the mean"""
    if maptype== 'Polar_Stereo': fig = plt.figure(fignr, figsize= (8,6) )
    elif maptype== 'PlateCarree': fig = plt.figure(fignr, figsize= (8,6) )
        
    plt.clf()
    
    if maptype== 'Polar_Stereo': ax= Plot_Polar_Stereo(fig, central_longitude= 20, extent= [-10, 50, 60, 80])
    elif maptype== 'PlateCarree': ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-40, 85, 50, 85])
    scale_bar(ax, 500, location= (0.06, 0.04))
    
    
    
     
    dsmean= ds[var].mean(dim='time')    
    dscontmean = ds[cont_var].mean(dim='time')  /9.81  
    
    cf= ax.contourf(ds.lon, ds.lat, dsmean, transform= ccrs.PlateCarree(), levels= levels, cmap= cmap, extend='both')
    cs= ax.contour(ds.lon, ds.lat, dscontmean, transform= ccrs.PlateCarree(), colors='k', linewidths= 1, levels= np.arange( 0, 5E3, 20))
    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')
    
    
    
    
    
    fig.tight_layout()
    
    if save: fig.savefig(savedir + var+'_'+cont_var+'_'+str(plevel)+'_'+ typ+ '_'+ str(node) , bbox_inches='tight')
    
    
    
  

"""with colorbar"""
cb= fig.colorbar(cf, ax= ax, shrink=0.7,
                     orientation= 'horizontal')
    
cb.set_label(ds[var].long_name + ' '+ str(plevel)+ 'hPa ['+ ds[var].units + ']', size=12)  

if save: fig.savefig(savedir + var+'_'+cont_var+'_'+str(plevel)+'_'+ str(node)+'_bar' , bbox_inches='tight')



