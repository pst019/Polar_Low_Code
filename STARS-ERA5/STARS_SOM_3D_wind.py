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

homedir=Mediadir+'home/'


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
save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM2/'
#
fignr= 7

imp_list=True
#imp_list=False

imp_ds=True
#imp_ds=False


"""SOM var"""
Obs= 'allObs_SOMprep' #the other versions were deleted

x=3 #3
y=3 #x+1

mirror_top_bottom=False
rotate_clock=False


SOM=9


if x== 3 and y == 4: mirror_top_bottom=True
if x== 3 and y == 3: rotate_clock=True

PLCG_type, smooth_param= 'track_smth', '1E-3'

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


"""evolve var"""
#Obs_evol=1
#Obs_evol= 'mature'
#mat_evol_step= 100 #in percent how much of time the PL has overcome towards the mature phase (0- initial stage, 100- mature stage)
#mat_evol_step= 25
Obs_evol='allObs'

if Obs_evol== 'mature': Obs_evol_str= Obs_evol +str(mat_evol_step)
elif Obs_evol =='allObs': Obs_evol_str= Obs_evol
else: Obs_evol_str= 'Obsnr'+str(Obs_evol)


"""shading variable"""
#ano_var=True
ano_var=False

imp_vars= ['u', 'v', 'w']
#var='t'
#plevel_var= [1000, 850, 700, 500]
plevel_var= [1000, 925, 850, 700]







 
    
    
"""import Stoll list"""
if imp_list:
#    SOM_filedir= Mediadir + "ERA5_STARS/PL_centred_fields/SOM/"+Svar_full+"_"+str(Splevel)+'_'+Obs+"_x"+str(x)+"_y"+str(y)

    SOM_filedir= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/SOM/"+Svar_full+"_"+str(Splevel)+'_'+Obs+ '_track-smth-'+smooth_param+'_dur'+str(lifelim)+"_x"+ str(x)+"_y"+str(y)

    ds_init= xr.open_dataset(SOM_filedir+"_cluster.nc")

    Stoll, Stoll_ind= imp_standard_Stoll() #SOM_filedir)#, Obs_evol)
    df_nodes= pd.read_csv(SOM_filedir+"_node_nr.txt", sep=" ")

    
  
    df_nodes['time']= pd.to_datetime(df_nodes.date.apply(str), format='%Y%m%d%H')
    df_nodes['ID']= [str(int(df_nodes.PLnr[n]))[:-3]+'_'+str(int(df_nodes.PLnr[n])%1000) for n in range(len(df_nodes.PLnr))]
#        df_nodes= df_nodes.drop(columns=['date', 'PLnr'])
    Stoll= Stoll.set_index(['ID', 'time'])
    df_nodes= df_nodes.set_index(['ID', 'time'])       
#        Stoll= pd.concat([Stoll, df_nodes], axis= 1)

    """to rotate/mirror the nodes"""
    df_nodes_orig= df_nodes.copy()
    if mirror_top_bottom:
        df_nodes['node']= (y-1 -(df_nodes.node.values-1)//x)*x + (df_nodes.node.values -1)%x +1
    if rotate_clock:
        x_n= (df_nodes['node'].values-1)%x +1
        y_n= (df_nodes['node'].values-1)//x +1
        y_n= y+1 - y_n #reverse the y_nodes
        
        df_nodes['node']= (x_n -1)*x + y_n








if imp_ds:
    print('import from all time step list')

    """file for shading variable"""
    for ni, i_var in enumerate(imp_vars):
        file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +i_var +'_allObs_track-smth-'+smooth_param+'.nc'
        if os.path.isfile(file):
            ds1= xr.open_dataset(file)
        else:
            ds1= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +i_var + '_all_levs_allObs_track-smth-'+smooth_param+'.nc')
            ds1= ds1.sel(plev= plevel_var)
            
        if ni == 0: ds= ds1
        else: ds= xr.merge([ds, ds1])     
    



    
    """add the node to ds """
    ds_node_var= np.zeros(len(ds.time))
    
    for s in range(len(df_nodes)):
        df_nodes_s= df_nodes.iloc[s]
        ds_index= np.argwhere(np.logical_and(ds.PLnr.values== df_nodes_s.PLnr, ds.time.values== np.datetime64(df_nodes_s.name[1])))[0][0]
        ds_node_var[ds_index]= df_nodes_s.node
    
    ds['node']= (('time'), ds_node_var)
    
    
    """select only the specific SOM node"""
    ds= ds.where(ds.node == SOM, drop=True) #
    
 
    
    """compute"""
    U= np.sqrt(ds.u**2 + ds.v**2)
    wind_beering= UV2Direction(ds.u, ds.v)
    track_beering = ds.beering.values
    
    rot_beering= (track_beering[:, np.newaxis, np.newaxis, np.newaxis]- wind_beering)%360
    v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='') #maybe u_zonal and v_meridional are exchanged in WindSpeedDirection2UV, but this seem to work"""

    ds['u_r']= (('time', 'plev', 'x', 'y'), u_r )
    ds['u_r'].attrs['units'] = 'm/s'
    ds['u_r'].attrs['long_name']='Tangential wind'

    
    ds['v_r']= (('time', 'plev', 'x', 'y'), v_r )
    ds['v_r'].attrs['units'] = 'm/s'
    ds['v_r'].attrs['long_name']='Azimuthal wind'
    




"""plot 3D"""
fig= plt.figure(fignr) #, figsize= (2.5*x,2.5*y +0) )

fignr+=1
plt.clf()
    
from mpl_toolkits.mplot3d import Axes3D    

ax = fig.add_subplot(111, projection='3d')

ds= ds.transpose('time', 'x', 'y', 'plev')

en= 2 #every nth

u= np.mean(ds['u_r'], axis= 0).values[::en, ::en] *3.6
v= np.mean(ds['v_r'], axis= 0).values[::en, ::en] *3.6
w= np.mean(ds['w'], axis= 0).values[::en, ::en]  *36 #3600/100


x,y,z= np.meshgrid(ds.x.values[::en], ds.y.values[::en], plevel_var)
ax.quiver(x, y, z, u, v, w, length=1 )#, normalize=True)

#cs = ax.contour(ds.x.values, ds.y.values, np.mean(ds[cont2_full], axis= 0).values , 
#               offset= 500, levels= cont_levels2, colors='g', zorder= 100) 

ax.set_xlim(-300, 300)
ax.set_ylim(-300, 300)
ax.set_zlim(1010, 650)
 
#, zdir='z', offset=np.min(Z), cmap=cm.ocean)




# contourf(self, X, Y, Z, *args, zdir='z', offset=None, **kwargs)


#"""make the arrow"""
#arrow_ax = fig.add_axes([0.75, 0.03, 0.2, 0.04])
#import matplotlib.patches as mpatches
#
#arrow = mpatches.FancyArrowPatch((0.2, 0.7), (.7, 0.7), mutation_scale=20)
#arrow_ax.add_patch(arrow)
#arrow_ax.text(0,0 , 'Propagation direction', fontweight= 'bold')#, fontsize=13)
#arrow_ax.set_frame_on(False)
#arrow_ax.get_xaxis().set_visible(False)
#arrow_ax.get_yaxis().set_visible(False)





