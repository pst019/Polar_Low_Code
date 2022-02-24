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
import matplotlib as mpl


plt.rcParams.update({'font.size': 18})


import numpy as np
from f_useful import *
from f_STARS import *

#
save= True
save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM_shear/'
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





#typ='Propagation'
typ='Shear'
    




"""import Stoll list"""
Stoll_imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


Stoll= pd.read_csv(Stoll_imp_dir + 'Stoll_list_noRojo.csv')
Stoll= Stoll.set_index(['ID', 'time'])

Stoll_dir_speed= pd.read_csv(Stoll_imp_dir + 'Stoll_list_dir-speed_smth1E-3.csv')
Stoll_dir_speed= Stoll_dir_speed.set_index(['ID', 'time'])

S_nodes= pd.read_csv(Stoll_imp_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')
S_nodes= S_nodes.set_index(['ID', 'time'])

S_ERA5= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5_1.csv')
S_ERA5= S_ERA5.set_index(['ID', 'time'])

S_ERA5_2= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5_2.csv')
S_ERA5_2= S_ERA5_2.set_index(['ID', 'time'])

Stoll= pd.concat([Stoll, Stoll_dir_speed, S_nodes, S_ERA5, S_ERA5_2], axis=1)


Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps






if typ== 'Propagation':
    speed_var= 'prop_speed'
    dir_var= 'prop_dir'
    bins= np.arange(0, 21, 5)
    title= "Speed [km$\cdot$h$^{-1}$]"
    typname= typ

elif typ== 'Shear':
    l1= 925
    l2= 700 #700
    mean= 250
    speed_var= 'vert_shear_strength_vec'+str(l1)+'-'+str(l2)+'_mean-'+str(mean)
#    speed_var= 'vert_shear_strength_vec'+str(l1)+'-500_mean-'+str(mean)
    
#    dir_var= 'vert_shear_angle_vec'+str(l1)+'-'+str(l2)+'_mean-'+str(mean)
    dir_var= 'vert_shear_angle_vec3_'+str(l1)+'-925-700-'+str(l2)+'_mean-'+str(mean)

#    dir_var= 'vert_shear_angle925-700_mean-400'

    Stoll[speed_var]*= 1E3
    Stoll[dir_var] = (Stoll[dir_var]+90)%360 #to rotate to the right
    bins= np.arange(0, 5, 1)
    title= "Shear strength\n[m$\cdot$s$^{-1}$ (km)$^{-1}$]"

    typname= typ+ str(l1)+'-'+str(l2)+'_mean-'+str(mean)



print('plot the SOMs')
#fig= plt.figure(fignr, figsize= (2.7*x,2.7*y) )
fig= plt.figure(fignr, figsize= (3,3)) # +0) )

#fignr+=1
plt.clf()



#good solution to put all in same: https://stackoverflow.com/questions/57027970/subplots-in-windrose-diagram

for iSOM in range(1, x*y +1):
    
    print('SOM ', iSOM)
    
    S_SOM= Stoll[Stoll.node== iSOM]

    fig= plt.figure(fignr, figsize= (3.5,3.5)) # +0) )
    plt.clf()
    
    if typ =='Propagation':
        ax = WindroseAxes.from_ax(fig= fig)
    elif typ== 'Shear':
        ax = WindroseAxes.from_ax(fig= fig)#, theta_labels= ["FS", "", "Left", "", "RS","",  "Right", ""]) #fig, rect)
#        ax.set_xticklabels(['N', '', 'E', '', 'S', '', 'W', ''])
        ax.set_xticklabels(["FS", "", "Left", "", "RS","",  "Right", ""])


    axb= ax.bar(S_SOM[dir_var], S_SOM[speed_var], normed=True, bins=bins,  opening=1 , edgecolor='k', cmap=mpl.cm.viridis_r)
    ax.set_rgrids([]) #remove the circular axis
    
    ax.set_title('SOM '+str(iSOM), position=(0.1, 1.01), color= 'r')


    if iSOM== 1:
        fig_legend= plt.figure(fignr+1, figsize= (10,3.5)) # +0) )
        plt.clf()

        ax2 = WindroseAxes.from_ax(fig= fig_legend)
        ax2.bar(S_SOM[dir_var], S_SOM[speed_var], normed=True, bins=bins,  opening=1 , edgecolor='k', cmap=mpl.cm.viridis_r)

        ax2.set_legend()
        ax2.legend(title=title, loc= (-1.1, 0), decimal_places=0)
#        fig_legend.tight_layout()
#        fig_legend= plt.figure(fignr+1, figsize= (3.5,3.5)) # +0) )
#        fig_legend.legend(axb)#, title="Speed [km$\cdot$h$^{-1}$]", decimal_places=0)
    
        if save:
            save_name=savedir+ typ+'_SOM-legend'
           
            print(save_name)
            fig_legend.savefig(save_name , bbox_inches='tight', dpi= 150)    
    

#
    if save:
        save_name=savedir+ typname+'_SOM-'+str(iSOM)
       
        print(save_name)
        fig.savefig(save_name , bbox_inches='tight', dpi= 150)
    

  





