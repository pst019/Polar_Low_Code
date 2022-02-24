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

Stoll_dir_speed= pd.read_csv(Stoll_imp_dir + 'Stoll_list_dir-speed_smth1E-3.csv')
Stoll_dir_speed= Stoll_dir_speed.set_index(['ID', 'time'])

S_nodes= pd.read_csv(Stoll_imp_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')
S_nodes= S_nodes.set_index(['ID', 'time'])



Stoll= pd.concat([Stoll, Stoll_dir_speed, S_nodes], axis=1)


Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps





print('plot the SOMs')
#fig= plt.figure(fignr, figsize= (2.7*x,2.7*y) )
fig= plt.figure(fignr, figsize= (3,3)) # +0) )

#fignr+=1
plt.clf()



for iSOM in [1]: #range(1, x*y +1):
    
    print('SOM ', iSOM)
    
    S_SOM= Stoll[Stoll.node== iSOM]

    fig= plt.figure(fignr, figsize= (3.5,3.5)) # +0) )
    plt.clf()
    
    ax = WindroseAxes.from_ax(fig= fig) #fig, rect)
#    

    axb= ax.bar(S_SOM.prop_dir, S_SOM.prop_speed, normed=True, bins=np.arange(0, 21, 5),  opening=1 , edgecolor='k', cmap=mpl.cm.viridis_r)
    ax.set_rgrids([]) #remove the circular axis
    
    ax.set_title('SOM '+str(iSOM), position=(0.1, 1.01), color= 'r')


    if iSOM== 1:
        fig_legend= plt.figure(fignr+1, figsize= (8,3.5)) # +0) )
        plt.clf()

        ax2 = WindroseAxes.from_ax(fig= fig_legend)
        ax2.bar(S_SOM.prop_dir, S_SOM.prop_speed, normed=True, bins=np.arange(0, 21, 5),  opening=1 , edgecolor='k', cmap=mpl.cm.viridis_r)

        ax2.set_legend()
        ax2.legend(title="Speed [km$\cdot$h$^{-1}$]", loc= (-.8, 0), decimal_places=0)
#        fig_legend.tight_layout()
#        fig_legend= plt.figure(fignr+1, figsize= (3.5,3.5)) # +0) )
#        fig_legend.legend(axb)#, title="Speed [km$\cdot$h$^{-1}$]", decimal_places=0)
    
#        if save:
#            save_name=savedir+ 'Wind_speed-dir_SOM-legend'
#           
#            print(save_name)
#            fig_legend.savefig(save_name , bbox_inches='tight', dpi= 150)    
    

#
    if save:
        save_name=savedir+ 'Wind_speed-dir_SOM-'+str(iSOM)
       
        print(save_name)
        fig.savefig(save_name , bbox_inches='tight', dpi= 150)
    

  





