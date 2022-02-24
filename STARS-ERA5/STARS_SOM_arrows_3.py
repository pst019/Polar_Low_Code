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


plt.rcParams.update({'font.size': 13})

#
save= True
save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM5/'
#
fignr= 6


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





dist=250
 
    
    
"""import Stoll list"""
Stoll_imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


Stoll= pd.read_csv(Stoll_imp_dir + 'Stoll_list_noRojo.csv')
Stoll= Stoll.set_index(['ID', 'time'])

Stoll_dir_speed= pd.read_csv(Stoll_imp_dir + 'Stoll_list_dir-speed_smth1E-3.csv')
Stoll_dir_speed= Stoll_dir_speed.set_index(['ID', 'time'])

S_nodes= pd.read_csv(Stoll_imp_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')
S_nodes= S_nodes.set_index(['ID', 'time'])

S_ERA5= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5.csv')
S_ERA5= S_ERA5.set_index(['ID', 'time'])

S_ERA5_2= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5_2.csv')
S_ERA5_2= S_ERA5_2.set_index(['ID', 'time'])

Stoll= pd.concat([Stoll, Stoll_dir_speed, S_nodes, S_ERA5, S_ERA5_2], axis=1)


Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps




print('plot the SOMs')
fig= plt.figure(fignr, figsize= (5.6,7) )

#axs = fig.subplots(2, 5)#, sharex=True, sharey=True)

fignr+=1
plt.clf()

import matplotlib.patches as mpatches
space= 3

for iSOM in range(1, x*y+1):
    
    print('SOM ', iSOM)
    
#    ax1= plt.subplot(y, x, iSOM)
    S_SOM= Stoll[Stoll.node== iSOM]

    
#    ax1= axs[(iSOM-1)//5, (iSOM-1)%5]
    ax1=plt.subplot(3, 3, iSOM)
    
    
    Prop_speed= S_SOM.prop_speed.mean()
    
#    arrow = mpatches.FancyArrowPatch((0, 0), (Prop_speed, 0),  mutation_scale=10)
#    ax1.add_patch(arrow)
   
    ax1.add_patch(mpatches.FancyArrowPatch((0, 0), (Prop_speed, 0),  mutation_scale=10, color='r', label='Prop') )
#    ax1.text(-1, 12, 'Propag.',  horizontalalignment='right', verticalalignment='center')
#    ax1.annotate('Propagation', xy=(Prop_speed, 16), xytext=(0, 16),  horizontalalignment='right', verticalalignment='center',  arrowprops=dict(facecolor='r')) #, shrink=0.05))

    
#    ax.annotate('', xy=(0, 0), xytext=(10, 4),
#                   arrowprops=dict(facecolor='black', shrink=0.05))
#    ax1.quiver(Prop_speed, 1)
#    
    
    wind_beering= UV2Direction(S_SOM['10u_mean-'+str(dist)], S_SOM['10v_mean-'+str(dist)])
    rot_beering = wind_beering - S_SOM.prop_dir
#    
    U= np.sqrt(S_SOM['10u_mean-'+str(dist)]**2+ S_SOM['10v_mean-'+str(dist)]**2)
#    
    rot_u, rot_v = WindSpeedDirection2UV(U, rot_beering, orientation='' )
        
    ax1.add_patch(mpatches.FancyArrowPatch((0, -3*space), (rot_v.mean(), -3*space -rot_u.mean()), color='k',  mutation_scale=10, label= '10m' ))
#    ax1.text(-1, -10, '10 m',  horizontalalignment='right', verticalalignment='center')

    


    
    
    for iplev, plevel in enumerate([925, 850, 700, 500]):#, 300]):
        if plevel < 850: iplev += 1
        wind_beering= UV2Direction(S_SOM['u'+str(plevel)+'_mean-'+str(dist)], S_SOM['v'+str(plevel)+'_mean-'+str(dist)])
        rot_beering = wind_beering - S_SOM.prop_dir
    #    rot_beering%=360
    #    
        U= np.sqrt(S_SOM['u'+str(plevel)+'_mean-'+str(dist)]**2+ S_SOM['v'+str(plevel)+'_mean-'+str(dist)]**2)
    #    
        rot_u, rot_v = WindSpeedDirection2UV(U, rot_beering, orientation='' )
        
        if plevel== 700: u_700, v_700 = rot_u, rot_v
        if plevel== 850: u_850, v_850 = rot_u, rot_v
        
        ax1.add_patch(mpatches.FancyArrowPatch((0, -2*space+ iplev*space), (rot_v.mean(), -2*space+ iplev*space -rot_u.mean()), color='k', mutation_scale=10, label= str(plevel) ))
#        ax1.text(-1, -6.5+ iplev*3.5, str(plevel)+' hPa',  horizontalalignment='right', verticalalignment='center')

        
    """yaxis left"""
    ax1.set_ylim(-3*space -3, 3*space+3)

    ax1.yaxis.set_ticks(np.arange(-3*space, 2.1*space, space))
    
    y_right = ax1.twinx()
    y_right.set_ylim(-3*space -3, 3*space+3)
    y_right.set_yticks(np.arange(-4*space, 4.1*space, space))
    y_right.set_yticklabels('')


    if iSOM in [1,4,7]:
        a=['10m', '925 hPa', '850 hPa', 'Propag.', '700 hPa', '500 hPa', '300 hPa']
        ax1.set_yticklabels(a)
        ax1.get_yticklabels()[3].set_color('red') 
    else:
        ax1.set_yticklabels('')
        
    

    """xaxis"""
    ax1.set_xlim(-.5,17)

    ax1.set_xticks(np.arange(0,17, space))
    if iSOM <7:     ax1.set_xticklabels('')
    else:
        ax1.set_xticklabels(['0', '', '6', '', '12', ''])
#    ax1.yaxis.tick_right()

#    ax1.set_yticklabels('')
    
    ax1.set_title('SOM '+str(iSOM), color= 'red', size= 13)
#        ax1.annotate(str(plevel), xy=(rot_v.mean(), -10 + iplev*5  -rot_u.mean()), xytext=(0, -10 + iplev*5 ),  horizontalalignment='right', verticalalignment='center',  arrowprops=dict(facecolor='k')) #, shrink=0.05))

    
    
    #    wind_beering= UV2Direction(S_SOM['u500_mean-'+str(dist)], S_SOM['v500_mean-'+str(dist)])
    #    rot_beering = wind_beering - S_SOM.prop_dir
    #    rot_beering%=360
    ##    
    #    U= np.sqrt(S_SOM['u500_mean-'+str(dist)]**2+ S_SOM['v500_mean-'+str(dist)]**2)
    #    
#        rot_u, rot_v = WindSpeedDirection2UV(U, rot_beering, orientation=''  )
#            
#        ax1.add_patch(mpatches.FancyArrowPatch((0, 0), (rot_v.mean(), -rot_u.mean()),  mutation_scale=10) )

#    ax1.legend()
    
#    ax1.set_xlim(-10, 16.5)
#    ax1.set_xticks([0, 5, 10, 15])
#    ax1.set_ylim(-13, 13)
#    ax1.set_yticklabels('')
#    
#
fig.text(0.5, 0.04, r"Horizontal velocity in propagation direction [m$\cdot$s$^{-1}$]", ha='center')
fig.text(0.92, 0.5, 'Horizontal velocity perpendicular to propagation [m$\cdot$s$^{-1}$]', va='center', rotation='vertical')

#fig.tight_layout()


#
if save:
    save_name=savedir+ 'SOM_Wind_vectors-'+str(dist)+'_reduced'
   
    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')
    

  





