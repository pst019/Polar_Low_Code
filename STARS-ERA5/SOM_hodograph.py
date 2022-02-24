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
from windrose import WindroseAxes

import numpy as np
from f_useful import *
from f_STARS import *


plt.rcParams.update({'font.size': 13})

#
save= True
#save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM_shear/'
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
#fig= plt.figure(fignr, figsize= (5.6,7) )

#axs = fig.subplots(2, 5)#, sharex=True, sharey=True)

#fignr+=1
#plt.clf()



for iSOM in range(1, x*y+1):
    
    print('SOM ', iSOM)
    fig= plt.figure(fignr, figsize= (3.,2.3)) # +0) )
    plt.clf()
    ax = fig.add_subplot(111)
    
#    plt.clf()

#    ax1= plt.subplot(y, x, iSOM)
    S_SOM= Stoll[Stoll.node== iSOM]

    
#    ax1= axs[(iSOM-1)//5, (iSOM-1)%5]
#    ax1=plt.subplot(3, 3, iSOM)
    
    
    Prop_speed= S_SOM.prop_speed.mean()
    
   
    ax.scatter(0,0, marker='x', color='k')
    
    theta = np.linspace(-np.pi/2, np.pi/2, 60)
    circlex = np.cos(theta)
    circley = np.sin(theta)
    plt.plot(5*circlex, 5*circley, c= 'k')
    plt.plot(10*circlex, 10*circley, c= 'k')
    plt.plot(15*circlex, 15*circley, c= 'k')
    
#    ax1.add_patch(mpatches.FancyArrowPatch((0, 0), (Prop_speed, 0),  mutation_scale=10, color='r', label='Prop') )
    ax.scatter(Prop_speed, 0, marker='x', c='r', s= 100, label= 'Propag.') #(-1, 12, 'Propag.',  horizontalalignment='right', verticalalignment='center')

##    
# 
    hod_u, hod_v=[], []
    
    """10m"""
#    wind_beering= UV2Direction(S_SOM['10u_mean-'+str(dist)], S_SOM['10v_mean-'+str(dist)])
#    rot_beering = wind_beering - S_SOM.prop_dir
#    
#    U= np.sqrt(S_SOM['10u_mean-'+str(dist)]**2+ S_SOM['10v_mean-'+str(dist)]**2)
#    
#    rot_u, rot_v = WindSpeedDirection2UV(U, rot_beering, orientation='' )
#        
#    hod_u += [rot_u.mean()]
#    hod_v += [rot_v.mean()]
#    
#    ax1.scatter(rot_v.mean(), -rot_u.mean(), marker= '^', color='k' )


    """plevels"""
    for iplev, plevel in enumerate([925, 850, 700, 500]):#, 300]):
        wind_beering= UV2Direction(S_SOM['u'+str(plevel)+'_mean-'+str(dist)], S_SOM['v'+str(plevel)+'_mean-'+str(dist)])
        rot_beering = wind_beering - S_SOM.prop_dir
 
        U= np.sqrt(S_SOM['u'+str(plevel)+'_mean-'+str(dist)]**2+ S_SOM['v'+str(plevel)+'_mean-'+str(dist)]**2)
        rot_u, rot_v = WindSpeedDirection2UV(U, rot_beering, orientation='' )
        hod_u += [rot_u.mean()]
        hod_v += [rot_v.mean()]

        if iplev== 0:
            ax.scatter(rot_v.mean(), -rot_u.mean(), marker= 's', color='k', s= 100 , label= '925 hPa')        
        elif iplev== 3:
            ax.scatter(rot_v.mean(), -rot_u.mean(), marker= '^', color='k', s= 100 , label= '500 hPa')
        elif iplev== 1:
            ax.scatter(rot_v.mean(), -rot_u.mean(), marker ='d', color='k' , label= '850 hPa')
        else:
            ax.scatter(rot_v.mean(), -rot_u.mean(), color='k', label= '700 hPa' )            

    plt.plot(hod_v, -1*np.array(hod_u), color='k')

    plt.xlim(-.5,18)
#    plt.xlabel('Wind in Prop. dir')
    plt.xticks([0, 5, 10, 15])
    plt.ylim(-7, 7)
#

    if save:
        save_name=savedir+ 'Hodo_SOM_'+str(iSOM)
       
        print(save_name)
        plt.savefig(save_name , bbox_inches='tight')


    """make the legend"""
    if iSOM == 1:
        figlegend = plt.figure(fignr+1, figsize=(8,.55))
        plt.clf()
        
        plt.figlegend(*ax.get_legend_handles_labels(), ncol= 5 )

        if save:
            save_name=savedir+ 'Hodo_legend'
           
            print(save_name)
            plt.savefig(save_name , bbox_inches='tight')
#    
#
#  
#
#
#
#
#
