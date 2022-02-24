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
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM_shear/'
#
fignr= 1




#"""SOM var"""
#Obs= 'allObs_SOMprep' #the other versions were deleted
#
#x=3 #3
#y=3 #x+1
#
#mirror_top_bottom=False
#rotate_clock=False
#
#if x== 3 and y == 4: mirror_top_bottom=True
#if x== 3 and y == 3: rotate_clock=True
#
#PLCG_type, smooth_param= 'track_smth', '1E-3'
#
#lifelim= 6 #lifetime limit
#
#vo_tend_excl= False
#vo_tend_thresh= 0
#
#
#Splevel=850
#Sano=True
##ano= False
#Svar= 't'
##Svar='z'
##Svar='vo'
##Svar='U'
##Svar='q'
#
#if Sano==True: Svar_full= Svar+ '_ano'
#else: Svar_full= Svar





dist=250
 
    
    
"""import Stoll list"""
Stoll_imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


Stoll= pd.read_csv(Stoll_imp_dir + 'Stoll_list_noRojo.csv')
Stoll= Stoll.set_index(['ID', 'time'])

Stoll_dir_speed= pd.read_csv(Stoll_imp_dir + 'Stoll_list_dir-speed_smth1E-3.csv')
Stoll_dir_speed= Stoll_dir_speed.set_index(['ID', 'time'])

#S_nodes= pd.read_csv(Stoll_imp_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')
#S_nodes= S_nodes.set_index(['ID', 'time'])

S_ERA5= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5.csv')
S_ERA5= S_ERA5.set_index(['ID', 'time'])

S_ERA5_2= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5_2.csv')
S_ERA5_2= S_ERA5_2.set_index(['ID', 'time'])

Stoll= pd.concat([Stoll, Stoll_dir_speed,  S_ERA5, S_ERA5_2], axis=1) #S_nodes,


Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps








fig= plt.figure(fignr) #, figsize= (5.6,7) )
fignr+=1
plt.clf()
#axs = fig.subplots(2, 5)#, sharex=True, sharey=True)





#Stoll['u']= Stoll['u700_mean-250'] 
#Stoll['v']= Stoll['v700_mean-250'] 

#Stoll['u']= 0.5* (Stoll['u925_mean-250'] + Stoll['u700_mean-250'])
#Stoll['v']= 0.5* (Stoll['v925_mean-250'] + Stoll['v700_mean-250'])

Stoll['u']= 1/3* (Stoll['u850_mean-250'] + 2* Stoll['u700_mean-250'])
Stoll['v']= 1/3* (Stoll['v850_mean-250'] + 2* Stoll['v700_mean-250'])

Stoll['u']= 0.5* (Stoll['u925_mean-250'] + Stoll['u500_mean-250'])
Stoll['v']= 0.5* (Stoll['v925_mean-250'] + Stoll['v500_mean-250'])

#Stoll['u']= 1/5* (3*Stoll['u925_mean-250'] + 2*Stoll['u500_mean-250'])
#Stoll['v']= 1/5* (3*Stoll['v925_mean-250'] + 2*Stoll['v500_mean-250'])

Stoll['U_dir'] = UV2Direction(Stoll['u'], Stoll['v'])
Stoll['U'] = np.sqrt(Stoll['u']**2 + Stoll['v']**2)

plt.scatter(Stoll['prop_speed'], Stoll['U'])
plt.plot(np.arange(0, 22), np.arange(0, 22), c= 'r')

print('Correlation: ', np.corrcoef(Stoll['prop_speed'], Stoll['U'])[0, 1])

fig= plt.figure(fignr) #, figsize= (5.6,7) )
fignr+=1
plt.clf()


Stoll= Stoll[Stoll['prop_speed'] > 3]


plt.scatter((Stoll['prop_dir'])%360, Stoll['U_dir']%360)
plt.plot(np.arange(0, 360), np.arange(0, 360), c= 'r')

print('Correlation: ', np.corrcoef((Stoll['prop_dir'])%360, Stoll['U_dir']%360)[0, 1])




#        

if save:
    save_name=savedir+ 'Est_prop_'
   
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
