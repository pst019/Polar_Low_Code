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
    Mediadir= '/run/media/pst019/PatsOrange/'

homedir=Mediadir+'home/'


import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import xarray as xr
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from windrose import WindroseAxes
import matplotlib as mpl
from scipy import stats

plt.rcParams.update({'font.size': 12})


import numpy as np
from f_useful import *
from f_STARS import *

#
save= True
save= False
#savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM_shear/'
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/new_class_3/'

#
fignr= 6

load_ds=True

all_figs= False #displays a lot of other figures for one node

with_cats=True #display the lines for the categorisation

#typ='Bias_cor'
typ=''
    


shear_type = 'wind_shear_2prop'
#shear_type = 'wind_shear_2mean-wind'
#shear_type = 'thermal_grad_2mean-wind'

"""the shear angle and speed"""
#l1, l2, mean= 850, 700, 500
#l1, l2, mean= 850, 700, 250
l1, l2, mean= 925, 500, 500
#l1, l2, mean= 925, 500, 250
#l1, l2, mean= 925, 700, 250


strength_level= 1.5



xSOM=3 #3
ySOM=3 #x+1


lifelim= 6 #lifetime limit





"""import Stoll list"""
Stoll_imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


Stoll= pd.read_csv(Stoll_imp_dir + 'Stoll_list_noRojo.csv')
Stoll= Stoll.set_index(['ID', 'time'])

Stoll_dir_speed= pd.read_csv(Stoll_imp_dir + 'Stoll_list_dir-speed_smth1E-3.csv')
Stoll_dir_speed= Stoll_dir_speed.set_index(['ID', 'time'])

S_nodes= pd.read_csv(Stoll_imp_dir + 'Stoll_nodes_x'+str(xSOM)+'_y'+str(ySOM)+'.csv')
S_nodes= S_nodes.set_index(['ID', 'time'])

S_ERA5= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5.csv')
S_ERA5= S_ERA5.set_index(['ID', 'time'])

S_ERA5_2= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5_2.csv')
S_ERA5_2= S_ERA5_2.set_index(['ID', 'time'])

S_ERA5_3= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5_3.csv')
S_ERA5_3= S_ERA5_3.set_index(['ID', 'time'])

Stoll= pd.concat([Stoll, Stoll_dir_speed, S_nodes, S_ERA5, S_ERA5_2, S_ERA5_3], axis=1)


Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps








if shear_type == 'wind_shear_2prop':
    speed_var= 'vert_shear_strength_vec'+str(l1)+'-'+str(l2)+'_mean-'+str(mean)
    
    #this gives the shear angle calculated from the differential wind vector in compass direction
    dir_var= 'vert_shear_angle_compass'+str(l1)+'-'+str(l2)+'_mean-'+str(mean)
    


#this gets the shear angle by the thickness gradient
#dir_var= 'vert_shear_angle_vec'+str(l1)+'-'+str(l2)+'_mean-'+str(mean)


if shear_type == 'wind_shear_2mean-wind':
    speed_var= 'vert_shear_strength_vec'+str(l1)+'-'+str(l2)+'_mean-'+str(mean)
    
    #this gives the shear angle calculated from the differential wind vector
#    dir_var= 'vert_shear_angle_vec3_'+str(l1)+'-850-700-'+str(l2)+'_mean-'+str(mean)
#    dir_var= 'vert_shear_angle_vec3_'+str(l1)+'-850-700-'+str(l2)+'_mean-'+str(mean)   
    dir_var= 'vert_shear_angle_vec3_'+str(l1)+'-'+str(l1)+'-'+str(l2)+'-'+str(l2)+'_mean-'+str(mean)



if shear_type == 'thermal_grad_2mean-wind':
    speed_var= 'vert_shear_strength_vec'+str(l1)+'-'+str(l2)+'_mean-'+str(mean)
#    speed_var= 'vert_shear_strength_vec'+str(l1)+'-'+str(l2)+'_mean-'+str(250)

    #this gets the shear angle by the thickness gradient
    dir_var= 'vert_shear_angle_vec'+str(l1)+'-'+str(l2)+'_mean-'+str(mean)


if 'compass' in dir_var:
    Stoll[dir_var] = (Stoll[dir_var] - Stoll['prop_dir'])%360   


#Stoll= Stoll[Stoll['prop_speed'] > 3]


Stoll[speed_var]*= 1E3
Stoll[dir_var] = (Stoll[dir_var]+90)%360 #to rotate to the right
#Stoll[dir_var] = (Stoll[dir_var])%360 #to rotate to the right
Stoll[dir_var] = np.deg2rad(Stoll[dir_var]) #to rotate to the right



"""new classification"""
if typ == '':
    Stoll['shear_category']= 0
    Stoll.loc[np.logical_and(Stoll[dir_var] < 3*np.pi/4, Stoll[dir_var] >= 1*np.pi/4), 'shear_category'] ='forward'
    Stoll.loc[np.logical_and(Stoll[dir_var] < 5*np.pi/4, Stoll[dir_var] >= 3*np.pi/4), 'shear_category'] ='right'
    Stoll.loc[np.logical_and(Stoll[dir_var] < 7*np.pi/4, Stoll[dir_var] >= 5*np.pi/4), 'shear_category'] ='reverse'
    Stoll.loc[np.logical_or(Stoll[dir_var] < 1*np.pi/4, Stoll[dir_var] >= 7*np.pi/4), 'shear_category'] ='left'
    Stoll.loc[Stoll[speed_var] <= strength_level, 'shear_category'] = 'weak'


S= Stoll[Stoll[speed_var] > strength_level ]



a= S[dir_var].groupby('ID').max() - S[dir_var].groupby('ID').min()
b= ((S[dir_var]+np.pi)%(2*np.pi)).groupby('ID').max() - ((S[dir_var]+np.pi)%(2*np.pi)).groupby('ID').min()

max_angle= np.min([a, b], axis= 0)


x= np.linspace(0, 2*np.pi, 100)

kernel= stats.gaussian_kde(max_angle)

fig= plt.figure(fignr, figsize= (6,6))
#fignr+= 1
#
plt.clf()
plt.plot(x, kernel(x), lw= 3 , color= 'k' )

print('Change shear angle more than 90 degree:', len(max_angle[max_angle > np.pi/2]/len(max_angle)), ' of ', len(max_angle) )



S_strong= remove_dublicate(Stoll[Stoll[speed_var] > strength_level ].reset_index()['ID'])
S_weak= remove_dublicate(Stoll[Stoll[speed_var] < strength_level ].reset_index()['ID'])

print(len(set(S_strong).intersection(set(S_weak)) ) /len(remove_dublicate(Stoll.reset_index()['ID']) ) ,
len(set(S_strong).intersection(set(S_weak)) ) ,len(remove_dublicate(Stoll.reset_index()['ID']) ) )