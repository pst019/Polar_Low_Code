#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 18:05:34 2019

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

import pandas as pd 
import numpy as np
from f_useful import *
from f_STARS import *


save= False
#save= True
#save_dir= Mediadir + '/ERA5_STARS/Stoll_list_other-landexcl/'
save_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


#Stype='system'
Stype='step'

list_type='Stoll'


system_char= 'mean'
#system_char= 'max'
#system_char='min'
#system_char='med'
#system_char='initial'
#system_char='final'




  
"""Stoll systems"""
test= 'version4'
dist= 150

Stoll_STARS_file="ERA5_STARS/STARS-PMC-merge/"+test+"_dist_"+str(dist)+"_handfix.csv"
Stoll = pd.read_csv(Mediadir+ Stoll_STARS_file)
Stoll['time']= pd.DatetimeIndex(Stoll.time)

Stoll= Stoll_Obs_nr(Stoll) #creates Stoll['Obs']
#Stoll= Stoll_excl(Stoll, landfraction= 0.5)   
Stoll= Stoll_excl(Stoll)   

#Stoll = merge_ERA5(Stoll, Mediadir, "PL/STARS/Stoll_ERA5_handfix_3.csv") #the ERA5 merge has to happen with the same stoll list as was used for the production
#Stoll= Stoll_excl_lifetime(Stoll, lifetime_threshold = 5)

Stoll= Stoll.rename(columns={'Stoll nr': 'ID', 'Rojo nr': 'Rojo_ID', 'Rojo nr old': 'Rojo_ID_orig',
                             'STARS Obs nr': 'Rojo_Obs', 'STARS lat': 'Rojo_lat', 'STARS lon': 'Rojo_lon',
                             'Press_min': 'Rojo_SLP'})
Stoll= Stoll.drop(columns=['Unnamed: 21']) # 'STARS lat', 'STARS lon', 'STARS Obs nr', 'dist'])



"""put the Obs as the first column"""
cols = list(Stoll.columns)
cols = [cols[-1]] + cols[:-1]
Stoll = Stoll[cols]


Stoll= Stoll.set_index(['ID', 'time'])
if save: Stoll.to_csv(save_dir + 'Stoll_list_prior_calcs.csv')


"""some calcs"""
Stoll= Stoll.reset_index()
Stoll, Stoll_ind= Stoll_individual_systems(Stoll, ID_name='ID', lat_name= 'Rojo_lat', lon_name= 'Rojo_lon',
                                           system_char=system_char, update_S= True)
Stoll_ind= calc_system_duration(Stoll, Stoll_ind, ID_name= 'ID', Obs_name='Obs')

Stoll= Stoll.join(Stoll_ind['Duration'], on= 'ID') #adding the duration also the Stoll
#Stoll['Lifetime fraction']= (Stoll['Obs']-1)/Stoll['Duration'] #not so good


#make the better Lifetime fraction
PLnr= remove_dublicate(Stoll.ID)
Stoll['Lifetime fraction']= -1* np.ones(len(Stoll))

for ID in PLnr:

    S_now= Stoll[Stoll.ID == ID]
    times= pd.to_datetime(S_now.index)
    dt= (times- times[0])/np.timedelta64(1,'h')
    Lifetime_frac= dt/dt[-1]    
    Stoll.loc[Stoll.ID == ID, 'Lifetime fraction']= Lifetime_frac



Stoll= Stoll_param_tend(Stoll, parameter= 'vo', outname= 'vo_tendency')


Stoll= Stoll.set_index(['ID', 'time'])


print('In Stoll list are:')

Stollnrs= remove_dublicate(Stoll.reset_index().ID)
print('Stoll centers:' , len(Stollnrs) ) 


Rojonrs= remove_dublicate(Stoll['Rojo_ID'])
print('of different Rojo centers:', len(Rojonrs) ) 

Rojoevents= remove_dublicate([ID.split('.')[0] for ID in Stoll['Rojo_ID'] ])
print('of different Rojo events:', len(Rojoevents) ) 

Rojo_primaries= [nr for nr in Rojonrs if '.' not in nr]
print('of different Rojo primaries:', len(Rojo_primaries))

if save: Stoll.to_csv(save_dir + 'Stoll_list.csv')


"""Stoll no Rojo"""
Stoll= Stoll.drop(columns=['Comment', 'track file', 'Rojo_ID', 'Rojo_ID_orig', 'Rojo_SLP', 'U_10min_knots',
                           'row_idx', 'track_idx', 'Rojo_Obs', 'Rojo_lat', 'Rojo_lon',
                           'Cloud diameter', 'Morphology', 'Morph_red', 'dist'])

if save: Stoll.to_csv(save_dir + 'Stoll_list_noRojo.csv')



"""Stoll matched with ERA_5"""
ERA5_file= "PL/STARS/Stoll_ERA5_handfix_5.csv"
S_ERA5 = pd.read_csv(Mediadir+ERA5_file, sep=',', header= 0)
S_ERA5= pd.concat([Stoll.reset_index()[['ID', 'time']],S_ERA5], axis=1)
S_ERA5= S_ERA5.set_index(['ID', 'time'])



if save: S_ERA5.to_csv(save_dir + 'Stoll_ERA5_3.csv')



"""nodes"""
x=3 #3
y=3 #x+1

Svar_full= 't_ano'
Splevel= 850
Obs= 'allObs_SOMprep'
smooth_param=  '1E-3'
lifelim= 6 

mirror_top_bottom=False
rotate_clock=False

if x== 3 and y == 4: mirror_top_bottom=True
if x== 3 and y == 3: rotate_clock=True


SOM_filedir= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/SOM/"+Svar_full+"_"+str(Splevel)+'_'+Obs+ '_track-smth-'+smooth_param+'_dur'+str(lifelim)+"_x"+ str(x)+"_y"+str(y)
df_nodes= pd.read_csv(SOM_filedir+"_node_nr.txt", sep=" ")

df_nodes['time']= pd.to_datetime(df_nodes.date.apply(str), format='%Y%m%d%H')
df_nodes['ID']= [str(int(df_nodes.PLnr[n]))[:-3]+'_'+str(int(df_nodes.PLnr[n])%1000) for n in range(len(df_nodes.PLnr))]
df_nodes= df_nodes.drop(columns=['date', 'PLnr'])
df_nodes= df_nodes.set_index(['ID', 'time'])


if mirror_top_bottom:
    df_nodes['node']= (y-1 -(df_nodes.node.values-1)//x)*x + (df_nodes.node.values -1)%x +1
if rotate_clock:
    x_n= (df_nodes['node'].values-1)%x +1
    y_n= (df_nodes['node'].values-1)//x +1
    y_n= y+1 - y_n #reverse the y_nodes
    
    df_nodes['node']= (x_n -1)*x + y_n




if save: df_nodes.to_csv(save_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')




"""evol_nodes"""
Stoll= pd.concat([Stoll, df_nodes], axis= 1)
Stoll= Stoll.reset_index()


for ind,ID in enumerate(Stoll_ind.ID):
    S= Stoll[Stoll['ID'] == ID]
    Snode_rr= [int(k) for k,_g in itertools.groupby(S.node.fillna(-1) )] #dr - repetitive remove , removes if several instances after each other are the same    

    if ind== 0:
        node_evol= [Snode_rr]
    else:
        node_evol+= [Snode_rr]

    
#    if hopp_remove: #remove a jump back and forth
    for i,node in enumerate(Snode_rr):
        if i == 0:#the first
            Snode_rr_rem_hopp= [node]
        elif i < len(Snode_rr)-1: #all middle ones
            if Snode_rr[i-1]!= Snode_rr[i+1]:
                if Snode_rr_rem_hopp[-1] != node:
                    Snode_rr_rem_hopp += [node]
            elif Snode_rr[i-1] == Snode_rr[i+1] and Snode_rr[i-1]== -1:
                #nodes between nans should not be removed (but nan between same nodes is removed by step before)
                Snode_rr_rem_hopp += [node]
        elif i == len(Snode_rr)-1: #the last
            if Snode_rr_rem_hopp[-1] != node:
                Snode_rr_rem_hopp += [node]        
#    Snode_rr= Snode_rr_rem_hopp
    
    if ind== 0:
        node_evol_hopp_rm= [Snode_rr_rem_hopp]
    else:
        node_evol_hopp_rm+= [Snode_rr_rem_hopp]
    
#    remove the -1 values
#    if m1_rem:        
    Snode_rr= [x for x in Snode_rr_rem_hopp if x != -1]    

    if ind== 0:
        node_evol_hopp_rm_nan_rm= [Snode_rr]
    else:
        node_evol_hopp_rm_nan_rm+= [Snode_rr]


Stoll_ind['node_evol']=node_evol
Stoll_ind['node_evol_hopp_rm']= node_evol_hopp_rm
Stoll_ind['node_evol_hopp_rm_nan_rm']= node_evol_hopp_rm_nan_rm


Stoll= Stoll.set_index(['ID', 'time'])

node_full= Stoll['node'].to_frame().join(Stoll_ind[['node_evol', 'node_evol_hopp_rm', 'node_evol_hopp_rm_nan_rm']], on= 'ID')


if save: node_full.to_csv(save_dir + 'Stoll_nodes_full_x'+str(x)+'_y'+str(y)+'.csv')

#to print how many
Stoll_ind['node_evol_hopp_rm_nan_rm']= Stoll_ind['node_evol_hopp_rm_nan_rm'].astype(str)
stats= Stoll_ind['node_evol_hopp_rm_nan_rm'].value_counts()#
if save:  stats.to_csv(save_dir + 'Stoll_nodes_evolution_stat_x'+str(x)+'_y'+str(y)+'.csv')
#
#
#pd.concat([Stoll, df_nodes], axis=1)

"""PCs"""
PC_dir= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/PCs/"
PC_file= "t850_allObs_track-smth1E-3_PCs.csv"
df=pd.read_csv(PC_dir+ PC_file)
df=df.drop(columns='t_PC4')
#        df= df.rename(columns={'t_PC1':'t_mat_PC1', 't_PC2':'t_mat_PC2', 't_PC3':'t_mat_PC3'})

df= df.set_index(['ID', 'Obs'])
df= pd.concat([df, Stoll.reset_index().set_index(['ID', 'Obs'])['time'] ], axis= 1)

df= df.reset_index().set_index(['ID', 'time'])

if save: df.to_csv(save_dir + 'Stoll_PCs.csv')
        




"""shear categories"""
Stoll_imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


Stoll= pd.read_csv(Stoll_imp_dir + 'Stoll_list_noRojo.csv')
Stoll= Stoll.set_index(['ID', 'time'])

Stoll_dir_speed= pd.read_csv(Stoll_imp_dir + 'Stoll_list_dir-speed_smth1E-3.csv')
Stoll_dir_speed= Stoll_dir_speed.drop(columns='Obs')
Stoll_dir_speed= Stoll_dir_speed.set_index(['ID', 'time'])

S_nodes= pd.read_csv(Stoll_imp_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')
S_nodes= S_nodes.set_index(['ID', 'time'])

S_ERA5= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5.csv')
S_ERA5= S_ERA5.set_index(['ID', 'time'])
#
S_ERA5_2= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5_2.csv')
S_ERA5_2= S_ERA5_2.set_index(['ID', 'time'])

Stoll= pd.concat([Stoll, S_nodes, Stoll_dir_speed, S_ERA5, S_ERA5_2], axis=1)


l1= 925
l2= 500
mean= 500
strength_level= 1.5


speed_var= 'vert_shear_strength_vec'+str(l1)+'-'+str(l2)+'_mean-'+str(mean)
#speed_var= 'vert_shear_strength_vec'+str(l1)+'-'+str(l2)+'_mean-'+str(250)

#this gets the shear angle by the thickness gradient
#dir_var= 'vert_shear_angle_vec'+str(l1)+'-'+str(l2)+'_mean-'+str(mean)

#this gives the shear angle calculated from the differential wind vector
#dir_var= 'vert_shear_angle_vec3_'+str(l1)+'-850-700-'+str(l2)+'_mean-'+str(mean)

#this gives the shear angle calculated from the differential wind vector in compass direction
dir_var= 'vert_shear_angle_compass'+str(l1)+'-'+str(l2)+'_mean-'+str(mean)


if 'compass' in dir_var:
    Stoll[dir_var] = (Stoll[dir_var] - Stoll['prop_dir'])%360   


#Stoll= Stoll[Stoll['prop_speed'] > 3]


Stoll[speed_var]*= 1E3
Stoll[dir_var] = (Stoll[dir_var]+90)%360 #to rotate to the right
#Stoll[dir_var] = (Stoll[dir_var])%360 #to rotate to the right
Stoll[dir_var] = np.deg2rad(Stoll[dir_var]) #to rotate to the right


Stoll['shear_category']= 0
Stoll.loc[np.logical_and(Stoll[dir_var] < 3*np.pi/4, Stoll[dir_var] >= 1*np.pi/4), 'shear_category'] = 1#'forward'
Stoll.loc[np.logical_and(Stoll[dir_var] < 5*np.pi/4, Stoll[dir_var] >= 3*np.pi/4), 'shear_category'] = 2 #'right'
Stoll.loc[np.logical_and(Stoll[dir_var] < 7*np.pi/4, Stoll[dir_var] >= 5*np.pi/4), 'shear_category'] = 3 #'reverse'
Stoll.loc[np.logical_or(Stoll[dir_var] < 1*np.pi/4, Stoll[dir_var] >= 7*np.pi/4), 'shear_category'] = 4 #'left'
Stoll.loc[Stoll[speed_var] <= strength_level, 'shear_category'] = 5 #'low'


Stoll['shear_category_corectbias']= 0
Stoll.loc[np.logical_and(Stoll[dir_var] < 5*np.pi/8, Stoll[dir_var] >= 1*np.pi/8), 'shear_category_corectbias'] = 1#'forward'
Stoll.loc[np.logical_and(Stoll[dir_var] < 9*np.pi/8, Stoll[dir_var] >= 5*np.pi/8), 'shear_category_corectbias'] = 2 #'right'
Stoll.loc[np.logical_and(Stoll[dir_var] < 13*np.pi/8, Stoll[dir_var] >= 9*np.pi/8), 'shear_category_corectbias'] = 3 #'reverse'
Stoll.loc[np.logical_or(Stoll[dir_var] < 1*np.pi/8, Stoll[dir_var] >= 13*np.pi/8), 'shear_category_corectbias'] = 4 #'left'
Stoll.loc[Stoll[speed_var] <= strength_level, 'shear_category_corectbias'] = 5 #'low'

if save: Stoll[['shear_category', 'shear_category_corectbias']].to_csv(save_dir + 'shear_categories.csv')
