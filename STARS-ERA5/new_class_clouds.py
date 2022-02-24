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

#plt.rcParams.update({'font.size': 10})

#
save= True
#save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/new_class_3/'
#
fignr= 9

imp_list=True
#imp_list=False


cloud_hist= 'percentage' #[percentage, number] - determines if the cloud morphology histogram is displayed in number of time steps or in percentage within the SOM



#PLCG_type, proplim= 'stearing_flow', 3 #propagation speed limit

lifelim= 6 #lifetime limit








"""import Stoll list"""
imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


Stoll= pd.read_csv(imp_dir + 'Stoll_list.csv')
Stoll= Stoll.set_index(['ID', 'time'])

#S_nodes= pd.read_csv(imp_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')
#S_nodes= S_nodes.set_index(['ID', 'time'])

S_shear_cats= pd.read_csv(imp_dir + 'shear_categories.csv')
S_shear_cats= S_shear_cats.set_index(['ID', 'time'])


Stoll= pd.concat([Stoll, S_shear_cats], axis=1)


Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps






x, y= 3, 2


"""cloud distributions"""
fig= plt.figure(fignr, figsize= (2.7*x,2.7*y) )
fignr+=1
plt.clf()


Sdist= 100
#only the steps where the cloud morphology is observed and the distance is lower than Sdist
Stoll_cloud= Stoll.loc[Stoll['Rojo_lat'] > 1]
Stoll_cloud= Stoll_cloud.loc[distance((Stoll_cloud.lat, Stoll_cloud.lon), (Stoll_cloud['Rojo_lat'], Stoll_cloud['Rojo_lon'])) < Sdist]
#    Stoll_cloud= remove_morphs(Stoll_cloud, 'Morphology', count= 25, how='remove', excludelist=['-', 'U'])

#should maybe adapt this
morphlist= list(Stoll_cloud.groupby('Morph_red').size().sort_values(ascending=False).index)[:7]
if 'other' in morphlist: morphlist.sort(key = 'other'.__eq__) #this moves other to the end of the list




for iSOM in range(1, 6):
    ax1= plt.subplot(y, x, iSOM)
    
#    Morph_node= Stoll_cloud[Stoll_cloud.node== iSOM]['Morph_red']
    Morph_node= Stoll_cloud[Stoll_cloud.shear_category== iSOM]['Morph_red']


    for morph in morphlist:
#        plt.hist(Sx[Sx['Morphology'] == morph], bins= bins , label=morph, bottom= bot, edgecolor= 'k', lw= 1)
        if cloud_hist== 'number':
            ax1.bar(morph, len(Morph_node[Morph_node == morph]), edgecolor= 'k', lw= 1)
        
        if cloud_hist== 'percentage':
            ax1.bar(morph, len(Morph_node[Morph_node == morph])/len(Morph_node), edgecolor= 'k', lw= 1)
            ax1.set_ylim(0, 0.65)

#    perc_of_SOM= len(np.where(Stoll.shear_category == iSOM)[0])/len(Stoll)
    plt.title('Cat '+str(iSOM) ) #+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')
#    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')

#    ax1.axis('off')
#
#
#
plt.tight_layout()

if save:
    save_name=savedir+ 'Cats_clouds_x'+str(x)+"_y"+str(y)
    print(save_name)

    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')




"""reduced cloud distribution"""
fig= plt.figure(fignr, figsize= (6,2.7) )
fignr+=1
plt.clf()

Stoll_c2= Stoll_cloud
Stoll_c2.loc[Stoll_c2.shear_category.isin([1,2,3,4]), 'shear_category']= 1
#Stoll_c2.loc[Stoll_c2.node.isin([5,6,8]), 'shear_cate']= 5

SOM_label_list=['strong shear', 'weak shear']

for iS, iSOM in enumerate([1,5]):
    print(iSOM)
    ax1= plt.subplot(1, 2, iS+1)
    
#    Morph_node= Stoll_c2[Stoll_c2.node== iSOM]['Morph_red']
    Morph_node= Stoll_c2[Stoll_c2.shear_category== iSOM]['Morph_red']


    for morph in morphlist:
#        plt.hist(Sx[Sx['Morphology'] == morph], bins= bins , label=morph, bottom= bot, edgecolor= 'k', lw= 1)
        if cloud_hist== 'number':
            ax1.bar(morph, len(Morph_node[Morph_node == morph]), edgecolor= 'k', lw= 1)
        
        if cloud_hist== 'percentage':
            cloud_share= len(Morph_node[Morph_node == morph])/len(Morph_node)
            print(morph, cloud_share)
            ax1.bar(morph, cloud_share, edgecolor= 'k', lw= 1)
            ax1.set_ylim(0, 0.6)

#    perc_of_SOM= len(np.where(Stoll_c2.node == iSOM)[0])/len(Stoll_c2.node[Stoll_c2.node > 0])
    plt.title(SOM_label_list[iS]) #+ ' ('+str(int(np.round(perc_of_SOM*100, 0)))+'%)')
#    plt.title(SOM_label_list[iS]+ ' ('+str(int(np.round(perc_of_SOM*100, 0)))+'%)')
#    ax1.axis('off')
#
#
#
    if iS== 0: plt.ylabel('Frequency')    
    
plt.tight_layout()    



if save:
    save_name=savedir+ 'Cats_clouds_reduced'
    print(save_name)

    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')

