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
    Mediadir= '/media/'+user+'/1692A00D929FEF8B/'

else:
    Mediadir= '/run/media/pst019/1692A00D929FEF8B/'
'/home/'+user+'/home/'

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


save= True
#save= False
#savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM/'

save_extra=''

##
fignr= 3

#imp_list=True
imp_list=False

imp_rot_nc=True
#imp_rot_nc=False


hopp_remove=True #remove a jump back and forth from the nodes, e.g. if the order of the nodes for one PL ('8_1') is [3, 2, -1, 4, 7, 10, 7, -1, 7] it is reduced to [3, 2, -1, 4, 7]
if hopp_remove: save_extra+= '_hopp_rem'

m1_rem=True #remove -1 nodes (nan nodes), due to slow motion or leaving the domain
if m1_rem: save_extra+= '_m1_rem'

"""SOM var"""
Obs= 'allObs_SOMprep'

plevel=850
ano=True
#ano= False
var= 't'
#var='z'
#var='vo'
#var='U'
#var='q'

if ano==True: var_full= var+ '_ano'
else: var_full= var

"""evolve var"""
Obs_evol='allObs'

if Obs_evol== 'mature': Obs_evol_str= Obs_evol +str(mat_evol_step)
elif Obs_evol =='allObs': Obs_evol_str= Obs_evol
else: Obs_evol_str= 'Obsnr'+str(Obs_evol)

ano_evol=True
#ano_evol=False
var_evol='t'

if ano_evol==True: var_evol_full= var_evol+ '_ano'
else: var_evol_full= var_evol

ano_evol_cont=True
var_evol_cont='z'

if ano_evol_cont==True: var_evol_cont_full= var_evol_cont+ '_ano'
else: var_evol_cont_full= var_evol_cont


if var_evol_full == 't_ano': 
    levels = np.linspace(-7, 7, 15)
if var_evol_cont_full == 'z_ano':
    cont_levels = np.linspace(-80, 80, 17)


x=3 #3
y=x+1




sym= False
#cmap='jet'
#cmap= 'viridis'
#cmap= 'Greys'
#cmap= 'Greys_r'
#cmap= 'Blues'
#cmap= 'Reds'

if "ano" in var_full: sym= True
cmap= 'RdBu_r'
#
#
filedir= Mediadir + "ERA5_STARS/PL_centred_fields/SOM/"+var_full+"_"+str(plevel)+'_'+Obs+"_x"+str(x)+"_y"+str(y)
ds_init= xr.open_dataset(filedir+"_cluster.nc")

#txt = pd.read_csv(filedir+"_node_nr.txt", sep=" ")



"""import Stoll list"""
if imp_list:
    """Stoll systems"""
    test= 'version4'
    dist= 150
    system_char='initial' # for the Stoll_ind
    
    Stoll_STARS_file="ERA5_STARS/STARS-PMC-merge/"+test+"_dist_"+str(dist)+"_handfix.csv"
    Stoll = pd.read_csv(Mediadir+ Stoll_STARS_file)
    Stoll['time']= pd.DatetimeIndex(Stoll.time)
    
    Stoll= Stoll_Obs_nr(Stoll) #creates Stoll['Obs']
    Stoll= Stoll_excl(Stoll) #, lifetime_threshold=6) - does not work since it should be excluded before SOMs are done
    Stoll= Stoll.rename(columns={'Stoll nr': 'ID'})
    
    Stoll = merge_ERA5(Stoll, Mediadir, "PL/STARS/Stoll_ERA5_handfix.csv")
    Stoll= Stoll.rename(columns={
        'U_max_200': 'Wind Speed 10m (max)',
        'vo': 'Vorticity$_{850}$ (centre)',
        'slp': 'Sea-level pressure (min)',
        'blh_max_200': 'Boundary layer height (max)',
        'cape_max_200': 'CAPE (max)', 
        'skt_med_200': 'Skin temperature (med)',
       'skt-t500_max_200': 'SST- T$_{500}$ (max)', 
       'skt-t700_max_200': 'SST -T$_{700}$ (max)',
        'tp_mean_200': 'Total precip. (mean) ',
        'cp_mean_200': 'Convective precip. (mean)',
        'sf_mean_200': 'Snow fall (mean)',
        'lsp_mean_200': 'Large-scale precip. (mean)', 
        'sshf_mean_200': 'Sensible heat flux (mean)',
        'slhf_mean_200': 'Latent heat flux (mean)',
       'grad_t850_max_200': 'Grad T$_{850}$ (max)',
       'baroclinic_gr_filter4_925-700_max_200': 'Baroclinic growth (max)' ,
       'barotropic_gr_filter4_850_max_200': 'Barotropic growth (max)',
       'vert_shear_angle925-700_mean_200': 'Vertical shear angle (mean)'
       })
    
    
    Stoll= Stoll.drop(columns=['Comment', 'track file', 'Rojo nr', 'Rojo nr old', 'Press_min', 'U_10min_knots', 'row_idx', 'track_idx'])
       
    
    """Stoll_individual_systems"""
    Stoll_ind= Stoll_individual_systems(Stoll, ID_name='ID', system_char=system_char)
#    Stoll_ind= calc_system_duration(Stoll, Stoll_ind, ID_name= 'ID', Obs_name='Obs')
    
#    Stoll_ind= Stoll_ind.drop(columns=[ 'STARS lat', 'STARS lon', 'STARS Obs nr'])
#    Stoll= Stoll.drop(columns=[ 'STARS lat', 'STARS lon', 'STARS Obs nr', 'dist'])
    Stoll_ind= calc_Obs_mature(Stoll, Stoll_ind, intensity_var='Vorticity$_{850}$ (centre)')

    df_nodes= pd.read_csv(filedir+"_node_nr.txt", sep=" ")

    
#    if Obs_evol=='allObs':
    df_nodes['time']= pd.to_datetime(df_nodes.date.apply(str), format='%Y%m%d%H')
    df_nodes['ID']= [str(int(df_nodes.PLnr[n]))[:-3]+'_'+str(int(df_nodes.PLnr[n])%1000) for n in range(len(df_nodes.PLnr))]

#        df_nodes= df_nodes.drop(columns=['date', 'PLnr'])
    Stoll= Stoll.set_index(['ID', 'time'])
    df_nodes= df_nodes.set_index(['ID', 'time'])
    
    Stoll= pd.concat([Stoll, df_nodes], axis= 1)




"""make the figure on which the original SOMs are displayed"""
fig= plt.figure(fignr, figsize= (2*x,2*y) )
fignr+=1
plt.clf()

if sym: vextr= np.max([np.max(ds_init.field), -np.min(ds_init.field)])
else: vmax, vmin= float(np.max(ds_init.field)), float(np.min(ds_init.field))
    
for iSOM in range(1, x*y +1):
    ax1= plt.subplot(y, x, iSOM)
    if sym:
        cf= plt.contourf(ds_init.x, ds_init.y, ds_init.field.sel(SOM=iSOM), cmap= cmap, vmin= -vextr, vmax= vextr)
        cs= plt.contour(ds_init.x, ds_init.y, ds_init.field.sel(SOM=iSOM), colors='k', vmin= -vextr, vmax= vextr, linewidth= 1)
    else:
        cf= plt.contourf(ds_init.x, ds_init.y, ds_init.field.sel(SOM=iSOM), cmap= cmap, vmin= vmin, vmax= vmax)
        cs= plt.contour(ds_init.x, ds_init.y, ds_init.field.sel(SOM=iSOM), colors='k', vmin= vmin, vmax= vmax, linewidth= 1)

    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.1f')

    perc_of_SOM= len(np.where(df_nodes.node == iSOM)[0])/len(df_nodes.node)
    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')



plt.tight_layout()

fig.subplots_adjust(bottom=0.1)
cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.015])
cb= fig.colorbar(cf, cax=cbar_ax, orientation="horizontal")

if var_full== 'z_ano': labelvar= 'Geopotential height anomaly [m]'
elif var_full== 'z': labelvar= 'Geopotential height [m]'

elif var_full== 't_ano': labelvar= 'Temperature anomaly [K]'
elif var_full== 't': labelvar= 'Temperature [K]'

elif var_full== 'U_ano': labelvar= 'Wind speed anomaly [m/s]'
elif var_full== 'U': labelvar= 'Wind speed [m/s]'

elif var_full== 'q_ano': labelvar= 'Specific humidity anomaly [g(kg)]'
elif var_full== 'q': labelvar= 'Specific humidity [g/kg]'
else: labelvar= var_full

cb.set_label(labelvar, size=14) 

if save:
    save_name=savedir+ 'SOM_'+var_full+"_"+str(plevel)+"_"+Obs+"_x"+str(x)+"_y"+str(y)
    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')



"""evolution"""
fig= plt.figure(fignr, figsize= (2*x,2*y + 1))
plt.clf()
fignr+=1


axs = fig.subplots(y, x)
axs= np.ravel(axs)

fig.subplots_adjust(wspace=0.25, hspace=0.20, left=0.1 -0.01*x, right=0.99, top=0.93+ 0.01*x, bottom=0.1)# 1.2/(1.5+ 2*y))



print('import from all time step list')
ds= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields/" +var_evol + '_'+ str(plevel) + '_allObs.nc')
ds2= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields/" +var_evol_cont + '_'+ str(plevel) + '_allObs.nc')
ds= xr.merge([ds, ds2])


"""add the node to ds - only for allObs"""
#if Obs_evol == 'allObs':
ds_node_var= np.zeros(len(ds.time))

for s in range(len(df_nodes)):
    df_nodes_s= df_nodes.iloc[s]
    ds_index= np.argwhere(np.logical_and(ds.PLnr.values== df_nodes_s.PLnr, ds.time.values== np.datetime64(df_nodes_s.name[1])))[0][0]
    ds_node_var[ds_index]= df_nodes_s.node

ds['node']= (('time'), ds_node_var)
ds_Obs_evol= ds
    
#for s in range(len(Stoll)):
#    Stoll_s= Stoll.iloc[s]
#    ds_index= np.argwhere(np.logical_and(ds.PLnr.values== Stoll_s.PLnr, ds.Obs.values== Stoll_s.Obs))[0][0]
#    ds_node[ds_index]= Stoll_s.node

if ano_evol == True:
    ds_Obs_evol[var_evol_full]= (('time', 'x', 'y'), ds_Obs_evol[var_evol].values- np.mean(ds_Obs_evol[var_evol].values, axis= (1,2))[:, np.newaxis, np.newaxis])
if ano_evol_cont == True:
    ds_Obs_evol[var_evol_cont_full]= (('time', 'x', 'y'), ds_Obs_evol[var_evol_cont].values- np.mean(ds_Obs_evol[var_evol_cont].values, axis= (1,2))[:, np.newaxis, np.newaxis])
    



"""plot the SOMs"""
for iSOM in range(1, x*y +1):
    
    
    if Obs_evol != 'allObs':
        node_list= Stoll_ind[Stoll_ind.node== iSOM].index.values
#        weight_list= 1/Stoll_ind[Stoll_ind.node== iSOM]['K_SOM.distances'].values
        ds_node= ds_Obs_evol.where(ds_Obs_evol.ID.isin(node_list), drop=True)
    
    if Obs_evol == 'allObs':
        ds_node= ds_Obs_evol.where(ds_Obs_evol.node == iSOM, drop=True)
    
    ax= axs[iSOM-1]
    if sym:
#        cf= plt.contourf(ds_node.x, ds_node.y, np.average(ds_node[var_evol_full], axis= 0, weights=weight_list), cmap= cmap, vmin= -vextr, vmax= vextr)
#        cf= plt.pcolor(ds_node.x, ds_node.y, np.average(ds_node[var_evol_full], axis= 0, weights=weight_list), cmap= cmap, vmin= -vextr, vmax= vextr)

        cf= ax.contourf(ds_node.x, ds_node.y, np.average(ds_node[var_evol_full], axis= 0), levels= levels, cmap= cmap, extend='both')
        
        cs= ax.contour(ds_node.x, ds_node.y, np.mean(ds_node[var_evol_cont_full], axis= 0), levels= cont_levels, colors='k')
    else:
        cf= plt.contourf(ds_node.x, ds_node.y, np.mean(ds_node[var_evol_full], axis= 0), cmap= cmap, vmin= vmin, vmax= vmax)
        cs= plt.contour(ds_node.x, ds_node.y, np.mean(ds_node[var_evol_full], axis= 0), colors='k', vmin= vmin, vmax= vmax)

    ax.clabel(cs, cont_levels[::2], fontsize=10, fmt='%1.0f')#, inline=1)

    perc_of_SOM= len(np.where(df_nodes.node == iSOM)[0])/len(df_nodes.node)
    ax.set_title(str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')

    ax.set_xticklabels('')
    ax.set_yticklabels('')
    ax.set_yticks(np.arange(-500,501, 250))
    ax.set_xticks(np.arange(-500,501, 250))
    



"""make the arrows, and circles"""

n_nodes= x*y
index= ['initial', -1]+ list(np.arange(1, n_nodes+1)) #-1 is for the nan values
columns= ['lysis', -1]+ list(np.arange(1, n_nodes+1))

node_evolv_df= pd.DataFrame(np.zeros((n_nodes+2,n_nodes+2)), index= index, columns=columns)
node_mature_df= pd.DataFrame(np.zeros((n_nodes+1, 1)), index= index[1:], columns=['mature_node'])

for ind in Stoll_ind.index.values:
    S= Stoll.loc[ind]
    Snode_rr= [int(k) for k,_g in itertools.groupby(S.node.fillna(-1) )] #dr - repetitive remove , removes if several instances after each other are the same    

    
    if hopp_remove: #remove a jump back and forth
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
        Snode_rr= Snode_rr_rem_hopp
    
    #remove the -1 values
    if m1_rem:        
        Snode_rr= [x for x in Snode_rr if x != -1]        
    
    if len(Snode_rr) > 0: #it can happen that the only node is a nan value, which is excluded by the line before
        
        node_evolv_df.loc['initial'][Snode_rr[0]] +=1
        node_evolv_df.loc[Snode_rr[-1]]['lysis'] +=1
        
        for ni in range(len(Snode_rr) -1):
            node_evolv_df.loc[Snode_rr[ni]][Snode_rr[ni+1]] +=1
    
        mature_node= Stoll[Stoll.Obs == Stoll_ind.loc[ind].Obs_mature].loc[ind].node.values[0] #find the node of the mature stage for this PL
        if np.isnan(mature_node):
            if m1_rem: 
                S_dropna= S.dropna(subset=['node'])
                mature_node= S_dropna.iloc[S_dropna['Vorticity$_{850}$ (centre)'].values.argmax()].node
            else: mature_node = -1
        node_mature_df.loc[mature_node] += 1

##to get it in percentage
#for ind in node_evolv_df.index:
#    node_evolv_df.loc[ind]/= node_evolv_df.loc[ind].sum()


bbox_right = dict(boxstyle="rarrow", fc=(0.8, 0.9, 0.9), ec="b", lw=2) #fc - background color, ec - edge color
bbox_left = dict(boxstyle="larrow", fc=(0.8, 0.9, 0.9), ec="b", lw=2)
bbox_gen = dict(boxstyle="circle", fc="lightgreen", ec="g", lw=2)
bbox_mat = dict(boxstyle="circle", fc="lightsalmon", ec="r", lw=2)
bbox_lysis = dict(boxstyle="circle", fc="wheat", ec="orange", lw=2)


arrow_limit=5
gen_limit= 3

def text_size(value, minsize=5, maxsize= 13, minval=arrow_limit, maxval=node_evolv_df[2:].max()[2:].max()): #the [2:] is done to exclude the nan, genesis and lysis
    """calculates a linear interpolation of the text size"""
    if value >= minval:
        size= minsize + (maxsize-minsize)* (value-minval)/(maxval -minval)
    else: size= 0
    return size

for node in range(1, x*y +1):
    goes_to_node= node_evolv_df.loc[node]
    ax= axs[node-1]
    
    #if goes_to_node -1 =left
    if (node-1)%x != 0 and goes_to_node[node-1] > arrow_limit:
        value= goes_to_node[node-1]
        ax.text(-500, -150, str(int(value)) , ha="right", va="center", size=text_size(value), fontweight='bold', bbox=bbox_left)
    
    #if goes_to_node +1 = right
    if (node-1)%x != (x-1) and goes_to_node[node+1] > arrow_limit:
        value= goes_to_node[node+1]
        ax.text(500, -325, str(int(value)) , ha="left", va="center", size=text_size(value), fontweight='bold', bbox=bbox_right)
    
    #if goes_to_node row+1 (node+x) = down
    if node <= x*(y-1) and goes_to_node[node+x] > arrow_limit:
        value= goes_to_node[node+x]
        ax.text(-350, -500, str(int(value)) , ha="center", va="top", rotation= 90, size=text_size(value), fontweight='bold', bbox=bbox_left, zorder=2)
    #if goes_to_node row-1 (node-x) = up
    if node > x and goes_to_node[node-x] > arrow_limit:
        value= goes_to_node[node-x]
        ax.text(350, 500, str(int(value)) , ha="center", va="bottom", rotation= 90, size=text_size(value), fontweight='bold', bbox=bbox_right)
    
    #if goes_to_node (-1 row-1) (node-x-1) =diagonal up left
    if node > x and (node-1)%x != 0 and goes_to_node[node-x-1] > arrow_limit:
        value= goes_to_node[node-x-1]
        ax.text(-500, 500, str(int(value)) , ha="center", va="center", rotation= 315, size=text_size(value), fontweight='bold', bbox=bbox_left)
        
#        ax.text(-420, 450, str(int(value)) , ha="right", va="bottom", rotation= 315, size=text_size(value), fontweight='bold', bbox=bbox_left)
    #if goes_to_node (+1 row-1) (node-x+1) =diagonal up right
    if node > x and (node-1)%x != (x-1) and goes_to_node[node-x+1] > arrow_limit:
        value=goes_to_node[node-x+1]
        ax.text(500, 500, str(int(value)) , ha="center", va="center", rotation= 45, size=text_size(value), fontweight='bold', bbox=bbox_right)
    #if goes_to_node (-1 row+1) (node+x-1) =diagonal down left
    if node <= x*(y-1) and (node-1)%x != 0 and goes_to_node[node+x-1] > arrow_limit:
        value= goes_to_node[node+x-1]
        ax.text(-500, -500, str(int(value)) , ha="center", va="center", rotation= 45, size=text_size(value), fontweight='bold', bbox=bbox_left)
    #if goes_to_node (-1 row+1) (node+x+1) =diagonal down right
    if node <= x*(y-1) and (node-1)%x != (x-1) and goes_to_node[node+x+1] > arrow_limit:
        value= goes_to_node[node+x+1]
        ax.text(500, -500, str(int(value)) , ha="center", va="center", rotation= 315, size=text_size(value), fontweight='bold', bbox=bbox_right)

    #genesis
    value= node_evolv_df.loc['initial'][node]
    if value > gen_limit:
        ax.text(-600, 400, str(int(value)), ha="center", va="center", size=text_size(value, minsize=5, maxsize= 13, minval= 3, maxval= node_mature_df.max().max()), fontweight='bold', bbox=bbox_gen)

    #mature
    value= node_mature_df.loc[node].values[0]
    if value > gen_limit:
        ax.text(-600, 225, str(int(value)), ha="center", va="center", size=text_size(value, minsize=5, maxsize= 13, minval= 3, maxval= node_mature_df.max().max()), fontweight='bold', bbox=bbox_mat)

    #lysis
    value= goes_to_node['lysis']
    if value > gen_limit:
        ax.text(-600, 50, str(int(value)), ha="center", va="center", size=text_size(value, minsize=5, maxsize= 13, minval= 3, maxval= node_mature_df.max().max()), fontweight='bold', bbox=bbox_lysis)





#plt.tight_layout()

cbar_ax = fig.add_axes([0.05, 0.06, 0.4, 0.015])
cb= fig.colorbar(cf, cax=cbar_ax, orientation="horizontal")

"""make the circle legend"""
circle_ax = fig.add_axes([0.52, 0.03, 0.25, 0.04])
import matplotlib.patches as mpatches

#arrow = mpatches.FancyArrowPatch((0.2, 0.7), (.7, 0.7), mutation_scale=20)
#arrow_ax.add_patch(arrow)
circle_ax.text(0, 0.5, 'evo-\nlution', ha="center", va="center", size=2*y, fontweight='bold', bbox=bbox_right)
circle_ax.text(0.4, 0.5, 'gen-\nesis', ha="center", va="center", size=2*y, fontweight='bold', bbox=bbox_gen)
circle_ax.text(0.7, 0.5, 'mat-\nure', ha="center", va="center", size=2*y, fontweight='bold', bbox=bbox_mat)
circle_ax.text(1, 0.5, 'lysis', ha="center", va="center", size=2*y, fontweight='bold', bbox=bbox_lysis)

#arrow_ax.text(0,0 , 'Propagation direction', fontweight= 'bold')#, fontsize=13)
circle_ax.set_frame_on(False)
circle_ax.get_xaxis().set_visible(False)
circle_ax.get_yaxis().set_visible(False)


"""make the arrow"""
arrow_ax = fig.add_axes([0.82, 0.03, 0.12, 0.04])
import matplotlib.patches as mpatches

arrow = mpatches.FancyArrowPatch((0.2, .8), (.7, .8), mutation_scale=20)
arrow_ax.add_patch(arrow)
arrow_ax.text(0,-0.2 , 'Propagation\ndirection', fontweight= 'bold', fontsize=8)
arrow_ax.set_frame_on(False)
arrow_ax.get_xaxis().set_visible(False)
arrow_ax.get_yaxis().set_visible(False)



if var_evol_full== 'z_ano': labelvar= 'Geopotential height anomaly [m]'
elif var_evol_full== 'z': labelvar= 'Geopotential height [m]'

elif var_evol_full== 't_ano': labelvar= 'Temperature anomaly [K]'
elif var_evol_full== 't': labelvar= 'Temperature [K]'

elif var_evol_full== 'U_ano': labelvar= 'Wind speed anomaly [m/s]'
elif var_evol_full== 'U': labelvar= 'Wind speed [m/s]'

elif var_evol_full== 'q_ano': labelvar= 'Specific humidity anomaly [g(kg)]'
elif var_evol_full== 'q': labelvar= 'Specific humidity [g/kg]'
else: labelvar= var_evol_full

cb.set_label(labelvar, size=13) 



if save:
    save_name=savedir+ 'SOM_evolution_'+var_full+"_"+str(plevel)+"_"+Obs+'--mean_'+var_evol_full+"_"+var_evol_cont_full+"_"+str(plevel)+"_"+Obs_evol_str+ save_extra+"_x"+str(x)+"_y"+str(y)
    print(save_name)
    plt.savefig(save_name , dpi= 140)#, bbox_inches='tight')


print('could include steps towards nan - so steps to no propagation speed')