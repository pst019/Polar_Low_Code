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
import numpy as np
import itertools
from f_useful import *
from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs


#save= True
save= False
savedir= homedir+ 'Polar_Low/STARS-ana/Figs/SOM2/'

save_extra=''

##
fignr= 3

imp_list=True
#imp_list=False

imp_rot_nc=True
#imp_rot_nc=False


hopp_remove=True #remove a jump back and forth from the nodes, e.g. if the order of the nodes for one PL ('8_1') is [3, 2, -1, 4, 7, 10, 7, -1, 7] it is reduced to [3, 2, -1, 4, 7]
if hopp_remove: save_extra+= '_hopp_rem'

m1_rem=True #remove -1 nodes (nan nodes), due to slow motion or leaving the domain
if m1_rem: save_extra+= '_m1_rem'

"""SOM var"""
Obs= 'allObs_SOMprep'

x=3 #3
y=3 #x+1

mirror_top_bottom=False
rotate_clock=False

if x== 3 and y == 4: mirror_top_bottom=True
if x== 3 and y == 3: rotate_clock=True


#PLCG_type, proplim= 'stearing_flow', 3 #propagation speed limit
PLCG_type, smooth_param= 'track_smth', '1E-3'

lifelim= 6 #lifetime limit

Splevel=850
Sano=True
#Sano= False
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
Obs_evol='allObs'

Obs_evol_str= Obs_evol


ano_var=True
#ano_var=False
var='t'
plevel_var=850
cmap= 'RdBu_r'


var_full= var
if ano_var==True:  var_full= var_full+ '_ano'
if type(plevel_var) == int: var_full= var_full+ '_'+str(plevel_var)


if 't_ano' in var_full:  levels = np.linspace(-7, 7, 15)



"""contour variable"""
ano_cont=True
cont='z'
plevel_cont=850


cont_full= cont
if ano_cont==True:  cont_full += '_ano'
if type(plevel_cont) == int: cont_full += '_'+str(plevel_cont)

if 'z_ano' in cont_full: cont_levels = np.arange(-200, 201, 10)




"""import Stoll list"""
if imp_list:
#    SOM_filedir= Mediadir + "ERA5_STARS/PL_centred_fields/SOM/"+Svar_full+"_"+str(Splevel)+'_'+Obs+"_x"+str(x)+"_y"+str(y)
    if PLCG_type== 'stearing_flow':
        SOM_filedir= Mediadir + "ERA5_STARS/PL_centred_fields/SOM/"+Svar_full+"_"+str(Splevel)+'_'+Obs+'_vel'+str(proplim)+'_dur'+str(lifelim)+"_x"+ str(x)+"_y"+str(y)
    elif PLCG_type== 'track_smth':
        SOM_filedir= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/SOM/"+Svar_full+"_"+str(Splevel)+'_'+Obs+ '_track-smth-'+smooth_param+'_dur'+str(lifelim)+"_x"+ str(x)+"_y"+str(y)

    ds_init= xr.open_dataset(SOM_filedir+"_cluster.nc")

    Stoll, Stoll_ind= imp_standard_Stoll() #SOM_filedir)#, Obs_evol)
    df_nodes= pd.read_csv(SOM_filedir+"_node_nr.txt", sep=" ")

#    df_nodes_orig= df_nodes
#    if mirror_top_bottom:
#        df_nodes['node']= (y-1 -(df_nodes.node.values-1)//x)*x + (df_nodes.node.values -1)%x +1
 
    if Obs_evol=='allObs':
        df_nodes['time']= pd.to_datetime(df_nodes.date.apply(str), format='%Y%m%d%H')
        df_nodes['ID']= [str(int(df_nodes.PLnr[n]))[:-3]+'_'+str(int(df_nodes.PLnr[n])%1000) for n in range(len(df_nodes.PLnr))]
#        df_nodes= df_nodes.drop(columns=['date', 'PLnr'])
        Stoll= Stoll.set_index(['ID', 'time'])
        df_nodes= df_nodes.set_index(['ID', 'time'])   
        
        """to rotate/mirror the nodes"""
        df_nodes_orig= df_nodes.copy()
        if mirror_top_bottom:
            df_nodes['node']= (y-1 -(df_nodes.node.values-1)//x)*x + (df_nodes.node.values -1)%x +1
        if rotate_clock:
            x_n= (df_nodes['node'].values-1)%x +1
            y_n= (df_nodes['node'].values-1)//x +1
            y_n= y+1 - y_n #reverse the y_nodes
            
            df_nodes['node']= (x_n -1)*x + y_n
            
            
        Stoll= pd.concat([Stoll, df_nodes], axis= 1)

    




"""make the figure on which the original SOMs are displayed"""
fig= plt.figure(fignr, figsize= (2.5*x,2.5*y +0) )
fignr+=1
plt.clf()

if S_sym: vextr= np.max([np.max(ds_init.field), -np.min(ds_init.field)])
else: vmax, vmin= float(np.max(ds_init.field)), float(np.min(ds_init.field))
    
for iSOM in range(1, x*y +1):
    ax1= plt.subplot(y, x, iSOM)
    if S_sym:
        cf= plt.contourf(ds_init.x, ds_init.y, ds_init.field.sel(SOM=iSOM), cmap= S_cmap, vmin= -vextr, vmax= vextr)
        cs= plt.contour(ds_init.x, ds_init.y, ds_init.field.sel(SOM=iSOM), colors='k', vmin= -vextr, vmax= vextr, linewidth= 1)
    else:
        cf= plt.contourf(ds_init.x, ds_init.y, ds_init.field.sel(SOM=iSOM), cmap= S_cmap, vmin= vmin, vmax= vmax)
        cs= plt.contour(ds_init.x, ds_init.y, ds_init.field.sel(SOM=iSOM), colors='k', vmin= vmin, vmax= vmax, linewidth= 1)

    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.1f')

    perc_of_SOM= len(np.where(df_nodes_orig.node == iSOM)[0])/len(df_nodes_orig.node)
    plt.title('SOM '+str(iSOM)+ ' ('+str(np.round(perc_of_SOM*100, 1))+'%)')


plt.tight_layout()

fig.subplots_adjust(bottom=0.1)
cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.015])
cb= fig.colorbar(cf, cax=cbar_ax, orientation="horizontal")

if Svar_full== 'z_ano': Slabelvar= 'Geopotential height anomaly [m]'
elif Svar_full== 'z': Slabelvar= 'Geopotential height [m]'

elif Svar_full== 't_ano': Slabelvar= 'Temperature anomaly [K]'
elif Svar_full== 't': Slabelvar= 'Temperature [K]'

elif Svar_full== 'U_ano': Slabelvar= 'Wind speed anomaly [m/s]'
elif Svar_full== 'U': Slabelvar= 'Wind speed [m/s]'

elif Svar_full== 'q_ano': Slabelvar= 'Specific humidity anomaly [g(kg)]'
elif Svar_full== 'q': Slabelvar= 'Specific humidity [g/kg]'
else: Slabelvar= Svar_full

cb.set_label(Slabelvar, size=14) 

if save:
    if PLCG_type== 'stearing_flow': save_name=savedir+ 'SOM_'+Svar_full+"_"+str(Splevel)+"_"+Obs+"_x"+str(x)+"_y"+str(y)
    elif PLCG_type== 'track_smth': save_name=savedir+ 'SOM_track-smth'+smooth_param+'_'+Svar_full+"_"+str(Splevel)+"_"+Obs+"_x"+str(x)+"_y"+str(y)

    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')



"""evolution"""
fig= plt.figure(fignr, figsize= (2*x+.5,2*y + 1))
plt.clf()
fignr+=1


axs = fig.subplots(y, x)
axs= np.ravel(axs)

fig.subplots_adjust(wspace=0.25, hspace=0.20, left=0.1 -0.01*x, right=0.99, top=0.93+ 0.01*x, bottom=0.1)# 1.2/(1.5+ 2*y))


if imp_rot_nc:

    print('import from all time step list')
    
    """file for shading variable"""
    if plevel_var != '': plev_str= '_'+ str(plevel_var)
    else: plev_str= plevel_var
    
    if PLCG_type== 'stearing_flow': file= Mediadir + "ERA5_STARS/PL_centred_fields/" +var + plev_str + '_allObs.nc'
    elif PLCG_type== 'track_smth': file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +var + plev_str + '_allObs_track-smth-'+smooth_param+'.nc'


    
    if os.path.isfile(file): ds= xr.open_dataset(file)
    else:
        if PLCG_type== 'stearing_flow': ds= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields/" +var + '_all_levs_allObs.nc')
        elif PLCG_type== 'track_smth': ds= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +var + '_all_levs_allObs_track-smth-'+smooth_param+'.nc')
            
        ds= ds.sel(plev= plevel_var)

    ds= ds.rename({var: var_full})

    """file for cont variable"""
    if plevel_cont != '': plev_str= '_'+ str(plevel_cont)
    else: plev_str= plevel_cont
    
    if PLCG_type== 'stearing_flow': file= Mediadir + "ERA5_STARS/PL_centred_fields/" +cont + plev_str + '_allObs.nc'
    elif PLCG_type== 'track_smth': file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +cont + plev_str + '_allObs_track-smth-'+smooth_param+'.nc'

    if os.path.isfile(file): ds2= xr.open_dataset(file)
    else:
        if PLCG_type== 'stearing_flow': ds2= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields/" +cont + '_all_levs_allObs.nc')
        elif PLCG_type== 'track_smth': ds2= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +cont + '_all_levs_allObs_track-smth-'+smooth_param+'.nc')

        ds2= ds2.sel(plev= plevel_cont)

    ds2= ds2.rename({cont: cont_full})

    ds= xr.merge([ds, ds2])
    
    
    """add the node to ds - only for allObs"""
    ds_node_var= np.zeros(len(ds.time))
    ds_node_var[:]= np.NaN #makes an empty nan array- so locations without a node are placed with nans
    
    for s in range(len(df_nodes)):
        df_nodes_s= df_nodes.iloc[s]
        ds_index= np.argwhere(np.logical_and(ds.PLnr.values== df_nodes_s.PLnr, ds.time.values== np.datetime64(df_nodes_s.name[1])))[0][0]
        ds_node_var[ds_index]= df_nodes_s.node
    
    ds['node']= (('time'), ds_node_var)
    ds_Obs_evol= ds
        

    if ano_var == True:
        ds_Obs_evol[var_full]= (('time', 'x', 'y'), ds_Obs_evol[var_full].values- np.mean(ds_Obs_evol[var_full].values, axis= (1,2))[:, np.newaxis, np.newaxis])
    if ano_cont == True:
        ds_Obs_evol[cont_full]= (('time', 'x', 'y'), ds_Obs_evol[cont_full].values- np.mean(ds_Obs_evol[cont_full].values, axis= (1,2))[:, np.newaxis, np.newaxis])
        


#if mirror_top_bottom:
#    ds_Obs_evol['node']= (('time'), (y-1 - ds_Obs_evol.node.values//x) *x + ds_Obs_evol.node.values%x )


"""plot the SOMs"""
for iSOM in range(1, x*y +1):
    
    ds_node= ds_Obs_evol.where(ds_Obs_evol.node == iSOM, drop=True)
    
    ax= axs[iSOM-1]
    
    cf= ax.contourf(ds_node.x, ds_node.y, np.mean(ds_node[var_full], axis= 0), levels= levels, cmap= cmap, extend='both')      
    
    cs= ax.contour(ds_node.x, ds_node.y, np.mean(ds_node[cont_full], axis= 0), levels= cont_levels, colors='k')#, linestyles='solid')

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

cbar_ax = fig.add_axes([0.05, 0.065, 0.4, 0.015])
cb= fig.colorbar(cf, cax=cbar_ax, orientation="horizontal")

if 't_ano' in var_full: labelvar= 'Temperature anomaly [K]'
elif ano_var== False: labelvar= ds[var_full].long_name + ' ['+ ds[var_full].units+ ']'
else: labelvar= ''

cb.set_label(labelvar, size=10) 


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


"""make the propagation arrow"""
arrow_ax = fig.add_axes([0.82, 0.03, 0.12, 0.045])
import matplotlib.patches as mpatches

arrow = mpatches.FancyArrowPatch((0.2, .8), (.7, .75), mutation_scale=20)
arrow_ax.add_patch(arrow)
arrow_ax.text(0,-0.2 , 'Propagation\ndirection', fontweight= 'bold', fontsize=8)
arrow_ax.set_frame_on(False)
arrow_ax.get_xaxis().set_visible(False)
arrow_ax.get_yaxis().set_visible(False)






if save:
    if PLCG_type== 'stearing_flow':
        save_name=savedir+ 'SOM_evolution_'+Svar_full+"_"+str(Splevel)+"_"+Obs+'--mean_'+var_full+"_"+cont_full+"_"+Obs_evol_str+ save_extra+"_x"+str(x)+"_y"+str(y)
    elif PLCG_type== 'track_smth':
        save_name=savedir+ 'SOM_evolution_track-smth'+smooth_param+'_'+Svar_full+"_"+str(Splevel)+"_"+Obs+'--mean_'+var_full+"_"+cont_full+"_"+Obs_evol_str+ save_extra+"_x"+str(x)+"_y"+str(y)

    print(save_name)
    plt.savefig(save_name , dpi= 140)#, bbox_inches='tight')


print('could include steps towards nan - so steps to no propagation speed')