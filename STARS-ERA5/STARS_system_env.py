#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 10:01:45 2019

@author: pst019
"""

save=True

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

import pandas as pd
import xarray as xr

import matplotlib.pyplot as plt
import numpy as np
from f_useful import *
from f_STARS import *

fignr= 1


#imp_lists=True
imp_lists=False

if imp_lists:
    STARS = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv")
    STARS = merge_ERA5(STARS, Mediadir, "PL/STARS/STARS_ERA5.csv")
    
    
    STARS_ind= STARS_individual_systems(STARS)
    
    """get individual systems"""
    
    test= 'version4'
    dist= 150
    
    Stoll_STARS_file="ERA5_STARS/STARS-PMC-merge/"+test+"_dist_"+str(dist)+"_handfix.csv"
    Stoll = pd.read_csv(Mediadir+ Stoll_STARS_file)
    Stoll['time']= pd.DatetimeIndex(Stoll.time)
    
    Stoll= Stoll_excl(Stoll)
    Stoll['ID']= Stoll['Stoll nr']
    Stoll['Obs']= Stoll['row_idx'] +1
    
    """Stoll_individual_systems"""
    Stoll_ind= Stoll_individual_systems(Stoll)




"""-----------------------------------------------------------------------------"""
"""set the dataset"""
#Sx, Stype = STARS_ind, 'system'
#S_total=STARS

Sx, Stype = Stoll_ind, 'system'
S_total= Stoll


#Sx, Stype= STARS, 'step'
#Sx, Stype= Stoll, 'Stoll step'


"""set variables"""
save=True
#save=False



#var, file_var, vmin, vmax= 't', 't', 255, 262
#var, file_var, vmin, vmax= 't_ano', 't', -3, 3

var, file_var, vmin, vmax= 'z', 'z', 1180, 1320
#var, file_var, vmin, vmax= 'z_ano', 'z', -100, 100

#var, file_var, vmin, vmax= 'q', 'q', 0.6, 1.8
#var, file_var, vmin, vmax= 'vo', 'vo', 0, 30
var, file_var, vmin, vmax= 'U', 'U', 4, 16


plevel= 850
ano=False


sym= False
#cmap='jet'
#cmap= 'viridis'
#cmap= 'Greys'
#cmap= 'Greys_r'
#cmap= 'Blues'
#cmap= 'Reds'

if "ano" in var: sym= True
cmap= 'RdBu_r'


"""-----------------------------------------------------------------------------"""
"""preprocess STARS"""
#S= remove_morphs(S, 'Morphology', count= 50, how='replace', excludelist=['-', 'U'])
if Stype=='step':
    Sx= remove_morphs(Sx, 'Morphology', count= 100, how='remove', excludelist=['-', 'U'])

if Stype=='Stoll step':
    Sdist= 100
    #only the steps where the cloud morphology is observed and the distance is lower than Sdist
    Sx= Sx.loc[Sx['STARS lat'] > 1]
    Sx= Sx.loc[distance((Sx.lat, Sx.lon), (Sx['STARS lat'], Sx['STARS lon'])) < Sdist]
    Sx= remove_morphs(Sx, 'Morphology', count= 25, how='remove', excludelist=['-', 'U'])



"""plot morphology"""
fig= plt.figure(fignr)
ax= fig.add_subplot(111)

fignr+=1
plt.clf()


if Stype =='system':
    Sx['Morphology']= Sx['Morphology'].replace({'': 'others'})

    morphlist= list(Sx.groupby('Morphology').size().sort_values(ascending=False).index)

    
#    for morph in pd.crosstab(S_ind['Month'], S_ind.PL_Morph).columns:
    for morph in morphlist:
        print(morph)
        color = next(ax._get_lines.prop_cycler)['color']
        S_morph= Sx[Sx['Morphology'] == morph]
    #    print('Morph:', morph)
    
        S_morph_count= S_morph['PL_Morph_full'].value_counts()
        
        l= len(S_morph_count)
        if l == 1:
            plt.bar(morph, len(S_morph), color= color, edgecolor= 'k', lw= 1)
    
        if l > 1:
            count= 0
            for li in range(l):
    #            print(S_morph_count.index[li])
                b= plt.bar(morph, S_morph_count[li], bottom= count,  color= color, edgecolor= 'k', lw= 1)#hatch= patterns[li],
                
                            
                if S_morph_count[li] > 4:
                    plt.text(morph, count+ S_morph_count[li]/2 -.3,  S_morph_count.index[li], ha='center', va='center')
    
                count+= S_morph_count[li]
    
    plt.ylabel('Number of polar lows')


if Stype == 'step':
#    morphs= S.groupby('Morphology').size()
    morphs= Sx.groupby('Morphology').size().sort_values(ascending=False)
    print(morphs)
 
    morphs.plot(kind= 'bar')
    #S = Spiraliform, C = Comma shaped, MGR = Merry-go-round, W = Wave system, U = Undefined, T = Transition between different forms, H = Hybrid, - = Systems don't appear completely on imagery
    plt.ylabel('Number of timesteps')
    plt.tight_layout()



if save:
    savedir= homedir+ 'Polar_Low/STARS-ana/Figs/STARS-Stoll/Sub-Morph_'+Stype
    print(savedir)
    plt.savefig(savedir) 


"""-----------------------------------------------------------------------------"""
"""import the PL centred field"""

Obsnr=1
file= Mediadir + "ERA5_STARS/PL_centred_fields/"+file_var+"_"+str(plevel)+"_Obsnr"+str(Obsnr)+".nc"
ds= xr.open_dataset(file)
#dsPLnr= [str(int(ds.PLnr[n].values))[:-3]+'_'+str(int(ds.PLnr[n].values))[-1] for n in range(len(ds.PLnr))]
ds['PL_nr_str']= (('time'), [str(int(ds.PLnr[n].values))[:-3]+'_'+str(int(ds.PLnr[n].values))[-1] for n in range(len(ds.PLnr))] )


fig= plt.figure(fignr, figsize=(9, 6.5))
fignr+=1
plt.clf()

for mi, morph in enumerate(morphlist):
    if Stype=='system':
        Sx_morph_now= Sx[Sx['Morphology'] == morph].index.values

        common= common_list(ds.PL_nr_str.values, Sx_morph_now) #gets the common list
        ds_now= ds.where(ds.PL_nr_str.isin(common), drop=True )
        
    ax1= plt.subplot(2, 3, mi+1)    

#    if sym:
#        cf= plt.contourf(ds.x, ds.y, ds_now.mean(dim='time')[var], cmap= cmap)#, vmin= -vextr, vmax= vextr)
#        cs= plt.contour(ds.x, ds.y, ds_now.mean(dim='time')[var], colors='k')#, vmin= -vextr, vmax= vextr, linewidth= 1)
#    else:
    cf= plt.contourf(ds.x, ds.y, ds_now.mean(dim='time')[var], cmap= cmap, vmin= vmin, vmax= vmax)
    cs= plt.contour(ds.x, ds.y, ds_now.mean(dim='time')[var], colors='k', vmin= vmin, vmax= vmax, linewidth= 1)

    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.1f')


    plt.title(morph)
        

plt.tight_layout()
fig.subplots_adjust(bottom=0.15)
cbar_ax = fig.add_axes([0.2, 0.08, 0.6, 0.02])
cb= fig.colorbar(cf, cax=cbar_ax, orientation="horizontal")   



if var== 'z_ano': labelvar= 'Geopotential height anomaly [m]'
elif var== 'z': labelvar= 'Geopotential height [m]'

elif var== 't_ano': labelvar= 'Temperature anomaly [K]'
elif var== 't': labelvar= 'Temperature [K]'

elif var== 'U_ano': labelvar= 'Wind speed anomaly [m/s]'
elif var== 'U': labelvar= 'Wind speed [m/s]'

elif var== 'q_ano': labelvar= 'Specific humidity anomaly [g(kg)]'
elif var== 'q': labelvar= 'Specific humidity [g/kg]'
else: labelvar= var

cb.set_label(labelvar, size=14) 

if save:
    savedir= homedir+ 'Polar_Low/STARS-ana/Figs/STARS-Stoll/Sub-Morph_'+Stype+'_env_'+var
    print(savedir)
    plt.savefig(savedir) 


#"""------------ plot the figure------------"""   
#if var != 'Prop_dir':
#    plt.figure(fignr)
#    fignr+=1
#    plt.clf() 
#     
#    if var not in ['Season', 'Month']: ax1= plt.subplot(211)
#      
#    bot=False
#    
#    morphlist= list(Sx.groupby('Morphology').size().sort_values(ascending=False).index) #the morphologies sorted by size
#    for morph in morphlist:
#        botn,b_dn,p_dn= plt.hist(Sx[Sx['Morphology'] == morph][var], bins= bins , label=morph, bottom= bot, edgecolor= 'k', lw= 1)
#        bot += botn
#    #
#    plt.xlim(bins[0], bins[-1])
#    if Stype == 'step':
#        plt.ylabel('Number of timesteps')
#    else:
#        plt.ylabel('Number of polar lows')
#        
#    plt.legend()
#    
#    
#    
#    if var== 'Season': plt.xticks(bins[1::2]+.5)
#    
#    
#    if var not in ['Season', 'Month']: ax1.xaxis.set_ticks_position('both')
#    
#    
#    """-------------- plot density distributions ------------------"""
#    if var not in ['Season', 'Month']:
#        ax2= plt.subplot(212, sharex= ax1)
#        
#        from scipy import stats
#        
#        xplot= np.linspace(bins[0], bins[-1], 100)
#        
##        kernel= stats.gaussian_kde(Sx[var].dropna())
##        plt.plot(xplot, kernel(xplot), label='all', lw= 2)
#        
#        for morph in morphlist:
#            kernel= stats.gaussian_kde(Sx[Sx['Morphology'] == morph][var].dropna())
#            plt.plot(xplot, kernel(xplot), label=morph)
#        
#        plt.legend()       
##        plt.legend(ncol= 2)
#        plt.ylabel('Normalized density')
#        ax2.xaxis.set_ticks_position('both')
#        
#    
#        
#    plt.xlabel(var)
#    
#    
#    if var== 'multiple': plt.xlabel('Number of simultaneous systems')
#    elif var=='Duration': plt.xlabel('Life time [h]')
#    elif var=='Distance': plt.xlabel('Distance between start and end point [km]')
#    elif var=='Prop_speed': plt.xlabel('Propagation speed [km/h]')
#    elif var=='Diameter': plt.xlabel('Diameter of the polar low [km]')
#    elif var=='Diameter_max': plt.xlabel('Maximum diameter of the polar low [km]')
#    elif var=='Diameter_min': plt.xlabel('Minimum diameter in mature stage of the polar low [km]')
#    elif var=='U_400': plt.xlabel('Maximum wind speed within 400 km radius [m/s]')
#    elif 'skt-t500' in var: plt.xlabel('Max skin temperature - T$_{500}$ within 200 km [K]')
#    elif 'skt-t700' in var: plt.xlabel('Max skin temperature - T$_{700}$ within 200 km [K]')
#    elif 'skt' in var: plt.xlabel('Median skin temperature within 200 km [K]')
#    elif 'grad_t' in var: plt.xlabel('Max temperature gradient within 300 km [K/100km]')
#    elif 'msl' in var: plt.xlabel('Minimum sea level pressure within 50 km [hPa]')
#    elif 'tp' in var: plt.xlabel('Mean precipitation within 300 km [mm/h]')
#    elif 'cp' in var: plt.xlabel('Mean convective precipitation within 300 km [mm/h]')
#    elif 'sf' in var: plt.xlabel('Mean snowfall within 300 km [mm/h]')
#    elif 'sshf' in var: plt.xlabel('Mean surface sensible heat flux within 300 km [W/m$^2$]')
#    elif 'slhf' in var: plt.xlabel('Mean surface latent heat flux within 300 km [W/m$^2$]')
#    elif 'barotropic' in var: plt.xlabel('Barotropic e-folding time [h]')
#    elif 'baroclinic' in var: plt.xlabel('Baroclinic e-folding time [h]')
#    
#    elif var[:2] == 'vo' and var != 'vortex_type': plt.xlabel('Max vorticity at 850hPa within 200 km [10$^{-5}$ 1/s]')
#    
#    
#    plt.tight_layout()
#
#
#"""prop direction - wind rose"""
#if var== 'Prop_dir':
#    fignr+=1
#
#    fig= plt.figure(fignr, figsize=(9.5, 6))
#    fignr+=1
#    plt.clf()
#    
#    from f_meteo import *
#    from windrose import WindroseAxes #pip install git+https://github.com/python-windrose/windrose
#    import matplotlib.gridspec as gridspec
#    import matplotlib.ticker as tkr
#    
#    gs = gridspec.GridSpec(2, 3)     # (nblines, nbcol)
#    # return lists of bottom and top position of rows, left and right positions of columns.
#
#    bottom, top, left, right = gs.get_grid_positions(fig)  # [bottom, top, left, right]
#
#    Sx['Prop_dir']= [calculate_initial_compass_bearing((Sx['start_lat'][i], Sx['start_lon'][i]), (Sx['end_lat'][i], Sx['end_lon'][i]) ) for i in range(len(Sx))]
#
#    col, row= 0,0
#    
#    rect = [left[col],  bottom[row],  right[col]-left[col],  0.9*(top[row]-bottom[row])]
#    
#    ax = WindroseAxes(fig, rect)
#    fig.add_axes(ax)
#
#    ax.bar(Sx['Prop_dir'], Sx['Prop_speed'], normed=False, bins=np.arange(0, 61, 20), opening=0.9, edgecolor='white')
#    ax.set_title('All', position=(0.1, 1.01), color= 'r', fontsize= 15)
##    ax.set_yticks(np.arange(0, 15.5, 5))
#    ax.yaxis.set_major_formatter(tkr.FormatStrFormatter('%2.0f'))
#    
#    ax.set_legend()
#    ax.legend(title="Speed [km/h]", loc= (-.3,0), decimal_places=0)
#
#    
#    for mi, morph in enumerate(morphlist[1:]):
#        col, row= (mi+1)%3, (mi+1)//3
##        print(mi, morph, col, row)
#        rect = [left[col],  bottom[row],  right[col]-left[col],  0.9*(top[row]-bottom[row])]
#
#        ax = WindroseAxes(fig, rect)
#        fig.add_axes(ax)
#    
#        ax.bar(Sx[Sx['Morphology'] == morph]['Prop_dir'], Sx[Sx['Morphology'] == morph]['Prop_speed'], normed=False, bins=np.arange(0, 61, 20),  opening=0.9, edgecolor='white')
#        ax.set_title(morph, position=(0.1, 1.01), color= 'r', fontsize= 15)
##        ax.set_yticks(np.arange(0, 15.5, 5))
#        ax.yaxis.set_major_formatter(tkr.FormatStrFormatter('%2.0f'))    
#    
#
#
#
#    
#
#
#
#if save:
#    if Stype== 'system': Stype += '_'+system_char
#    
#    savedir= homedir+ 'Polar_Low/STARS-ana/Figs/2/'+'Morph_'+Stype+'_'+var
#    print(savedir)
#    plt.savefig(savedir, bbox_inches='tight')
#
#
#
#


