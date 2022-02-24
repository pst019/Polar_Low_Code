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

elif typ == 'Bias_cor':
    Stoll.loc[np.logical_and(Stoll[dir_var] < 5*np.pi/8, Stoll[dir_var] >= 1*np.pi/8), 'shear_category'] = 'forward'
    Stoll.loc[np.logical_and(Stoll[dir_var] < 9*np.pi/8, Stoll[dir_var] >= 5*np.pi/8), 'shear_category'] = 'right'
    Stoll.loc[np.logical_and(Stoll[dir_var] < 13*np.pi/8, Stoll[dir_var] >= 9*np.pi/8), 'shear_category'] = 'reverse'
    Stoll.loc[np.logical_or(Stoll[dir_var] < 1*np.pi/8, Stoll[dir_var] >= 13*np.pi/8), 'shear_category'] = 'left'
    Stoll.loc[Stoll[speed_var] <= strength_level, 'shear_category'] = 'weak'


bins= np.arange(0, 5, 1)
title= "Shear strength\n[m$\cdot$s$^{-1}$ (km)$^{-1}$]"

typname= typ+ str(l1)+'-'+str(l2)+'_mean-'+str(mean)

#S_shear= Stoll[Stoll[speed_var] > 1.5]
#S_shear= S_shear[np.logical_and(S_shear[dir_var] < 7*np.pi/4, S_shear[dir_var] >= 5*np.pi/4)]





#good solution to put all in same: https://stackoverflow.com/questions/57027970/subplots-in-windrose-diagram


"""plot all points and distribution"""
if all_figs:
    plt.figure(fignr, figsize= (6,6) )
    fignr+= 1
    plt.clf()
    
    xyall = np.vstack([np.sin(Stoll[dir_var])*Stoll[speed_var], np.cos(Stoll[dir_var])*Stoll[speed_var]]) 
    xall, yall= xyall
    plt.scatter(xall,yall, zorder= 0, s= .2)
    plt.title('all PL points')
    
    
    nbins= 100
    k = stats.gaussian_kde(xyall) #, bw_method= 'silverman')
    xi, yi = np.mgrid[xall.min():xall.max():nbins*1j, yall.min():yall.max():nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    
    lev= np.percentile(zi, 95) #this makes 5% of the points with highest rate 
    #- maybe better to find a level such that 70% of points are included
    #    print(sum(zi[zi> lev])/sum(zi))
    
    zi= zi.reshape(xi.shape)
    
    cset= plt.contour(xi, yi, zi, zorder =1) #, levels=[lev], colors= color)#,  label=str(iSOM))
    plt.contour(xi, yi, zi, levels=[lev], colors= 'r', zorder= 2)#,  label=str(iSOM))
    plt.clabel(cset, inline=1, fontsize=10)
    
    
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    theta = np.linspace(-np.pi, np.pi, 120)
    circlex = np.cos(theta)
    circley = np.sin(theta)
    plt.plot(1*circlex, 1*circley, c= 'white')
    plt.plot(2*circlex, 2*circley, c= 'white')
    plt.plot(0,0, 'o', color= 'white')
    
    
    
    
    print('plot the SOMs')
    plt.figure(fignr, figsize= (6,6) )
    fignr+= 1
    plt.clf()
    ax = plt.subplot(111)#, polar=True)
    
    #fignr+=1
    plt.clf()
    for iSOM in [9]: 
    #for iSOM in range(1, xSOM*ySOM +1):
        color = next(ax._get_lines.prop_cycler)['color']
        print('SOM ', iSOM)
        
        S_SOM= Stoll[Stoll.node== iSOM]
    
       
        xy = np.vstack([np.sin(S_SOM[dir_var])*S_SOM[speed_var], np.cos(S_SOM[dir_var])*S_SOM[speed_var]]) 
        x, y= xy
    
        """make the cartesian plot"""
    #    plt.plot(xall,yall, 'x', zorder= -1, alpha= 0.4, c= 'k') #plot all points in the background    
        plt.scatter(xall,yall, zorder= -1, s= .5, c='k')
    
        plt.scatter(x,y, zorder= 0, s= .5)
    #    plt.plot(x,y, 'x', zorder= 0) #points of this SOM
    
    
        plt.xlim(-10, 10)
        plt.ylim(-10, 10) 
        theta = np.linspace(-np.pi, np.pi, 120)
        circlex = np.cos(theta)
        circley = np.sin(theta)
        plt.plot(1*circlex, 1*circley, c= 'white')
        plt.plot(2*circlex, 2*circley, c= 'white')
        plt.plot(0,0, 'o', color= 'white')
        
    
        plt.title('SOM ' +str(iSOM) )
    
        
        """make the contours for the density"""
        nbins= 100
    
        k = stats.gaussian_kde(xy) #, bw_method= 'silverman')
        xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
        zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    
        lev= np.percentile(zi, 95) #this makes 5% of the points with highest rate 
        #- maybe better to find a level such that 70% of points are included
    #    print(sum(zi[zi> lev])/sum(zi))
        
        zi= zi.reshape(xi.shape)
    
        cset= plt.contour(xi, yi, zi, zorder =1) #, levels=[lev], colors= color)#,  label=str(iSOM))
        plt.contour(xi, yi, zi, levels=[lev], colors= 'r', zorder= 2)#,  label=str(iSOM))
    
        plt.clabel(cset, inline=1, fontsize=10)
    
    #    ctf = ax.contourf(S_SOM[dir_var], S_SOM[speed_var], z)
    #    plt.colorbar(ctf)
    
        totalcount= len(x)
        count= 0
        for point in range(totalcount):
            dist= (x[point]- xi)**2 + (y[point]- yi)**2
            if zi[np.unravel_index(dist.argmin(), dist.shape)] >= lev: count += 1
        
        print(count/totalcount)
    
    
    
        """new classification"""
        plt.figure(fignr, figsize= (6,6) )
        fignr+= 1
        plt.clf()
    
    
        plt.plot(xall,yall, 'x', zorder= -1, alpha= 0.4, c= 'k') #plot all points in the background    
    
        S_shear= Stoll.loc[Stoll['shear_category'] == 'left']
    
    #    S_shear= Stoll[Stoll[speed_var] > 1.5]
    #    S_shear= S_shear[np.logical_and(S_shear[dir_var] < 7*np.pi/4, S_shear[dir_var] >= 5*np.pi/4)]
    
        xy = np.vstack([np.sin(S_shear[dir_var])*S_shear[speed_var], np.cos(S_shear[dir_var])*S_shear[speed_var]]) 
        x, y= xy
        plt.plot(x,y, 'x', zorder= 0) #points of this SOM
    
    
        plt.xlim(-10, 10)
        plt.ylim(-10, 10) 
        theta = np.linspace(-np.pi, np.pi, 120)
        circlex = np.cos(theta)
        circley = np.sin(theta)
        plt.plot(1*circlex, 1*circley, c= 'white')
        plt.plot(2*circlex, 2*circley, c= 'white')
        plt.plot(0,0, 'o', color= 'white')
    
    
        """make the polar coordinate plot"""
        plt.figure(fignr)
        fignr+= 1
    
        plt.clf()
    #    plt.plot(S_SOM[dir_var], S_SOM[speed_var], 'x')
        plt.polar(-1* (S_SOM[dir_var] - np.pi/2) , S_SOM[speed_var], 'x')    
        plt.title('SOM ' +str(iSOM) )
        
        fig2= plt.figure(fignr)
        fignr+= 1
    
        plt.clf()
    
    
    
        """make the distribution of the shear vector in wind rose"""
        ax = WindroseAxes.from_ax(fig= fig2)
        ax.set_xticklabels(["FS", "", "Left", "", "RS","",  "Right", ""])
    
    
        axb= ax.bar(np.rad2deg(S_SOM[dir_var]), S_SOM[speed_var], normed=True, bins=bins,  opening=1 , edgecolor='k', cmap=mpl.cm.viridis_r)
        ax.set_rgrids([]) #remove the circular axis
    #    
        ax.set_title('SOM '+str(iSOM), position=(0.1, 1.01), color= 'r')
    #
    #
    #    if iSOM== 1:
    #        fig_legend= plt.figure(fignr+1, figsize= (10,3.5)) # +0) )
    #        plt.clf()
    #
    #        ax2 = WindroseAxes.from_ax(fig= fig_legend)
    #        ax2.bar(S_SOM[dir_var], S_SOM[speed_var], normed=True, bins=bins,  opening=1 , edgecolor='k', cmap=mpl.cm.viridis_r)
    #
    #        ax2.set_legend()
    #        ax2.legend(title=title, loc= (-1.1, 0), decimal_places=0)
    ##        fig_legend.tight_layout()
    ##        fig_legend= plt.figure(fignr+1, figsize= (3.5,3.5)) # +0) )
    ##        fig_legend.legend(axb)#, title="Speed [km$\cdot$h$^{-1}$]", decimal_places=0)
    #    
    #        if save:
    #            save_name=savedir+ typ+'_SOM-legend'
    #           
    #            print(save_name)
    #            fig_legend.savefig(save_name , bbox_inches='tight', dpi= 150)    
    #    
    #
        
    



"""the contour plot for all SOMS"""

fig= plt.figure(fignr, figsize= (6,6))
fignr+= 1

plt.clf()
ax = plt.subplot(111)

#plt.scatter(xall,yall, zorder= -1, s= .5, c='k')


for iSOM in range(1, xSOM*ySOM +1):
    color = next(ax._get_lines.prop_cycler)['color']
    print('SOM ', iSOM)
    
    S_SOM= Stoll[Stoll.node== iSOM]

    
    xy = np.vstack([np.sin(S_SOM[dir_var])*S_SOM[speed_var], np.cos(S_SOM[dir_var])*S_SOM[speed_var]]) 
    x, y= xy
    
    plt.scatter(x,y, zorder= 0, s= .1, c=color)
#    plt.plot(np.median(x), np.median(y), 'o', color= color,  label=str(iSOM))
#    plt.plot(np.mean(x), np.mean(y), 'o', color= color,  label=str(iSOM))
#    plt.scatter(np.mean(x), np.mean(y), color= color, edgecolor='k', label=str(iSOM), zorder=10, s= 80)
    plt.scatter(np.mean(x), np.mean(y), color= color, edgecolor='k', zorder=3, s= 180)
    plt.text(np.mean(x), np.mean(y) -0.05, str(iSOM), horizontalalignment= 'center', verticalalignment= 'center', fontsize=12, fontweight='bold', zorder= 4)

    """make the contour that includes 50% of points"""
#    fraction = 0
#
#    totalcount= len(x)
#    nbins= 100
#
#    k = stats.gaussian_kde(xy) #, bw_method= 'silverman')
#    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
#    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
#    zi= zi.reshape(xi.shape)
#
#    percentile= 97
#    lev= np.percentile(zi, percentile) #this makes 5% of the points with highest rate 
#    
#    while fraction < 0.5:
#        """a method to calculate how many points are included by the contour"""
#        count= 0
#        for point in range(totalcount):
#            dist= (x[point]- xi)**2 + (y[point]- yi)**2
#            if zi[np.unravel_index(dist.argmin(), dist.shape)] >= lev: count += 1
#        
#        fraction= count/totalcount
#        print('points included in contour:', fraction)
#        
#        if fraction < 0.5:
#            percentile =percentile - (0.5 - fraction)*10
#        
#        lev= np.percentile(zi, percentile)
#    
#    print('percentile', percentile)
#    
#
#    plt.contour(xi, yi, zi, levels=[lev], colors= color)



lim= 6
plt.xlim(-lim, lim)
plt.ylim(-lim, lim)


if with_cats:
    theta = np.linspace(-np.pi, np.pi, 120)
    circlex = np.cos(theta)
    circley = np.sin(theta)
    plt.plot(strength_level*circlex, strength_level*circley, c= 'k')
    #plt.plot(2*circlex, 2*circley, c= 'k')  
    
    line_x= np.linspace(np.sqrt(0.5)*strength_level, lim, 2)

    
    if typ == '':
        plt.plot(line_x, line_x, c='k')
        plt.plot(line_x, -line_x, c='k')
        plt.plot(-line_x, -line_x, c='k')
        plt.plot(-line_x, line_x, c='k')
    
    
    elif typ == 'Bias_cor':
        plt.plot(line_x, 2* line_x, c='k')
        plt.plot(2*line_x, -line_x, c='k')
        plt.plot(-line_x, -2* line_x, c='k')
        plt.plot(-2* line_x, line_x, c='k')

    plt.text(4.5, 0, 'forward \n'+ str(np.round(100*len(Stoll[Stoll['shear_category']=='forward'])/len(Stoll), 1))+'%'
             , horizontalalignment= 'center', verticalalignment= 'center', fontsize=14)
    plt.text(0, -5, 'right \n'+ str(np.round(100*len(Stoll[Stoll['shear_category']=='right'])/len(Stoll), 1))+'%'
             , horizontalalignment= 'center', verticalalignment= 'center', fontsize=14)
    plt.text(-4, 0, 'reverse \n'+ str(np.round(100*len(Stoll[Stoll['shear_category']=='reverse'])/len(Stoll), 1))+'%'
             , horizontalalignment= 'center', verticalalignment= 'center', fontsize=14)
    plt.text(0, 5, 'left \n'+ str(np.round(100*len(Stoll[Stoll['shear_category']=='left'])/len(Stoll), 1))+'%'
             , horizontalalignment= 'center', verticalalignment= 'center', fontsize=14)
    plt.text(0, 0.2, 'weak \n'+ str(np.round(100*len(Stoll[Stoll['shear_category']=='weak'])/len(Stoll), 1))+'%'
             , horizontalalignment= 'center', verticalalignment= 'center', fontsize=14)

#plt.legend(ncol=3, loc= 2)

#plt.xlabel('Shear in propagation direction [m/s per km]')
#plt.ylabel('Shear across propagation [m/s per km]')

plt.xlabel('Shear in propagation direction [10$^{-3}$ s$^{-1}$]')
plt.ylabel('Shear across propagation [10$^{-3}$ s$^{-1}$]')



label_shear= shear_type+ '_'+ str(l1) +'-'+str(l2) + 'mean-'+ str(mean)

if save:
    
    save_name=savedir+ 'classification2-'+label_shear
    if with_cats== False:
        save_name+= '_no_categories'
        
    print(save_name)
    fig.savefig(save_name , bbox_inches='tight', dpi= 150)


  
plt.title(label_shear)




