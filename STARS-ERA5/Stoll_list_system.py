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
    
    Stoll_STARS_file="ERA5_STARS/STARS-PMC-merge/"+test+"_dist_"+str(dist)+".csv"
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

#system_char= 'mean'
#system_char= 'max'
#system_char='min'
#system_char='med'
#system_char='initial'
#system_char='final'




#Sx, Stype= STARS, 'step'
#Sx, Stype= Stoll, 'Stoll step'


"""set variables"""
#save=True
save=False

#var, bins='Season', np.arange(1999, 2020,1)-.5
#var, bins='Month', np.arange(.5, 12.6, 1)
#var, bins= 'multiple', np.arange(10)+.5
#var, bins= 'Duration', np.arange(0, 116, 5)
#var, bins= 'Distance', np.arange(0,2101, 100)  #"""travelled distance - start-end"""
#var, bins= 'Prop_speed', np.arange(0, 100.1, 4)
#var, bins= 'Diameter', np.arange(0, 851, 50)
#var, bins= 'Diameter_max', np.arange(0, 851, 50) #one is 1200
#var, bins= 'Diameter_min', np.arange(0, 651, 50)
##
#var, bins= 'U_400', np.arange(5, 30.1, 1)
#var, bins= 'skt_med_200', np.arange(258, 286, 1)
#var, bins= 'msl_min_50', np.arange(950, 1030.1, 2)
#var, bins= 'blh_max_250', np.arange(1000, 3201.1, 100)
#var, bins= 'cape_max_250', np.arange(0, 301, 10)
#var, bins= 'sshf_mean_300', np.arange(0, 301, 10)
#var, bins= 'slhf_mean_300', np.arange(0, 301, 10)
#var, bins= 'tp_mean_300', np.arange(0, .65, .02) #total precip
#var, bins= 'cp_mean_300', np.arange(0, .65, .02) #convective precip
#var, bins= 'sf_mean_300', np.arange(0, .65, .02) #snow fall
#
#var, bins= 'vo850_max_200', np.arange(0, 120.1, 4)

#var, bins= 'skt-t500_max_200', np.arange(34, 60)
#var, bins= 'skt-t700_max_200', np.arange(20, 40.1)
#var, bins= 'grad_t850_max_300', np.arange(0, 20.1)

#var, bins= 'vert_shear_angle[925, 700]_mean_200', np.arange(0, 180.1, 10)
#
#var, bins= 'barotropic850_max_200', np.arange(0, 25.1, 1) #it is the barotropic growth rate, which is converted to the efolding time later in the code
#var, bins= 'barotropic_efolding_filter4_850_min_200', np.arange(0, 30.1, 1)
#var, bins= 'barotropic_efolding_filter10_850_min_200', np.arange(0, 80.1, 1)

#var, bins= 'baroclinic_efolding_filter4_925-700_min_200', np.arange(0, 20.1, 1)


#var, bins= 'lat', np.arange(51, 82.1, 1) #55, 82
#var, bins= 'start_lat', np.arange(51, 82.1, 1) #58, 82
#var, bins= 'end_lat', np.arange(51, 82.1, 1)
#var, bins= 'lon', np.arange(-30, 70.1, 5) #-30, 55.1
#var, bins= 'start_lon', np.arange(-30, 70.1, 5) #-35, 55.1
#var, bins= 'end_lon', np.arange(-30, 70.1, 5) #-30, 85.1

#var= 'Prop_dir'


#variable for stoll
#var, bins= 'vo', np.arange(0, 70.1, 1)
#var, bins= 'Vort_diameter', np.arange(0, 1101, 50) #from vorticity area
#var, bins= 'slp', np.arange(950, 1030.1, 2)
var, bins= 'vortex_type', np.arange(0, 4.1, 1)




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
    savedir= homedir+ 'Polar_Low/STARS-ana/Figs/2/Sub-Morph_'+Stype
    print(savedir)
    plt.savefig(savedir) 

"""-----------------------------------------------------------------------------"""
"""create variable"""


if var=='Season':
    Sx['Season']= Sx['Season'].astype(int)

if var=='Month':
    Sx['Month']= Sx['Month'].round().astype(int)

if var=='multiple':
    Sx['PL_ID']= Sx['ID'].str.replace(r'\D', '').astype(int)
    Sx['multiple']= Sx.groupby('PL_ID')['PL_ID'].transform('count')

if var in ['Duration', 'Prop_speed', 'Prop_dir']:
    #Sx['start_time']
    starttime= S_total.loc[S_total['Obs'] == 1].groupby(['ID']).first()['time']
    Sx= pd.concat([Sx, starttime], axis= 1, sort=True)
    Sx= Sx.rename(columns={"time":"starttime"})
    
    #endtime
    endobs= Sx['Obs']*2 -1
    endobs.loc[(endobs - endobs.astype(int)) > .01]= np.nan
    
    #endobs.index[0], endobs[0]
    #S[S['ID'] == endobs.index[0]]
    endtime= S_total[np.logical_and(S_total['ID'] == endobs.index[0], S_total['Obs'] == endobs[0])][['ID', 'time']]
    for i in range(1, len(endobs)):
        a= S_total[np.logical_and(S_total['ID'] == endobs.index[i], S_total['Obs'] == endobs[i])][['ID', 'time']]
        endtime= pd.concat([endtime, a])
    
    #endtime.rename(columns={"time":"endtime"})
    endtime= endtime.set_index('ID')
    
    Sx= pd.concat([Sx, endtime], axis= 1, sort=True)
    Sx= Sx.rename(columns={"time":"endtime"})   
    
    Sx['Duration']= np.round((Sx['endtime'] - Sx['starttime'])/np.timedelta64(1,'h'))

if var in ['Distance', 'Prop_speed', 'Prop_dir', 'start_lat', 'start_lon', 'end_lat', 'end_lon']:
    Sx= Sx.rename(columns={"lat":"avg_lat", "lon":"avg_lon"})
    
    startloc= S_total.loc[S_total['Obs'] == 1].groupby(['ID']).first()[['lat', 'lon']]
    Sx= pd.concat([Sx, startloc], axis= 1, sort=True)
    Sx= Sx.rename(columns={"lat":"start_lat", "lon":"start_lon"})
    
    endobs= Sx['Obs']*2 -1
    endobs.loc[(endobs - endobs.astype(int)) > .01]= np.nan
    
    endloc= S_total[np.logical_and(S_total['ID'] == endobs.index[0], S_total['Obs'] == endobs[0])][['ID', 'lat', 'lon']]
    for i in range(1, len(endobs)):
        a= S_total[np.logical_and(S_total['ID'] == endobs.index[i], S_total['Obs'] == endobs[i])][['ID', 'lat', 'lon']]
        endloc= pd.concat([endloc, a])
    
    #endtime.rename(columns={"time":"endtime"})
    endloc= endloc.set_index('ID')
    
    Sx= pd.concat([Sx, endloc], axis= 1, sort=True)
    Sx= Sx.rename(columns={"lat":"end_lat", "lon":"end_lon"})
    
    Sx['Distance']= 110* np.sqrt( (Sx['start_lat']- Sx['end_lat'])**2+ (np.cos(np.deg2rad(0.5*(Sx['start_lat']+ Sx['end_lat'])))* (Sx['start_lon']- Sx['end_lon']))**2)


if var in ['Prop_speed', 'Prop_dir']: Sx['Prop_speed']= Sx['Distance']/Sx['Duration']
    
    
if var == 'Vort_diameter':
    Sx['Vort_diameter']= 2* np.sqrt(Sx['area']) *4/np.pi
    S_total['Vort_diameter']= 2* np.sqrt(S_total['area']) *4/np.pi


if Stype == 'system':
#    Sx.reset_index(inplace=True, drop=True)
    if system_char== 'max': Sx[var]= S_total.groupby(['ID']).max()[var]
    elif system_char== 'min': Sx[var]= S_total.groupby(['ID']).min()[var]
    elif system_char== 'median': Sx[var]= S_total.groupby(['ID']).median()[var]
    elif system_char== 'initial': Sx[var]=S_total.groupby(['ID']).first()[var]
    elif system_char== 'final': Sx[var]=S_total.groupby(['ID']).last()[var]


        
#if Stype=='system': #to remove systems without category - this is done for "steps" earlier
#    Sx= remove_morphs(Sx, 'Morphology', count= 1, how='remove', excludelist=[''])
#S.loc[S['Obs'] == 1].groupby(['ID']).first()['time']




if 'msl_' in var: Sx[var]/= 100
if var[:2]=='vo' and var != 'vortex_type': Sx[var]*= 1E5
if var in ['barotropic850_max_200']: Sx[var]= 1/Sx[var]* 1/(60**2)




"""------------ plot the figure------------"""   
if var != 'Prop_dir':
    plt.figure(fignr)
    fignr+=1
    plt.clf() 
     
    if var not in ['Season', 'Month']: ax1= plt.subplot(211)
      
    bot=False
    
    morphlist= list(Sx.groupby('Morphology').size().sort_values(ascending=False).index) #the morphologies sorted by size
    for morph in morphlist:
        botn,b_dn,p_dn= plt.hist(Sx[Sx['Morphology'] == morph][var], bins= bins , label=morph, bottom= bot, edgecolor= 'k', lw= 1)
        bot += botn
    #
    plt.xlim(bins[0], bins[-1])
    if Stype == 'step':
        plt.ylabel('Number of timesteps')
    else:
        plt.ylabel('Number of polar lows')
        
    plt.legend()
    
    
    
    if var== 'Season': plt.xticks(bins[1::2]+.5)
    
    
    if var not in ['Season', 'Month']: ax1.xaxis.set_ticks_position('both')
    
    
    """-------------- plot density distributions ------------------"""
    if var not in ['Season', 'Month']:
        ax2= plt.subplot(212, sharex= ax1)
        
        from scipy import stats
        
        xplot= np.linspace(bins[0], bins[-1], 100)
        
#        kernel= stats.gaussian_kde(Sx[var].dropna())
#        plt.plot(xplot, kernel(xplot), label='all', lw= 2)
        
        for morph in morphlist:
            kernel= stats.gaussian_kde(Sx[Sx['Morphology'] == morph][var].dropna())
            plt.plot(xplot, kernel(xplot), label=morph)
        
        plt.legend()       
#        plt.legend(ncol= 2)
        plt.ylabel('Normalized density')
        ax2.xaxis.set_ticks_position('both')
        
    
        
    plt.xlabel(var)
    
    
    if var== 'multiple': plt.xlabel('Number of simultaneous systems')
    elif var=='Duration': plt.xlabel('Life time [h]')
    elif var=='Distance': plt.xlabel('Distance between start and end point [km]')
    elif var=='Prop_speed': plt.xlabel('Propagation speed [km/h]')
    elif var=='Diameter': plt.xlabel('Diameter of the polar low [km]')
    elif var=='Diameter_max': plt.xlabel('Maximum diameter of the polar low [km]')
    elif var=='Diameter_min': plt.xlabel('Minimum diameter in mature stage of the polar low [km]')
    elif var=='U_400': plt.xlabel('Maximum wind speed within 400 km radius [m/s]')
    elif 'skt-t500' in var: plt.xlabel('Max skin temperature - T$_{500}$ within 200 km [K]')
    elif 'skt-t700' in var: plt.xlabel('Max skin temperature - T$_{700}$ within 200 km [K]')
    elif 'skt' in var: plt.xlabel('Median skin temperature within 200 km [K]')
    elif 'grad_t' in var: plt.xlabel('Max temperature gradient within 300 km [K/100km]')
    elif 'msl' in var: plt.xlabel('Minimum sea level pressure within 50 km [hPa]')
    elif 'tp' in var: plt.xlabel('Mean precipitation within 300 km [mm/h]')
    elif 'cp' in var: plt.xlabel('Mean convective precipitation within 300 km [mm/h]')
    elif 'sf' in var: plt.xlabel('Mean snowfall within 300 km [mm/h]')
    elif 'sshf' in var: plt.xlabel('Mean surface sensible heat flux within 300 km [W/m$^2$]')
    elif 'slhf' in var: plt.xlabel('Mean surface latent heat flux within 300 km [W/m$^2$]')
    elif 'barotropic' in var: plt.xlabel('Barotropic e-folding time [h]')
    elif 'baroclinic' in var: plt.xlabel('Baroclinic e-folding time [h]')
    
    elif var[:2] == 'vo' and var != 'vortex_type': plt.xlabel('Max vorticity at 850hPa within 200 km [10$^{-5}$ 1/s]')
    
    
    plt.tight_layout()


"""prop direction - wind rose"""
if var== 'Prop_dir':
    fignr+=1

    fig= plt.figure(fignr, figsize=(9.5, 6))
    fignr+=1
    plt.clf()
    
    from f_meteo import *
    from windrose import WindroseAxes #pip install git+https://github.com/python-windrose/windrose
    import matplotlib.gridspec as gridspec
    import matplotlib.ticker as tkr
    
    gs = gridspec.GridSpec(2, 3)     # (nblines, nbcol)
    # return lists of bottom and top position of rows, left and right positions of columns.

    bottom, top, left, right = gs.get_grid_positions(fig)  # [bottom, top, left, right]

    Sx['Prop_dir']= [calculate_initial_compass_bearing((Sx['start_lat'][i], Sx['start_lon'][i]), (Sx['end_lat'][i], Sx['end_lon'][i]) ) for i in range(len(Sx))]

    col, row= 0,0
    
    rect = [left[col],  bottom[row],  right[col]-left[col],  0.9*(top[row]-bottom[row])]
    
    ax = WindroseAxes(fig, rect)
    fig.add_axes(ax)

    ax.bar(Sx['Prop_dir'], Sx['Prop_speed'], normed=False, bins=np.arange(0, 61, 20), opening=0.9, edgecolor='white')
    ax.set_title('All', position=(0.1, 1.01), color= 'r', fontsize= 15)
#    ax.set_yticks(np.arange(0, 15.5, 5))
    ax.yaxis.set_major_formatter(tkr.FormatStrFormatter('%2.0f'))
    
    ax.set_legend()
    ax.legend(title="Speed [km/h]", loc= (-.3,0), decimal_places=0)

    
    for mi, morph in enumerate(morphlist[1:]):
        col, row= (mi+1)%3, (mi+1)//3
#        print(mi, morph, col, row)
        rect = [left[col],  bottom[row],  right[col]-left[col],  0.9*(top[row]-bottom[row])]

        ax = WindroseAxes(fig, rect)
        fig.add_axes(ax)
    
        ax.bar(Sx[Sx['Morphology'] == morph]['Prop_dir'], Sx[Sx['Morphology'] == morph]['Prop_speed'], normed=False, bins=np.arange(0, 61, 20),  opening=0.9, edgecolor='white')
        ax.set_title(morph, position=(0.1, 1.01), color= 'r', fontsize= 15)
#        ax.set_yticks(np.arange(0, 15.5, 5))
        ax.yaxis.set_major_formatter(tkr.FormatStrFormatter('%2.0f'))    
    



    



if save:
    if Stype== 'system': Stype += '_'+system_char
    
    savedir= homedir+ 'Polar_Low/STARS-ana/Figs/2/'+'Morph_'+Stype+'_'+var
    print(savedir)
    plt.savefig(savedir, bbox_inches='tight')






