#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 10:39:19 2019

@author: pst019
"""


import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/PatsOrange/'

homedir=Mediadir+'home/'


import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/test_Denis/mc_era5-master/code') #to get obs_tracks_api


import pandas as pd 
import numpy as np
import xarray as xr
import scipy.ndimage.filters as filters

from f_useful import *
from f_meteo import *
from obs_tracks_api import prepare_tracks  #conda install -c dennissergeev octant


import os
user = os.getcwd().split('/')[2]




def import_STARS(Mediadir= '/media/pst019/1692A00D929FEF8B/', file= "PL/STARS/Rojo-etal_2019.csv",
                 droplist= ['Optional_diameter', 'Comment', 'Comment_visual', 'Comment_uncertainties', 'Press_min', 'U_10min_knots', 'U_3sec_kntots']):

    S = pd.read_csv(Mediadir+file, sep=',', header= 27)


    new_cols= list(S.columns)
    new_cols[1] = 'Time'
    new_cols[6:8] = ['Diameter', 'Optional_diameter']
    new_cols[10:16] = ['Comment', 'Comment_visual', 'Comment_uncertainties', 'Press_min', 'U_10min_knots', 'U_3sec_kntots']
    S.columns= new_cols
    
    S.drop(droplist, axis=1, inplace=True )
    
    S['time'] = pd.to_datetime(S['Time'])
    S.drop(['Time'], axis=1, inplace=True)
    
    S= S.rename(columns={"Latitude": "lat", "Longitude": "lon"})
    
#    S['Season']= S['Season'].str.slice(0, 4).astype(int)
    S['Month']= S['time'].dt.month
    #S = S.set_index('time')
    S['ID']= S['ID'].replace({'98': '97'}) #98 is according to Rojo a continuation of 97
    S= S.replace({'comma': 'C', 'undefined': 'U'}) # some wrong morphologies
    S['Morphology']= S['Morphology'].replace({'  ': ' ', 'H - C': 'C - H', 'MGR - C': 'C - MGR', 'T - C': 'C - T', 'U - C': 'C - U', 'H - S':'S - H', 'T - S': 'S - T', 'U - S':'S - U', 'S - W': 'W - S', 'S - MGR': 'MGR - S' }, regex= True)

    S= S[~np.isnan(S.lat)] #removes two rows where lat, lon, time is nan
    return S
    


def read_stars_file(fname):
    """Read data into a `pandas.DataFrame` from the standard file."""

    def _date_parser(*x):
        return pd.datetime.strptime(" ".join(x), "%Y %m %d %H %M")

    dtype_tuple = (int,) + 5 * (str,) + 4 * (float,)
    dtypes = {k: v for k, v in enumerate(dtype_tuple)}

    df = pd.read_csv(
        fname,
        dtype=dtypes,
        sep=r"\s+",
        skiprows=0,
        date_parser=_date_parser,
        parse_dates={"time": [1, 2, 3, 4, 5]},
    )
    return df


def read_Gunnar_stars(STARSdir):
    """Read both North and South subsets of STARS."""
    df_n = read_stars_file(STARSdir + "STARS_TRACKS.csv")
    df_s = read_stars_file(STARSdir + "STARS_TRACKS_south.csv")

    df_s.N += df_n.N.values[-1]

    return df_n.append(df_s).reset_index(drop=True)



def merge_ERA5(S, Mediadir= '/media/pst019/1692A00D929FEF8B/', file= "PL/STARS/STARS_ERA5.csv"):
    S_ERA5 = pd.read_csv(Mediadir+file, sep=',', header= 0)
    S= pd.concat([S, S_ERA5], axis= 1)

    return S



def remove_morphs(S, var, count= 50, how='replace', excludelist= None):
    """remove or replace variable var that occur less than "count" times
    how= ['replace', 'remove'] - replace is putting 'other' instead
    more items specified in excludelist can be removed/replaced
    """
    
    varcounts= S.groupby(var).size()
    varcounts= varcounts[varcounts.values > count]
    varlist= list(varcounts.index)
    
    if excludelist is not None:
        varlist= [x for x in varlist if x not in excludelist]
#    
    if how =='replace':
        boollist= ~S[var].isin(varlist)
#        S[boollist][var]='other'
        S.loc[boollist, var] = 'other'
    
    if how== 'remove':
        S= S.loc[S[var].isin(varlist)]
        
    return S


def STARS_individual_systems(S, remove_list= ['-', 'T', 'U', 'H'], system_char='mean'):
    """make a list of the Systems from the list of all time steps
    system_char=['mean', 'max', 'min', median', 'initial'] - expresses how the system variables are calculated"""
    from f_useful import remove_repeater, split_items, remove_from_list, remove_dublicate

    if system_char== 'mean':  S_ind= S.groupby(['ID']).mean()
    elif system_char== 'max': S_ind= S.groupby(['ID']).max()
    elif system_char== 'min': S_ind= S.groupby(['ID']).min()
    elif system_char== 'median': S_ind= S.groupby(['ID']).median()
    elif system_char== 'initial': S_ind= S.groupby(['ID']).first()
    elif system_char== 'final': S_ind= S.groupby(['ID']).last()
    
    S_ind['ID']= S_ind.index
    
    
    """system morphology"""
    ##S = Spiraliform, C = Comma shaped, MGR = Merry-go-round, W = Wave system, U = Undefined, T = Transition between different forms, H = Hybrid, - = Systems don't appear completely on imagery
    
      
    S['Morph_red']= S['Morphology'].replace({'T - ': '', ' - T': '', '  ': ' ', 'U - ': '', ' - U': '', ' - H': '', 'H - ': ''}, regex= True)
    S_ind['PL_Morph_full']= [" ".join(remove_repeater(split_items(remove_repeater(remove_from_list(S[S['ID'] == ID]['Morph_red'], remove_list)), ' - '))) for ID in S_ind['ID']]
    
    #splits transition, remove dublicates, sort by alphabet
    S['Morph_red']= S['Morphology'].replace({'  ': ' '}, regex= True)
    S_ind['Morphology']= [" ".join(sorted(remove_dublicate(remove_from_list(split_items(S[S['ID']== ID]['Morph_red'], ' - '), remove_list)))) for ID in S_ind['ID']]
    
    S_ind.loc[S_ind['Morphology'].str.contains('MGR'), 'Morphology']= 'MGR+'
    S_ind.loc[S_ind['Morphology'].str.contains('W'), 'Morphology']= 'W+'
    S_ind['Morphology']= S_ind['Morphology'].replace({'C S': 'C-S'})
    
    
    """full morphology"""
    
    #S_ind.loc[S_ind['PL_Morph_full'].str.contains(' MGR '), 'PL_Morph_full']= '..MGR..'
    S_ind.loc[S_ind['PL_Morph_full'].str.contains(' MGR'), 'PL_Morph_full']= '.. MGR'
    S_ind.loc[S_ind['PL_Morph_full'].str.contains('MGR '), 'PL_Morph_full']= 'MGR ..'
    
    #S_ind.loc[np.logical_and.reduce((S_ind['PL_Morph_full'] != 'C W C' , S_ind['PL_Morph_full'] != 'C W S' , S_ind['PL_Morph_full'].str.contains(' W '))), 'PL_Morph_full']= '..W..'
    #S_ind.loc[np.logical_and(~S_ind['PL_Morph_full'].str.contains(' W ')  , S_ind['PL_Morph_full'].str.contains(' W')), 'PL_Morph_full']= '..W'
    
    #S_ind.loc[S_ind['PL_Morph_full'].str.contains(' W '), 'PL_Morph_full']= '..W..'
    #S_ind.loc[S_ind['PL_Morph_full'].str.contains('W '), 'PL_Morph_full']= 'W..'
    #S_ind.loc[S_ind['PL_Morph_full'].str.contains(' W'), 'PL_Morph_full']= '..W'
    
    #S_ind.loc[S_ind['PL_Morph_full'].str.contains('C S C'), 'PL_Morph_full']= 'CSC+'
    #S_ind.loc[np.logical_or(S_ind['PL_Morph_full'] == 'S C S' , S_ind['PL_Morph_full'].str.contains('S C S C')), 'PL_Morph_full']= '. S C S ..'
    #S_ind.loc[S_ind['PL_Morph_full'].str.contains('S C S'), 'PL_Morph_full']= 'CSC+'

#    S_ind['Morphology']= S_ind['PL_Morph']

    
    return S_ind




def Stoll_Obs_nr(Stoll):
    """create the Obs nr of each Stoll PL"""
    
    Stoll['Obs']= ""

    for Stoll_nr in remove_dublicate(Stoll['Stoll nr']):
        Stolltime= Stoll.loc[Stoll['Stoll nr'] == Stoll_nr, 'time']
        timelist= remove_dublicate(np.sort(Stolltime))
        grid= np.meshgrid(timelist, Stolltime)
        
        labellist= np.argwhere(grid[0] == grid[1])
        Stoll.loc[Stoll['Stoll nr'] == Stoll_nr, 'Obs']= labellist[:, 1] +1
    return Stoll



def Stoll_Vort_tend(Stoll):
    """Stoll vorticity tendency
    could be removed since param_tend is more general"""
    import scipy.signal

    Stoll['Vort_tend']= ""
    
    for ID in remove_dublicate(Stoll['ID']):
        vo= Stoll.loc[Stoll['ID'] == ID, 'Vorticity$_{850}$ (centre)'].values
        if len(vo)>= 11: window= 11
        else: window= len(vo) -(len(vo)+1)%2 #to get the largest uneven number smaller than the length of vorticity
        
        vodf= scipy.signal.savgol_filter(vo, window_length=window, polyorder=2, deriv=1)

        Stoll.loc[Stoll['ID'] == ID, 'Vort_tend']= vodf
        
    return Stoll


def Stoll_param_tend(Stoll, parameter='Vorticity$_{850}$ (centre)', outname= 'Vort_tend'):
    """Stoll calculate the tendency in the parameter (first derivative with Obs - like time derivative - if all no gaps would be present
    this is just the tendency from one to the next timestep - inaccurate if gaps occur
    """
    import scipy.signal

    Stoll[outname]= ""
    
    for ID in remove_dublicate(Stoll['ID']):
        vo= Stoll.loc[Stoll['ID'] == ID, parameter].values
        if len(vo)>= 11: window= 11
        else: window= len(vo) -(len(vo)+1)%2 #to get the largest uneven number smaller than the length of vorticity
        
        vodf= scipy.signal.savgol_filter(vo, window_length=window, polyorder=2, deriv=1)

        Stoll.loc[Stoll['ID'] == ID, outname]= vodf
        
    return Stoll


def Stoll_individual_systems(S, remove_list= ['-', 'T', 'U', 'H'], Sdist= 100,
                             ID_name= 'Stoll nr', system_char='mean',
                             lat_name= 'STARS lat', lon_name= 'STARS lon',
                             update_S= False):
    
    """make a list of the Systems from the list of all time steps
    Sdist is the maximum distance the Rojo track and the matched Stoll track has to have in order to included for the cloud morphology
    ID_name is the name of the variable that labels the Stoll systems 
    system_char=['mean', 'max', 'min', median', 'initial'] - expresses how the system variables are calculated
    S[Morphology] is the original
    S[Morph_red] - gives a reduction of the Morphology by sorting transition situations to clear situations
    update_S - specifies if the cloud morphology of S should also be updated
    """


    if system_char== 'mean':  S_ind= S.groupby([ID_name]).mean()
    elif system_char== 'max': S_ind= S.groupby([ID_name]).max()
    elif system_char== 'min': S_ind= S.groupby([ID_name]).min()
    elif system_char== 'median': S_ind= S.groupby([ID_name]).median()
    elif system_char== 'initial': S_ind= S.groupby([ID_name]).first()
    elif system_char== 'final': S_ind= S.groupby([ID_name]).last()
    
    S_ind[ID_name]= S_ind.index

    #trivial reduction
    S['Morph_red']= S['Morphology'].replace({'T - ': '', ' - T': '', '  ': ' ', 'U - ': '', ' - U': '', ' - H': '', 'H - ': ''}, regex= True)
    #stronger reduction
    S['Morph_red']= S['Morph_red'].replace({'C - ': '', ' - S': '', '-': 'other', 'H': 'other', 'U': 'other', 'T': 'other'}, regex= True)

    S['dist']= distance((S.lat, S.lon), (S[lat_name], S[lon_name])) 
    
    S_ind['PL_Morph_full']= [" ".join(
            remove_repeater(split_items(remove_repeater(remove_from_list(
                    S[np.logical_and(S[ID_name] == ID, S['dist'] <= Sdist)]['Morph_red'], remove_list)), ' - '))) 
        for ID in S_ind.index]
    
#    S['Morph_red']= S['Morphology'].replace({'  ': ' '}, regex= True)
    S_ind['Morphology']= [" ".join(
            sorted(remove_dublicate(remove_from_list(split_items(
                    S[np.logical_and(S[ID_name] == ID, S['dist'] <= Sdist)]['Morph_red'], ' - '), remove_list))))
         for ID in S_ind.index]
    
    S_ind.loc[S_ind['Morphology'].str.contains('MGR'), 'Morphology']= 'MGR+'
    S_ind.loc[S_ind['Morphology'].str.contains('W'), 'Morphology']= 'W+'
    S_ind['Morphology']= S_ind['Morphology'].replace({'C S': 'C-S'})
    
    """full morphology"""
    #S_ind.loc[S_ind['PL_Morph_full'].str.contains(' MGR '), 'PL_Morph_full']= '..MGR..'
    S_ind.loc[S_ind['PL_Morph_full'].str.contains(' MGR'), 'PL_Morph_full']= '.. MGR'
    S_ind.loc[S_ind['PL_Morph_full'].str.contains('MGR '), 'PL_Morph_full']= 'MGR ..'
    
    #S_ind.loc[np.logical_and.reduce((S_ind['PL_Morph_full'] != 'C W C' , S_ind['PL_Morph_full'] != 'C W S' , S_ind['PL_Morph_full'].str.contains(' W '))), 'PL_Morph_full']= '..W..'
    #S_ind.loc[np.logical_and(~S_ind['PL_Morph_full'].str.contains(' W ')  , S_ind['PL_Morph_full'].str.contains(' W')), 'PL_Morph_full']= '..W'
    
    #S_ind.loc[S_ind['PL_Morph_full'].str.contains(' W '), 'PL_Morph_full']= '..W..'
    #S_ind.loc[S_ind['PL_Morph_full'].str.contains('W '), 'PL_Morph_full']= 'W..'
    #S_ind.loc[S_ind['PL_Morph_full'].str.contains(' W'), 'PL_Morph_full']= '..W'
    
    #S_ind.loc[S_ind['PL_Morph_full'].str.contains('C S C'), 'PL_Morph_full']= 'CSC+'
    #S_ind.loc[np.logical_or(S_ind['PL_Morph_full'] == 'S C S' , S_ind['PL_Morph_full'].str.contains('S C S C')), 'PL_Morph_full']= '. S C S ..'
    #S_ind.loc[S_ind['PL_Morph_full'].str.contains('S C S'), 'PL_Morph_full']= 'CSC+'

#    S_ind['Morphology']= S_ind['PL_Morph']
    
    
    if update_S:
        return S, S_ind
    else:
        return S_ind


def Stoll_excl(Stoll, double=True, lifetime_threshold= 2, landfraction= False, land=True, give_print=True):
    """excludes tracks from Stoll list
    - if double= True, exclude systems that occur multiple times
    - it excludes systems that have a lifetime shorter than the lifetime threshold
    - landfraction can take a value between 0 and 1, then this gives the fraction of time steps above which the system is excluded
    - if land= True, exclude systems that are on land at start, middle and end
    - lifetime_threshold is in hours (so number of timesteps +1)"""
    
    excl_list=[]
    if double == True:
        """ exclude the dublicate - sometimes the Stoll systems matched to different Rojo PLs are the same """
        excl_list_double=[] #this is the list of the Stoll matches that are double in the Stoll list
        excl_dict_double={} #the Stoll nr of excluded system and the Rojo number that was assigned to it
        
        for i_track_file in np.arange(np.min(Stoll['track file'].values), np.max(Stoll['track file'].values)+1):
            Stoll_trnr= Stoll[Stoll['track file'] == i_track_file]    
        
            Stoll_nrs_in_trnr= remove_dublicate(Stoll_trnr['Stoll nr'])
            
            d= {} #dictionary in which the Stoll nr and the corresponding tr_id_nr are saved
            for Stoll_nr in Stoll_nrs_in_trnr:
                d[Stoll_nr]= remove_dublicate(Stoll_trnr[Stoll_trnr['Stoll nr'] == Stoll_nr]['track_idx'])
        
            for i_S, Stoll_nr in enumerate(list(d.keys())):
                #if two Stoll_nr have a common member in the tr_id_nr it is added to the excl_list
                if common_member(list(d.values())[i_S], flat_list(list(d.values())[:i_S])):
                    excl_list_double += [Stoll_nr]
                    excl_dict_double[Stoll_nr]= remove_dublicate(Stoll[Stoll['Stoll nr'] == Stoll_nr]['Rojo nr'])[0]

            #for some reason the order matters?!
#            for i_S, Stoll_nr in enumerate(list(d.keys())):
#                #if two Stoll_nr have a common member in the tr_id_nr it is added to the excl_list
#                if common_member(list(d.values())[i_S], flat_list(list(d.values())[(i_S+1):])):
#                    excl_list_double += [Stoll_nr]
#                    excl_dict_double[Stoll_nr]= remove_dublicate(Stoll[Stoll['Stoll nr'] == Stoll_nr]['Rojo nr'])[0]

         
        if give_print: print('Nr PLs excluded due to double occurrance: ', len(excl_list_double))
        excl_list += excl_dict_double

           
    """exclude systems with short lifetime""" 
    Stoll_tracks= prepare_tracks(Stoll, system_id_label="Stoll nr")


    excl_list_lifetime= [track['Stoll nr'][0] for track in Stoll_tracks if track.lifetime_h < lifetime_threshold]
    
    if give_print: print('Nr PLs excluded by lifetime: ', len(excl_list_lifetime))
    excl_list += excl_list_lifetime
    
    """exclude systems occuring mainly on land"""
    if landfraction != False:
        print('have to write this')
        lsm = xr.open_dataarray(Mediadir+ "ERA5_STARS/data/lsm.nc")
        lsm = lsm.squeeze()  # remove singular time dimension
        
        lsm.values[lsm > 0]= 1

        excl_list_land=[]
        for track in Stoll_tracks:
            land_list= np.array([int(lsm.sel(latitude= track.lat.values[i], longitude=track.lon.values[i])) for i in range(len(track.lat))])
            

            if sum(land_list)/len(land_list) > landfraction: excl_list_land += [track['Stoll nr'][0]]
    
        if give_print: print('Nr PLs excluded by land: ', len(excl_list_land))
        excl_list += excl_list_land    
    
    elif land == True:
        lsm = xr.open_dataarray(Mediadir+ "ERA5_STARS/data/lsm.nc")
        lsm = lsm.squeeze()  # remove singular time dimension
        
        lsm.values[lsm > 0]= 1
        #extremely unefficient
        #Stoll['land']= [int(lsm.sel(latitude= Stoll.lat.values[i], longitude=Stoll.lon.values[i]) ) for i in range(len(Stoll))]
        
        excl_list_land=[]
        for track in Stoll_tracks:
            land_start=int(lsm.sel(latitude= track.lat.values[0], longitude=track.lon.values[0]) )
            land_final=int(lsm.sel(latitude= track.lat.values[-1], longitude=track.lon.values[-1]) )
            i_middle= int(track.lifetime_h//2)
            land_middle=int(lsm.sel(latitude= track.lat.values[i_middle], longitude=track.lon.values[i_middle]) )
        
            if land_start+ land_middle+ land_final== 3: excl_list_land += [track['Stoll nr'][0]]
    
        if give_print: print('Nr PLs excluded by land: ', len(excl_list_land))
        excl_list += excl_list_land


    Stoll= Stoll[[Stoll_nr not in excl_list for Stoll_nr in Stoll['Stoll nr']]]
    Stoll.index = np.arange(len(Stoll))
    return Stoll



def Stoll_excl_lifetime(Stoll, lifetime_threshold= 5):
    """exclude systems with short lifetime
    lifetime_threshold of 5 is the same as 6 hourly time steps""" 
    Stoll_tracks= prepare_tracks(Stoll, system_id_label="Stoll nr")


    excl_list_lifetime= [track['Stoll nr'][0] for track in Stoll_tracks if track.lifetime_h < lifetime_threshold]
    
    print('Nr PLs excluded by lifetime: ', len(excl_list_lifetime))

    Stoll= Stoll[[Stoll_nr not in excl_list_lifetime for Stoll_nr in Stoll['Stoll nr']]]
    Stoll.index = np.arange(len(Stoll))
    return Stoll


def imp_standard_Stoll(drop=True, rename=True,
                       ERA5_file= "PL/STARS/Stoll_ERA5_handfix.csv"): #SOM_filedir):#, Obs_evol):
    """drop specifies if variables in Stoll and Stoll_ind are dropped"""
    """maybe split in two"""
    """Stoll systems"""
    test= 'version4'
    dist= 150
    system_char='initial' # for the Stoll_ind
    
    Stoll_STARS_file="ERA5_STARS/STARS-PMC-merge/"+test+"_dist_"+str(dist)+"_handfix.csv"
    Stoll = pd.read_csv(Mediadir+ Stoll_STARS_file)
    Stoll['time']= pd.DatetimeIndex(Stoll.time)
    
    Stoll= Stoll_Obs_nr(Stoll) #creates Stoll['Obs']
    Stoll= Stoll_excl(Stoll)
    Stoll= Stoll.rename(columns={'Stoll nr': 'ID'})
    
    Stoll = merge_ERA5(Stoll, Mediadir, ERA5_file)
    if rename: Stoll= Stoll.rename(columns={
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
    
    if drop: Stoll= Stoll.drop(columns=['Comment', 'track file', 'Rojo nr', 'Rojo nr old', 'Press_min', 'U_10min_knots', 'row_idx', 'track_idx'])
       
    
    """Stoll_individual_systems"""
    Stoll, Stoll_ind= Stoll_individual_systems(Stoll, ID_name='ID', system_char=system_char, update_S= True)
    Stoll_ind= calc_system_duration(Stoll, Stoll_ind, ID_name= 'ID', Obs_name='Obs')
    
    if drop:
        Stoll_ind= Stoll_ind.drop(columns=[ 'STARS lat', 'STARS lon', 'STARS Obs nr'])
        Stoll= Stoll.drop(columns=[ 'STARS lat', 'STARS lon', 'STARS Obs nr', 'dist'])
    
    if rename:
        Stoll_ind= calc_Obs_mature(Stoll, Stoll_ind, intensity_var='Vorticity$_{850}$ (centre)')

#    df_nodes= pd.read_csv(SOM_filedir+"_node_nr.txt", sep=" ")
#    df_dist= pd.read_csv(filedir+"_distances.txt", sep=" ")
#    df_nodes= pd.concat([df_nodes, df_dist], axis= 1)


##    if Obs_evol!='allObs':
#    df_nodes['ID']= [str(int(df_nodes.PLnr[n]))[:-3]+'_'+str(int(df_nodes.PLnr[n])%1000) for n in range(len(df_nodes.PLnr))]
#    df_nodes= df_nodes.set_index('ID')
##        df_nodes= df_nodes.drop(columns='as.character(dates)')
#
#    Stoll_ind=Stoll_ind.join(df_nodes['node'], how='outer')
##        Stoll_ind=Stoll_ind.join(df_nodes['K_SOM.distances'], how='outer')
#
#    
##    if Obs_evol=='allObs':
#    df_nodes= pd.read_csv(SOM_filedir+"_node_nr.txt", sep=" ")
#
#    df_nodes['time']= pd.to_datetime(df_nodes.date.apply(str), format='%Y%m%d%H')
#    df_nodes['ID']= [str(int(df_nodes.PLnr[n]))[:-3]+'_'+str(int(df_nodes.PLnr[n])%1000) for n in range(len(df_nodes.PLnr))]
#
##        df_nodes= df_nodes.drop(columns=['date', 'PLnr'])
#    Stoll= Stoll.set_index(['ID', 'time'])
#    df_nodes= df_nodes.set_index(['ID', 'time'])
#    
#    Stoll= pd.concat([Stoll, df_nodes], axis= 1)
    
    return Stoll, Stoll_ind



def calc_Obs_mature(S, S_ind, intensity_var='vo'):
    """calculate the observation number of the mature phase, defined as maximum in the intensity variable """
    
    if 'ID' not in S.columns:
        S= S.reset_index()

    
    Obs_mature= []

    for ID_now in S_ind.index.values:
        S_now= S.loc[S['ID'] == ID_now]
        Obs_mature+= [S_now.iloc[np.argmax(S_now[intensity_var].values)].Obs]
    
#    S_ind.reindex([remove_dublicate(S.ID)])
    S_ind['Obs_mature']= Obs_mature
    return S_ind


def calc_Obs_mature_2(S, S_ind, intensity_var='vo'):
    """calculate the observation number of the mature phase, defined as maximum in the intensity variable
    This version also adds the mature stage with a 1 to the S dataset (0 elsewhere)"""

    if 'ID' not in S.columns:
        reset_index= True
        S= S.reset_index()
    else: reset_index= False

    Obs_mature_list= []
    S['Obs_mature']= 0
    
    for ID_now in S_ind.index.values:
#        print(ID_now)
        S_now= S.loc[S['ID'] == ID_now]
        obs_mature= S_now.iloc[np.argmax(S_now[intensity_var].values)].Obs
        
#        S[np.logical_and(S['ID']== ID_now, S['Obs'] == obs_mature)]['Obs_mature']= 1
        Obs_mature_list += [obs_mature]
        
#    S_ind.reindex([remove_dublicate(S.ID)])
    S_ind['Obs_mature']= Obs_mature_list
    
    for i in range(len(S_ind)):
        S.loc[np.logical_and(S['ID']== S_ind['Obs_mature'].index[i], S['Obs']== S_ind['Obs_mature'][i]), 'Obs_mature']= 1
    
    if reset_index:
        S= S.set_index(['ID', 'time'])
    
    return S, S_ind


def calc_system_duration(S, S_ind, ID_name= 'ID', Obs_name='Obs'):
    """calculates the duration of the STARS/Stoll PLs from the timestel list (STARS) and the individual system list (STARS_ind)"""
    starttime= S.loc[S[Obs_name] == 1].groupby([ID_name]).first()['time']


    endobs= S.groupby([ID_name]).last()[Obs_name]
    
    endtime= S[np.logical_and(S[ID_name] == endobs.index[0], S[Obs_name] == endobs[0])][[ID_name, 'time']]
    for i in range(1, len(endobs)):
        a= S[np.logical_and(S[ID_name] == endobs.index[i], S[Obs_name] == endobs[i])][[ID_name, 'time']]
        if len(a) > 1: print(a)
        endtime= pd.concat([endtime, a])
    
    endtime= endtime.set_index(ID_name)['time']
    
    duration= endtime- starttime
    
    S_ind['Duration']= duration/np.timedelta64(1,'h')
    return S_ind




def imp_ds(var, plevel_var, var_full, lsm_mask=False, compute=False, imp_var=None, data_gausfilter='', smooth_param= '1E-3'):
    """import one variable of the polar low centred grid"""
    
    if compute== False:
        if plevel_var != '': plev_str= '_all_levs'
        else: plev_str= ''
        
        file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +var + plev_str + '_allObs_track-smth-'+smooth_param+'.nc'
        ds= xr.open_dataset(file)
        if plevel_var != '':
            ds= ds.sel(plev= plevel_var)
            ds = ds.drop('plev')
        ds= ds.rename({var: var_full})



    elif compute:
        for ni, i_var in enumerate(imp_var):
            file= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +i_var +'_allObs_track-smth-'+smooth_param+'.nc'
            if os.path.isfile(file):
                ds1= xr.open_dataset(file)
            else:
                ds1= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +i_var + '_all_levs_allObs_track-smth-'+smooth_param+'.nc')
                ds1= ds1.sel(plev= plevel_var)
                
            if ni == 0: ds2= ds1
            else: ds2= xr.merge([ds2, ds1]) 


        if var== 'adv_t':
            if data_gausfilter:
                ds2['t'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.t, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
                ds2['t'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.t, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )

            grad_t_vec= np.gradient(ds2.t, 25E3) #two components of the temperature gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km

            U= np.sqrt(ds2.u**2 + ds2.v**2)
            wind_beering= UV2Direction(ds2.u, ds2.v)
            track_beering = ds2.beering.values
            
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='')
            
            ds2[var_full]= (('time', 'x', 'y'), -1* (u_r * grad_t_vec[2] + v_r * grad_t_vec[1] ) * 3600 )
            ds2[var_full].attrs['units'] = 'K/h'
            ds2[var_full].attrs['long_name']='Horizontal temperature advection' 

  
        if var== 'tp':                  
            ds2[var_full]= (('time', 'x', 'y'), (ds2['cp']+ ds2['lsp'] ))
#            ds2[var_full] *= 24
            ds2[var_full].attrs['units'] = 'mm $\cdot$ h^{-1}'
            ds2[var_full].attrs['long_name']='Total precipitation' 



        if var=='t_diff':
            ds2[var_full]= (('time', 'x', 'y'), ds2['t'].isel(plev=0) - ds2['t'].isel(plev=1) )
            ds2[var_full].attrs['units'] = 'K'

        if var=='sst-t':
            ds2[var_full]= (('time', 'x', 'y'), ds2['sst'] - ds2['t'])
            ds2[var_full].attrs['units'] = 'K'            

        if var=='sst-theta':
            ds2[var_full]= (('time', 'x', 'y'), ds2['sst'] - PotTemp(ds2.t, plevel_var) )
            ds2[var_full].attrs['units'] = 'K'     
            ds2[var_full].attrs['long_name']='SST - theta'+str(plevel_var)

        if var=='theta_sst-theta':
            theta= PotTemp(ds2.t, plevel_var)
            theta_sst= PotTemp(ds2.sst, ds2.msl/100)

            ds2[var_full]= (('time', 'x', 'y'), theta_sst - theta )
            ds2[var_full].attrs['units'] = 'K'     
            ds2[var_full].attrs['long_name']='theta_sst-theta'+str(plevel_var)


        if var=='theta_e_sst-theta_e':
            theta_e= EquiPotTemp(ds2.t, ds2.q/1E3, plevel_var)
#            theta_e_sst= EquiPotTemp(ds2.sst, RH2SH(100, 1000, ds2.sst)/1E3, 1000)
            theta_e_sst= EquiPotTemp(ds2.sst, RH2SH(100, ds2.msl/100, ds2.sst)/1E3, ds2.msl/100)

            ds2[var_full]= (('time', 'x', 'y'), theta_e_sst - theta_e )
            ds2[var_full].attrs['units'] = 'K'     
            ds2[var_full].attrs['long_name']='theta_e_sst-theta_e'+str(plevel_var)


        if var=='theta_e_2m-theta_e':
            theta_e= EquiPotTemp(ds2.t, ds2.q/1E3, plevel_var)
#            theta_e_sst= EquiPotTemp(ds.sst, RH2SH(100, 1000, ds.sst)/1E3, 1000)
            theta_e_2m= EquiPotTemp(ds2['2t'], RH2SH(100, ds2.msl/100, ds2['2t'])/1E3, ds2.msl/100)

            ds2[var_full]= (('time', 'x', 'y'), theta_e_2m - theta_e )
            ds2[var_full].attrs['units'] = 'K'     
            ds2[var_full].attrs['long_name']='theta_e_2m-theta_e'+str(plevel_var)


        if var=='theta_diff':
            theta_1, theta_0= PotTemp(ds2.t.isel(plev=1), plevel_var[1]), PotTemp(ds2.t.isel(plev=0), plevel_var[0])

            ds2[var_full]= (('time', 'x', 'y'), theta_0 - theta_1)
            ds2[var_full].attrs['units'] = 'K'

        if var=='theta_e_diff':
            theta_e_1, theta_e_0= EquiPotTemp(ds2.t.isel(plev=1), ds2.q.isel(plev=1)/1E3, plevel_var[1]), EquiPotTemp(ds2.t.isel(plev=0), ds2.q.isel(plev=0)/1E3, plevel_var[0])
            
            print('here')
            ds2[var_full]= (('time', 'x', 'y'), theta_e_0 - theta_e_1)
            ds2[var_full].attrs['units'] = 'K'

        if var=='rh':
            rh= SH2RH(ds2.q, plevel_var, ds2.t)

            ds2[var_full]= (('time', 'x', 'y'), rh)
            ds2[var_full].attrs['units'] = '%'
            ds2[var_full].attrs['long_name']='Relative humidity'
            
            
        if var=='N':
            h_diff=  (ds2['z'].isel(plev=1) - ds2['z'].isel(plev=0) ) ###later: why not /9.81?

            theta_1, theta_0= PotTemp(ds2.t.isel(plev=1), plevel_var[1]), PotTemp(ds2.t.isel(plev=0), plevel_var[0])
            N= np.sqrt(9.81/(theta_1+theta_0)/2 *  (theta_1- theta_0)/h_diff)        

            ds2[var_full]= (('time', 'x', 'y'),  N)
            ds2[var_full].attrs['units']= '1/s' 

        if var== 'grad_t':
            if data_gausfilter:
                ds2['t'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.t, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
                ds2['t'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.t, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )

            grad_t_vec= np.gradient(ds2.t, 0.25) #two components of the temperature gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km
            ds2[var_full]= (('time', 'x', 'y'), np.sqrt(grad_t_vec[0]**2 + grad_t_vec[1]**2) )
            ds2[var_full].attrs['units']= 'K/100km'

        if var== 'u_g':
#            if data_gausfilter:
#                ds2['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.z, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
#                ds2['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.z, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )
#
            f= np.sin(np.deg2rad(Stoll.lat.values)) * 2*2*np.pi/(24*60**2)

            grad_z_vec = np.gradient(ds2.z *9.81 , 25E3)
            u_g= -1/f[:, np.newaxis, np.newaxis]* grad_z_vec[1]
            
            ds2[var_full]= (('time', 'x', 'y'), u_g )
            ds2[var_full].attrs['units'] = 'm/s'
            ds2[var_full].attrs['long_name']='Tangential geostrophic wind'
            
        if var== 'v_g':
#            if data_gausfilter:
#                ds2['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.z, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
#                ds2['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.z, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )

            f= np.sin(np.deg2rad(Stoll.lat.values)) * 2*2*np.pi/(24*60**2)
            grad_z_vec = np.gradient(ds2.z *9.81 , 25E3)
            v_g= 1/f[:, np.newaxis, np.newaxis]* grad_z_vec[2] 
            
            ds2[var_full]= (('time', 'x', 'y'), v_g )
            ds2[var_full].attrs['units'] = 'm/s'
            ds2[var_full].attrs['long_name']='Azimuthal geostrophic wind'            

        if var== 'vo_g':
#            if data_gausfilter:
#                ds2['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.z, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
#                ds2['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.z, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )

            f= np.sin(np.deg2rad(Stoll.lat.values)) * 2*2*np.pi/(24*60**2)
            grad_z_vec = np.gradient(ds2.z *9.81 , 25E3)
            u_g= -1/f[:, np.newaxis, np.newaxis]* grad_z_vec[1]
            v_g= 1/f[:, np.newaxis, np.newaxis]* grad_z_vec[2] 
            
            vo_g = np.gradient(v_g, 25E3, axis= 2) - np.gradient(u_g, 25E3, axis= 1)
            
            ds2[var_full]= (('time', 'x', 'y'), vo_g*1E5 )
            ds2[var_full].attrs['units'] = '1/s'
            ds2[var_full].attrs['long_name']='Geostrophic vorticity'  




        if var== 'wind_r': #"""rotate the wind with propagation direction"""
            uvar, vvar= imp_var[0], imp_var[1]
            uvar_full, vvar_full= var_full[0], var_full[1]
            
            U= np.sqrt(ds2[uvar]**2 + ds2[vvar]**2)
            wind_beering= UV2Direction(ds2[uvar], ds2[vvar])
            track_beering = ds2.beering.values
            
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='') #maybe u_zonal and v_meridional are exchanged in WindSpeedDirection2UV, but this seem to work"""
        
            ds2[uvar_full]= (('time', 'x', 'y'), u_r )
            ds2[vvar_full]= (('time', 'x', 'y'), v_r )
        
            ds2[uvar_full].attrs['units'] = 'm/s'
            ds2[uvar_full].attrs['long_name']='x wind'
        
            ds2[vvar_full].attrs['units'] = 'm/s'
            ds2[vvar_full].attrs['long_name']='y wind'

        if var== 'u_r':
            U= np.sqrt(ds2.u**2 + ds2.v**2)
            wind_beering= UV2Direction(ds2.u, ds2.v)
            track_beering = ds2.beering.values
            
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindspeedDirection2UV(U, rot_beering, orientation='') #maybe u_zonal and v_meridional are exchanged in Winds2peedDirection2UV, but this seem to work"""

            ds2[var_full]= (('time', 'x', 'y'), u_r )
            ds2[var_full].attrs['units'] = 'm/s'
            ds2[var_full].attrs['long_name']='Tangential wind'

        if var== 'v_r':
            U= np.sqrt(ds2.u**2 + ds2.v**2)
            wind_beering= UV2Direction(ds2.u, ds2.v)
            track_beering = ds2.beering.values
            
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindspeedDirection2UV(U, rot_beering, orientation='')
            ds2[var_full]= (('time', 'x', 'y'), v_r )
            ds2[var_full].attrs['units'] = 'm/s'
            ds2[var_full].attrs['long_name']='Azimuthal wind'

        if var== 'u_a':
            U= np.sqrt(ds2.u**2 + ds2.v**2)
            wind_beering= UV2Direction(ds2.u, ds2.v)
            track_beering = ds2.beering.values           
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindspeedDirection2UV(U, rot_beering, orientation='') #maybe u_zonal and v_meridional are exchanged in Winds2peedDirection2UV, but this seem to work"""

            f= np.sin(np.deg2rad(Stoll.lat.values)) * 2*2*np.pi/(24*60**2)
            grad_z_vec = np.gradient(ds2.z *9.81 , 25E3)
            u_g= -1/f[:, np.newaxis, np.newaxis]* grad_z_vec[1]
#            v_g= 1/f[:, np.newaxis, np.newaxis]* grad_z_vec[2] 
            
            u_a= u_r - u_g
#            v_a= v_r - v_g
            
            ds2[var_full]= (('time', 'x', 'y'), u_a )
            ds2[var_full].attrs['units'] = 'm/s'
            ds2[var_full].attrs['long_name']='Ageostrophic tangential wind' 

        if var== 'v_a':
            U= np.sqrt(ds2.u**2 + ds2.v**2)
            wind_beering= UV2Direction(ds2.u, ds2.v)
            track_beering = ds2.beering.values           
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='') #maybe u_zonal and v_meridional are exchanged in WindSpeedDirection2UV, but this seem to work"""

            f= np.sin(np.deg2rad(Stoll.lat.values)) * 2*2*np.pi/(24*60**2)
            grad_z_vec = np.gradient(ds2.z *9.81 , 25E3)
            v_g= 1/f[:, np.newaxis, np.newaxis]* grad_z_vec[2] 
            
            v_a= v_r - v_g
            
            ds2[var_full]= (('time', 'x', 'y'), v_a )
            ds2[var_full].attrs['units'] = 'm/s'
            ds2[var_full].attrs['long_name']='Ageostrophic azimuthal wind' 
            
        if var== 'U':
            U= np.sqrt(ds2.u**2 + ds2.v**2)

            ds2[var_full]= (('time', 'x', 'y'), U )
            ds2[var_full].attrs['units'] = 'm/s'
            ds2[var_full].attrs['long_name']= 'Wind speed'

        if var== '10U':
            U= np.sqrt(ds2['10u']**2 + ds2['10v']**2)

            ds2[var_full]= (('time', 'x', 'y'), U )
            ds2[var_full].attrs['units'] = r"m $\cdot$ s$^{-1}$"
            ds2[var_full].attrs['long_name']= 'Wind speed at 10m'

#        if var== 'U_r':
#            U= np.sqrt(ds2.u**2 + ds2.v**2)
#            wind_beering= UV2Direction(ds2.u, ds2.v)
#            track_beering = ds2.beering.values
#            
#            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
#            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='') #maybe u_zonal and v_meridional are exchanged in WindSpeedDirection2UV, but this seem to work"""
#
#            ds2[var_full]= (('time', 'x', 'y'), np.sqrt(u_r**2 + v_r**2) )
#            ds2[var_full].attrs['units'] = 'm/s'
#            ds2[var_full].attrs['long_name']='Rotated wind speed'

        if var== 'U_g':
            f= np.sin(np.deg2rad(Stoll.lat.values)) * 2*2*np.pi/(24*60**2)
            grad_z_vec = np.gradient(ds2.z *9.81 , 25E3)
            u_g= -1/f[:, np.newaxis, np.newaxis]* grad_z_vec[1]
            v_g= 1/f[:, np.newaxis, np.newaxis]* grad_z_vec[2] 
            
            ds2[var_full]= (('time', 'x', 'y'), np.sqrt(u_g**2 + v_g**2) )
            ds2[var_full].attrs['units'] = 'm/s'
            ds2[var_full].attrs['long_name']='Geostrophic wind speed' 

        if var== 'U_a':
            U= np.sqrt(ds2.u**2 + ds2.v**2)
            wind_beering= UV2Direction(ds2.u, ds2.v)
            track_beering = ds2.beering.values           
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='') #maybe u_zonal and v_meridional are exchanged in WindSpeedDirection2UV, but this seem to work"""

            f= np.sin(np.deg2rad(Stoll.lat.values)) * 2*2*np.pi/(24*60**2)
            grad_z_vec = np.gradient(ds2.z *9.81 , 25E3)
            u_g= -1/f[:, np.newaxis, np.newaxis]* grad_z_vec[1]
            v_g= 1/f[:, np.newaxis, np.newaxis]* grad_z_vec[2] 
            
            u_a= u_r - u_g
            v_a= v_r - v_g
            
            ds2[var_full]= (('time', 'x', 'y'), np.sqrt(u_a**2 + v_a**2) )
            ds2[var_full].attrs['units'] = 'm/s'
            ds2[var_full].attrs['long_name']='Ageostrophic wind speed'  

           
        if var== 'adv_t_g':
            if data_gausfilter:
                ds2['t'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.t, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
                ds2['t'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.t, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )
#                ds2['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.z, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
#                ds2['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.z, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )

            grad_t_vec= np.gradient(ds2.t, 25E3) #two components of the temperature gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km
            f= np.sin(np.deg2rad(Stoll.lat.values)) * 2*2*np.pi/(24*60**2)
            grad_z_vec = np.gradient(ds2.z *9.81 , 25E3)
            u_g= -1/f[:, np.newaxis, np.newaxis]* grad_z_vec[1]
            v_g= 1/f[:, np.newaxis, np.newaxis]* grad_z_vec[2] 
            
            ds2[var_full]= (('time', 'x', 'y'), -1* (u_g * grad_t_vec[2] + v_g * grad_t_vec[1] ) * 3600 )
            ds2[var_full].attrs['units'] = 'K/h'
            ds2[var_full].attrs['long_name']='Horizontal geostrophic temperature advection' 


        if var== 'adv_t':
            if data_gausfilter:
                ds2['t'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.t, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
                ds2['t'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.t, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )

            grad_t_vec= np.gradient(ds2.t, 25E3) #two components of the temperature gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km

            U= np.sqrt(ds2.u**2 + ds2.v**2)
            wind_beering= UV2Direction(ds2.u, ds2.v)
            track_beering = ds2.beering.values
            
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='')
            
            ds2[var_full]= (('time', 'x', 'y'), -1* (u_r * grad_t_vec[2] + v_r * grad_t_vec[1] ) * 3600 )
            ds2[var_full].attrs['units'] = 'K/h'
            ds2[var_full].attrs['long_name']='Horizontal temperature advection' 


      
        if var== 'adv_q':
            if data_gausfilter:
                ds2['q'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.q, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
                ds2['q'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.q, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )
#                ds2['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.z, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
#                ds2['z'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.z, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )

            grad_q_vec= np.gradient(ds2.q, 25E3) #two components of the temperature gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km
            
            U= np.sqrt(ds2.u**2 + ds2.v**2)
            wind_beering= UV2Direction(ds2.u, ds2.v)
            track_beering = ds2.beering.values
            
            rot_beering= (track_beering[:, np.newaxis, np.newaxis]- wind_beering)%360
            v_r, u_r = WindSpeedDirection2UV(U, rot_beering, orientation='')
            
            ds2[var_full]= (('time', 'x', 'y'), -1* (u_r * grad_q_vec[2] + v_r * grad_q_vec[1] ) * 3600 )
            ds2[var_full].attrs['units'] = '(g/kg)/h'
            ds2[var_full].attrs['long_name']='Horizontal humidity advection' 



        if var== 'adv_q_g':
            if data_gausfilter:
                ds2['q'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.q, sigma= data_gausfilter, axis= -1, mode='nearest', truncate= 1.) )
                ds2['q'] = (('time', 'x', 'y'), filters.gaussian_filter1d(ds2.q, sigma= data_gausfilter, axis= -2, mode='nearest', truncate= 1.) )

            grad_q_vec= np.gradient(ds2.q, 25E3) #two components of the temperature gradient. Since the grid spacing is 25km, 0.25 is the grid spacing such that the gradient is K/100km
            f= np.sin(np.deg2rad(Stoll.lat.values)) * 2*2*np.pi/(24*60**2)
            grad_z_vec = np.gradient(ds2.z *9.81 , 25E3)
            u_g= -1/f[:, np.newaxis, np.newaxis]* grad_z_vec[1]
            v_g= 1/f[:, np.newaxis, np.newaxis]* grad_z_vec[2] 
            
            ds2[var_full]= (('time', 'x', 'y'), -1* (u_g * grad_q_vec[2] + v_g * grad_q_vec[1] ) * 3600 )
            ds2[var_full].attrs['units'] = '(g/kg)/h'
            ds2[var_full].attrs['long_name']='Horizontal humidity advection' 

        if var== 'flux':
            ds2[var_full]= (('time', 'x', 'y'), ds2['sshf'] + ds2['slhf'] )
            ds2[var_full].attrs['long_name']='Turbulent heat flux'
#            ds2= ds2[var_full]

#        if var== 'Bowen':
#            
#            ds2[var_full]= (('time', 'x', 'y'), ds2['sshf']/ds2['slhf'] )
#            ds2[var_full].attrs['units'] = 'sensible/latent'
#            ds2[var_full].attrs['long_name']='Bowen ratio' 
#            ds2= ds2[var_full]
      

        ds=  ds2 #[var_full]


#    print(ds)

    """do land-sea ice mask""" 
    if lsm_mask:
        ds3= xr.open_dataset( Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/ci_allObs_track-smth-"+smooth_param+'.nc')
        ds3['ci']= ds3['ci'].fillna(value= 1)
           
        ds[var_full]= ds[var_full].where(ds3['ci'] < 0.3)


    if var_full == 'msl':
        ds['msl']/= 1E2
        ds['msl'].attrs['units']= 'hPa'

    if var_full == 'z_pv':
        ds['z_pv']/= 9.81
        ds['z_pv'].attrs['units']= 'm'

    if var == 'd':
        ds[var_full]*= 1E5
        ds[var_full].attrs['units']= '1e-5 '+ ds[var_full].attrs['units']
        
    if var_full in ['slhf', 'sshf', 'flux'] :
        ds[var_full]/= -60**2
        ds[var_full].attrs['units'] = r"W m$^{-2}$"  
        
    if var_full in ['tp', 'cp', 'lsp']:
        ds[var_full]*= 1E3
        ds[var_full].attrs['units']= 'mm'




    return ds

    

def OctantTrack_to_df(track):
    """make OctantTrack to pandas dataframe"""

    tr_nr= track.index.labels[0][0]
    max_row_index= track.index.labels[1][-1]
    
    xr= track.to_xarray()
    xr= xr.sel(track_idx = tr_nr)
    xr= xr.where(xr.row_idx <= max_row_index, drop=True)
    return xr.to_dataframe()