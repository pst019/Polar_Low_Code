#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 10:39:19 2019

@author: pst019
"""

import pandas as pd 
import numpy as np


def import_STARS(Mediadir= '/media/pst019/1692A00D929FEF8B/', file= "PL/STARS/Rojo-etal_2019.csv"):

    S = pd.read_csv(Mediadir+file, sep=',', header= 27)


    new_cols= list(S.columns)
    new_cols[1] = 'Time'
    new_cols[6:8] = ['Diameter', 'Optional_diameter']
    new_cols[10:16] = ['Comment', 'Comment_visual', 'Comment_uncertainties', 'Press_min', 'U_10min_knots', 'U_3sec_kntots']
    S.columns= new_cols
    
    S.drop(['Optional_diameter', 'Comment', 'Comment_visual', 'Comment_uncertainties', 'Press_min', 'U_10min_knots', 'U_3sec_kntots'],axis=1, inplace=True )
    
    S['datetime'] = pd.to_datetime(S['Time'])
    S.drop(['Time'], axis=1, inplace=True)
    
#    S['Season']= S['Season'].str.slice(0, 4).astype(int)
    S['Month']= S['datetime'].dt.month
    #S = S.set_index('datetime')
    S['ID']= S['ID'].replace({'98': '97'}) #98 is according to Rojo a continuation of 97
    S= S.replace({'comma': 'C', 'undefined': 'U'}) # some wrong morphologies
    S['Morphology']= S['Morphology'].replace({'  ': ' ', 'H - C': 'C - H', 'MGR - C': 'C - MGR', 'T - C': 'C - T', 'U - C': 'C - U', 'H - S':'S - H', 'T - S': 'S - T', 'U - S':'S - U', 'S - W': 'W - S', 'S - MGR': 'MGR - S' }, regex= True)

    
    return S
    

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


def STARS_individual_systems(S, remove_list= ['-', 'T', 'U', 'H']):
    from f_useful import remove_repeater, split_items, remove_from_list, remove_dublicate
    
    S_ind= S.groupby(['ID']).mean()
    S_ind['ID']= S_ind.index
    
    
    """system morphology"""
    ##S = Spiraliform, C = Comma shaped, MGR = Merry-go-round, W = Wave system, U = Undefined, T = Transition between different forms, H = Hybrid, - = Systems don't appear completely on imagery
    
      
    S['Morph_red']= S['Morphology'].replace({'T - ': '', ' - T': '', '  ': ' ', 'U - ': '', ' - U': '', ' - H': '', 'H - ': ''}, regex= True)
    S_ind['PL_Morph_full']= [" ".join(remove_repeater(split_items(remove_repeater(remove_from_list(S[S['ID'] == ID]['Morph_red'], remove_list)), ' - '))) for ID in S_ind['ID']]
    
    #splits transition, remove dublicates, sort by alphabet
    S['Morph_red']= S['Morphology'].replace({'  ': ' '}, regex= True)
    S_ind['PL_Morph']= [" ".join(sorted(remove_dublicate(remove_from_list(split_items(S[S['ID']== ID]['Morph_red'], ' - '), remove_list)))) for ID in S_ind['ID']]
    
    S_ind.loc[S_ind['PL_Morph'].str.contains('MGR'), 'PL_Morph']= 'MGR+'
    S_ind.loc[S_ind['PL_Morph'].str.contains('W'), 'PL_Morph']= 'W+'
    S_ind['PL_Morph']= S_ind['PL_Morph'].replace({'C S': 'C-S'})
    
    
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

    S_ind['Morphology']= S_ind['PL_Morph']

    
    return S_ind



def grad_x(var, lat, lon):
    """calculates the gradient in x-direction of the variable with dimensions(lat, lon) with lat starting from the pole. (The ERA-5 data that I have downloaded)"""
    import xarray as xr
    if type(lat) == xr.core.dataarray.DataArray:
        lat= lat.values
        lon= lon.values
    
    diff_deg= lon[1] - lon[0]
    dx= diff_deg *110E3 * np.cos(np.deg2rad(lat))
    dx= np.tile(dx, (len(lon), 1)).T
    return np.gradient(var, axis= 1)/dx
    
    
def grad_y(var, lat, lon):
    import xarray as xr
    if type(lat) == xr.core.dataarray.DataArray:
        lat= lat.values
        lon= lon.values
        
    dy= (lat[1] - lat[0]) *110E3 
    return np.gradient(var, dy, axis= 0)




def value_in_rad(field, lat, lon, center_lat, center_lon, distance, type= 'mean', intype='pd'):
    """calculuates the 'type'=[mean, max, min, median]' around the center_lat, center_lon, with radius= distance for the field with coordinates (lat, lon)
    intype gives the type of field which can be 'pd' or 'np'"""
    
    if len(np.shape(lat))== 1: #create 2d lat, lon
        lon, lat= np.meshgrid(lon, lat)
    
#    print( center_lat )
#    if type(center_lat) == pd.core.series.Series: #convert the center location from pandas series to a value
#        center_lat, center_lon= center_lat.values[0], center_lon.values[0]
        
    
    dist_now= 110* np.sqrt( (center_lat- lat)**2+ (np.cos(np.deg2rad(center_lat))* (center_lon- lon))**2)
    
   
#    if type(field) is not np.ndarray: field= field.values
#    if not isinstance(type(field), np.ndarray): field= field.values
    if intype=='pd': field= field.values
    
    if type== 'mean': value= np.mean(field[dist_now < distance])
    elif type== 'max': value= np.max(field[dist_now < distance])
    elif type== 'min': value= np.min(field[dist_now < distance])
    elif type== 'med': value= np.median(field[dist_now < distance])

    return value