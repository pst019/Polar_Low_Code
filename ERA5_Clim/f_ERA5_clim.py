#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 18:10:28 2021

@author: pst019
"""

import pandas as pd
import numpy as np
from datetime import datetime

from f_useful import *


def import_Rojo_STARS(Mediadir= '/media/pst019/1692A00D929FEF8B/', file= "PL/STARS/Rojo-etal_2019.csv",
                 droplist= ['Optional_diameter', 'Comment', 'Comment_visual', 'Comment_uncertainties', 'Press_min', 'U_10min_knots', 'U_3sec_kntots']):

    """in f_STARS this was just called import_STARS"""
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

    # def _date_parser(*x):
    #     return pd.datetime.strptime(" ".join(x), "%Y %m %d %H %M")

    def _date_parser(*x):
        return datetime.strptime(" ".join(x), "%Y %m %d %H %M")

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


def import_Gunnar_STARS(Mediadir= '/media/pst019/1692A00D929FEF8B/', file_path= 'PL/PLclim/STARS/'):
    """Read both North and South subsets of STARS."""
    STARSdir= Mediadir + file_path
    df_n = read_stars_file(STARSdir + "STARS_TRACKS.csv")
    df_s = read_stars_file(STARSdir + "STARS_TRACKS_south.csv")

    df_s.N += df_n.N.values[-1]

    S= df_n.append(df_s).reset_index(drop=True)
    S= S.rename({'N':'ID'}, axis= 1)
    
    return S




def S_Obs_nr(S):
    """create the Obs nr of each PL"""
    
    S['Obs']= ""

    for ID in remove_dublicate(S['ID']):
        Stime= S.loc[S['ID'] == ID, 'time']
        timelist= remove_dublicate(np.sort(Stime))
        grid= np.meshgrid(timelist, Stime)
        
        labellist= np.argwhere(grid[0] == grid[1])
        S.loc[S['ID'] == ID, 'Obs']= labellist[:, 1] +1
    return S

def S_step(S):
    """create the Obs nr of each PL"""
    
    S['step']= ""

    for ID in remove_dublicate(S['ID']):
        Stime= S.loc[S['ID'] == ID, 'time']
        timelist= remove_dublicate(np.sort(Stime))
        grid= np.meshgrid(timelist, Stime)
        
        labellist= np.argwhere(grid[0] == grid[1])
        S.loc[S['ID'] == ID, 'step']= labellist[:, 1]
    return S