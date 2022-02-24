#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 18:01:07 2020

@author: pst019
"""

Obs_mature= []

for ID_now in remove_dublicate(S.ID):
    S_now= S.loc[S['ID'] == ID_now]
    Obs_mature+= [S_now.iloc[np.argmax(S_now.vo.values)].Obs]