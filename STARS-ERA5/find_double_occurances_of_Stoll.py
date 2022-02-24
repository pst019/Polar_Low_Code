#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 10:18:41 2020

@author: pst019
"""

#import the list from somewhere else
S= Stoll

time_list= remove_dublicate(S.time)

for t in time_list:
    S_now= S.loc[S['time']==t]
    if len(S_now['Stoll nr']) != len(remove_dublicate(S_now['Stoll nr'])):
        print(t, S_now['Stoll nr'])
        