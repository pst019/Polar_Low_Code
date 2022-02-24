#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 13:04:06 2020

@author: pst019
"""

#from xarray get specific index
ds.sel(time=np.datetime64('2001-01-02T17:00:00'))

#from xarray variable
ds.where(ds.PLnr==8001, drop=True).node

#from pandas df
df_nodes.loc['8_1']


S= Stoll.loc[ind]