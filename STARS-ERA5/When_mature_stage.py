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
#    homedir= '/home/'+user+'/home/'
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/PatsOrange/'

homedir=Mediadir+'home/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from f_useful import *
from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs
import xarray as xr
from scipy import stats


save= True
#save= False
savedir= homedir + 'Polar_Low/STARS-ana/Figs/SOM_geographic_composites/'

fignr= 1

imp_data=True
#imp_data=False

#maptype= 'PlateCarree'
maptype= 'Polar_Stereo'



x, y= 3, 3
lifelim= 6


node= 1

"""import Stoll list"""
Stoll_imp_dir= Mediadir + '/ERA5_STARS/Stoll_list/'


Stoll= pd.read_csv(Stoll_imp_dir + 'Stoll_list_noRojo.csv')
Stoll= Stoll.set_index(['ID', 'time'])

S_nodes= pd.read_csv(Stoll_imp_dir + 'Stoll_nodes_x'+str(x)+'_y'+str(y)+'.csv')
S_nodes= S_nodes.set_index(['ID', 'time'])

S_ERA5= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5.csv')
S_ERA5= S_ERA5.set_index(['ID', 'time'])

S_ERA5_2= pd.read_csv(Stoll_imp_dir + 'Stoll_ERA5_2.csv')
S_ERA5_2= S_ERA5_2.set_index(['ID', 'time'])

Stoll= pd.concat([Stoll, S_nodes, S_ERA5, S_ERA5_2], axis=1)


Stoll= Stoll.loc[Stoll['Duration'] >= lifelim-1] #Duration gives actually the number of time steps

Stoll= Stoll.reset_index()
Stoll['PLnr']= [int(ID.split('_')[0])*1E3 + int(ID.split('_')[1]) for ID in Stoll.ID.values]



Stoll_ind= Stoll.groupby(['ID']).mean()

Stoll_ind= calc_Obs_mature(Stoll, Stoll_ind, intensity_var='vo')



Stoll_ind['Obs_mature_fraction']= (Stoll_ind['Obs_mature'] -1) / (Stoll_ind['Obs']*2 -1)

plt.figure(fignr)
plt.clf()
plt.scatter(Stoll_ind['Obs_mature_fraction'], len(Stoll_ind)*[0], s= 3)

x= np.linspace(0, 1, 100)
kernel= stats.gaussian_kde(Stoll_ind['Obs_mature_fraction'])
plt.plot(x, kernel(x) )
plt.xlabel('Lifetime fraction of the mature stage')
plt.ylabel('Frequency')



plt.figure(fignr+1)
plt.clf()
var= 'vo_tendency'

plt.scatter(Stoll[var], len(Stoll)*[0], s= 3)
x= np.linspace(Stoll[var].min(), Stoll[var].max(), 100)


kernel= stats.gaussian_kde(Stoll[var])
plt.plot(x, kernel(x) )
plt.xlabel(var)
plt.xlim(x[0], x[-1])
plt.ylabel('Frequency')