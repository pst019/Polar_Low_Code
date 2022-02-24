#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 15:06:41 2017

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/1692A00D929FEF8B/'


from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import scipy.ndimage.filters as filters

import sys  #to import the functions from a different directory

sys.path.insert(0, homedir +'Polar_Low/Code2/Functions/')

from f_imp_ERA2 import *
from f_imp_ASR import *
from f_impLists import *
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)


fignr= 1

syear= 2000
eyear= 2012 #this is not included
smonths= 1 #including this
emonths= 4 #not including this

#sea of japan
latmax= 48
latmin= 35
lonmin= 130
lonmax= 150 #should rather be 140


#mediterenean sea
latmax= 50
latmin= 30
lonmin= 0
lonmax= 60 

#nordic sea -wilhelmsen
latmax= 80
latmin= 64
lonmin= -12
lonmax= 60 

#STARS list South
#syear= 2006
#eyear= 2012 
#latmax= 64
#latmin= 52
#lonmin= -30
#lonmax= 10 

#smirnova
syear= 2000#1979
eyear= 2009
smonths= 10 #including this
emonths= 4 #not including this
latmax= 80
latmin= 64
lonmin= -10
lonmax= 60 

d= dataASR(syear, month=1, sday=1, level='surf')

#d= data(['SST'], year=syear, month= 4)
nrPLpoints=0 #the number of polar low points
nrPLs= 0

for y, year in enumerate(range(syear, eyear)):
    TPL=TRACK_PLlist_Year(year, name='7', model='ASR')
    
    #ERA
    
    #just take the PL points
    PLlist= TPL.PLlist[:, TPL.PLlist[-1] == 1]
    
    #exclude lat lon region
    PLlist=PLlist[:, PLlist[4] <lonmax]
    PLlist=PLlist[:, PLlist[4] >lonmin]
    PLlist=PLlist[:, PLlist[5] <latmax]
    PLlist=PLlist[:, PLlist[5] >latmin]

    if year==syear: PLlist=PLlist[:, PLlist[1] >= smonths]
    if year== eyear-1: PLlist=PLlist[:, PLlist[1] <emonths]

    PLmonthnumber= PLlist[1]*1E5+PLlist[2] #combine month of the PL and the TRACK number
    nrPLs += len(remove_dublicate(PLmonthnumber))
    
#    print('PLnrs in that area:', remove_dublicate(PLlist[2]) )
    
    for i in range(len(PLlist[0])):
#        print(PLlist[2, i], 'at time:', np.int(PLlist[0, i]), np.int(PLlist[1,i]), np.int(PLlist[3, i]))
#        print(int(PLlist[2, i]), 'at time:', np.int(PLlist[0, i]), np.int(PLlist[1,i]), np.int(PLlist[3, i])//8+1, np.int(PLlist[3, i])%8*3 ,
#              'lon', int(PLlist[4,i]), 'lat', int(PLlist[5,i]) )#, 'Vorticity', PLlist[-5, i])
        nrPLpoints +=1

                      
print('number of PL points:', nrPLpoints)
print('number of PLs:', nrPLs)




