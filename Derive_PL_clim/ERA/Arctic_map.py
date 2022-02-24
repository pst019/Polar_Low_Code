#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 13:25:51 2017

@author: pst019
from ASR-ERA_Inv_random_3 of folder ASR
"""

import sys
import os
user = os.getcwd().split('/')[2]
#print('user: ', user)

homedir= '/home/'+user+'/home/'

sys.path.insert(0, code/')
from f_meteo import *

#sys.path.insert(0, code/Derive_PL_clim/ERA/')
#from f_imp_ERA2 import *

sys.path.insert(0, homedir+'/MaPro/Code_office/paper/')
from funclass5 import *

sys.path.insert(0, code/AROME/')
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)

fignr= 1


hemi= 'NH'
#hemi='SH'


plt.figure(fignr)
fignr +=1
plt.clf()
#plt.subplot(2, 3, 2)


if hemi== 'NH':
    arc= region(latmin= 40, latmax= 90, lonmin= -180, lonmax= 180, model= 'ERA', name= 'ARC')
elif hemi== 'SH':
    arc= region(latmin= -90, latmax= -40, lonmin= -180, lonmax= 180, model= 'ERA', name= 'ANT')

arc.data(origdata=True, dataset= 'SIC', save=False) #, regresstype= 'quad')

if hemi== 'NH':
    map= Polar_map(latsouth= arc.lat[0], fillcontinents= True)
elif hemi =='SH':
    map= Polar_map(latsouth= arc.lat[-1], fillcontinents= True, hemisp= hemi)

#map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
grid= np.meshgrid(arc.lon, arc.lat) #pcolormesh needs corner points for projection
Lon, Lat= map(grid[0], grid[1])

#de= data(['SST'], year, month)
#de.SST[de.SST<272]= np.nan #this sets the area with ice to land
#de.SST[de.SST>=272]= -.3 #this sets the area with ice to land

#de.SST= ma.array(de.SST, mask= isnan(de.SST))

if hemi== 'NH':
    SICaverage= np.mean(np.concatenate((arc.origdata['SIC'][:, :4], arc.origdata['SIC'][:, 10:]), axis= 1), axis= (0,1))
elif hemi =='SH':
    SICaverage= np.mean(arc.origdata['SIC'][:, 4:11], axis= (0,1))


SICaverage[SICaverage >=0.2]= 1
SICaverage[SICaverage <0.2]= -1
          
PlotIce(Lon, Lat, SICaverage, map, icelimit= 0)


if hemi == 'NH':
    x,y = map(-90,85)
    plt.text(x,y,'Sea Ice')
    
    x,y = map(-40,58)
    plt.text(x,y,'Denmark \n Strait')
    
    #x,y = map(0,71)
    #plt.text(x,y,'Tromsø \n Flake')
    
    x,y = map(-3,66)
    plt.text(x,y,'Nordic \n Seas')

    x,y = map(18,73)
    plt.text(x,y,'Barents \n Sea')
    
    x,y = map(-55,55)
    plt.text(x,y,'Labrador \n Sea')
    
    x,y = map(-130,55)
    plt.text(x,y,'Gulf of \nAlaska')
    
    x,y = map(140,50)
    plt.text(x,y,'Sea of \n Japan')
    
    x,y = map(-42,42)
    plt.text(x,y,'North Atlantic')
    
    x,y = map(-169,45)
    plt.text(x,y,'North Pacific')
    
    x,y = map(160,60)
    plt.text(x,y,'Sea of\nOkhotsk')
    
    x,y = map(-172,60)
    plt.text(x,y,'Bering \n Sea')  
#    x,y = map(-120,70)
#    plt.text(x,y,'Beaufort \n Sea')
    
if hemi == 'SH':
    x,y = map(-90, -82)
    plt.text(x,y,'Antarctica')  
    
    x,y = map(-35, -62)
    plt.text(x,y,'Sea Ice')  

    x,y = map(-35, -42)
    plt.text(x,y,'South Atlantic')  

    x,y = map(-135, -42)
    plt.text(x,y,'South Pacific')  

    x,y = map(53, -60)
    plt.text(x,y,'Indian Ocean')  

    x,y = map(-88, -45)
    plt.text(x,y,'Bellingshausen \n Sea') 

    x,y = map(-108, -47)
    plt.text(x,y,'Amundsen \n Sea')

    x,y = map(122, -57)
    plt.text(x,y,'Mawson \n Sea')

    x,y = map(95, -57)
    plt.text(x,y,'Davis \n Sea')

    x,y = map(178, -43)
    plt.text(x,y,'New \nZealand')
    
    x,y = map(-74, -42)
    plt.text(x,y,'South \nAmerica')
    

#plt.tight_layout()
