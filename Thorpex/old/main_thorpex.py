#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 12:57:58 2016

@author: pst019

Plot the thorpex data

"""


import os
user = os.getcwd().split('/')[2]

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.mlab import griddata #"""interpolate the data"""

import sys  #to import the functions from a different directory
sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_mapplot import * #make basemap plots


Mediadir= '/media/'+user+'/1692A00D929FEF8B/'
#homedir= "/home/"+user+"/home/"

"""global variables"""
fignr= 0
#"""Lambert coordinates"""
lllon, lllat, urlon, urlat= -3, 70, 18, 75
lat0, lon0= 75, 0 #(lat0, lon0) = center point

#lllon, lllat, urlon, urlat= -5, 65, 28, 76
#lat0, lon0= 70, 0 #(lat0, lon0) = center point

"""b) times"""
year, month = 2008, 3
day, hour= 3, 12                 
                    
"""c) plot characteristics"""
PlotPressWind= True

#number_plevels= 20
pleveldist= 1                   

PlotTheta700= False


PlotTheta_e850= False
"""end global variables"""


"""first run import_thorpex.py"""
level='surf'
exec(open("/home/"+user+"/polar_low_code/Thorpex/import_thorpex.py").read())


plt.figure(fignr)
plt.clf()
fignr+= 1

map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)


"""plot the observation points and their time"""
xpt, ypt= map(lon, lat)
map.plot(xpt,ypt,'bo')

for i in range(len(UTC)):
    plt.text(xpt[i]+100, ypt[i]+100, str(i)+'\n'+UTC[i])

  
#from scipy.interpolate import interp2d

"""make lat lon grid"""
resol= 30
loni= np.linspace(min(lon), max(lon), resol)
lati= np.linspace(min(lat), max(lat), resol)
long, latg= np.meshgrid(loni, lati)
lonm, latm= map(long, latg) #translate it to basemap


"""Plot the surface pressure and wind velocity"""
if PlotPressWind== True:
    presi= griddata(lon, lat, pres, long, latg, interp='linear') #interp='nn' is possibly better
    Ui =griddata(lon, lat, U, long, latg, interp='linear') 
    
    PlotContours(lonm, latm, presi, map, leveldist= pleveldist, numbers=False)
    PlotWindVelo(lonm, latm, Ui, map, Umax= 25)
#    PlotLocalMax(presi, threshold=1010, distance=20, map= map, lon=loni, lat=lati, data2=U, threshold2=18, distance2=80/AAres)    

    plt.title('Thorpex')
    plt.tight_layout()


if PlotTheta_e850 == True:
    level=850    #"""import_thorpex.py for a different level"""
    exec(open("/home/"+user+"/polar_low_code/Thorpex/import_thorpex.py").read())

    
    plt.figure(fignr)
    plt.clf()
    fignr+= 1
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
    
    
    """plot the observation points and their time"""
    xpt, ypt= map(lon, lat)
    map.plot(xpt,ypt,'bo')
    
    for i in range(len(UTC)):
        plt.text(xpt[i]+100, ypt[i]+100, UTC[i])


    thetai= griddata(lon, lat, theta, long, latg)#, interp='linear')
    theta_ei= griddata(lon, lat, theta_e, long, latg)#, interp='linear')

    RHi =griddata(lon, lat, RH, long, latg) ###donot know if this is good

    """Plot the surface pressure"""
    PlotContours(lonm, latm, presi, map, leveldist= pleveldist)

    """Plot the Theta"""
    #PlotColorMap(lonm, latm, thetai, map,  variable='850 hPa theta')
    PlotColorMap(lonm, latm, theta_ei, map,  variable=r"$\theta_{e,850}$ [K]")
    plt.title('Thorpex')
