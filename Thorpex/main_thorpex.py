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
from f_imp_thorpex import * #import the thorpex data

Mediadir= '/media/'+user+'/PatsOrange/'
#homedir= "/home/"+user+"/home/"

"""global variables"""
fignr= 1
#"""Lambert coordinates"""
#lllon, lllat, urlon, urlat= -3, 70, 18, 75
#lat0, lon0= 75, 0 #(lat0, lon0) = center point

#lllon, lllat, urlon, urlat= -5, 65, 28, 76
#lat0, lon0= 70, 0 #(lat0, lon0) = center point

"""b) times"""
year, month = 2008, 3
day, hour= 4, 12                 


if year == 2008 and month== 3 and day== 3 and hour<15:
    lllon, lllat, urlon, urlat= -3, 70, 18, 75
    lat0, lon0= 75, 0 #(lat0, lon0) = center point

elif year == 2008 and month== 3 and day== 3 and hour>=15:
    lllon, lllat, urlon, urlat= -5, 69.5, 10, 74.5
    lat0, lon0= 75, 0 #(lat0, lon0) = center point

elif year == 2008 and month== 3 and day== 4:
    lllon, lllat, urlon, urlat= -1, 64, 12, 70
    lat0, lon0= 70, 0 #(lat0, lon0) = center point   
                    
"""c) plot characteristics"""
PlotPressWind= True

#number_plevels= 20
pleveldist= 1                   

PlotTheta700= False


PlotTheta_e850= True
pot_tempe850_bounds= np.arange(264, 286, 1)


"""d) define vertical cross section"""
startdropnr= 2 #this is the number of the dropsonde (starting from one)
enddropnr= 4 #this is not included

PlotVert= False



"""end global variables"""


"""first run import_thorpex.py"""
thor= data(year, month, day, hour, plevels= np.arange(1000, 450, -10))

startdrop= int(np.where(thor.dropnr== startdropnr)[0]) #this is the index of the dropsonde
enddrop= int(np.where(thor.dropnr== enddropnr)[0]) #this is the index of the dropsonde

plt.figure(fignr)
plt.clf()
fignr+= 1

map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)


"""plot the observation points and their time"""
xpt, ypt= map(thor.lon10, thor.lat10)
map.plot(xpt,ypt,'bo')

for i in range(len(thor.UTC10)):
    plt.text(xpt[i]+100, ypt[i]+100, str(thor.dropnr10[i])+'\n'+thor.UTC10[i])

  
#from scipy.interpolate import interp2d

"""make lat lon grid"""
resol= 30
loni= np.linspace(min(thor.lon10), max(thor.lon10), resol)
lati= np.linspace(min(thor.lat10), max(thor.lat10), resol)
long, latg= np.meshgrid(loni, lati)
lonm, latm= map(long, latg) #translate it to basemap


"""Plot the surface pressure and wind velocity"""
if PlotPressWind== True:
    presi= griddata(thor.lon10, thor.lat10, thor.pres0, long, latg, interp='linear') #interp='nn' is possibly better
    Ui =griddata(thor.lon10, thor.lat10, thor.U10, long, latg, interp='linear') 
    
    PlotContours(lonm, latm, presi, map, leveldist= pleveldist, numbers=False)
    PlotWindVelo(lonm, latm, Ui, map, Umax= 25)
#    PlotLocalMax(presi, threshold=1010, distance=20, map= map, lon=loni, lat=lati, data2=U, threshold2=18, distance2=80/AAres)    

    plt.title('Thorpex')
    plt.tight_layout()


if PlotTheta_e850 == True:
#    thor= data(year, month, day, hour, level=850)
    pressure= 850
    
    
    
    plt.figure(fignr)
    plt.clf()
    fignr+= 1
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)
    
    
    """plot the observation points and their time"""
    latpr= thor.lat[thor.pres== pressure]        
    lonpr= thor.lon[thor.pres== pressure] 
    
    xpt, ypt= map(lonpr, latpr)
    map.plot(xpt,ypt,'bo')
    
    for i in range(len(xpt)):
        plt.text(xpt[i]+100, ypt[i]+100, str(thor.dropnr[np.where(thor.pres== pressure)[0]][i])+'\n'+ str(thor.datetime[thor.pres== pressure][i])[11:16])
        
        
    loni= np.linspace(min(lonpr), max(lonpr), resol)
    lati= np.linspace(min(latpr), max(latpr), resol)
    long, latg= np.meshgrid(loni, lati)
    lonm, latm= map(long, latg) #translate it to basemap

    thetai= griddata(thor.lon[thor.pres== pressure], latpr, thor.theta[thor.pres== pressure], long, latg, interp='linear')
    theta_ei= griddata(lonpr, latpr, thor.theta_e[thor.pres== pressure], long, latg, interp='linear')

    RHi =griddata(lonpr, latpr, thor.RH[thor.pres== pressure], long, latg, interp='linear') ###donot know if this is good

    """Plot the surface pressure"""
#    PlotContours(lonm, latm, presi, map, leveldist= pleveldist)

    """Plot the Theta"""
    #PlotColorMap(lonm, latm, thetai, map,  variable='850 hPa theta')
    PlotColorMap3(lonm, latm, theta_ei, map, bounds= pot_tempe850_bounds, label=r"$\theta_{e,"+str(pressure)+"}$ [K]")
    plt.title('Thorpex')


if PlotVert == True:
    lon0drop= thor.lon[startdrop,-20]
    lat0drop= thor.lat[startdrop,-20]   

    dist= 110* np.sqrt((thor.lat[startdrop:enddrop,-20]- lat0drop)**2+
                       (np.cos(np.deg2rad(lat0drop))*(thor.lon[startdrop:enddrop,-20]- lon0drop))**2)
    plt.figure(fignr)
    plt.clf()
    fignr+= 1
#    PlotCross_sec_T(thor2D.lon[:, 15], thor2D.press[0])
#    plt.contourf(thor2D.theta.T)
    
    """potential temperature"""
    cs= plt.contour(dist, thor.plevels, thor.theta[startdrop:enddrop].T, levels= np.arange(260, 285, 1))
    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f', color= 'black')
    plt.title('Potential temperature')
    
#    plt.plot(dist, [460]*len(dist), 'o')
    for i in range(enddrop- startdrop):
        plt.text(dist[i], 454, str(thor.dropnr[startdrop+ i]), horizontalalignment='center')
    
    plt.ylim([1000, 450])
#    plt.gca().invert_xaxis()

#    """equivalent potential temperature"""
#    cs= plt.contour(np.arange(startdrop, enddrop), np.arange(1000, 250, -10), thor2D.theta_e[startdrop:enddrop].T, levels= np.arange(268, 285, 1))
#    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f', color= 'black')
#    plt.title('Equivalent potential temperature')

    plt.figure(fignr)
    plt.clf()
    fignr+= 1
    cs= plt.contour(dist, thor.plevels, thor.RH[startdrop:enddrop].T)#, levels= np.arange(0, 1, .1))
#    cs= plt.pcolor(dist, thor.plevels, thor.RH[startdrop:enddrop].T)#, levels= np.arange(0, 1, .1))

    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f', color= 'black')
    plt.title('Relative humidity')    

    for i in range(enddrop- startdrop):
        plt.text(dist[i], 454, str(thor.dropnr[startdrop+ i]), horizontalalignment='center')

    """wind speed"""
#    cs= plt.contour(np.arange(startdrop, enddrop), np.arange(1000, 250, -10), thor2D.U[startdrop:enddrop].T, levels= np.arange(0, 30, 2))
#    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f', color= 'black')

    plt.ylim([1000, 450])
#    plt.gca().invert_xaxis()
#  
#    
#    plt.figure(fignr)
#    plt.clf()
#    fignr+= 1
#
#    plt.plot(thor.RH[5])
    
    plt.tight_layout()

