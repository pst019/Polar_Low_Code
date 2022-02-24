#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

Here are functions to plot fields onto a basemap
"""

import numpy as np
#import os
#from pylab import *
#from datetime import date, datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib as mpl

#from f_plot_on_map import * # not needed in this script, but it was outsourced from here.

print('better use f_plot_fields.py')

def PlotContours(Lon, Lat, psurf, map, nrlevels=10, leveldist=None,levels=None, numbers=True, color= 'k'):
    """ contours for example the pressure
    nrlevels - gives the number of displayed levels
    leveldist - gives distance between levels, if specified the nlevels is ignored
    levels - can be an array that specifies the levels to display, if specified nrlevels and leveldist are ignored
    numbers - True if the contours are labeled
    color - color of the contours (None is s color map)"""
    if levels is not None:
        cs= map.contour(Lon, Lat, psurf, levels, linewidths= 1. , colors= color)
    elif leveldist is not None:
        levels= np.arange(np.round(np.min(psurf)- np.min(psurf)%leveldist), np.round(np.max(psurf)+ leveldist), leveldist)
        cs= map.contour(Lon, Lat, psurf, levels, linewidths= 1. , colors= color)        
    else:
        cs= map.contour(Lon, Lat, psurf, nrlevels, linewidths= 1. , colors= color)#, colors= 6*['b']+ 6*['r'],)
    if numbers == True: plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f')#, color= 'black')
    #plt.tight_layout()


def PlotColorMap(Lon, Lat, w, m, variable='w'):
    """ plot color map, version 3 below is much more flexible , e.g. vertical velocity"""
    boxnr= 21
    colors= [(plt.cm.RdBu_r(h)) for h in range(256)]
    new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=boxnr)
    cs= m.pcolormesh(Lon, Lat, w , cmap=new_map, alpha= 1.)
    cb = m.colorbar(cs,location='right',pad='3%')
    cb.ax.tick_params(labelsize=14)
    
    if variable == 'w': variable='vertical velocity [m/s]'
    cb.set_label(variable, size=14)
    #plt.tight_layout()



def PlotColorMap3(Lon, Lat, data, map, maxlevel= None, symetric=True, bounds=None, label='', color= 'RdBu', boxnr= 21):
    """ plot a color map, e.g. vertical velocity
    if symetric == True it is symetric around 0 and the maxlevel is calculated automatically
    best version of PlotColorMap"""
    
    if color== 'RdBu': colors= [(plt.cm.RdBu_r(h)) for h in range(256)]
    elif color== 'seismic': colors= [(plt.cm.seismic(h)) for h in range(256)]
    elif color== 'blue': colors= [(plt.cm.Blues(h)) for h in range(256)]
    elif color == 'red': colors= [(plt.cm.Reds(h)) for h in range(256)]
    else: print('wrong color')
#    if bounds != None: boxnr = len(bounds)
    new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors) #, N=boxnr)

    if symetric is True and maxlevel is None and bounds is None:
        maxlevel= np.max(np.abs(data))

    if maxlevel is None and bounds is None:
        cs= map.pcolormesh(Lon, Lat, data , cmap=new_map, alpha= 1.)
        cb = map.colorbar(cs,location='right',pad='3%')
    else:
        if bounds is None: bounds= np.round(np.linspace(-maxlevel, maxlevel, boxnr+1), int(np.log10(85/maxlevel)))
#        bounds= np.round(list(np.linspace(-maxlevel, 0, boxnr//2+1))+list(np.linspace(0, maxlevel, boxnr//2+1)), int(np.log10(85/maxlevel)))        
#        print(maxlevel)
        
        norm= mpl.colors.BoundaryNorm(bounds, new_map.N)
        cs= map.pcolormesh(Lon, Lat, data, norm= norm, cmap=new_map, alpha= 1.)
        cb = map.colorbar(cs, boundaries= bounds, norm= norm, location='right',pad='3%')
        
    if label == 'w': label='vertical velocity [m/s]'
    cb.set_label(label, size=14)    
    cb.ax.tick_params(labelsize=14)

    #plt.tight_layout()


def PlotColorMap4(Lon, Lat, data, map, maxlevel= None, symetric=True, bounds=None, label='', color= 'RdBu', boxnr= 21):
    """ plot a color map, e.g. vertical velocity
    if symetric == True it is symetric around 0 and the maxlevel is calculated automatically
    best version of PlotColorMap"""
    
    if color== 'RdBu': colors= [(plt.cm.RdBu_r(h)) for h in range(256)]
    elif color== 'seismic': colors= [(plt.cm.seismic(h)) for h in range(256)]
    elif color== 'blue': colors= [(plt.cm.Blues(h)) for h in range(256)]
    elif color== 'inverse_blue': colors= [(plt.cm.Blues(h)) for h in range(255, 0, -1)]

    elif color == 'red': colors= ['azure']+[(plt.cm.Reds(h)) for h in range(256)]
    else: print('wrong color')
#    if bounds != None: boxnr = len(bounds)
    new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors) #, N=boxnr)

    if bounds is None:
        if maxlevel is not None: minlevel= maxlevel
        if maxlevel is None and bounds is None:
            if symetric is True:
                maxlevel, minlevel= np.max(np.abs(data)), -np.max(np.abs(data))
            else:
                maxlevel, minlevel= np.max(data), np.min(data)       
        
        bounds= np.round(np.linspace(minlevel, maxlevel, boxnr+1), int(np.log10(85/maxlevel)))
#        bounds= np.round(list(np.linspace(-maxlevel, 0, boxnr//2+1))+list(np.linspace(0, maxlevel, boxnr//2+1)), int(np.log10(85/maxlevel)))        
#        print(maxlevel)

    Lon= 0.5* (Lon[1:, 1:]+ Lon[:-1, :-1])
    Lat= 0.5* (Lat[1:, 1:]+ Lat[:-1, :-1])
        
    norm= mpl.colors.BoundaryNorm(bounds, new_map.N)
    cs= map.pcolormesh(Lon, Lat, data[1:, 1:], norm= norm, cmap=new_map, alpha= 1.)
    cb = map.colorbar(cs, boundaries= bounds, norm= norm, location='right',pad='3%')
        
    cb.set_label(label, size=14)    
    cb.ax.tick_params(labelsize=14)



def PlotSurfTemp(Lon, Lat, Tsurf, m, maxlevel= 20, save= False):
    boxnr= 17
    colors= [(plt.cm.RdBu_r(h)) for h in range(256)]
    new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=boxnr)

    Lon= 0.5* (Lon[1:, 1:]+ Lon[:-1, :-1])
    Lat= 0.5* (Lat[1:, 1:]+ Lat[:-1, :-1])

    cs= m.pcolormesh(Lon, Lat, Tsurf[1:, 1:] , cmap=new_map, vmin = -maxlevel, vmax =maxlevel, alpha= 1.) 

    #version2
#    new_map= plt.cm.RdBu_r
#    cs= m.contourf(Lon, Lat, Tsurf , cmap=new_map)
    cb = m.colorbar(cs,location='right',pad='3%')
    cb.ax.tick_params(labelsize=14)
    cb.set_label('Surface Temperature [$^{\circ}$C]', size=14)
    #cb.set_label('Temperature at 2m (C)')
    
    #plt.tight_layout()
    #plt.savefig(savedir+ time +'_Tsurf', dpi = 50)


def PlotIce(Lon, Lat, ice , map, icelimit= 271.5):
    """ plot ice
    icemin - minimum value for which ice is plotted"""
    boxnr= 2
    colors= ['skyblue', 'white'] #[(plt.cm.RdBu_r(h)) for h in range(256)]
    new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=boxnr)
#    cs= map.pcolor(Lon, Lat, ice , cmap=new_map, vmin= icelimit-10, vmax= icelimit+10, alpha= .1) #the "on-off" limit is the avg of vmin and vmax .

    Lon= 0.5* (Lon[1:, 1:]+ Lon[:-1, :-1])
    Lat= 0.5* (Lat[1:, 1:]+ Lat[:-1, :-1])

    cs= map.pcolormesh(Lon, Lat, ice[1:, 1:] , cmap=new_map, vmin= icelimit-5, vmax= icelimit+5, alpha= .5)

#    cb = map.colorbar(cs,location='right',pad='3%')
    
    #plt.tight_layout()
    
    
def PlotPrecip(Lon, Lat, prec, m):
    """ the precipiation of the last hour"""
    boxnr= 5
    colors= [(plt.cm.Blues(h)) for h in range(256)]
    new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors,  N=boxnr)
    
    cs= m.pcolormesh(Lon, Lat, prec , cmap=new_map, vmax= boxnr, alpha= 1.)
    cb = m.colorbar(cs,location='right',pad='3%')
    cb.ax.tick_params(labelsize=14)
    cb.set_ticks(np.arange(boxnr+1))
    cb.set_label('Precipitation [mm/h]', size=14)
    #plt.tight_layout()

def PlotHeatFlux(Lon, Lat, heat, m, label= 'Latent'):
    boxnr= 13
    colors= [(plt.cm.RdBu_r(h)) for h in range(256)]
    new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=boxnr)
    maxlevel= np.max(heat)
    cs= m.pcolormesh(Lon, Lat, heat , cmap=new_map, vmin= -maxlevel, vmax= maxlevel, alpha= 1.)
    cb = m.colorbar(cs,location='right',pad='3%')
    cb.ax.tick_params(labelsize=14)
    
    cb.set_label(label+' heat flux [W/m2]', size=14)
    #plt.tight_layout()
    
def PlotPV(Lon, Lat, PV, m):
    boxnr= 21
    colors= [(plt.cm.seismic(h)) for h in range(256)]
    new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=boxnr)

    maxlevel= np.max(np.abs(PV* 1E6))
    cs= m.pcolormesh(Lon, Lat, PV* 1E6 , cmap=new_map, vmin = -maxlevel, vmax =maxlevel) #, vmin = -maxlevel, vmax =maxlevel, alpha= 1.) 
    cb = m.colorbar(cs,location='right',pad='3%')
    cb.ax.tick_params(labelsize=14)
    cb.set_label('PV [10$^{-6}$ K m$^2$/ kg s]', size=14)
    #plt.tight_layout()

    
def PlotVort(Lon, Lat, vort, map, boxnr= None, maxlevel= None):
#    colors= [(plt.cm.RdBu_r(h)) for h in range(256)]
    
    if maxlevel== None: maxlevel= np.round(np.max(np.abs(vort))-.5)
#    print(maxlevel)
    if boxnr == None: boxnr= 6
#    print(boxnr)
    colors= [(plt.cm.RdBu_r(h)) for h in list(np.arange(0,128, 128//(boxnr -1)))+2*[128]] + [(plt.cm.RdBu_r(h)) for h in list(np.arange(256, 128, -128//(boxnr -1))[::-1])]
    new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=2*boxnr)

    cs= map.pcolormesh(Lon, Lat, vort , cmap=new_map, vmin = -maxlevel, vmax =maxlevel, alpha= 1.) 
    
    cb = map.colorbar(cs,location='right',pad='3%')
    cb.ax.tick_params(labelsize=14)
    cb.set_label('Vorticity [10E-5 1/s]', size=14)
    
    #plt.tight_layout()


def PlotWind(lon, lat, u, v, map, nx= 20, ny=20, rot=True, alen=None, option=None, color=False, scale= 400, colorbarb= 'k'):
    """ winds with barbs, plot this one last, because it is overwritten by other color plots
    alen= arrow length: the lenght of the arrows can be specified by alen
    nx, ny, specifies how many arrows there are in the x and y direction
    for option='ASR' the wind direction has to be opposed, this is done by option ASR"""
#    print('plot winds')
    if len(np.shape(lon))==1: #for ERA data lon and lat are put into 2D array
        lon, lat= np.meshgrid(lon, lat)
        
    everyx= np.max([np.shape(u)[0]//nx, 1]) #put a barb in every nth box, such that you have 20 barbs in the longer dim
    everyy= np.max([np.shape(u)[1]//ny, 1]) #20 for a quadratical map , 50 for era
    u = u[0::everyx, 0::everyy]
    v = v[0::everyx, 0::everyy]
    lon= lon[0::everyx, 0::everyy]
    lat= lat[0::everyx, 0::everyy]
    
    if option=='ASR':
#        print('change sign')
        u, v = -u, -v
    
    if rot:
        u,v, Lon, Lat = map.rotate_vector(u,v,lon,lat ,returnxy=True)
    else:
        Lon,Lat = map(lon,lat)
    #barbs = map.barbs(Lon,Lat,uproj,vproj,length=5,barbcolor='k',flagcolor='r',linewidth=0.5)

    if alen != None:
        print('plot', u.shape)
        Q = map.quiver(Lon,Lat,u,v, pivot= 'mid', scale=scale, color= colorbarb)
        qk = plt.quiverkey(Q, 0.12, 0.97, alen, str(alen)+' m/s', labelpos='W')
    
    else:        
        Umax= np.max(np.sqrt(u**2+v**2)) #//5*5 #rounds to nearest 5
#        print(Umax)
        if Umax < 40:
#            if color==True:
#                Q = map.quiver(Lon,Lat,u,v, color= np.sqrt(u**2+v**2)/Umax, pivot= 'mid', scale=500)
#                qk = plt.quiverkey(Q, 0.22, 0.97, 20, '20 m/s', labelpos='W')         
#            else:
            Q = map.quiver(Lon,Lat,u,v, pivot= 'mid', scale=scale, color= colorbarb)
            qk = plt.quiverkey(Q, 0.12, 0.97, 20, '20 m/s', labelpos='W')
        else:
            Q = map.quiver(Lon,Lat,u,v, pivot= 'mid', scale=1000, color= colorbarb)
            qk = plt.quiverkey(Q, 0.12, 0.97, 75, '75 m/s', labelpos='W')
            
#        Q = map.quiver(Lon,Lat,uproj,vproj, pivot= 'mid', scale=25*Umax)
#        qk = plt.quiverkey(Q, 0.22, 0.97, Umax, str(int(Umax))+' m/s', labelpos='W')
    #plt.tight_layout()


def PlotWindVelo(Lon, Lat, U, map, Umax=None):
    if Umax== None: Umax= np.round(np.max(U))
    PlotColorMap4(Lon, Lat, U, map, bounds= np.arange(0, Umax +1, 2), color='blue', label='Wind velocity [m/s]')


#def PlotWindVelo(Lon, Lat, U, m, Umax=None, bar=True, plotbar=False):
#    """wind speed with Teresas blue colors
#    bar - specifies if the colorbar should be given in plot
#    plotbar, if bar=False - specifies if the bar should be plotted in an extra graph"""
#    if Umax != None:
#        boxnr= Umax//5
#        if boxnr >=10: boxnr //= 2
#        colors= [(plt.cm.YlGnBu(h)) for h in range(256)]
#        new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=boxnr)
#    else:
#        new_map= plt.cm.YlGnBu
#    
#    Lon= 0.5* (Lon[1:, 1:]+ Lon[:-1, :-1])
#    Lat= 0.5* (Lat[1:, 1:]+ Lat[:-1, :-1])
#    
#    cwind = m.pcolormesh(Lon, Lat, U[1:, 1:], cmap=new_map, vmin= 0, vmax= Umax)
#    #plt.tight_layout()
#    if bar== True:
#        cb = m.colorbar(cwind,"right", size="5%", pad="2%")
#        cb.set_label('Wind speed (m/s)')
#    else:
#        if plotbar==True:
#            plt.figure(10)
#            plt.clf()
#            cax = plt.axes([0.1, 0.2, 0.8, 0.05])
#            cb = plt.colorbar(cwind, orientation='horizontal', cax= cax)
#            cb.set_label('Wind speed (m/s)', size= 20)        
#            cb.ax.tick_params(labelsize=18)
##        #plt.tight_layout()

    
def PlotStaticStability(Lon, Lat, Tdiff, m, difftype='SST-T500'):
    """Tdiff = Tsurf - Th
    difftype =['SST-T500', thetaSST-theta500']"""
    
    if difftype== 'thetaSST-theta500': lev= [-13, -9, -6, -3, 0, 3]
    else: lev= [40, 43, 46, 49, 52]
    cs= m.contour(Lon, Lat, Tdiff, levels= lev, linewidths= 1. , colors= 'red' )#, colors= 6*['b']+ 6*['r'],)

#    col= [(plt.cm.YlOrBr(h)) for h in np.linspace(0, 256, 9)]
#    col= [(plt.cm.YlOrBr(h)) for h in range(0, 256, 256//10)]
#    cs= m.contour(Lon, Lat, Tdiff, levels= [38, 40, 42, 44, 46, 48, 50, 52, 54], colors=col,  linewidths= 1. )#, colors= 'k' )#, colors= 6*['b']+ 6*['r'],)

    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f', color= 'black')
    #plt.tight_layout()
    
    

def PlotContourStatStabTheta_e(Lon, Lat, Tdiff, m):
    """Tdiff = Tsurf - Th"""
    levlist=[-16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8]
    levlist=[-15, -12, -9, -6, -3, 0, 3, 6, 9]
    cs= m.contour(Lon, Lat, Tdiff, levels= levlist, linewidths= 1. )#, colors= 'k' )#, colors= 6*['b']+ 6*['r'],)
    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f', color= 'black')
    #plt.tight_layout() 
    
def PlotColor_from0(Lon, Lat, data, map, boxnr=11, maxlevel= None, log=False, label= None, col='YlOrBr'):
    """plot the count density of polar lows - have to work on the box nr
    maxlevel can specify the maximum level in the plot - if not given it is calculated by itself by the max of the data
    log - if the count density should be plotted logaritmic or not
    col= [YlOrBr , BlueRed] defines the used colors
    mainly replaced by PlotColorMap_3"""
#    boxnr= np.max(data*fact)+1
    
    if col=='BlueRed': colors= [(plt.cm.RdBu_r(h)) for h in range(256)]
    else: colors= ['white']+[(plt.cm.YlOrBr(h)) for h in range(256)]

    new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=boxnr)
    
    if maxlevel== None: maxlevel= np.max(np.abs(data))
    
    if log==False:    
        cs= map.pcolormesh(Lon, Lat, data, cmap=new_map, vmin= 0, vmax= maxlevel, alpha= 1.)
    else:
        from matplotlib.colors import LogNorm
        cs= map.pcolor(Lon, Lat, data, norm= LogNorm(vmin=0.005, vmax= maxlevel), cmap=new_map)
    
    cb = map.colorbar(cs,location='right',pad='3%')   
    
    if log == False: cb.set_ticks(np.linspace(0, maxlevel, boxnr+1))
#    tick_locator = ticker.MaxNLocator(nbins=boxnr-3)
#    cb.locator = tick_locator
#    cb.update_ticks()
    if label != None: cb.set_label(label)



def PlotColor_setbounds(Lon, Lat, data, map, bounds, label= None, col='YlOrBr'):
    """plot with setting your own boundaries, for example: bounds= [0, 0.1, 0.3, 1, 5]
    col= [YlOrBr , BlueRed] defines the used colors
    mainly replaced by PlotColorMap_3"""
#    boxnr= np.max(data*fact)+1
    
    if col=='BlueRed': colors= [(plt.cm.RdBu_r(h)) for h in range(256)]
    elif col=='Blues': colors= [(plt.cm.YlGnBu(h)) for h in range(256)]
    else: colors= ['white']+[(plt.cm.YlOrBr(h)) for h in range(256)]

    new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=len(bounds))
    norm= mpl.colors.BoundaryNorm(bounds, new_map.N)
    
    cs= map.pcolormesh(Lon, Lat, data, norm= norm, cmap=new_map, alpha= 1.)
    
    cb = map.colorbar(cs, boundaries= bounds, norm= norm, location='right',pad='3%')   
    
    if label != None: cb.set_label(label, size=14)
    cb.ax.tick_params(labelsize=14)



    

    

    

def PlotCross_sec_T(press, lon, field):
    """now for temp"""
    nrlevels= 10 #np.arange(960,1080,4)
    cs= plt.contour(lon, press , field, nrlevels, linewidths= 1. , colors= 'k' )
    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f', color= 'black')
    plt.gca().invert_yaxis()
    plt.ylabel('Pressure [hPa]')
    plt.xlabel('Longitude [deg east]')
    #plt.tight_layout()
     
def PlotCross_sec_hum(press, lon, field):
    "for rel hum"""
#    cs= plt.contour(lon, press , field, linewidths= 1. , colors= 'k' )
    cs= plt.contourf(lon, press , field, linewidths= 1. , cmap=plt.cm.YlGnBu, levels=np.linspace(0,1, 11))
    cb = plt.colorbar()#cs,location='right',pad='3%')
    cb.set_label('relative humidity', size=14)
#    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.1f', color= 'black')
#    plt.gca().invert_yaxis()
    plt.ylabel('Pressure [hPa]')
    plt.xlabel('Longitude [deg east]')
    #plt.tight_layout()

def PlotCross_sec_vervel(press, lon, w):
    """plot the vertical velocity"""
    ze= np.zeros(np.shape(w))
    
    Q= plt.quiver(lon, press, ze, w, scale= 10)
    plt.quiverkey(Q, 0.15, 0.97, 1, '1 m/s', labelpos='W')
    plt.title('vertical velocity')
    
    plt.gca().invert_yaxis()
    #plt.tight_layout()
    
def PCross_sec_horvel(press, lon, u, v):
    """plot the horizontal velocity with north south (up), east west (right) wind on a cross section"""
    """some areas are white because of dicontinuities"""
    U = np.sqrt(u**2 + v**2)
    cs= plt.contourf(lon, press , U, cmap=plt.cm.YlOrRd)
    cb = plt.colorbar()#cs,location='right',pad='3%')
    cb.set_label('velocity', size=14)
    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f', color= 'black')
    
    Q= plt.quiver(lon, press, u, v)
    plt.quiverkey(Q, 0.15, 0.97, 20, '20 m/s', labelpos='W')
    plt.title('horizontal velocity')
       
    plt.ylabel('Pressure [hPa]')
    plt.xlabel('Longitude [deg east]')
        
    plt.gca().invert_yaxis()
    plt.ylim([1030, 0])
    #plt.tight_layout()       
    
def PCross_sec_horvel2(press, lon, u, v):
    """plot the horizontal velocity with north south (up), east west (right) wind on a cross section"""
    """as PCross_sec_horvel but with plot for every box so that nothing gets white"""
    boxnr= 7
    colors= [(plt.cm.YlOrRd(h)) for h in range(256)]
    new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors,  N=boxnr)

    U = np.sqrt(u**2 + v**2)

    #cs= plt.contourf(d.lon[y,x], d.press, U, cmap=new_map) #to discontinous for the data
    cs= plt.pcolormesh(lon, np.append(press, 0), U, cmap=new_map)
    
    cb = plt.colorbar()#cs,location='right',pad='3%')
    cb.set_label('velocity', size=14)
    #plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f', color= 'black') #together with contourf
    
    Q= plt.quiver(lon, press, u, v)
    plt.quiverkey(Q, 0.15, 0.97, 20, '20 m/s', labelpos='W')
    #plt.barbs(d.lon[y,x], d.press, ze, d.w*100)
    plt.title('horizontal velocity')
    
    plt.ylabel('Pressure [hPa]')
    plt.xlabel('Longitude [deg east]')
    
    plt.gca().invert_yaxis()
    plt.ylim([1030, 0])
    ##plt.tight_layout()
    
    
    

    