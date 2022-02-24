#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 22:46:08 2018

@author: pst019
"""

import numpy as np
import matplotlib.pyplot as plt
from f_impLists import *

def mapcross_section(dlon, dlat, x_start, y_start, x_end, y_end, map, extent= 0, coordinates='xy', color='k', markersize= 2):
    """plot a cross section onto the map
    extent- specifies an extention in % on each side of the start and end point.
    coordinates=['xy' or 'latlon']- specifies if the start and endpoints are given in pixel (xy) or in latitude-longitude (latlon)
    xy: xstart - startpoint in x direction
    latlon: xstart - startpoint in lon direction
    it gives out the x, y points that are needed for the imp_cross_section"""
    
    if coordinates== 'xy':
        lon= dlon[[x_start, x_end], [y_start, y_end]]
        lat= dlat[[x_start, x_end], [y_start, y_end]]
    elif coordinates== 'latlon':
        lon_start, lon_end, lat_start, lat_end= x_start, x_end, y_start, y_end
        lon= np.array([lon_start, lon_end])
        lat= np.array([lat_start, lat_end])
        
        dist= np.sqrt((np.cos(np.deg2rad(lat_start))*(dlon - lon_start))**2+ (dlat - lat_start)**2)
        x_start, y_start= np.where(dist== np.min(dist))
        
        dist= np.sqrt((np.cos(np.deg2rad(lat_end))*(dlon - lon_end))**2+ (dlat - lat_end)**2)
        x_end, y_end= np.where(dist== np.min(dist))
        
        x_start, y_start, x_end, y_end= int(x_start), int(y_start), int(x_end), int(y_end)
        
    
    xpt, ypt= map(lon, lat)
#    map.plot(xpt,ypt,'o') #plot the start and end point
    #plt.text(xpt[i]+100, ypt[i]+100, UTC[i])
    
    
    res= int(np.sqrt((x_start- x_end)**2 + (y_start - y_end)**2))
    print(res, x_start, x_end)
    
    #this is the way with extent= 0
#    x= [int(round(x_start+ (x_end - x_start)*i/res )) for i in range(res+1)]
#    y= [int(round(y_start+ (y_end - y_start)*i/res )) for i in range(res+1)]

    x= [int(round(x_start+ (x_end - x_start)*i/res )) for i in range(- int(extent*res/100), res+ int(extent*res/100) +1)]
    y= [int(round(y_start+ (y_end - y_start)*i/res )) for i in range(- int(extent*res/100), res+ int(extent*res/100) +1)]
    
    lon= dlon[x, y]
    lat= dlat[x, y]
    xpt, ypt= map(lon, lat)
    map.plot(xpt,ypt,'o', markersize= markersize, color=color)
    return x, y
    

def PointsOnMap(maxpoints, lon, lat, map, PLnumber=[]):
    """plot points onto the map
    maxpoints is a tuple/array with latcoord in first row and loncoord in second row
    write the PLnumber next to the point if given"""
    lonp= lon[maxpoints[1]]
    latp= lat[maxpoints[0]]
    xpt, ypt= map(lonp, latp)
    map.plot(xpt,ypt,'bo')
    
    if len(PLnumber) > 0:
        for n in range(len(lonp)):
            plt.text(xpt[n], ypt[n], str(int(PLnumber[n])), fontsize=10, ha='center',va='top',color='black')
    
       
    
def CurrentTrack(t, PLnumber, PLlistraw, map, d):
    """plot the polar tracks at time t on map"""
    PLnrs= PLnumber[PLlistraw[0]== t]
    print(PLnrs)
    c= ['b', 'g', 'r', 'y', 'm', 'c', 'k']*len(PLnumber)*2
    for n in PLnrs:
        lonp= d.lon[PLlistraw[2][PLnumber==n]]
        latp= d.lat[PLlistraw[1][PLnumber==n]]
        xpt, ypt= map(lonp, latp)
        map.plot(xpt,ypt,'x', color=c[int(n-1)])
        map.plot(xpt,ypt,'-', color=c[int(n-1)])
        map.plot(xpt[0],ypt[0],'o', label= str(int(n)), color=c[int(n-1)])#, markersize= 10)
        
        PLdurationtotime= np.where(PLlistraw[0][PLnumber==n]==t)[0][0]  #the amount of time steps the pl has already lived at this point of time
        map.plot(xpt[PLdurationtotime],ypt[PLdurationtotime],'^', color=c[int(n-1)])#, markersize= 10)
    
        plt.text(xpt[0],ypt[0], str(int(n)), fontsize=10, #fontweight= 'bold',
                            ha='center',va='top',color='black')

        
def PlotMonthsPLtracks(PLnumber, PLlistraw, map, d):
    c= ['b', 'g', 'r', 'y', 'm', 'c', 'k']*len(PLnumber)*2
    
    for n in remove_dublicate(PLnumber)[:50]:
        lonp= d.lon[PLlistraw[2][PLnumber==n]]
        latp= d.lat[PLlistraw[1][PLnumber==n]]
        xpt, ypt= map(lonp, latp)
        map.plot(xpt,ypt,'x', color=c[int(n-1)])
        map.plot(xpt,ypt,'-', color=c[int(n-1)])
        map.plot(xpt[0],ypt[0],'o', label= str(int(n)), color=c[int(n-1)])#, markersize= 10)
        plt.text(xpt[0],ypt[0], str(int(n)), fontsize=10, #fontweight= 'bold',
                            ha='center',va='top',color='black')
   
 
    
def PlotTRACK_PL(TPL, t, map, track=True, nowPL=True, number=True, othert=False, color='r'):
    """TRACK polar lows -now
    TPL - list of TRACK PLs, t - considered 6hourly timepoint
    used by Inv_timeTandS_3
    if track==False, than only the actual PL point is plottet on the map
    if nowPL= False (best with track=True), shows tracks of cyclones, that get PLs at other points in time are shown
    if othert is a time value, than it plots the point where the cyclone is at othert that is a PL at t"""
    PLt= TPL.PLlist[1]== t
    
    if nowPL== True:                   
        PL1= TPL.PLlist[-1]==1
        PLnow= np.array([a and b for a, b in zip(PLt, PL1)])
    else:
        PLnow=PLt
        
    PLnrs= TPL.PLlist[0, PLnow]
    
    TPLxpt, TPLypt= map(TPL.PLlist[2][PLnow], TPL.PLlist[3][PLnow])
    if othert== False:
        map.plot(TPLxpt,TPLypt,'o', color=color)
    
    if number is True:
        for n in range(len(TPLxpt)):
            plt.text(TPLxpt[n],TPLypt[n], str(int(TPL.PLlist[0][PLnow][n])), fontsize=10, ha='center',va='top',color='green')

    if track is True:        
        for n in PLnrs:
            lonp= TPL.PLlist[2][TPL.PLlist[0]==n]
            latp= TPL.PLlist[3][TPL.PLlist[0]==n]
            xpt, ypt= map(lonp, latp)
            map.plot(xpt,ypt,'gx') 
            map.plot(xpt,ypt,'g-') 
            map.plot(xpt[0],ypt[0],'g^', label= str(int(n))) 
        
            lonp=lonp[TPL.PLlist[-1][TPL.PLlist[0]==n]==1]
            latp=latp[TPL.PLlist[-1][TPL.PLlist[0]==n]==1]
            xpt, ypt= map(lonp, latp)
            map.plot(xpt,ypt,'yx')
            
            if othert != False:
                lonp= TPL.PLlist[2][np.logical_and(TPL.PLlist[0]==n, TPL.PLlist[1]==othert)]
                latp= TPL.PLlist[3][np.logical_and(TPL.PLlist[0]==n, TPL.PLlist[1]==othert)]
                xpt, ypt= map(lonp, latp)
                map.plot(xpt,ypt,'o', color=color)                 
#            map.plot(xpt,ypt,'y-') 


def PlotallTRACKpoints(T, tcomb, map, text=False, color='b'):
    """Plotting all TRACK cyclones points happening now
    if text=True it displays the number of the TRACK cyclone"""
    Txpt, Typt= map(T.PLlist[2][T.PLlist[1]== tcomb], T.PLlist[3][T.PLlist[1]== tcomb])
    print(Txpt)
    map.plot(Txpt,Typt,'o', color=color)
    if text==True:
        for n in range(len(Txpt)):
            plt.text(Txpt[n],Typt[n], str(int(T.PLlist[0][T.PLlist[1]== tcomb][n])), fontsize=10, ha='center',va='top',color='black')


def PlotallTRACK_PLpoints(TPL, t, map, color='y'):
    """TRACK polar lows points happening now in yellow"""
    PLnow= np.where(np.logical_and(TPL.PLlist[1]== t, TPL.PLlist[-1]==1))
    TPLxpt, TPLypt= map(TPL.PLlist[2][PLnow], TPL.PLlist[3][PLnow])
    map.plot(TPLxpt,TPLypt,'o', color=color)

def PlotGivenTRACKcyclones(T, TPLnrs, map, color='b'):
    """Plot the TRACKs with numbers TPLnrs (list) in blue"""
    for n in TPLnrs:
        lonp= T.PLlist[2][T.PLlist[0]==n]
        latp= T.PLlist[3][T.PLlist[0]==n]
        xpt, ypt= map(lonp, latp)
        map.plot(xpt,ypt,'x', color=color) 
        map.plot(xpt,ypt,'-', color=color) 
        map.plot(xpt[0],ypt[0],'^', label= str(int(n)), color=color) 
        plt.text(xpt[0], ypt[0], str(int(n)), fontsize=10, ha='center',va='top',color=color)

def PlotPLpointsinTRACKcyclone(TPL, TPLnrs, map, color='y'):
    """mark the PL points of the TRACKS cyclones with numbers TPLnrs yellow"""
    for TPLnr in TPLnrs:
        index= np.where(np.logical_and(TPL.PLlist[0]== TPLnr, TPL.PLlist[-1]==1))[0]
        TPLxpt, TPLypt= map(TPL.PLlist[2][index], TPL.PLlist[3][index])
        map.plot(TPLxpt,TPLypt,'x', color=color)

    
def PlotSTARS_PL(S, tcomb, map):
    """Plot all STARS polar lows happening now with their tracks - return the number of the STARS PLs
    used by Inv_timeTandS_3"""
    Sxpt, Sypt= map(S.lon[S.tPL== tcomb], S.lat[S.tPL== tcomb])
    map.plot(Sxpt,Sypt,'ro')
    
    for n in range(len(Sxpt)):
        plt.text(Sxpt[n],Sypt[n], str(int(S.PLnumber[S.tPL== tcomb][n])), fontsize=10, ha='center',va='top',color='red')
    
    
    PLnrs= S.PLnumber[S.tPL== tcomb]
    for n in PLnrs:
        lonp= S.lon[S.PLnumber==n]
        latp= S.lat[S.PLnumber==n]
        xpt, ypt= map(lonp, latp)
        map.plot(xpt,ypt,'rx') 
        map.plot(xpt,ypt,'r-') 
        map.plot(xpt[0],ypt[0],'r^', label= str(int(n)))
        
    return PLnrs

def PlotSTARSandMatchingTRACK(S, T, TPL, t, tcomb, map):
    """Plot the STARS PL, the corresponding TRACK cyclone, 
    all TRACK and TRACKPL points now and PLpoints during the corresponding TRACK cyclone
    used by Inv_timeTandS_4"""
    SPLnrs= PlotSTARS_PL(S, tcomb, map)
    
    """Plotting all TRACK cyclones points happening now"""
    PlotallTRACKpoints(T, tcomb, map, text=False)
    
    """TRACK polar lows points happening now"""
    PlotallTRACK_PLpoints(TPL, t, map)
    
    TPLnrs= MatchingTRACKtoSTARS(SPLnrs) #find the TRACK cyclones corresponding to the STARS PLs
    
    """Plot the TRACKs of the corresponding TRACKS cyclones"""
    PlotGivenTRACKcyclones(T, TPLnrs, map)
    
    """mark the PL points of the corresponding TRACKS cyclone yellow"""
    PlotPLpointsinTRACKcyclone(TPL, TPLnrs, map)
    
    
    
def PlotLocalMax(data_in, threshold, distance, map, lon, lat, typ='min',
                 data2=None, threshold2=None, distance2=None, typ2='max',
                 color='r', dot=True, value=True, yoffset= None, roundorder= 1,
                 latbound= None, lonbound= None):
    """ plot the local max/min (typ) of data, below/above threshold (scalar)
    with a minimum distance of distance (scalar or tuple of shape(data)) to the closest next local max/min
    on map with lon/lat data
    data2 can give a second criteria to be fullfill threshold2 within distance2
    dot - specifies if a dot is plotted in map for the location of the max/min
    yoffset - displacement of the dot in y direction
    value - specifies if the value should be plotted next to the dot
    roundorder - gives how many digits behind the comma are displayed
    latbound -[latmin, latmax] - picks only local maxima  with latmin < lat < latmax
    """   

#    print('if there is a problem with tight layout, specify latbound, lonbound')

    import scipy.ndimage.filters as filters

    data = np.copy(data_in)

    if data2 is not None and threshold2 is None:
        if typ2== 'max':
            threshold2 = np.min(data2)
        elif typ2== 'min':
            threshold2= np.max(data2)
    if data2 is not None and distance2 is None:
        distance2= distance

    if lon.ndim==1 and data.ndim==2:
        lon, lat= np.meshgrid(lon, lat)
        
    if latbound is not None:
        data[lat< latbound[0]]=0
        data[lat> latbound[1]]=0
        
    if lonbound is not None:
        data[lon< lonbound[0]]=0
        data[lon> lonbound[1]]=0

    if typ== 'min': data_max = filters.minimum_filter(data, distance)
    elif typ == 'max':
        
        def maximum_filter_ignore_nan(array, *args, **kwargs):
            nans = np.isnan(array)
            replaced = np.where(nans, 0 , array)
            return filters.maximum_filter(replaced, *args, **kwargs)

        data_max = maximum_filter_ignore_nan(data, distance)

#        data_max = filters.maximum_filter(data, distance)

    maxima = (data == data_max)
    if typ== 'min':
        maxima[data_max > threshold]=0
    elif typ== 'max':
        maxima[data_max < threshold]=0
              
    if data2 is not None:
        data2_max = filters.maximum_filter(data2, distance2)
        if typ2== 'max': maxima[data2_max < threshold2]= 0
        elif typ2== 'max': maxima[data2_max > threshold2]= 0

#    print(maxima)    
#    print(np.where(maxima==True))
#    print(lon[maxima])
    xp, yp= map(lon[maxima], lat[maxima])    
    if dot==True and yoffset is None: yoffset = 0.05*(map.ymax-map.ymin)
    else: yoffset= 0
    
    for i in range(len(xp)):
        if dot==True: map.plot(xp[i], yp[i], 'o', color=color)
        if value== True:
            plt.text(xp[i], yp[i]+yoffset, ("{:."+str(roundorder)+"f}").format(round(data[maxima][i], roundorder)), fontsize=13, fontweight= 'bold',
                                ha='center',va='top',color=color)#bbox = dict(boxstyle="square",ec='None',fc=(1,1,1,0.5)))
    
    return data[maxima], lon[maxima], lat[maxima]
#    maxima_wind = (U == data_max_wind)
#    maxima_wind[data_max > threshold]=0
##    maxima[d.lon < -10] = 0
#    maxima_wind[data_max_wind < Uthreshold]= 0
#          
#    xp, yp= map(d.lon[maxima_wind], d.lat[maxima_wind])    
#    for i in range(len(xp)):
##        map.plot(xp[i], yp[i], 'o', color='black')
#        plt.text(xp[i], yp[i], str(round(U[maxima_wind][i], 1)), fontsize=13, fontweight= 'bold',
#                            ha='center',va='top',color='b',
#                            )#bbox = dict(boxstyle="square",ec='None',fc=(1,1,1,0.5)))

    