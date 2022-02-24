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

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'

import matplotlib.pyplot as plt

sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')

from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_meteo import *
import f_imp_thorpex as thorpex


"""global variables"""

fignr= 5

#maptype='AA'
maptype='AA_half'
#maptype='Lambert'
#maptype='polar'

if maptype=='AA': plt.figure(fignr, figsize= (6, 6))
elif maptype=='AA_half': plt.figure(fignr, figsize= (6, 4.7))
else: plt.figure(fignr, figsize= (6, 4.5))
fignr += 1
plt.clf()



"""time"""
year=2008
month=3
day= 3
shour= 0

lacktime=30 #in hours

t= lacktime//3


"""Plot the map"""
if maptype== 'AA': map = AA_map()
elif maptype== 'AA_half': map = AA_map_half()
elif maptype=='polar': map= Polar_map(latsouth= 50) #d.lat[-1])
else: 
    plotday= day+ (shour+lacktime)//24
    plothour= (shour+lacktime)%24
    #to match with Thorpex                    
    if year == 2008 and month== 3 and plotday== 3 and plothour<15:
        lllon, lllat, urlon, urlat= -4, 70, 19, 75
        lat0, lon0= 75, 0 #(lat0, lon0) = center point
    
    elif year == 2008 and month== 3 and plotday== 3 and plothour>=15:
        lllon, lllat, urlon, urlat= -6, 69.5, 13, 74.5
        lat0, lon0= 75, 0 #(lat0, lon0) = center point
    
    elif year == 2008 and month== 3 and plotday== 4:
        lllon, lllat, urlon, urlat= -3, 63.5, 15, 69.5
        lat0, lon0= 70, 0 #(lat0, lon0) = center point 
    
    map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)



#var= 'PressWind'
var= 'PressWind_advanced'
#var= 'RH_lev'
#var= 'SH_lev'
#var= 'Theta_lev'
#var= 'Geop_wind_lev'
#var= 'OLW'
#var= 'Vort_lev'

ilevel= 1 #0 -950, 1- 850, 2-700, 3-500


Thorpex= False
Plotbox= False
#boxlon, boxlat= [-4, 13, 18.5, -4, -4], [70, 69.3, 75, 76, 70] #to display the first lambert map
#boxlon, boxlat= [-5, 8, 10, -6.5, -5], [69.5, 69.5, 74.5, 74.5, 69.5]   # lllon, lllat, urlon, urlat= -5, 69.5, 10, 74.5 to display the second lamber map
#boxlon, boxlat= [-3, 12, 15, -3, -3], [63.5, 63.5, 70, 70, 63.5] #lllon, lllat, urlon, urlat= -3, 63.5, 15, 69.5 #to display the third lambert map
boxlon, boxlat= [0, 10, 10, 0, 0], [67.9, 67.9, 74, 74, 68] #to display the area without fluxes


title_extra='' #jet'

save= True
#save=False

"""import the data"""
from netCDF4 import Dataset, num2date
#if 'PressWind' in var: #but includes also SLP


if var=='OLW':
    file= Mediadir + 'ECMWF/HRES/oper_SFC_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'_'+str(shour).zfill(2)+'_OLR.nc'
    nc= Dataset(file)
    print(nc.variables.keys())
    varlist= list(nc.variables.keys())

#    tim= nc.variables['time'][t]
#    datetime= num2date(tim, nc.variables['time'].units)
#    lat= nc.variables[varlist[2]][:]
#    lon= nc.variables[varlist[1]][:]
    OLW= nc.variables['ttr'][:][t-1:t+1]
    OLW= (OLW[1]- OLW[0])/(-3*60**2)

if 'lev' not in var:
    file= Mediadir + 'ECMWF/HRES/oper_SFC_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'_'+str(shour).zfill(2)+'.nc'
    nc= Dataset(file)
    print(nc.variables.keys())
    varlist= list(nc.variables.keys())
    tim= nc.variables['time'][t]
    datetime= num2date(tim, nc.variables['time'].units)
    SLP= nc.variables['msl'][t]/100
    
#    if file== '/media/pst019/1692A00D929FEF8B/ECMWF/HRES/oper_SFC_20080302_12.nc':
#        lat= nc.variables['lat'][:]
#        lon= nc.variables['lon'][:]
#        T= nc.variables['2t'][t]                                 
#        u= nc.variables['10u'][t]
#        v= nc.variables['10v'][t]  
    
#    else:
#    lat= nc.variables[varlist[2]][:]
#    lon= nc.variables[varlist[1]][:]
    lat= nc.variables[varlist[1]][:]
    lon= nc.variables[varlist[0]][:]
    T= nc.variables['t2m'][t]                                 
    u= nc.variables['u10'][t]
    v= nc.variables['v10'][t]
       
    U= np.sqrt(u**2+ v**2) 



elif 'lev' in var:
    file= Mediadir + 'ECMWF/HRES/oper_SFC_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'_'+str(shour).zfill(2)+'.nc'
    nc= Dataset(file)
    SLP= nc.variables['msl'][t]/100    
    
    print('if SLP is missing you have to get it from the SFC file...')
    file= Mediadir + 'ECMWF/HRES/oper_PL_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'_'+str(shour).zfill(2)+'.nc'
    
    nc= Dataset(file)
    print(nc.variables.keys())
    varlist= list(nc.variables.keys())
    
                        
    tim= nc.variables['time'][t]
    datetime= num2date(tim, nc.variables['time'].units)
    
    if 'latitude' in varlist:
        lat= nc.variables['latitude'][:]
        lon= nc.variables['longitude'][:]
        level= int(nc.variables['level'][ilevel])
    
#    if file== '/media/pst019/1692A00D929FEF8B/ECMWF/HRES/oper_PL_20080302_12.nc' or file== '/media/pst019/1692A00D929FEF8B/ECMWF/HRES/oper_PL_20080304_12.nc':
    else:
        lat= nc.variables['lat'][:]
        lon= nc.variables['lon'][:]
        level= int(nc.variables['plev'][ilevel]//100)
        
    
    SH= nc.variables['q'][t, ilevel]
    T= nc.variables['t'][t, ilevel]                                 
    u= nc.variables['u'][t, ilevel]
    v= nc.variables['v'][t, ilevel]
    U= np.sqrt(u**2+ v**2) 
    geop= nc.variables['z'][t, ilevel]/9.81

#    
print(datetime)

grid= np.meshgrid(lon, lat) #pcolormesh needs corner points for projection
Lon, Lat= map(grid[0], grid[1])

if var == 'PressWind_advanced':
    PlotContours(Lon, Lat, SLP, map, leveldist= 1, numbers=False)
    
    PlotWindVelo(Lon, Lat, U, map, Umax= 25, color='YlBu')
    PlotLocalMax(SLP, threshold=1000, distance=60, map= map, lon=lon, lat=lat,
                     data2=U, threshold2=10, distance2= 30) 
#
    PlotLocalMax(U, threshold=18, distance=50, map= map, lon=lon, lat=lat, typ='max',
                         color='orange', dot=False, roundorder= 0, latbound=[68.5, 80])
    
    
if var== 'RH_lev':
    PlotContours(Lon, Lat, SLP, map, leveldist= 1, numbers=False)

    import metpy.calc as mpcalc
    from metpy.units import units    
    dRH= mpcalc.relative_humidity_from_specific_humidity(SH, T* units.K, level* units.hPa)
    
    bounds, colors= PlotColorMap4(Lon, Lat, dRH , map, label= 'Relative humidity at '+str(level)+' hPa', color='blue', bounds= np.array([0, 0.4, 0.7, 0.8, 0.9, 0.95, 1]) )


if var== 'SH_lev':
    PlotContours(Lon, Lat, SLP, map, leveldist= 1, numbers=False)

    bounds, colors= PlotColorMap4(Lon, Lat, SH*1E3 , map, label= 'Specific humidity [g/kg] at '+str(level)+' hPa',
                  color='jet_r', bounds= np.arange(0, 2.3, .2))


    import metpy.calc as mpcalc
    from metpy.units import units    
    dRH= mpcalc.relative_humidity_from_specific_humidity(SH, T* units.K, level* units.hPa)

    PlotContours(Lon, Lat, dRH, map, levels= [0.9], numbers= False, color='white')
    


if var== 'Theta_lev':
    PlotContours(Lon, Lat, SLP, map, leveldist= 1, numbers= False)
#    PlotContours(Lon, Lat, d.geop, map, leveldist= 10)
    if level == 500: bounds= np.arange(278, 288, 1)
    elif level == 700: bounds= np.arange(270, 282, 1) #for pot temp 700
    elif level == 850: bounds= np.arange(264, 280, 1) #for pot temp 850
    elif level <= 950: bounds= np.arange(260, 278, 1) #for pot temp 950 
    else: bounds= None
    
#        bounds= np.arange(267, 290.5, 1)
#        bounds= np.arange(260, 291, 2)
    
    theta= PotTemp(T, level)
    bounds, colors=PlotColorMap4(Lon, Lat, theta, map, bounds= bounds, color='RdBu', label=r"$\theta_{"+str(level)+"}$ [K]")

if var == 'Geop_wind_lev':
    PlotContours(Lon, Lat, geop, map, leveldist= 10, alpha= .7, numbers= False)
#    PlotWindVelo(Lon, Lat, np.sqrt(d.u**2+ d.v**2), map)
    bounds, colors= PlotColorMap4(Lon, Lat, np.sqrt(u**2+ v**2), map, color='YlBu', label='Wind velocity [m/s] at '+str(level)+'hPa', bounds= np.arange(0, 31, 3))


    uzonal, vmeri= u,v #Wind2ZonalMeridional(d.lon, d.lat, d.u, d.v) #this rotates the wind in model coordinate direction to zonal and meridional winds
    PlotWind2(lon, lat, uzonal, vmeri, map, everyx= 3, everyy= 13, rot=False, color='white', quiverkey= False)


if var == 'OLW':
#        levels= np.arange(940, 1040, 3)
    PlotContours(Lon, Lat, SLP, map, leveldist= 3, numbers=False, color= 'r', alpha= .6)
    
    PlotColorMap4(Lon, Lat, OLW , map, label= 'Outgoing longwave radiation [W/m$^2$]', color='grey', bounds= np.arange(130, 241, 10))



if var== 'Vort_lev':
    PlotContours(Lon, Lat, geop, map, leveldist= 10, alpha= .7, numbers= False)
 
    dy= (lat[0]- lat[1])*110E3 #distance along axis 1
    dx= dy* np.cos(np.deg2rad(lat)) #distance along axis 0

    #axis 0: latitude = y-axis  -   axis 1: longitude = x-axis
    vort= np.gradient(v, axis= 1)/dx[:, np.newaxis] - np.gradient(u, -dy, axis= 0) #-dy since the latitudes start from the north  
#    gausfilterdist= 100
    vort= filters.gaussian_filter(vort, sigma= 10, mode='nearest', truncate= 4) #truncate = 1   
    
    PlotColorMap4(Lon, Lat, vort*1E5, map, color='RdBu', symetric=True)

    PlotLocalMax(vort*1E5, threshold=12, distance=200, map= map, lon=lon, lat=lat, typ='max',
                         color='orange', dot=False, roundorder= 0, latbound=[67, 80])


#
#"""ERA wind and surface pressure"""
#if var == 'PressWind':   
#    d.imp_u10()
#    U10= np.sqrt(d.u10[0]**2+d.v10[0]**2)
#
#    PlotWindVelo(Lon, Lat, U10, map, Umax= 25)
#   
#    PlotWind(d.lon, d.lat, d.u10[0], d.v10[0], map)
#    PlotContours(Lon, Lat, d.MSLP[0], map, leveldist= pleveldist)
#
#
#if var == 'PressWind_advanced':
#    PlotContours(Lon, Lat, d.MSLP[0], map, leveldist= 1, numbers=False)
#
#    d.imp_u10()
#    U10= np.sqrt(d.u10[0]**2+d.v10[0]**2)
#    PlotWindVelo(Lon, Lat, U10, map, Umax= 25, color='YlBu')
#
#    PlotLocalMax(d.MSLP[0], threshold=1010, distance=10, map= map, lon=d.lon, lat=d.lat,
#                     data2=U10, threshold2=10, distance2= 7) 
#
#    PlotLocalMax(U10, threshold=15, distance=10, map= map, lon=d.lon, lat=d.lat, typ='max',
#                 color='orange', dot=False, roundorder= 0)        
#
#
#
#
#if var == 'RH850':
#    PlotContours(Lon, Lat, d.MSLP[0], map, leveldist= pleveldist)
#    d.impvar('sHum850')
#    d.impvar('T850')
#    import metpy.calc as mpcalc
#    from metpy.units import units    
#    dRH= mpcalc.relative_humidity_from_specific_humidity(d.sHum850[0], d.T850[0]* units.K, 850* units.hPa)
#    
#    PlotColorMap4(Lon, Lat, dRH , map, 'Relative humidity at 850 hPa', color='blue', bounds= np.array([0, 0.4, 0.7, 0.8, 0.9, 0.95, 1]) )
#
#



"""Plot Thorpex points"""
if Thorpex== True:
    if year == 2008 and month== 3 and plotday== 3 and plothour<15:
        excl= [3, 5, 12] #for the first flight       
    elif year == 2008 and month== 3 and plotday== 3 and plothour>=15:
        excl=[1] #for the second flight       
#        excl=[1,9,11] #for the second flight
    elif year == 2008 and month== 3 and plotday== 4:
        excl= [1, 5, 13]# for the third flight

    
    thor= thorpex.data(year, month, plotday, plothour, plevels= [level ], exclude=excl)
#    thor2= thorpex.data(year, month, day, hour, level= level)
    
    xpt, ypt= map(thor.lon, thor.lat)
    #map.plot(xpt,ypt,'bo')
    
    for i in np.arange(len(xpt)):
        scattercolor= 'black'
        if var=='Geop_wind_lev': scattercolor= 'r'
#            if var is not 'Geop_wind_lev':
        plt.text(xpt[i]-10000, ypt[i]+20000, str(thor.dropnr[i]), color=scattercolor, fontsize=14)               
#                plt.text(xpt[i]+100, ypt[i]+100, str(thor.dropnr[i])+'\n'+ str(thor.datetime[i, 0])[11:16], color='r')
                
        
        if var== 'Theta_lev' or var=='Theta_newbounds_lev':
            value= thor.theta[i]
        elif 'RH_lev' in var:
            value= thor.RH[i]
        elif var== 'SH_lev':
            value= thor.SH[i]*1E3
        elif var== 'Geop_wind_lev':
            value= thor.U[i]                
        else:
            break

        if np.isnan(value[0]): print('Exclude ', thor.dropnr[i])   #in case the value that should be plotted is none it exits the loop
        else:
#                a= value - bounds
#                ind= np.where(a== min(i for i in a if i >= 0))[0][0]+1 #this goes wrong if the thorpex value is outside the bounds
            a= np.abs(value - bounds) # in order to find the color of the point, it calculates the distance of the value to the bounds of the colormap
            ind= np.where(a== np.min(a))[0][0] #finds the index of bounds that should be used as color for the point

            colornow= colors[ind]
    
    #        plt.text(xpt[i]+100, ypt[i]+100, str(np.round(value[0],2) ) )
            scattercolor= 'black'
            if var=='Geop_wind_lev': scattercolor= 'r'
            
            map.scatter(xpt[i],ypt[i], color= colornow, s= 150, edgecolors= scattercolor)
    
    
            if var== 'Geop_wind_lev':
                lon2d, lat2d= np.meshgrid(lon, lat)
                dist= np.sqrt( (thor.lat[i,0]- lat2d)**2+ (np.cos(np.deg2rad(thor.lat[i,0]))* (thor.lon[i,0]- lon2d))**2)
                xpos, ypos= np.where(dist== np.min(dist))
#                    plt.text(xpt[i], ypt[i]+10000, str(int(thor.alt[i,0]- d.geop[xpos[0], ypos[0]])), color='r')  #depature from the geopotential height  
                
                print(thor.dropnr[i], str(int(thor.alt[i,0]- geop[xpos[0], ypos[0]])))        
                u,v, Lon, Lat = map.rotate_vector(thor.u[i], thor.v[i] , thor.lon[i], thor.lat[i] ,returnxy=True)
                plt.quiver(Lon, Lat, u, v, color= 'r', scale = 500, width= .004, headwidth= 5)



if Plotbox== True:
    x,y= map(boxlon, boxlat) #longitude and latitude coordinates of the corner points
    map.plot(x, y, color='b', linewidth= 2)
    title_extra +='_box'


if save== False:
    plt.title('HRES_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+str(shour).zfill(2)+'+'+str(lacktime).zfill(2)+'_'+var+title_extra)

#plt.tight_layout()

if save== True:
    if Thorpex==True:
        savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/DropComp/'
    else:
        savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/Fields/'  
    if '_lev' in var:
        var += str(level)
    savename= savedir + 'HRES_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+str(shour).zfill(2)+'+'+str(lacktime).zfill(2)+'_'+var+title_extra
    print(savename)
    plt.tight_layout()
    plt.savefig(savename, bbox_inches= 'tight')
