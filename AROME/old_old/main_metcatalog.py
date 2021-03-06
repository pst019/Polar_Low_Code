# -*- coding: utf-8 -*-

#!/usr/bin/env python
# Read data from an opendap server


from f_imp_metcatalog import *


"""plot and basemap configurations"""
import numpy as np
import matplotlib.pyplot as plt
fignr= 1
fig = plt.figure(fignr)
plt.clf()


nup= 949
nright = 739 


""" times """
year = 2015
month = 11
day= 19

#year = 2015
#month = 12
#day= 13
#
#year = 2015
#month = 12
#day= 26

#year = 2016
#month = 4
#day= 11

hour= 18 #in 3 hour steps


""" savedirectory"""
savedir= r'/home/'+user+'/home/Polar_Low/Graphs/winter15_16/'
time= str(year)+'_'+str(month).zfill(2)+'_'+str(day).zfill(2)+'_'+str(hour).zfill(2)

#import os
#if not os.path.exists(savedir): os.makedirs(savedir)


from mpl_toolkits.basemap import Basemap

#map= Basemap(projection= 'nplaea',resolution='c',lon_0=0,boundinglat= 50 ,area_thresh=10000, round=1)

#map = Basemap(width=3000000,height=3200000,
#            resolution='l',projection='aea',\
#            lat_1=60.,lat_2=87,lon_0=22,lat_0=76)     

# Lambert Conformal Conic map.
map = Basemap(llcrnrlon=0.,llcrnrlat=60.,urcrnrlon=110.,urcrnrlat=80.,
            projection='lcc',lat_1=67.,lat_2=85.,lon_0=30,
            resolution ='l',area_thresh=2000.)
            
map.drawcoastlines(color='black')
map.drawmeridians(np.arange(0,360,10))#, labels=[0,0,0,1])
map.drawparallels(np.arange(-90,90,10))#, labels=[1,0,0,0])

boxnr= 31
colors= [(plt.cm.seismic(h)) for h in range(256)]
new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=boxnr)





"""a) surface situation """

d= data(year, month, day, hour)
d.imp(usurf= True, vsurf= True, Tsurf= True, psurf= True)

Lon, Lat= map(d.lon, d.lat)

maxlevel = 20
cs= map.pcolormesh(Lon, Lat, d.Tsurf , cmap=new_map, vmin = -maxlevel, vmax =maxlevel, alpha= 1.) 
cb = map.colorbar(cs,location='right',pad='3%')
cb.ax.tick_params(labelsize=16)
cb.set_label('Surface Temperature [$^{\circ}$C]', size= 16)
#plt.savefig(savedir+ time +'_Tsurf', dpi = 50)

""" contours - the pressure"""
cs= map.contour(Lon, Lat, d.psurf, 15, linewidths= 1. , colors= 'k' )#, colors= 6*['b']+ 6*['r'],)
plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f', color= 'black')

""" winds """
u = d.usurf[0::10, 0::10]
v = d.vsurf[0::10, 0::10]
lon= d.lon[0::10, 0::10]
lat= d.lat[0::10, 0::10]

uproj,vproj, Lon, Lat = map.rotate_vector(u,v,lon,lat ,returnxy=True)
#barbs = map.barbs(Lon,Lat,uproj,vproj,length=5,barbcolor='k',flagcolor='r',linewidth=0.5)
Q = map.quiver(Lon,Lat,uproj,vproj,scale=500)
qk = plt.quiverkey(Q, 0.14, 0.07, 20, '20 m/s', labelpos='W')

plt.savefig(savedir+ time +'_Tsurf_psurf_uvsurf', dpi = 50)


"""b) situation at one pressure level"""
fignr += 1
fig = plt.figure(fignr)
plt.clf()

d= data(year, month, day, hour)
plevel= 7
d.imp(uh= plevel, vh= plevel, geoh= plevel, wh= plevel)
print('height', d.press[plevel])

Lon, Lat= map(d.lon, d.lat)
cs= map.contour(Lon, Lat, d.geoh, 15, linewidths= 1. , colors= 'k' )#, colors= 6*['b']+ 6*['r'],)
plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f', color= 'black')

""" vertical velocity"""
maxlevel= 1
cs= map.pcolormesh(Lon, Lat, d.wh , cmap=new_map, vmin = -maxlevel, vmax =maxlevel, alpha= 1.)
cb = map.colorbar(cs,location='right',pad='3%')
cb.ax.tick_params(labelsize=16)
cb.set_label('vertical velocity [m/s]', size= 16)

""" winds """
u = d.uh[0::10, 0::10]
v = d.vh[0::10, 0::10]
lon= d.lon[0::10, 0::10]
lat= d.lat[0::10, 0::10]

uproj,vproj, Lon, Lat = map.rotate_vector(u,v,lon,lat ,returnxy=True)
#barbs = map.barbs(Lon,Lat,uproj,vproj,length=5,barbcolor='k',flagcolor='r',linewidth=0.5)
Q = map.quiver(Lon,Lat,uproj,vproj,scale=500)
qk = plt.quiverkey(Q, 0.14, 0.07, 20, '20 m/s', labelpos='W')



"""c) instability"""
fignr += 1
fig = plt.figure(fignr)
plt.clf()

d= data(year, month, day, hour)
plevel= 6
d.imp(Tsurf= True, Th= 5, PVh= plevel, APE= True)
print('height', d.press[plevel])
                 
Lon, Lat= map(d.lon, d.lat)

""" vertical instability"""
cs= map.contour(Lon, Lat, d.Tsurf-d.Th, levels= [38, 40, 42, 44, 46, 48, 50, 52, 54], linewidths= 1. )#, colors= 'k' )#, colors= 6*['b']+ 6*['r'],)
plt.clabel(cs, fontsize=10, inline=1, fmt='%1.0f', color= 'black')


"""PV """
cs= map.pcolormesh(Lon, Lat, d.PVh* 1E6 , cmap=new_map) #, vmin = -maxlevel, vmax =maxlevel, alpha= 1.) 
cb = map.colorbar(cs,location='right',pad='3%')
cb.ax.tick_params(labelsize=16)
cb.set_label('PV [10$^{-6}$ K m$^2$/ kg s', size= 16)


""" APE """
#cs= map.pcolormesh(Lon, Lat, d.APE , cmap=new_map) #, vmin = -maxlevel, vmax =maxlevel, alpha= 1.) 
#cb = map.colorbar(cs,location='right',pad='3%')
#cb.ax.tick_params(labelsize=16)
#cb.set_label('APE [J/kg]', size= 16)

#plt.savefig(savedir+ time +'_Tsurf', dpi = 50)





#fig.tight_layout()
