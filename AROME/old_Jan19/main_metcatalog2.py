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


""" import data for a) surface"""
d= data(year, month, day, hour)
d.imp(usurf= True, vsurf= True, Tsurf= True, psurf= True)


"""make the map"""
from f_mapplot import *

ncfile= d.dataset
AA_map()
Lon, Lat = m(d.lon, d.lat)


"""a) surface situation """
from f_plot_fields import *

PlotSurfTemp(Lon, Lat, d.Tsurf, m)
PlotContours(Lon, Lat, d.psurf, m)
PlotWind(lon= d.lon, lat= d.lat, u= d.usurf, v= d.vsurf, m=m)

#plt.savefig(savedir+ time +'_Tsurf_psurf_uvsurf', dpi = 50)


"""b) situation at one pressure level"""
fignr += 1
fig = plt.figure(fignr)
plt.clf()

d= data(year, month, day, hour)

plevel= 7
d.imp(uh= plevel, vh= plevel, geoh= plevel, wh= plevel)
print('height', d.press[plevel])

PlotContours(Lon, Lat, d.geoh, m)

""" vertical velocity"""
PlotColorMap(Lon, Lat, d.wh, m, maxlevel= 1)

""" winds """
PlotWind(lon= d.lon, lat= d.lat, u= d.uh, v= d.vh, m=m)


"""c) instability"""
fignr += 1
fig = plt.figure(fignr)
plt.clf()

d= data(year, month, day, hour)
plevel= 6
d.imp(Tsurf= True, Th= 5, PVh= plevel, APE= True)
print('height', d.press[plevel])
                 

""" vertical instability"""
PlotStaticStability(Lon, Lat, d.Tsurf-d.Th, m)

"""PV """
PlotPV(Lon, Lat, d.PVh, m)


""" APE """
#cs= map.pcolormesh(Lon, Lat, d.APE , cmap=new_map) #, vmin = -maxlevel, vmax =maxlevel, alpha= 1.) 
#cb = map.colorbar(cs,location='right',pad='3%')
#cb.ax.tick_params(labelsize=16)
#cb.set_label('APE [J/kg]', size= 16)
#
##plt.savefig(savedir+ time +'_Tsurf', dpi = 50)
#
#
#
#
#
##fig.tight_layout()
