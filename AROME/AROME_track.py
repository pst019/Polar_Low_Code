#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

This is from plot_comparison
"""

import time
start = time.time()


import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'



from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import numpy as np
from scipy import stats

# import own modules
import sys  #to import the functions from a different directory
sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_imp_AROME import *  # Read in netcdf file
import f_imp_thorpex as thorpex
from f_meteo import *
from f_useful import *


"""global variables"""
fignr= 13

plt.figure(fignr, figsize= (6, 6))
fignr+=1
plt.clf()
                  
                   
"""b) times"""
year, month = 2008, 3
day, hour= 4, 1                   


"""Lambert coordinates"""
map = AA_map()


fileday, filehour= 2, 0        
t= (day- fileday)*24 + (hour- filehour)  # -1 
 
    
"""name of the experiment"""
#exp_name= '080303_warmctr'
#exp_name='080303_cold_sice'
#exp_name= '080303_warmsens_noQH'
#exp_name='080303_warmsens_2FLX'
#exp_name='080303_warmsens_nocondens'
#exp_name='870226_cold_ERA_can'
#exp_name='080303_cold_pseudo2'
#exp_name='080303_coldsens_noFLX_AREA'
exp_name= 'DA_080303_CTR'              
exp_name= 'DA_080301_cycling'
exp_name= '080302_cold_pseudo'

"""c) plot characteristics"""
pleveldist= 3                   

#var= 'PressWind_advanced'
var= 'Vort_lev'


""" presure level (if the variable ends with '_var' """
pn = 4
#1-950, 4- 850hPa, 6- 700hPa, 8 - 500hPa

save=False

savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/IntensityFieldsTracks/'

title_extra=''

"""end global variables"""



"""start module"""
AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
                                                 
AAres=2 #2- every second datapoint is taken
d= data(filename= AAfilename, res= AAres)

Lon,Lat = map(d.lon,d.lat)
    
    
d.imp_uv_ft(pn= pn)
d.imp_surf(tn= t)

threshold= 7
propdistance= 130E3 #in m   #130  70
distance= int(180E3/(AAres * 2.5E3)) #180  100   #make distance reasonable larger than propdistance to avoid "doublematches" twice is on the safe side



"""calculate vorticity"""
dx= 2500* AAres   
vort_full= (np.gradient(d.v_ft, dx, axis= -1)- np.gradient(d.u_ft, dx, axis= -2))*1E5     

#    vort= (np.gradient(d.v[t], dx, axis= 1)- np.gradient(d.u[t], dx, axis= 0))*1E5         
import scipy.ndimage.filters as filters
gausfilterdist= 100
#    vortfilter= filters.gaussian_filter(vort, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 4.)

vortfilter_full= filters.gaussian_filter1d(vort_full, axis= -1, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 4.)
vortfilter_full= filters.gaussian_filter1d(vortfilter_full, axis= -2, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 4.)
          
#print(data[t][maxima[t]], d.lon[maxima[t]], d.lat[maxima[t]])

bounds, colors= PlotColorMap4(Lon, Lat, vortfilter_full[t], map, label= r"$\xi_{"+str(d.chosen_plevel)+"}$ [10$^{-5}$ 1/s] Gauss filter of std "+str(gausfilterdist)+"km", bounds= np.arange(-27, 27.1, 1), color='RdBu', symetric=True)
PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
title_extra= '_GaussFilter'+str(gausfilterdist)+'km'

maxdata, maxlon, maxlat= PlotLocalMax(vortfilter_full[t], threshold=threshold, distance=distance, map= map, lon=d.lon, lat=d.lat, typ='max', color='g')

    






data= vortfilter_full

latbound= [60, 71]
lonbound= [-5, 10]
timebound= [20, 30]
timebound= [(3- fileday)* 24 + (24- filehour), (3- fileday)* 24 + (34- filehour)] # the time in which the PL has to propage through the box

#data[np.tile(d.lat, (vortfilter_full.shape[0],1,1) ) > latbound[1]] = 0

""" a tracking algorithm"""
vortlist, lonlist, latlist, labellist = EasyTrack(vortfilter_full, d.lon, d.lat, distance, threshold, propdistance, latbound= latbound, lonbound= lonbound)

vortlist, lonlist, latlist, labellist= Tracks_in_lat_lon(vortlist, lonlist, latlist, labellist, latbound, lonbound, timebound)
            
""" change the formate - get one cyclone"""
cyclnratt= labellist[t]
cyclnr= 2



for cyclnr in cyclnratt:
    tcycl, vortcycl, loncycl, latcycl= Data_One_Track(vortlist, lonlist, latlist, labellist, cyclnr)
    
#    plt.figure(fignr-2)
    xpt, ypt= map(loncycl, latcycl)
#    map.plot(xpt,ypt,'bx')
    for i in range(len(loncycl)):
        value= vortcycl[i]
            
#            colors= [(plt.cm.RdBu_r(h)) for h in range(0, 256)]
#            colornow= colors[ int((value-bounds[0])/ (bounds[-1] - bounds[0]) * 256)]

        a= value - bounds
        ind= np.where(a== min(i for i in a if i >= 0))[0][0]+1
        colornow= colors[ind]
        
        
        map.scatter(xpt[i],ypt[i], color= colornow, s= 100, edgecolors= 'r')        
#        plt.text(xpt[i]+100, ypt[i]+100, str(int(vortcyclnr[i])), color='black')
        if i == 0:
            plt.text(xpt[i], ypt[i], str(int(cyclnr)), color='black')


plt.figure(fignr)
fignr+=1
plt.clf()

for cyclnr in remove_dublicate2D(labellist): 
    tcycl, vortcycl, loncycl, latcycl= Data_One_Track(vortlist, lonlist, latlist, labellist, cyclnr)
    
#    plt.figure(fignr-1)
    if len(vortcycl) >0:
        plt.plot(d.datetime[tcycl], vortcycl, label= str(cyclnr))


plt.legend()
plt.ylabel('Vorticity [10$^{-5}$ 1/s]')

plt.gca().xaxis.set_major_locator(mdates.DayLocator())
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d.%m'))

plt.gca().xaxis.set_minor_formatter(mdates.DateFormatter('%H'))
plt.gca().xaxis.set_minor_locator(mdates.HourLocator(byhour=[6,12,18]))


if save== False:
    plt.figure(fignr -2)
    plt.title('Arome '+str(d.datetime[t])[:-3]+ title_extra)

plt.tight_layout()

if save==True:
    if '_lev' in var:
        var += str(d.chosen_plevel)
    plt.figure(fignr -2)
    plt.savefig(savedir+'Arome_'+exp_name+'_'+str(t).zfill(2)+'_'+ var +title_extra +'track')

    plt.figure(fignr -1)
    plt.savefig(savedir+'Arome_'+exp_name+'_'+str(t).zfill(2)+'_'+ var +title_extra +'trackplot')


end = time.time()
print('execution time: ', end - start)


"""find the most important cyclone, in terms of vorticity, that also exists at time tn"""
#for l in labellist[t]:
#    tcycl, datacycl, loncycl, latcycl= Data_One_Track(vortlist, lonlist, latlist, labellist, l)
#    print(l, sum(datacycl))