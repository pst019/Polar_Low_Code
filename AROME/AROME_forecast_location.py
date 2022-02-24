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

import scipy.ndimage.filters as filters

"""global variables"""
fignr= 6

plt.figure(fignr, figsize= (4.5, 4.5))
fignr+=1
plt.clf()
                  
                   
"""b) times"""
year, month = 2008, 3
day, hour= 4, 4                   


"""Lambert coordinates"""
#map = AA_map_half()
map = Lambert_map(lllon= -5, lllat= 68, urlon=30, urlat= 70, lat0= 77.5, lon0= -25, res='i', fill=False, coastline=True, latdist=5, londist=10)


#t= (day- fileday)*24 + (hour- filehour)  # -1 
 
    
"""name of the experiment"""
exp_name= 'DA_080301_cycling'
fileday, filehour= 3, 18        

#
#exp_name= 'DA_080303_CTR'              
#fileday, filehour= 3, 0        

     
"""c) plot characteristics"""
pleveldist= 3                   

#var= 'PressWind_advanced'
var= 'Vort_lev'


""" presure level (if the variable ends with '_var' """
pn = 4
#1-950, 4- 850hPa, 6- 700hPa, 8 - 500hPa

save=False
#save= True

savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/DifferentTracks/'

title_extra=''

threshold= 6
propdistance= 130E3 #in m   #130  70
AAres=2 #2- every second datapoint is taken

distance= int(180E3/(AAres * 2.5E3)) #180  100   #make distance reasonable larger than propdistance to avoid "doublematches" twice is on the safe side

latbound= [60, 69] #60, 69
lonbound= [-5, 9] #-5, 9


gausfilterdist= 100


"""end global variables"""
exptype=''

#exp_list= ['DA_080303_p2SST']
#exp_list=['DA_080301_cycling']
#exp_list=['DA_08030212']
#exp_list=['080302_cold_pseudo']
#fileday_list= [2]
#filehour_list= [0]


#exp_list= ['080302_cold_pseudo', 'DA_08030212', 'DA_080303_CTR', 'DA_080301_cycling', 'DA_080301_cycling', 'DA_080301_cycling']
#fileday_list= [2, 2, 3, 3, 4, 4]
#filehour_list= [0, 12, 0, 12, 0, 6]
#
#exp_list= ['DA_08030212', 'DA_080303_CTR', 'DA_080301_cycling', 'DA_080301_cycling', 'DA_080301_cycling', 'DA_080301_cycling']
#fileday_list= [2, 3, 3, 4, 4, 4]
#filehour_list= [12, 0, 12, 0, 6, 12]


#exptype='_sens_SST'
#exp_list= ['DA_080303_p6SST','DA_080303_p4SST','DA_080303_p2SST', 'DA_080303_CTR', 'DA_080303_m2SST','DA_080303_m4SST','DA_080303_m6SST']
#legendlabel=[ '+6SST', '+4SST', '+2SST', 'CTR', '-2SST', '-4SST', '-6SST']
#fileday_list= [3]*len(exp_list)
#filehour_list= [0]*len(exp_list)

exptype='_sens_pSST'
exp_list= ['DA_080303_p6SST','DA_080303_p4SST','DA_080303_p2SST', 'DA_080303_CTR']
legendlabel=[ '+6SST', '+4SST', '+2SST', 'CTR']
fileday_list= [3]*len(exp_list)
filehour_list= [0]*len(exp_list)


#exptype='_sens'
#exp_list= ['DA_080303_CTR', 'DA_080303_noQH', 'DA_080303_noCondens', 'DA_080303_noTH', 'DA_080303_noFLX', 'DA_080303_noFLXarea', 'DA_080303_2FLX']
#legendlabel= ['CTR', 'noQH', 'noCond', 'noTH', 'noFLX', 'noFLX-A', '2FLX']
#fileday_list= [3]*len(exp_list)
#filehour_list= [0]*len(exp_list)


point_day, point_hour= 4, 6
point_day2, point_hour2= 3, 12

"""start module"""
for ie, exp_name in enumerate(exp_list):
    color = next(plt.gca()._get_lines.prop_cycler)['color']
    print(exp_name)
    
    timebound= [(3- fileday_list[ie])* 24 + (20- filehour_list[ie]), (3- fileday_list[ie])* 24 + (37- filehour_list[ie])] # the time in which the PL has to propage through the box

    AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday_list[ie]).zfill(2)+str(filehour_list[ie]).zfill(2)+'_fp_extract.nc'
                                                     
    d= data(filename= AAfilename, res= AAres)
    
    Lon,Lat = map(d.lon,d.lat)
        
        
    d.imp_uv_ft(pn= pn)
    #d.imp_surf(tn= t)
    
    
    
    
    """calculate vorticity"""
    dx= 2500* AAres   
    vort_full= (np.gradient(d.v_ft, dx, axis= -1)- np.gradient(d.u_ft, dx, axis= -2))*1E5     
    
    vortfilter_full= filters.gaussian_filter1d(vort_full, axis= -1, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 4.)
    vortfilter_full= filters.gaussian_filter1d(vortfilter_full, axis= -2, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 4.)
              
       
    
    """ a tracking algorithm"""
    vortlist, lonlist, latlist, labellist = EasyTrack(vortfilter_full, d.lon, d.lat, distance, threshold, propdistance, latbound= latbound, lonbound= lonbound)
    vortlist, lonlist, latlist, labellist= Tracks_in_lat_lon(vortlist, lonlist, latlist, labellist, latbound, lonbound, timebound)
                
    """ change the formate - get one cyclone"""
#    t= np.max([timebound[0], 0])
#    cyclnratt= labellist[t]
#    cyclnratt= labellist[tpoint]
    cyclnratt= remove_dublicate2D(list(labellist))

    print(cyclnratt)
    for ci, cyclnr in enumerate(cyclnratt):
        tcycl, vortcycl, loncycl, latcycl= Data_One_Track(vortlist, lonlist, latlist, labellist, cyclnr)
        
        xpt, ypt= map(loncycl, latcycl)
#        map.plot(xpt[0],ypt[0],'o', color=color, label= 'SIM-'+str(fileday_list[ie]).zfill(2)+'-'+str(filehour_list[ie] ).zfill(2) )

        tpoint= (point_day- fileday_list[ie])* 24 + (point_hour - filehour_list[ie])
        tpoint2= (point_day2- fileday_list[ie])* 24 + (point_hour2 - filehour_list[ie])

        if ci== 0:
            if exptype !='':
                map.plot(xpt[tcycl==tpoint],ypt[tcycl==tpoint],'s', color=color, label= legendlabel[ie], zorder= 5, markersize=7)
            else:
                map.plot(xpt[tcycl==tpoint],ypt[tcycl==tpoint],'s', color=color, label= 'SIM-'+str(fileday_list[ie]).zfill(2)+'-'+str(filehour_list[ie] ).zfill(2), zorder= 5, markersize=7)
        else: map.plot(xpt[tcycl==tpoint],ypt[tcycl==tpoint],'s', color=color, zorder= 5, markersize=7)

        if tpoint2 >= 0: map.plot(xpt[tcycl==tpoint2],ypt[tcycl==tpoint2],'o', color=color, zorder= 5, markersize=8)

        #plot the location of analysis time steps
        map.plot(xpt[tcycl==0],ypt[tcycl==0],'x', color='k', zorder= 5)
        
#        map.plot(xpt[tcycl==tpoint],ypt[tcycl==tpoint],'s', color=color, label= '+ '+str(tpoint).zfill(2), zorder= 5)

        map.plot(xpt,ypt, color=color)



plt.legend(loc=1)

plt.tight_layout()

    
if save==True:
    print(savedir+'Arome_tracks_analysis_mark'+exptype)
    plt.savefig(savedir+'Arome_tracks_analysis_mark'+exptype)

