#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 16:38:42 2017

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'


import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})

from scipy import stats

import sys

#this is the path where the functions are saved
sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions/')
from f_imp_ERA2 import *
from f_impLists import *
from f_meteo import *



import time
start_time = time.time()


vortlim= 5
latcircle= 70
circlewidth= 2

syear= 1979
eyear= 2017

CYlist= np.zeros((7, 0))

for year in np.arange(syear, eyear):
    for month in [1, 2, 12]:
       
        T=TRACK_list_month2(year, month)
        T.local(minlon= -180, maxlon= 180, minlat= latcircle- circlewidth, maxlat= latcircle+ circlewidth) #has to be within ERA data
#        CYlist #nr, tCY, lon, lat, vort
        
        T.PLlist= T.PLlist[:, T.PLlist[4] > vortlim]

        T.PLlist= np.vstack((T.PLlist, year* np.ones((1, T.PLlist.shape[1])), month* np.ones((1, T.PLlist.shape[1]))))
        CYlist= np.hstack((CYlist, T.PLlist))
    


nrCYs, percentile97 , percentile90, percentile70, yearvec= [], [], [], [], []

for year in np.arange(syear, eyear):
    for i, month in enumerate([1, 2, 12]):


        CYlistmonths=  CYlist[:,np.logical_and(CYlist[5] == year, CYlist[6] == month)]
        
        nrCYs += [len(CYlistmonths[4])]
        percentile97 += [np.percentile(CYlistmonths[4], 97.5)]
        percentile90 += [np.percentile(CYlistmonths[4], 90)]
        percentile70 += [np.percentile(CYlistmonths[4], 75)]
        
        yearvec += [year+ i/3]



plt.figure(1, figsize= (8, 4))
plt.clf()
plt.plot(yearvec, nrCYs, 'bx')
plt.ylabel('Number of 6-hourly cyclone timesteps')# in latitude band')
plt.xlim([syear, eyear])


slope, intercept, r_value, p_value, std_err = stats.linregress(yearvec, nrCYs)
plt.plot(yearvec,  intercept +np.array(yearvec)*slope, label= 'trend: '+str(np.round(slope*10 ,3))+ ' cycl/decade \np-value: '+str(np.round(p_value, 3)), color='b')
plt.legend()
plt.tight_layout()

plt.savefig(homedir + 'Polar_Low/ERA_cycl/Number_vortlim'+str(vortlim)+'_lat'+str(latcircle)+'_width'+str(circlewidth))


plt.figure(2, figsize= (8, 4))
plt.clf()
plt.plot(yearvec, percentile97, 'bx')
plt.plot(yearvec, percentile90, 'gx')
plt.plot(yearvec, percentile70, 'rx')

plt.ylabel('Percentile in Vorticity [10$^{-5}$1/s]')
plt.xlim([syear, eyear])

slope, intercept, r_value, p_value, std_err = stats.linregress(yearvec, percentile97)
plt.plot(yearvec,  intercept +np.array(yearvec)*slope, label= 'trend: '+str(np.round(slope*10 ,3))+ ' vort/decade \np-value: '+str(np.round(p_value, 3)), color= 'b')
plt.text(1980, 11.5, '97.5 percentile', color='b')


slope, intercept, r_value, p_value, std_err = stats.linregress(yearvec, percentile90)
plt.plot(yearvec,  intercept +np.array(yearvec)*slope, label= 'trend: '+str(np.round(slope*10 ,3))+ ' vort/decade \np-value: '+str(np.round(p_value, 3)), color= 'g')
plt.text(1980, 9, '90 percentile', color='g')


slope, intercept, r_value, p_value, std_err = stats.linregress(yearvec, percentile70)
plt.plot(yearvec,  intercept +np.array(yearvec)*slope, label= 'trend: '+str(np.round(slope*10 ,3))+ ' vort/decade \np-value: '+str(np.round(p_value, 3)), color= 'r')
plt.text(1980, 6, '70 percentile', color='r')

plt.legend()
plt.tight_layout()
plt.savefig(homedir + 'Polar_Low/ERA_cycl/Percentile_vortlim'+str(vortlim)+'_lat'+str(latcircle)+'_width'+str(circlewidth))



np.savetxt(homedir + 'Polar_Low/ERA_cycl/Cyclone.txt', CYlist)



print("--- %s seconds ---" % (time.time() - start_time))