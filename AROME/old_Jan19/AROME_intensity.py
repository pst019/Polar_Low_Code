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
Mediadir= '/media/'+user+'/1692A00D929FEF8B/'



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
import matplotlib.dates as mdates

#import itertools

marker = ['-', '--', ':', '-.', '>', 's', '8', 'p']

"""global variables"""
fignr= 7

plt.figure(fignr)
fignr+=1
plt.clf()
                  
                   
"""b) times"""
year, month = 2008, 3
day, hour= 3, 3                   


"""Lambert coordinates"""



fileday, filehour= 3, 0        
t= (day- fileday)*24 + (hour- filehour)  # -1 
 
    
"""name of the experiment"""
#exp_name= '080303_warmctr'
exp_name='080303_cold_sice'
#exp_name= '080303_warmsens_noQH'
#exp_name='080303_warmsens_noFLX'
#exp_name='080303_warmsens_nocondens'
#exp_name='870226_cold_ERA_can'
#exp_name='080303_cold_pseudo2'
                   
"""c) plot characteristics"""

#var= 'PressWind_advanced'
var= 'Vort_lev'


""" presure level (if the variable ends with '_var' """
pn = 4
#1-950, 4- 850hPa, 6- 700hPa, 8 - 500hPa

save=True
savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/IntensityFieldsTracks/'

title_extra=''

"""end global variables"""

threshold= 10
propdistance= 130E3 #in m   #70
distance= int(180E3/(AAres * 2.5E3)) #100   #make distance reasonable larger than propdistance to avoid "doublematches" twice is on the safe side


latbound= [60, 71]
lonbound= [-5, 10]
timebound= [20, 30]

for exp_name in ['080303_warmctr', '080303_warmsens_noFLX', '080303_coldsens_noFLX_AREA', '080303_warmsens_2FLX', '080303_warmsens_noQH',
                 '080303_warmsens_noTH', '080303_warmsens_nocondens', '080303_cold_pseudo2',  '080303_cold_sice',
                 '080303_cold_ERA']:

    """import data"""
    AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
                                                     
    AAres=2 #2- every second datapoint is taken
    d= data(filename= AAfilename, res= AAres)     
        
    d.imp_uv_ft(pn= pn)
         
    """calculate vorticity"""
    if var== 'Vort_lev':
        dx= 2500* AAres
        
        vort_full= (np.gradient(d.v_ft, dx, axis= -1)- np.gradient(d.u_ft, dx, axis= -2))*1E5     
        import scipy.ndimage.filters as filters
        gausfilterdist= 60
        
        vortfilter_full= filters.gaussian_filter1d(vort_full, axis= -1, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 4.)
        vortfilter_full= filters.gaussian_filter1d(vortfilter_full, axis= -2, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 4.)
     
        title_extra= '_GaussFilter'+str(gausfilterdist)+'km'
        
             
    
    
    """ a tracking algorithm"""
    vortlist, lonlist, latlist, labellist = EasyTrack(vortfilter_full, d.lon, d.lat, distance, threshold, propdistance)
    vortlist, lonlist, latlist, labellist= Tracks_in_lat_lon(vortlist, lonlist, latlist, labellist, latbound, lonbound, timebound)

    
    
    """ change the formate - get one cyclone"""
    cyclnratt= labellist[t]
    cyclnr= 2
    
    color = next(plt.gca()._get_lines.prop_cycler)['color']
    
    for i, cyclnr in enumerate(remove_dublicate2D(labellist)): # [2]: #cyclnratt:
        tcycl, vortcycl, loncycl, latcycl= Data_One_Track(vortlist, lonlist, latlist, labellist, cyclnr)
        
        if i== 0:
            plt.plot(d.datetime[tcycl], vortcycl, color= color, linestyle=marker[i], label= exp_name[7:]) # +'_'+ str(cyclnr))
        else:
            plt.plot(d.datetime[tcycl], vortcycl, color= color, linestyle=marker[i])
    
plt.legend()

plt.ylabel('Vorticity [10$^{-5}$ 1/s]')

plt.gca().xaxis.set_major_locator(mdates.DayLocator())
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d.%m'))

plt.gca().xaxis.set_minor_formatter(mdates.DateFormatter('%H'))
plt.gca().xaxis.set_minor_locator(mdates.HourLocator(byhour=[6,12,18]))



#if save== False:
#    plt.title('Arome '+str(d.datetime[t])[:-3]+ title_extra)

plt.tight_layout()

if save==True:
    if '_lev' in var:
        var += str(d.chosen_plevel)
    plt.savefig(savedir+'Arome_allexp_trackplot')



end = time.time()
print('execution time: ', end - start)


"""find the most important cyclone, in terms of vorticity, that also exists at time tn"""
#for l in labellist[t]:
#    tcycl, datacycl, loncycl, latcycl= Data_One_Track(maxvortlist, lonlist, latlist, labellist, l)
#    print(l, sum(datacycl))