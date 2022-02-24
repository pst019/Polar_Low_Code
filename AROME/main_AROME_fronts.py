#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

This is from plot_comparison
"""

import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'



from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import numpy as np
import time
#from datetime import date, datetime, timedelta
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
fignr= 2

#maptype='AA'
maptype='AA_half'
#maptype='Lambert'

if maptype=='AA': plt.figure(fignr, figsize= (6, 6))
elif maptype=='AA_half': plt.figure(fignr, figsize= (6, 4.7))
else: plt.figure(fignr, figsize= (6, 4.5))

fignr+=1
 
"""time specification: time of the plot"""
year, month = 2008, 3
day, hour= 3, 18   

"""time of the file"""
fileday, filehour= 3, 0        
t= (day- fileday)*24 + (hour- filehour)  # -1 
#
#year, month = 1987, 2
#day, hour= 26, 12    
#
#fileday, filehour= 25, 12        
#t= (day- fileday)*24 + (hour- filehour)  # -1 

tlist= [t]
#tlist= np.arange(0, 49, 1)                                            

"""Lambert coordinates"""
#lllon, lllat, urlon, urlat= -15, 65, 50, 75
#lat0, lon0= 75, 0 #(lat0, lon0) = center point

#lllon, lllat, urlon, urlat= 0, 67, 60, 75
#lat0, lon0= 70, 0 #(lat0, lon0) = center point

#to match with Thorpex                    
if year == 2008 and month== 3 and day== 3 and hour<15:
    lllon, lllat, urlon, urlat= -4, 70, 19, 75
    lat0, lon0= 75, 0 #(lat0, lon0) = center point

elif year == 2008 and month== 3 and day== 3 and hour>=15:
    lllon, lllat, urlon, urlat= -6, 69.5, 13, 74.5
    lat0, lon0= 75, 0 #(lat0, lon0) = center point

elif year == 2008 and month== 3 and day== 4:
    lllon, lllat, urlon, urlat= -3, 63.5, 15, 69.5
    lat0, lon0= 70, 0 #(lat0, lon0) = center point    
    


    
"""name of the experiment"""
exp_name= '080303_warmctr'
#exp_name='080303_cold_ERA'
#exp_name= '080303_warmsens_noQH'
#exp_name='080303_warmsens_noFLX'
#exp_name='080303_warmsens_nocondens'
#exp_name='870226_cold_ERA_can'
#exp_name='080303_cold_pseudo2'
#exp_name='080304_cold_pseudo'
#exp_name='080303_coldsens_noFLX_AREA'
#exp_name='08030312_cycling'

  
"""c) plot characteristics"""
pleveldist= 1                   


var= 'Front_lev'
#bounds= np.arange(264, 286, 1) #for pottempe850


""" presure level (if the variable ends with '_var' """
pn = 4
#1-950, 4- 850hPa, 6- 700hPa, 8 - 500hPa



        
save=False

#savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/PseudoSat_movie/'
if var in ['CTT', 'WVTC', 'WVT', 'CWR']:
    savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/PseudoSat2/'
elif Thorpex==True:
    savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/DropComp/'

else:
    savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/Fields/'
#    savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/IntensityFields/'
                 
                 
"""end global variables"""




for t in tlist:
    print(t)
    plt.clf()
    title_extra=''

    if maptype== 'AA': map = AA_map()
    elif maptype== 'AA_half': map = AA_map_half()
    else: map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)    
    
    
    AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
                                                     
    AAres=1 #2- every second datapoint is taken
    d= data(filename= AAfilename, res= AAres)
    
    Lon,Lat = map(d.lon,d.lat)
        
        
    d.imp_surf(tn= t)
    
    if '_lev' in var:
        d.imp_level(pn= pn, tn= t)
 



    if var== 'Front_lev':
        gausfilterdist = 20 #in km
        gausfilterdist_th= 20
        
               
        spechum= d.RH/(.263 *d.pres[pn])* np.exp(17.67*(d.T-273.15) / (d.T-29.65))
        theta_e= EquiPotTemp(d.T, spechum, d.chosen_plevel)
        theta_e_f= filters.gaussian_filter(theta_e, sigma= gausfilterdist_th/(AAres*2.5), truncate= 4.)

        
#        PlotColorMap3(Lon, Lat, theta_e, map, bounds= bounds, label=r"$\theta_{e,"+str(d.chosen_plevel)+"}$ [K]")
        PlotContours(Lon, Lat, theta_e_f, map, leveldist= pleveldist)
#        PlotContours(Lon, Lat, theta_e, map, leveldist= pleveldist, color='blue')
        
        
#        plt.figure(fignr, figsize= (6, 4.7), clear=True)
#        map = AA_map_half()
#        fignr +=1
        
        dx= 2500* AAres
        grad_theta_e= np.gradient(theta_e_f, dx, edge_order=1)
        abs_grad_theta_e= np.sqrt(grad_theta_e[0]**2+ grad_theta_e[1]**2)
        abs_grad_theta_e_filter= filters.gaussian_filter(abs_grad_theta_e, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 4.)

        PlotColorMap4(Lon, Lat, abs_grad_theta_e_filter *1E3*1E2, map, color= 'red', bounds=np.arange(0,10.1,1), label=r"$\nabla \theta_{e,"+str(d.chosen_plevel)+"}$ [K/100km]")
        
        TFP= - np.sum(np.array(np.gradient(abs_grad_theta_e, dx, edge_order=1)) * grad_theta_e/abs_grad_theta_e  , axis= 0)
        TFP_filter= filters.gaussian_filter(TFP, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 4.)

        plt.figure(fignr, figsize= (6, 4.7), clear=True)
        map = AA_map_half()
        fignr +=1


        PlotColorMap4(Lon, Lat, TFP_filter, map, color= 'RdBu', symetric=True, label=r"TFP$_{"+str(d.chosen_plevel)+"}$")

        front= np.zeros(np.shape(theta_e))
        front[np.logical_and(np.abs(TFP_filter)< .1* 1E-9, abs_grad_theta_e_filter >3*1E-5)]= 1

        plt.figure(fignr, figsize= (6, 4.7), clear=True)
        map = AA_map_half()
        fignr +=1
             
        PlotColorMap4(Lon, Lat, front, map, color= 'red')#, label=r"$\nabla \theta_{e,"+str(d.chosen_plevel)+"}$ [K/100km]")
             

        

    
    if save== False:
        plt.title('Arome '+str(d.datetime[t])[:-3]+ title_extra)
    
    plt.tight_layout()
    
    if save==True:
        if '_lev' in var:
            var += str(d.chosen_plevel)
        plt.savefig(savedir+'Arome_'+exp_name+'_'+str(t).zfill(2)+'_'+ var +title_extra)
        print(savedir+'Arome_'+exp_name+'_'+str(t).zfill(2)+'_'+ var +title_extra)



