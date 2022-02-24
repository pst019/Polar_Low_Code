#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  5 14:52:50 2018

@author: pst019

conda install -c conda-forge metpy

https://github.com/Unidata/MetPy/tree/master/staticdata
"""

import os
user = os.getcwd().split('/')[2]

import sys
sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_mapplot import * #make basemap plots
from f_imp_thorpex import data as Tdata #import the thorpex data

from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from metpy.cbook import get_test_data
import metpy.calc as mpcalc
from metpy.plots import SkewT, Hodograph
from metpy.units import units, concatenate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from metpy.calc.thermo import *

Mediadir= '/media/'+user+'/1692A00D929FEF8B/'
savedir= '/home/'+user+'/home/Polar_Low/Graphs/2008_03_03/Soundings/'

year, month = 2008, 3
day, hour= 3, 12  


dropid=1
#exclude= 3, 5, 6,7( above 500hPa)
for dropid in [x for x in range(8, 9) if x not in [3, 5]]:  #for the first flight
#for dropid in [x for x in range(1, 10) if x not in [1, 9]]: #for the second flight
#for dropid in [x for x in range(19, 20) if x not in [5, 13]]: #for the third flight

    
    fignr= 8


    
    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dictated by the typical meteorological plot
    
    
    
    
    thor= Tdata(year, month, day, hour, level='one_profile', dropid= dropid)
    
    
    thorV= Tdata(year, month, day, hour, level='vertical', plevels= np.array([1000, 950, 925, 900, 850, 800, 700, 600, 500, 400]))
    idrop= int(np.where(thorV.dropnr == dropid)[0])
    
    p= thorV.pres[idrop]* units.hPa
    T= thorV.T[idrop]* units.degC
    Td= dewpoint_rh(T, thorV.RH[idrop])
    
    
    """get the AA data and also plot it"""
#    hour= thor.datetime[0].hour #the hour where the dropsonde reached the ground
#    
#    from f_imp_AROME import data as Adata
##    exp_name= '080303_warmctr'
#    exp_name= '080303_cold_pseudo2'
##    exp_name= '080303_warmsens_noTH'
#    fileday, filehour= 3, 0        
#    t= (day- fileday)*24 + (hour- filehour)  # -1
#    t= 12
#    AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
#    AA= Adata(filename= AAfilename, res=1)
#    
#    dist= (AA.lat - thor.lat[0])**2 + (np.cos(np.deg2rad(thor.lat[0])) * (AA.lon - thor.lon[0])**2)
#    x,y= np.where(dist == np.min(dist))
#    
#    AA.imp_cross_sec_full(xn= x[0], yn= y[0], tn= t)
    
    import metpy.calc as mpcalc
    from metpy.units import units
    AARH_fromSH= mpcalc.relative_humidity_from_specific_humidity(AA.SH, AA.T* units.K, AA.pres* units.hPa)

#    thorVRH_fromSH= mpcalc.relative_humidity_from_specific_humidity(thorV.SH[idrop], T, p)
    
   
    """profiles"""
    plt.figure(fignr)
    fignr +=1
    plt.clf()
    
   
    plt.plot(thor.RH, thor.pres, 'g', label='Thor')
    plt.plot(thorV.RH[idrop], thorV.pres[idrop], 'go')
#    plt.plot(thorVRH_fromSH, thorV.pres[idrop], 'ro')
    
    
    plt.plot(AA.RH, AA.pres, 'o', label='AA')
    plt.plot(AARH_fromSH, AA.pres, 'o', label='AA from SH')

    plt.xlabel('Relative humidity [%]')
    plt.ylabel('Altitude [hPa]')
    plt.ylim([1000, 400])
#    plt.xlim([0.4,1.05])
    plt.title('Sounding Dropsonde '+str(dropid))

    plt.legend()
    