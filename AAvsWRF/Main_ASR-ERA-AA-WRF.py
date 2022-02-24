#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 12:27:21 2017

@author: pst019
"""

import numpy as np

import os
user = os.getcwd().split('/')[2]
homedir= '/home/'+user+'/home/'



"""global variables"""
fignr= 1
savedir='Polar_Low/AromeArctic-vs-WRF/Graphs/2008_03_03/'
savefigs= False
                   
"""b) time"""
year=2008
month=3           
day= 4
hour= 00

"""a) Lambert coordinates"""
lllon, lllat, urlon, urlat= -15, 65, 50, 75
lat0, lon0= 75, 0 #(lat0, lon0) = center point

lllon, lllat, urlon, urlat= -5, 67, 28, 76
lat0, lon0= 70, 0 #(lat0, lon0) = center point
                   
#if day == 3:                                      
#    lllon, lllat, urlon, urlat= -5, 67, 28, 76
#    lat0, lon0= 70, 0 #(lat0, lon0) = center point
#
#elif day == 4:          
#    lllon, lllat, urlon, urlat= 0, 65, 18, 70
#    lat0, lon0= 75, 0 #(lat0, lon0) = center point

"""c) plot characteristics"""
PlotPressWind= True
#number_plevels= 20
pleveldist= 1                   

PlotTheta700= False
pot_temp_bounds= np.arange(264, 290, 1)

PlotTheta_e850= False

PlotVorticity=False

"""end global variables"""


"""execute ASR and ERA"""
exec(open("main_ASR-ERA.py").read())
          
#"""execute ECoperational"""
#exec(open("main_ECo.py").read())
               
"""execute WRF"""
exec(open(homedir+"Polar_Low/Code2/WRF/main_WRF.py").read())

"""execute AA"""
exec(open(homedir+"Polar_Low/Code2/AROME/main_AROME.py").read())

"""execute Thorpex"""
exec(open(homedir+"Polar_Low/Code2/Thorpex/main_thorpex.py").read())