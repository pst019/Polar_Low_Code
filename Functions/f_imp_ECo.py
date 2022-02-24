#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

Import data from the EC operational
"""

from netCDF4 import Dataset
import numpy as np
import time

import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'

class dataECo:
    def __init__(self, fieldname='multiplefields', timeperiod='2008.03.02-04', filedir=Mediadir +'/ERA/ECoperational/',
                 tstart=0, tend=None):
        """import data from the ECoperational
        tstart and tend gives the start and the end point in the 6hourly timsteps starting 6 hours after the start date
        """

        self.tstart= tstart
        

        nc= Dataset(filedir+fieldname+'.'+timeperiod+'.nc')
        
#        variables= nc.variables.keys()
        #this one is even more interesting:
#        print(nc.variables.values())
        
        """reduce resolution"""
        if tend ==None:
            self.tim = nc.variables['time'][tstart:]
            self.tend= len(self.tim)
        else:
            self.tim = nc.variables['time'][tstart:tend]
            self.tend=tend
#            
        self.lat = nc.variables['lat'][:]
        self.lon = nc.variables['lon'][:]
        self.level= nc.variables['lev'][:]/100 #hPa
        



        self.T = nc.variables['var130'][self.tstart:self.tend, :, :, :] # (time, lat, lon)
        self.sHum= nc.variables['var133'][self.tstart:self.tend,:, :, :]   #kg H2O / kg air
        self.vort= nc.variables['var138'][self.tstart:self.tend, :, :, :]*1E5 
        self.uPL= nc.variables['var131'][self.tstart:self.tend,:, :, :]
        self.vPL= nc.variables['var132'][self.tstart:self.tend,:, :, :]
        self.Geop= nc.variables['var129'][self.tstart:self.tend,:, :, :]
