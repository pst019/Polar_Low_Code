#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

This is for importing HRES data. It is also a test to see how gribdata actually looks like
"""

from netCDF4 import Dataset
import numpy as np
import time

import os
user = os.getcwd().split('/')[2]

from f_useful import *

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'

filedir= Mediadir+'ECMWF/HRES/oper_sfc20180223_0000_av.nc'
#year
#month
#day
#tstart

#        print(nc.variables.keys())
#this one is even more interesting:
#        print(nc.variables.values())
#"""should remove height above msl"""

#nc= Dataset(filedir)
#tim = nc.variables['time'][:]
#
#T2= nc.variables['var167_2'][:]


filedir= Mediadir+'ECMWF/HRES/oper_sfc20180223_0000_av.grib'
import pygrib

grbs= pygrib.open(filedir)

# some information on https://jswhit.github.io/pygrib/docs/

#some commands
#grbs.select()
#grbdata= grbs.message(2)
#
#grbdata.latlons()[0][150:180,0]