#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 12:00:39 2017

@author: pst019
With wrf_python a lot of things can be done.
"""

from __future__ import print_function

from netCDF4 import Dataset
from wrf import getvar

#datadir= "/media/'+user+'/1692A00D929FEF8B/ASR/2000/subasr15km.anl.3D.20001231.nc"
datadir= "/media/'+user+'/1692A00D929FEF8B/ASR/asr30km.anl.3D.20000101.nc"
datadir= "/media/'+user+'/1692A00D929FEF8B/PL/WRF/wrfout_d01_2008-03-03_00:00:00"

ncfile = Dataset(datadir)

vort = getvar(ncfile, "avo", timeidx= 1, meta=False)
#
#print(p)