#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 15:06:41 2017

@author: pst019
from test_sphericalharmonics2.py
"""


import os
user = os.getcwd().split('/')[2]


import numpy as np
import matplotlib.pyplot as plt

# fast spherical harmonic lib from https://bitbucket.org/nschaeff/shtns
import shtns

Dropboxdir= '/home/'+user+'/Dropbox/'

import sys
sys.path.insert(0, '/home/'+user+'/Dropbox/Polar_Low/polar_low_code/Functions')
from f_imp_ERA2 import *
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
#from f_useful import *
from f_plot_on_map import *


fignr= 1

class Spharmt(object):
    """
    wrapper class for commonly used spectral transform operations in
    atmospheric models
    Jeffrey S. Whitaker <jeffrey.s.whitaker@noaa.gov>
    """
    def __init__(self,nlons,nlats,ntrunc,gridtype='gaussian'):
        """initialize
        nlons:  number of longitudes
        nlats:  number of latitudes
        ntrunc: spectral truncation
        gridtype: 'gaussian' (default) or 'regular'"""
        self._shtns = shtns.sht(ntrunc, ntrunc, 1,\
                shtns.sht_fourpi|shtns.SHT_NO_CS_PHASE)
        if gridtype == 'gaussian':
            self._shtns.set_grid(nlats,nlons,shtns.sht_quick_init|shtns.SHT_PHI_CONTIGUOUS,1.e-8)
#        elif gridtype == 'regular':
#            self._shtns.set_grid(nlats,nlons,shtns.sht_reg_dct|shtns.SHT_PHI_CONTIGUOUS,1.e-8)
#        self._shtns.print_info()
        self.nlons = nlons
        self.nlats = nlats
        self.nlm = self._shtns.nlm

    def grdtospec(self,data):
        """compute spectral coefficients from gridded data"""
        data = np.ascontiguousarray(data, dtype=np.float)
        if data.ndim == 2:
            dataspec = np.empty(self.nlm, dtype=np.complex)
            self._shtns.spat_to_SH(data, dataspec)
        elif data.ndim == 3:
            dataspec = np.empty((data.shape[0],self.nlm), dtype=np.complex)
            for k,d in enumerate(data):
                self._shtns.spat_to_SH(d, dataspec[k])
        else:
            raise IndexError('data must be 2d or 3d')
        return dataspec
    def spectogrd(self,dataspec):
        """compute gridded data from spectral coefficients"""
        dataspec = np.ascontiguousarray(dataspec, dtype=np.complex)
        if dataspec.ndim == 1:
            data = np.empty((self.nlats,self.nlons), dtype=np.float)
            self._shtns.SH_to_spat(dataspec, data)
        elif dataspec.ndim == 2:
            data = np.empty((dataspec.shape[0],self.nlats,self.nlons), dtype=np.float)
            for k,d in enumerate(dataspec):
                self._shtns.SH_to_spat(d, data[k])
        else:
            raise IndexError('dataspec must be 1d or 2d')
        return data





year= 2008
month, day, hour= 3, 4, 12
t= ((day-1)*24+hour)//6 #hours since beginning of month

"""import data"""
d= data('Vort', year, month)


nlons = 720#len(d.lon)  # number of longitudes
nlats = nlons/2  #len(d.lat)   # for gaussian grid.
ntim= 124

ntrunc = 100 #int(nlons/10) #int(nlons/2)-1  # spectral truncation (for alias-free computations)
ntrunc2= 40

vort = np.zeros((ntim, nlats, nlons))
vort[:, 10:120, :]= d.vort    
 
x = Spharmt(nlons,nlats,ntrunc) #,gridtype='regular')   
vortspec= x.grdtospec(vort)
vortback= x.spectogrd(vortspec)[t]


x = Spharmt(nlons,nlats,ntrunc2) #,gridtype='regular')
vort2spec= x.grdtospec(vort)
vort2back= x.spectogrd(vort2spec)[t]

boxnr= 17
colors= [(plt.cm.RdBu_r(h)) for h in range(256)]
new_map = plt.matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=boxnr)



#plt.figure(fignr)
#fignr += 1
#plt.clf()
#map= Polar_map(latsouth= d.lat[-1])
#grid= np.meshgrid(d.lon, d.lat) #pcolormesh needs corner points for projection
#Lon, Lat= map(grid[0], grid[1])
#
#maxlevel= np.max(np.abs(d.vort[t]))
#PlotVort(Lon, Lat, d.vort[t], map, maxlevel= maxlevel)
#
#
#plt.title('unfiltered')
#
#plt.figure(fignr)
#fignr += 1
#plt.clf()
#map= Polar_map(latsouth= d.lat[-1])
#maxlevel= np.max(np.abs( vortback[10:120]))
#
#PlotVort(Lon, Lat, vortback[10:120], map, maxlevel= maxlevel)
#
#plt.title('truncation'+str(ntrunc))
#
#    
#plt.figure(fignr)
#fignr += 1
#plt.clf()
#map= Polar_map(latsouth= d.lat[-1])
#maxlevel= np.max(np.abs( vort2back[10:120]))
#
#PlotVort(Lon, Lat, vort2back[10:120], map, maxlevel= maxlevel)
#
#plt.title('truncation'+str(ntrunc2))


plt.figure(fignr)
fignr += 1
plt.clf()
map= Polar_map(latsouth= d.lat[-1])
maxlevel= np.max(np.abs(vortback[10:120]- vort2back[10:120]))

PlotVort(Lon, Lat, vortback[10:120]- vort2back[10:120], map, maxlevel= maxlevel)

plt.title('truncation'+str(ntrunc)+'-'+str(ntrunc2))

PlotLocalMax(vortback[10:120]- vort2back[10:120], 2, 10, map, d.lon, d.lat, typ='max', value=False)
