#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 15:06:41 2017

@author: pst019
"""

from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

from f_imp_ERA2 import *
import sys  #to import the functions from a different directory
sys.path.insert(0, '/home/'+user+'/codeAROME/')
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
import scipy.ndimage.filters as filters

year= 2008
#month= 3
#day= 3
#hour = 240#3 3 -6
month, day, hour= 3, 4, 12
t= ((day-1)*24+hour)//6 #hours since beginning of month

"""import data"""
#d= data('Vort', year, month)
#vortfilter2= FourierFilter3D(d.lon, d.lat, d.vort)



lon= d.lon
lat= d.lat
vort= d.vort
T_low= 20
T_up= 150

from scipy.fftpack import fftfreq, fft2, ifft2
transf= fft2(vort[t])

dlat, dlon= lat[0]-lat[1], lon[1]-lon[0]
ltim= np.shape(vort)[0]

wny= fftfreq(len(lat), dlat/360) #wavenr in y
wny= np.tile(wny, (len(lon), 1)).T #dublication of the wavenr in y
#wny= np.tile(wny, (ltim, 1, 1))
#wnx= np.array([list(fftfreq(len(lon), dlon*np.cos(np.deg2rad(la))/360)) for la in range(len(lat))]) # wavenr in x
wnx= np.array([list(fftfreq(len(lon), dlon*np.cos(np.deg2rad(la))/360)) for la in lat]) # wavenr in x

#wnx= np.tile(wnx, (ltim, 1, 1))

transf[np.sqrt(wnx**2) <T_low] = 0
transf[np.sqrt(wnx**2) >T_up] = 0

#transf[np.sqrt(wny**2) <T_low] = 0
#transf[np.sqrt(wny**2) >T_up] = 0

#transf[np.abs(wnx+wny) <T_low] = 0
#transf[np.abs(wnx+wny) >T_up] = 0

#transf[np.sqrt(wnx**2+wny**2) <T_low] = 0
#transf[np.sqrt(wnx**2+wny**2) >T_up] = 0

#transf[np.sqrt(wnx**2+wny**2) <0.8* T_low] = 0
#transf[np.sqrt(wnx**2+wny**2) >1.2* T_up] = 0
#
#transf[np.sqrt(wnx**2+wny**2) <1.2* T_low] /= 2
#transf[np.sqrt(wnx**2+wny**2) >0.8* T_up] /= 2

#transf[np.sqrt(wnx**2+wny**2) <0.8* T_low] = 0
#transf[np.sqrt(wnx**2+wny**2) <0.9* T_low] *= .32
#transf[np.sqrt(wnx**2+wny**2) <1* T_low] *= .625
#transf[np.sqrt(wnx**2+wny**2) <1.2* T_low] *= .8
#       
#transf[np.sqrt(wnx**2+wny**2) >1.2* T_up] = 0
#transf[np.sqrt(wnx**2+wny**2) >1.1* T_up] *= 0.32
#transf[np.sqrt(wnx**2+wny**2) >1* T_up] *= 0.625
#transf[np.sqrt(wnx**2+wny**2) >0.8* T_up] *= 0.8  
       
vortfilter= np.real(ifft2(transf))       


"""real fft"""
#from numpy.fft import rfftfreq, rfft2, irfft2
#transf= rfft2(vort[t])
#
#dlat, dlon= lat[0]-lat[1], lon[1]-lon[0]
#ltim= np.shape(vort)[0]
#
#wny= rfftfreq(len(lat), dlat/360) #wavenr in y
#wny= np.tile(wny, (len(lon), 1)).T #dublication of the wavenr in y
##wny= np.tile(wny, (ltim, 1, 1))
##wnx= np.array([list(fftfreq(len(lon), dlon*np.cos(np.deg2rad(la))/360)) for la in range(len(lat))]) # wavenr in x
#wnx= np.array([list(rfftfreq(len(lon), dlon*np.cos(np.deg2rad(la))/360)) for la in lat]) # wavenr in x
#
#
#transf[np.sqrt(wnx**2+wny**2) <T_low] = 0
#transf[np.sqrt(wnx**2+wny**2) >T_up] = 0
#
#       
#vortfilter= np.real(irfft2(transf))       


plt.figure(6)
plt.clf()
l= -13
print('lat', lat[l])

plt.plot(lon, vort[t, l, :])
plt.plot(lon, vortfilter[l, :])


#"""test"""
transfx= fft(vort[t, l, :])
wnxx= fftfreq(len(lon), dlon*np.cos(np.deg2rad(lat[l]))/360) # wavenr in x
transfx[np.sqrt(wnxx**2) <T_low] = 0
transfx[np.sqrt(wnxx**2) >T_up] = 0

vortfilterx= np.real(ifft(transfx)) 
plt.plot(lon, vortfilterx[:])


transfx= rfft(vort[t, l, :])
wnxx= rfftfreq(len(lon), dlon*np.cos(np.deg2rad(lat[l]))/360) # wavenr in x
transfx[np.sqrt(wnxx**2) <T_low] = 0
transfx[np.sqrt(wnxx**2) >T_up] = 0

vortfilterx= irfft(transfx)
plt.plot(lon, vortfilterx[:])


plt.figure(7)
plt.clf()

plt.plot(wnxx, transfx)


plt.figure(8)
plt.clf()

plt.plot(wnx[l], transf[l])

#plt.figure(7)
#plt.clf()
#l= -50 #-290
#print('lon', lon[l])
#
#plt.plot(lat, vort[t, :, l])
#plt.plot(lat, vortfilter[:, l])



"""PL detection"""
#maxima= LocalMax(vortfilter, neighborhood_size= (1,10,15), threshold= 5)
#maxima= LocalMax_dist(vortfilter, d.lat, distlat= 5, threshold= 2)

#maxima[d.SST < 272]= 0
#maxima[d.SST- d.T500 < 43]= 0



"""Plot"""

#plt.figure(1)
#plt.clf()
#map= Polar_map(latsouth= d.lat[-1])
#grid= np.meshgrid(d.lon, d.lat) #pcolormesh needs corner points for projection
#Lon, Lat= map(grid[0], grid[1])
#
#PlotVort(Lon, Lat, d.vort[t], map)
##maxpoints3D= np.where(maxima == True)
##maxpoints= (maxpoints3D[1][np.where(maxpoints3D[0]==t)[0]], maxpoints3D[2][np.where(maxpoints3D[0]==t)[0]])
##PointsOnMap(maxpoints, d.lon, d.lat, map)
#plt.title('unfiltered vorticity')
#

plt.figure(2)
plt.clf()
map= Polar_map(latsouth= d.lat[-1])

PlotVort(Lon, Lat, vortfilter, map)
#maxpoints3D= np.where(maxima == True)
#maxpoints= (maxpoints3D[1][np.where(maxpoints3D[0]==t)[0]], maxpoints3D[2][np.where(maxpoints3D[0]==t)[0]])
#PointsOnMap(maxpoints, d.lon, d.lat, map)
plt.title('filtered vorticity')


#plt.figure(2)
#plt.clf()
#map= Polar_map(latsouth= d.lat[-1])
#
##maxima= LocalMax(vortfilter2, neighborhood_size= (1,5,7), threshold= 4)
#
#
#PlotVort(Lon, Lat, vortfilter2[t], map)
##maxpoints3D= np.where(maxima == True)
##maxpoints= (maxpoints3D[1][np.where(maxpoints3D[0]==t)[0]], maxpoints3D[2][np.where(maxpoints3D[0]==t)[0]])
##PointsOnMap(maxpoints, d.lon, d.lat, map)
#plt.title('filtered vorticity 2')

