#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 10:01:45 2019

@author: pst019
"""


from scipy.interpolate import griddata
import great_circle_calculator.great_circle_calculator as gcc # pip install great-circle-calculator
import numpy as np

        
        
        
"""variables to specify"""
ncells = 10 # number of grid cells in each direction 
grid_dist= 25 #grid distance in km
radius= grid_dist * ncells


#lat- lon of your original grid as 2D array, this should be the CARA lat lon variable, 
longrid, latgrid=  ds.lat, ds.lon # or for ERA5: np.meshgrid(ds.lon, ds.lat), ds is my dataset

var = ds[variable] # the field variable you want to interpolate to the PL centred grid on the same grid as (longrid, latgrid)

        
center_lon= 50 #center location
center_lat= 70

#

bearing = 0 # direction of the track at a given time step as an angle compared to North.
"""could be calculated by the following: (center_lon is the longitude at the time you are looking at, center_lon_+1 at the next time step)


if trackpoint == 0:
    bearing= gcc.bearing_at_p1( (center_lon, center_lat, center_lon_+1, center_lat_+1) )
elif trackpoint == "last":
    bearing= gcc.bearing_at_p2( (center_lon_-1, center_lat_-1, center_lon, center_lat) ) 
else:
    bearing= gcc.bearing_at_p2( (center_lon_-1, center_lat_-1, center_lon_1, center_lat_1) )
"""




#create coordinates of the the tangental axis in the propogation
tanax= [list(gcc.point_given_start_and_bearing((center_lon, center_lat), bearing, n*grid_dist*1E3)) for n in np.arange(-ncells, ncells+1)]
tanax= np.array(tanax)

#the bearing along the tangential axis
point0= gcc.point_given_start_and_bearing((center_lon, center_lat), bearing, -(radius+grid_dist)*1E3) #point one before the start of the tangential axis, used for calculation of bearing_axis
bearing_axis= [gcc.bearing_at_p2((point0), (tanax[m,0], tanax[m,1])) for m in range(len(tanax))] #the wind direction along the tangential axis

#creation of the PL centred grid
center_grid= [[list(gcc.point_given_start_and_bearing((tuple(tanax[m])), (bearing_axis[m]-90)%360, n*grid_dist*1E3)) for m in range(len(tanax))] for n in np.arange(-ncells, ncells+1)]
center_grid= np.array(center_grid)




#this gives you the interpolated variable on the PL centred field
var_interp= griddata(tuple([np.ravel(longrid), np.ravel(latgrid)]) , np.ravel(var), (center_grid[:,:,0], center_grid[:,:,1]), method='linear')

