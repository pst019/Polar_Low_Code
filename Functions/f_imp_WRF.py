#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

This is from plot_comparison to import data
"""

from netCDF4 import Dataset
import numpy as np
import time


class WRFdata:
    def __init__(self, filename, res=1):
        """res- to coarse the resolution: take just every res datapoint"""
        self.filename= filename
        print(filename)
        
        nc= Dataset(self.filename)
#        print(nc.variables.keys())
#        this one is even more interesting:
#        print(nc.variables.values())
        """should remove height above msl"""
        
        """reduce resolution"""
        self.res = res
        
        self.tim = nc.variables['XTIME'][:]/60 #in hours after the start of the simulation
        self.lat = nc.variables['XLAT'][0, :] #[::self.res, ::self.res]
        self.lon = nc.variables['XLONG'][0,:] #[::self.res, ::self.res]
        self.lev = nc.variables['LEV'][:]/100 # pressure level in hPa
        self.lat0= nc.TRUELAT2 #a latitude for the map projection (maybe take MOAD_CEN_LAT)- for map projection
        self.lon0= nc.STAND_LON #standard longitude - this is needed to make the map projection
#        self.press= nc.variables['bottom_top'][:]#[::-1]  ???
        """be careful the first variable in lat and lon is the y- variable and the second x: lat = (y, x)"""
#        print('levels:'+self.lev)
        
    def imp_standard(self):
        """import all levels for the most important data"""
        nc= Dataset(self.filename)

        # Variables
        self.u10= nc.variables['U10'][:]
        self.v10= nc.variables['V10'][:]
#        self.u= nc.variables['UU'][:]
#        self.v= nc.variables['VV'][:]
#        self.geo= nc.variables['GHT'][:]
        self.T= nc.variables['TT'][:]
        self.psfc = nc.variables['PSFC'][:]/100#[:, ::self.res, ::self.res]/100 #in hPa
        self.SST = nc.variables['SST'][:]
#        self.PV= nc.variables['PVO'][:]

        #u= nc.variables['U'] #this is actually in a different grid
#        self.u = nc.variables['x_wind_pl'][:][:, ::-1, ::self.res, ::self.res]
#        self.v = nc.variables['y_wind_pl'][:][:, ::-1, ::self.res, ::self.res]
    
        
#    def imp_level(self, pn, tn):
#        """import just the data of one level pn for one point of time tn"""     
#        print('presslevel: '+ str(self.press[pn])+ ', time: '+ str(time.localtime(self.tim[tn])[:4]))
#        pn = len(self.press) - pn -1 #since I inversed the pressure
#        
##        from datetime import datetime
##        chosentim= datetime(year, month, day, hour)
##        chosentim_s = (t-datetime(1970,1,1)).total_seconds()
##        chosent = np.where(chosentim_s = self.tim)
#        
#        nc= Dataset(self.filename)
#
#        # Variables
#        self.Temp = nc.variables['air_temperature_pl'][tn, pn][::self.res, ::self.res]
#        self.mslp = nc.variables['air_pressure_at_sea_level'][tn, 0][::self.res, ::self.res]/100 #in hPa
#        self.u = nc.variables['x_wind_pl'][tn, pn][::self.res, ::self.res]
#        self.v = nc.variables['y_wind_pl'][tn, pn][::self.res, ::self.res]        
#        self.PV = nc.variables['ertel_potential_vorticity_pl'][tn, pn][::self.res, ::self.res]
#        self.BL = nc.variables['atmosphere_boundary_layer_thickness'][tn, 0][::self.res, ::self.res]
#        self.geop = nc.variables['geopotential_pl'][tn, pn][::self.res, ::self.res]/9.81
#        self.relhum = nc.variables['relative_humidity_pl'][tn, pn][::self.res, ::self.res]
#        self.w = nc.variables['upward_air_velocity_pl'][tn, pn][::self.res, ::self.res]
#        if tn >= 1: self.prec = (nc.variables['precipitation_amount_acc'][tn+1, 0][::self.res, ::self.res] - nc.variables['precipitation_amount_acc'][tn-1, 0][::self.res, ::self.res] )/2       
#        if tn >= 1: self.LH = (-nc.variables['integral_of_surface_downward_latent_heat_evaporation_flux_wrt_time'][tn+1, 0][::self.res, ::self.res] + nc.variables['integral_of_surface_downward_latent_heat_evaporation_flux_wrt_time'][tn-1, 0][::self.res, ::self.res] +
#                               -nc.variables['integral_of_surface_downward_latent_heat_sublimation_flux_wrt_time'][tn+1, 0][::self.res, ::self.res] + nc.variables['integral_of_surface_downward_latent_heat_sublimation_flux_wrt_time'][tn-1, 0][::self.res, ::self.res] )/(2*60*60) #upward W/m**2       
#        if tn >= 1: self.SH = (-nc.variables['integral_of_surface_downward_sensible_heat_flux_wrt_time'][tn+1, 0][::self.res, ::self.res] + nc.variables['integral_of_surface_downward_sensible_heat_flux_wrt_time'][tn-1, 0][::self.res, ::self.res] )/(2*60*60) #upward W/m**2      
#
#        
#        
##    def imp_cross_sec(self, xn, tn):
##        nc= Dataset(self.filename)
##
##        # Variables
##        self.T = nc.variables['air_temperature_pl'][tn, :, xn][::-1, ::self.res]
###        self.mslp = nc.variables['air_pressure_at_sea_level'][tn, 0][::self.res, ::self.res]/100 #in hPa
##        self.u = nc.variables['x_wind_pl'][tn, :, xn][::-1, ::self.res]
##        self.v = nc.variables['y_wind_pl'][tn, :, xn][::-1, ::self.res]       
##        self.PV = nc.variables['ertel_potential_vorticity_pl'][tn, :, xn][::-1, ::self.res]
##        self.BL = nc.variables['atmosphere_boundary_layer_thickness'][tn, :, xn][::-1, ::self.res]
##        self.geop = nc.variables['geopotential_pl'][tn, :, xn][::-1, ::self.res]
##        self.relhum = nc.variables['relative_humidity_pl'][tn, :, xn][::-1, ::self.res]
##        self.w = nc.variables['upward_air_velocity_pl'][tn, :, xn][::-1, ::self.res]
#        
#        
#    def imp_cross_sec(self, xn, yn, tn):
#        """import data for a cross section following xn and yn at time tn
#        xn and yn from d.imp_cross_sec"""
#        nc= Dataset(self.filename)
#
#        # Variables
#        self.Temp = nc.variables['air_temperature_pl'][tn][::-1, ::self.res, ::self.res][:, yn, xn]
##        self.mslp = nc.variables['air_pressure_at_sea_level'][tn, 0][::self.res, ::self.res]/100 #in hPa
#        self.u = nc.variables['x_wind_pl'][tn][::-1, ::self.res, ::self.res][:, yn, xn]
#        self.v = nc.variables['y_wind_pl'][tn][::-1, ::self.res, ::self.res][:, yn, xn]       
#        self.PV = nc.variables['ertel_potential_vorticity_pl'][tn][::-1, ::self.res, ::self.res][:, yn, xn]
#        self.BL = nc.variables['atmosphere_boundary_layer_thickness'][tn, 0][::self.res, ::self.res][yn, xn]
#        self.geop = nc.variables['geopotential_pl'][tn][::-1, ::self.res, ::self.res][:, yn, xn]/9.81
#        self.relhum = nc.variables['relative_humidity_pl'][tn][::-1, ::self.res, ::self.res][:, yn, xn]
#        self.w = nc.variables['upward_air_velocity_pl'][tn][::-1, ::self.res, ::self.res][:, yn, xn]       
        
        

def mslp(geo, T):
    """calculate the mean sea level pressure from the 1000hPa geopotential height and the 1000hPa temperature"""
    return 1000* np.exp(9.81* geo/(287* T))