# -*- coding: utf-8 -*-
"""
Created on Sat Apr 16 10:46:59 2016

@author: pst019
"""

import numpy as np
import netCDF4

class data:
    def __init__(self, year, month, day, hour):
        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        
        
    def imp(self, Tsurf= False, Th= False, psurf= False, geoh= False, uh= False, vh= False, usurf= False, vsurf= False,
            PVh= False, APE= False, wh= False ):
        """import data from an online met.no catalogue
        Tsurf - air temperature 0m
        Th - temperature at pressure levels
        psurf - air pressure at sea level
        geoh - geopotential height at pressure levels (pl)
        uh - east - west wind at pl
        vh - north south wind at pl
        wh - upward air velocity at pl
        usurf - east - west 10m wind
        vsurf - north south 10m wind
        PVh - Ertel potential vorticity at pl
        APE - specific convective available potential energy
        """
    
        res='4'
        res2= '4'
        
        url = ('http://thredds.met.no/thredds/dodsC/aromearcticraw/'+str(self.year)+'/'+str(self.month).zfill(2)+'/'+str(self.day).zfill(2)+
        '/arome_arctic_extracted_2_5km_'+str(self.year)+str(self.month).zfill(2)+str(self.day).zfill(2)+'T'+str(self.hour).zfill(2)+
        'Z.nc?time[0:1:66],pressure[0:1:10],longitude[0:'+res+':948][0:'+res+':738],latitude[0:'+res+':948][0:'+res+':738]'
        + (',air_temperature_0m[0:1:0][0:1:0][0:'+res+':948][0:'+res+':738]' if Tsurf== True else '') #this is hopefully the surface temperature
        + (',air_temperature_pl[0:1:0]['+str(Th)+':1:'+str(Th)+'][0:'+res+':948][0:'+res+':738]' if type(Th)== int else'') #this should be the air temperature at a given height
        + (',air_pressure_at_sea_level[0:1:0][0:1:0][0:'+res+':948][0:'+res+':738]' if psurf== True else '') #surface pressure
        + (',geopotential_pl[0:1:0]['+str(geoh)+':1:'+str(geoh)+'][0:'+res+':948][0:'+res+':738]' if type(geoh)== int else'') # at 
        + (',x_wind_10m[0:1:0][0:1:0][0:'+res2+':948][0:'+res2+':738]' if usurf== True else '')
        + (',y_wind_10m[0:1:0][0:1:0][0:'+res2+':948][0:'+res2+':738]' if vsurf== True else '')
        + (',x_wind_pl[0:1:0]['+str(uh)+':1:'+str(uh)+'][0:'+res2+':948][0:'+res2+':738]' if type(uh)== int else'')
        + (',y_wind_pl[0:1:0]['+str(vh)+':1:'+str(vh)+'][0:'+res2+':948][0:'+res2+':738]' if type(vh)== int else'')
        + (',upward_air_velocity_pl[0:1:0]['+str(wh)+':1:'+str(wh)+'][0:'+res+':948][0:'+res+':738]'  if type(wh)== int else'') #this is relative strong for the last PL
#        + (',surface_geopotential[0:1:0][0:1:0][0:'+res+':948][0:'+res+':738]' if geosurf== True else '')
        + (',ertel_potential_vorticity_pl[0:1:0]['+str(PVh)+':1:'+str(PVh)+'][0:'+res+':948][0:'+res+':738]' if type(PVh)== int else'')
        + (',specific_convective_available_potential_energy[0:1:0][0:1:0][0:'+res+':948][0:'+res+':738]' if APE== True else '')
        + (',projection_lambert')
        )
        
        
        print('get data from '+url)
        
        # create a dataset object
        dataset = netCDF4.Dataset(url)
         
        # lookup a variable
        self.time = dataset.variables['time'][:]
        self.press = dataset.variables['pressure'][:]
        self.lon= dataset.variables['longitude'][:]
        self.lat= dataset.variables['latitude'][:]
        
        if Tsurf == True: self.Tsurf= dataset.variables['air_temperature_0m'][0,0] - 273.15
        if type(Th)== int: self.Th= dataset.variables['air_temperature_pl'][0,0] - 273.15
        if psurf == True: self.psurf= dataset.variables['air_pressure_at_sea_level'][0,0] /100
        if type(geoh) == int: self.geoh= dataset.variables['geopotential_pl'][0,0]/9.81

        if usurf == True: self.usurf =dataset.variables['x_wind_10m'][:][0,0]
        if vsurf == True: self.vsurf =dataset.variables['y_wind_10m'][:][0,0]
        if type(uh) == int: self.uh =dataset.variables['x_wind_pl'][:][0,0]
        if type(vh)== int: self.vh =dataset.variables['y_wind_pl'][:][0,0]
        if type(wh)== int: self.wh =dataset.variables['upward_air_velocity_pl'][0,0]
        #sgeo= dataset.variables['surface_geopotential'][0,0]/9.81
        if type(PVh)== int: self.PVh = dataset.variables['ertel_potential_vorticity_pl'][0,0]
        if APE== True: self.APE= dataset.variables['specific_convective_available_potential_energy'][0,0]
    
        self.dataset=dataset
#        print('height', self.press[int(height)])