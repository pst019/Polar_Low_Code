#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

This is from plot_comparison to import data
"""

from netCDF4 import Dataset, num2date
import numpy as np
import time


class data:
    def __init__(self, filename, res= 5):
        """res - gives how much the resolution is decreased"""
        self.filename= filename
        
        nc= Dataset(self.filename)
        variables= nc.variables.keys()
#        print(nc.variables.keys())
        #this one is even more interesting:
#        print(nc.variables.values())
        """should remove height above msl"""
        
        """reduce resolution"""
        self.res = res
        
        self.tim = nc.variables['time'][:]
        self.datetime = num2date(nc.variables['time'][:], nc.variables['time'].units)       
        self.lat = nc.variables['latitude'][:][::self.res, ::self.res]
        self.lon = nc.variables['longitude'][:][::self.res, ::self.res]
        
#        if '_fp' in filename:
        if 'pressure' in variables:  
            self.pres= nc.variables['pressure'][:][::-1]
        elif 'pressure0' in variables:
            self.pres= nc.variables['pressure0'][:][::-1]
#            print('pressure0', self.pres)
#        else:
            #imp_standard, imp_surf, imp_atm ... do not work
        
        
        print(' relhum has been changed to RH and Temp to T and press to pres, SH for sensible heat should be changed to (SenH) otherwise it means specific humidity (SH)')
        """be careful the first variable in lat and lon is the y- variable and the second x: lat = (y, x)"""
        
    def imp_standard(self):
        """import all levels for the most important data -for fp"""
        nc= Dataset(self.filename)

        # Variables
        self.T = nc.variables['air_temperature_pl'][:][:, ::-1, ::self.res, ::self.res]
        self.mslp = nc.variables['air_pressure_at_sea_level'][:, 0][:, ::self.res, ::self.res]/100 #in hPa
        self.u = nc.variables['x_wind_pl'][:][:, ::-1, ::self.res, ::self.res] #
        self.v = nc.variables['y_wind_pl'][:][:, ::-1, ::self.res, ::self.res]
    


    def imp_surf(self, tn):
        """import just the data of the surface for one point of time tn -for fp"""     
#        print('presslevel: '+ str(self.pres[pn])+ ', time: '+ str(time.localtime(self.tim[tn])[:4]))
#        pn = len(self.pres) - pn -1 #since I inversed the pressure
        
        nc= Dataset(self.filename)

        # Variables
        self.T0m = nc.variables['air_temperature_0m'][tn, 0][::self.res, ::self.res]
        self.T2m = nc.variables['air_temperature_2m'][tn, 0][::self.res, ::self.res]        
        self.mslp = nc.variables['air_pressure_at_sea_level'][tn, 0][::self.res, ::self.res]/100 #in hPa
        self.u10m = nc.variables['x_wind_10m'][tn, 0][::self.res, ::self.res]
        self.v10m = nc.variables['y_wind_10m'][tn, 0][::self.res, ::self.res]        

    def imp_geop1000(self, tn):
        self.geop1000 = nc.variables['geopotential_pl'][tn, 0][::self.res, ::self.res]/9.81
     

    def imp_atm(self, tn):
        """import data that is only available for one level for the atmosphere -for fp"""
        nc= Dataset(self.filename)
        print(nc.variables.keys())

        self.cloud_conv = nc.variables['convective_cloud_area_fraction'][tn, 0][::self.res, ::self.res]
        self.cloud_high = nc.variables['high_type_cloud_area_fraction'][tn, 0][::self.res, ::self.res]
        self.cloud_med = nc.variables['medium_type_cloud_area_fraction'][tn, 0][::self.res, ::self.res]
        self.cloud_low = nc.variables['low_type_cloud_area_fraction'][tn, 0][::self.res, ::self.res]

        self.BL = nc.variables['atmosphere_boundary_layer_thickness'][tn, 0][::self.res, ::self.res]
        self.CAPE = nc.variables['specific_convective_available_potential_energy'][tn, 0][::self.res, ::self.res]  #J/kg
        self.CIN = nc.variables['atmosphere_convective_inhibition'][tn, 0][::self.res, ::self.res] #J/kg
        self.LCL = nc.variables['lifting_condensation_level'][tn, 0][::self.res, ::self.res]
        self.LFC = nc.variables['atmosphere_level_of_free_convection'][tn, 0][::self.res, ::self.res]

        
        #it is the flux between this and the next hour
        self.mslp = nc.variables['air_pressure_at_sea_level'][tn, 0][::self.res, ::self.res]/100 #in hPa
          
        self.prec = (nc.variables['precipitation_amount_acc'][tn+1, 0][::self.res, ::self.res] - nc.variables['precipitation_amount_acc'][tn, 0][::self.res, ::self.res] )       
        self.LatH = (-nc.variables['integral_of_surface_downward_latent_heat_evaporation_flux_wrt_time'][tn+1, 0][::self.res, ::self.res] + nc.variables['integral_of_surface_downward_latent_heat_evaporation_flux_wrt_time'][tn, 0][::self.res, ::self.res] +
                               -nc.variables['integral_of_surface_downward_latent_heat_sublimation_flux_wrt_time'][tn+1, 0][::self.res, ::self.res] + nc.variables['integral_of_surface_downward_latent_heat_sublimation_flux_wrt_time'][tn, 0][::self.res, ::self.res] )/(60*60) #upward W/m**2       
        self.SenH = (-nc.variables['integral_of_surface_downward_sensible_heat_flux_wrt_time'][tn+1, 0][::self.res, ::self.res] + nc.variables['integral_of_surface_downward_sensible_heat_flux_wrt_time'][tn, 0][::self.res, ::self.res] )/(60*60) #upward W/m**2      
        self.OLW = (-nc.variables['integral_of_toa_outgoing_longwave_flux_wrt_time'][tn+1, 0][::self.res, ::self.res] + nc.variables['integral_of_toa_outgoing_longwave_flux_wrt_time'][tn, 0][::self.res, ::self.res] )/(60*60) #upward W/m**2      
       
       
    def imp_level(self, pn, tn):
        """import just the data of one level pn for one point of time tn -for fp"""     
        print('presslevel: '+ str(self.pres[pn])+ ', time: '+ str(self.datetime[tn]))
        self.chosen_plevel = int(self.pres[pn])
        
        pn = len(self.pres) - pn -1 #since I inversed the pressure
        
        nc= Dataset(self.filename)

        # Variables
        self.T = nc.variables['air_temperature_pl'][tn, pn][::self.res, ::self.res]
        ### for older codes this one should be implemented again
#        self.mslp = nc.variables['air_pressure_at_sea_level'][tn, 0][::self.res, ::self.res]/100 #in hPa
        self.u = nc.variables['x_wind_pl'][tn, pn][::self.res, ::self.res]
        self.v = nc.variables['y_wind_pl'][tn, pn][::self.res, ::self.res]        
        self.PV = nc.variables['ertel_potential_vorticity_pl'][tn, pn][::self.res, ::self.res] #K m2 kg-1 s-1
        self.geop = nc.variables['geopotential_pl'][tn, pn][::self.res, ::self.res]/9.81
        self.RH = nc.variables['relative_humidity_pl'][tn, pn][::self.res, ::self.res]
        self.w = nc.variables['upward_air_velocity_pl'][tn, pn][::self.res, ::self.res]
        self.cloud = nc.variables['cloud_area_fraction_pl'][tn, pn][::self.res, ::self.res]
        
        if 'specific_humidity_pl' in nc.variables.keys():
            self.SH = nc.variables['specific_humidity_pl'][tn, pn][::self.res, ::self.res]
       

        
        
    def imp_cross_sec(self, xn, yn, tn):
        """import data for a cross section following xn and yn at time tn -for fp 
        xn and yn from d.imp_cross_sec
        I rotated xn and yn - maybe this introduced an error somewhere"""
        nc= Dataset(self.filename)

        # Variables
        self.T = nc.variables['air_temperature_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]
#        self.mslp = nc.variables['air_pressure_at_sea_level'][tn, 0][::self.res, ::self.res]/100 #in hPa
        self.u = nc.variables['x_wind_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]
        self.v = nc.variables['y_wind_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]       
        self.PV = nc.variables['ertel_potential_vorticity_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]
        self.BL = nc.variables['atmosphere_boundary_layer_thickness'][tn, 0][::self.res, ::self.res][xn, yn]
        self.geop = nc.variables['geopotential_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]/9.81
        self.RH = nc.variables['relative_humidity_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]
        self.w = nc.variables['upward_air_velocity_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]       

    def imp_cross_sec_reduced(self, xn, yn, tn):
        """import data for a cross section following xn and yn at time tn -for fp
        xn and yn from d.imp_cross_sec
        I rotated xn and yn - maybe this introduced an error somewhere"""
        nc= Dataset(self.filename)

        # Variables
        self.T = nc.variables['air_temperature_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]
        self.u = nc.variables['x_wind_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]
        self.v = nc.variables['y_wind_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]       
        self.RH = nc.variables['relative_humidity_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]
        if 'specific_humidity_pl' in nc.variables.keys():
            self.SH = nc.variables['specific_humidity_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]
        
    def imp_cross_sec_full(self, xn, yn, tn):
        """import data for a cross section following xn and yn at time tn -for fp
        xn and yn from d.imp_cross_sec
        I rotated xn and yn - maybe this introduced an error somewhere"""
        nc= Dataset(self.filename)

        # Variables
        self.T = nc.variables['air_temperature_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]
#        self.mslp = nc.variables['air_pressure_at_sea_level'][tn, 0][::self.res, ::self.res]/100 #in hPa
        self.u = nc.variables['x_wind_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]
        self.v = nc.variables['y_wind_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]       
        self.PV = nc.variables['ertel_potential_vorticity_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]
        self.BL = nc.variables['atmosphere_boundary_layer_thickness'][tn, 0][::self.res, ::self.res][xn, yn]
        self.geop = nc.variables['geopotential_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]/9.81
        self.RH = nc.variables['relative_humidity_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]
        if 'specific_humidity_pl' in nc.variables.keys():
            self.SH = nc.variables['specific_humidity_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]
        self.w = nc.variables['upward_air_velocity_pl'][tn][::-1, ::self.res, ::self.res][:, xn, yn]

        self.T0m = nc.variables['air_temperature_0m'][tn, 0][::self.res, ::self.res][xn, yn]
        self.T2m = nc.variables['air_temperature_2m'][tn, 0][::self.res, ::self.res] [xn, yn]       
        self.mslp = nc.variables['air_pressure_at_sea_level'][tn, 0][::self.res, ::self.res][xn, yn]/100 #in hPa
        self.u10m = nc.variables['x_wind_10m'][tn, 0][::self.res, ::self.res][xn, yn]
        self.v10m = nc.variables['y_wind_10m'][tn, 0][::self.res, ::self.res][xn, yn]
        
#        self.cloud_conv = nc.variables['convective_cloud_area_fraction'][tn, 0][::self.res, ::self.res][xn, yn]  #it is always 0
        self.cloud_high = nc.variables['high_type_cloud_area_fraction'][tn, 0][::self.res, ::self.res][xn, yn]
        self.cloud_med = nc.variables['medium_type_cloud_area_fraction'][tn, 0][::self.res, ::self.res][xn, yn]
        self.cloud_low = nc.variables['low_type_cloud_area_fraction'][tn, 0][::self.res, ::self.res][xn, yn]

        self.BL = nc.variables['atmosphere_boundary_layer_thickness'][tn, 0][::self.res, ::self.res][xn, yn]
        self.CAPE = nc.variables['specific_convective_available_potential_energy'][tn, 0][::self.res, ::self.res][xn, yn]  #J/kg
        self.CIN = nc.variables['atmosphere_convective_inhibition'][tn, 0][::self.res, ::self.res][xn, yn] #J/kg
        self.LCL = nc.variables['lifting_condensation_level'][tn, 0][::self.res, ::self.res][xn, yn]
        self.LFC = nc.variables['atmosphere_level_of_free_convection'][tn, 0][::self.res, ::self.res][xn, yn]

          
        if tn >= 1: self.prec = (nc.variables['precipitation_amount_acc'][tn+1, 0][::self.res, ::self.res][xn, yn] - nc.variables['precipitation_amount_acc'][tn-1, 0][::self.res, ::self.res][xn, yn] )/2  #mm/h
        if tn >= 1: self.LatH = (-nc.variables['integral_of_surface_downward_latent_heat_evaporation_flux_wrt_time'][tn+1, 0][::self.res, ::self.res][xn, yn] + nc.variables['integral_of_surface_downward_latent_heat_evaporation_flux_wrt_time'][tn-1, 0][::self.res, ::self.res][xn, yn] +
                               -nc.variables['integral_of_surface_downward_latent_heat_sublimation_flux_wrt_time'][tn+1, 0][::self.res, ::self.res][xn, yn] + nc.variables['integral_of_surface_downward_latent_heat_sublimation_flux_wrt_time'][tn-1, 0][::self.res, ::self.res][xn, yn] )/(2*60*60) #upward W/m**2       
        if tn >= 1: self.SenH = (-nc.variables['integral_of_surface_downward_sensible_heat_flux_wrt_time'][tn+1, 0][::self.res, ::self.res][xn, yn] + nc.variables['integral_of_surface_downward_sensible_heat_flux_wrt_time'][tn-1, 0][::self.res, ::self.res][xn, yn] )/(2*60*60) #upward W/m**2      
        if tn >= 1: self.LWsurf = (-nc.variables['integral_of_surface_net_downward_longwave_flux_wrt_time'][tn+1, 0][::self.res, ::self.res][xn, yn] + nc.variables['integral_of_surface_net_downward_longwave_flux_wrt_time'][tn-1, 0][::self.res, ::self.res][xn, yn] )/(2*60*60) #upward W/m**2      
        if tn >= 1: self.SWsurf = (-nc.variables['integral_of_surface_net_downward_shortwave_flux_wrt_time'][tn+1, 0][::self.res, ::self.res][xn, yn] + nc.variables['integral_of_surface_net_downward_shortwave_flux_wrt_time'][tn-1, 0][::self.res, ::self.res][xn, yn] )/(2*60*60) #upward W/m**2      
        if tn >= 1: self.LWTOA = (-nc.variables['integral_of_toa_outgoing_longwave_flux_wrt_time'][tn+1, 0][::self.res, ::self.res][xn, yn] + nc.variables['integral_of_toa_outgoing_longwave_flux_wrt_time'][tn-1, 0][::self.res, ::self.res][xn, yn] )/(2*60*60) #upward W/m**2      
        if tn >= 1: self.SWTOA = (-nc.variables['integral_of_toa_net_downward_shortwave_flux_wrt_time'][tn+1, 0][::self.res, ::self.res][xn, yn] + nc.variables['integral_of_toa_net_downward_shortwave_flux_wrt_time'][tn-1, 0][::self.res, ::self.res][xn, yn] )/(2*60*60) #upward W/m**2      




    def imp_pseudo(self, tn):
        """import data for a cross section following xn and yn at time tn -for the full files"""    
        filename= self.filename.replace('_fp', '')
        print(filename)
        nc= Dataset(filename)
#        print(nc.variables.values())
        
        self.tim_ps = nc.variables['time'][:]
        self.datetime_ps = num2date(nc.variables['time'][:], nc.variables['time'].units) 
        self.CTT = nc.variables['cloud_top_temperature'][tn, 0][::self.res, ::self.res]
        self.WVT = nc.variables['water_vapor_temperature'][tn, 0][::self.res, ::self.res]
        self.WVTC = nc.variables['water_vapor_temperature_with_clouds'][tn, 0][::self.res, ::self.res]      
        self.CWR = nc.variables['cloud_water_reflectivity'][tn, 0][::self.res, ::self.res]


    def imp_cross_sec_grib(self, xn, yn, tn):
        """import data for a cross section following xn and yn at time tn -for fp
        xn and yn from d.imp_cross_sec
        I rotated xn and yn - maybe this introduced an error somewhere"""
        filename= self.filename.replace('_fp', '')
        nc= Dataset(filename)

        # Variables
        self.T = nc.variables['air_temperature_ml'][tn][::-1, ::self.res, ::self.res][:, xn, yn]
        self.u = nc.variables['x_wind_ml'][tn][::-1, ::self.res, ::self.res][:, xn, yn]
        self.v = nc.variables['y_wind_ml'][tn][::-1, ::self.res, ::self.res][:, xn, yn]       
        self.SH = nc.variables['specific_humidity_ml'][tn][::-1, ::self.res, ::self.res][:, xn, yn]

        ap= nc.variables['ap'][:]
        b= nc.variables['b'][:]       
        ps= nc.variables['surface_air_pressure'][tn, :][::self.res,::self.res][:, xn, yn]
        print(ps)
#        self.pres= (ap.reshape([65, 1, 1]) + b.reshape([65, 1,1])* ps )[::-1] /100 #hPa
        self.pres= (ap.reshape([65]) + b.reshape([65])* ps )[::-1] /100 #hPa



    def imp_uv_ft(self, pn):
        """import the vorticity of level pn for the full time period"""     
#        print('presslevel: '+ str(self.pres[pn]) )
        self.chosen_plevel = int(self.pres[pn])
        
        pn = len(self.pres) - pn -1 #since I inversed the pressure       
        nc= Dataset(self.filename)

        # Variables
        self.u_ft = nc.variables['x_wind_pl'][:, pn][:, ::self.res, ::self.res]
        self.v_ft = nc.variables['y_wind_pl'][:, pn][:, ::self.res, ::self.res]                        

    def imp_surf_ft(self):
        """import the 10m winds and the mslp for the full time period"""
        nc= Dataset(self.filename)

        # Variables
        self.mslp_ft = nc.variables['air_pressure_at_sea_level'][:, 0][:, ::self.res, ::self.res]/100 #in hPa
        self.u10m_ft = nc.variables['x_wind_10m'][:, 0][:, ::self.res, ::self.res]
        self.v10m_ft = nc.variables['y_wind_10m'][:, 0][:, ::self.res, ::self.res]                         
#import from the internet

       
        #import numpy as np
#import netCDF4
#
#class data:
#    def __init__(self, year, month, day, hour):
#        self.year = year
#        self.month = month
#        self.day = day
#        self.hour = hour
#        
#        
#    def imp(self, Tsurf= False, Th= False, psurf= False, geoh= False, uh= False, vh= False, usurf= False, vsurf= False,
#            PVh= False, APE= False, wh= False ):
#        """import data from an online met.no catalogue
#        Tsurf - air temperature 0m
#        Th - temperature at pressure levels
#        psurf - air pressure at sea level
#        geoh - geopotential height at pressure levels (pl)
#        uh - east - west wind at pl
#        vh - north south wind at pl
#        wh - upward air velocity at pl
#        usurf - east - west 10m wind
#        vsurf - north south 10m wind
#        PVh - Ertel potential vorticity at pl
#        APE - specific convective available potential energy
#        """
#    
#        res='4'
#        res2= '4'
#        
#        url = ('http://thredds.met.no/thredds/dodsC/aromearcticraw/'+str(self.year)+'/'+str(self.month).zfill(2)+'/'+str(self.day).zfill(2)+
#        '/arome_arctic_extracted_2_5km_'+str(self.year)+str(self.month).zfill(2)+str(self.day).zfill(2)+'T'+str(self.hour).zfill(2)+
#        'Z.nc?time[0:1:66],pressure[0:1:10],longitude[0:'+res+':948][0:'+res+':738],latitude[0:'+res+':948][0:'+res+':738]'
#        + (',air_temperature_0m[0:1:0][0:1:0][0:'+res+':948][0:'+res+':738]' if Tsurf== True else '') #this is hopefully the surface temperature
#        + (',air_temperature_pl[0:1:0]['+str(Th)+':1:'+str(Th)+'][0:'+res+':948][0:'+res+':738]' if type(Th)== int else'') #this should be the air temperature at a given height
#        + (',air_pressure_at_sea_level[0:1:0][0:1:0][0:'+res+':948][0:'+res+':738]' if psurf== True else '') #surface pressure
#        + (',geopotential_pl[0:1:0]['+str(geoh)+':1:'+str(geoh)+'][0:'+res+':948][0:'+res+':738]' if type(geoh)== int else'') # at 
#        + (',x_wind_10m[0:1:0][0:1:0][0:'+res2+':948][0:'+res2+':738]' if usurf== True else '')
#        + (',y_wind_10m[0:1:0][0:1:0][0:'+res2+':948][0:'+res2+':738]' if vsurf== True else '')
#        + (',x_wind_pl[0:1:0]['+str(uh)+':1:'+str(uh)+'][0:'+res2+':948][0:'+res2+':738]' if type(uh)== int else'')
#        + (',y_wind_pl[0:1:0]['+str(vh)+':1:'+str(vh)+'][0:'+res2+':948][0:'+res2+':738]' if type(vh)== int else'')
#        + (',upward_air_velocity_pl[0:1:0]['+str(wh)+':1:'+str(wh)+'][0:'+res+':948][0:'+res+':738]'  if type(wh)== int else'') #this is relative strong for the last PL
##        + (',surface_geopotential[0:1:0][0:1:0][0:'+res+':948][0:'+res+':738]' if geosurf== True else '')
#        + (',ertel_potential_vorticity_pl[0:1:0]['+str(PVh)+':1:'+str(PVh)+'][0:'+res+':948][0:'+res+':738]' if type(PVh)== int else'')
#        + (',specific_convective_available_potential_energy[0:1:0][0:1:0][0:'+res+':948][0:'+res+':738]' if APE== True else '')
#        + (',projection_lambert')
#        )
#        
#        
#        print('get data from '+url)
#        
#        # create a dataset object
#        dataset = netCDF4.Dataset(url)
#         
#        # lookup a variable
#        self.time = dataset.variables['time'][:]
#        self.pres = dataset.variables['pressure'][:]
#        self.lon= dataset.variables['longitude'][:]
#        self.lat= dataset.variables['latitude'][:]
#        
#        if Tsurf == True: self.Tsurf= dataset.variables['air_temperature_0m'][0,0] - 273.15
#        if type(Th)== int: self.Th= dataset.variables['air_temperature_pl'][0,0] - 273.15
#        if psurf == True: self.psurf= dataset.variables['air_pressure_at_sea_level'][0,0] /100
#        if type(geoh) == int: self.geoh= dataset.variables['geopotential_pl'][0,0]/9.81
#
#        if usurf == True: self.usurf =dataset.variables['x_wind_10m'][:][0,0]
#        if vsurf == True: self.vsurf =dataset.variables['y_wind_10m'][:][0,0]
#        if type(uh) == int: self.uh =dataset.variables['x_wind_pl'][:][0,0]
#        if type(vh)== int: self.vh =dataset.variables['y_wind_pl'][:][0,0]
#        if type(wh)== int: self.wh =dataset.variables['upward_air_velocity_pl'][0,0]
#        #sgeo= dataset.variables['surface_geopotential'][0,0]/9.81
#        if type(PVh)== int: self.PVh = dataset.variables['ertel_potential_vorticity_pl'][0,0]
#        if APE== True: self.APE= dataset.variables['specific_convective_available_potential_energy'][0,0]
#    
#        self.dataset=dataset
##        print('height', self.pres[int(height)])