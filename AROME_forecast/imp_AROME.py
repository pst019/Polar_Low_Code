#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 15:07:08 2019

@author: pst019
"""

import os

user = os.getcwd().split('/')[2]

homedir = '/home/' + user + '/home/'
Mediadir = '/media/' + user + '/1692A00D929FEF8B/'

import numpy as np
from netCDF4 import Dataset
import pickle

import xarray as xr
import matplotlib.pyplot as plt


""" Specify the directory where you want to save the AROME-Arctic data."""
AROMEdir = Mediadir + 'PL/AA/operational/'

year = 2018
month = 11
day = 3
hour = 0

#url = ('http://thredds.met.no/thredds/dodsC/aromearcticarchive/' + str(year) + '/' + str(month).zfill(2) + '/' + str(
#    day).zfill(2) +
#       '/arome_arctic_full_2_5km_' + str(year) + str(month).zfill(2) + str(day).zfill(2) + 'T' + str(hour).zfill(
#            2) + 'Z.nc')

url = ('http://thredds.met.no/thredds/dodsC/aromearcticarchive/2018/02/10/arome_arctic_pp_2_5km_20180210T21Z.nc')
# url = ('http://thredds.met.no/thredds/fileServer/metusers/eivinds/AA_Patrick/AS500_2018021400.nc')
# url = (r"http://thredds.met.no/thredds/dodsC/aromearcticarchive/2018/12/31/arome_arctic_pp_2_5km_20181231T21Z.nc")

# url= ('http://thredds.met.no/thredds/fileServer/aromearcticarchive/2018/04/06/arome_arctic_full_2_5km_20180406T21Z.nc')
# url = ('http://thredds.met.no/thredds/dodsC/aromearcticraw/'+str(year)+'/'+str(month).zfill(2)+'/'+str(day).zfill(2)+
#        '/arome_arctic_extracted_2_5km_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'T'+str(hour).zfill(2)+
#        'Z.nc')
#url= ('http://thredds.met.no/thredds/dodsC/met.no/observations/stations/SN99938.nc')

#url = ('http://thredds.met.no/thredds/dodsC/aromearcticarchive/2018/12/31/arome_arctic_pp_2_5km_20181231T00Z.nc')

#dataset = Dataset(url) #+'?air_temperature_2m,time,longitude,latitude,height1')  #can specify the variables we want to import

ds = xr.open_dataset(url)  # , decode_times=False)


plt.figure(1)


ds.close()

# """ To get different datasets. Day and hour specify the start time of the AROME simulation. Data is obtained for the whole 66h? of the AROME simulation. """
##for day in range(18, 19):
##    for hour in [0]: #, 6, 12, 18]:
#        print(day, hour)
#
#        """ To get one of the horizontal fields """
#        #extracted
##        url = ('http://thredds.met.no/thredds/dodsC/aromearcticarchive/2018/02/'+str(day).zfill(2)+'/arome_arctic_extracted_2_5km_201802'+str(day).zfill(2)+'T'+str(hour).zfill(2)+'Z.nc')
#        #full - has more variables, but is slower
#
#        dataset = Dataset(url) #+'?air_temperature_2m,time,longitude,latitude,height1')  #can specify the variables we want to import
#
#        """if you want to know which fields are available and a bit more desciption"""
#        #print(dataset.variables.keys())
#        #this one is even more interesting:
#        #print(dataset.variables.values())
#        #lookup a variable
#        #time = dataset.variables['time'][:]
#
#        filename= AROMEdir + 'AROME_'+str(day).zfill(2)+'_02_2018_'+str(hour).zfill(2)
##
##        t = num2date(dataset.variables['time'][:], dataset.variables['time'].units)
##        pickle.dump(t, open(filename+'_t.p', 'wb'))
##
##        T2= dataset.variables['air_temperature_2m'][:] -273.15
##        pickle.dump(T2, open(filename+'_T2.p', 'wb'))
##
#        q2= dataset.variables['relative_humidity_2m'][:]
#        pickle.dump(q2, open(filename+'_RH.p', 'wb'))
##
##        SP=dataset.variables['surface_air_pressure'][:]
##        pickle.dump(SP, open(filename+'_SP.p', 'wb'))
##
##        u10= dataset.variables['x_wind_10m'][:]
##        pickle.dump(u10, open(filename+'_u10.p', 'wb'))
##
##        v10= dataset.variables['y_wind_10m'][:]
##        pickle.dump(v10, open(filename+'_v10.p', 'wb'))
##
##        LWR= dataset.variables['integral_of_surface_net_downward_longwave_flux_wrt_time'][:]
##        pickle.dump(LWR, open(filename+'_LWR.p', 'wb'))
#
##        LWR_down= dataset.variables['integral_of_surface_downwelling_longwave_flux_in_air_wrt_time'][:]
##        pickle.dump(LWR_down, open(filename+'_LWR_down.p', 'wb'))
#
##        SWR= dataset.variables['integral_of_surface_net_downward_shortwave_flux_wrt_time'][:]
##        pickle.dump(SWR, open(filename+'_SWR.p', 'wb'))
##
##        SWR_down= dataset.variables['integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time'][:]
##        pickle.dump(SWR_down, open(filename+'_SWR_down.p', 'wb'))
#
##        cloud= dataset.variables['cloud_area_fraction'][:]
##        pickle.dump(cloud, open(filename+'_cloud.p', 'wb'))
#
##        SH= dataset.variables['integral_of_surface_downward_sensible_heat_flux_wrt_time'][:]
##        pickle.dump(SH, open(filename+'_SH.p', 'wb'))
##
##        LH_ev= dataset.variables['integral_of_surface_downward_latent_heat_evaporation_flux_wrt_time'][:]
##        pickle.dump(LH_ev, open(filename+'_LH_ev.p', 'wb'))
##
##        LH_su= dataset.variables['integral_of_surface_downward_latent_heat_sublimation_flux_wrt_time'][:]
##        pickle.dump(LH_su, open(filename+'_LH_su.p', 'wb'))
#
##        lat= dataset.variables['latitude'][:]
##        lon= dataset.variables['longitude'][:]
##        pickle.dump([lat, lon], open(AROMEdir + 'AROME_lat-lon.p', 'wb'))
#
#
#        """get improved vertical profile for one grid cell / station (for Thea)"""
##        url = ('http://thredds.met.no/thredds/dodsC/aromearcticarchive/2018/02/'+str(day).zfill(2)+'/arome_arctic_full_2_5km_201802'+str(day).zfill(2)+'T'+str(hour).zfill(2)+'Z.nc')
##        dataset = Dataset(url) #+'?air_temperature_2m,time,longitude,latitude,height1')  #can specify the variables we want to import
#
##        lat= dataset.variables['latitude'][:]
##        lon= dataset.variables['longitude'][:]
##
##        statname='NyAle'
##        latstat= 78.9232
##        lonstat= 11.9231
#
##        statname='Jan'
##        latstat= 70.9396
##        lonstat= -8.6679
#
##        statname='Bjorn'
##        latstat= 74.5037
##        lonstat= 19.0011
#
##        statname='And'
##        latstat= 69.3152
##        lonstat= 16.1309
##
##        dist_lat= (np.cos(np.deg2rad(latstat))*(lon-lonstat))**2 + (lat-latstat)**2
##        xpos, ypos= np.where(dist_lat == np.min(dist_lat))
##
##        ap= dataset.variables['ap'][:]
##        b= dataset.variables['b'][:]
##
##        ps= dataset.variables['surface_air_pressure'][:,:, xpos, ypos]
##        T_m= dataset.variables['air_temperature_ml'][:,:, xpos, ypos]
##
##        p= (ap.reshape([1, 65, 1, 1]) + b.reshape([1, 65, 1,1])* ps ) /100 #hPa
##
##        filename= AROMEdir + 'AROME_'+str(day).zfill(2)+'_02_2018_'+str(hour).zfill(2)+'_Tp_'+statname+'.p'
##        pickle.dump([p, T_m], open(filename, 'wb'))
#
#
#        """get improved profile for several grid cells in the sourounding of one station"""
##        url = ('http://thredds.met.no/thredds/dodsC/aromearcticarchive/2018/02/'+str(day).zfill(2)+'/arome_arctic_full_2_5km_201802'+str(day).zfill(2)+'T'+str(hour).zfill(2)+'Z.nc')
##        dataset = Dataset(url) #+'?air_temperature_2m,time,longitude,latitude,height1')  #can specify the variables we want to import
##
##        lat= dataset.variables['latitude'][:]
##        lon= dataset.variables['longitude'][:]
##
##        statname='Adv'
##        latstat= 78+12.10/60
##        lonstat= 15+49.41/60
##
##        dist_lat= (np.cos(np.deg2rad(latstat))*(lon-lonstat))**2 + (lat-latstat)**2
##        xpos, ypos= np.where(dist_lat == np.min(dist_lat))
##        xpos, ypos= xpos[0], ypos[0]
##
##        ap= dataset.variables['ap'][:]
##        b= dataset.variables['b'][:]
##
##        ps= dataset.variables['surface_air_pressure'][:24,:, xpos: xpos+3, ypos-5: ypos+1 ]
##        T_m= dataset.variables['air_temperature_ml'][:24,:, xpos: xpos+3, ypos-5: ypos+1]
##        q_m= dataset.variables['specific_humidity_ml'][:24,:, xpos: xpos+3, ypos-5: ypos+1]
##        x_wind_m= dataset.variables['x_wind_ml'][:24,:, xpos: xpos+3, ypos-5: ypos+1]
##        y_wind_m= dataset.variables['y_wind_ml'][:24,:, xpos: xpos+3, ypos-5: ypos+1]
##        lon= lon[xpos: xpos+3, ypos-5: ypos+1]
##        lat= lat[xpos: xpos+3, ypos-5: ypos+1]
##
##        p= (ap.reshape([1, 65, 1, 1]) + b.reshape([1, 65, 1,1])* ps ) /100 #hPa
##
##        filename= AROMEdir + 'AROME_'+str(day).zfill(2)+'_02_2018_'+str(hour).zfill(2)+'_around_'+statname+'.p'
##        pickle.dump([p, lon, lat, T_m, q_m, x_wind_m, y_wind_m], open(filename, 'wb'))
##
##        filename= AROMEdir + 'AROME_'+str(day).zfill(2)+'_02_2018_'+str(hour).zfill(2)+'_around_'+statname+'2.p'
##        pickle.dump([p, lon, lat, T_m, q_m, x_wind_m, y_wind_m], open(filename, 'wb'), protocol=2)
#
#        dataset.close()


#    """ get temp profile for one grid cells - this can be deleted """
#    latNyAle= 78.9232
#    lonNyAle= 11.9231
##
##    #latJan= 70.9396
##    #lonJan= -8.6679
##    #
##    #latBjorn= 74.5037
##    #lonBjorn= 19.0011
##    #
##    #LatAnd= 69.3152
##    #lonAnd= 16.1309
##
#    dist_lat= (np.cos(np.deg2rad(latNyAle))*(lon-lonNyAle))**2 + (lat-latNyAle)**2
#    xpos, ypos= np.where(dist_lat == np.min(dist_lat))
##
#    p= dataset.variables['pressure'][:]
#    T_p= dataset.variables['air_temperature_pl'][:,:, xpos, ypos]
##    #
##    print('write')
##
#    filename= AROMEdir + 'AROME_'+str(day).zfill(2)+'_02_2018_'+str(hour).zfill(2)+'_Tp_NyAle.p'
#    pickle.dump([p, T_p], open(filename, 'wb'))


""" To get general_inf - information on the altitude and the land-area fraction"""
# url = ('http://thredds.met.no/thredds/dodsC/aromearcticarchive/2018/02/10/arome_arctic_pp_2_5km_20180210T21Z.nc')
# dataset = Dataset(url)  #can specify the variables we want to import
##z= dataset.variables['altitude'][:]
##
##lon= dataset.variables['longitude'][:]
##lat= dataset.variables['latitude'][:]
##
##pickle.dump([lon, lat, z], open(AROMEdir + 'AROME_general_inf.p', 'wb'))
#
# land= dataset.variables['land_area_fraction'][:]
# pickle.dump(land, open(AROMEdir + 'AROME_land.p', 'wb'))
# dataset.close()
