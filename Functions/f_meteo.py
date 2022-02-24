#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 12:45:26 2017

@author: pst019
"""

import numpy as np

def PotTemp(T, plev):
    """calculate the potential temperature
    both T and plev can be scalars or arrays (of same dimension)
    T given in K
    plev given in hPa"""
    RCp= 2/7
    return T* (1000/plev)**RCp


def EquiPotTemp(T, sHum, plev, approx=True):
    """calculate the equivalent potential temperature
    all variables can be scalars or arrays (of same dimension)
    T given in K
    sHum given in ERA unit [kg/m**3]? (also called water vapor mixing ratio)
    plev given in hPa
    approx - True gives an approximation good around 0C
    approx - False - not tested - is a more accurate formula from https://en.wikipedia.org/wiki/Latent_heat"""
    if approx: Lc, Cp= 2501, 1.006
    else:
        Lc= 2500.8 - 2.36*(T-273.15) + 0.0016*(T-273.15)**2 - 0.00006*(T-273.15)**3
        Cp= 1.006 #actually this also varies with T and p
        
    return PotTemp(T, plev)* np.exp(sHum* Lc / (T * Cp))


def RH2SH(RH, p, T):
    """translate from relative humidity (RH) in % to specific humditiy (SH) in g/kg
    at pressure p in hPa and temperature T in K
    rel humidity to specific humidity: http://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
    spechum= relhum[0to1]*100/(.263 * press[Pa]) * exp(17.67*(T[k]-273.15)/(T[K]-29.65))"""
#    return RH/(.263 *p)* np.exp(17.67*(T-273.15) / (T-29.65))
    SH= RH/(.263 * (p*100) )* np.exp(17.67*(T-273.15) / (T-29.65))
    return 1000* SH

#def RHfromSH(SH, p, T ):
#    """using metpy"""
#    from metpy.units import units
#    from metpy.calc.thermo import thermo
#    RH = thermo.relative_humidity_from_specific_humidity(SH, T* units.K, p*units.hPa)
#    return RH


def SH2RH(SH, p, T):
    """translate from specific humidity (SH) in g/kg to relative humditiy (RH) in %
    at pressure p in hPa and temperature T in K 
    https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity     
    
    This is approx equivalent to the metpy function:
    AA_var= relative_humidity_from_specific_humidity(AA_SH, T_m* units.K, p*units.Pa) *100    
    AA_var= SH2RH(AA_SH*1E3, p/100, T_m)
    p *100 #convert hPa to Pa
    SH = 1E3 # convert g/kg to kg/kg"""
    
    RH=  0.263 * (p*100)* (SH/1E3) * 1/(np.exp(17.67*(T- 273.15) / (T-29.65) ) )
    return RH


def Wind2ZonalMeridional(lon, lat, u, v):
    """rotate the x,y wind vector (u,v) for the coordinates lon, lat to zonal and meridional wind
    tested for AROME only"""
    
    del_lon= np.gradient(lon, axis= 0)* np.cos(np.deg2rad(lat)) #find the difference in longitude between grid cells along x-axis
    del_lat= np.gradient(lat, axis= 0)
    
    rot= np.rad2deg(np.arctan2(del_lon, del_lat)) #how much to rotate the gridcell anticlockwise to be oriented N-S    
    dir_raw= np.rad2deg(np.arctan2(u, v)) #the wind direction with respect to the grid (0 being upwards, 90 rightwards, -90 leftw) 
#    print(np.shape(dir_raw), np.shape(rot))
    dir_rot= (dir_raw + rot)%360 #the wind direction after the rotation. Now it is the direction with respect to the North pole (0 Northwards, 90 E, 270 W)

    wind_vel= np.sqrt(u**2 +v**2)
    u_zonal= wind_vel* np.sin(np.deg2rad(dir_rot))
    v_meridional= wind_vel* np.cos(np.deg2rad(dir_rot))
  
    return u_zonal, v_meridional


def WindSpeedDirection2UV(wind_vel, wind_dir, orientation='comming from'):
    """translate the wind speed and the direction [in degree] to the zonal (u) and meridional (v) wind
    if orientation= 'comming from', the wind_dir is measuring from where the wind is coming from - otherwise where it is "going to" """
    """for STARS_loc_SOM_single_var_2 it appears that u_zonal and v_meridional are exchanged!!"""
    
    if orientation== 'comming from':  wind_dir = (wind_dir+ 180)%360
    
    u_zonal= wind_vel* np.sin(np.deg2rad(wind_dir))
    v_meridional= wind_vel* np.cos(np.deg2rad(wind_dir))
    return u_zonal, v_meridional


def UV2Direction(u, v):
    """calculates the wind direction [in degree] from the u and v wind vector """
    dir_rad= np.arctan2(u, v)
    return np.rad2deg(dir_rad)



def CalculateThermalWind(thick, lat, lon):
    """calculates the thermal wind from the geopotential thickness (2D or 3D)
    lat and lon are 1D arrays"""
    
    a= 6370E3
    omega=2*np.pi/(24*60**2)
    f= 2*omega* np.sin(np.deg2rad(lat))
    
    if len(np.shape(thick)) ==3:
        thick_dlat= -np.gradient(thick, axis=1)
        thick_dlon= np.gradient(thick, axis=2)
    elif len(np.shape(thick)) ==2:
        thick_dlat= -np.gradient(thick, axis=0)
        thick_dlon= np.gradient(thick, axis=1)        
    
    dlat=np.deg2rad(lat[0]- lat[1])
    dlon=np.deg2rad(lon[1]- lon[0])
    
    u_T= -1/(f[:,None]) *1/(a*dlat)* thick_dlat
    v_T= 1/(f[:,None]) *1/(a*dlon)* 1/(np.cos(np.deg2rad(lat))[:,None])* thick_dlon 
    return u_T, v_T


def LatentHeat_fromSnow(snow):
    """calculates the release of latent heat [W/m**2] from precipiation in form of snow
    snow [mm/(m**2 h)= kg/(m**2 h)]"""
    L=2834E3 #J/kg - latent heat release from gas to solid at 0 C - does not vary much with T
    
    return L*snow/60**2 #J/kg * kg/m**2 h * h/s = J/s /m**2 = W/m**2

def LatentHeat_fromRain(rain):
    """calculates the release of latent heat [W/m**2] from precipiation in form of snow
    snow [mm/(m**2 h)= kg/(m**2 h)]"""
    L=2501E3 #J/kg - latent heat release from gas to water at 0 C - could be more precise if T is included: see EquiPotTemp
    
    return L*rain/60**2 #J/kg * kg/m**2 h * h/s = J/s /m**2 = W/m**2



def CalcVort(u, v, dx= 2.5, filter=True, filterdist= 60, gaustruncate= 1):
    """calculate the vorticity from the horizontal wind vectors - tested for AA
    The vorticity can also be filtered
    u,v - can be 2D or 3D. If 3D, the first dimension is time 
    dx- grid spacing [km]
    filterdist - distance of the filter [km]"""
    
    vort= (np.gradient(v, dx*1E3, axis= -1)- np.gradient(u, dx*1E3, axis= -2))*1E5         #axis=1 is derivative in x direction, axis= 0 in y direction

    if filter:
        import scipy.ndimage.filters as filters        
#        vort= filters.gaussian_filter(vort, sigma= filterdist/dx, mode='nearest', truncate= gaustruncate) - for 2D
        vort= filters.gaussian_filter1d(vort, axis= -1, sigma= filterdist/dx, mode='constant', truncate= gaustruncate)
        vort= filters.gaussian_filter1d(vort, axis= -2, sigma= filterdist/dx, mode='constant', truncate= gaustruncate)

    return vort



def calculate_initial_compass_bearing(pointA, pointB):
    import math
    """
    From: https://gist.github.com/jeromer/2005586
    Calculates the bearing between two points.
    The formulae used is the following:
        θ = atan2(sin(Δlong).cos(lat2),
                  cos(lat1).sin(lat2) − sin(lat1).cos(lat2).cos(Δlong))
    :Parameters:
      - `pointA: The tuple representing the latitude/longitude for the
        first point. Latitude and longitude must be in decimal degrees
      - `pointB: The tuple representing the latitude/longitude for the
        second point. Latitude and longitude must be in decimal degrees
    :Returns:
      The bearing in degrees
    :Returns Type:
      float
    """
    if (type(pointA) != tuple) or (type(pointB) != tuple):
        raise TypeError("Only tuples are supported as arguments")

    lat1 = math.radians(pointA[0])
    lat2 = math.radians(pointB[0])

    diffLong = math.radians(pointB[1] - pointA[1])

    x = math.sin(diffLong) * math.cos(lat2)
    y = math.cos(lat1) * math.sin(lat2) - (math.sin(lat1)
            * math.cos(lat2) * math.cos(diffLong))

    initial_bearing = math.atan2(x, y)

    # Now we have the initial bearing but math.atan2 return values
    # from -180° to + 180° which is not what we want for a compass bearing
    # The solution is to normalize the initial bearing as shown below
    initial_bearing = math.degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing