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
    T given in K (check)
    plev given in hPa"""
    RCp= 2/7
    return T* (1000/plev)**RCp


def EquiPotTemp(T, sHum, plev):
    """calculate the equivalent potential temperature
    all variables can be scalars or arrays (of same dimension)
    sHum given in ERA unit [kg/m**3]? (also called water vapor mixing ratio)"""
    Lc, Cp= 2501, 1.006
    return PotTemp(T, plev)* np.exp(sHum* Lc / (T * Cp))


def RH2SH(RH, p, T):
    """translate from relative humidity (RH) in % to specific humditiy (SH) in ???kg/m**3 g/kg
    at pressure p in hPa and temperature T in K
    rel humidity to specific humidity: http://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
    spechum= relhum[0to1]*100/(26.3 * press[hPa]) * exp(17.67*(T[k]-273.15)/(T[K]-29.65))"""
    return RH/(.263 *p)* np.exp(17.67*(T-273.15) / (T-29.65))



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
