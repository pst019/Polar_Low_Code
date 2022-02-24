#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 18:04:12 2018

@author: pst019
retrieve data from ECMWF
see examples: https://software.ecmwf.int/wiki/display/WEBAPI/Python+ERA-interim+examples
"""

from ecmwfapi import ECMWFDataServer
from ecmwfapi import ECMWFService

    
mediadir= "/media/pst019/1692A00D929FEF8B/ECMWF/"

leveltype = "pl" #(pl, sfc)
area= "85/-30/50/60"

year= 2008
month= 3
day= 4

#var, param = 'T', "130.128"
#var, param = 'U', "131/132"
#var, param = 'SH', "133.128"
#var, param, leveltype = 'SLP', "151.128", "sfc"
var, param, leveltype= 'OLW', "179.128", "sfc"  #top of atmosphere outgoing longwave radiation - as pseudo sat image

dataset='ERA5' #(HRES, ERA5, ERA-I)
#dataset='HRES'

"""for ERA5"""
if dataset == 'ERA5':
    #plevels= "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000"
    
    plevels = "200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000"
#    times = "00:00:00/03:00:00/06:00:00/09:00:00/12:00:00/15:00:00/18:00:00/21:00:00"
    times = "00/01/02/03/04/05/06/07/08/09/10/11/12/13/14/15/16/17/18/19/20/21/22/23" #"00/03/06/09/12/15/18/21" #"00:00:00/01:00:00/02:00:00"
#    times = "0/to/24/by/3"
    
    server = ECMWFDataServer()
    server.retrieve({
        "class"     : "ea",
        "stream"    : "oper",        
        "dataset"   : "era5",
        "type"      : "an",
        "expver"    : "1",
        
        "levtype"   : leveltype,
        "levelist"  : plevels,        
        
        "grid"      : "0.25/0.25",
        "area"      : area,        
               
        "param"     : param,
        
        "date"      :  str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2) ,
        "time"      : times,
        "format"    : "netcdf",    
        "target"    : mediadir +dataset +"/"+ str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+var+"_2.nc",
    })


"""for ERA-I"""
if dataset == 'ERA-I':
    server = ECMWFDataServer()
    server.retrieve({
        'class'     : "ei",            
        'stream'    : "oper",
        'dataset'   : "interim",
        'type'      : "an",
        
        'levtype'   : leveltype,
#        "levelist": plevels,        
        'grid'      : "0.75/0.75",
        "area"      : area,        
        
        'param'     : param,

        'date'      : str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2),
        'time'      : "00/06/12/18",

        'step'      : "0",
        'format'    : "netcdf",
        'target'    : mediadir + dataset +"/"+ str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+var+'.nc',
     })    
    




    
"""for ECMWF operational"""
if dataset == 'HRES':
    """to obtain analysis"""
#    plevels= "200/250/300/400/500/600/700/800/850/900/925/950/1000"
#    times = "00:00:00/06:00:00/12:00:00/18:00:00"
#    
#    server.retrieve({
#        "class": "od",
#        "dataset": "Oper",
#        "date":  str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2) ,               
#        "expver": "1",
#        "levtype": leveltype,
#        "levelist": plevels,
#        "param": param,
#        "stream": "oper",
#        "time": times,
#        "type": "an",
#        "format": "netcdf",    
#        "target": mediadir+ dataset +"/"+str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+var+'.nc',
#    })

    """to obtain the surface situation"""
#    server = ECMWFService("mars")
#    server.execute(
#        {
#        "class": "od",
#        "date": str(year)+str(month).zfill(2)+str(day).zfill(2),
#        "expver": "1",
#        "levtype": "sfc",
#        "param": "165.128/166.128/167.128/168.128/228.128",
#        "step": "0/to/3/by/1",
#        "stream": "oper",
#        "time": "12",
#        "type": "fc", 
#        "grid": "0.125/0.125",
#        "area": "60/-10/50/10"
#        },
#        mediadir+ dataset +"/"+str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+var+'.grib')

    """to obtain vertical profiles for a small area"""
    plevels= "300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000"
    times = "12:00:00/18:00:00"
    
    server.retrieve({
        "class": "od",
        "dataset": "Oper",
        "date":  str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2) ,               
        "expver": "1",
        "levtype": leveltype,
        "levelist": plevels,
        "param": "165.128/166.128/167.128/168.128/228.128",
        "stream": "oper",
        "time": times,
        "type": "fc",
        "format": "netcdf", 
        "grid": "0.125/0.125",
        "area": "76/-7/68/10",        
        "target": mediadir+ dataset +"/"+str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_vertical.nc',
    })
    