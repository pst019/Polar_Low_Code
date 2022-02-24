#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 18:04:12 2018

@author: pst019
see examples: https://software.ecmwf.int/wiki/display/WEBAPI/Python+ERA-interim+examples
"""

from ecmwfapi import ECMWFDataServer
    
server = ECMWFDataServer()
    
mediadir= "/media/pst019/1692A00D929FEF8B/ECMWF/"

leveltype = "pl"
times = "00:00:00" #/03:00:00/06:00:00/09:00:00/12:00:00/15:00:00/18:00:00/21:00:00"
year= 2018
month= 2
day= 23

var, param = 'T', "130.128"
#var, param = 'U', "131/132"
#var, param = 'SH', "133.128"


"""for ECMWF operational"""

plevels= "200/250/300/400/500/600/700/800/850/900/925/950/1000"
times = "00:00:00" #/06:00:00/12:00:00/18:00:00"

#server.retrieve({
#    "class": "od",
#    "dataset": "Oper",
#    "date":  str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2) ,               
#    "expver": "1",
#    "levtype": leveltype,
#    "levelist": plevels,
#    "param": param,
#    "stream": "oper",
#    "time": times,
#    "type": "an",
#    "format": "netcdf",    
#    "target": mediadir + "HRES/"+str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+var+'.nc',
#})

server.retrieve({
    "class": "od",
    "date":  str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2) ,               
    "expver": "1",
    "levtype": "sfc",
    "param": 167.128,
    "step":0/1,
    "stream": "oper",
    "time": times,
    "type": "fc",
    "format": "netcdf",    
    "target": mediadir + "HRES/"+str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_'+var+'.nc',
})

