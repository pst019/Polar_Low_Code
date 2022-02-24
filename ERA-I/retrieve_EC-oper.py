#!/usr/bin/env python
from ecmwfapi import ECMWFService

mediadir= "/media/pst019/1692A00D929FEF8B/ECMWF/HRES/"

  
#server = ECMWFService("mars")
#server.execute(
#    {
#    "class": "od",
#    "date": "20161106",
#    "expver": "1",
#    "levtype": "sfc",
#    "param": "165.128/166.128/167.128/168.128/228.128",
#    "step": "0/to/3/by/1",
#    "stream": "oper",
#    "time": "12",
#    "type": "fc", 
#    "grid": "0.125/0.125",
#    "area": "60/-10/50/10"
#    },
#    "EC_oper_2016110612.grib")

server = ECMWFService("mars")
server.execute(
    {
    "class": "od",
    "date": "20080223",
    "expver": "1",
    "levtype": "sfc",
    "param": "165.128/166.128/167.128/172.128",
    "step": "0/to/24/by/1",
    "stream": "oper",
    "time": "0",
    "type": "fc", 
    "grid": "av",
    "area": "82/10/73/20",
    "format":"netcdf",
    },
    mediadir+"20080223.nc")

#Model terrain height	mterh	m	260183