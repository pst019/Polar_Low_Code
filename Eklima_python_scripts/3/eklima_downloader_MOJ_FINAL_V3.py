#!/usr/bin/python

import numpy as np
import eklima_MOJ_FINAL_V3 as e

from datetime import datetime as dt
from copy import copy

# ['rr_1','dd','ta','sa']
# elems = ['sa'] #['ff', 'fg_1', 'fm']
# elems = ['ta','ff','po']
elems = ['ta','ff','po','uu','dd','rr_1','rr_12','tax','tan','tss','sa']
years = np.arange(2017,2018)

#filelist=set(['99735']) #,'99740','99752','99927','99935','99938','99710','99720','99754','99760','99765','99790','99820','99840','99880','99910','99928','99870']);
filelist=set(['99870','99840']);

test={}

for year in years:
    supported = 0
    for elem in elems:
        
        stations = e._get_stations_from_tstype_elems('hourly', elem)
        supported = set(stations.keys())
        supported_stations = stations
        supported = supported.intersection(set(stations.keys()))
        
        supported2=set(supported).intersection(filelist)
        print(supported2)
        
        for stnr in list(supported_stations.keys()): 
            if stnr not in supported2:
                del supported_stations[stnr]
        
        test[elem]=supported_stations
        
        print(test[elem])
        
        # Downloading matching data	
        e.download_hourly_by_elems(elem, dt(year,1,1), dt(year,12,31,23,59,59), stations=copy(test[elem]), savemat=True)
