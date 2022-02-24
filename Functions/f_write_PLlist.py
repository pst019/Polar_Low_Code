#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 10:36:55 2018

@author: pst019

This was before in f_imp_ERA(2)
It is writing PLlists in the TRACK2PL code
"""

        
def WritePLlistComp2S(PLlist, year, hemi='', name='', model='ERA', rowname= None):
    """write the PLlist into an excel file for each year - used by TRACK2PL_7 and ASR_TRACK2PL_3 (the newest versions)
    complete PLlist (year, month, TrackPLnr, point in time 6 hourly starting with 0 of month (ERA-reanalysis - timepoint),
    lon, lat, filtered vorticity, SST- theta_e500, maximum tropopause wind north of point (UPVnorth), PLpoint(yes= 1, no= 0))
    including shear"""
    
    import csv
    harddisk= Mediadir +'/PL/PL_file_'+model+'/'
    if hemi== '_SH_' and model=='ERA':
        with open(harddisk+'PL_file_complete'+hemi+name+str(year)+'.csv', 'w') as csvfile:
            if rowname== None: rowname= ['year', 'month', 'TPLnr', 'time pnt', 'lon', 'lat', 'vort_filt', 'thetaSST-theta500', 'UPVsouth', 'PLpoint']
            writer = csv.writer(csvfile, delimiter="\t")
            writer.writerow(rowname)
            [writer.writerow(r) for r in PLlist.T]
    elif model=='ERA':
        with open(harddisk+'PL_file_complete'+name+str(year)+'.csv', 'w') as csvfile:
            if rowname== None: rowname= ['year', 'month', 'TPLnr', 'time pnt', 'lon', 'lat', 'vort_filt', 'thetaSST-theta500', 'UPVnorth', 'shear', 'PLpoint']
            writer = csv.writer(csvfile, delimiter="\t")
            writer.writerow(rowname)
            [writer.writerow(r) for r in PLlist.T]
    elif model=='ASR':
        with open(harddisk+'PL_file_complete'+name+str(year)+'.csv', 'w') as csvfile:
            if rowname== None: rowname= ['year', 'month', 'TPLnr', 'time pnt (3hourly)', 'lon', 'lat', 'vort_filt', 'U10', 'thetaSST-theta500', 'U500north', 'PLpoint']
            writer = csv.writer(csvfile, delimiter="\t")
            writer.writerow(rowname)
            [writer.writerow(r) for r in PLlist.T]
