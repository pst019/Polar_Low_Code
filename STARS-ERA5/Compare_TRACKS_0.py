#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 16:29:52 2019

@author: pst019
"""


#import pandas as pd
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from octant.core import TrackRun, OctantTrack
from f_carto import *


from pathlib import Path


track_dir= "/media/pst019/1692A00D929FEF8B/ERA5_STARS/tracks"
t_dir= Path(".") / track_dir
#t_dir= Path(".") / "/media/pst019/1692A00D929FEF8B/ERA5_STARS/pmctrack-master/results/test6"

#tr= TrackRun(t_dir)
tr = TrackRun(t_dir, columns=['lon', 'lat', 'vo', 'time', 'area', 'vortex_type', 'slp'])
print(tr)
print(tr.columns)


fig = plt.figure()



def lcc_map(fig, subplot_grd=111, clon=None, clat=None, coast=None, extent=None):
    """Create cartopy axes in Lambert Conformal projection."""
    proj = ccrs.LambertConformal(central_longitude=clon, central_latitude=clat)

    # Draw a set of axes with coastlines:
    ax = fig.add_subplot(subplot_grd, projection=proj)
    if isinstance(extent, list):
        ax.set_extent(extent, crs=ccrs.PlateCarree())
    if isinstance(coast, dict):
        feature = cartopy.feature.NaturalEarthFeature(
            name="coastline", category="physical", **coast
        )
        ax.add_feature(feature)
    return ax




plt.rcParams["figure.figsize"] = (15, 10)  # change the default figure size
COAST = dict(scale="50m", alpha=0.5, edgecolor="#333333", facecolor="#AAAAAA")
clon = 10
clat = 75
extent = [-20, 50, 65, 85]
LCC_KW = dict(clon=clon, clat=clat, coast=COAST, extent=extent)

mapkw = dict(transform=ccrs.PlateCarree())





conditions = [
    ("long_lived", [lambda ot: ot.lifetime_h >= 6]),
    (
        "far_travelled_and_very_long_lived",
        [lambda ot: ot.lifetime_h >= 36, lambda ot: ot.gen_lys_dist_km > 300.0],
    ),
    ("strong", [lambda x: x.max_vort > 1e-3]),
]

tr.classify(conditions)


single= tr["long_lived"]
print(single.gen_lys_dist_km)


#ax= Plot_Polar_Stereo(fig, central_longitude= 30, extent= [-10, 50, 60, 80], subplot= (1,1,1))
ax= Plot_PlateCarree(fig, central_longitude= 0, extent= [-45, 45, 50, 80], subplot= (1,1,1))


#ax = lcc_map(fig, **LCC_KW)
single.plot_track(ax=ax, color=None)




"""file from Denis"""
#filename = "era5_run000_2000_2018__matched_to_stars_6h__bs2000_beta100.h5"
#my_track_run = TrackRun.from_archive(filename)
#print(my_track_run.data)