# -*- coding: utf-8 -*-
"""
Common objects for PMC climatology code
"""
import calendar

import numpy as np

from scipy.ndimage.filters import gaussian_filter


# Categorisation
CAT = "pmc"

# Smoothing
SMOOTH_FUNC = gaussian_filter
SMOOTH_KW = {"sigma": (1.2, 4.5)}

# Columns of vortrack text files
columns = ["lon", "lat", "vo", "time", "area", "vortex_type", "slp"]

# For legend keys
aliases = {"era5": "ERA5", "interim": "ERA-Interim"}

# Shorthand list
datasets = ["era5", "interim"]

dset_names = (
    ("era5_run000", "ERA5, CTRL"),
    ("interim_run100", "ERA-Interim, CTRL"),
    #     ('interim_run100', 'ERA-Interim, LVT')
)
dset_ctrl_tstep_h = {"era5": 1, "interim": 3}

# Time period parameters
START_YEAR = 2000
nyr = 18
winters = [f"{START_YEAR+i}_{START_YEAR+i+1}" for i in range(nyr)]
period = f"{winters[0][:4]}_{winters[-1][-4:]}"
winter_dates = {k: (f"{k.split('_')[0]}-10-01", f"{k.split('_')[1]}-04-30") for k in winters}

ndays_per_month_total = np.zeros((12), dtype=int)
for yr in range(START_YEAR, START_YEAR + nyr):
    ndays_per_month_total += np.array(calendar.mdays[1:])
    if calendar.isleap(yr + 1):
        # +1 because we start from autumn
        ndays_per_month_total[1] += 1

month_weights = 30 / ndays_per_month_total

# Bounding box of study area
# Set to 1 degree bigger than TrackRun.conf.extent
# because otherwise land mask passed to check_by_mask will have no effect
bbox = [-21, 51, 64, 86]
# Additional bbox, equivalent to TrackRun.conf.extent
inner_bbox = [-20, 50, 65, 85]

# Geographical names to put on map plots
toponyms = [
    dict(name="Svalbard", lon=20, lat=79),
    dict(name="Greenland", lon=-35, lat=80),
    dict(name="Norway", lon=23, lat=70),
    dict(name="Fram\nStrait", lon=0, lat=80),
    dict(name="Greenland\nSea", lon=-5, lat=73),
    dict(name="Barents\nSea", lon=40, lat=75),
    dict(name="Norwegian\nSea", lon=5, lat=66),
    dict(name="Franz\nJosef\nLand", lon=53, lat=80.5),
]

# For axis labels
conf_key_typeset = {
    "zeta_max0": r"$\zeta_{max}$",
    "zeta_min0": r"$\zeta_{min}$",
    "r_steering": "$r_{steer}$",
    "smth_type": "Smoothing type",
    "r_smth": "$r_{smooth}$",
    "del_r": "$r_{link}$",
    "merge_opt": "Merging option",
    "halo_r": "$r_{halo}$",
    "vor_lvl": r"$\zeta_{level}$",
}


def _runs_grid_formatter(run_dict):
    txt = ""
    for k, v in run_dict.items():
        key = conf_key_typeset.get(k, k)
        #         val = misc.unit_format(v)
        #         if val == '':
        #             val = '1'
        if v < 1:
            val = f"{v:5.5f}"
        else:
            val = v
        s = f"{key} = {val}"
        txt += f" {s:}\n"
    if not txt:
        txt = "CTRL"
    else:
        txt = txt.strip("\n")
    return txt
