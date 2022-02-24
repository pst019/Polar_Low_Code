# -*- coding: utf-8 -*-
"""Paths to data."""
from os import getenv
from pathlib import Path


# Top-level directory containing code and data (one level up)
topdir = Path(__file__).absolute().parent.parent

# Output directories
plotdir = topdir / "figures"

# Local data
accdir = topdir / "data" / "tracks" / "accacia"
acctracks = accdir / "pmc_loc_time_ch4_20Mar-02Apr.txt"
starsdir = topdir / "data" / "tracks" / "stars"
starstracks = starsdir / "PolarLow_tracks_North_2002_2011"

# Satellite data
ascat_file = topdir / "data"/ "ascat" / "ascat_20130326_123601_metopa_33384_eps_o_125_2101_ovw.l2.nc"
avhrr_file = topdir / "data"/ "avhrr" / "201303261220.ch4.r8.nps.tif"

# External data
if getenv("TRUSTZONE") is not None:
    # if on MONSooN
    datadir = Path(getenv("DATADIR"))
    trackresdir = datadir / "pmctrack" / "output"
    procdir = datadir / "pmctrack" / "processed_data"
else:
    datadir = Path.home() / "phd"
    trackresdir = datadir / "pmc_tracking" / "results"
    procdir = datadir / "pmc_tracking" / "results" / "processed_data"
runsgridfile = trackresdir / "runs_grid.json"

# Reanalyses
ra_dir = datadir / "reanalysis"
era5_dir = ra_dir / "era5"
interim_dir = ra_dir / "interim"
