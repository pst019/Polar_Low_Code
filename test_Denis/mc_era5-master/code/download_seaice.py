#!/usr/bin/env python3
import cdsapi

import mypaths

AREA = "90/-30/60/60"
TOPDIR = mypaths.ra_dir
PRODUCT_NAME = "era5"


if __name__ == "__main__":
    outdir = TOPDIR / PRODUCT_NAME
    outdir.mkdir(exist_ok=True)
    fname = f"{PRODUCT_NAME}.an.sfc.2000-2018.sea_ice_cover.nc"
    target = str(outdir / fname)

    c = cdsapi.Client()

    c.retrieve(
        "reanalysis-era5-single-levels",
        {
            "variable": "sea_ice_cover",
            "product_type": "reanalysis",
            "year": [
                "2000",
                "2001",
                "2002",
                "2003",
                "2004",
                "2005",
                "2006",
                "2007",
                "2008",
                "2009",
                "2010",
                "2011",
                "2012",
                "2013",
                "2014",
                "2015",
                "2016",
                "2017",
                "2018",
            ],
            "month": ["01", "02", "03", "04", "10", "11", "12"],
            "day": [
                "01",
                "02",
                "03",
                "04",
                "05",
                "06",
                "07",
                "08",
                "09",
                "10",
                "11",
                "12",
                "13",
                "14",
                "15",
                "16",
                "17",
                "18",
                "19",
                "20",
                "21",
                "22",
                "23",
                "24",
                "25",
                "26",
                "27",
                "28",
                "29",
                "30",
                "31",
            ],
            "time": "12:00",
            "area": AREA,
            "format": "netcdf",
        },
        fname,
    )
