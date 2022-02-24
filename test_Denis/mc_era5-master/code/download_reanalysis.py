#!/usr/bin/env python3
"""Download ERA5 data using CDS API."""
import argparse
from calendar import monthrange
from cdsapi import Client
from datetime import datetime
from loguru import logger as L
from pathlib import Path
import sys
import xarray as xr

import mypaths


# Top-level directory
# The files are downloaded to TOPDIR / ['era5', 'interim']
PRODUCT_NAME = "era5"
PRODUCT_TYPE = "reanalysis"
TOPDIR = mypaths.ra_dir
LEVELS = ["500", "600", "700", "800", "850", "900", "925", "950", "975", "1000"]
# VARNAMES = ['msl', 'u/v/vo']
LEVTYPES = {
    "sst": "single",
    "ci": "single",
    "msl": "single",
    "lsm": "single",
    "u": "pressure",
    "v": "pressure",
    "vo": "pressure",
}
LEV_ABBR = {"single": "sfc", "pressure": "pl"}
# Subset (clip) to an area. Specify as N/W/S/E in Geographic lat/long degrees.
AREA = "90/-30/60/60"
NAME = "reanalysis-{product_name}-{levtype}-levels"
FORMAT = "netcdf"

SCRIPT = Path(Path(__file__).name).with_suffix("")
# Logging set up
LOGPATH = Path(__file__).parent / "logs"
LOGPATH.mkdir(exist_ok=True)
LOGFILE = LOGPATH / "{}_{:%Y%m%d%H%M}.log".format(SCRIPT, datetime.now())
L.add(LOGFILE)
NC_SAVE_KW = dict(
    encoding={"time": {"units": "hours since 1900-01-01 00:00:0.0"}}, unlimited_dims=("time",)
)


def download(client, req, varnames):
    """Retrieve data from CDS server."""
    outdir = TOPDIR / req["dataset"]
    outdir.mkdir(exist_ok=True)
    dt0, dt1 = req["date"].split("/to/")
    dt0, dt1 = [datetime.strptime(i, "%Y-%m-%d") for i in [dt0, dt1]]
    fnames = []
    for varname in varnames:
        req["levtype"] = LEVTYPES[varname]
        if req["levtype"] == "pl":
            req["levelist"] = LEVELS
        fname = f"{req['dataset']}.{req['type']}.{req['levtype']}" f".{dt0:%Y%m%d}-{dt1:%Y%m%d}.nc"
        req["target"] = str(outdir / fname)
        req["param"] = varname
        client.retrieve(req)
        fnames.append(req["target"])
    return fnames


def slice_files(fnames, req):
    ds = xr.open_mfdataset(fnames)
    outdir = Path(fnames[0]).parent
    for varname in ds.data_vars:
        for m in set(ds.time.dt.month.values):
            month_subset = ds[varname].loc[ds.time.dt.month == m]
            month_subset = month_subset.sortby("time")
            fname_out = (
                f"{req['dataset']}.{req['type']}.{LEVTYPES[varname]}"
                f".{month_subset.time.dt.year.values[0]}"
                f".{m:02d}.{varname}.nc"
            )
            L.info(fname_out)
            month_subset.to_netcdf(str(outdir / fname_out), **NC_SAVE_KW)


def parse_args(args=None):
    example = """
./get_era5_data.py -vvv -y 2010-2012 -m 1,2,3,4,10,11,12 -t 1 -n msl,u,v,vo
"""
    ap = argparse.ArgumentParser(
        SCRIPT,
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=f"Example of use:\n{example}",
    )

    ap.add_argument(
        "-q", "--quiet", action="store_true", help="Suppress log messages (write to log file only)"
    )
    ap.add_argument(
        "-y", "--years", required=True, type=str, help="Years list (sep: ,) or span (sep: -)"
    )
    ap.add_argument(
        "-m", "--months", required=True, type=str, help="Months list (sep: ,) or span (sep: -)"
    )
    ap.add_argument(
        "-n", "--names", required=True, type=str, help="Comma-separated list of variables"
    )
    ap.add_argument(
        "-t", "--timestep", type=int, required=True, choices=[1, 3, 6], help="Time step (h)"
    )
    ap.add_argument(
        "--hours",
        type=str,
        required=False,
        choices=["all", "midday"],
        default="all",
        help="All times of the day or 12:00 (e.g. for sea ice)",
    )

    ap.add_argument(
        "--slice",
        action="store_true",
        default=False,
        help="Slice the downloaded files month- and variable-wise",
    )

    return ap.parse_args(args)


def main(args=None):
    """Main entry point of the script."""
    args = parse_args(args)
    if args.quiet:
        # write messages to log file only
        L.remove(0)

    # fnames = []
    if "-" in args.years:
        years = args.years.split("-")
        years = list(range(int(years[0]), int(years[1]) + 1))
    elif "," in args.years:
        years = [int(i) for i in args.years.split(",")]
    else:
        try:
            years = [int(args.years)]
        except TypeError:
            raise argparse.ArgumentTypeError("Wrong format of --years argument")
    L.info(f"years = {years}")

    if args.months:
        if "-" in args.months:
            months = args.months.split("-")
            months = list(range(int(months[0]), int(months[1]) + 1))
        elif "," in args.months:
            months = [int(i) for i in args.months.split(",")]
        else:
            try:
                months = [int(args.months)]
            except TypeError:
                raise argparse.ArgumentTypeError("Wrong format of --months argument")
    else:
        months = list(range(1, 13))
    L.info(f"months = {months}")

    if "," in args.names:
        varnames = args.names.split(",")
    else:
        varnames = [args.names]
    L.info(f"varnames = {varnames}")

    tstep = args.timestep
    if args.hours == "all":
        times = "/".join(["{:02d}:00".format(i) for i in range(0, 24, tstep)])
    elif args.hours == "midday":
        times = "12:00:00"

    base_req = {
        "product_type": PRODUCT_TYPE,
        "time": times,
        "area": AREA,
        # Optional: output in NetCDF format.
        # If you want data in GRIB format, remove this line.
        "format": FORMAT,
    }

    C = Client()
    outdir = TOPDIR / PRODUCT_NAME
    outdir.mkdir(exist_ok=True)

    for year in years:
        for month in months:
            for varname in varnames:
                days = [f"{i:02d}" for i in range(1, monthrange(year, month)[1] + 1)]
                req = base_req.copy()
                req.update(variable=varname, year=str(year), month=f"{month:02d}", day=days)
                levtype = LEVTYPES[varname]
                if levtype == "pressure":
                    req.update(pressure_level=LEVELS)

                L.info(f"Full request:\n{req}")

                fname = f"{PRODUCT_NAME}.an.{LEV_ABBR[levtype]}" f".{year}.{month:02d}.{varname}.nc"
                target = str(outdir / fname)
                L.info(f"target = {target}")
                C.retrieve(NAME.format(product_name=PRODUCT_NAME, levtype=levtype), req, target)

    # fnames += download(C, req, varnames)

    # L.info('Retrieved files:')
    # for fname in fnames:
    #     L.info(fname)

    # if args.slice:
    #     L.info('Slicing data...')
    #     slice_files(fnames, req)


if __name__ == "__main__":
    sys.exit(main())
