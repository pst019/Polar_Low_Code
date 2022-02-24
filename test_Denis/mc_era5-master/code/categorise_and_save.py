#!/usr/bin/env python3
"""Categorise selected runs and serialise TrackRun object for later use."""
import argparse
import sys
from pathlib import Path
from textwrap import dedent

from loguru import logger

import numpy as np

from octant.core import TrackRun
from octant.decor import get_pbar
from octant.misc import add_domain_bounds_to_mask, check_by_arr_thresh, check_by_mask

import xarray as xr

from common_defs import bbox, columns, period, winters, SMOOTH_FUNC, SMOOTH_KW
import mypaths

# Select runs
# runs2process = dict(era5=[0])  # , interim=[100, 106])
SCRIPT = Path(__file__).name
lsm_paths = {"era5": mypaths.era5_dir / "lsm.nc", "interim": mypaths.interim_dir / "lsm.nc"}


def parse_args(args=None):
    """Parse command line arguments."""
    epilog = dedent(
        f"""Example of use:
    ./{SCRIPT} -n era5 --runs 0,3,10,11 -ll 10,20,65,85
    """
    )
    ap = argparse.ArgumentParser(
        SCRIPT,
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=epilog,
    )

    ap.add_argument(
        "-n",
        "--name",
        type=str,
        required=True,
        choices=["era5", "interim"],
        help="Name of the dataset",
    )
    ap.add_argument(
        "-r",
        "--runs",
        type=str,
        required=True,
        help="Run numbers, comma- or dash-separated (inclusive range)",
    )
    ap.add_argument(
        "--betterlandmask",
        action="store_true",
        help=("Use smoothed ERA5 land mask for both reanalyses"),
    )

    ag_sub = ap.add_argument_group(title="Subset")
    ag_sub.add_argument(
        "-ll",
        "--lonlat",
        type=str,
        default=",".join([str(i) for i in bbox]),
        help=("Lon-lat bounding box (lon0,lon1,lat0,lat1)"),
    )

    ag_etc = ap.add_argument_group(title="Other")
    ag_etc.add_argument(
        "--progressbar", action="store_true", help=("Show progress bar if available")
    )

    return ap.parse_args(args)


def get_lsm(path_to_file, bbox=None, shift=False):
    """Load land-sea mask from a file and crop a region defined by `bbox`."""
    # Load land-sea mask
    lsm = xr.open_dataarray(path_to_file).squeeze()
    if shift:
        # Shift longitude coordinate from (0, 359) to (-180, 179)
        lsm = lsm.assign_coords(longitude=(((lsm.longitude + 180) % 360) - 180)).sortby("longitude")
    if bbox is not None:
        # Select a subset of longitudes and latitudes to speed up categorisation
        lsm = lsm.sel(
            longitude=(lsm.longitude >= bbox[0]) & (lsm.longitude <= bbox[1]),
            latitude=(lsm.latitude >= bbox[2]) & (lsm.latitude <= bbox[3]),
        )
    # Set all non-zero values to 1 to make it the mask binary (0 or 1)
    # lsm = lsm.where(lsm == 0, 1)

    lsm.attrs["units"] = 1

    return lsm


def main(args=None):
    """Loop over track runs and categorise them according `cat_kw`."""
    args = parse_args(args)

    if args.progressbar:
        from octant import RUNTIME

        RUNTIME.enable_progress_bar = True
    pbar = get_pbar()

    if "-" in args.runs:
        _start, _end = args.runs.split("-")
        runs = [*range(int(_start), int(_end) + 1)]
    else:  # if ',' in args.runs:
        runs = [int(i) for i in args.runs.split(",")]
    runs2process = {args.name: runs}

    outer_box = [int(i) for i in args.lonlat.split(",")]
    inner_box = [outer_box[0] + 1, outer_box[1] - 1, outer_box[2] + 1, outer_box[3] - 1]
    gen_box = [outer_box[0] + 1, outer_box[1] - 1, outer_box[2] + 1, outer_box[3] - 3]

    if args.betterlandmask:
        mask = get_lsm(lsm_paths["era5"], bbox=outer_box, shift=True)
        mask = xr.apply_ufunc(SMOOTH_FUNC, mask, kwargs=SMOOTH_KW)
        mask = add_domain_bounds_to_mask(mask, inner_box)
        # Additional constraint on genesis over sea ice covered area
        gen_mask = add_domain_bounds_to_mask(mask, gen_box)
    else:
        mask = get_lsm(lsm_paths[args.name], bbox=outer_box, shift=True)
    lon2d, lat2d = np.meshgrid(mask.longitude, mask.latitude)

    # Additional arguments for land-mask function
    # mask_func_kw = dict(lsm=mask, lmask_thresh=0.5, dist=100.0)
    conditions = [
        (
            "pmc",
            [
                # Mesoscale
                lambda ot: ((ot.vortex_type != 0).sum() / ot.shape[0] < 0.2),
                # Sensible speed
                lambda ot: (ot.average_speed / 3.6) <= 30,
                # Non-stationary
                lambda ot: ot.total_dist_km >= 100.0 and ot.lifetime_h >= 3.0,
                # Far from orography and domain boundaries
                lambda ot: check_by_mask(
                    ot,
                    None,  # Do not pass TrackRun because domain boundaries are already in `mask`
                    mask,
                    lmask_thresh=0.2,
                    time_frac=0.2,
                    dist=60.0,
                    check_domain_bounds=False,
                ),
                # Maritime genesis
                lambda ot: check_by_arr_thresh(
                    ot.xs(0, level="row_idx"), arr=gen_mask, arr_thresh=0.2, oper="le", dist=120.0
                ),
            ],
        )
    ]

    for dset, run_nums in pbar(runs2process.items()):  # , desc="dset"):
        for run_num in pbar(run_nums):  # , leave=False, desc="run_num"):
            logger.info(run_num)
            full_tr = TrackRun()
            for winter in pbar(winters):  # , desc="winter", leave=False):
                logger.info(f"winter: {winter}")
                track_res_dir = mypaths.trackresdir / dset / f"run{run_num:03d}" / winter
                _tr = TrackRun(track_res_dir, columns=columns)
                logger.debug(f"TrackRun size: {len(_tr)}")
                if len(_tr) > 0:
                    logger.info("Begin classification")
                    _tr.classify(conditions, True)
                full_tr += _tr

            full_tr.to_archive(mypaths.procdir / f"{dset}_run{run_num:03d}_{period}.h5")


if __name__ == "__main__":
    sys.exit(main())
