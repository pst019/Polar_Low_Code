# coding: utf-8
"""Match cyclone tracks from ERA5 and ERA-Interim to a reference list of polar lows."""
import json
from loguru import logger as L
from pathlib import Path

import octant
from octant.core import TrackRun
from octant.decor import get_pbar

from common_defs import CAT, bbox, datasets, period, winters
import mypaths
from obs_tracks_api import read_all_accacia, read_all_stars, prepare_tracks


NAME = "stars"
RUN_GROUP = "vort_thresh"

RUN_GROUPS = {
    "vort_thresh": {
        "start": 0,
        "paths": {
            "era5": mypaths.procdir / "runs_grid_vort_thresh_era5.json",
            "interim": mypaths.procdir / "runs_grid_vort_thresh_interim.json",
        },
    },
    "tfreq": {
        "start": 100,
        "paths": {
            "era5": mypaths.procdir / "runs_grid_tfreq_era5.json",
            "interim": mypaths.procdir / "runs_grid_tfreq_interim.json",
        },
    },
    "extra": {
        "start": 200,
        "paths": {"interim": mypaths.procdir / "runs_grid_low_thresh_interim.json"},
    },
}

match_options = [
    # dict(method="bs2000", beta=25.0),
    dict(method="bs2000", beta=50.0),
    # dict(method="bs2000", beta=75.0),
    dict(method="bs2000", beta=100.0),
    # dict(method="simple", thresh_dist=150.0),
    # dict(method="simple", thresh_dist=200.0),
    # dict(method="simple", thresh_dist=250.0),
    # dict(method="simple", thresh_dist=300.0),
]

REF_DATASETS = {
    "stars": {
        "time_dict": {
            k: (f"{k.split('_')[0]}-10-01", f"{k.split('_')[1]}-04-30") for k in winters[1:11]
        },
        "load_func": read_all_stars,
        "filter_func": [lambda ot: ot.within_rectangle(*bbox) and ot.lifetime_h > 0],
    },
    "accacia": {
        "time_dict": {"accacia": ("2013-03-15", "2013-04-05")},
        "load_func": read_all_accacia,
        "filter_func": [lambda ot: ot.within_rectangle(*bbox) and ot.lifetime_h > 0],
    },
}


def _make_match_label(match_kwargs, delim="_"):
    match_kwargs_label = []
    for k, v in match_kwargs.items():
        try:
            vv = int(v)
        except ValueError:
            vv = v
        match_kwargs_label.append(f"{k}={vv}")
    return delim.join(match_kwargs_label)


@L.catch
def main():
    LOGPATH = Path(__file__).parent / "logs"
    LOGPATH.mkdir(exist_ok=True)
    # L.remove(0)
    L.add(LOGPATH / f"log_match_to_{NAME}_{{time}}.log")
    octant.RUNTIME.enable_progress_bar = True
    pbar = get_pbar(use="tqdm")
    octant.RUNTIME.enable_progress_bar = False

    obs_tracks = prepare_tracks(
        REF_DATASETS[NAME]["load_func"](), filter_funcs=REF_DATASETS[NAME]["filter_func"]
    )
    n_ref = len(obs_tracks)
    L.debug(f"Number of suitable tracks: {n_ref}")

    # Define an output directory and create it if it doesn't exist
    output_dir = mypaths.procdir / "matches"
    output_dir.mkdir(exist_ok=True)

    # Loop over datasets, runs, subsets, matching methods
    for dset in pbar(datasets[:]):  # , desc="dataset"):
        with RUN_GROUPS[RUN_GROUP]["paths"][dset].open("r") as fp:
            runs_grid = json.load(fp)
        for run_id, _ in pbar(enumerate(runs_grid, RUN_GROUPS[RUN_GROUP]["start"])):
            TR = TrackRun.from_archive(mypaths.procdir / f"{dset}_run{run_id:03d}_{period}.h5")
            L.debug(mypaths.procdir / f"{dset}_run{run_id:03d}_{period}.h5")
            L.debug(TR)
            for match_kwargs in pbar(match_options):  # , desc="match options"):
                match_pairs_abs = []
                for winter, w_dates in pbar(REF_DATASETS[NAME]["time_dict"].items()):
                    tr = TR.time_slice(*w_dates)
                    L.debug(tr)
                    L.debug(run_id)
                    L.debug(match_kwargs)
                    L.debug(winter)
                    match_pairs = tr.match_tracks(obs_tracks, subset=CAT, **match_kwargs)
                    for match_pair in match_pairs:
                        match_pairs_abs.append(
                            (match_pair[0], obs_tracks[match_pair[1]].N.unique()[0])
                        )
                match_kwargs_label = _make_match_label(match_kwargs)

                # Save matching pairs to a text file
                fname = f"{dset}_run{run_id:03d}_{period}_{NAME}_{match_kwargs_label}.txt"
                with (output_dir / fname).open("w") as fout:
                    fout.write(
                        f"""# {dset}
# {run_id:03d}
# {period}
# {match_kwargs_label}
"""
                    )
                    for match_pair in match_pairs_abs:
                        fout.write("{:d},{:d}\n".format(*match_pair))


if __name__ == "__main__":
    main()
