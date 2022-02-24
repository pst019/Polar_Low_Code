# -*- coding: utf-8 -*-
"""Functions for loading STARS and ACCACIA datasets of PMCs."""
from octant.core import OctantTrack

import pandas as pd

import mypaths


def read_stars_file(fname=mypaths.starsdir / "PolarLow_tracks_North_2002_2011"):
    """Read data into a `pandas.DataFrame` from the standard file."""

    def _date_parser(*x):
        return pd.datetime.strptime(" ".join(x), "%Y %m %d %H %M")

    dtype_tuple = (int,) + 5 * (str,) + 4 * (float,)
    dtypes = {k: v for k, v in enumerate(dtype_tuple)}

    df = pd.read_csv(
        fname,
        dtype=dtypes,
        sep=r"\s+",
        skiprows=5,
        date_parser=_date_parser,
        parse_dates={"time": [1, 2, 3, 4, 5]},
    )
    return df


def read_all_stars():
    """Read both North and South subsets of STARS."""
    df_n = read_stars_file(fname=mypaths.starsdir / "PolarLow_tracks_North_2002_2011")
    df_s = read_stars_file(fname=mypaths.starsdir / "PolarLow_tracks_South_2002_2011")

    df_s.N += df_n.N.values[-1]

    return df_n.append(df_s).reset_index(drop=True)


def read_all_accacia():
    """Load ACCACIA tracks as `pandas.DataFrame`"""

    def _date_parser(x):
        return pd.datetime.strptime(x, "%Y%m%d%H%M")

    df = pd.read_csv(
        mypaths.acctracks,
        delimiter="\t",
        names=["N", "time", "lon", "lat"],
        parse_dates=["time"],
        date_parser=_date_parser,
    )
    return df


def prepare_tracks(obs_df, filter_funcs=[], system_id_label="N"):
    """Make a list of those tracks that satisfy the list of conditions."""

    selected = []
    for i, df in obs_df.groupby(system_id_label):
        ot = OctantTrack.from_df(df)

        flag = True
        for func in filter_funcs:
            flag &= func(ot)

        if flag:
            selected.append(ot)
    return selected
