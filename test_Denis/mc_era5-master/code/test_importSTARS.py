#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 10:22:36 2019

@author: pst019
"""

from datetime import datetime
import json
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import string


from common_defs import (
    aliases,
    bbox,
    conf_key_typeset,
    datasets,
    dset_ctrl_tstep_h,
    nyr,
    winters,
    winter_dates,
    _runs_grid_formatter,
)
from plot_utils import cc, use_style
from match_to_ref import match_options, _make_match_label, REF_DATASETS
from obs_tracks_api import prepare_tracks, read_all_accacia, read_all_stars
import mypaths


REF_SET = "stars"

STARS= read_all_stars()
stars_tracks = prepare_tracks(STARS) #, REF_DATASETS[REF_SET]["filter_func"])