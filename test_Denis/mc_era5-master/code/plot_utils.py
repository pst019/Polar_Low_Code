# -*- coding: utf-8 -*-
"""
Various plotting functions and useful parameters
"""
from itertools import cycle
import string

import cartopy.crs as ccrs
import matplotlib as mpl

import matplotlib.pyplot as plt
import numpy as np


cc = plt.rcParams["axes.prop_cycle"]
icc = cycle(cc)

iletters = iter(string.ascii_lowercase)

# Common plotting settings
CBARKW = dict(orientation="vertical")
AXGR_KW = dict(
    axes_pad=0.2, cbar_location="right", cbar_mode="single", cbar_pad=0.2, cbar_size="3%"
)

# Common cartopy settings
trans = dict(transform=ccrs.PlateCarree())
COAST = dict(scale="50m", alpha=0.5, edgecolor="#333333", facecolor="#AAAAAA")
clon = 15
clat = 75
extent = [-20, 50, 65, 85]
LCC_KW = dict(clon=clon, clat=clat, coast=COAST, extent=extent, ticks=None)

mapkey = dict(transform=ccrs.PlateCarree())

clev101 = list(np.linspace(-1, 1, 9))
clev101.remove(0)
clev101 = np.array(clev101)

# Density map plot settings
abs_plt_kw = dict(cmap="Oranges", extend="max", **trans)


def div_cmap(
    numcolors=11,
    name="bwr_div_cmap",
    mincol="blue",
    midcol="white",
    maxcol="red",
    under=None,
    midcol_alpha=0,
    over=None,
):
    """
    Create a custom diverging colormap with three colors

    Default is blue to transparent white to red with 11 colors.
    Colors can be specified in any way understandable
    by matplotlib.colors.ColorConverter.to_rgb()
    """
    c_max = mpl.colors.colorConverter.to_rgba(maxcol)
    c_min = mpl.colors.colorConverter.to_rgba(mincol)
    c_mid = mpl.colors.colorConverter.to_rgba(midcol, alpha=midcol_alpha)
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        name=name, colors=[c_min, c_mid, c_max], N=numcolors
    )
    if under is not None:
        cmap.set_under(under)
    if over is not None:
        cmap.set_over(over)
    return cmap


def use_style():
    plt.style.use("paperfig.mplstyle")


def heatmap(data, row_labels, col_labels, ax=None, cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Arguments:
        data       : A 2D numpy array of shape (N,M)
        row_labels : A list or array of length N with the labels
                     for the rows
        col_labels : A list or array of length M with the labels
                     for the columns
    Optional arguments:
        ax         : A matplotlib.axes.Axes instance to which the heatmap
                     is plotted. If not provided, use current axes or
                     create a new one.
        cbar_kw    : A dictionary with arguments to
                     :meth:`matplotlib.Figure.colorbar`.
        cbarlabel  : The label for the colorbar
    All other arguments are directly passed on to the imshow call.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-45, ha="right", rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="w", linestyle="-", linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(
    im, data=None, valfmt="{x:.2f}", textcolors=["black", "white"], threshold=None, **textkw
):
    """
    A function to annotate a heatmap.

    Arguments:
        im         : The AxesImage to be labeled.
    Optional arguments:
        data       : Data used to annotate. If None, the image's data is used.
        valfmt     : The format of the annotations inside the heatmap.
                     This should either use the string format method, e.g.
                     "$ {x:.2f}", or be a :class:`matplotlib.ticker.Formatter`.
        textcolors : A list or array of two color specifications. The first is
                     used for values below a threshold, the second for those
                     above.
        threshold  : Value in data units according to which the colors from
                     textcolors are applied. If None (the default) uses the
                     middle of the colormap as separation.

    Further arguments are passed on to the created text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max()) / 2.0

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center", verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = mpl.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[(im.norm(data[i, j]) > threshold).astype(int)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts
