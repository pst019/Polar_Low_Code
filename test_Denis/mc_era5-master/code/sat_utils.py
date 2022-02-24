# -*- coding: utf-8 -*-
"""Functions to load satellite data."""
import cartopy.crs as ccrs
import rasterio

import mypaths


def read_raster_stereo(filename=mypaths.avhrr_file):
    """
    Read the image and essential metadata from a GeoTIFF file.

    Used for reading satellite imagery in polar stereographic projection
    from NERC Satellite Receiving Station, Dundee University, Scotland
    (http://www.sat.dundee.ac.uk/).

    Parameters
    ----------
    filename: str or path-like
        path to the GeoTIFF file

    Returns
    -------
    im: numpy array
        image data
    extent: list
        image extent in units of the original projection
    crs: cartopy.crs.Projection
        Stereographic projection of the image
    """
    with rasterio.open(filename, "r") as src:
        proj = src.crs.to_dict()
        crs = ccrs.Stereographic(
            central_latitude=proj["lat_0"],
            central_longitude=proj["lon_0"],
            false_easting=proj["x_0"],
            false_northing=proj["y_0"],
            true_scale_latitude=proj.get("lat_ts"),
            globe=ccrs.Globe(datum=proj.get("datum")),
        )

        # read image into ndarray
        im = src.read().squeeze()

        # calculate extent of raster
        # --------------------------
        # note that the order of transform parameters changed since rasterio v1.0
        # See https://rasterio.readthedocs.io/en/stable/topics/migrating-to-v1.html for details
        #
        # One of the biggest API changes on the road to Rasterio 1.0 is the full deprecation of
        # GDAL-style geotransforms in favor of the affine library.
        # For reference, an affine.Affine() looks like:
        # affine.Affine(a, b, c,
        #               d, e, f)
        # and a GDAL geotransform looks like:
        # (c, a, b, f, d, e)
        xmin = src.transform[2]
        xmax = src.transform[2] + src.transform[0] * src.width
        ymin = src.transform[5] + src.transform[4] * src.height
        ymax = src.transform[5]
        extent = [xmin, xmax, ymin, ymax]

    return im, extent, crs
