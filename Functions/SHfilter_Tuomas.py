import numpy as np
import pyshtools as sh
from pyshtools.expand import SHExpandDH


def SHfilter(dat,lmax = 7):
    """Smoothing filter for data on an NxN or 2NxN grid. Truncated at
    spherical harmonic order lmax (default lmax = 5). """
    shapes = dat.shape
    if shapes[0] == shapes[1]:
        s = 1
    elif 2*shapes[0] == shapes[1]:
        s = 2
    coeffs = SHExpandDH(dat,sampling = s)
    clm = sh.SHCoeffs.from_array(coeffs)
    clm.coeffs[:, lmax:, :] = 0
    grid = clm.expand(grid = 'DH2')
    d = grid.data
    return d

#
# def SHfilter(dat,lmax = 5):
#     """Smoothing filter for data on an NxN or 2NxN grid. Truncated at
#     spherical harmonic order lmax (default lmax = 5). """
#     shapes = dat.shape
#     if shapes[0] == shapes[1]:
#         s = 1
#     elif 2*shapes[0] == shapes[1]:
#         s = 2
#     coeffs = SHExpandDH(dat,sampling = s)
#     clm = sh.SHCoeffs.from_array(coeffs)
#     clm.coeffs[:, lmax, :] = 0
#     grid = clm.expand(grid = 'DH2')
#     d = grid.data
#     return d
