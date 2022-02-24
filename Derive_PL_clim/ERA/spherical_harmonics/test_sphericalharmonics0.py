# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy            # numpy for arrays
import shtns            # shtns module
                                         #   compiled and installe using ./configure --enable-python && make && make install
 
lmax = 7                        # maximum degree of spherical harmonic representation.
mmax = 3                        # maximum order of spherical harmonic representation.

sh = shtns.sht(lmax, mmax)              # create sht object with given lmax and mmax (orthonormalized)
# sh = shtns.sht(lmax, mmax, mres=2, norm=shtns.sht_schmidt | shtns.SHT_NO_CS_PHASE)    # use schmidt semi-normalized harmonics
 
nlat, nphi = sh.set_grid()              # build default grid (gauss grid, phi-contiguous)
print(sh.nlat, sh.nphi)                 # displays the latitudinal and longitudinal grid sizes.
 
cost = sh.cos_theta                     # latitudinal coordinates of the grid as cos(theta)
el = sh.l                                       # array of size sh.nlm giving the spherical harmonic degree l for any sh coefficient
l2 = el*(el+1)                          # array l(l+1) that is useful for computing laplacian

### use advanced options to create a regular grid, theta-contiguous, and with south-pole comming first.
# nlat = lmax*2
# nphi = mmax*3
# grid_typ = shtns.sht_reg_fast | shtns.SHT_THETA_CONTIGUOUS | shtns.SHT_SOUTH_POLE_FIRST
# polar_opt_threshold = 1.0e-10

ylm = sh.spec_array()           # a spherical harmonic spectral array, same as numpy.zeros(sh.nlm, dtype=complex)
vr = sh.spat_array()            # a spatial array, same as numpy.zeros(sh.spat_shape)
 
ylm[sh.idx(1,0)] = 1.0          # set sh coefficient l=1, m=0 to value 1

print(ylm[sh.l == 1])                   # print all l=1 coefficients
print(ylm[sh.m == 1])                   # print all m=1 coefficients

ylm = ylm * l2                          # multiply by l(l+1)

y = sh.synth(ylm)               # transform sh description ylm into spatial representation y (scalar transform)
print(y)                                # display spatial field
i_theta = sh.nlat/2
i_phi = 1
print(y[i_theta, i_phi])        # spatial element of coordinate i_theta, i_phi

zlm = sh.analys(y)              # transform the spatial field back to spectral
print(zlm)

### compute gradients :
dy_dt, dy_dp = sh.synth_grad(ylm)       # compute gradients of ylm on a spatial grid.
 
#   71 ### vector transforms :
#   72 #slm, tlm = sh.analys(v_theta,v_phi)            # with two arguments, sh.analys() performs a 2D vector transform (spherical coordinates without radial)
#   73 #qlm, slm, tlm = sh.analys(v_r,v_theta,v_phi)   # with three arguments, sh.alanys() performs a 3D vector tansform (spherical coordinates)
#   74 #
#   75 #v_theta, v_phi = sh.synth(slm, tlm)    # same hold for synthesis, using spheroidal/toroidal vector spherical harmonics.
#   76 #v_r = sh.synth(qlm)
#   77 
