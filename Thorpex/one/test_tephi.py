#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  5 14:52:50 2018

@author: pst019

conda install -c conda-forge metpy

https://github.com/Unidata/MetPy/tree/master/staticdata
"""

from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from metpy.cbook import get_test_data
import metpy.calc as mpcalc
from metpy.io import get_upper_air_data
from metpy.plots import SkewT, Hodograph
from metpy.units import units, concatenate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from metpy.io.upperair import UseSampleData
with UseSampleData(): # Only needed to use our local sample data
# Download and parse the data
    dataset = get_upper_air_data(datetime(1999, 5, 4, 0), 'OUN')
p = dataset.variables['pressure'][:]
T = dataset.variables['temperature'][:]
Td = dataset.variables['dewpoint'][:]
u = dataset.variables['u_wind'][:]
v = dataset.variables['v_wind'][:]

# Create a new figure. The dimensions here give a good aspect ratio
fig = plt.figure(1, figsize=(10, 8))
plt.clf()
skew = SkewT(fig, rotation=45, subplot= (1, 2, 1))
# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot
skew.plot(p, T, 'r')
skew.plot(p, Td, 'g')
skew.plot_barbs(p, u, v, xloc= .2)
skew.ax.set_ylim(1000, 100)
skew.ax.set_xlim(-40, 60)

# Calculate LCL height and plot as black dot
l = mpcalc.lcl(p[0], T[0], Td[0])
lcl_temp = mpcalc.dry_lapse(concatenate((p[0], l[0])), T[0])[-1].to('degC')
#skew.plot(l, lcl_temp, 'ko', markerfacecolor='black')
# Calculate full parcel profile and add to plot as black line
prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')
skew.plot(p, prof, 'k', linewidth=2)
# Example of coloring area between profiles
skew.ax.fill_betweenx(p, T, prof, where=T>=prof, facecolor='blue', alpha=0.4)
skew.ax.fill_betweenx(p, T, prof, where=T<prof, facecolor='red', alpha=0.4)
# An example of a slanted line at constant T -- in this case the 0
# isotherm
l = skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)
# Add the relevant special lines
skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
skew.plot_mixing_lines()
# Show the plot




#ax_hod = inset_axes(skew.ax, '40%', '40%', loc=1)
#h = Hodograph(ax_hod, component_range=80.)
#h.add_grid(increment=20)
#h.plot_colormapped(u, v, np.hypot(u, v))

#plt.show()

