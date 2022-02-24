#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  5 14:52:50 2018

@author: pst019

conda install -c conda-forge metpy

https://github.com/Unidata/MetPy/tree/master/staticdata
"""

import os
user = os.getcwd().split('/')[2]

import sys
sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_mapplot import * #make basemap plots
from f_imp_thorpex import * #import the thorpex data

Mediadir= '/media/'+user+'/1692A00D929FEF8B/'
#homedir= "/home/"+user+"/home/"

year, month = 2008, 3
day, hour= 3, 12  
thor= data(year, month, day, hour, plevels= np.arange(950, 450, -10))


from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from metpy.cbook import get_test_data
import metpy.calc as mpcalc
from metpy.plots import SkewT, Hodograph
from metpy.units import units, concatenate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# Create a new figure. The dimensions here give a good aspect ratio
fig = plt.figure(1, figsize=(8, 8))
plt.clf()
skew = SkewT(fig, rotation=45)
# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot

dropid=6

p= thor.pres[dropid]* units.hPa
T= thor.T[dropid]* units.degC
RH= thor.RH[dropid]/100

from metpy.calc.thermo import *
Td= dewpoint_rh(T, RH)

skew.plot(p, T, 'r')
skew.plot(p, Td, 'g')


thor2= data(year, month, day, hour, level='one_profile', dropid= dropid+2)
p= thor2.pres* units.hPa
T= thor2.T* units.degC

Td= thor2.Td * units.degC

#from metpy.io import get_upper_air_data
#from metpy.io.upperair import UseSampleData
#with UseSampleData(): # Only needed to use our local sample data
## Download and parse the data
#    dataset = get_upper_air_data(datetime(1999, 5, 4, 0), 'OUN')
#    
#p = dataset.variables['pressure'][:]
#T = dataset.variables['temperature'][:]
#Td = dataset.variables['dewpoint'][:]
#u = dataset.variables['u_wind'][:]
#v = dataset.variables['v_wind'][:]



skew.plot(p, T, 'y')
skew.plot(p, Td, 'b')
#skew.plot_barbs(p, u, v)
skew.ax.set_ylim(1000, 400)
skew.ax.set_xlim(-30, 10)

# Calculate LCL height and plot as black dot
#l = mpcalc.lcl(p[0], T[0], Td[0])
#lcl_temp = mpcalc.dry_lapse(concatenate((p[0], l[0])), T[0])[-1].to('degC')
#skew.plot(l, lcl_temp, 'ko', markerfacecolor='black')
# Calculate full parcel profile and add to plot as black line
#prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')
#skew.plot(p, prof, 'k', linewidth=2)
# Example of coloring area between profiles
#skew.ax.fill_betweenx(p, T, prof, where=T>=prof, facecolor='blue', alpha=0.4)
#skew.ax.fill_betweenx(p, T, prof, where=T<prof, facecolor='red', alpha=0.4)
# An example of a slanted line at constant T -- in this case the 0
# isotherm
#l = skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)
# Add the relevant special lines
skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
skew.plot_mixing_lines()
# Show the plot


"""make the hodograph"""
#ax_hod = inset_axes(skew.ax, '40%', '40%', loc=1)
#h = Hodograph(ax_hod, component_range=80.)
#h.add_grid(increment=20)
#h.plot_colormapped(u, v, np.hypot(u, v))

#plt.show()

