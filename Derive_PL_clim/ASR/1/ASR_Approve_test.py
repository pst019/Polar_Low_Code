#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 13:25:51 2017

@author: pst019
"""

import pickle
import sys
sys.path.insert(0, '/home/'+user+'/codeDerive_PL_clim/ERA/')
from f_imp_ERA2 import *
from f_imp_ASR import *
from f_impLists import *
sys.path.insert(0, '/home/'+user+'/codeAROME/')
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)


print('make filedir = "/media/'+user+'/1692A00D929FEF8B/PL/TRACKoutput/ERA_VOR850_" in line 426 in f_impLists')

fignr= 4


year, month= 2002, 1
day, h3=  12, 2 #hour3 is hour*3
t = (day-1)*8 + h3 #hours since beginning of month
#t= 3 #52

print(year, month, day,  h3*3)
  
Sdurationlim= 2 #6d #in hours
   
"""import data"""
T=TRACK_list_month2(year, month, duration_limit= 2, timestep='threehours', model='ASR', approved= False)
 
"""SLP and wind speed"""
#plt.figure(fignr)
#fignr+=1
#plt.clf()
###plt.subplot(2, 3, 1)
#
#d= dataASR(year, month, day, sh3= h3, eh3= h3+2, level='surf', reso= 1)
#d.impvar('SLP', level='surf')
#d.impvar('U', level='surf')
#d.impvar('V', level='surf')
#d.impvar('U', level=500)
#d.impvar('V', level=500)
#
##d= dataASR(year, month, day, eday= day+1, sh3= h3, eh3= h3-1, level='all')
#
#map= ASR_map(res='c')
#Lon, Lat= map(d.lon, d.lat)
#
#PlotContours(Lon, Lat, d.SLP[0], map)
#PlotWindVelo(Lon, Lat, np.sqrt(d.u10[0]**2+ d.v10[0]**2), map, Umax= 25)
#
#"""TRACK cyclones"""
#"""the position now"""
#Txpt, Typt= map(T.PLlist[2][T.PLlist[1]== t], T.PLlist[3][T.PLlist[1]== t])
#map.plot(Txpt,Typt,'bo')


"""SLP and wind speed"""
#plt.figure(fignr)
#fignr+=1
#plt.clf()
#map= ASR_map(res='c')
#
#PlotWindVelo(Lon, Lat, np.sqrt(d.u500[0]**2+ d.v500[0]**2), map, Umax= 70)
#
#
#"""TRACK cyclones"""
#"""the position now"""
#Txpt, Typt= map(T.PLlist[2][T.PLlist[1]== t], T.PLlist[3][T.PLlist[1]== t])
#map.plot(Txpt,Typt,'bo')
#for n in range(len(Txpt)):
#    plt.text(Txpt[n],Typt[n], str(int(T.PLlist[0][T.PLlist[1]== t][n])), fontsize=10, ha='center',va='top',color='black')
#
#
#"""make a "north of mask" """
#indexpt= np.where(np.logical_and(T.PLlist[0] == 2146, T.PLlist[1]==t))[0]
#
#lonpt= T.PLlist[2][np.logical_and(T.PLlist[0] == 2146, T.PLlist[1]==t)]
#latpt= T.PLlist[3][np.logical_and(T.PLlist[0] == 2146, T.PLlist[1]==t)]
#
#masknorth = np.logical_and.reduce((d.lat > latpt, d.lon > lonpt -1, d.lon < lonpt +1))
#
#U500north= np.max(np.sqrt(d.u500[0][masknorth]**2 + d.v500[0][masknorth]**2))
#
#plt.figure(fignr)
#fignr+=1
#plt.clf()
#map= ASR_map(res='c')
#PlotColorMap3(Lon, Lat, masknorth, map, maxlevel= None, bounds=None, color= 'RdBu', boxnr= 2)


""" SST data"""
plt.figure(fignr)
fignr += 1
plt.clf()

map= ASR_map(res='c')

d= dataASR(year, month, day, eday= day+1, sh3= h3, eh3= h3-1, level='surf')
d.impvar('SST', level='surf')
d.SST[d.SST < 271.5] = 0
PlotSurfTemp(Lon, Lat, d.SST[0]-273, map)

#"""truncated data"""
#
#map= Polar_map(latsouth= d.lat[-1])
#ntrunc, ntrunc2 = 100, 40
#vorttrunc_dir=Mediadir+'ERA/pickle/Vort_trunc/Vort_trunc'
#vortfilter_comb= pickle.load(open(vorttrunc_dir+str(year)+'_'+str(month).zfill(2)+'_T'+str(ntrunc)+'-'+str(ntrunc2), 'rb'), encoding='latin1').astype('float')
#
#PlotVort(Lon, Lat, vortfilter_comb[tcomb], map)
#plt.title('truncated vorticity '+str(year)+' '+str(month))
##
#
##
#"""TRACK cyclones"""
#"""the position now"""
#Txpt, Typt= map(T.lon[T.tPL== tcomb], T.lat[T.tPL== tcomb])
#map.plot(Txpt,Typt,'bo')
#for n in range(len(Txpt)):
#    plt.text(Txpt[n],Typt[n], str(int(T.PLnumber[T.tPL== tcomb][n])), fontsize=10, ha='center',va='top',color='black')
#
