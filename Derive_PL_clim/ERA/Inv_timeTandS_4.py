#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 13:25:51 2017

@author: pst019
"""

import os
user = os.getcwd().split('/')[2]

import pickle
import sys
sys.path.insert(0, '/home/'+user+'/polar_low_code/Functions/')
from f_imp_ERA2 import *
from f_impLists import *
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)

Mediadir= '/media/'+user+'/1692A00D929FEF8B/'

fignr= 1

year, month= 2009, 4

#day, hour=  22, 18
#t= ((day-1)*24+hour)//6 #hours since beginning of month
t= 18
tcomb= t*2
#tcomb = 134
#t= tcomb//2
print(year, month, (t//4)+1,  t%4*6)
  
minlon= 0 #-45
maxlon= 100 #65
minlat= 30 #54
maxlat= 70 #81

Sdurationlim= 2 #6 #in hours
   
"""import data"""
T=TRACK_list_month2(year, month, duration_limit= 1, timestep='threehours')
#T.local(minlon, maxlon, minlat, maxlat) #for combined
ASRT=TRACK_list_month2(year, month, duration_limit= 1, timestep='threehours', model='ASR')

        
TPL=TRACK_PLlist_2(year, month, name='3')
ASRTPL=TRACK_PLlist_2(year, month, model='ASR', name='4')


d= data('SST', year, month, tstart= t, tend= t+1)          

S= STARSList(year, smonth= month, filename='STARS_TRACKS.csv', durationlim= Sdurationlim, timestep= 'threehours')
S2= STARSList(syear= year, smonth= month, filename='STARS_TRACKS_south.csv', durationlim= Sdurationlim, timestep= 'threehours')
S.PLnumber= np.concatenate((S.PLnumber, S2.PLnumber))
S.lat= np.concatenate((S.lat, S2.lat))
S.lon= np.concatenate((S.lon, S2.lon))
S.tPL= np.concatenate((S.tPL, S2.tPL))


#map= Lat_Lon_map(lllon= minlon, lllat= minlat, urlon= maxlon, urlat= maxlat)
map= AA_map()
#map= Polar_map(latsouth= d.lat[-1])

grid= np.meshgrid(d.lon, d.lat) #pcolormesh needs corner points for projection
Lon, Lat= map(grid[0], grid[1])
#

"""vorticity truncated data"""
plt.figure(fignr)
fignr += 1
plt.clf()

map= AA_map()


ntrunc, ntrunc2 = 100, 40
vorttrunc_dir=Mediadir+'ERA/pickle/Vort_trunc/Vort_trunc'
vortfilter_comb= pickle.load(open(vorttrunc_dir+str(year)+'_'+str(month).zfill(2)+'_T'+str(ntrunc)+'-'+str(ntrunc2), 'rb'), encoding='latin1').astype('float')


PlotVort(Lon, Lat, vortfilter_comb[tcomb], map, maxlevel= 12)
plt.title('truncated vorticity')
#
#PlotTRACK_PL(TPL, t, map)
PlotSTARSandMatchingTRACK(S, T, TPL, t, tcomb, map)


"""wind and surface pressure"""
plt.figure(fignr)
fignr +=1
plt.clf()
#plt.subplot(2, 3, 2)

d= data(['MSLP'], year, month, tstart= t, tend= t+1)
d.imp_u10()
U10= np.sqrt(d.u10[0]**2+d.v10[0]**2)

#map= Polar_map(latsouth= d.lat[-1])
map= AA_map()

PlotWindVelo(Lon, Lat, U10, map, Umax= 25)
PlotWind(d.lon, d.lat, d.u10[0], d.v10[0], map)
PlotContours(Lon, Lat, d.MSLP[0], map)

plt.title('surface wind and pressure')
#PlotTRACK_PL(TPL, t, map)
#PlotTRACK_PL(ASRTPL, tcomb, map)
#
#PlotSTARS_PL(S, tcomb, map)




"""static stability (SST-T500)"""    
#plt.figure(fignr)
#fignr +=1
#plt.clf()
##plt.subplot(2, 3, 5)
#
#map= Polar_map(latsouth= d.lat[-1])
#d= data(['SST', 'T'], year, month, tstart= t, tend= t+1)
#PlotStaticStability(Lon, Lat, (d.SST[0]- d.T[0,2]) , map)
#
#
#d.impvar('uPL')
#d.impvar('vPL')
#
#uPLavg= np.mean(d.uPL[0], axis=0)
#vPLavg= np.mean(d.vPL[0], axis=0)
#
#PlotWindVelo(Lon, Lat, np.sqrt(uPLavg**2+ vPLavg**2), map, Umax= 40)
#PlotWind(d.lon, d.lat, uPLavg, vPLavg, map)
#
#
#plt.title('static stability and mean troposphere wind ')
#PlotTRACK_PL(TPL, t, map)
#



"""tropopause wind"""
#plt.figure(fignr)
#fignr +=1
#plt.clf()
##plt.subplot(2, 3, 3)
#
#d= data(['uPV', 'vPV'], year, month, tstart= t, tend= t+1)
#
#PlotWindVelo(Lon, Lat, np.sqrt(d.uPV[0]**2+ d.vPV[0]**2), map, Umax= 80)
#PlotWind(d.lon, d.lat, d.uPV[0], d.vPV[0], map)
#
##map= Polar_map(latsouth= d.lat[-1])
#map= AA_map()
#
#
#plt.title('tropopause wind ')
#PlotTRACK_PL(TPL, t, map)

"""Plot all STARS polar lows happening now with their tracks """
#Sxpt, Sypt= map(S.lon[S.tPL== tcomb], S.lat[S.tPL== tcomb])
#map.plot(Sxpt,Sypt,'ro')
#
#for n in range(len(Sxpt)): #write PLnumber
#    plt.text(Sxpt[n]*1.01,Sypt[n]*1.01, str(int(S.PLnumber[S.tPL== tcomb][n])), fontsize=10, ha='center',va='top',color='red')
#
#
#SPLnrs= S.PLnumber[S.tPL== tcomb]
#for n in SPLnrs:
#    lonp= S.lon[S.PLnumber==n]
#    latp= S.lat[S.PLnumber==n]
#    xpt, ypt= map(lonp, latp)
#    map.plot(xpt,ypt,'rx') 
#    map.plot(xpt,ypt,'r-') 
#    map.plot(xpt[0],ypt[0],'r^', label= str(int(n)))


"""plot some TRACK properties"""
PlotSTARSandMatchingTRACK(S, T, TPL, t, tcomb, map)
"""Plot the STARS PL, the corresponding TRACK cyclone, 
#all TRACK and TRACKPL points now and PLpoints during the corresponding TRACK cyclone"""
#SPLnrs= PlotSTARS_PL(S, tcomb, map)
#PlotallTRACKpoints(ASRT, tcomb, map, color='g')
#
#PlotallTRACK_PLpoints(ASRTPL, tcomb, map, color='cyan')
#TPLnrs= MatchingTRACKtoSTARS(SPLnrs, model='ASR')
#PlotGivenTRACKcyclones(ASRT, TPLnrs, map, color= 'g')
#PlotPLpointsinTRACKcyclone(ASRTPL, TPLnrs, map, color='cyan')
#    
    
"""calculate the thermal wind"""
#d.impvar('Geop')
#thick=d.Geop[0,1]- d.Geop[0,0]
#u_T, v_T= CalculateThermalWind(thick, d.lat, d.lon)
#
#        
#plt.figure(fignr)
#fignr +=1
#plt.clf()
##plt.subplot(2, 3, 6)
#
#map= Polar_map(latsouth= d.lat[-1])
#PlotContours(Lon, Lat, thick, map)
#PlotColorMap(Lon, Lat, thick, map, variable='geopotential thickness')
#
#PlotWind(d.lon, d.lat, u_T, v_T, map, alen= 30)
#PlotTRACK_PL(TPL, t, map)
#
#plt.title('thermal wind')
#
#
#"""calculate the angle between thermal wind and troposphere wind"""
#alpha= np.rad2deg(np.arccos((u_T*uPLavg + v_T*vPLavg)/(np.sqrt(u_T**2+v_T**2)* np.sqrt(uPLavg**2+vPLavg**2))))
#
#latdist= 4
#
#
##PLt= TPL.PLlist[1]== t
##PL1= TPL.PLlist[-1]==1
##PLnow= np.array([a and b for a, b in zip(PLt, PL1)])
##PLlistnow= TPL.PLlist[:, PLnow]
##
##alpham= []
##for i, t in enumerate(PLlistnow[1]):  #loop through every PL point            
##    londist= int(latdist//np.cos(np.deg2rad(PLlistnow[3][i])))  
##    Tlatraw= np.array(np.round(-(PLlistnow[3][i]-d.lat[0])*2), dtype= int)
##    Tlonraw= np.array(np.round((PLlistnow[2][i] - d.lon[0])*2), dtype= int)%720         
##    mask= radMask2D((len(d.lat), len(d.lon)), (Tlatraw, Tlonraw), radiusx= latdist, radiusy= londist) #lat lon mask
##
##    alpham += [np.mean(alpha[mask])]
##
##print('PLnr',PLlistnow[0], 'wind angle', alpham)
#
#fignr= 6
#"""plot the angle of the thermal wind"""
#plt.figure(fignr)
#fignr +=1
#plt.clf()
#map= Polar_map(latsouth= d.lat[-1])
#
#PlotColor_from0(Lon, Lat, alpha, map, boxnr= 9, maxlevel= 180, col='BlueRed', label='Angle of thermal wind vs avg wind')
#
#
#
#"""static stability (SST-T500)"""    
##plt.figure(fignr)
##fignr +=1
##plt.clf()
###plt.subplot(2, 3, 5)
##
##map= Polar_map(latsouth= d.lat[-1])
##d= data(['SST', 'T'], year, month, tstart= t, tend= t+1)
##PlotStaticStability(Lon, Lat, (d.SST[0]- d.T[0,2]) , map)
##plt.title('static stability (SST-T500)')
##
##PlotTRACK_PL(TPL, t, map)
#
#"""Stability (SST-theta_e500)"""
##plt.figure(fignr)
##fignr+=1
##plt.clf()
####plt.subplot(2, 3, 1)
##
##d= data(['SST', 'T', 'sHum', 'MSLP'], year, month, tstart= t, tend= t+1)
##RCp, Lc, Cp= 2/7, 2501, 1#R/Cp
##
##theta500= d.T[0,2]*(d.MSLP[0]/500)**RCp
##theta_e500= theta500* np.exp(d.sHum[0,2]* Lc/(d.T[0,2]* Cp))
##stab= d.SST[0] - theta_e500
##
##map= Polar_map(latsouth= d.lat[-1])
##
###PlotStatStabTheta_e(Lon, Lat, stab, map) #does not exist anymore
##PlotContourStatStabTheta_e(Lon, Lat, stab, map)
##plt.title('static stability (SST-theta_e500)')
#
#
#"""mean wind of pressure level 925 and 700hPa"""
##fignr= 5
##
##plt.figure(fignr)
##fignr +=1
##plt.clf()
#
##d.impvar('uPL')
##d.impvar('vPL')
##
##uPLavg= np.mean(d.uPL[0], axis=0)
##vPLavg= np.mean(d.vPL[0], axis=0)
##
##PlotWindVelo(Lon, Lat, np.sqrt(uPLavg**2+ vPLavg**2), map, Umax= np.max(np.abs(uPLavg)))
##PlotWind(d.lon, d.lat, uPLavg, vPLavg, map)
##
##map= Polar_map(latsouth= d.lat[-1])
##
##plt.title('mean troposphere wind ')
##PlotTRACK_PL(TPL, t, map)