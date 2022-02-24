#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 13:25:51 2017

@author: pst019
from ASR-ERA_Inv_random_3 of folder ASR
"""


import pickle
import sys
import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'


import matplotlib.pyplot as plt

sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')
from f_imp_ERA2 import *
from f_imp_ASR import *
from f_impLists import *
from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_meteo import *

from random import randint


"""global variables"""

fignr= 0

#maptype='AA'
maptype='AA_half'
#maptype='Lambert'
#maptype='polar'

if maptype=='AA': plt.figure(fignr, figsize= (6, 6))
elif maptype=='AA_half': plt.figure(fignr, figsize= (6, 4.7))
else: plt.figure(fignr, figsize= (6, 4.5))
fignr += 1
plt.clf()

"""Lambert coordinates"""
lllon, lllat, urlon, urlat= -15, 63, 60, 75
lat0, lon0= 75, 0 #(lat0, lon0) = center point

                   
if maptype== 'AA': map = AA_map()
elif maptype== 'AA_half': map = AA_map_half()
elif maptype=='polar': map= Polar_map(latsouth= 50) #d.lat[-1])
else: map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)





                    
"""time"""
year=2008
month=3
day= 3
hour= 18

t= (day-1)*4 + hour//6


#var= 'PressWind'
var= 'PressWind_advanced'
#var= 'PressWind_Geop'
#var= 'Theta700'
#var= 'Vort'
#var= 'Vort_trunc'
#var= 'Theta_e850'
#var = 'U_trop'
#var= 'front'
#var = 'Cloud_med'
#var = 'OW'
#var= 'RH850'
var='SH850'


"""c) plot characteristics"""
pleveldist= 4

"""end global variables """

title_extra='jet' #''
"""module starts"""
   

save= True
savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/Fields/'
 
   
    
print(year, month, day, hour, 't: ', t)



d= data(['MSLP'], year, month, tstart= t, tend= t+1) #just to have d.lat
grid= np.meshgrid(d.lon, d.lat) #pcolormesh needs corner points for projection
Lon, Lat= map(grid[0], grid[1])


"""Plot location of TRACK and STARS PLs"""
Track = False #plot the TRACK PL
STARS = False # plot the STARS PL


"""ERA-vorticity truncated data"""
if var == 'Vort_trunc':
    ntrunc, ntrunc2 = 100, 40
    vorttrunc_dir=Mediadir+'ERA/pickle/Vort_trunc/Vort_trunc'
    vortfilter_comb= pickle.load(open(vorttrunc_dir+str(year)+'_'+str(month).zfill(2)+'_T'+str(ntrunc)+'-'+str(ntrunc2), 'rb'), encoding='latin1').astype('float')

    PlotVort(Lon, Lat, vortfilter_comb[t], map, maxlevel= 12)
    plt.title('ERA truncated vorticity')
#    PlotTRACK_PL(TPLera, t, map, track=False)


if var == 'Vort':  
    d= data(['Vort'], year, month, tstart= t, tend= t+1) #just to have d.lat

    PlotVort(Lon, Lat, d.vort[0], map)#, maxlevel= 12)
    plt.title('ERA vorticity')

fignr= 16

if var == 'OW':
    d.impvar('uPL') #wind at 925 and 700hPa
    d.impvar('vPL')
    
    PlotContours(Lon, Lat, d.MSLP[0], map, leveldist= pleveldist, numbers=False )
#    PlotWindVelo(Lon, Lat, np.sqrt(d.uPL[0,0]**2 + d.vPL[0,0]**2), map, Umax= 25)
    
    dy= (d.lat[0]- d.lat[1])*110E3 #distance along axis 1
    dx= dy* np.tile(np.cos(np.deg2rad(d.lat[::-1])), (len(d.lon),1)) #distance along axis 0
    v= d.vPL[0,0, ::-1].T # transformation to (east- west, south - north) - grid. since the latitudes go from north to south
    u= d.uPL[0,0, ::-1] .T       

#    vort= (np.gradient(v, dx , axis= 0)- np.gradient(u, dy, axis= 1))*1E5   
#    PlotVort(Lon, Lat, vort[:,::-1].T, map)

#    shear= (np.gradient(v, dx , axis= 0)+ np.gradient(u, dy, axis= 1))*1E5   
#    PlotVort(Lon, Lat, shear[:,::-1].T, map, label='Shear [10$^{-5}$ 1/s]')

    stretching= (np.gradient(u, dx , axis= 0)- np.gradient(v, dy, axis= 1))*1E5                                              
    PlotVort(Lon, Lat, stretching[:,::-1].T, map, label='Stretching [10$^{-5}$ 1/s]')

#    OW= np.sqrt(vort**2 - shear**2 - stretching**2)
#    PlotVort(Lon, Lat, OW[:,::-1].T, map, label='Sqrt OW [10$^{-5}$ 1/s]')   

#    PlotWind(d.lon, d.lat, d.uPL[0, 0], d.vPL[0, 0], map, nx= 60, ny= 100)


"""ERA wind and surface pressure"""
if var == 'PressWind':   
    d.imp_u10()
    U10= np.sqrt(d.u10[0]**2+d.v10[0]**2)

    PlotWindVelo(Lon, Lat, U10, map, Umax= 25)
   
    PlotWind(d.lon, d.lat, d.u10[0], d.v10[0], map)
    PlotContours(Lon, Lat, d.MSLP[0], map, leveldist= pleveldist)


if var == 'PressWind_advanced':
    PlotContours(Lon, Lat, d.MSLP[0], map, leveldist= 1, numbers=False)

    d.imp_u10()
    U10= np.sqrt(d.u10[0]**2+d.v10[0]**2)
    PlotWindVelo(Lon, Lat, U10, map, Umax= 25, color='YlBu')

    PlotLocalMax(d.MSLP[0], threshold=1010, distance=15, map= map, lon=d.lon, lat=d.lat,
                     data2=U10, threshold2=10, distance2= 7) 

    PlotLocalMax(U10, threshold=15, distance=10, map= map, lon=d.lon, lat=d.lat, typ='max',
                 color='orange', dot=False, roundorder= 0)        



if var == 'PressWind_Geop':   

    d.SST= np.ma.array(d.SST, mask= np.isnan(d.SST))
    
    d.imp_u10()
    U10= np.sqrt(d.u10[0]**2+d.v10[0]**2)

    PlotWindVelo(Lon, Lat, U10, map, Umax= 25)
    
    PlotWind(d.lon, d.lat, d.u10[0], d.v10[0], map)
    PlotContours(Lon, Lat, d.MSLP[0], map, leveldist= pleveldist)
    
    d.impvar('Geop500')
    PlotContours(Lon, Lat, d.Geop500[0]/9.81 , map, nrlevels= 10, color='r')


"""potential temperature at 700 hPa  + static stability"""
if var == 'Theta700':  
    PlotContours(Lon, Lat, d.MSLP[0], map, leveldist= pleveldist)
    
    Theta= PotTemp(d.T[0, 1], plev= 700)
    bounds= np.arange(264, 290, 1)
    PlotColorMap3(Lon, Lat, Theta, map, symetric=False, color='RdBu', bounds= bounds, label= r"Potential temperature at 700hPa [K]")
    plt.title(r"ERA $\theta_{700}$")
    
    theta500= d.T[0,2]*(1000/500)**RCp
    thetaSST= d.SST[0]*(1000/d.MSLP[0])**RCp
    PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')
    

"""equivalent potential temperature at 850 hPa"""
if var == 'Theta_e850':
    d.impvar('T850')
    d.impvar('sHum850')
    Theta_e= EquiPotTemp(d.T850[0], d.sHum850[0], plev= 850)
    
    PlotContours(Lon, Lat, d.MSLP[0], map, leveldist= pleveldist)   
    PlotColorMap3(Lon, Lat, Theta_e, map, symetric=False, color='RdBu', bounds= np.arange(264, 292, 1), label= r"$\theta_{e, 850}$ [K]")
    

    
"""ERA - tropopause wind"""
if var == 'U_trop':
    d.impvar('uPV')
    d.impvar('vPV')
    d.impvar('Geop500')
    
    PlotWindVelo(Lon, Lat, np.sqrt(d.uPV[0]**2+ d.vPV[0]**2), map, Umax= 80)
    #PlotColor_setbounds(Lon, Lat, np.sqrt(d.uPV[0]**2+ d.vPV[0]**2), map, bounds= [0, 30.5, 50, 80], col='Blues')
    
    PlotWind(d.lon, d.lat, d.uPV[0], d.vPV[0], map)
    PlotContours(Lon, Lat, d.Geop500[0]/9.81 , map, nrlevels= 10, color='r')
    

"""frontal zones"""
if var == 'front':
    d.impvar('T850')
    d.impvar('sHum850')
    Theta_e= EquiPotTemp(d.T850[0], d.sHum850[0], plev= 850)
        
    #PlotColorMap3(Lon, Lat, Theta_e, map, symetric=False) #, color='blue') #, boxnr=15)#, maxlevel= 12)
    PlotContours(Lon, Lat, Theta_e, map, nrlevels= 20, color='r')
    PlotContours(Lon, Lat, d.MSLP[0], map, nrlevels= 20)
    #    PlotContours(Lon, Lat, d.Geop500[0]/9.81 , map, nrlevels= 10, color='r')


    latdist= 55.
    londist= latdist*np.cos(np.deg2rad(d.lat))
    londist= np.tile(londist, (720,1)).T
    g= np.gradient(Theta_e, latdist, londist) 
    absg= np.sqrt(g[0]**2+ g[1]**2) *100 #change in Theta_e /100km
    #    
    PlotColorMap3(Lon, Lat, absg, map, symetric=False, color='blue', bounds= np.linspace(0, 15, 6), label= r"$|\nabla T|$ [K/100km]")
    #    PlotColorMap3(Lon, Lat, Theta_e, map, symetric=False, color='RdBu', bounds= np.arange(264, 292, 1), label= r"Equivalent potential temperature [K]")
    

if var == 'RH850':
    PlotContours(Lon, Lat, d.MSLP[0], map, leveldist= 1, numbers=False)
    d.impvar('sHum850')
    d.impvar('T850')
    import metpy.calc as mpcalc
    from metpy.units import units    
    dRH= mpcalc.relative_humidity_from_specific_humidity(d.sHum850[0], d.T850[0]* units.K, 850* units.hPa)
    
    PlotColorMap4(Lon, Lat, dRH , map, label='Relative humidity at 850 hPa', color='blue', bounds= np.array([0, 0.4, 0.7, 0.8, 0.9, 0.95, 1]) )


if var == 'SH850':
    PlotContours(Lon, Lat, d.MSLP[0], map, leveldist= 1, numbers=False)
    d.impvar('sHum850')
    d.impvar('T850')
    import metpy.calc as mpcalc
    from metpy.units import units    
    dRH= mpcalc.relative_humidity_from_specific_humidity(d.sHum850[0], d.T850[0]* units.K, 850* units.hPa)

    PlotColorMap4(Lon, Lat, d.sHum850[0]*1E3 , map, label= 'Specific humidity [g/kg] at 850 hPa',
                  color='jet_r', bounds= np.arange(0, 2.3, .2))    
    PlotContours(Lon, Lat, dRH, map, levels= [0.9], numbers= False, color='white')


"""ERA -medium cloud"""
if var == 'Cloud_med':
    if year == 2000:
        d.impvar('MedCloud')
        
        PlotColorMap3(Lon, Lat, d.MedCloud[0], map, symetric=False, color='blue') #, boxnr=15)#, maxlevel= 12)
        plt.title('ERA Medium Cloud')
        
    else: print('data not downloaded for this year')



if Track== True:
    TPLera=TRACK_PLlist_2(year, month, name='7')
    PlotTRACK_PL(TPLera, t, map, track=False)
#    PlotTRACK_PL(TPLera, t, map, nowPL=False, track=True)

if STARS == True:
    S= STARSList(year, smonth= month, filename='STARS_TRACKS.csv', durationlim= 3, timestep= 'threehours')
    PlotSTARS_PL(S, t, map)#, track=False, color= 'b')



if save== False:
    plt.title('ERA-I '+str(d.datetime[0])[:-3]+ title_extra)

plt.tight_layout()

if save== True:
    print('save')
    savename= savedir + 'ERA_I_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+str(hour).zfill(2)+'_'+var+title_extra
    plt.savefig(savename, pad_inches=0)
    print(savename)