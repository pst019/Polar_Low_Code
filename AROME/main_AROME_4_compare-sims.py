#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

This is from plot_comparison
"""

import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'



from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import numpy as np
import time
#from datetime import date, datetime, timedelta
from scipy import stats

# import own modules
import sys  #to import the functions from a different directory
sys.path.insert(0, '/home/'+user+'/home/Polar_Low/polar_low_code/Functions')

from f_mapplot import * #make basemap plots
from f_plot_fields import * #plot fields (onto basemap or cross section)
from f_imp_AROME import *  # Read in netcdf file
import f_imp_thorpex as thorpex
from f_meteo import *
from f_useful import *


"""global variables"""
fignr= 1

#maptype='AA'
maptype='AA_half'
#maptype='Lambert'

 
"""time specification: time of the plot"""
year, month = 2008, 3
day, hour= 4, 0
    
"""name of the experiment"""
#exp_name= '080303_warmctr'
#exp_name='080303_cold_ERA'
#exp_name= '080303_warmsens_noQH'
#exp_name='080303_warmsens_noFLX'
#exp_name='080303_warmsens_nocondens'
#exp_name='870226_cold_ERA_can'
#exp_name, fileday, filehour='080303_cold_pseudo2', 3, 0
#exp_name, fileday, filehour='080304_cold_pseudo', 4, 0
#exp_name='080303_coldsens_noFLX_AREA'
#exp_name, fileday, filehour='08030312_cycling', 3, 12

#exp_name, fileday, filehour='DA_080301_cycling', 4, 6
exp_name, fileday, filehour='DA_080303_CTR', 3, 0
#exp_name, fileday, filehour='DA_080303_noHARATU', 3, 0

#exp_name, fileday, filehour='DA_080303_noFLXarea', 3, 0
#exp_name, fileday, filehour='DA_08030212', 2, 12
#exp_name, fileday, filehour='DA_080303_p6SST', 3, 0
#exp_name, fileday, filehour='DA_080303_2FLX', 3, 0

#exp_name_2='DA_080303_noOCND2'
#exp_name_2='DA_080303_noHARATU'
exp_name_2='DA_080303_m6SST'


exp_list= [exp_name, exp_name_2]
nexp= len(exp_list)


"""time of the file"""
#fileday, filehour= 4, 0        
t= (day- fileday)*24 + (hour- filehour)  # -1 
#
#year, month = 1987, 2
#day, hour= 26, 12    
#
#fileday, filehour= 25, 12        
#t= (day- fileday)*24 + (hour- filehour)  # -1 

tlist= [t]
#tlist= [19] #np.arange(0, 3)#42, 3)                                            

"""Lambert coordinates"""
#lllon, lllat, urlon, urlat= -15, 65, 50, 75
#lat0, lon0= 75, 0 #(lat0, lon0) = center point

#lllon, lllat, urlon, urlat= 0, 67, 60, 75
#lat0, lon0= 70, 0 #(lat0, lon0) = center point

#to match with Thorpex                    
if year == 2008 and month== 3 and day== 3 and hour<15:
    lllon, lllat, urlon, urlat= -4, 70, 19, 75
    lat0, lon0= 75, 0 #(lat0, lon0) = center point

elif year == 2008 and month== 3 and day== 3 and hour>=15:
    lllon, lllat, urlon, urlat= -6, 69.5, 13, 74.5
    lat0, lon0= 75, 0 #(lat0, lon0) = center point

elif year == 2008 and month== 3 and day== 4:
    #for AA domain moved south
#    lllon, lllat, urlon, urlat= -3, 63.5, 15, 69.5
#    lat0, lon0= 70, 0 #(lat0, lon0) = center point    

    #for original AA domain
    lllon, lllat, urlon, urlat= -0.5, 67.5, 12.5, 68.8
    lat0, lon0= 68, -25 #(lat0, lon0) = center point    

  
"""c) plot characteristics"""
pleveldist=1

#var= 'PressWind'
#var= 'PressWindarrow'
#var= 'PressWind_advanced'
#var= 'PressWindVortmax'
#var= 'SLP_laplacian'
#var= 'SLP_filtered'
#var= 'MeanSLPminSLP'
#var= 'Theta_lev'

#var= 'Theta_e_lev'
#var= 'Front_lev'
#var= 'Grad_Theta_e_lev'
#bounds= np.arange(264, 286, 1) #for pottempe850
#var= 'Vort_lev'
#var='Shear_lev'
#var= 'Stretching_lev'
#var= 'OW_lev'
#var= 'Vort_lev_T'
#var= 'Thickness'
#var= 'test_rot'
#var= 'Geop_lev'
#var= 'Geop_wind_lev'

#var= 'Geop_laplacian_lev'
#var= 'Precip'
#var= 'Precip_acc'
#var= 'RH_lev'
#var= 'RH_lev_fromSH'
#var= 'SH_lev'

#var= 'TH'
#var= 'TH_thor' #for comparison with thorpex
#var= 'LH'
#var= 'LH_thor'
var= 'FLX'
#var= 'PBH'
#var= 'W_lev'
#var= 'PV_lev'
#var= 'Tsurf'
#var= 'T2m'
#var= 'Cloud_conv' # is somehow always 0.
#var= 'Cloud_high'
#var= 'Cloud_med'    
#var= 'Cloud_low'
#var= 'CAPE'
#var= 'CIN' #not very helpful
#var= 'LCL'  #not very helpful
#var= 'LFC' #not very helpful
#var= 'CTT' #best to compare against sat images
#var= 'WVTC'
#var= 'CWR'
#var='OLW' #top of atmosphere outgoing longwave radiation

""" presure level (if the variable ends with '_var' """
pn = 4
#1-950, 4- 850hPa, 6- 700hPa, 8 - 500hPa



Thorpex= False
Plotbox= False
boxlon, boxlat= [-4, 13, 18.5, -4, -4], [70, 69.3, 75, 76, 70] #to display the first lambert map
#boxlon, boxlat= [-5, 8, 10, -6.5, -5], [69.5, 69.5, 74.5, 74.5, 69.5]   # lllon, lllat, urlon, urlat= -5, 69.5, 10, 74.5 to display the second lamber map
#boxlon, boxlat= [-3, 12, 15, -3, -3], [63.5, 63.5, 70, 70, 63.5] #lllon, lllat, urlon, urlat= -3, 63.5, 15, 69.5 #to display the third lambert map
#boxlon, boxlat= [0, 10, 10, 0, 0], [67.9, 67.9, 74, 74, 68] #to display the area without fluxes

                 
save=False

#savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/PseudoSat_movie/'
if var in ['CTT', 'WVTC', 'WVT', 'CWR']:
    savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/PseudoSat3/'
elif Thorpex==True:
    savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/DropComp2/'
else:
    savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/Fields2/'
#    savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/IntensityFields/'
                 
                 
"""end global variables"""




for t in tlist:
#    for exp in explist:
    print('t', t)
    
    if maptype=='AA': plt.figure(fignr, figsize= (6, 6))
    elif maptype=='AA_half': plt.figure(fignr, figsize= (6, 4.7))
    else: plt.figure(fignr, figsize= (6, 4.5))    
    fignr+=1

    plt.clf()
    
    title_extra=''

    if maptype== 'AA': map = AA_map()
    elif maptype== 'AA_half': map = AA_map_half()
    else: map= Lambert_map(lllon, lllat, urlon, urlat, lat0, lon0)    
    
    
    AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'                                                     
    AAfilename_2= Mediadir+'PL/AA/ec/'+exp_name_2+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'                                                     

    AAres=2 #2- every second datapoint is taken
    d= data(filename= AAfilename, res= AAres)
    d2= data(filename= AAfilename_2, res= AAres)

    
    Lon,Lat = map(d.lon,d.lat)
        
        
    d.imp_surf(tn= t)
    d2.imp_surf(tn= t)

    
    if var in ['FLX', 'LH', 'LH_thor', 'TH', 'TH_thor', 'PBH', 'CAPE', 'CIN', 'LCL', 'LFC', 'OLW','Precip','Precip_acc'] or 'Cloud' in var:
        d.imp_atm(tn= t)
        d2.imp_atm(tn= t)
    
    if '_lev' in var:
        d.imp_level(pn= pn, tn= t)
        d2.imp_level(pn= pn, tn= t)
    
    if var in ['CTT', 'WVTC', 'WVT', 'CWR']:
        if exp_name in ['08030312_cycling']: t//=3 #for some experiments the pseudo sat is only retrieved every 3 hour
                       
        d.imp_pseudo(tn= t)
        d2.imp_pseudo(tn= t)
        
        if exp_name in ['08030312_cycling']: t*=3
    #var= 'PressWind'
    
    
    """MSLP, Surface temperature and surface winds"""
    if var== 'PressWind':
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotWindVelo(Lon, Lat, np.sqrt(d.u10m**2+ d.v10m**2), map, Umax= 25)
    
    if var== 'PressWindarrow':
        PlotContours(Lon, Lat, d.mslp, map, leveldist= 2)
    #    PlotWindVelo(Lon, Lat, np.sqrt(d.u10m**2+ d.v10m**2), map, Umax= 25)
        PlotWind(d.lon, d.lat, d.u10m, d.v10m, map, rot=False, arrowtype='barb') #put this after PlotColorMap, otherwise it disappears
                 
                 
    
    """MSLP, Surface temperature and surface winds -advanced"""
    if var== 'PressWind_advanced':
        PlotContours(Lon, Lat, d.mslp, map, leveldist= 1, numbers=False)
    
        U= np.sqrt(d.u10m**2+d.v10m**2)
        PlotWindVelo(Lon, Lat, U, map, Umax= 25, color='YlBu')    
        #find local pressure min and plot them in the map   
        PlotLocalMax(d.mslp, threshold=1010, distance=100/AAres, map= map, lon=d.lon, lat=d.lat,
                     data2=U, threshold2=18, distance2=150/AAres)    
      
        PlotLocalMax(U, threshold=20, distance=100/AAres, map= map, lon=d.lon, lat=d.lat, typ='max',
                     color='orange', dot=False, roundorder= 0)    
    
    
    if var== 'PressWindVortmax':
        PlotContours(Lon, Lat, d.mslp, map, leveldist= 2, numbers=False)
    
        U= np.sqrt(d.u10m**2+d.v10m**2)
    #    PlotWindVelo(Lon, Lat, U, map, Umax= 25)    
        #find local pressure min and plot them in the map   
    #    PlotLocalMax(d.mslp, threshold=1010, distance=100/AAres, map= map, lon=d.lon, lat=d.lat,
    #                 data2=U, threshold2=16, distance2=80/AAres)    
        PlotWind(d.lon, d.lat, d.u10m, d.v10m, map, rot=False, arrowtype='barb') #put this after PlotColorMap, otherwise it disappears
    
      
        PlotLocalMax(U, threshold=20, distance=100/AAres, map= map, lon=d.lon, lat=d.lat, typ='max',
                     color='r', dot=False, roundorder= 0)    
    
    
        d.imp_level(pn= 4, tn= t)
        dx= 2500* AAres
        vort= (np.gradient(d.v, dx, axis= 1)- np.gradient(d.u, dx, axis= 0))*1E5         
        #probably one has to rotate u and v first to proper x and y coordinates
        import scipy.ndimage.filters as filters
        gausfilterdist= 60
        vortfilter= filters.gaussian_filter(vort, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 1.)
    #    PlotColorMap4(Lon, Lat, vortfilter, map, label= r"$\xi_{"+str(d.chosen_plevel)+"}$ [10$^{-5}$ 1/s]" #Gauss filter of std "+str(gausfilterdist)+"km"
    #    , color='RdBu', symetric=True)
    
        PlotLocalMax(vortfilter, threshold=10, distance=250/(AAres*2.5), map= map, lon=d.lon, lat=d.lat, typ='max', color='y', roundorder= 0)    
    
    
    """potential temperature at 700 hPa and stability"""
    if var== 'Theta_lev':
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist, numbers= False)
    #    PlotContours(Lon, Lat, d.geop, map, leveldist= 10)
        if d.chosen_plevel == 500: bounds= np.arange(278, 288, 1)
        elif d.chosen_plevel == 700: bounds= np.arange(270, 282, 1) #for pot temp 700
        elif d.chosen_plevel == 850: bounds= np.arange(264, 280, 1) #for pot temp 850
        elif d.chosen_plevel <= 950: bounds= np.arange(260, 278, 1) #for pot temp 950 
        else: bounds= None
        
#        bounds= np.arange(267, 290.5, 1)
#        bounds= np.arange(260, 291, 2)
        
        theta= PotTemp(d.T, d.chosen_plevel)
        bounds, colors=PlotColorMap4(Lon, Lat, theta, map, bounds= bounds, color='RdBu', label=r"$\theta_{"+str(d.chosen_plevel)+"}$ [K]")
        
    #    title_extra= '_static_stability'
    #     
    #    #"""stability"""
    #    d.imp_level(pn= 8, tn= t) # 8 = 500hPa
    #    theta500= PotTemp(d.T, d.chosen_plevel)
    #    thetaSST= PotTemp(d.T0m, d.mslp)  
    #    PlotStaticStability(Lon, Lat, (thetaSST- theta500) , map, difftype='thetaSST-theta500')


    
    """ equivalent potential temperature at 850hPa"""
    if var== 'Theta_e_lev':
        bounds= None
        spechum= d.RH/(.263 *d.pres[pn])* np.exp(17.67*(d.T-273.15) / (d.T-29.65))
        theta_e= EquiPotTemp(d.T, spechum, d.chosen_plevel)
        
        PlotColorMap4(Lon, Lat, theta_e, map, bounds= bounds, label=r"$\theta_{e,"+str(d.chosen_plevel)+"}$ [K]")
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)



    if var== 'Front_lev':
        gausfilterdist = 20 #in km
        spechum= d.RH/(.263 *d.pres[pn])* np.exp(17.67*(d.T-273.15) / (d.T-29.65))
        theta_e= EquiPotTemp(d.T, spechum, d.chosen_plevel)
        
#        PlotColorMap3(Lon, Lat, theta_e, map, bounds= bounds, label=r"$\theta_{e,"+str(d.chosen_plevel)+"}$ [K]")
        PlotContours(Lon, Lat, theta_e, map, leveldist= pleveldist)
        
        dx= 2500* AAres
        grad_theta_e= np.gradient(theta_e, dx, edge_order=1)
        abs_grad_theta_e= np.sqrt(grad_theta_e[0]**2+ grad_theta_e[1]**2)
        abs_grad_theta_e_filter= filters.gaussian_filter(abs_grad_theta_e, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 4.)

#        PlotColorMap4(Lon, Lat, abs_grad_theta_e_filter *1E3*1E2, map, color= 'red', bounds=np.arange(0,10.1,1), label=r"$\nabla \theta_{e,"+str(d.chosen_plevel)+"}$ [K/100km]")
        
        TFP= - np.sum(np.array(np.gradient(abs_grad_theta_e, dx, edge_order=1)) * grad_theta_e/abs_grad_theta_e  , axis= 0)
        TFP_filter= filters.gaussian_filter(TFP, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 4.)

#        PlotColorMap4(Lon, Lat, TFP_filter, map, color= 'RdBu', symetric=True, label=r"TFP$_{"+str(d.chosen_plevel)+"}$")

        front= np.zeros(np.shape(theta_e))
        front[np.logical_and(np.abs(TFP_filter)< .3* 1E-9, abs_grad_theta_e_filter >5*1E-5)]= 1
             
        PlotColorMap4(Lon, Lat, front, map, color= 'red')#, label=r"$\nabla \theta_{e,"+str(d.chosen_plevel)+"}$ [K/100km]")


    if var== 'Grad_Theta_e_lev':
        gausfilterdist = 20 #in km
        spechum= d.RH/(.263 *d.pres[pn])* np.exp(17.67*(d.T-273.15) / (d.T-29.65))
        theta_e= EquiPotTemp(d.T, spechum, d.chosen_plevel)
        
#        PlotColorMap3(Lon, Lat, theta_e, map, bounds= bounds, label=r"$\theta_{e,"+str(d.chosen_plevel)+"}$ [K]")
        PlotContours(Lon, Lat, theta_e, map, leveldist= 2)
        
        dx= 2500* AAres
        grad_theta_e= np.gradient(theta_e, dx, edge_order=1)
        abs_grad_theta_e= np.sqrt(grad_theta_e[0]**2+ grad_theta_e[1]**2)
        abs_grad_theta_e_filter= filters.gaussian_filter(abs_grad_theta_e, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 4.)

        PlotColorMap4(Lon, Lat, abs_grad_theta_e_filter *1E3*1E2, map, color= 'red', bounds=np.arange(0,10.1,1), label=r"$\nabla \theta_{e,"+str(d.chosen_plevel)+"}$ [K/100km]")
        
       
    """mean(SLP) - SLP"""
    if var== 'MeanSLPminSLP':  # make AAres= 10 #otherwise it becomes really uneffecient
        radius= 300
        pixrad= radius/(2.5*AAres) #the pixelradius
        SLP_radmean= CircularFilter_samedist(d.mslp, radius= pixrad)
        diffslp= d.mslp - SLP_radmean
        diffslp[:int(.5*pixrad)], diffslp[:, :int(.5*pixrad)], diffslp[-int(.5*pixrad):], diffslp[:, -int(.5*pixrad):]= 0, 0, 0, 0  #cout of edges
    
        PlotColorMap4(Lon, Lat, diffslp, map, label=r'$\overline{SLP} - SLP$ [hPa] in '+str(radius)+'km radius', color='RdBu', symetric=True)
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        title_extra='_radius'+str(radius)+'km'
    
    """filtered SLP"""
    if var== 'SLP_filtered':
        T_low, T_up= 60, 120 
        MSLPfourierfilter= FourierFilter2d_equaldist(d.mslp, dist= 2.5* AAres, T_low= T_low, T_up= T_up)
    
        PlotColorMap4(Lon, Lat, MSLPfourierfilter, map, label= r"SLP filtered T"+str(T_low)+'-T'+str(T_up) , color='RdBu', symetric=True)
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        title_extra='_T'+str(T_low)+'-T'+str(T_up)
    
    
    """filtered Laplacian of SLP"""
    if var== 'SLP_laplacian':
        dx= 2500* AAres
        SLPlap= np.gradient(np.gradient(d.mslp*1E2, dx, axis= 0, edge_order=1), dx, axis= 0, edge_order=1) + np.gradient(np.gradient(d.mslp*1E2, dx, axis= 1, edge_order=1), dx, axis= 1, edge_order=1)
        SLPlap[:3], SLPlap[:, :3], SLPlap[-3:], SLPlap[:, -3:]= 0, 0, 0, 0  #cout of edges
    
        import scipy.ndimage.filters as filters
        gausfilterdist= 100  
        SLPlapfilter= filters.gaussian_filter(SLPlap, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 4.)
    
        PlotColorMap4(Lon, Lat, SLPlapfilter, map, color= 'RdBu', symetric= True, label= r"$\nabla^2$ SLP [Pa/m$^2$] Gauss filter of std "+str(gausfilterdist)+"km")
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        title_extra= '_GaussFilter'+str(gausfilterdist)+'km'
    
    
    """calculate vorticity"""
    if var== 'Vort_lev':
        dx= 2500* AAres
        vort= (np.gradient(d.v, dx, axis= 1)- np.gradient(d.u, dx, axis= 0))*1E5         
        #probably one has to rotate u and v first to proper x and y coordinates
        import scipy.ndimage.filters as filters
        gausfilterdist= 100# 40 #100 #60
        gaustrauncate=2
        vortfilter= filters.gaussian_filter(vort, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= gaustruncate) #truncate = 1
        PlotColorMap4(Lon, Lat, vortfilter, map, label= r"$\xi_{"+str(d.chosen_plevel)+"}$ [10$^{-5}$ 1/s]" #Gauss filter of std "+str(gausfilterdist)+"km"
        , color='RdBu', symetric=True)
        
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist, numbers=False)
        U= np.sqrt(d.u10m**2+d.v10m**2)
    
        PlotLocalMax(vortfilter, threshold=10, distance=250/(AAres*2.5), map= map, lon=d.lon, lat=d.lat, typ='max', color='black', roundorder= 0)    
    #    PlotLocalMax(d.mslp, threshold=1010, distance=100/AAres, map= map, lon=d.lon, lat=d.lat,
    #                 data2=U, threshold2=16, distance2=80/AAres, color='g')      
        title_extra= '_GaussFilter'+str(gausfilterdist)+'km'
    
    
    
    """calculate shear deformation"""
    if var== 'Shear_lev':
        dx= 2500* AAres
        vort= (np.gradient(d.v, dx, axis= 1)+ np.gradient(d.u, dx, axis= 0))*1E5 #it is not transformation invariant
    #    vort= np.sqrt(np.gradient(d.v, dx, axis= 1)**2+ np.gradient(d.u, dx, axis= 0)**2)*1E5
    
    
        import scipy.ndimage.filters as filters
        gausfilterdist= 60
        vortfilter= filters.gaussian_filter(vort, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 1.)
        PlotColorMap4(Lon, Lat, vortfilter, map, label= r"$(\partial_x v + \partial_y u)_{"+str(d.chosen_plevel)+"}$ [10$^{-5}$ 1/s] Gauss filter of std "+str(gausfilterdist)+"km", color='RdBu', symetric=True)
    #    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotContours(Lon, Lat, d.geop, map, leveldist= 10)
        PlotWind(d.lon, d.lat, d.u, d.v, map, rot=False)
    
        title_extra= '_GaussFilter'+str(gausfilterdist)+'km'
    
    
    """calculate stretching deformation - this looks more like like the shear"""
    if var== 'Stretching_lev':
        dx= 2500* AAres
        vort= (np.gradient(d.u, dx, axis= 1)- np.gradient(d.v, dx, axis= 0))*1E5
              
    #    vort= (np.gradient(d.u, dx, axis= 0)) *1E5 #- np.gradient(d.v, dx, axis= 0))*1E5 #to look only at the u shear
    #    vort= (np.gradient(d.v, dx, axis= 1)) *1E5 #- np.gradient(d.v, dx, axis= 0))*1E5
    
        import scipy.ndimage.filters as filters
        gausfilterdist= 60
        vortfilter= filters.gaussian_filter(vort, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 1.)
        PlotColorMap4(Lon, Lat, vortfilter, map, label= r"$(\partial_x u - \partial_y v)_{"+str(d.chosen_plevel)+"}$ [10$^{-5}$ 1/s] Gauss filter of std "+str(gausfilterdist)+"km", color='RdBu', symetric=True)
    #    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotContours(Lon, Lat, d.geop, map, leveldist= 10)
        PlotWind(d.lon, d.lat, d.u, d.v, map, rot=False)
        
    #    PlotWind(d.lon, d.lat, d.u, np.zeros(d.v.shape), map, rot=False, nx= 50) #to look only at the u shear
    #    PlotWind(d.lon, d.lat, np.zeros(d.v.shape), d.v, map, rot=False, ny= 50)
           
        title_extra= '_GaussFilter'+str(gausfilterdist)+'km'
    
    """calculate Okubo-Weiss parameter (vort^2 - (shear^2 + stretching^2)"""
    if var== 'OW_lev':
        dx= 2500* AAres
        vort= 1E10* ( (np.gradient(d.v, dx, axis= 1)- np.gradient(d.u, dx, axis= 0))**2 - (np.gradient(d.v, dx, axis= 1)+ np.gradient(d.u, dx, axis= 0))**2 - (np.gradient(d.u, dx, axis= 1)- np.gradient(d.v, dx, axis= 0))**2 )
        import scipy.ndimage.filters as filters
        gausfilterdist= 60
        vortfilter= filters.gaussian_filter(vort, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 1.)
        vortfilter[vortfilter < 0]= 0
    #    PlotColorMap4(Lon, Lat, vortfilter, map, label= r"OW$_{"+str(d.chosen_plevel)+"}$ [10$^{-8}$ 1/s$^2$] Gauss filter of std "+str(gausfilterdist)+"km", color='RdBu', symetric=True)
        PlotColorMap4(Lon, Lat, np.sqrt(vortfilter), map, label= r"$\sqrt{OW_{"+str(d.chosen_plevel)+"}}$ [10$^{-5}$ 1/s] Gauss filter of std "+str(gausfilterdist)+"km", color='red')
    
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        title_extra= '_GaussFilter'+str(gausfilterdist)+'km'
    
        
    """the fourier filtered vorticity"""
    if var== 'Vort_lev_T':
    #    d.imp_level(pn= 4, tn= t) #850hPa
        
        dx= 2500* AAres
        vort= (np.gradient(d.v, dx, axis= 1)- np.gradient(d.u, dx, axis= 0))*1E5
        
        T_low, T_up= 40, 100          
        vortfourierfilter= FourierFilter2d_equaldist(vort, dist= 2.5* AAres, T_low= T_low, T_up= T_up)
    
        PlotColorMap4(Lon, Lat, vortfourierfilter, map, label= r"$\xi_{"+str(d.chosen_plevel)+"}$ [10$^{-5}$ 1/s] filtered T"+str(T_low)+'-T'+str(T_up) , color='RdBu', symetric=True)
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        title_extra='_filtered_T'+str(T_low)+'-T'+str(T_up)
    
    """ different geopotential levels and mslp"""
    if var== 'Geop_lev':
    #    title_extra=' SLP and '+str(d.chosen_plevel)+'hPa geopotential height'
        PlotContours(Lon, Lat, d.mslp, map, leveldist= 2)
        PlotContours(Lon, Lat, d.geop, map, color='r', leveldist= 10)
    #    PlotContours(Lonm, Latm, geop500, m, color=[(plt.cm.RdBu_r(h)) for h in range(0,257,32)], leveldist= 20)
    
    
    
    if var== 'test_rot': #the rotation of the wind field
        del_lon= (d.lon[1:, :]- d.lon[:-1, :])* np.cos(np.deg2rad(d.lat[1:, :]))
        del_lat= d.lat[1:, :]- d.lat[:-1, :]
        
        rot= np.rad2deg(np.arctan2(del_lon, del_lat)) #how much to rotate the gridcell anticlockwise to be oriented N-S
    #    PlotColorMap4(Lon[1:], Lat[1:], rot, map)
        
        dir_raw= np.rad2deg(np.arctan2(d.u[1:], d.v[1:])) #the wind direction in respect to the grid (0 being upwards, 90 rightwards, -90 leftw)
    #    PlotColorMap4(Lon[1:], Lat[1:], dir_raw, map)
    
        dir_rot= (dir_raw + rot)%360 #the wind direction with respect to the North pole (0 Northwards, 90 E, 270 W)
    #    PlotColorMap4(Lon[1:], Lat[1:], dir_rot, map)
    
        wind_vel= np.sqrt(d.u[1:]**2 +d.v[1:]**2)
        u= wind_vel* np.sin(np.deg2rad(dir_rot))
        v= wind_vel* np.cos(np.deg2rad(dir_rot))
        
        PlotWind(d.lon[1:], d.lat[1:], u, v, map, rot=True, nx= 50, ny= 50) #this should give the same as only the wind field
    
    
    if var == 'Geop_wind_lev':
        PlotContours(Lon, Lat, d.geop, map, leveldist= 10, alpha= .7, numbers= False)
    #    PlotWindVelo(Lon, Lat, np.sqrt(d.u**2+ d.v**2), map)
        bounds, colors= PlotColorMap4(Lon, Lat, np.sqrt(d.u**2+ d.v**2), map, color='YlBu', label='Wind velocity [m/s] at '+str(d.chosen_plevel)+'hPa', bounds= np.arange(0, 31, 3))


        uzonal, vmeri= Wind2ZonalMeridional(d.lon, d.lat, d.u, d.v) #this rotates the wind in model coordinate direction to zonal and meridional winds
        PlotWind2(d.lon, d.lat, uzonal, vmeri, map, everyx= 13, everyy= 13, rot=False, color='white', quiverkey= False)
    
#        PlotWind2(d.lon, d.lat, d.u, d.v, map, rot=False, everyx= 13, everyy= 13, color='white') #this was a wrong rotation
    
    
    """filtered Laplacian of the geopotential height"""
    if var== 'Geop_laplacian_lev':
        dx= 2500* AAres
        geoplap= np.gradient(np.gradient(d.geop, dx, axis= 0, edge_order=1), dx, axis= 0, edge_order=1) + np.gradient(np.gradient(d.geop, dx, axis= 1, edge_order=1), dx, axis= 1, edge_order=1)
        geoplap[:3], geoplap[:, :3], geoplap[-3:], geoplap[:, -3:]= 0, 0, 0, 0  #cout of edges
    
        import scipy.ndimage.filters as filters
        gausfilterdist= 40  
        geoplapfilter= filters.gaussian_filter(geoplap, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 4.)
    
        PlotColorMap4(Lon, Lat, geoplapfilter, map, color= 'RdBu', symetric= True, label= r"$\nabla^2$ Z$_{"+str(d.chosen_plevel)+"}$ [m/m$^2$] Gauss filter of std "+str(gausfilterdist)+"km")
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist, numbers=False)
        PlotContours(Lon, Lat, d.geop, map, color='r', leveldist= 10, numbers=False)
    
        title_extra= '_GaussFilter'+str(gausfilterdist)+'km'
     
    """geopotential 1000 -500 hPa thickness"""
    if var== 'Thickness':
        d.imp_geop1000(tn= t)
        d.imp_level(pn= 8, tn= t)
    
        geop_diff= d.geop - d.geop1000
        PlotContours(Lon, Lat, geop_diff, map)
        PlotColorMap3(Lon, Lat, geop_diff, map, symetric=False, label='1000-500 hPa thickness')
        title_extra='_1000-500 hPa thickness'
                                      
       
    """ Precip """
    if var== 'Precip': 
        precipbounds= [0, 0.1, 0.2, 0.5, 1, 2, 5, 10]
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.prec, map, bounds= precipbounds, color='blue', label='Precipitation [mm/h]')

    if var== 'Precip_acc':
#        precipbounds= [0, 0.1, 0.2, 0.5, 1, 2, 5, 10]
#        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
#        PlotColorMap4(Lon, Lat, d2.prec_acc - d.prec_acc, map, color='RdBu', symetric=True, label='Precipitation in '+str(t)+'h [mm]')

        gausfilterdist= 30  
        prec_diff_filter= filters.gaussian_filter(d2.prec_acc - d.prec_acc, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= 4.)
        PlotColorMap4(Lon, Lat, prec_diff_filter, map, color='RdBu', symetric=True, label='Precipitation in '+str(t)+'h [mm]')

    
    """ Sensible Heat Flux"""
    if var== 'TH':
        fluxbounds= [-800,  -500, -300, -200 , -100, -50, 0, 50, 100, 200, 300, 500, 800]
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.SenH, map, color='RdBuWhite', bounds= fluxbounds, label=r"Sensible heat flux [W/m$^2$]")

    if var== 'TH_thor':
        fluxbounds= np.arange(0, 280, 2)
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.SenH, map, color='jet', bounds= fluxbounds, label=r"Sensible heat flux [W/m$^2$]")
       
    
    if var== 'LH':
        fluxbounds= [-800,  -500, -300, -200 , -100, -50, 50, 100, 200, 300, 500, 800]
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.LatH, map, color='RdBu', bounds= fluxbounds, label=r"Latent heat flux [W/m$^2$]")

    if var== 'LH_thor':
        fluxbounds= np.arange(0, 280, 2)
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.LatH, map, color='jet', bounds= fluxbounds, label=r"Latent heat flux [W/m$^2$]")
       
    
    if var== 'FLX':
#        fluxbounds= [0,  100, 200, 300, 500, 800, 1200]
#        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d2.SenH+ d2.LatH - (d.SenH+ d.LatH), map, symetric=True, color= 'RdBuWhite', label=r"Turbulent heat flux [W/m$^2$]")
       
    
    if var== 'PBH':
        BLbounds= [0,  500, 1000, 2000, 3000, 4000, 5000]
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.BL, map, bounds=BLbounds,  color= 'blue', label=r"Boundary layer height [m]")
       
    
    if var== 'W_lev':
    #    d.imp_level(pn= 6, tn= t)
#        wbounds= [-5, -1, -0.5, -0.2, -0.1, 0.1, 0.2, 0.5, 1, 5]
        wbounds= [-5, -2, -1, -0.5, -0.2, 0.2, 0.5, 1, 2, 5]
        print('max updraft: ', str(np.max(d.w)))
        
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.w, map, bounds=wbounds, color= 'RdBuWhite', label= r"vertical velocity at "+str(d.chosen_plevel)+" [m/s]")
    
    
    
    if var== 'PV_lev':
        PVbounds= np.array([-5, -3, -2, -1, -0.5, 0,  0.5, 1, 2, 3, 5])
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.PV *1E6, map, bounds=PVbounds, color= 'RdBu', label= str(d.chosen_plevel) + r"hPa PV [$10^{-6}$K m$^2$ kg$^{-1}$ s$^{-1}$]")
    
    
    if var== 'Tsurf':
        Tsurfbounds= np.array([-40, -20, -10, -5, -2, 0,  2,  5, 10, 20, 40])  
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.T0m -273.15 , map, bounds= Tsurfbounds, color= 'RdBu', label=r"Surface temperature [K]")
            
    
    if var== 'T2m':
        Tsurfbounds= np.array([-40, -20, -10, -5, -2, 0,  2,  5, 10, 20, 40])  
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.T2m -273.15 , map, bounds= Tsurfbounds, color= 'RdBu', label=r"2m temperature [K]")
          
    
    if var== 'RH_lev':
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        bounds, colors= PlotColorMap4(Lon, Lat, d.RH, map, symetric=False, label='Relative humidity at ' +str(d.chosen_plevel)+'hPa', color='blue', bounds= np.array([0, 0.4, 0.7, 0.8, 0.9, 0.95, 1]) )
    
    if var== 'SH_lev':
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist, numbers= False)
        bounds, colors= PlotColorMap4(Lon, Lat, d.SH*1E3, map, symetric=False, label='Specific humidity [g/kg] at ' +str(d.chosen_plevel)+'hPa', color='blue', bounds= np.arange(0, np.max(d.SH)*1E4, 5)*1E-1)
    
    
    if var== 'RH_lev_fromSH':
        import metpy.calc as mpcalc
        from metpy.units import units
#        RH= mpcalc.relative_humidity_from_specific_humidity(d.SH, d.T* units.K, d.chosen_plevel* units.hPa)
        RH= SH2RH(d.SH*1E3, d.chosen_plevel, d.T)/100
  
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist, numbers= False)
        bounds, colors= PlotColorMap4(Lon, Lat, RH, map, symetric=False, label='Relative humidity at ' +str(d.chosen_plevel)+'hPa', color='blue', bounds= np.array([0, 0.4, 0.7, 0.8, 0.9, 0.95, 1]) )
    
    
    if var== 'Cloud_conv':
    #    Cloud_convbounds= np.array([0,  100,  200, 300, 500, 800])
        Cloud_bounds= np.array([0, 0.2, 0.4, 0.6, 0.8, 1])
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.cloud_conv , map, bounds=Cloud_bounds, color= 'brown', label=r"Convective Clouds")
    
    if var== 'Cloud_high':
        Cloud_bounds= np.array([0, 0.2, 0.4, 0.6, 0.8, 1])    
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.cloud_high , map, bounds=Cloud_bounds, color= 'grey', label=r"High cloud cover")
    
    if var== 'Cloud_med':
        Cloud_bounds= np.array([0, 0.2, 0.4, 0.6, 0.8, 1])    
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.cloud_med , map, bounds=Cloud_bounds, color= 'grey', label=r"Medium cloud cover") 
    
    
    if var== 'Cloud_low':
        Cloud_bounds= np.array([0, 0.2, 0.4, 0.6, 0.8, 1])    
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.cloud_low , map, bounds=Cloud_bounds, color= 'grey', label=r"Low cloud cover")
    
    if var== 'CAPE':
        CAPEbounds= np.array([0,  100,  200, 300, 500, 800])    
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.CAPE , map, bounds=CAPEbounds, color= 'brown', label=r"CAPE [J/kg]")
    
    if var== 'CIN':
        CINbounds= np.array([-100, -10, -1, -0.1, 0])    
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.CIN , map, bounds=CINbounds, color= 'brown', label=r"CIN [J/kg]")
     
        
    if var== 'LFC':
        LFCbounds= np.array([0, 10, 100, 200, 500, 1000, 10000, 30000])    
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.LFC , map, bounds=LFCbounds, color= 'brown', label=r"LFC [m]")
    
    if var== 'LCL':
        LCLbounds= np.array([0, 10, 100, 1000, 3000, 5000, 10000, 30000])  
        PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.LCL , map, bounds=LCLbounds, color= 'brown', label=r"LCL [m]")
    
    if var== 'CTT':
        levels= np.arange(940, 1040, 3)
        PlotContours(Lon, Lat, d.mslp, map, levels= levels, numbers=False, color= 'r', alpha= .6)
        bounds= np.arange(210, 281, 2)   
    
        PlotColorMap4(Lon, Lat, d.CTT , map, color= 'grey', bounds= bounds, label=r"Cloud top temperature [K]")
    #    PlotColorMap4(Lon, Lat, d.CTT , map, label=r"Cloud top temperature [K]", boxnr= 50)
    
    if var== 'WVTC':
    #    Cloud_bounds= np.array([0, 0.2, 0.4, 0.6, 0.8, 1])    
    #    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.WVTC , map, color= 'grey', label=r"Water vapor temperature with clouds [K]", boxnr= 128)
    
    if var== 'WVT':
    #    Cloud_bounds= np.array([0, 0.2, 0.4, 0.6, 0.8, 1])    
    #    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.WVT , map, color= 'grey', label=r"Water vapor temperature [K]", boxnr= 100)
    
    if var== 'CWR':
    #    Cloud_bounds= np.array([0, 0.2, 0.4, 0.6, 0.8, 1])    
    #    PlotContours(Lon, Lat, d.mslp, map, leveldist= pleveldist)
        PlotColorMap4(Lon, Lat, d.CWR , map, color= 'grey', label=r"Cloud water reflectivity", bounds= np.linspace(0,1., 101))

    if var =='OLW':
        levels= np.arange(940, 1040, 3)
        PlotContours(Lon, Lat, d.mslp, map, levels= levels, numbers=False, color= 'r', alpha= .6)
        
        PlotColorMap4(Lon, Lat, d.OLW , map, label= 'Outgoing longwave radiation [W/m$^2$]', color='grey')
        
    
    """Plot Thorpex points"""
    if Thorpex== True:
        if year == 2008 and month== 3 and day== 3 and hour<15:
            excl= [3, 5, 12] #for the first flight       
        elif year == 2008 and month== 3 and day== 3 and hour>=15:
            excl=[1,9,11] #for the second flight
        elif year == 2008 and month== 3 and day== 4:
#            excl= [1, 5, 13]# for the third flight - for AA domain moved southwards
            excl= [1, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]# for the third flight - for original AA domain

        
        thor= thorpex.data(year, month, day, hour, plevels= [d.chosen_plevel], exclude=excl)
    #    thor2= thorpex.data(year, month, day, hour, level= d.chosen_plevel)
        
        xpt, ypt= map(thor.lon, thor.lat)
        #map.plot(xpt,ypt,'bo')
        
        for i in np.arange(len(xpt)):
            scattercolor= 'black'
            if var=='Geop_wind_lev': scattercolor= 'r'
#            if var is not 'Geop_wind_lev':
            plt.text(xpt[i]-10000, ypt[i]+20000, str(thor.dropnr[i]), color=scattercolor, fontsize=14)               
#                plt.text(xpt[i]+100, ypt[i]+100, str(thor.dropnr[i])+'\n'+ str(thor.datetime[i, 0])[11:16], color='r')
                    
            
            if var== 'Theta_lev' or var=='Theta_newbounds_lev':
                value= thor.theta[i]
            elif 'RH_lev' in var:
                value= thor.RH[i]
            elif var== 'SH_lev':
                value= thor.SH[i]*1E3
            elif var== 'Geop_wind_lev':
                value= thor.U[i]                
            else:
                break

            if np.isnan(value[0]): print('Exclude ', thor.dropnr[i])   #in case the value that should be plotted is none it exits the loop
            else:
#                a= value - bounds
#                ind= np.where(a== min(i for i in a if i >= 0))[0][0]+1 #this goes wrong if the thorpex value is outside the bounds
                a= np.abs(value - bounds) # in order to find the color of the point, it calculates the distance of the value to the bounds of the colormap
                ind= np.where(a== np.min(a))[0][0] #finds the index of bounds that should be used as color for the point

                colornow= colors[ind]
        
        #        plt.text(xpt[i]+100, ypt[i]+100, str(np.round(value[0],2) ) )
                scattercolor= 'black'
                if var=='Geop_wind_lev': scattercolor= 'r'
                
                map.scatter(xpt[i],ypt[i], color= colornow, s= 150, edgecolors= scattercolor)
        
        
                if var== 'Geop_wind_lev':
                    dist= np.sqrt( (thor.lat[i,0]- d.lat)**2+ (np.cos(np.deg2rad(thor.lat[i,0]))* (thor.lon[i,0]- d.lon))**2)
                    xpos, ypos= np.where(dist== np.min(dist))
#                    plt.text(xpt[i], ypt[i]+10000, str(int(thor.alt[i,0]- d.geop[xpos[0], ypos[0]])), color='r')  #depature from the geopotential height  
                    
                    print(thor.dropnr[i], str(int(thor.alt[i,0]- d.geop[xpos[0], ypos[0]])))        
                    u,v, Lon, Lat = map.rotate_vector(thor.u[i], thor.v[i] , thor.lon[i], thor.lat[i] ,returnxy=True)
                    plt.quiver(Lon, Lat, u, v, color= 'r', scale = 500, width= .004, headwidth= 5)
    
    
    
    if Plotbox== True:
        x,y= map(boxlon, boxlat) #longitude and latitude coordinates of the corner points
        map.plot(x, y, color='b', linewidth= 2)
        title_extra +='_box'
    
    
    if save== False:
        plt.title('Arome '+str(d.datetime[t])[:-3]+ title_extra +'\n' + exp_name_2 +' - '+ exp_name)
    
    plt.tight_layout()
    
    if save==True:
        if '_lev' in var:
            var += str(d.chosen_plevel)
            
        print('exp_name:', exp_name)
        if exp_name== 'DA_080301_cycling': save_exp_name= exp_name+'_'+str(fileday).zfill(2)+str(filehour).zfill(2)
        else: save_exp_name= exp_name
        
        save_exp_name= exp_name_2+'-'+save_exp_name
        
        savename= savedir+'Arome_'+save_exp_name+'_'+str(t).zfill(2)+'_'+ var +title_extra
        plt.savefig(savename)
        print(savename)







    


    

"""for cross section definition"""
#x_start, y_start= 10, 30
#x_end, y_end= 110, 20
#
#x, y= mapcross_section(d.lon, d.lat, x_start, y_start, x_end, y_end, m)
