#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:07:08 2017

@author: pst019

This is from plot_comparison
"""

import time
start = time.time()


import os
user = os.getcwd().split('/')[2]

homedir= '/home/'+user+'/home/'
Mediadir= '/media/'+user+'/PatsOrange/'



from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})

import numpy as np
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
import matplotlib.dates as mdates
import scipy.ndimage.filters as filters

#import itertools

plt.rcParams.update({'font.size': 12})




"""global variables"""
fignr= 6
                 

"""b) times"""
year, month = 2008, 3
day, hour= 3, 6 #22                   


"""Lambert coordinates"""
#maptype='AA'
maptype='AA_half'
#maptype='AA_south'


fileday, filehour= 3, 0        
t= (day- fileday)*24 + (hour- filehour)  # -1 
 
    
"""name of the experiment"""
exp_name, label_name= '080303_warmctr', 'CTR'
#exp_name='080303_cold_sice'
#exp_name= '080303_warmsens_noQH'
#exp_name, label_name='080303_warmsens_2FLX', '2FLX'
#exp_name='080303_warmsens_nocondens'
#exp_name='870226_cold_ERA_can'
#exp_name, label_name= '080303_cold_ERA', 'COLD-ERA'
#exp_name, label_name ='080303_cold_pseudo2', 'Cold'
#exp_name, label_name='080303_coldsens_noFLX_AREA', 'noFLX-Area'
#exp_name, label_name= '080303_cold', 'COLD'
#exp_name= '080303_warmsens_notopo'
#exp_name, label_name='080304_cold_pseudo', 'Cold_04'


exp_name, label_name= 'DA_080303_CTR', 'CTR'
#exp_name, label_name= 'DA_080303_noCondens', 'noCondens'
#exp_name, label_name= 'DA_080303_2FLX', '2FLX'
#exp_name, label_name= 'DA_080303_noQH', 'noQH'


#exp_name, label_name= 'DA_080303_noFLXarea', 'noFLXarea'
#exp_name, label_name= 'DA_080303_p6SST', '+6SST'

#exp_name, label_name= 'DA_080301_cycling', 'cycling'
#exp_name, label_name= 'DA_08030212', 'early_start'



explist=[exp_name]
legendlabel= [label_name]


colortype= 1 #use the python default colors - might be specified differently below


#explist= ['DA_080303_CTR', 'DA_080303_2FLX']
#legendlabel= ['CTR', '2FLX']

#explist= ['DA_080303_CTR', 'DA_080303_p6SST']
#legendlabel= ['CTR', '+6SST']

#explist= ['DA_080303_CTR', 'DA_080303_noFLXarea', 'DA_080303_noFLX', 'DA_080303_2FLX']
#legendlabel= ['CTR', 'noFLXarea', 'noFLX', '2FLX']
#
#explist= ['DA_080303_CTR', 'DA_080303_noFLXarea']
#legendlabel= ['CTR', 'noFLXarea']


#explist= ['DA_080303_p6SST', 'DA_080303_p4SST', 'DA_080303_p2SST', 'DA_080303_CTR', 'DA_080303_m2SST', 'DA_080303_m4SST', 'DA_080303_m6SST']
#legendlabel= ['+6SST', '+4SST', '+2SST', 'CTR', '-2SST', '-4SST', '-6SST']
#colortype= 2

#explist= ['DA_080303_CTR', 'DA_080303_noQH', 'DA_080303_noCondens', 'DA_080303_noTH', 'DA_080303_noFLX']
#legendlabel= ['CTR', 'noQH', 'noCondens', 'noTH', 'noFLX']

explist= ['DA_080303_CTR', 'DA_080303_noQH', 'DA_080303_noCondens', 'DA_080303_noTH', 'DA_080303_noFLX', 'DA_080303_noFLXarea', 'DA_080303_2FLX']
legendlabel= ['CTR', 'noQH', 'noCond', 'noTH', 'noFLX', 'noFLX-A', '2FLX']
colortype= 3

#explist= ['DA_080303_CTR', 'DA_080303_noHARATU', 'DA_080303_noOCND2']
#legendlabel= ['CTR', 'noHaratu', 'noOCND2']

#explist= ['DA_080303_CTR', 'DA_080303_noQH', 'DA_080303_noCondens']
#legendlabel= ['CTR', 'noQH', 'noCondens']

#explist= ['DA_080303_m6SST', 'DA_080303_noFLX']
#legendlabel= ['m6SST', 'noFLX']


"""variabel for the maplot. Heat, Baro have to be 1 if they are chosen here"""

mapvar= 'PressWindVort'
#mapvar= 'Baroclin_lev'
#mapvar='Heat'


""" colors and markers"""
marker = ['-', '--', ':', '-.', (0, (2,2)), (0, (3, 2, 1, 2, 1, 2)), (0, (4, 2, 4, 2, 1, 2))]

alldiffmarker= False
maxnri = 2 #maximum number of cyclones for one experiment in the intensity graph



if colortype==2:
    colorlist= [plt.cm.Spectral(h) for h in [0, 40, 80]] + ['k'] + [plt.cm.Spectral(h) for h in [170, 210, 250]]     
#    colorlist= [plt.cm.Spectral(h) for h in range(0, 256, 256//8)] #for specified colors (colortype= 2); 8 is chosen for 7 lines since yellow is skipped

if colortype==3:
    colorlist= ['k','#17becf','#1f77b4','#2ca02c','#d62728','#ff7f0e','#8c564b']

"""variables shown in the evolution - main plot"""
Baro= 1 #[0,1]
Baroclin= 'dry' #[dry, moist]

Heat= 0
Heatvar= 'Condens'

number_subplots= 3+ Baro + Heat

mainfig= plt.figure(fignr, figsize= (7,13))
fignr+=1
plt.clf()
                  
ax1= plt.subplot(number_subplots, 1, 1)    
ax2= plt.subplot(number_subplots, 1, 2)  
ax3= plt.subplot(number_subplots, 1, 3)
axlist= [ax1, ax2, ax3]
                     
if Baro== 1:
    ax4= plt.subplot(number_subplots, 1, 4)
    axlist += [ax4]
if Heat== 1:
    ax5= plt.subplot(number_subplots, 1, 5)   
    axlist += [ax5]
    



""" presure level (if the variable ends with '_var' """
pn = 4
#1-950, 4- 850hPa, 6- 700hPa, 8 - 500hPa

save=True
#save=False
savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/IntensityFieldsTracks4/'

title_extra=''

"""end global variables"""
AAres=5 #2- every second datapoint is taken

gausfilterdist= 100 #60 #40
gaustruncate= 2
threshold= 8 #10 #13
propdistance= 180E3 #in km   #70
distance= int(200E3/(AAres * 2.5E3)) #100   #make distance reasonable larger than propdistance to avoid "doublematches" twice is on the safe side


latbound= [60, 71]
lonbound= [-5, 10]
timebound= [20, 30]

#for 4 March:
if fileday== 4: timebound= [0, 6]

winddist= 400/(AAres *2.5) #150
slpdist= 100/(AAres *2.5)


grad_thetadist= 400/(AAres *2.5)
theta_gausfilterdist= 100




for iexp,exp_name in enumerate(explist): 
    if 'noFLXarea' in exp_name:
        Plotbox=True
        boxlon, boxlat= [0, 10, 10, 0, 0], [67.9, 67.9, 74, 74, 68] #to display the area without fluxes
    else: Plotbox=False
    
    """import data"""
    AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
                                                     
    d= data(filename= AAfilename, res= AAres)     
        
    d.imp_uv_ft(pn= pn, additional='all')

    d.imp_surf_ft()

    U= np.sqrt(d.u10m_ft**2+d.v10m_ft**2)

         
    """calculate vorticity"""
    dx= 2500* AAres
    
    vort_full= (np.gradient(d.v_ft, dx, axis= -1)- np.gradient(d.u_ft, dx, axis= -2))*1E5     
    
    vortfilter_full= filters.gaussian_filter1d(vort_full, axis= -1, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= gaustruncate)
    vortfilter_full= filters.gaussian_filter1d(vortfilter_full, axis= -2, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= gaustruncate)
 

    """calculate Grad Theta e"""
    if Baroclin=='moist':
        theta_e= EquiPotTemp(d.T_ft, d.SH_ft, d.chosen_plevel)
        bounds, Baro_label= np.arange(0,10.1,1), r"$\nabla \theta_{e,"+str(d.chosen_plevel)+"}$ [K/100km]"
    elif Baroclin=='dry':
        theta_e= PotTemp(d.T_ft, d.chosen_plevel)
        bounds, Baro_label= np.arange(0,4.6,.5), r"$\nabla \theta_{"+str(d.chosen_plevel)+"}$ [K/100km]"
    
    baroclin_filter= filters.gaussian_filter1d(theta_e, axis=-1, sigma= theta_gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)
    baroclin_filter= filters.gaussian_filter1d(baroclin_filter, axis=-2, sigma= theta_gausfilterdist/(AAres*2.5), mode='nearest', truncate= 1.)

    dx= 2500* AAres
    grad_baroclin_filter= np.gradient(baroclin_filter, dx, edge_order=1)
    abs_grad_baroclin_filter= np.sqrt(grad_baroclin_filter[-1]**2+ grad_baroclin_filter[-2]**2) *1E5 #[K/100km]



    """Plot the data for one experiment"""
    if maptype=='AA': fig= plt.figure(fignr, figsize= (6, 6))  
    elif maptype in ['AA_half', 'AA_south']: fig=plt.figure(fignr, figsize= (6, 4.8))
    fignr+=1
    plt.clf()    

    if maptype=='AA': map = AA_map()
    elif maptype=='AA_half': map = AA_map_half()
    elif maptype== 'AA_south': map = AA_South_map()

    Lon,Lat = map(d.lon,d.lat)

    U= np.sqrt(d.u10m_ft**2+d.v10m_ft**2)

    if mapvar == 'PressWindVort':
        PlotContours(Lon, Lat, d.mslp_ft[t], map, leveldist= 2, numbers=False)   
        PlotWindVelo(Lon, Lat, U[t], map, Umax= 25, color='YlBu')
        
        PlotLocalMax(d.mslp_ft[t], threshold=1010, distance=400/AAres, map= map, lon=d.lon, lat=d.lat,
                     data2=U[t], threshold2=15, distance2=500/AAres, color='k', roundorder= 0, dot= False)         

        d.imp_level(pn= 8, tn= t)
        PlotContours(Lon, Lat, d.geop, map, leveldist= 20, numbers= True, color='red')

    
    elif mapvar == 'Baroclin_lev':#if baroclinicty should be displayed
        PlotContours(Lon, Lat, baroclin_filter[t], map, leveldist= 2)
        PlotColorMap4(Lon, Lat, abs_grad_baroclin_filter[t], map, color= 'red', bounds=np.arange(0,7.1,1), label= Baro_label)

    elif mapvar == 'Heat':
        fluxbounds= [-800,  -500, -300, -200 , -100, -50, 0, 50, 100, 200, 300, 500, 800]       
        PlotContours(Lon, Lat, d.mslp_ft[t], map, leveldist= 2)
#        PlotColorMap4(Lon, Lat, d.SenH_ft[t], map, color='RdBuWhite', bounds=fluxbounds, label=r"Sensible heat flux [W/m$^2$]")
#        PlotColorMap4(Lon, Lat, d.SenH_ft[t], map,  color='blue', label='Sensible heat flux [W/m$^2$]')
        
    if Heat== 1:
        if Heatvar== 'Condens':
            Condens= LatentHeat_fromSnow(d.prec_ft)
            gausfilterdist= 10
            Condens= filters.gaussian_filter1d(Condens, axis=-1, sigma= gausfilterdist/(d.res*2.5), mode='nearest', truncate= 1.)
            Condens= filters.gaussian_filter1d(Condens, axis=-2, sigma= gausfilterdist/(d.res*2.5), mode='nearest', truncate= 1.)
    
            PlotColorMap4(Lon, Lat, Condens[t], map,  color='blue', label='Latent heat release [W/m$^2$]')


    if save== False: plt.title('Arome '+ legendlabel[iexp] + ' '+ str(d.datetime[t])[:-3]+ title_extra)
    
    #find local pressure min and plot them in the map 
#    maxdata, maxlon, maxlat= PlotLocalMax(vortfilter_full[t], threshold=threshold, distance=distance, map= map, lon=d.lon, lat=d.lat, typ='max', color='purple')   
#    PlotLocalMax(d.mslp_ft[t], threshold=1010, distance=100/AAres, map= map, lon=d.lon, lat=d.lat,
#                 data2=U[t], threshold2=16, distance2=80/AAres)         
#    PlotLocalMax(U[t], threshold=20, distance=100/AAres, map= map, lon=d.lon, lat=d.lat, typ='max', color='orange', dot=True, yoffset= -200)    


    if Plotbox== True:
        x,y= map(boxlon, boxlat) #longitude and latitude coordinates of the corner points
        map.plot(x, y, color='b', linewidth= 2)
    
    
    """ a tracking algorithm"""
    vortlist, lonlist, latlist, labellist = EasyTrack(vortfilter_full, d.lon, d.lat, distance, threshold, propdistance)
    vortlist, lonlist, latlist, labellist= Tracks_in_lat_lon(vortlist, lonlist, latlist, labellist, latbound, lonbound, timebound)

    
    
    """ change the formate - get one cyclone"""
    cyclnratt= labellist[t]
    
    if colortype==1:
        color = next(ax1._get_lines.prop_cycler)['color']
        if alldiffmarker==True: mark= marker[iexp]
    
    if colortype>=2:
        icolorlabel= iexp
#        if iexp >= 4 and colortype==2: icolorlabel+=1 #to skip a yellow color since otherwise the experiments can not be seperated easily
        color= colorlist[icolorlabel]
        if alldiffmarker==True: mark= marker[icolorlabel]
        
    for i, cyclnr in enumerate(remove_dublicate2D(labellist)): # [2]: #cyclnratt:
        tcycl, vortcycl, loncycl, latcycl= Data_One_Track(vortlist, lonlist, latlist, labellist, cyclnr)


        
        if i== 0:
            if alldiffmarker==False: mark= marker[i]
            
#            label= exp_name[7:]
#            if label== 'cold_pseudo2': label= 'cold'
            ax1.plot(d.datetime[tcycl], vortcycl, color= color, linestyle=mark, label= legendlabel[iexp], lw= 2) # +'_'+ str(cyclnr))
        elif i < maxnri:
            if alldiffmarker==False: mark= marker[i]
            ax1.plot(d.datetime[tcycl], vortcycl, color= color, linestyle=mark, lw= 2)


        edgedist =15//AAres #distance from boundary of AA domain in order that the max is not taken at the boundary where disturbances can occure
        windcycl, windmaxlon, windmaxlat= Track_OtherMax(U[:, edgedist: -edgedist, edgedist: -edgedist], winddist, tcycl, latcycl, d.lat[edgedist: -edgedist, edgedist: -edgedist], d.lon[edgedist: -edgedist, edgedist: -edgedist])      #the max should not be on the boundary itself 
        slpcycl, slpmaxlon, slpmaxlat= Track_OtherMax(d.mslp_ft, slpdist, tcycl, latcycl, d.lat, d.lon, local='min')

        edgedist = int(1.5* theta_gausfilterdist/(AAres*2.5) ) #distance from boundary of AA domain in order that the max is not taken at the boundary where disturbances can occure
        
        tlim= 30 #until when the baroclinicity is plotted
        grad_thetacycl, grad_thetamaxlon, grad_thetamaxlat= Track_OtherMax(abs_grad_baroclin_filter[:, edgedist: -edgedist, edgedist: -edgedist], grad_thetadist, tcycl[tcycl<=tlim], latcycl, d.lat[edgedist: -edgedist, edgedist: -edgedist], d.lon[edgedist: -edgedist, edgedist: -edgedist], local='max')
        
        
        if Heat== 1: #calculate mean heat
    #        Heatdistkm= 400
            Heatdist= 300/(d.res *2.5)
            edgedist= 15//AAres
            
            tcycl_red= tcycl[tcycl < len(d.tim)-1] #since the heat fluxes have one shorter array
            
            #the flux maxes, but is not soo good since the condensational part is very dependent on the filter then
    #        SenHcycl, SenHlon, SenHlat= Track_OtherMax(d.SenH_ft[:, edgedist: -edgedist, edgedist: -edgedist], Heatdist, tcycl_red, latcycl, d.lat[edgedist: -edgedist, edgedist: -edgedist], d.lon[edgedist: -edgedist, edgedist: -edgedist], local='max')
    #        LatHcycl, LatHlon, LatHlat= Track_OtherMax(d.LatH_ft[:, edgedist: -edgedist, edgedist: -edgedist], Heatdist, tcycl_red, latcycl, d.lat[edgedist: -edgedist, edgedist: -edgedist], d.lon[edgedist: -edgedist, edgedist: -edgedist], local='max')
    #        Condenscycl, Condenslon, Condenslat= Track_OtherMax(Condens[:, edgedist: -edgedist, edgedist: -edgedist], Heatdist, tcycl_red, latcycl, d.lat[edgedist: -edgedist, edgedist: -edgedist], d.lon[edgedist: -edgedist, edgedist: -edgedist], local='max')
        

            cpt= [np.where(d.lon== loncycl[t]) for t in tcycl]
        #        cpt= [np.where(d.lon== windmaxlon[t]) for t in tcycl]
            
            SenHMean= MeanAroundPoint(cpt, Heatdist, d.SenH_ft[tcycl_red])
            LatHMean= MeanAroundPoint(cpt, Heatdist, d.LatH_ft[tcycl_red])
            CondensMean= MeanAroundPoint(cpt, Heatdist, Condens[tcycl_red])
        


        if i < maxnri:
            if alldiffmarker==False: mark= marker[i]           
            ax2.plot(d.datetime[tcycl], windcycl, color= color, linestyle=mark, lw= 2)#, label= str(cyclnr))       
            ax3.plot(d.datetime[tcycl], slpcycl, color= color, linestyle=mark, lw= 2)#, label= str(cyclnr))
            if Baro== 1: ax4.plot(d.datetime[tcycl[tcycl<=tlim]], grad_thetacycl, color= color, linestyle=mark, lw= 2)#, label= str(cyclnr))

            if Heat== 1:
                ax5.plot(d.datetime[tcycl_red], SenHMean, color= color, linestyle=marker[i], lw= 2)#, label= str(cyclnr))
                ax5.plot(d.datetime[tcycl_red], LatHMean, color= color, linestyle=marker[i+1], lw= 2)#, label= str(cyclnr))
                ax5.plot(d.datetime[tcycl_red], CondensMean, color= color, linestyle=marker[i+2], lw= 2)#, label= str(cyclnr))
            
            #the flux maxes
#            ax5.plot(d.datetime[tcycl_red], SenHcycl, color= color, linestyle=marker[i])#, label= str(cyclnr))
#            ax5.plot(d.datetime[tcycl_red], LatHcycl, color= color, linestyle=marker[i+1])#, label= str(cyclnr))
#            ax5.plot(d.datetime[tcycl_red], Condenscycl, color= color, linestyle=marker[i+2])#, label= str(cyclnr))



        """plot some of the data in map"""
        ind= np.where(tcycl == t)[0]
        if len(ind)> 0: #only if the cyclone exists for this point in time
            xpt, ypt= map(loncycl[ind], latcycl[ind])
    
            map.plot(xpt,ypt, 'o', color= 'r')        
            plt.text(xpt, ypt+ 0.01*(map.ymax-map.ymin), str(int(np.round(vortcycl[ind]))), color='r', fontsize=15, fontweight= 'bold')
    
            xpt, ypt= map(windmaxlon[ind], windmaxlat[ind])
#            map.plot(xpt,ypt, 'o', color= 'orange')
            plt.text(xpt- 0.04*(map.xmax-map.xmin), ypt- 0.02*(map.ymax-map.ymin), str(int(np.round(windcycl[ind]))), color='orange', fontsize=15, fontweight= 'bold')

            """if plot a circle"""
#            dist= 110* np.sqrt((np.cos(np.deg2rad(d.lat))*(d.lon - lonlist[ind[0]]))**2 + (d.lat - latlist[ind[0]])**2)
#            dist= 110* np.sqrt((np.cos(np.deg2rad(d.lat))*(d.lon - windmaxlon[ind[0]]))**2 + (d.lat - windmaxlat[ind[0]])**2)


#            datanow= U[t]
#            datanow[dist > 300]= np.nan
#            PlotWindVelo(Lon, Lat, datanow, map, Umax= 25, color='YlBu')
#            datanow= d.SenH_ft[t]
#            datanow= Condens[t]
#            datanow[dist> Heatdist*(2.5*d.res)]= np.nan
            
#            PlotColorMap4(Lon, Lat, datanow, map, color='blue', bounds= np.arange(0, 601, 100), label=r"Sensible heat flux [W/m$^2$]")


            if mapvar == 'Baroclin_lev':
                xpt, ypt= map(grad_thetamaxlon[ind], grad_thetamaxlat[ind])
#                map.plot(xpt,ypt, 'o', color= 'green')
                plt.text(xpt- 0.04*(map.xmax-map.xmin), ypt- 0.01*(map.ymax-map.ymin), str(np.round(np.round(grad_thetacycl[ind])[0], 1)), color='green', fontsize=15, fontweight= 'bold')


        plt.tight_layout()
        if save==True:
            savename= (savedir+'Arome_'+legendlabel[iexp]+'_'+str(t).zfill(2)+'_'+ mapvar+ # str(d.chosen_plevel) +
                        title_extra)
            if Plotbox: savename +='_box'
            plt.savefig(savename)
            print(savename)
    

ax1.set_ylabel('Vorticity [10$^{-5}$ s$^{-1}$]', fontsize= 13)
ax2.set_ylabel('Wind speed [ms$^{-1}$]', fontsize= 13)
ax3.set_ylabel('Sea level pressure [hPa]', fontsize= 13)
if Baro== 1: ax4.set_ylabel(Baro_label, fontsize= 13 ) #r"$\nabla \theta_{e,"+str(d.chosen_plevel)+"}$ [K/100km]", fontsize= 12)
if Heat== 1: ax5.set_ylabel("Heat fluxes [W m$^{-2}$]", fontsize= 13)


for ax in axlist:
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d.%m'))
    ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=[6,12,18]))
    ax.tick_params(labelsize=13)
    ax.set_xbound([d.datetime[0], d.datetime[-8]])
    
# to change the yaxis for the SLP
ax3.set_yticks(np.arange(980, 1001, 5))
#from matplotlib.ticker import FormatStrFormatter
#ax3.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))


ax1.legend(fontsize= 12, ncol= 4)
if Heat== 1: ax5.legend(['Sensible flux', 'Latent flux', 'Latent release'], fontsize= 12, ncol= 1) #be careful with the order

plt.tight_layout()

if save==True:
    savename=''
    for legendlabeli in legendlabel:
        savename += '_'+ legendlabeli #i[11:]
    if Baro==1: savename += '_Baro_'+Baroclin
    
    if Heat==1: savename += '_Heat'
    
    mainfig.savefig(savedir+'Arome_intensity'+savename, bbox_inches='tight' )



end = time.time()
print('execution time: ', end - start)


"""find the most important cyclone, in terms of vorticity, that also exists at time tn"""
#for l in labellist[t]:
#    tcycl, datacycl, loncycl, latcycl= Data_One_Track(maxvortlist, lonlist, latlist, labellist, l)
#    print(l, sum(datacycl))