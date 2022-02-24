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

marker = ['-', '--', ':', '-.', '>', 's', '8', 'p']
maxnri = 1 #maximum number of cyclones for one experiment in the intensity graph

colorlist= [plt.cm.Spectral(h) for h in range(0, 256, 256//8)] #for specified colors (colortype= 2)
colortype= 1 #use the python default colors - might be specified differently below

"""global variables"""
fignr= 5

mainfig= plt.figure(fignr, figsize= (7,8))
fignr+=1
plt.clf()
                  
ax1= plt.subplot(311)    
ax2= plt.subplot(312)
ax3= plt.subplot(313)                   

"""b) times"""
year, month = 2008, 3
day, hour= 4, 0 #22                   


"""Lambert coordinates"""
#maptype=='AA'
maptype='AA_half'

fileday, filehour= 3, 0        
t= (day- fileday)*24 + (hour- filehour)  # -1 
 
    
"""name of the experiment"""
#exp_name, label_name= '080303_warmctr', 'CTR'
#exp_name='080303_cold_sice'
#exp_name= '080303_warmsens_noQH'
#exp_name, label_name='080303_warmsens_2FLX', '2FLX'
#exp_name='080303_warmsens_nocondens'
#exp_name='870226_cold_ERA_can'
#exp_name, label_name= '080303_cold_ERA', 'COLD-ERA'
#exp_name='080303_cold_pseudo2'
#exp_name='080304_cold_pseudo'

#exp_name, label_name='080303_coldsens_noFLX_AREA', 'noFLX-Area'
#exp_name, label_name= '080303_cold', 'COLD'
#exp_name= '080303_warmsens_notopo'

exp_name, label_name= 'DA_080303_CTR', 'CTR'
#exp_name, label_name= 'DA_080303_noFLX', 'noFLX'
#exp_name, label_name= 'DA_080303_noFLXarea', 'noFLXarea'
#exp_name, label_name= 'DA_080303_p2SST', 'p2SST'

#exp_name, label_name= 'DA_080301_cycling', 'cycling'
#exp_name, label_name= 'DA_08030212', 'early_start'

"""c) plot characteristics"""

var= 'PressWindVort'
#var= 'Vort_lev'


""" presure level (if the variable ends with '_var' """
pn = 4
#1-950, 4- 850hPa, 6- 700hPa, 8 - 500hPa

save=False
savedir= homedir+'Polar_Low/AromeArctic-vs-WRF/2008_03_03/IntensityFieldsTracks4/'

title_extra=''

"""end global variables"""
AAres=2 #2- every second datapoint is taken


threshold= 8 #10
propdistance= 180E3 #in km   #70
distance= int(200E3/(AAres * 2.5E3)) #100   #make distance reasonable larger than propdistance to avoid "doublematches" twice is on the safe side


latbound= [60, 71]
lonbound= [-5, 10]
timebound= [20, 30]

winddist= 400/(AAres *2.5) #150
slpdist= 100/(AAres *2.5)
gausfilterdist= 100 #60
gaustruncate= 2


Plotbox=False
boxlon, boxlat= [0, 10, 10, 0, 0], [67.9, 67.9, 74, 74, 68] #to display the area without fluxes

explist=[exp_name]
legendlabel= [label_name]


#explist= ['DA_080303_CTR', 'DA_080303_2FLX']
#legendlabel= ['CTR', '2FLX']

#explist= ['DA_080303_CTR', 'DA_080303_noFLXarea', 'DA_080303_noFLX', 'DA_080303_2FLX']
#legendlabel= ['CTR', 'noFLXarea', 'noFLX', '2FLX']
#
#explist= ['DA_080303_CTR', 'DA_080303_noFLXarea']
#legendlabel= ['CTR', 'noFLXarea', 'noFLX', '2FLX']


#explist= ['DA_080303_p6SST', 'DA_080303_p4SST', 'DA_080303_p2SST', 'DA_080303_CTR', 'DA_080303_m2SST', 'DA_080303_m4SST', 'DA_080303_m6SST']
#legendlabel= ['p6SST', 'p4SST', 'p2SST', 'CTR', 'm2SST', 'm4SST', 'm6SST']
#colortype= 2

#explist= ['DA_080303_CTR', 'DA_080303_noQH', 'DA_080303_noCondens', 'DA_080303_noTH', 'DA_080303_noFLX']
#legendlabel= ['CTR', 'noQH', 'noCondens', 'noTH', 'noFLX']
#
#explist= ['DA_080303_CTR', 'DA_080303_noHARATU', 'DA_080303_noOCND2']
#legendlabel= ['CTR', 'noHaratu', 'noOCND2']

#explist= ['DA_080303_CTR', 'DA_080303_noQH', 'DA_080303_noCondens']
#legendlabel= ['CTR', 'noQH', 'noCondens']

#explist= ['DA_080303_m6SST', 'DA_080303_noFLX']
#legendlabel= ['m6SST', 'noFLX']


for iexp,exp_name in enumerate(explist): 
    
    """import data"""
    AAfilename= Mediadir+'PL/AA/ec/'+exp_name+'_'+str(year)+str(month).zfill(2)+str(fileday).zfill(2)+str(filehour).zfill(2)+'_fp_extract.nc'
                                                     
    d= data(filename= AAfilename, res= AAres)     
        
    d.imp_uv_ft(pn= pn)

    d.imp_surf_ft()

    U= np.sqrt(d.u10m_ft**2+d.v10m_ft**2)

         
    """calculate vorticity"""
    dx= 2500* AAres
    
    vort_full= (np.gradient(d.v_ft, dx, axis= -1)- np.gradient(d.u_ft, dx, axis= -2))*1E5     
    
    vortfilter_full= filters.gaussian_filter1d(vort_full, axis= -1, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= gaustruncate)
    vortfilter_full= filters.gaussian_filter1d(vortfilter_full, axis= -2, sigma= gausfilterdist/(AAres*2.5), mode='constant', truncate= gaustruncate)
        

    """Plot the data for one experiment"""
    if maptype=='AA': fig= plt.figure(fignr, figsize= (6, 6))  
    elif maptype=='AA_half': fig=plt.figure(fignr, figsize= (6, 4.8))
    fignr+=1
    plt.clf()    

    if maptype=='AA': map = AA_map()
    elif maptype=='AA_half': map = AA_map_half()

    Lon,Lat = map(d.lon,d.lat)

    PlotContours(Lon, Lat, d.mslp_ft[t], map, leveldist= 2, numbers=False)   
    
    U= np.sqrt(d.u10m_ft**2+d.v10m_ft**2)
    PlotWindVelo(Lon, Lat, U[t], map, Umax= 25, color='YlBu')

    if save== False: plt.title('Arome '+ legendlabel[iexp] + ' '+ str(d.datetime[t])[:-3]+ title_extra)
    
    #find local pressure min and plot them in the map 
#    maxdata, maxlon, maxlat= PlotLocalMax(vortfilter_full[t], threshold=threshold, distance=distance, map= map, lon=d.lon, lat=d.lat, typ='max', color='purple')   
#    PlotLocalMax(d.mslp_ft[t], threshold=1010, distance=100/AAres, map= map, lon=d.lon, lat=d.lat,
#                 data2=U[t], threshold2=16, distance2=80/AAres)         
#    PlotLocalMax(U[t], threshold=20, distance=100/AAres, map= map, lon=d.lon, lat=d.lat, typ='max', color='orange', dot=True, yoffset= -200)    


    if Plotbox== True:
        x,y= map(boxlon, boxlat) #longitude and latitude coordinates of the corner points
        map.plot(x, y, color='b', linewidth= 2)
        title_extra +='_box'           
    
    
    """ a tracking algorithm"""
    vortlist, lonlist, latlist, labellist = EasyTrack(vortfilter_full, d.lon, d.lat, distance, threshold, propdistance)
    vortlist, lonlist, latlist, labellist= Tracks_in_lat_lon(vortlist, lonlist, latlist, labellist, latbound, lonbound, timebound)

    
    
    """ change the formate - get one cyclone"""
    cyclnratt= labellist[t]
    
    if colortype==1: color = next(ax1._get_lines.prop_cycler)['color']
    if colortype==2:
        icolorlabel= iexp
        if iexp >= 4: icolorlabel+=1 #to skip a yellow color since otherwise the experiments can not be seperated easily
        color= colorlist[icolorlabel]
    
    for i, cyclnr in enumerate(remove_dublicate2D(labellist)): # [2]: #cyclnratt:
        tcycl, vortcycl, loncycl, latcycl= Data_One_Track(vortlist, lonlist, latlist, labellist, cyclnr)

        
        if i== 0:
#            label= exp_name[7:]
#            if label== 'cold_pseudo2': label= 'cold'
            ax1.plot(d.datetime[tcycl], vortcycl, color= color, linestyle=marker[i], label= legendlabel[iexp]) # +'_'+ str(cyclnr))
        elif i < maxnri:
            ax1.plot(d.datetime[tcycl], vortcycl, color= color, linestyle=marker[i])


        Udist =10//AAres #distance from boundary of AA domain in order that the max is not taken at the boundary where disturbances can occure
        windcycl, windmaxlon, windmaxlat= Track_OtherMax(U[:, Udist: -Udist, Udist: -Udist], winddist, tcycl, latcycl, d.lat[Udist: -Udist, Udist: -Udist], d.lon[Udist: -Udist, Udist: -Udist])      #the max should not be on the boundary itself 
        slpcycl, slpmaxlon, slpmaxlat= Track_OtherMax(d.mslp_ft, slpdist, tcycl, latcycl, d.lat, d.lon, local='min')

        if i < maxnri:
            ax2.plot(d.datetime[tcycl], windcycl, color= color, linestyle=marker[i])#, label= str(cyclnr))       
            ax3.plot(d.datetime[tcycl], slpcycl, color= color, linestyle=marker[i])#, label= str(cyclnr))



        """plot some of the data in map"""
        ind= np.where(tcycl == t)[0]
        if len(ind)> 0: #only if the cyclone exists for this point in time
            xpt, ypt= map(loncycl[ind], latcycl[ind])
    
            map.plot(xpt,ypt, 'o', color= 'r')        
            plt.text(xpt, ypt+ 0.01*(map.ymax-map.ymin), str(int(np.round(vortcycl[ind]))), color='r', fontsize=15, fontweight= 'bold')
    
            xpt, ypt= map(windmaxlon[ind], windmaxlat[ind])
    #        map.plot(xpt,ypt, 'o', color= 'orange')        
            
            plt.text(xpt- 0.04*(map.xmax-map.xmin), ypt- 0.01*(map.ymax-map.ymin), str(int(np.round(windcycl[ind]))), color='orange', fontsize=15, fontweight= 'bold')

        plt.tight_layout()
        if save==True:
            if '_lev' in var:
                var += str(d.chosen_plevel)
            plt.savefig(savedir+'Arome_'+legendlabel[iexp]+'_'+str(t).zfill(2)+'_'+ var +title_extra)
            print(savedir+'Arome_'+legendlabel[iexp]+'_'+str(t).zfill(2)+'_'+ var +title_extra)
    
ax1.legend(fontsize= 12, ncol= 2)

ax1.set_ylabel('Vorticity [10$^{-5}$ s$^{-1}$]', fontsize= 12)
ax2.set_ylabel('Wind speed [ms$^{-1}$]', fontsize= 12)
ax3.set_ylabel('Sea level pressure [hPa]', fontsize= 12)

for ax in [ax1, ax2, ax3]:
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d.%m'))
    ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=[6,12,18]))
    ax.tick_params(labelsize=11)
    ax.set_xbound([d.datetime[0], d.datetime[-6]])
    
# to change the yaxis for the SLP
ax3.set_yticks(np.arange(975, 1001, 5))
#from matplotlib.ticker import FormatStrFormatter
#ax3.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))


plt.tight_layout()

if save==True:
    savename=''
    for legendlabeli in legendlabel:
        savename += '+'+ legendlabeli #i[11:]
        
    mainfig.savefig(savedir+'Arome_intensity'+savename, bbox_inches='tight' )



end = time.time()
print('execution time: ', end - start)


"""find the most important cyclone, in terms of vorticity, that also exists at time tn"""
#for l in labellist[t]:
#    tcycl, datacycl, loncycl, latcycl= Data_One_Track(maxvortlist, lonlist, latlist, labellist, l)
#    print(l, sum(datacycl))