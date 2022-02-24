#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 10:01:45 2019

@author: pst019
"""


import time
start = time.time()


import os
user = os.getcwd().split('/')[2]

if user=='pst019':
    homedir= '/home/'+user+'/home/'
    Mediadir= '/media/'+user+'/PatsOrange/'

else:
    Mediadir= '/run/media/pst019/1692A00D929FEF8B/'
'/home/'+user+'/home/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from eofs.standard import Eof

from f_useful import *
from f_meteo import *
from f_STARS import *
from f_carto import *
import cartopy.crs as ccrs
from scipy.interpolate import griddata
import netCDF4 as nc
from datetime import datetime
import great_circle_calculator.great_circle_calculator as gcc # pip install great-circle-calculator
import csv
import xarray as xr


savedir= homedir+ 'Polar_Low/STARS-ana/Figs/EOF2/'
save= True
#save= False

fignr= 7

imp_lists=True
#imp_lists=False

##import from the raw ERA-5 data
##imp_raw_data=True
#imp_raw_data=False
#
##save the PL-CG in .nc file
##save_local_field= True
#save_local_field= False

#import from the the PL-CG .nc file
imp_rot_nc=True
#imp_rot_nc=False

##only if imp_rot_nc is True
#impallObs=True #'import from all time step list'
##impallObs=False #import from datasets of one single timestep per PL


#write the PC to a csv file
write=True #write the PC to a csv file
#write=False


"""specify system ID, Obsnr"""
#Obs= 1
#Obs='mature'
Obs='allObs'


#PLCG_type='track_smth'
smooth_param= '1E-3'

lifetime= 6 #can be set to None if no systems should be excluded by the lifetime
#track_type='Rojo'
track_type='Stoll'




radius= 500 #radius of the local grid km
grid_dist= 25 #in km
ncells= radius//grid_dist
dist_grid_axis= np.arange(-ncells, ncells+1)*grid_dist

steering_rad= 200



plevel=None
vmin,vmax=None,None


"""specify field"""
var_type=''
#var_type='ano'

filetype='tracking'
plevel= 850
#var='vo'
#var='z'
#var='U'
#var='q'
var='t'




sym= False
#cmap='jet'
#cmap= 'viridis'
#cmap= 'Greys'
#cmap= 'Greys_r'
#cmap= 'Blues'
#cmap= 'Reds'

#sym= True
cmap= 'RdBu_r'


if imp_lists:
    if track_type=='Rojo':
        S = import_STARS(Mediadir, "PL/STARS/Rojo-etal_2019.csv")
        S_ind= STARS_individual_systems(S)
    
    elif track_type== 'Stoll':
        test= 'version4'
        dist= 150
        
        Stoll_STARS_file="ERA5_STARS/STARS-PMC-merge/"+test+"_dist_"+str(dist)+"_handfix.csv"
        
        Stoll = pd.read_csv(Mediadir+ Stoll_STARS_file)
        Stoll['time']= pd.DatetimeIndex(Stoll.time)
        
        Stoll= Stoll_excl(Stoll)#, lifetime_threshold= lifetime)
        Stoll= Stoll_Obs_nr(Stoll) #get the Obs nr
    
        Stoll= Stoll.rename(columns={"Stoll nr": "ID"})
        
        
        """Stoll_individual_systems"""
        Stoll_ind= Stoll_individual_systems(Stoll, ID_name= 'ID')    
    
        S= Stoll
        S_ind= Stoll_ind

    if Obs =='mature':
        print('have to fix this:')
        #    S_ind= calc_Obs_mature(S, S_ind)


"""get the rotated nc file"""
if imp_rot_nc:
#    if impallObs:
    print('import from all time step list')
    ds= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/" +var + '_'+ str(plevel) + '_allObs_track-smth-'+smooth_param+'.nc')

        
    ds['ID']= ('time', [str(int(ds.PLnr[n]))[:-3]+'_'+str(int(ds.PLnr[n])%1000) for n in range(len(ds.PLnr))] )
    
    if type(Obs)== int: ds_now= ds.where(ds.Obs== Obs, drop=True)
    elif Obs =='mature':
        mature_index= [np.where(np.logical_and(ds.ID== S_ind.ID[n], ds.Obs== S_ind.Obs_mature[n]))[0][0] for n in range(len(S_ind.ID)) ]
        ds_now= ds.isel(time= mature_index)
    elif Obs =='allObs':
        ds_now= ds

          
        
    ds_now= ds_now.dropna(dim='time') #to exclude nan values for touching boundary
#    print('excluded PLs due to nan values:', set(S_ind.ID.values) - set(ds_now.ID.values) )
    
    if lifetime: #exclusion due to short lifetime
        ds_now['duration'] = ('time', np.zeros(len(ds_now['PLnr'])) )
        duration_ds=ds_now.groupby(ds_now['PLnr']).last()['Obs']
        for PLnow in duration_ds.PLnr.values:
            ds_now['duration'][ds_now['PLnr']== PLnow] =duration_ds.sel(PLnr= PLnow).values
        
        ds_now= ds_now.where(ds_now.duration >= lifetime, drop=True)
    
    
    var_interp= ds_now[var].values
    PLID= ds_now['ID'].values
        
  


### have to do something if only single time step is of interest
#    else: #import from datasets of one single timestep per PL
#        if PLCG_type == 'stearing_flow':
#            if type(Obs)== int: ds= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields/"+var+ '_'+ str(plevel)+ '_Obsnr'+str(Obs)+'.nc')  
#            else: ds= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields/" +var + '_'+ str(plevel) +"_" + Obs +'.nc')
#            
#            var_interp= ds[var].values   
#        
#            PLID= [str(int(ds.PLnr[n]))[:-3]+'_'+str(int(ds.PLnr[n])%1000) for n in range(len(ds.PLnr))]
#        elif PLCG_type == 'track_smth': print('not included')
#            




var_ano= var_interp - np.mean(var_interp, axis= (1,2))[:, np.newaxis, np.newaxis]    
if var_type=='ano': var_now= var_ano
else: var_now= var_interp

       




print('shape of the utilized variable:', np.shape(var_now))



"""plot the mean"""
plt.figure(fignr, figsize= (5.5,8) )
fignr+=1
plt.clf()

ax1= plt.subplot(211)




 
varmean= np.nanmean(var_now, axis= 0)    
    
    
if sym== False:
    cf= plt.contourf(dist_grid_axis, dist_grid_axis, varmean,  cmap= cmap)
    cs= plt.contour(dist_grid_axis, dist_grid_axis, varmean,  colors='k', linewidths= 1)
    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.1f')


elif sym== True:
    vextr= np.max([np.max(varmean), -np.min(varmean)])
    cf= plt.pcolormesh(dist_grid_axis, dist_grid_axis, varmean,  cmap= cmap, vmin= -vextr, vmax= vextr)


cb= plt.colorbar(cf, label= 'Mean ' + ds[var].long_name + ' at '+str(plevel)+ 'hPa [' + ds[var].units + ']' )

varlabel= var
if plevel != None: varlabel = var+'_'+str(plevel)
#cb.set_label(varlabel + ' ['+ ds.units + ']', size=14)  
#plt.title('Mean')



"""plot the standard deviation"""
ax2= plt.subplot(212, sharex= ax1)

varstd= np.nanstd(var_now, axis= 0)    

cf= plt.contourf(dist_grid_axis, dist_grid_axis, varstd, cmap= 'Reds')
cb= plt.colorbar(cf, label= ds[var].long_name + ' standard deviation [' + ds[var].units + ']' )

cs= plt.contour(dist_grid_axis, dist_grid_axis, varstd, colors='k', linewidths= 1)
plt.clabel(cs, fontsize=10, inline=1, fmt='%1.1f')

#cb.set_label(varlabel + ' ['+ ds.units + ']', size=14)  

#plt.title('Standard deviation')


plt.tight_layout()

if save:
    if type(Obs)== int:save_name=savedir+ 'EOF_track-smth'+smooth_param+'_'+ var +var_type+ '_'+ str(plevel) + '_Obsnr'+ str(Obs)+  '_mean'
    else: save_name=savedir+ 'EOF_track-smth'+smooth_param+'_'+ var +var_type+ '_'+ str(plevel) + '_'+Obs+ '_mean'

    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')


"""calculate the EOFs"""
from eofs.standard import Eof

#dsano= ds - dsmean

fig = plt.figure(fignr, figsize= (8,8) )
fignr+=1
plt.clf()



solver = Eof(var_now)

EOFnr= 4


var_frac= np.round(solver.varianceFraction(neigs=EOFnr)*100, 1)
print('Variance fraction [%]: ', var_frac)

eigenvalue_error= solver.northTest(neigs=EOFnr, vfscaled=True)
print('Eigenvalue error [%]: ', np.round(eigenvalue_error*100, 1))


Sign_adapt= np.array([1, 1, -1, 1])


eofs = solver.eofsAsCovariance(neofs=EOFnr, pcscaling= 1)#Empirical orthogonal functions (EOFs) expressed as the covariance between the PCs and the time series of the Eof input dataset at each grid point.  pcscaling= 1 - scaling to unit variance
#eofs = solver.eofsAsCorrelation(neofs=EOFnr) #Empirical orthogonal functions (EOFs) expressed as the correlation between PCs and the time series of the Eof input dataset at each grid point.
#eofs = solver.eofs(neofs=EOFnr)#, eofscaling= 2) #this is the same as "eofsAsCovariance with pcscaling=1


eofs *= Sign_adapt[:, np.newaxis, np.newaxis]

"""find which PC is most important for each PL"""
PC= solver.pcs(npcs= EOFnr, pcscaling= 1)

PC *= Sign_adapt[np.newaxis, :]


#PCabs= np.abs(PC)
#MostImpPCofPL= np.where( PCabs== np.max(PCabs, axis= 1))[1] #finds the most important PC for each PL
#PCSign= np.sign(PC.values[(PCabs== np.max(PCabs, axis= 1)).values]) #the sign of the most important PC

for iEOFnr in range(EOFnr):

    plt.subplot(2,2,iEOFnr+1)
    
   
    vextr= np.max([np.max(eofs), -np.min(eofs)])
    cf= plt.contourf(dist_grid_axis, dist_grid_axis, eofs[iEOFnr], cmap= cmap, vmin= -vextr, vmax= vextr)
    if iEOFnr== 1: cf_legend= cf
#    cb= fig.colorbar(cf, ax= ax, shrink=0.7)
    
    cs= plt.contour(dist_grid_axis, dist_grid_axis, eofs[iEOFnr], colors='k', linewidths= 1)
#    roundlevel= int(np.log10(vextr) + 1)
    plt.clabel(cs, fontsize=10, inline=1, fmt='%1.1f')
#    varlabel= var
#    if plevel != None: varlabel = var+'_'+str(plevel)
#    #cb.set_label(varlabel + ' ['+ ds.units + ']', size=14)  
    plt.title('EOF '+str(iEOFnr+1)+' ('+str(int(np.round(var_frac[iEOFnr])))+'%)')

#    #plot the points where the EOF is in positive phase in yellow
#    PLsThisPC= S[S['Obs'] == Obs ][np.logical_and(MostImpPCofPL == iEOFnr, PCSign >= 1)] #get all PLs that are connected to this PC
#    ax.scatter(PLsThisPC.lon, PLsThisPC.lat, transform=ccrs.PlateCarree(), s= 1, color='y')
#
#    #plot the tracks of the most dominant PC
#    S_ID_list= list(S[S['Obs'] == Obs ][np.logical_and(MostImpPCofPL == iEOFnr, PCSign >= 1)]['ID'])
#    S_now= S[ [ID in S_ID_list for ID in S['ID']] ]
#    S_tracks= prepare_tracks(S_now, system_id_label="ID")
#    for track in S_tracks:
#        track.plot_track(ax=ax, color='y')#, label='Stoll excluded')    
#
#    
#    PLsThisPC= S[S['Obs'] == Obs ][np.logical_and(MostImpPCofPL == iEOFnr, PCSign <= -1)] #get all PLs that are connected to this PC
#    ax.scatter(PLsThisPC.lon, PLsThisPC.lat, transform=ccrs.PlateCarree(), s= 1, color= 'g')
#
#    S_ID_list= list(S[S['Obs'] == Obs ][np.logical_and(MostImpPCofPL == iEOFnr, PCSign <= -1)]['ID'])
#    S_now= S[ [ID in S_ID_list for ID in S['ID']] ]
#    S_tracks= prepare_tracks(S_now, system_id_label="ID")
#    for track in S_tracks:
#        track.plot_track(ax=ax, color='g')#, label='Stoll excluded')   


plt.tight_layout()

fig.subplots_adjust(bottom=0.11)
cbar_ax = fig.add_axes([0.09, 0.06, 0.6, 0.015])
#cb= fig.colorbar(cf, cax=cbar_ax, orientation="horizontal")
cb= plt.colorbar(cf_legend, cax=cbar_ax, orientation="horizontal",
                 label= ds[var].long_name + ' [' + ds[var].units + ']' )


"""make the arrow"""
arrow_ax = fig.add_axes([0.75, 0.03, 0.2, 0.04])
import matplotlib.patches as mpatches

arrow = mpatches.FancyArrowPatch((0.2, 0.7), (.7, 0.7), mutation_scale=20)
arrow_ax.add_patch(arrow)
arrow_ax.text(0,0 , 'Propagation direction', fontweight= 'bold')#, fontsize=13)
arrow_ax.set_frame_on(False)
arrow_ax.get_xaxis().set_visible(False)
arrow_ax.get_yaxis().set_visible(False)

if save:
    if type(Obs)== int:save_name=savedir+ 'EOF_track-smth'+smooth_param+'_'+ var +var_type+ '_'+ str(plevel) + '_Obsnr'+ str(Obs)+  '_EOF'
    else: save_name=savedir+ 'EOF_track-smth'+smooth_param+'_'+ var +var_type+ '_'+ str(plevel) + '_'+Obs+ '_EOF'

    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')





if write: 
    print('write the PCs into csv file')

    if type(Obs)== int: csv_name= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/PCs/"+var+ var_type+ str(plevel) +'_Obsnr'+ str(Obs)+ '_track-smth'+smooth_param+ "_PCs.csv"
    else: csv_name= Mediadir + "ERA5_STARS/PL_centred_fields_smooth-tracks/PCs/"+var+ var_type + str(plevel)+ '_'+Obs+ '_track-smth'+smooth_param+ "_PCs.csv"
    
    print(csv_name)
    
    col_names=[var+ var_type +'_PC'+str(iEOFnr+1) for iEOFnr in range(EOFnr) ]
    df= pd.DataFrame(PC, columns= col_names)
    df['ID']= ds_now.ID
    df['Obs']= ds_now.Obs
    df= df.set_index('ID')
    
    exist= os.path.isfile(csv_name) 
    df.to_csv(csv_name, sep=',' )
        

    

