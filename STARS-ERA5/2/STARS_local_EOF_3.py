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
    Mediadir= '/media/'+user+'/1692A00D929FEF8B/'

else:
    Mediadir= '/run/media/pst019/1692A00D929FEF8B/'
'/home/'+user+'/home/'

import sys  #to import the functions from a different directory
sys.path.insert(0, homedir+ 'Polar_Low/polar_low_code/Functions')

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
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


savedir= homedir+ 'Polar_Low/STARS-ana/Figs/EOF/'
#save= True
save= False

fignr= 7

imp_lists=True
#imp_lists=False

#import from the raw ERA-5 data
imp_raw_data=True
#imp_raw_data=False

#save the PL-CG in .nc file
#save_local_field= True
save_local_field= False

#import from the the PL-CG .nc file
#imp_rot_nc=True
imp_rot_nc=False

#only if imp_rot_nc is True
impallObs=True #'import from all time step list'
#impallObs=False #import from datasets of one single timestep per PL


#write the PC to a csv file
write=True #write the PC to a csv file
#write=False


"""specify system ID, Obsnr"""
Obs= 1
#Obs='mature'
#Obs='allObs'


prop_speed_excl=3 #to exclude systems below this propagation speed [m/s]

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
        
        Stoll= Stoll_excl(Stoll, lifetime_threshold= 6)
        Stoll= Stoll_Obs_nr(Stoll) #get the Obs nr
    
        Stoll= Stoll.rename(columns={"Stoll nr": "ID"})
        
        
        """Stoll_individual_systems"""
        Stoll_ind= Stoll_individual_systems(Stoll, ID_name= 'ID')    
    
        S= Stoll
        S_ind= Stoll_ind


    S_ind= calc_Obs_mature(S, S_ind)


"""get the rotated nc file"""
if imp_rot_nc:
    if impallObs:
        print('import from all time step list')
        ds= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields/" +var + '_'+ str(plevel) + '_allObs.nc')
        ds['ID']= ('time', [str(int(ds.PLnr[n]))[:-3]+'_'+str(int(ds.PLnr[n])%1000) for n in range(len(ds.PLnr))] )
        
        if type(Obs)== int: ds_now= ds.where(ds.Obs== Obs, drop=True)
        elif Obs =='mature':
            mature_index= [np.where(np.logical_and(ds.ID== S_ind.ID[n], ds.Obs== S_ind.Obs_mature[n]))[0][0] for n in range(len(S_ind.ID)) ]
            ds_now= ds.isel(time= mature_index)
        elif Obs =='allObs':
            ds_now= ds
            
        ds_now= ds_now.where(ds_now.steer_vel > prop_speed_excl, drop=True)            
            
        ds_now= ds_now.dropna(dim='time') #to exclude nan values
        print('excluded PLs due to nan values:', set(S_ind.ID.values) - set(ds_now.ID.values) )
        var_interp= ds_now[var].values
        PLID= ds_now['ID'].values
        
  


    else: #import from datasets of one single timestep per PL
        if type(Obs)== int: ds= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields/"+var+ '_'+ str(plevel)+ '_Obsnr'+str(Obs)+'.nc')  
        else: ds= xr.open_dataset(Mediadir + "ERA5_STARS/PL_centred_fields/" +var + '_'+ str(plevel) +"_" + Obs +'.nc')
        
        var_interp= ds[var].values   
    
        PLID= [str(int(ds.PLnr[n]))[:-3]+'_'+str(int(ds.PLnr[n])%1000) for n in range(len(ds.PLnr))]

        

"""get and process the raw ERA-5 data"""
if imp_raw_data:
    
    for i, ID in enumerate(remove_dublicate(S['ID'])):
        print(i, ID)

        if type(Obs)== int:    
            S_now= S[np.logical_and(S['ID'] == ID, S['Obs'] == Obs ) ]
        if Obs== 'mature':
            S_now= S[np.logical_and(S['ID'] == ID, S['Obs'] == S_ind.loc[ID]['Obs_mature']) ]
        
        datetime_now= S_now.time.dt.round("H")      
        
        filetime= str(datetime_now.dt.year.values[0])+'_'+str(datetime_now.dt.month.values[0]).zfill(2)+'_'+str(datetime_now.dt.day.values[0]).zfill(2)
        hour= datetime_now.dt.hour.values[0]
    
        print(datetime_now)
    
    
    
        """general"""      
        ds0= xr.open_dataset(Mediadir + "ERA5_STARS/data/"+filetype+ "_era5_"+ filetime + '.nc')
        ds0= ds0.isel(time= hour)
        ds0['plev']/= 100

        
        """get the center position, the steering vector"""
        center_lon, center_lat= S_now.lon.values[0], S_now.lat.values[0]
        
        u_steer= value_in_rad(ds0['u'].sel(plev= [1000, 700]).mean(axis= 0), ds0.lat, ds0.lon, center_lat, center_lon, steering_rad)
        v_steer= value_in_rad(ds0['v'].sel(plev= [1000, 700]).mean(axis= 0), ds0.lat, ds0.lon, center_lat, center_lon, steering_rad)

        steer_vel= np.sqrt(u_steer**2 + v_steer**2)
        if steer_vel < 3:
            print('Should be considered if the rotation works with such a steering wind of {} km/h'.format(np.round(steer_vel*3.6)) )


        beering= UV2Direction(u_steer, v_steer) #wind direction

        #the tangental axis in the strearing wind direction
        tanax= [list(gcc.point_given_start_and_bearing((center_lon, center_lat), beering, n*grid_dist*1E3)) for n in np.arange(-ncells, ncells+1)]
        tanax= np.array(tanax)
        
        #the beering along the tangential axis
        point0= gcc.point_given_start_and_bearing((center_lon, center_lat), beering, -(radius+grid_dist)*1E3) #point one before the start of the tangential axis, used for calculation of beering_axis
        beering_axis= [gcc.bearing_at_p2((point0), (tanax[m,0], tanax[m,1])) for m in range(len(tanax))] #the wind direction along the tangential axis
        
        #creation of the PL centred grid
        center_grid= [[list(gcc.point_given_start_and_bearing((tuple(tanax[m])), (beering_axis[m]-90)%360, n*grid_dist*1E3)) for m in range(len(tanax))] for n in np.arange(-ncells, ncells+1)]
        center_grid= np.array(center_grid)



        """ interpolate to local grid"""
        ds0= ds0.sel(plev= plevel)
        if var=='z': ds0[var]/=9.81
        if var=='vo': ds0[var] *=1E5
        if var=='q': ds0[var] *=1E3
        if var=='U': ds0[var]= np.sqrt(ds0['u']**2+ ds0['v']**2)
        
        longrid, latgrid= np.meshgrid(ds0.lon, ds0.lat)
        var_interp0= griddata(tuple([np.ravel(longrid), np.ravel(latgrid)]) , np.ravel(ds0[var]), (center_grid[:,:,0], center_grid[:,:,1]), method='linear')


    
        if i == 0:
            var_interp=np.zeros((0, np.shape(var_interp0)[0],np.shape(var_interp0)[1]) )
            datetime_vec= datetime_now.values #time vector for the output netcdf file
            PLID= [ID]
#            PLnr_list_float= [int(ID.split('_')[0])*1E3 + int(ID.split('_')[1])]
            steer_vel_vec= [steer_vel]
        if len(np.argwhere(np.isnan(var_interp0))) == 0:    #no nan values in data array, can be due to boundaries of domain
            var_interp= np.vstack((var_interp, var_interp0[np.newaxis]))
            if i != 0:
                datetime_vec= np.append(datetime_vec, datetime_now.values)
                PLID += [ID]
#                PLnr_list_float += [int(ID.split('_')[0])*1E3 + int(ID.split('_')[1])]
                steer_vel_vec += [steer_vel]
        else: print('array contains nan values')




var_ano= var_interp - np.mean(var_interp, axis= (1,2))[:, np.newaxis, np.newaxis]    
if var_type=='ano': var_now= var_ano
else: var_now= var_interp

       
if save_local_field:
    ncdstime= nc.Dataset(Mediadir + "ERA5_STARS/data/"+filetype+ "_era5_"+ filetime + '.nc')
    calendar= ncdstime['time'].calendar
    timeunits= ncdstime['time'].units
#    cdftime = nc.netcdftime.utime(timeunits,calendar=calendar)
    ncdstime.close()
    
    
    if type(Obs)== int:
        nc_out_name=Mediadir + "ERA5_STARS/PL_centred_fields/"+var + '_'+ str(plevel) +'_Obsnr'+ str(Obs) +'.nc'
    else:
        nc_out_name= Mediadir + "ERA5_STARS/PL_centred_fields/"+var + '_'+ str(plevel)+"_" + Obs +'.nc'
    

    if os.path.isfile(nc_out_name): print('nc file exists')
#    else:
           
    ncds = nc.Dataset(nc_out_name,'w',format='NETCDF3_CLASSIC')

        
#    ncds = nc.Dataset(Mediadir + "ERA5_STARS/PL_centred_fields/"
#                      +var + '_'+ str(plevel) +'_Obsnr'+ str(Obs) +'.nc','w')
#                lon, lat = m(xeag,yeag,inverse=True)
#        ncds = nc.Dataset(self.bn+'_eag_sca.nc','w') 
    ncds.createDimension('time',var_interp.shape[0])        
    ncds.createDimension('x',var_interp.shape[1])
    ncds.createDimension('y',var_interp.shape[2])
    #ncds.createDimension('xy',field.shape[0]*field.shape[1])
#        ncds.createDimension('time',None)
    x = ncds.createVariable('x','f',('x',))
    x.setncattr('name','Distance in propagation direction')
    x.setncattr('units','km')
    x.setncattr('axis','X')
    x[:] = dist_grid_axis

    y = ncds.createVariable('y','f',('y',))
    y.setncattr('name','Distance perpendicular to propagation direction')
    y.setncattr('units','km')
    y.setncattr('axis','Y')
    y[:] = dist_grid_axis
    
    """add the datetime_vec also add the PL number"""
    nctime = ncds.createVariable('time','f',('time',))
    nctime.setncattr('name','time')
    nctime.setncattr('calendar', calendar)
    nctime.setncattr('axis','T')
    nctime.setncattr('units', timeunits)
    
    time_vec= pd.to_datetime(datetime_vec).values.astype('datetime64[h]')
    time_vec= time_vec.astype(datetime)
    time_vec= nc.date2num(time_vec, timeunits, calendar=calendar) #check why it does not work
    nctime[:]= time_vec

    steer_vel = ncds.createVariable('steer_vel', 'f' ,('time',)) #some PLs have to be excluded.
    steer_vel.setncattr('name','Velocity of the steering wind [m/s]')
    steer_vel[:]= steer_vel_vec


    PLnr_list_float = [int(PLID[n].split('_')[0])*1E3 + int(PLID[n].split('_')[1]) for n in range(len(PLID))]
    ncPLnr = ncds.createVariable('PLnr', 'f' ,('time',)) #some PLs have to be excluded.
    ncPLnr.setncattr('name','Polar low number from Stoll list')
    ncPLnr[:]= PLnr_list_float
#    ncPLnr[:]= nc.stringtochar(np.array(PLnr_list, 'S'))

    vari = ncds.createVariable(var,'f',('time','x','y',)) #,fill_value=self.FillValue)
#    vari.setncattr('coordinates','time x y')
    if var != 'U':
        vari.setncattr('standard_name', ds0[var].standard_name)
        vari.setncattr('long_name', ds0[var].long_name)
        vari.setncattr('units', ds0[var].units)
    vari[:,:,:] = var_interp
    
    varano = ncds.createVariable(var +'_ano', 'f',('time','x','y',)) #,fill_value=self.FillValue)
#    varano.setncattr('coordinates','time x y')
    if var != 'U':
        varano.setncattr('standard_name', ds0[var].standard_name + ' anomaly')
        varano.setncattr('long_name', ds0[var].long_name+ ' anomaly')
        varano.setncattr('units', ds0[var].units)

    varano[:,:,:] = var_ano

    ncds.setncattr('history',"Created by Patrick Stoll on %s." % \
                 (datetime.today().strftime("%Y-%m-%d") ) )
    ncds.close()        





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


cb= plt.colorbar(cf, label= 'Mean ' + ds[var].long_name + 'at '+str(plevel)+ ' [' + ds[var].units + ']' )

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
    if type(Obs)== int: save_name=savedir+ track_type + '_'+ var +var_type+ '_'+ str(plevel)+'_Obsnr'+ str(Obs)+ '_propexcl'+str(prop_speed_excl) + '_mean'
    else: save_name=savedir+ track_type + '_'+ var +var_type+ '_'+ str(plevel) + '_'+Obs+ '_propexcl'+str(prop_speed_excl)+ '_mean'

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


Sign_adapt= np.array([-1, -1, -1, 1])


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
    if type(Obs)== int: save_name=savedir+ track_type + '_'+ var +var_type+ '_'+ str(plevel)+'_Obsnr'+ str(Obs)+ '_propexcl'+str(prop_speed_excl)+ '_EOF'
    else: save_name=savedir+ track_type + '_'+ var +var_type+ '_'+ str(plevel) + '_'+Obs+ '_propexcl'+str(prop_speed_excl)+ '_EOF'

    print(save_name)
    plt.savefig(save_name , bbox_inches='tight')





if write: 
    print('write the PCs into csv file')

    if type(Obs)== int: csv_name= Mediadir + "ERA5_STARS/PL_centred_fields/PCs/"+var+ var_type+ str(plevel) +'_Obsnr'+ str(Obs)+ '_propexcl'+str(prop_speed_excl) +"_PCs.csv"
    else: csv_name= Mediadir + "ERA5_STARS/PL_centred_fields/PCs/"+var+ var_type + str(plevel)+ '_'+Obs+ '_propexcl'+str(prop_speed_excl)+"_PCs.csv"
    
    print(csv_name)
    
    col_names=[var+ var_type +'_PC'+str(iEOFnr+1) for iEOFnr in range(EOFnr) ]
    df= pd.DataFrame(PC, columns= col_names)
    df['ID']= ds_now.ID
    df['Obs']= ds_now.Obs
    df= df.set_index('ID')
    
    exist= os.path.isfile(csv_name) 
    df.to_csv(csv_name, sep=',' )
        

    

