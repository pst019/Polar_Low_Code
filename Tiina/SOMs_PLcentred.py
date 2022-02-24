
import numpy as np
import somoclu
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
## Defining projections for plots
#projm = ccrs.Stereographic(central_latitude=70, central_longitude= -40)
#proj = ccrs.PlateCarree()
#
#datadir = '/run/media/the038/external/data/racmo/SMB/'
##fname = 'smb_rec.1979-2017.nc'
#fname = 'smb_rec_mm.1979-2017.nc'
#
## Reading data
#d = xr.open_dataset(f'{datadir}{fname}')
## Selecting only one year for test purposes
##d = d.sel(time = (d.time.dt.year < 2002)*(d.time.dt.year > 2000))
#
## Reading data into memory. This will not work for large datasets.
## Need to find a workaround for that.
#d = d['SMB_rec'].data
#
#
#
## Reshaping data
#orig_x = d.shape[1]
#orig_y = d.shape[2]
#d = np.reshape(d, (d.shape[0], orig_x*orig_y))
#print(d.shape)


# Defining size of SOM
nrows = 3
ncols = 3


datadir= "/media/pst019/PatsOrange/ERA5_STARS/PL_centred_fields_smooth-tracks/"
var= 't'
level='_850'
ano= True
if ano== True: var_full= var+ '_ano'
else: var_full= var

file_ending= "_allObs_SOMprep_track-smth-1E-3_dur6.nc"

d = xr.open_dataset(f'{datadir}{var}{level}{file_ending}')



max_lev = 8
#max_lev = np.max(d)
min_lev = -8
#min_lev = np.min(d)
nlevels = 17
levels = np.linspace(min_lev,max_lev,nlevels)


d= d[var_full]

dim_time= len(d.time)
dim_x= len(d.x)
dim_y= len(d.y)
 
d_flat = np.reshape(d.data, (dim_time, dim_x*dim_y ) )
print(d_flat.shape)


# Creating instance of SOM
som = somoclu.Somoclu(ncols,nrows ) #,compactsupport= False)

# Training SOM on data.
som.train(d_flat, epochs=300)



#fig= plt.figure(1, figsize= (2.5*nrows,2.5*ncols +1) )        
fig, axs = plt.subplots(nrows,ncols, figsize = (10,10))
#d = np.reshape(d,(dim_time, dim_x, dim_y))



for row in range(nrows):
    temp_row = 0
    for col in range(ncols):
        bool_array = [np.array_equal(b,[col,row]) for b in som.bmus]
        print(sum(bool_array))
        t_dat = np.mean(d[bool_array,:,:],axis = 0)
#        ax1= plt.subplot(ncols, nrows, (col*ncols + row+1 ) )
        
#        cf = ax1.contourf(d.x,d.y,t_dat, cmap = 'seismic',extend = 'both')
        cf = axs[row,col].contourf(d.x,d.y,t_dat, levels= levels, cmap = 'RdBu_r',extend = 'both')
        
#        axs[row,col].coastlines(resolution = '10m')
#        plt.colorbar(cf,ax = axs[row,col])

#plt.contourf(d.x, d.y, t_dat)      
        
fig.subplots_adjust(bottom=0.11)
cbar_ax = fig.add_axes([0.09, 0.06, 0.4, 0.015])
cb= fig.colorbar(cf, cax=cbar_ax, orientation="horizontal")

#plt.show()



