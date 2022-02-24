import matplotlib
matplotlib.use("Agg")
import numpy as np
import somoclu
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
# Defining projections for plots
projm = ccrs.Stereographic(central_latitude=70, central_longitude= -40)
proj = ccrs.PlateCarree()

datadir = '/run/media/the038/external/data/racmo/SMB/'
#fname = 'smb_rec.1979-2017.nc'
fname = 'smb_rec_mm.1979-2017.nc'

# Reading data
d = xr.open_dataset(f'{datadir}{fname}')
# Selecting only one year for test purposes
#d = d.sel(time = (d.time.dt.year < 2002)*(d.time.dt.year > 2000))

# Reading data into memory. This will not work for large datasets.
# Need to find a workaround for that.
d = d['SMB_rec'].data



# Reshaping data
orig_x = d.shape[1]
orig_y = d.shape[2]
d = np.reshape(d, (d.shape[0], orig_x*orig_y))
print(d.shape)
# Defining size of SOM
nrows = 3
ncols = 4
# Creating instance of SOM
som = somoclu.Somoclu(ncols,nrows,compactsupport= False)

# Training SOM on data.
som.train(d,epochs=10)
# Viewing best matches
labels = range(d.shape[0])
#som.view_umatrix(bestmatches=True,  labels=labels)
#print(som.bmus.shape)

# Reading grid for plotting
gpath = '/run/media/the038/external/data/racmo/Icemask_Topo_Iceclasses\
_lon_lat_average_1km.nc'
grid = xr.open_dataset(gpath)
lat = grid['LAT']
lon = grid['LON']



fig, axs = plt.subplots(nrows,ncols, subplot_kw = dict(projection = projm),figsize = (10,10))
d = np.reshape(d,(d.shape[0],orig_x,orig_y))
max_lev = 20
#max_lev = np.max(d)
min_lev = -20
#min_lev = np.min(d)
nlevels = 41
levels = np.linspace(min_lev,max_lev,nlevels)

for row in range(nrows):
    temp_row = 0
    for col in range(ncols):
        bool_array = [np.array_equal(b,[col,row]) for b in som.bmus]
        print(sum(bool_array))
        t_dat = np.mean(d[bool_array,:,:],axis = 0)
        cf = axs[row,col].contourf(lon,lat,t_dat,transform = proj,levels = levels,
                cmap = 'seismic',extend = 'both')
        axs[row,col].coastlines(resolution = '10m')
        plt.colorbar(cf,ax = axs[row,col])

        


plt.show()

