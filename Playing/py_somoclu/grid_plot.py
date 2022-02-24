import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
# Defining projections for plots
# Projection for map
projm = ccrs.Stereographic(central_latitude=70, central_longitude= -40)
# Projection to use for transform
proj = ccrs.PlateCarree()

# Data paths
datadir = '/run/media/the038/external/data/racmo/SMB/'
fname = 'smb_rec.1979-2017.nc'
# Path to gridfile
gpath = '/run/media/the038/external/data/racmo/Icemask_Topo_Iceclasses\
_lon_lat_average_1km.nc'

# Reading grid
grid = xr.open_dataset(gpath)
lat = grid['LAT']
lon = grid['LON']
# Reading data
d = xr.open_dataset(f'{datadir}{fname}')
# Selecting one year from the data
year = 1979
d = d.sel(time = (d.time.dt.year == year))
# Reading data into memory
d = d['SMB_rec'].data 

# Dimensions of grid
nrows = 4
ncols = 5

fig, axs = plt.subplots(nrows,ncols,subplot_kw = dict(projection = projm))
i = 0
for row in range(nrows):
    for col in range(ncols):
        axs[row,col].contourf(lon,lat,d[i,:,:],transform = proj)
        axs[row,col].coastlines(resolution = '10m')
        i += 1
plt.show()
