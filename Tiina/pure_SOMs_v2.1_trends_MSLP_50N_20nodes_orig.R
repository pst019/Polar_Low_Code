# set up USER choices

out_dir <- getwd()
in_dir  <- "/media/sf_VBox_Ubuntu_Shared_Folder/SOM/Vertical_SOM/"
f_name  <- "msl_erai_0.75deg_6h_20090101-20181231_eag.nc"

setwd(out_dir)

seas   <- "ALL"
Nnodes <- 20
nx=round(sqrt(Nnodes),0)
ny=nx+1

seas
Nnodes

# ==================================================================================== #
# SET UP LIBRARY

library("fpc")
library("tools")
library("ncdf4")
library("kohonen")
library("clValid")
packageVersion("kohonen")

# ==================================================================================== #
# read in netcdf file to do clustering on
setwd(in_dir)
files <- system(paste("ls ",f_name,sep=""), intern = TRUE)
print(files)

var_name <- files
var_name_S=file_path_sans_ext(var_name)  # remove the extension using package "tools"
print(var_name_S)

nc  <- nc_open( var_name, readunlim=FALSE )

msl    <- ncvar_get(nc, varid="MSL")
t      <- ncvar_get(nc, "time")
lon    <- ncvar_get(nc, "nav_lon")
lat    <- ncvar_get(nc, "nav_lat")
tunits <- ncatt_get(nc, "time", "units")
tunits_ref        <- strsplit(tunits$value, " ")
tunits_ref_string <- paste(unlist(tunits_ref)[3], unlist(tunits_ref)[4])
paste("time unit is ....",unlist(tunits_ref)[1],unlist(tunits_ref)[2])

posixct_datetime <- as.POSIXct(t*3600, origin=tunits_ref_string, tz = "GMT") # note the *3600 is specific to "hours since", use 3600/60 for minutes since (MERRA2)
datetime         <- as.numeric(format(posixct_datetime, "%Y%m%d%H"))
mslp_dates       <- as.numeric(as.character(datetime))
mslp_dates[1:3]

nlat <- dim(msl)[1]
nlon <- dim(msl)[2]
nt   <- dim(t)

# this was the problem previously - PETER
# we were not transposing properly
tmp.vec.long <- as.vector(msl)
tmp.mat      <- matrix(tmp.vec.long, nrow = nlon * nlat, ncol = nt)
mslp         <- data.frame(t(tmp.mat))
#########################################
print(paste("obs dataset has dimensions .....", dim(mslp), sep=" "))


# ==================================================================================== #
#### run single #####
set.seed(5) #7
setwd(out_dir)

print("running clustering algorithm ........")
# note should increase rlen to 1000+ but small here (100) for quick testing
K_SOM <- som(X=data.matrix(mslp), grid = somgrid(nx, ny, "hex"), rlen=5000, alpha = c(0.05,0.01), radius = c(5,1), keep.data=T)
print("FIN")

K_SOM_d    <- data.frame(K_SOM$distances)        # get errors
K_SOM_SOM  <- data.frame(K_SOM$codes)            # get clusters gird
temp       <- data.frame(t(K_SOM_SOM))           # transpose
K_SOM_SOMc <- data.frame(unlist(temp))           # concatenate codebook vectors to 1 row
K_SOM_win  <- data.frame(K_SOM$unit.classif)     # get clusters win

# write some node-numbers into a file
winning_nodes <- cbind(as.character(mslp_dates),K_SOM_win)            
colnames(winning_nodes) <- c("date","node")
write.table(winning_nodes, file=paste("SOM_",seas,"_",nx,"_",ny,"_master.txt", sep=""),row.names=F, col.names=T, quote=F)
print("file written as ....", quote=F)
print(paste("SOM_",seas,"_",nx,"_",ny,"_5050.txt", sep=""))

# write distances (errors) into a file
error_distance <-  cbind(as.character(mslp_dates),K_SOM_d)
colnames(winning_nodes) <- c("date","distances")
write.table(error_distance, file=paste("SOM_",seas,"_",nx,"_",ny,"_distances.txt", sep=""),row.names=F, col.names=T, quote=F)
print("file written as ....", quote=F)
print(paste("SOM_",seas,"_",nx,"_",ny,"_distances.txt", sep=""))


# this may fail on EASE grid because lon and lat are 2d (expects 1d)
lon_1d <- seq(1,nlon,1)
lat_1d <- seq(1,nlat,1)
#node_nc_LAT  <- ncdim_def( "Lat", "degreesN", lat_1d)
node_nc_LON  <- ncdim_def( "Lon", "testing", lon_1d)
node_nc_LAT  <- ncdim_def( "Lat", "testing", lat_1d)
node_nc_N    <- ncdim_def( "nodeN", "node number", seq(1,nx*ny))
node_nc_data <- ncvar_def( "node1", "mslp,", list(node_nc_LON, node_nc_LAT, node_nc_N))
nc_new <- nc_create( "node_nc.nc", node_nc_data)
# put values into nc file and close
ncvar_put(nc=nc_new, varid=node_nc_data, K_SOM_SOMc[,1])
nc_close(nc_new)
print("nc file written")

