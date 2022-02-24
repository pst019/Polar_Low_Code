# set up USER choices
# on Backup

out_dir <- "/media/pst019/Backup/ERA5_STARS/PL_centred_fields_smooth-tracks/SOM/"
#in_dir  <- "/media/pst019/1692A00D929FEF8B/ERA5_STARS/PL_centred_fields/"
in_dir  <- "/media/pst019/Backup/ERA5_STARS/PL_centred_fields_smooth-tracks/"

var <- "t"
ano <- TRUE
#ano <- FALSE

#f_name  <- paste(var, "_850_Obsnr1.nc", sep="")
#f_name  <- paste(var, "_850_mature_prep.nc", sep="")
#f_name  <- paste(var, "_850_allObs_SOMprep.nc", sep="")
#f_name  <- paste(var, "_850_allObs_SOMprep_vel3_dur6.nc", sep="")
#f_name  <- paste(var, "_700_allObs_SOMprep_track-smth-1E-3_dur6.nc", sep="")

#tests:
#f_name  <- paste(var, "_850_Obs1_SOMprep_track-smth-1E-3_dur6.nc", sep="")
f_name  <- paste(var, "_850_Obsmature_SOMprep_track-smth-1E-3_dur6.nc", sep="")
#f_name  <- paste(var, "_850_Obslast_SOMprep_track-smth-1E-3_dur6.nc", sep="")

#f_name  <- paste(var, "_850_allObs_primary_SOMprep_track-smth-1E-3_dur12_Umax20.nc", sep="") 
#f_name  <- paste(var, "_850_allObs_SOMprep_track-smth-1E-3_dur6_matchdist75_Nmatches5.nc", sep="") 


setwd(out_dir)

seas   <- "ALL"
#Nnodes <- 4
x_nodes=3 # round(sqrt(Nnodes),0)
y_nodes=3 #x_nodes+1

seas
#Nnodes

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

nc  <- nc_open(f_name, readunlim=FALSE )

if (ano == TRUE){
dat    <- ncvar_get(nc, varid=paste(var, "_ano", sep="") )
} else {
dat    <- ncvar_get(nc, varid=var)
}

tim      <- ncvar_get(nc, "time")
x    <- ncvar_get(nc, "x")
y    <- ncvar_get(nc, "y")
PLnr <- ncvar_get(nc, "PLnr")
tunits <- ncatt_get(nc, "time", "units")
tunits_ref        <- strsplit(tunits$value, " ")
tunits_ref_string <- paste(unlist(tunits_ref)[3], unlist(tunits_ref)[4])
paste("time unit is ....",unlist(tunits_ref)[1],unlist(tunits_ref)[2])

posixct_datetime <- as.POSIXct(tim*3600, origin=tunits_ref_string, tz = "GMT") # note the *3600 is specific to "hours since", use 3600/60 for minutes since (MERRA2)
datetime         <- as.numeric(format(posixct_datetime, "%Y%m%d%H"))
dates       <- as.numeric(as.character(datetime))
dates[1:3]

nx <- dim(x)
ny <- dim(y)
ntim   <- dim(tim)

# this was the problem previously - PETER
# we were not transposing properly
tmp.vec.long <- as.vector(dat)
tmp.mat      <- matrix(tmp.vec.long, nrow = nx * ny, ncol = ntim)
dat_mat         <- data.frame(t(tmp.mat))
#########################################
print(paste("obs dataset has dimensions .....", dim(dat_mat), sep=" "))


# ==================================================================================== #
#### run single #####
#set.seed(5) # random number generation that can be reproduced
set.seed(1) # random number generation that can be reproduced

setwd(out_dir)
f_name_S=file_path_sans_ext(f_name) #file name short, without .nc
if (ano == TRUE){
f_name_S= paste(var,"_ano_", sub("^.*?_", "", f_name_S), sep="") }

print("running clustering algorithm ........")
# note should increase rlen to 1000+ but small here (100) for quick testing
K_SOM <- som(X=data.matrix(dat_mat), grid = somgrid(x_nodes, y_nodes, "hex"), rlen=10000, alpha = c(0.05,0.01), radius = c(5,1), keep.data=T)

print("FIN")

K_SOM_d    <- data.frame(K_SOM$distances)        # get errors
K_SOM_SOM  <- data.frame(K_SOM$codes)            # get clusters gird
temp       <- data.frame(t(K_SOM_SOM))           # transpose
K_SOM_SOMc <- data.frame(unlist(temp))           # concatenate codebook vectors to 1 row
K_SOM_win  <- data.frame(K_SOM$unit.classif)     # get clusters win

# write some node-numbers into a file
#winning_nodes <- cbind(as.character(dates),K_SOM_win)            
#colnames(winning_nodes) <- c("date","node")

#f_name_node_nr= paste(f_name_S,"_x",x_nodes,"_y",y_nodes,"_node_nr.txt", sep="")
#write.table(winning_nodes, file=f_name_node_nr, row.names=F, col.names=T, quote=F)
#print("file written as ....", quote=F)
#print(f_name_node_nr)

# write some node-numbers into a file
winning_nodes <- cbind(as.character(PLnr), as.character(dates), K_SOM_win)            
colnames(winning_nodes) <- c("PLnr", "date","node")

f_name_node_nr= paste(f_name_S,"_x",x_nodes,"_y",y_nodes,"_node_nr.txt", sep="")
write.table(winning_nodes, file=f_name_node_nr, row.names=F, col.names=T, quote=F)
print("file written as ....", quote=F)
print(f_name_node_nr)

# write distances (errors) into a file
error_distance <-  cbind(as.character(dates),K_SOM_d)
colnames(winning_nodes) <- c("date","distances")
f_name_dist= paste(f_name_S,"_x",x_nodes,"_y",y_nodes,"_distances.txt", sep="")

write.table(error_distance, file=f_name_dist, row.names=F, col.names=T, quote=F)
print("file written as ....", quote=F)
print(f_name_dist)


# this may fail on EASE grid because lon and lat are 2d (expects 1d)
x_vec <- seq(1,nx,1)
y_vec <- seq(1,ny,1)
node_nc_x  <- ncdim_def( "x", "Prop dir [km]", x)
node_nc_y  <- ncdim_def( "y", "Orth dir [km]", y)
node_nc_N    <- ncdim_def( "SOM", "SOM number", seq(1,x_nodes*y_nodes))
node_nc_data <- ncvar_def( "field", var, list(node_nc_x, node_nc_y, node_nc_N))

f_name_nc= paste(f_name_S,"_x",x_nodes,"_y",y_nodes,"_cluster.nc", sep="")
nc_new <- nc_create( f_name_nc, node_nc_data)
# put values into nc file and close
ncvar_put(nc=nc_new, varid=node_nc_data, K_SOM_SOMc[,1])
nc_close(nc_new)
print("nc file written")

