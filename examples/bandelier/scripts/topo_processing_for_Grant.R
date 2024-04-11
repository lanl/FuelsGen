library(terra)
library(here) #optional, I like this instead of setwd()

## Read in data
DEM <- rast(here("examples/bandelier/data/shape_files/DEM.tif")) #replace here(...) with your file location
aoi <- vect(here("examples/bandelier/data/shape_files/dome_aoi_shape.shp"))
# plot(aoi)
# aoi <- project(aoi,DEM)
plot(DEM)
# plot(aoi) # we don't really need this now, we are going to generate fuels on the entire domain
# put a 50m buffer around the plot
# ext_aoi = ext(aoi)
# ext_aoi[1] <- ext_aoi[1] - 50
# ext_aoi[2] <- ext_aoi[2] + 50
# ext_aoi[3] <- ext_aoi[3] - 50
# ext_aoi[4] <- ext_aoi[4] + 50
# DEM <- crop(DEM,ext_aoi)

cover_20m = rast(here('examples/bandelier/data/shape_files/Shrub_Cover_20m.tif'))
# ext(cover_20m)
# cover_20m = trim(cover_20m,value=NaN)
# ext(cover_20m)
# plot(cover_20m)

# DEM is on a much larger domain than Shrub_Cover_20m.tif so crop it
DEM = crop(DEM, ext(cover_20m))

# calculate slope and aspect for each 1m grid cell
slope_1m <- terra::terrain(DEM, v="slope", unit="degrees")
aspect_1m <- terra::terrain(DEM, v="aspect", unit="degrees")
aspect_1m <- (1 - cos((aspect_1m - 30) * pi/180))/2 #convert to southwestness

# aggregate to 20m res
slope_20m <- terra::aggregate(slope_1m, fact=20, fun="mean")
aspect_20m <- terra::aggregate(aspect_1m, fact=20, fun="mean")
elevation_20m = terra::aggregate(DEM, fact=20, fun="mean")

# make sure these all have the same cell dimensions
dim(cover_20m)
dim(slope_20m)
dim(aspect_20m)
dim(elevation_20m)

# get topo for plot locations
# NOTE 74 and 73 are out of order here and in allpts_buffer.shp we must account for that when reading from csv's
plot_topo = read.csv(here('examples/bandelier/data/plot_topo.csv'))

# plot(rast(here('shape_files/topo_rasters.tif')))
# plot(rast(here('shape_files/predicted_nirrgb.tif')))
# plot(rast(here('shape_files/nirrgb_1m_mask.tif')))

# get locations of plots
plots_buffer = vect(here('examples/bandelier/data/shape_files/allpts_buffer.shp'))
# convert to same coordinates
plots_buffer = project(plots_buffer,DEM)
plots_center = centroids(plots_buffer)
plots_center_crds = crds(plots_center)
plot_topo$elevation = extract(DEM,plots_center_crds)[,1]
plot_topo$x_crd = plots_center_crds[,1]
plot_topo$y_crd = plots_center_crds[,2]
# get elevation for plot locations from DEM
par(mfrow=c(2,2))
plot(slope_20m,main='slope 20m');points(plots_center,cex=.25)
plot(aspect_20m,main='aspect 20m');points(plots_center,cex=.25)
plot(elevation_20m,main='elevation 20m');points(plots_center,cex=.25)
plot(cover_20m,main='NAIP cover 20m');points(plots_center,cex=.25)

# what are the NAN values is that a boundary?
grid_topo = data.frame(slope=values(slope_20m),
                       aspect=values(aspect_20m),
                       elevation=values(elevation_20m)[,1], # these have 202 NAN
                       coverage = values(cover_20m)[,1],
                       x_crd=crds(elevation_20m)[,1], # these don't include the NAN locations
                       y_crd=crds(elevation_20m)[,2])
grid_topo = grid_topo[complete.cases(grid_topo),]

save(grid_topo,plot_topo,file=here('examples/bandelier/results/topology_df.RData'))
