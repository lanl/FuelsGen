library(terra)
library(here) #optional, I like this instead of setwd()

## Read in data
DEM <- rast(here("examples/bandelier/data/shape_files/DEM.tif")) #replace here(...) with your file location
# aoi <- vect(here("examples/bandelier/data/shape_files/dome_aoi_shape.shp"))

aoi = vect(here("examples/bandelier/data/shape_files/dome_aoi_new/dome_aoi_new.shp"))
bbox = vect(here("examples/bandelier/data/shape_files/bbox_for_modeling/bbox_for_modeling.shp"))

# put a 20m buffer around the bounding box for computing slope/aspect from DEM
ext_bbox = ext(bbox)
ext_bbox[1] <- ext_bbox[1] - 20
ext_bbox[2] <- ext_bbox[2] + 20
ext_bbox[3] <- ext_bbox[3] - 20
ext_bbox[4] <- ext_bbox[4] + 20
DEM = crop(DEM, ext_aoi)

par(mfrow=c(2,2))
plot(bbox)
# DEM = project(DEM,bbox)
plot(DEM,add=T)
plot(aoi,add=T)

cover_20m = rast(here('examples/bandelier/data/shape_files/Shrub_Cover_20m.tif'))
ext(cover_20m)
cover_20m = trim(cover_20m,value=NaN)
ext(cover_20m)
cover_20m = project(cover_20m,bbox)
plot(bbox)
plot(cover_20m,add=T)
plot(aoi,add=T)

# calculate slope and aspect for each 1m grid cell
slope_1m <- terra::terrain(DEM, v="slope", unit="degrees")
plot(bbox); plot(slope_1m,add=T); plot(aoi,add=T)
aspect_1m <- terra::terrain(DEM, v="aspect", unit="degrees")
aspect_1m <- (1 - cos((aspect_1m - 30) * pi/180))/2 #convert to southwestness
plot(bbox); plot(aspect_1m,add=T); plot(aoi,add=T)

# aggregate to 20m res
slope_20m <- terra::aggregate(slope_1m, fact=20, fun="mean")
slope_20m = crop(slope_20m,bbox)
aspect_20m <- terra::aggregate(aspect_1m, fact=20, fun="mean")
aspect_20m = crop(aspect_20m,bbox)
elevation_20m = terra::aggregate(DEM, fact=20, fun="mean")
elevation_20m = crop(elevation_20m,bbox)

# make sure these all have the same cell dimensions
dim(cover_20m)
dim(slope_20m)
dim(aspect_20m)
dim(elevation_20m)

par(mfrow=c(2,2))
plot(bbox);plot(slope_20m,main='slope 20m',add=T);points(plots_center,cex=.25);plot(aoi,add=T)
plot(bbox);plot(aspect_20m,main='aspect 20m',add=T);points(plots_center,cex=.25);plot(aoi,add=T)
plot(bbox);plot(elevation_20m,main='elevation 20m',add=T);points(plots_center,cex=.25);plot(aoi,add=T)
plot(bbox);plot(cover_20m,main='NAIP cover 20m',add=T);points(plots_center,cex=.25);plot(aoi,add=T)
all.equal(crds(elevation_20m),crds(slope_20m))
all.equal(crds(aspect_20m),crds(slope_20m))
all.equal(crds(aspect_20m),crds(cover_20m)) 

# get topo for plot locations
# NOTE 74 and 73 are out of order here and in allpts_buffer.shp we must account for that when reading from csv's
plot_topo = read.csv(here('examples/bandelier/data/plot_topo.csv'))
# NOTE we do not have scan information for plots 1 and 143, remove these
plot_topo = plot_topo %>% filter(!(Plot %in% c(1,143)))
load(here('examples/bandelier/results/plot_coverage.RData'))
# make sure that plot numbers match
all.equal(plot_topo$Plot,coverage$Plot)
plot_topo$coverage = coverage$coverage

# plot(rast(here('shape_files/topo_rasters.tif')))
# plot(rast(here('shape_files/predicted_nirrgb.tif')))
# plot(rast(here('shape_files/nirrgb_1m_mask.tif')))

# get locations of plots
plots_buffer = vect(here('examples/bandelier/data/shape_files/allpts_buffer.shp'))
# remove scans 1 and 143 from this because we don't have data for them
which.delete = which(values(plots_buffer)$Plot %in% c(1,143))
plots_buffer = plots_buffer[-which.delete]

# convert to same coordinates
plots_buffer = project(plots_buffer,DEM)
plots_center = centroids(plots_buffer)
plots_center_crds = crds(plots_center)
plot_topo$elevation = extract(DEM,plots_center_crds)[,1]
plot_topo$x_crd = plots_center_crds[,1]
plot_topo$y_crd = plots_center_crds[,2]
# get elevation for plot locations from DEM

# what are the NAN values is that a boundary?
grid_topo = data.frame(slope=values(slope_20m),
                       aspect=values(aspect_20m),
                       elevation=values(elevation_20m)[,1], # these have 202 NAN
                       coverage = values(cover_20m)[,1],
                       x_crd=crds(elevation_20m)[,1], # these don't include the NAN locations
                       y_crd=crds(elevation_20m)[,2])
grid_topo = grid_topo[complete.cases(grid_topo),]

save(grid_topo,plot_topo,file=here('examples/bandelier/results/topology_df.RData'))
