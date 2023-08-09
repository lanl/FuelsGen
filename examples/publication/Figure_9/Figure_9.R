# Canopy
chm = raster::raster(here::here("examples/CHM/LA_CHM.tif"))

# Road
chm_road = raster::raster(here::here("examples/CHM/LA_CHM_road.tif"))
raster::extent(chm_road) = raster::extent(chm)
raster::values(chm_road)[raster::values(chm_road) == 0] = -1
raster::values(chm_road)[raster::values(chm_road) != -1] = 0
chm_road = modify_raster(chm_road,type='zero')

chm = modify_raster(chm,type='regular')

# 9a
raster::plot(chm$raster)

# 9b
raster::plot(chm_road$raster,col=colorRampPalette(c("black","white"))(11))

fuel = gen_fuels(chm$dimX, chm$dimY,
                 density = .025,           # shrubs per unit area before downsampling
                 heterogeneity = 5,        # level of heterogeneity in shrub placement
                 heterogeneity.scale = 1,  # scale of the mean-zero GP realization
                 radius = 1, sd_radius = 0,# normal distribution parameters for shrub radius
                 height = NULL, sd_height = 0,# normal distribution parameters for shrub height
                 reps=1,                   # number of random maps to generate
                 GP.init.size=100,         # How densely to sample GP
                 seed=10,                  # random seed for reproducibility)
                 verbose = F)  
# 9c
plot_fuels(fuel)

B = c(1,10)
fuel = gen_fuels(chm$dimX, chm$dimY,
                 density = .025,           # shrubs per unit area before downsampling
                 heterogeneity = 5,        # level of heterogeneity in shrub placement
                 heterogeneity.scale = 1,  # scale of the mean-zero GP realization
                 radius = 1, sd_radius = 0,# normal distribution parameters for shrub radius
                 height = NULL, sd_height = 0,# normal distribution parameters for shrub height
                 sp.cov.locs = list(chm$locs[[1]],chm_road$locs[[1]]),
                 # locations where canopy height is known
                 sp.cov.vals = list(chm$vals[[1]],chm_road$vals[[1]]),     
                 # canopy height values
                 sp.cov.scale = B,         # Importance of chm for shrub placement
                 reps=1,                   # number of random maps to generate
                 GP.init.size=100,         # How densely to sample GP
                 seed=10,                  # random seed for reproducibility)
                 verbose = F)                 

# 9d
raster::plot(chm$raster,legend=F)
points(fuel$dat[[1]]$X,fuel$dat[[1]]$Y,cex=.5,pch=16)
