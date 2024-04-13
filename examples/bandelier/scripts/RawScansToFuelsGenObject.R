library(here)
library(dplyr)
files = list.files(here('examples/bandelier/data/Shrubs/'))
nfiles = length(files)
# NOTE 74 and 73 are out of order here and in allpts_buffer.shp and plot_topo.csv, 
# we must account for that when reading from csv's by putting them in the wrong order here as well
# ..., 72, 74, 73, 75, ...
tmp = files[50]
files[50] = files[51]
files[51] = tmp; rm(tmp)
files

dat = vector(mode='list',length=nfiles)
nshrub = scan_number = numeric(nfiles)
for(i in 1:nfiles){
  raw = read.csv(paste0(here('examples/bandelier/data/Shrubs/'),files[i]))
  dat[[i]] = data.frame(X=raw$X+10,Y=raw$Y+10,r=sqrt(raw$convhull_area/pi),h=raw$Z)
  
  # filter shrubs with r<.25
  dat[[i]] = dat[[i]] %>% filter(r>.25)
  
  # filter size to 20x20m square inside the 15m radius
  dat[[i]] = dat[[i]] %>% filter(between(X,0,20))
  dat[[i]] = dat[[i]] %>% filter(between(Y,0,20))
  nshrub[i] = nrow(dat[[i]])
  scan_number[i] = as.integer(stringr::str_sub(files[i],7,10))
}
save(nshrub,file=here('examples/bandelier/results/plot_nshrub.RData'))
scan_number
hist(nshrub)
boxplot(nshrub)
min(nshrub)

data = list(dat=dat,reps=length(dat),dimX=20,dimY=20)
mets = get_mets(fuels=data)
mets$mets[,2] = jitter(mets$mets[,2])
save(data,mets,nshrub,scan_number,file=here('examples/bandelier/results/plot_data.RData'))

# shrub coverage
coverage = numeric(length(data$dat))
for(k in 1:length(coverage)){
  cat(k,' ')
  dx = .1
  dy = .1
  stepsize = 0.01
  xstart = 0
  dimX = 20
  dimY = 20
  area_fractions = matrix(ncol=dimX / dx, nrow=dimY / dy)
  
  for (i in 1:(dimX / dx)) {
    ystart = 0
    for (j in 1:(dimY / dx)) {
      my_grid = pracma::meshgrid(seq(xstart, xstart + dx - stepsize, by=stepsize),
                                 seq(ystart, ystart + dy - stepsize, by=stepsize))
      mc_d = sqrt(plgp::distance(cbind(c(my_grid$X), c(my_grid$Y)), data$dat[[k]][,c("X", "Y")]))
      
      area_fractions[j,i] = mean(apply(mc_d < data$dat[[k]]$r, 1, any))
      ystart = ystart + dy
    }
    xstart = xstart + dx
  }
  total_area = mean(area_fractions) * dimX * dimY
  coverage[k] = total_area / dimX / dimY
  # area_with_overlap = sum(pi*data$dat[[k]]$r^2)
  # save(coverage,file=here('examples/bandelier/results/plot_coverage.RData'))
}
hist(coverage)
coverage[is.nan(coverage)] = 0
coverage = data.frame(Plot = scan_number,coverage = coverage)
save(coverage,file=here('examples/bandelier/results/plot_coverage.RData'))
