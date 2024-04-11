# This script contains user-facing functions should the fuels generation process be run
# from outside the shiny application.

#' @title Generate fuel layouts
#'
#' @description Generate fuel layouts on a rectangular domain using a handful of parameters to a generative model.
#' @param dimX (0,Inf): The size of the domain in X direction (meters)
#' @param dimY (0,Inf): The size of the domain in X direction (meters)
#' @param density (0,Inf): The mean number of fuels items per unit area (square meter)
#' @param radius (0,Inf): The Average radius of a circular fuel element
#' @param sd_radius (0,Inf): The standard deviation of fuel element radii 
#' @param height (0,Inf): The Average height of a fuel element. Default is NULL indicating we should not sample heights.
#' @param sd_height (0,Inf): The standard deviation of fuel element heights. Default is 0 indicating no variance in heights.
#' @param heterogeneity (0,Inf): How heterogeneous the fuel maps will be. Larger values result in more heterogeneity.
#' @param heterogeneity.scale (0,Inf): How big of an effect should the heterogeneity have? Scaling factor.
#' @param reps (integer): The number of fuel layouts to generate
#' @param seed (integer): Optional random seed for reproducibility
#' @export
#'
gen_fuels = function(dimX, dimY,
                     density, 
                     radius, sd_radius, 
                     height = NULL, sd_height = 0, 
                     heterogeneity, heterogeneity.scale = 1,
                     sp.cov.locs = NULL, sp.cov.vals = NULL, sp.cov.scale = NULL,
                     reps=1, GP.init.size = 1000, seed=NULL, logis.scale=.217622, parallel=F, verbose=T){
    if(verbose){
      cat('Generating',reps,'fuel maps.\n')
      if(is.null(seed))
        cat('For reproducible maps, set the random seed.\n')
      cat('Parameters:\n')
      cat('  density:       ',density,'\n')
      cat('  fuel radius:   ',radius,'\n')
      cat('  fuel radius sd:',sd_radius,'\n')
      cat('  fuel height:   ',height,'\n')
      cat('  fuel height sd:',sd_height,'\n')
      cat('  heterogeneity: ',heterogeneity,'\n')
      cat('  heterogeneity scale: ',heterogeneity.scale,'\n')
    }
    
    theta = c(density, heterogeneity, radius, sd_radius)
    if(!is.null(height))
        theta = c(theta, height, sd_height)
    data = gen_data(theta, dimX, dimY, heterogeneity.scale, sp.cov.locs, sp.cov.vals, sp.cov.scale, reps, GP.init.size, seed, logis.scale, parallel)
    return(data)
}

#' @title Plot fuel layouts
#'
#' @description Plot fuel layouts on a square tiled figure
#' @param data: fuels data object returned from gen_fuels
#' @export
#'
plot.fuelsgen = function(data){
    if(data$reps>0){
        plot.dim = min(5,ceiling(sqrt(data$reps)))
        plot.list = vector(mode='list',length=data$reps)
        hmin = 0; hmax = 0

        for(i in 1:data$reps){
            if(!is.null(data$dat[[i]]$h)){
              tmp = max(data$dat[[i]]$h)
              if(tmp>hmax)
                hmax = tmp 
              tmp = min(data$dat[[i]]$h)
              if(tmp<hmin)
                hmin = tmp 
                plot.list[[i]] = ggplot2::ggplot() +
                    ggforce::geom_circle(ggplot2::aes(x0 = X, y0 = Y, r = r, fill = h), data=data$dat[[i]]) + 
                    ggplot2::labs(x='X (m)',y='Y (m)',fill='height (m)') + 
                    ggplot2::theme_bw() + 
                    ggplot2::theme(plot.margin=ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                     axis.text = ggplot2::element_text(size = 15), 
                                     axis.title = ggplot2::element_text(size = 15),
                                     axis.ticks = ggplot2::element_line(linewidth = 1.5)) + 
                    ggplot2::coord_fixed(xlim = c(-1, data$dimX+1), ylim = c(-1, data$dimY+1))
            } else{
                hmax = 0
                plot.list[[i]] = ggplot2::ggplot() + 
                    ggforce::geom_circle(ggplot2::aes(x0 = X, y0 = Y, r = r), fill=adjustcolor('grey',alpha.f=.5), data=data$dat[[i]]) + 
                    ggplot2::labs(x='X (m)',y='Y (m)') + 
                    ggplot2::theme_bw() + 
                    ggplot2::theme(plot.margin=ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                   axis.text = ggplot2::element_text(size = 15),
                                   axis.title = ggplot2::element_text(size = 15),
                                   axis.ticks = ggplot2::element_line(linewidth = 1.5)) +
                    ggplot2::coord_fixed(xlim = c(-1, data$dimX+1), ylim = c(-1, data$dimY+1))
            }
        }
        
        if(!is.null(data$dat[[1]]$h)){
          return(patchwork::wrap_plots(plot.list,ncol=plot.dim,nrow=plot.dim, guides='collect') &
                   ggplot2::scale_fill_continuous(limits = c(hmin, hmax),breaks = seq(hmin,hmax,length.out=5),
                                                  labels = function(x) sprintf("%.2f", x)))
        } else{
          return(patchwork::wrap_plots(plot.list,ncol=plot.dim,nrow=plot.dim, guides='collect'))
        }
    } else{
        # single map
        if(!is.null(data$dat$h)){
            myplot = ggplot2::ggplot() +
                     ggforce::geom_circle(ggplot2::aes(x0 = X, y0 = Y, r = r, fill = h), data=data$dat[[1]]) + 
                     ggplot2::labs(x='X (m)',y='Y (m)',fill='height (m)') + 
                     ggplot2::theme_bw() + 
                     ggplot2::theme(aspect.ratio = data$dimY/data$dimX) +
                     ggplot2::coord_fixed(xlim = c(-1, data$dimX+1), ylim = c(-1, data$dimY+1))
        } else{
            myplot = ggplot2::ggplot() + 
                     ggforce::geom_circle(ggplot2::aes(x0 = X, y0 = Y, r = r), fill=adjustcolor('grey',alpha.f=.5), data=data$dat[[1]]) + 
                     ggplot2::labs(x='X (m)',y='Y (m)') + 
                     ggplot2::theme_bw() + 
                     ggplot2::theme(aspect.ratio = data$dimY/data$dimX) +
                     ggplot2::coord_fixed(xlim = c(-1, data$dimX+1), ylim = c(-1, data$dimY+1))
        }
        # plot code for IF we move to ellipses instead of circles
        #data$dat$angle = runif(n=nrow(data$dat),min=0,max=360)
        #myplot = ggplot2::ggplot() + 
        #         ggforce::geom_ellipse(ggplot2::aes(x0 = X, y0 = Y, 
        #                                            a = .8*r, b=1.2*r, 
        #                                            angle=angle, fill=h), data=data$dat) + 
        #         ggplot2::labs(x='X (m)',y='Y (m)',fill='height (m)')
        return(myplot)
    }
}

#' @title Fix raster
#'
#' @description Change raster extent and coordinates to work with the generative model. Also scale raster values to lie in [-1,1]. Return a list containing the raster, coordinates, and values all in the correct format for the model. 
#' @param data: raster object
#' @export
#'
modify_raster = function(data,type='regular'){
  dimX = data@extent@xmax-data@extent@xmin
  dimY = data@extent@ymax-data@extent@ymin
  raster::extent(data) = c(0,dimX,0,dimY)
  resolution = raster::res(data)
  coords = raster::coordinates(data)
  locs = list()
  locs$x = unique(coords[,1])
  locs$y = rev(unique(coords[,2])) # y's must be in increasing order for pracma::interp2
  vals = raster::values(data,format='matrix')
  if(type=='regular'){
    vals[vals<0] = 0
    vals = (vals - min(vals)) / (max(vals) - min(vals))
    vals = 2*vals - 1 # scale to [-1,1] so that low canopy reduces probability of shrub
  }
  raster::values(data) = vals
  # flip y dim to match locs$y
  vals = vals[nrow(vals):1,]
  
  return(list(raster=data,dimX=dimX,dimY=dimY,locs=list(locs),vals=list(vals)))
}