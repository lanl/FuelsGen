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
#' @param intensity.mean (-Inf,Inf): Mean of the intensity function
#' @param intensity.scale (0,Inf): Standard deviation of the intensity function.
#' @param reps (integer): The number of fuel layouts to generate
#' @param GP.init.size (integer): Number of grid cells in each direction for discrete representation of the intensity function.
#' @param seed (integer): Optional random seed for reproducibility
#' @param parallel (logical): Should the fuels be generated in parallel
#' @param verbose (logical): Print output related to generated fuels
#' @export
#'
gen_fuels = function(dimX, dimY, density=.1/dimX/dimY,
                     radius = 1, sd_radius = .1, 
                     height = NULL, sd_height = 0, 
                     heterogeneity = mean(c(dimX,dimY))/3, 
                     intensity.scale = 1,
                     sp.cov.locs = NULL, sp.cov.vals = NULL, sp.cov.scale = NULL,
                     reps=1, GP.init.size = 28, seed=NULL, sample.method='direct', verbose=T){
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
      # cat('  intensity mean:',intensity.mean,'\n')
      cat('  intensity sd: ',intensity.scale,'\n')
    }
    
    theta = c(density, heterogeneity, radius, sd_radius, 
              intensity.scale)
    if(!is.null(height))
      theta = c(theta, height, sd_height)
    data = gen_data(theta, dimX, dimY,
                    sp.cov.locs, sp.cov.vals, sp.cov.scale, reps, 
                    GP.init.size, seed, sample.method)
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

#' @title Diagnostic plots for mcmc calibration
#'
#' @description Plot a fuelsgen_mcmc object
#' @param mcmc: fuelsgen_mcmc object from fuelsgen::mcmc_MH_adaptive()
#' @param burn: Number of initial samples to discard in plots
#' @param type: 'trace' or 'density'. 'trace' gives trace plots for the parameters, while 'density' plots the KDE and prior over the prior range
#' @export
#'
plot.fuelsgen_mcmc = function(mcmc,burn=0,type='trace'){
  par(mfrow=c(2,3)) # 5 parameters in the calibration
  names = mcmc$par_names
  nsamp = nrow(mcmc$trace)
  lb = mcmc$prior$prior_params$lb
  ub = calib$prior$prior_params$ub
  # convert sigma to precision
  sig_id = which(names=='sigma')
  mcmc$trace[,sig_id] = 1/(mcmc$trace[,sig_id]^2)
  ub[sig_id] = 1.1*max(mcmc$trace[,sig_id])
  names[sig_id] = 'precision'
  if(type=='trace'){
    for(i in 1:length(names)){
      plot(mcmc$trace[(burn+1):nsamp,i],type='l',
           ylab=names[i],ylim=c(lb[i],ub[i]))
    }
  } else if(type=='density'){
    for(i in 1:length(names)){
      plot(density(mcmc$trace[,i],
                   from = lb[i],to = ub[i]),
           xlab=names[i],main='',ylab='density')
      
      # abline(v=mcmc$prior$theta_est[i])
      switch(mcmc$prior$prior_params$dist[i],
              gamma = {
                if(names[i]=='lambda'){
                  curve(dgamma(x,mcmc$prior$prior_params$lambda_shape,mcmc$prior$prior_params$lambda_rate),add=T,lty=2)
                } else{
                  curve(dgamma(x,mcmc$prior$prior_params$prec_shape,mcmc$prior$prior_params$prec_rate),add=T,lty=2)
                }
              },
              uniform = {
                curve(dunif(x,0,mcmc$prior$prior_params$rho_max),add=T,lty=2)
              },
              truncnorm = {
                curve(truncdist::dtrunc(x,'norm',0,3,mcmc$prior$prior_params$mu_mean,mcmc$prior$prior_params$mu_sd),add=T,lty=2)
              },
              hcauchy = {
                curve(LaplacesDemon::dhalfcauchy(x,mcmc$prior$prior_params$I_var_scale),add=T,lty=2)
              })
      # lines(priorfunc)
    }
  } else{
    stop('Currently supported types are trace and density')
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