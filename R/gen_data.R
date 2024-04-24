relu = function(x,m=0,bound=0){
  x[x<=bound] = m
  return(x)
}

# Generate 'reps' fuel maps with parameters theta over a domain of size [0,dimX]x[0,dimY]
gen_data = function(theta, dimX, dimY, heterogeneity.scale = 1, 
                    sp.cov.locs = NULL, sp.cov.vals = NULL, sp.cov.scale = NULL, 
                    reps = 1, GP.init.size=1000, seed = NULL, 
                    I.transform='logistic',
                    logis.scale=.217622, parallel=F)
{
  if(!is.null(seed))
    set.seed(seed)
  lambda = theta[4]
  rho = theta[1]
  mu_r = theta[2]
  sigma_r = theta[3]
  if(length(theta)>4){
    mu_h = theta[5]
    sigma_h = theta[6]
  } else{
    mu_h = NULL
  }
  
  W = spatstat.geom::owin(c(0,dimX),c(0,dimY), mask=matrix(TRUE, GP.init.size,GP.init.size))
  gamma = spatstat.random::rGRFgauss(W = W, mu = 0, var = heterogeneity.scale, scale = rho,nsim = reps)
  n_plus = rpois(reps, lambda * dimX * dimY)
  dat = vector(mode="list", length=reps)
  for(i in 1:reps){
    if(I.transform=='logistic'){
      gamma[[i]]$v = plogis(gamma[[i]]$v - mean(gamma[[i]]$v),0,logis.scale)
    } else if(I.transform=='exp'){
      gamma[[i]]$v = exp(gamma[[i]]$v - mean(gamma[[i]]$v))
    } else if(I.transform=='relu'){
      gamma[[i]]$v = relu(gamma[[i]]$v - mean(gamma[[i]]$v))
    }
    
    #### sample pixels directly ####
    nn = n_plus[i]
    if(!is.finite(nn))
      stop(paste("Unable to generate Poisson process with a mean of",
                 nn, "points"))
    if(nn>0){
      dx = gamma[[i]]$xstep/2
      dy = gamma[[i]]$ystep/2
      df = as.data.frame(gamma[[i]])
      npix = nrow(df)
      lpix = df$value
      ii = sample.int(npix, size=nn, replace=TRUE, prob=lpix)
      xx = df$x[ii] + runif(nn, -dx, dx)
      yy = df$y[ii] + runif(nn, -dy, dy)
      dat[[i]] = data.frame(X=xx,Y=yy)
      dat[[i]]$r = truncdist::rtrunc(nn,'norm',a=0,b=Inf,mu_r,sigma_r)
      if(!is.null(mu_h)){
        dat[[i]]$h = truncdist::rtrunc(nn,'norm',a=0,b=Inf,mu_h,sigma_h)
      }
    } else{
      dat[[i]] = data.frame(X=numeric(),Y=numeric(),r=numeric())
      if(!is.null(mu_h)){
        dat[[i]]$h = numeric()
      }
    }
  }
  
  realization = function(n_plus,gp_fit,XY,dimX,dimY,rho,heterogeneity.scale,mu_r,sigma_r,mu_h,sigma_h,
                         sp.cov.locs = NULL, sp.cov.vals = NULL, sp.cov.scale = NULL){
    if(n_plus>0){
      
      dat=data.frame(X=numeric(),Y=numeric())
      
      gp_mean = GP_sample(gp_fit,heterogeneity.scale)
      gp_mean = gp_mean - mean(gp_mean)
      
      n_cand = 50
      while(nrow(dat)<n_plus){
        # candidate points within half meter from boundary because sp.cov.locs are at .5 and dim-.5
        xy_cand = matrix(c(runif(n_cand,0.5,dimX-.5),runif(n_cand,0.5,dimY-.5)),nrow=n_cand)
        pred = GP_predict(gp_fit,XY,rho,xy_cand,gp_mean)
        if(!is.null(sp.cov.vals)){
          XB = matrix(0,nrow=length(pred),ncol=1)
          for(k in 1:length(sp.cov.vals)){
            XB = XB + sp.cov.scale[k]*pracma::interp2(x=sp.cov.locs[[k]]$x,y=sp.cov.locs[[k]]$y,Z=sp.cov.vals[[k]],xp=xy_cand[,1],yp=xy_cand[,2])
          }
          accept = runif(n_cand) < plogis(pred + as.numeric(XB),0,logis.scale)
        } else{
          accept = runif(n_cand) < plogis(pred,0,logis.scale)
          
        }
        dat = rbind(dat,data.frame(X=xy_cand[accept,1],Y=xy_cand[accept,2]))
      }
      # randomly sample back down to n_plus shrubs
      dat = dat[sample(seq_len(nrow(dat)),n_plus),]
      # give random radius and height to each shrub
      dat$r = truncdist::rtrunc(n_plus,'norm',a=0,b=Inf,mu_r,sigma_r)
      if(!is.null(mu_h)){
        dat$h = truncdist::rtrunc(n_plus,'norm',a=0,b=Inf,mu_h,sigma_h)
      }
      
    } else{
      dat=data.frame(X=numeric(),Y=numeric(),r=numeric())
      if(!is.null(mu_h))
        dat$h = numeric()
    }
    return(dat)
  }
  
  # if(parallel){
  #   cl = parallel::makeCluster(min(reps,parallel::detectCores()))
  #   doParallel::registerDoParallel(cl)
  #   dat = foreach::foreach(i=1:reps) %dopar% realization(n_plus[i],gp_fit,XY,dimX,dimY,rho,heterogeneity.scale,mu_r,sigma_r,mu_h,sigma_h,sp.cov.locs,sp.cov.vals,sp.cov.scale)
  #   parallel::stopCluster(cl)
  # } else{
  #   dat = vector(mode='list',length=reps)
  #   for(i in 1:reps){
  #     dat[[i]] = realization(n_plus[i],gp_fit,XY,dimX,dimY,rho,heterogeneity.scale,mu_r,sigma_r,mu_h,sigma_h,sp.cov.locs,sp.cov.vals,sp.cov.scale)
  #   }
  # }

  # return parameters as they may be needed to calculate prior
  ret = list(dat=dat,dimX=dimX,dimY=dimY,reps=reps,
             theta=theta,
             heterogeneity.scale = heterogeneity.scale, 
             sp.cov.locs = sp.cov.locs, sp.cov.vals = sp.cov.vals, sp.cov.scale = sp.cov.scale, 
             GP.init.size=GP.init.size, seed = seed, logis.scale=logis.scale, parallel=parallel,
             gamma=gamma)
  class(ret) = c('list','fuelsgen')
  return(ret)
}

# csv.filenames: vector of filenames for repeat observations
load_data = function(csv.filenames, dimX, dimY, heterogeneity.scale = 1, 
                    sp.cov.locs = NULL, sp.cov.vals = NULL, sp.cov.scale = NULL, 
                    GP.init.size=1000, logis.scale=.217622)
{
  reps = length(csv.filenames)
  dat = list(reps)
  for(i in 1:reps){
    # cols: X,Y,r,h
    dat[[i]] = data.frame(read.csv(csv.filenames[i]))
  }
  return(list(dat=dat,dimX=dimX,dimY=dimY,reps=reps,
              heterogeneity.scale = heterogeneity.scale, 
              sp.cov.locs = sp.cov.locs, sp.cov.vals = sp.cov.vals, sp.cov.scale = sp.cov.scale, 
              GP.init.size=GP.init.size, logis.scale=logis.scale))
}