# Given a matrix of xy locations and correlation parameters rho,sd
# build the covariance matrix and return the cholesky and inverse
# for use in drawing samples and making predictions.
GP_fit = function(xy,rho)
{
  n = nrow(xy)
  K = cov_mat(rho,xy) + diag(sqrt(.Machine$double.eps),n)
  # can we make chol or inv faster using kronecker since K is evaluated on a grid?
  cholK = chol(K)
  invK = chol2inv(cholK)
  return(list(cholK=cholK,invK=invK))
}

# Given the cholesky of the covariance matrix draw a random sample
# from the mean zero with stderr 'sd'
GP_sample = function(gp_fit,mu=0,sd=1)
{
  mu + sd^2*drop(rnorm(nrow(gp_fit$cholK)) %*% gp_fit$cholK)
}

# With covariance K defined on points xy predict from the GP
# at xy.pred using training points y.train
GP_predict = function(gp_fit,xy,rho,xy.pred,y.train)
{
  # cross-covariance of training locations and prediction locations
  K.xy.pred = cov_mat(rho,xy.pred,xy)
  # mean prediction
  K.xy.pred%*%gp_fit$invK%*%y.train
}

# compute squared exponential covariance matrix at points x,y
cov_mat = function(rho,x,y=NULL)
{
  if(is.null(y)){
    d = as.matrix(distances::distances(x))
  } else{
    d = sqrt(plgp::distance(x,y))
  }
  exp(-d^2 / (2*rho^2))
}

# Generate 'reps' fuel maps with parameters theta over a domain of size [0,dimX]x[0,dimY]
gen_data = function(theta, dimX, dimY, 
                    sp.cov.locs, sp.cov.vals, sp.cov.scale, 
                    reps, GP.init.size, seed, sample.method='direct')
{

  if(!is.null(seed))
    set.seed(seed)
  lambda = theta[1]
  rho = theta[2]
  mu_r = theta[3]
  sigma_r = theta[4]
  # intensity.mean = theta[5]
  intensity.var = theta[5]
  if(length(theta)>5){
    mu_h = theta[6]
    sigma_h = theta[7]
  } else{
    mu_h = NULL
  }

  W = owin(c(0,dimX),c(0,dimY), mask=matrix(TRUE, GP.init.size,GP.init.size))
  # for this formulation we don't care about the intensity mean because the intensity function
  # is going to get scaled so that it integrates to n
  gamma = exp(rGRFgauss(W = W,mu = 0,var = intensity.var,scale = rho,nsim = reps))

  if(sample.method=='direct'){
    #### sample pixels directly ####
    mu = sapply(1:reps,function(i) integral(gamma[[i]]))
    nn = rpois(reps,lambda*dimX*dimY)
    if(!all(is.finite(nn)))
      stop(paste("Unable to generate Poisson process with a mean of",
                 mu, "points"))
    Eta = lapply(1:reps, function(i) gamma[[i]] * nn[i] / mu[i])
    # mu = sapply(1:reps,function(i) integral(Eta[[i]]))
    dat <- vector(mode="list", length=reps)

    for(isim in seq_len(reps)){
      if(nn[isim]>0){
        dx <- Eta[[isim]]$xstep/2
        dy <- Eta[[isim]]$ystep/2
        df <- as.data.frame(Eta[[isim]])
        npix <- nrow(df)
        lpix <- df$value
        ni <- nn[isim]
        ii <- sample.int(npix, size=ni, replace=TRUE, prob=lpix)
        xx <- df$x[ii] + runif(ni, -dx, dx)
        yy <- df$y[ii] + runif(ni, -dy, dy)
        dat[[isim]] <- data.frame(X=xx,Y=yy)
        dat[[isim]]$r = truncdist::rtrunc(nn[isim],'norm',a=0,b=Inf,mu_r,sigma_r)
        if(!is.null(mu_h)){
          dat[[isim]]$h = truncdist::rtrunc(nn[isim],'norm',a=0,b=Inf,mu_h,sigma_h)
        }
      } else{
        dat[[isim]] = data.frame(X=numeric(),Y=numeric(),r=numeric())
        if(!is.null(mu_h)){
          dat[[isim]]$h = numeric()
        }
      }
    }
  } else{
    #### THINNING ####
    # maximum of the intensity is used for determining probability of point placement
    summ <- summary(gamma)
    gamma_max <- summ$max + 0.05 * diff(summ$range)
    # gamma_max = max(gamma)
    dat <- runifpoispp(gamma_max, W, nsim=1, drop=FALSE)
    for(isim in seq_len(nsim)) {
      X <- dat[[isim]]
      if(X$n > 0) {
        prob <- lambda[X, drop=FALSE]/gamma_max
        u <- runif(X$n)
        retain <- is.finite(prob) & (u <= prob)
        dat[[isim]] <- X[retain]
      }
    }
    
    
    cell_probs = gamma / gamma_max
    # generate uniform candidate points based on gamma_max
    n_cand = rpois(1,gamma_max)#*dimX*dimY)
    if(n_cand>0){
      # keep = logical(n_cand)
      xy_cand = matrix(c(runif(n_cand,0,dimX),
                         runif(n_cand,0,dimY)),
                       nrow=n_cand)
      # which cell are the candidates in
      xcell = floor(xy_cand[,1] / dx)
      ycell = floor(xy_cand[,2] / dy)
      cellnum = ycell*length(X) + xcell + 1
      # keep each point in cell i with p = gamma[i] / gamma_max
      u = runif(n_cand)
      keep = u < cell_probs[cellnum]
      # for(i in 1:nrow(xy_cand)){
      #   keep[i] = runif(1)<cell_probs[cellnum[i]]
      # }
      
      # loop over cells and 
      # for(i in 1:nrow(XY)){
      #   # find cell boundaries
      #   xb = XY[i,1] + c(-1,1)*dx
      #   yb = XY[i,2] + c(-1,1)*dy
      #   # find candidates within this cell
      #   cell_cand = which(between(xy_cand[,1],xb[1],xb[2]) & between(xy_cand[,2],yb[1],yb[2]))
      #   keep[cell_cand] = runif(length(cell_cand)) < cell_probs[i]
      # }
      dat = data.frame(X=xy_cand[keep,1],Y=xy_cand[keep,2])
    }
  }

  # give random radius and height to each shrub
  
  # else{
  #   dat = data.frame(X=numeric(),Y=numeric(),r=numeric(),h=numeric())
  # }
  # return(dat)
  # }

  # return parameters as they may be needed to calculate prior
  ret = list(dat=dat,dimX=dimX,dimY=dimY,reps=reps,
                theta=theta,
                #intensity.mean = intensity.mean,
                intensity.var = intensity.var, 
                sp.cov.locs = sp.cov.locs, sp.cov.vals = sp.cov.vals, sp.cov.scale = sp.cov.scale, 
                GP.init.size=GP.init.size, seed = seed)
  class(ret) = c('list','fuelsgen')
  return(ret)
}

# csv.filenames: vector of filenames for repeat observations
load_data = function(csv.filenames, dimX, dimY, intensity.var = 1, 
                    sp.cov.locs = NULL, sp.cov.vals = NULL, sp.cov.scale = NULL)
{
  reps = length(csv.filenames)
  dat = list(reps)
  for(i in 1:reps){
    # cols: X,Y,r,h
    dat[[i]] = data.frame(read.csv(csv.filenames[i]))
  }
  return(list(dat=dat,dimX=dimX,dimY=dimY,reps=reps,
              intensity.var = intensity.var, 
              sp.cov.locs = sp.cov.locs, sp.cov.vals = sp.cov.vals, sp.cov.scale = sp.cov.scale, 
              ))
}