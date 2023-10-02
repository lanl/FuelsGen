# Given a matrix of xy locations and correlation parameters rho,sd
# build the covariance matrix and return the cholesky and inverse
# for use in drawing samples and making predictions.
GP_fit = function(xy,rho)
{
  n = nrow(xy)
  K = cov_mat(rho,xy) + diag(sqrt(.Machine$double.eps),n)
  # can we make chol or inv faster using kronecker since K is evaluated on a grid?
  cholK = chol(K)
  return(list(cholK=cholK,invK=chol2inv(cholK)))
}

# Given the cholesky of the covariance matrix draw a random sample
# from the mean zero with stderr 'sd'
GP_sample = function(gp_fit,sd)
{
  sd^2*drop(rnorm(nrow(gp_fit$cholK)) %*% gp_fit$cholK)
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
gen_data = function(theta, dimX, dimY, heterogeneity.scale = 1, 
                    sp.cov.locs = NULL, sp.cov.vals = NULL, sp.cov.scale = NULL, 
                    reps = 1, GP.init.size=1000, seed = NULL, logis.scale=.217622, parallel=F)
{
  if(!is.null(seed))
    set.seed(seed)
  lambda = theta[1]
  rho = theta[2]
  mu_r = theta[3]
  sigma_r = theta[4]
  if(length(theta)>4){
    mu_h = theta[5]
    sigma_h = theta[6]
  } else{
    mu_h = NULL
  }
  
  grid.len = as.integer(sqrt(GP.init.size))

  X = seq(0,dimX,length.out=grid.len)
  Y = seq(0,dimY,length.out=grid.len)
  XY = expand.grid(X,Y)
  
  gp_fit = GP_fit(XY,rho)
  n_plus = rpois(reps, lambda * dimX * dimY)
  
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
  
  if(parallel){
    cl = parallel::makeCluster(min(reps,parallel::detectCores()))
    doParallel::registerDoParallel(cl)
    dat = foreach::foreach(i=1:reps) %dopar% realization(n_plus[i],gp_fit,XY,dimX,dimY,rho,heterogeneity.scale,mu_r,sigma_r,mu_h,sigma_h,sp.cov.locs,sp.cov.vals,sp.cov.scale)
    parallel::stopCluster(cl)
  } else{
    dat = vector(mode='list',length=reps)
    for(i in 1:reps){
      dat[[i]] = realization(n_plus[i],gp_fit,XY,dimX,dimY,rho,heterogeneity.scale,mu_r,sigma_r,mu_h,sigma_h,sp.cov.locs,sp.cov.vals,sp.cov.scale)
    }
  }

  # return parameters as they may be needed to calculate prior
  return(list(dat=dat,dimX=dimX,dimY=dimY,reps=reps,
              theta=theta,
              heterogeneity.scale = heterogeneity.scale, 
              sp.cov.locs = sp.cov.locs, sp.cov.vals = sp.cov.vals, sp.cov.scale = sp.cov.scale, 
              GP.init.size=GP.init.size, seed = seed, logis.scale=logis.scale, parallel=parallel))
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