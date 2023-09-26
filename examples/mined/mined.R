library(fuelsgen)
mined_logf_mu_only = function(mu,data=dat){
  
  theta = c(data$rho,mu,data$sigma,data$lambda)
  
  fuelsgen:::llh(data$y_obs, theta, data$dimX, data$dimY, data$prior$ldetS, data$prior$Sinv, data$gen_reps, data$avg, data$gen_parallel, data$mets_parallel) + 
      fuelsgen:::lprior(theta,data$prior$prior_params)
}
mined_logf_mu_sigma = function(theta,data=dat){
  
  theta = c(data$rho,theta,data$lambda)
  
  fuelsgen:::llh(data$y_obs, theta, data$dimX, data$dimY, data$prior$ldetS, data$prior$Sinv, data$gen_reps, data$avg, data$gen_parallel, data$mets_parallel) + 
    fuelsgen:::lprior(theta,data$prior$prior_params)
}
mined_logf_mu_lambda = function(theta,data=dat){
  
  theta = c(data$rho,theta[1],data$sigma,theta[2])
  
  fuelsgen:::llh(data$y_obs, theta, data$dimX, data$dimY, data$prior$ldetS, data$prior$Sinv, data$gen_reps, data$avg, data$gen_parallel, data$mets_parallel) + 
    fuelsgen:::lprior(theta,data$prior$prior_params)
}
mined_logf_sigma_lambda = function(theta,data=dat){
  
  theta = c(data$rho,data$mu,theta)
  
  fuelsgen:::llh(data$y_obs, theta, data$dimX, data$dimY, data$prior$ldetS, data$prior$Sinv, data$gen_reps, data$avg, data$gen_parallel, data$mets_parallel) + 
    fuelsgen:::lprior(theta,data$prior$prior_params)
}
mined_logf_rho_mu = function(theta,data=dat){
  theta = c(theta,data$sigma,data$lambda)
  if(any(theta<0)){
    return(-Inf)
  } else{
    return(fuelsgen:::llh(data$y_obs, theta, data$dimX, data$dimY, data$prior$ldetS, data$prior$Sinv, data$gen_reps, data$avg, data$gen_parallel, data$mets_parallel) + 
             fuelsgen:::lprior(theta,data$prior$prior_params))
  }
}

lambda = .25
rho = 3
mu = .5
sigma = .1
fuel = gen_fuels(dimX = 30, dimY = 30,
                 density = lambda,        # shrubs per unit area
                 heterogeneity = rho,     # level of heterogeneity in shrub placement
                 heterogeneity.scale = 1,      # scale of the mean-zero GP realization
                 radius = mu, sd_radius = sigma,  # normal distribution parameters for shrub radius
                 height = NULL, sd_height = 0, # normal distribution parameters for shrub height
                 reps=10,                      # number of random maps to generate
                 GP.init.size=100,             # How densely to sample GP (n^3 scaling, <=1000 is pretty fast)
                 seed=10,                      # random seed for reproducibility
                 parallel=F)                   # Parallel option will be faster for expensive generation only (parallel overhead), rng seed doesn't work with parallel

# compute spatial metrics on observed fuel maps
y_obs = get_mets(fuel, parallel = T, make.cluster = T)
prior = get_prior_info(fuel,y_obs,est_cov_obs = F,est_cov_samples = 25, est_cov_reps = 25, 
                       est_rho_prior = list(prior='gamma',params=c(1,.33)),gen_parallel = F,
                       mets_parallel = T, make.cluster=T)
dat = list(y_obs=y_obs,rho=rho,mu=mu,sigma=sigma,lambda=lambda,dimX=fuel$dimX,dimY=fuel$dimY,prior=prior,gen_reps=25,avg='llh',gen_parallel=F,mets_parallel=T)
save(rho,mu,sigma,lambda,fuel,y_obs,prior,dat,file='test_data.RData')

# Calibrate mu only with all other parameters fixed at the truth
init = matrix(seq(.3,.6,.01))
cl = parallel::makeCluster(8)
doParallel::registerDoParallel(cl)
samples = mined::mined(init,mined_logf_mu_only)
hist(samples$points,xlim=c(.3,.6))
plot(samples$cand,samples$candlf)
save(init,samples,file='mu_only_samples.RData')

# Calibrate mu,sigma jointly
init = as.matrix(expand.grid(seq(.45,.55,.025),
                   seq(.05,.15,.025)))
samples = mined::mined(init,mined_logf_mu_sigma)

# Calibrate mu,lambda jointly
init = as.matrix(expand.grid(seq(.45,.55,.05),
                             seq(.25,.3,.025)))
samples = mined::mined(init,mined_logf_mu_lambda)

# Calibrate sigma,lambda jointly
init = as.matrix(expand.grid(seq(.05,.15,.05),
                             seq(.2,.3,.05)))
samples = mined::mined(init,mined_logf_sigma_lambda)

# Calibrate mu,rho jointly
init = as.matrix(expand.grid(seq(.45,.55,.05),
                             seq(2.5,3.5,.5)))
samples = mined::mined(init,mined_logf_rho_mu)
