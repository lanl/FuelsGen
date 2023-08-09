
# This file shows how to calibrate the generative model to data generated from the model.

# generate data with known parameters
lambda.true = .25
rho.true = 3
mu.true = 1
sigma.true = .1
fuel = gen_fuels(dimX = 30, dimY = 30,
                 density = lambda.true,        # shrubs per unit area
                 heterogeneity = rho.true,     # level of heterogeneity in shrub placement
                 heterogeneity.scale = 1,      # scale of the mean-zero GP realization
                 radius = mu.true, sd_radius = sigma.true,  # normal distribution parameters for shrub radius
                 height = NULL, sd_height = 0, # normal distribution parameters for shrub height
                 reps=10,                      # number of random maps to generate
                 GP.init.size=100,             # How densely to sample GP (n^3 scaling, <=1000 is pretty fast)
                 seed=10,                      # random seed for reproducibility
                 parallel=F)                   # Parallel option will be faster for expensive generation only (parallel overhead), rng seed doesn't work with parallel

# compute spatial metrics on observed fuel maps
y_obs = get_mets(fuel, parallel = T, make.cluster = T)

# precompute important information for priors, adding additional simulation data if requested.
# list(prior='gamma',params=c(1,.33)) defines a Gamma(1,1/3) prior for the lengthscale
# in the paper we also explore a Uniform(0,10) defined as list(prior='unif',params=c(0,10))
prior = get_prior_info(fuel,y_obs,est_cov_obs = F,est_cov_samples = 8, est_cov_reps = 8, 
                                  est_rho_prior = list(prior='gamma',params=c(1,.33)),gen_parallel = F,
                                  mets_parallel = T, make.cluster=T)

# calibrate the model to y_obs
mcmc.out = mcmc(y_obs,fuel,prior,'./',n.samples=10000,n.burn=1000,gen_reps=25,avg='mets',
                           gen_parallel=F,mets_parallel=T,verbose=10,load.theta=F,update.every=1)
