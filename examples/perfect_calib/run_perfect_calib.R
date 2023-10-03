library(fuelsgen)
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
                 reps=1,                      # number of random maps to generate
                 GP.init.size=100,             # How densely to sample GP (n^3 scaling, <=1000 is pretty fast)
                 seed=10,                      # random seed for reproducibility
                 parallel=F)                   # Parallel option will be faster for expensive generation only (parallel overhead), rng seed doesn't work with parallel
plot_fuels(fuel)
y_obs = get_mets(fuel, parallel = T, make.cluster = T)
prior = get_prior_info(fuel,y_obs,est_cov_obs = T,est_cov_samples = 25, est_cov_reps = 25,
                       est_rho_prior = list(prior='unif',params=c(0,10)),gen_parallel = F,
                       mets_parallel = T, make.cluster=T)
save.image('examples/perfect_calib/perfect_ex_m_1_diag_relative_prior.RData')

mcmc.out = mcmc(y_obs,fuel,prior,'./examples/perfect_calib/perfect_ex_m_1_diag_relative_mcmc.RData',
                n.samples=100,n.burn=0,gen_reps=25,avg='mets',
                gen_parallel=F,mets_parallel=T,verbose=10,load.theta=F,update.every=1)
save.image('examples/perfect_calib/perfect_ex_m_1_diag_relative_prior.RData')
pairs(mcmc.out$theta)
par(mfrow=c(2,2))
for(i in 1:4){
  hist(mcmc.out$theta[,i])
}
for(i in 1:4){
  plot(mcmc.out$theta[,i],type='l')
}
