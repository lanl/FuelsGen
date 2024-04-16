library(spatstat)
library(fuelsgen)

# PROBLEM
# to get the uncertainty we want in n out of a lgcp
# we must reduce the variance sufficiently, which inhibits 
# heterogeneity even for large scale parameters. To get high heterogeneity,
# we must increase the variance of the GP, which makes the integral of the
# intensity function highly variable, which makes outcomes of n highly variable.
# cranking up the variance, which is needed for high heterogeneity, tends to also
# cause too much overlap in shrubs because the regions of high probability are clumpy
nsim = 100
tmp = rLGCP(mu=0,
            var = .1,
            scale=6,model='gauss',
            win = owin(c(0,20),c(0,20)),nsim=nsim)
plot(tmp[1:12])
n.sim = numeric(nsim)
for(i in 1:length(tmp)){
  n.sim[i] = tmp[[i]]$n
}
hist(n.sim)
summary(n.sim)

mean(n.sim)
var(n.sim)

################################################################################
# generate data from the model

dimX = 20
dimY = 20
reps = 1000
set.seed(2)
# Question: can we improve inference on rho by adding s2_rho to the list of calibration parameters
fuels = fuelsgen::gen_fuels(dimX,dimX,
                            density = .1,radius = .75,sd_radius = .2,
                            heterogeneity = 3,reps=reps, heterogeneity.scale = 1,
                            seed=1, GP.init.size = 32)
n.fg = numeric(reps)
for(i in 1:reps){n.fg[i]=nrow(fuels$dat[[i]])}
hist(n.fg)
mean(n.fg); var(n.fg) # n is drawn from a poisson, so we get poisson uncertainty
################################################################################
# calibration

# do a fuelsgen model
prior = get_prior_info(fuels,get_mets(fuels),est_cov_obs = T,GP.init.size = 64)
prior$theta_est
prop.sigma = .1^2*diag(prior$theta_est)
source('~/Desktop/R tools/Metro_Hastings_Stochastic.R')
profvis::profvis({
  # issues: proposal variance for rho and I_var are very high - low accenpt
  # correlation of intensity var and hetero are shockickly low
  # this is a bit slower than LGCP branch
  calib = fuelsgen:::mcmc_MH_adaptive(get_mets(fuels),fuels,prior,prop.sigma = prop.sigma,
                                      adapt.par = c(100,20,.5,.5),
                                      n.samples=2500,n.burn=0,gen_reps=25,avg='mets',
                                      GP.init.size = 32,
                                      mets_parallel=T,verbose=T)
})
plot(calib)
plot(calib,type='density')

# simulated number of trees

# I think what this is telling us is that we desire the variability in n that comes from a singly stochastic
# poisson process, but we desire the variability in simulations that is only possibly with a doubly stochastic
# poisson process. I'm actually a bit surprised that the uncertainty in our calibrated model is so spot on
# with the PPM uncertainty. I wonder what would happen if we fit the PPM in a Bayesian way. perhaps not much as
# this model is so simple. The ppm model here does estimate a Standard error on lamdba (probably just asymptotic N) and must use that for simulation.

# generate from the model