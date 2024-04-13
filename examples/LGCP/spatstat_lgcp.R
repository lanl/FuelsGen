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
nsim = 1000
tmp = rLGCP(mu=-3,
            var=2,
            scale=1,model='gauss',
            win = owin(c(0,20),c(0,20)),nsim=nsim)
plot(tmp[1:25],main='',arrange = T,mar.panel = c(1,1,1,1))
n.sim = numeric(nsim)
for(i in 1:length(tmp)){
  n.sim[i] = tmp[[i]]$n
}

x = 1:200
plot(dpois(x,median(n.sim)),type='l',lty=2,ylab='density',xlab='N')
lines(density(n.sim,from=0),lwd=1.5)
legend('topright',legend=c('n sim','pois'),lty=c(1,2),lwd=c(1.1,1))
################################################################################

dimX = 20
dimY = 20
set.seed(2)
# increasing intensity.scale from 1 to 2 generates way more points
# rho should be affecting it more
fuels = fuelsgen::gen_fuels(dimX,dimX,
                            .1,.75,.2,heterogeneity = 5,reps=10,
                            intensity.scale = .75,seed=1)
plot(fuels)

# n.obs = nrow(fuels$dat[[1]])

ppp.dat = ppp(x=fuels$dat[[1]]$X,
              y=fuels$dat[[1]]$Y,
              c(0,dimX),c(0,dimY))
quad = expand.grid(X=seq(0,dimX,.2),Y=seq(0,dimY,.2))
ux = sort(unique(quad$X)) 
uy = sort(unique(quad$Y)) 
nx = length(ux) 
ny = length(uy) 
col.ref = match(quad$X, ux) 
row.ref = match(quad$Y, uy) 
all.vec = rep(NA, max(row.ref)*max(col.ref)) 
vec.ref = (col.ref- 1)*max(row.ref) + row.ref 
all.vec[vec.ref] = 1 

quads = ppp(quad$X, quad$Y, window = owin(c(0,dimX),c(0,dimY))) 
Q = quadscheme(ppp.dat,dummy = quads, method = "grid", ntile = c(nx, ny), npix = c(nx, ny))

G = hyperframe(shrubs = ppp.dat)
X.des = polynom(quad$X,quad$Y,2)
int.list = list() 
for (i in 1:dim(X.des)[2]){ 
  all.vec = rep(NA, max(row.ref)*max(col.ref)) 
  vec.ref = (col.ref- 1)*max(row.ref) + row.ref 
  all.vec[vec.ref] = X.des[,i] 
  int.list[[i]] = im(matrix(all.vec, max(row.ref), max(col.ref), dimnames = list(uy, ux)), xcol = ux, yrow = uy) 
} 
names(int.list) = paste("V", 1:dim(X.des)[2], sep = "") 
int.form = as.formula(paste("~", paste(names(int.list), collapse = "+"))) 
mod = ppm(Q, trend = as.formula(int.form), covariates = int.list)
summary(mod)

# ppm does a great job at learning the kind of intensity we want (number of points)
# however, we need a stochastic part of the model (cox process i.e. doubly stochastic poisson process)
# because we want realizations of an intensity function that change w.r.t. space
mod = ppm(ppp.dat, interaction = Poisson())
summary(mod)

nsims = 1000
sims = simulate.ppm(mod,nsims)
# plot(sims)
n.sim.ppm = numeric(nsims)
for(i in 1:nsims){
  n.sim.ppm[i] = sims[[i]]$n
}
plot(sims[1:12]) # this is not what we want because we want more flexibility in the intensity function i.e. clumpyness, but not always in the same spot
plot(density(n.sim.ppm,from=0,to=200));abline(v=n.obs,col='red')

# what about a log-gaussian cox process?
nsims = 1000
mod.lgcpK = lgcp.estK(ppp.dat,covmodel = list(model='gauss'))
mod.lgcpK$modelpar
tmp = rLGCP(mu=mod.lgcpK$modelpar[3],
            var=mod.lgcpK$modelpar[1],
            scale=mod.lgcpK$modelpar[2],model='gauss',
            win = owin(c(0,dimX),c(0,dimY)),nsim=nsims)
n.sim.lgcp = numeric(nsims)
for(i in 1:nsims){n.sim.lgcp[i] = tmp[[i]]$n}
lines(density(n.sim.lgcp,from=0,to=200),col='cornflowerblue')#;abline(v=n.obs)
summary(n.sim.lgcp)
plot(tmp[1:12])
# this one produces realizations with ~ between 5 and 620 shrubs, way too variable

# mod.lgcp = lgcp.estpcf(ppp.dat,covmodel = list(model='gauss'))
# mod.lgcp
# tmp = rLGCP(mu=mod.lgcp$modelpar[3],
#             var=mod.lgcp$modelpar[1],
#             scale=mod.lgcp$modelpar[2],model = 'gauss',
#             win = owin(c(0,20),c(0,20)),nsim=nsims)
# plot(tmp[1:12])
# for(i in 1:nsims){n.sim[i] = tmp[[i]]$n}
# hist(n.sim);abline(v=n.obs)
# summary(n.sim) # again, far too variable



# do a fuelsgen model
prior = get_prior_info(fuels,get_mets(fuels),est_cov_obs = T,GP.init.size = 64)
prior$theta_est
prop.sigma = .1^2*diag(prior$theta_est)
source('~/Desktop/R tools/Metro_Hastings_Stochastic.R')
profvis::profvis({
  # issues: proposal variance for rho and I_var are very high - low accenpt
  # correlation of intensity var and hetero are shockickly low
calib = fuelsgen:::mcmc_MH_adaptive(get_mets(fuels),fuels,prior,prop.sigma = prop.sigma,
                                    adapt.par = c(100,20,.5,.5),
                                    n.samples=5000,n.burn=0,gen_reps=25,avg='mets',
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
fuels2 = fuelsgen::gen_fuels(dimX,dimX,1,.1,reps=1000,heterogeneity = 3,intensity.mean = 100, intensity.scale = .1)
n.obs2 = numeric(1000)
for(i in 1:1000){n.obs2[i] = nrow(fuels2$dat[[i]])}
# plot(density(rpois(1000,calib$trace[,4]*20*20),from=0,to=500),col='darkorange',ylim=c(0,.05),main='number of simulated points - lcgp vs fuelsgen',xlab='N');abline(v=n.obs,col='darkgreen')
plot(density(n.obs2,from=0,to=200),col='darkorange',ylim=c(0,.05),main='number of simulated points - lcgp vs fuelsgen',xlab='N');abline(v=n.obs,col='darkgreen')
lines(density(n.sim.lgcp,from=0,to=200),col='cornflowerblue')
lines(density(n.sim.ppm,from=0,to=200),col='maroon')
legend('topright',legend=c('calibrated fuelsgen','LGCP','Poisson PPM','observed n'),col=c('darkorange','cornflowerblue','maroon','darkgreen'),lty=1)
summary(rpois(1000,calib$trace[,4]*20*20))
qpois(.025,min(calib$trace[,4])*20*20)
qpois(.975,max(calib$trace[,4])*20*20)
