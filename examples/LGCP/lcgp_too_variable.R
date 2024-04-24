library(fuelsgen)
dimX = 20
dimY = 20
set.seed(2)
reps = 25
# Question: can we improve inference on rho by adding s2_rho to the list of calibration parameters
fuels = gen_fuels(dimX,dimX,
                  density = .1,radius = 1,sd_radius = .1,
                  heterogeneity = 4,reps=reps, heterogeneity.scale = 1,
                  seed=1, GP.init.size = 32,I.transform='logistic')
plot(fuels)
# fuelsplot = fuels
# fuelsplot$dat = fuelsplot$dat[1:25]
# fuelsplot$reps = 25
# plot(fuelsplot)
n.obs = numeric(reps)
for(i in 1:reps){n.obs[i] = nrow(fuels$dat[[i]])}
summary(n.obs)
mean(n.obs);var(n.obs)
# fit the model with a lgcp
W = spatstat.geom::owin(c(0,dimX),c(0,dimY))
ppp.list = vector(mode='list',length=reps)
for(i in 1:reps){
  ppp.list[[i]] = spatstat.geom::ppp(x = fuels$dat[[i]]$X,y = fuels$dat[[i]]$Y,window = W)
}

#Fit the trend part of the model to your replicated point pattern data using mppm, for example m <- mppm(Y ~ X, data=G)
G = spatstat.geom::hyperframe(shrubs = ppp.list)
trend = spatstat.model::mppm(shrubs ~ 1, G, interaction = spatstat.model::Poisson())
# trend_R = spatstat.model::mppm(shrubs ~ 1, G, interaction = spatstat.model::Strauss(1))
# plot(simulate(trend_R))

# Extract the fitted intensities for each point pattern using predict.mppm
lambda = predict(trend, type='lambda')

# For each point pattern, using the corresponding intensity obtained from the model, 
# compute the inhomogeneous K function using Kinhom (with argument ratio=TRUE)
Ki = lapply(1:reps, function(i) spatstat.explore::Kinhom(ppp.list[[i]],lambda$cif[[i]],ratio = T))

# Combine the K functions using pool
Kpool = do.call(spatstat.explore::pool,Ki)

# Estimate the cluster parameters of the LGCP by applying lgcp.estK to the pooled K function.
cox.mod = spatstat.model::lgcp.estK(Kpool, covmodel =  list('gauss'))
cox.mod$modelpar

# generate a number of realizations
cox.sim = spatstat.random::rLGCP(model='gauss',
                mu = coef(trend),var=cox.mod$modelpar[1],scale=cox.mod$modelpar[2],win = W,
                nsim = 500)
plot(cox.sim[1:25])
n.sim = numeric(500)
for(i in 1:500){n.sim[i] = cox.sim[[i]]$n}
hist(n.sim)
summary(n.sim)
mean(n.sim);var(n.sim)

par(mfrow=c(1,1),mar=c(4,4,4,4))
hist(n.sim,ylim=c(0,.08),freq = F,col=adjustcolor('darkgrey',alpha=.5),
     main='',ylab='density',xlab='n')
hist(n.obs,add=T,col=adjustcolor('darkgreen',alpha=.5),freq=F)
legend('topright',c('fitted LGCP realizations','training realizations'),
       col=c(adjustcolor('darkgrey',alpha=.5),adjustcolor('darkgreen',alpha=.5)),lty=1,lwd=8)

# calibrate with fuelsgen
library(spatstat)

# first use only the metrics defined in the paper
info = get_mets_info()
y_obs = get_mets(fuels,info)
prior = get_prior_info(fuels,y_obs,T)
prop = .1^2*diag(prior$theta_est)
prop[1,1] = .1*prior$theta_est[1]
fg.calib = mcmc_MH_adaptive(y_obs,fuels,prior,n.samples = 10000, n.burn = 0,I.transform = 'logistic',avg = 'mets',
                            prop.sigma = prop, adapt.par = c(100,50,.5,75), gen_reps = 22) # 12 core system so do 11*2
plot(fg.calib,burn=0) # using the set of summary statistics defined in the paper - 
plot(fg.calib,type = 'density',burn = 5000, truth = c(4,1,1/.1^2,.1))

# now use only the summary statistics defined in Moller book (much much faster)
info2 = get_mets_info(F,F,F,F,F,F,F,F,F,0,T,T,F)
y_obs2 = get_mets(fuels,info2)
prior2 = get_prior_info(fuels,y_obs2,T)
prop2 = .1^2*diag(prior2$theta_est)
prop2[1,1] = .1*prior2$theta_est[1]
fg.calib2 = mcmc_MH_adaptive(y_obs2,fuels,prior2,n.samples = 10000, n.burn = 0,I.transform = 'logistic',avg = 'mets',
                            prop.sigma = prop2, adapt.par = c(100,50,.5,75), gen_reps = 22) # 12 core system so do 11*2
plot(fg.calib2,burn=0) # using only summ stats from book - works decently well for constraining rho
plot(fg.calib2,type = 'density',burn = 6000, truth = c(4,1,1/.1^2,.1))

note = "comparison of calibration metrics. fg.calib contains a calibration to the fuels object using the metrics defined in the paper. fg.calib2 uses only the inhomogeneous K function and the pairwise correlation function evalutated at 3 points each"
save(note,fuels,fg.calib,fg.calib2,file='examples/LGCP/compare_calib_metrics.RData')

# possible n in simulations
burn = 5000
n = numeric(10000-burn)
for(i in 1:(10000-burn)){
  n[i] = rpois(1,dimX*dimY*fg.calib2$trace[i,4])
}
# possible n from model which generated training data
fuels_many = gen_fuels(dimX,dimX,
                       density = .1,radius = 1,sd_radius = .1,
                       heterogeneity = 4,reps=1000, heterogeneity.scale = 1,
                       seed=1, GP.init.size = 32,I.transform='logistic')
n.obs.poss = numeric(1000)
for(i in 1:1000){n.obs.poss[i]=nrow(fuels_many$dat[[i]])}
par(mfrow=c(1,1))
plot(density(n),ylim=c(0,.08),xlim=c(0,max(n.sim)),col='maroon')
lines(density(n.obs),col='darkgreen')
#lines(density(n.obs.poss),col='lightblue1')
lines(density(n.sim),col='darkorange')
lines(density(n.obs.poss),lty=2)
legend('topright',c('fitted FuelsGen Model','training realizations','fitted lgcp realizations'),
       col=c('maroon','darkgreen','darkorange'),lty=1,lwd=8)
  