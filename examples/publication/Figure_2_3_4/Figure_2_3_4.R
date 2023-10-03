load(here::here('examples/publication/Figure_2_3_4/mcmc_n_1_y_cov_avg_mets.RData'))
theta.n.1.avg.mets = theta
# load(here::here('examples/publication/Figure_2_3_4/mcmc_n_10_avg_mets.RData'))
# theta.n.10.avg.mets = theta
load(here::here('examples/publication/Figure_2_3_4/mcmc_n_25_y_cov_avg_mets.RData'))
theta.n.25.avg.mets = theta


# unif avg mets
# load(here::here('examples/publication/Figure_2_3_4/mcmc_n_1_rho_unif_0-10_avg_mets.RData'))
# theta.n.1.unif.0.10.avg.mets = theta
# load(here::here('examples/publication/Figure_2_3_4/mcmc_n_10_rho_unif_0-10_avg_mets.RData'))
# theta.n.10.unif.0.10.avg.mets = theta
# load(here::here('examples/publication/Figure_2_3_4/mcmc_n_25_rho_unif_0-10_avg_mets.RData'))
# theta.n.25.unif.0.10.avg.mets = theta

# gamma avg mets
load(here::here('examples/publication/Figure_2_3_4/mcmc_n_1_rho_gamma_3-1_avg_mets.RData'))
theta.n.1.gamma.3.1.avg.mets = theta
load(here::here('examples/publication/Figure_2_3_4/mcmc_n_10_rho_gamma_3-1_avg_mets.RData'))
theta.n.10.gamma.3.1.avg.mets = theta
load(here::here('examples/publication/Figure_2_3_4/mcmc_n_25_rho_gamma_3-1_avg_mets.RData'))
theta.n.25.gamma.3.1.avg.mets = theta

truth = c(3,1,.1^2,.25)

# Figure 2
par(usr=c(0,30,0,30))
set.seed(2)
png('examples/publication/Figure_2_3_4/perfect_prior_data_post.png',
    width=8,height=5,units='in',res=1200)
# Prior
n = 5
rho = runif(n,0,10)
mu = truncdist::rtrunc(n,'norm',a=0,b=3,mean=1.5,.5)
var = rgamma(n,1,scale=1e-3)
lambda = rgamma(n,1,1/.25)
prior_realization = vector(mode='list',length=n)
par(mfrow=c(3,5),mar=c(.5,.5,.5,.5))
for(i in 1:n){
  prior_realization[[i]] = gen_fuels(30,30,lambda[i],mu[i],sqrt(var[i]),NULL,NULL,rho[i],1,NULL,NULL,NULL,1,100,1)
  plot(x=c(),y=c(),xlim=c(0,30),ylim=c(0,30),xlab='X',ylab='Y',xaxt='n',yaxt='n',asp=1)
  for(j in 1:nrow(prior_realization[[i]]$dat[[1]])){
    plotrix::draw.circle(prior_realization[[i]]$dat[[1]]$X[j],prior_realization[[i]]$dat[[1]]$Y[j],prior_realization[[i]]$dat[[1]]$r[j],col=adjustcolor('grey',alpha.f=.5))
  }
}

# Observed Data
true_realization = gen_fuels(30,30,.25,1,.1,NULL,NULL,3,1,NULL,NULL,NULL,10,100,10)
for(i in 1:n){
  plot(x=c(),y=c(),xlim=c(0,30),ylim=c(0,30),xlab='X',ylab='Y',xaxt='n',yaxt='n',asp=1)
  for(j in 1:nrow(true_realization$dat[[i]])){
    plotrix::draw.circle(true_realization$dat[[i]]$X[j],true_realization$dat[[i]]$Y[j],true_realization$dat[[i]]$r[j],col=adjustcolor('grey',alpha.f=.5))
  }
}

# Posterior
rho = sample(theta[,1],n)
mu = sample(theta[,2],n)
sigma = sample(theta[,3],n)
lambda = sample(theta[,4],n)
post_realization = vector(mode='list',length=n)
for(i in 1:n){
  post_realization[[i]] = gen_fuels(30,30,lambda[i],mu[i],sigma[i],NULL,NULL,rho[i],1,NULL,NULL,NULL,1,100,i)
  plot(x=c(),y=c(),xlim=c(0,30),ylim=c(0,30),xlab='X',ylab='Y',xaxt='n',yaxt='n',asp=1)
  for(j in 1:nrow(post_realization[[i]]$dat[[1]])){
    plotrix::draw.circle(post_realization[[i]]$dat[[1]]$X[j],post_realization[[i]]$dat[[1]]$Y[j],post_realization[[i]]$dat[[1]]$r[j],col=adjustcolor('grey',alpha.f=.5))
  }
}
dev.off()

# Figure 3
png("examples/publication/Figure_2_3_4/n_10_vs_n_25.png",width=7,height=5,units="in",res=1200)
# contour_pairs(list(theta.n.10.avg.mets,theta.n.25.avg.mets),
#               labels = c('m=10','m=25'),
#               truth,mu.min=.7,mu.max=1.3)
contour_pairs(list(theta.n.25.avg.mets),
              labels = c('m=10','m=25'),
              truth,mu.min=.7,mu.max=1.3)
dev.off()

# Figure 4
png('examples/publication/Figure_2_3_4/n_prior_confidence.png')
labs = c('rho','mu','sigma^2','lambda')
lims = list(c(0,10),c(.8,1.5),c(0,.21),c(.1,.4))
truth = c(3,1,.1,.25)
# prior intervals rho, mu, sigma, lambda
# var.samp = rgamma(1e6,prior$prior_params$prec_shape,1/prior$prior_params$prec_rate)
# quant.prior = cbind(qunif(c(.025,.975),0,10),
#                     truncdist::qtrunc(c(.025,.975),'norm',a=0,b=3,mean=prior$prior_params$mu_mean,sd=prior$prior_params$mu_sd),
#                     quantile(sqrt(var.samp),c(.025,.975)),
#                     qgamma(c(.025,.975),prior$prior_params$lambda_shape,prior$prior_params$lambda_rate))
# lims = list(quant.prior[,1],quant.prior[,2],c(0,.25),quant.prior[,4])
par(mfrow=c(2,2),mar=c(2.8,2,1,2))
for(i in 1:4){
  plot(c(),c(),xlim=c(0,27),ylim=lims[[i]],xlab='',xaxt='n')
  title(xlab="m", line=1.75, cex.lab=1.2)
  
  axis(1,at=c(1,10,25))
  if(i==1){
    legend(inset=c(.1,.05),'topleft',legend=expression(rho),cex=1.5,bty='n')
  } else if(i==2){
    legend(inset=c(.1,.05),'topleft',legend=expression(mu),cex=1.5,bty='n')
  } else if(i==3){
    legend(inset=c(.1,.05),'topleft',legend=expression(sigma),cex=1.5,bty='n')
  } else{
    legend(inset=c(.1,.05),'topleft',legend=expression(lambda),cex=1.5,bty='n')
  }
  abline(h=truth[i],col='red',lty=2)
  # quant.unif = cbind(quantile(theta.n.1.unif.0.10.avg.mets[,i],c(.025,.975)),
  #                    quantile(theta.n.10.unif.0.10.avg.mets[,i],c(.025,.975)),
  #                    quantile(theta.n.25.unif.0.10.avg.mets[,i],c(.025,.975)))
  
  quant.obs = cbind(quantile(theta.n.1.avg.mets[,i],c(.025,.975)),
                    quantile(theta.n.25.avg.mets[,i],c(.025,.975)))
  
  quant.gam = cbind(quantile(theta.n.1.gamma.3.1.avg.mets[,i],c(.025,.975)),
                    quantile(theta.n.10.gamma.3.1.avg.mets[,i],c(.025,.975)),
                    quantile(theta.n.25.gamma.3.1.avg.mets[,i],c(.025,.975)))
  # quant.unif.2.4 = cbind(quantile(theta.n.10.unif.2.4[,i],c(.025,.975)))
  # arrows(x0=c(1,10,25),
  #        y0=quant.unif[1,],
  #        y1=quant.unif[2,],
  #        code=3, angle=90, length=0.05, col="darkgreen", lwd=1)
  # arrows(x0=c(10)+2,
  #        y0=quant.unif.2.4[1],
  #        y1=quant.unif.2.4[2],
  #        code=3, angle=90, length=0.05, col="maroon", lwd=1)
  arrows(x0=c(1,25)-1,
         y0=quant.obs[1,],
         y1=quant.obs[2,],
         code=3, angle=90, length=0.05, col="darkorange", lwd=1)
  arrows(x0=c(1,10,25)+1,
         y0=quant.gam[1,],
         y1=quant.gam[2,],
         code=3, angle=90, length=0.05, col="cornflowerblue", lwd=1)
  # arrows(x0=c(1,10,25)-1,
  #        y0=rep(quant.prior[1,i],3),
  #        y1=rep(quant.prior[2,i],3),
  #        code=3, angle=90, length=0.05, lwd=1)
  if(i==2)
    legend("topright", inset = c(0, 0),                   # Create legend outside of plot
           legend = c("Data","Data+Gam(3,1)","Data+Unif(0,10)"), # "Data+Unif(2,4)"),
           lty=1,
           col = c('darkorange','cornflowerblue','forestgreen'))#,'maroon'))
}
dev.off()
