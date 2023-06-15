#' @title Plot MCMC posteriors
#'
#' @description Pairs plot of posterior density estimates and prior distributions
#' @param theta.list list of matrices containing posterior samples
#' @param truth optional vector of known parameter values
#' @param cols colors for posterior density estimate plots
#' @param labels legend names for posterior density estimates
#' @export
#'
contour_pairs = function(theta.list,truth=NULL,cols=c('cornflowerblue','darkorange','forestgreen'),
                         rho.min=0,rho.max=10,rho.prior='unif',
                         mu.min=0,mu.max=3,mu.prior.mean=1.5,mu.prior.sd=.5,
                         var.min = 0, var.max=.05,var.prior.type='gamma',var.prior.params=c(1,1e-3),
                         lambda.max=.5,lambda.prior.type='gamma',lambda.prior.params=c(1,.25),
                         labels=c('n=10','n=25'))
{
  n.theta = length(theta.list)
  rho.x = seq(0,rho.max,length.out=1000)
  rho.prior = dunif(rho.x,rho.min,rho.max)
  # p(mu) ~ truncNormal
  mu.x = seq(0,mu.max,length.out=1000)
  mu.prior = truncdist::dtrunc(mu.x,'norm',a=0,b=3,mean=mu.prior.mean,sd=mu.prior.sd)
  # p(gamma) ~ Gamma
  var.x = seq(var.min,var.max,length.out=1000)
  var.prior = dgamma(var.x,var.prior.params[1],scale=var.prior.params[2])
  # p(lam) ~ Gamma - semi-informative gamma prior with mean at \hat{\lam}
  lambda.x = seq(0,lambda.max,length.out=1000)
  lambda.prior = dgamma(lambda.x,lambda.prior.params[1],scale=lambda.prior.params[2])
  
  par(mfrow=c(4,4),mar=c(2,2,.5,.5))
  
  # rho prior-posterior
  rho.dens = list(n.theta)
  max = 0
  for(i in 1:n.theta){
    rho.dens[[i]] = density(theta.list[[i]][,1],from=rho.min,to=rho.max)
    M = max(rho.dens[[i]]$y)
    if(M>max)
      max = M
  }
  plot(rho.x,rho.prior,type='l',ylim=c(0,max+.2),lty=2,xlab='',yaxt='n')
  for(i in 1:n.theta){
    lines(rho.dens[[i]],col=cols[i])
  }
  abline(v=truth[1],col='red',lty=2)
  legend(inset=c(.05,.05),'topright',legend=expression(rho),cex=1.5,bty='n')
  plot.new()
  if(!is.null(truth)){
    legend('center',legend=c('prior',labels,'truth'),lty=c(2,rep(1,n.theta),2),col=c('black',cols[1:n.theta],'red'),cex=1.2)
  } else{
    legend('center',legend=c('prior',labels),lty=c(2,rep(1,n.theta)),col=c('black',cols[1:n.theta]),cex=1.2)
  }
  plot.new()
  plot.new()
  
  # mu prior-posterior
  for(i in 1:n.theta){
    if(i==1){
      contour(MASS::kde2d(theta.list[[i]][,1], theta.list[[i]][,2]),nlevels = 10,col=cols[i],xlim=c(rho.min,rho.max),ylim=c(mu.min,mu.max))
    } else{
      contour(MASS::kde2d(theta.list[[i]][,1], theta.list[[i]][,2]),nlevels = 10,add=T,col=cols[i])
    }
  }
  points(truth[1],truth[2],col='red',pch=16)
  
  mu.dens = list(n.theta)
  max = 0
  for(i in 1:n.theta){
    mu.dens[[i]] = density(theta.list[[i]][,2],from=mu.min,to=mu.max)
    M = max(mu.dens[[i]]$y)
    if(M>max)
      max = M
  }
  plot(mu.x,mu.prior,type='l',xlim=c(mu.min,mu.max),ylim=c(0,max),lty=2,xlab='',yaxt='n')
  for(i in 1:n.theta){
    lines(mu.dens[[i]],col=cols[i])
  }
  abline(v=truth[2],col='red',lty=2)
  legend(inset=c(.05,.05),'topright',legend=expression(mu),cex=1.5,bty='n')
  plot.new();plot.new()

  # variance prior-posterior
  for(i in 1:n.theta){
    if(i==1){
      contour(MASS::kde2d(theta.list[[i]][,1], theta.list[[i]][,3]^2),nlevels = 10,col=cols[i],xlim=c(rho.min,rho.max),ylim=c(var.min,var.max))
    } else{
      contour(MASS::kde2d(theta.list[[i]][,1], theta.list[[i]][,3]^2),nlevels = 10,add=T,col=cols[i])
    }
  }
  points(truth[1],truth[3],col='red',pch=16)
  for(i in 1:n.theta){
    if(i==1){
      contour(MASS::kde2d(theta.list[[i]][,2], theta.list[[i]][,3]^2),nlevels = 10,col=cols[i],ylim=c(var.min,var.max),xlim=c(mu.min,mu.max))
    } else{
      contour(MASS::kde2d(theta.list[[i]][,2], theta.list[[i]][,3]^2),nlevels = 10,add=T,col=cols[i])
    }
  }
  points(truth[2],truth[3],col='red',pch=16)

  var.dens = list(n.theta)
  max = 0
  for(i in 1:n.theta){
    var.dens[[i]] = density(theta.list[[i]][,3]^2,from=0,to=var.max)
    M = max(var.dens[[i]]$y)
    if(M>max)
      max = M
  }
  plot(var.x,var.prior,type='l',xlim=c(var.min,var.max),ylim=c(0,max),lty=2,xlab='',yaxt='n')
  for(i in 1:n.theta){
    lines(var.dens[[i]],col=cols[i])
  }
  abline(v=truth[3],col='red',lty=2)
  legend(inset=c(.05,.05),'topright',legend=expression(sigma^2),cex=1.5,bty='n')
  plot.new()
  
  # lambda prior-posterior
  for(i in 1:n.theta){
    if(i==1){
      contour(MASS::kde2d(theta.list[[i]][,1], theta.list[[i]][,4]),nlevels = 10,col=cols[i],xlim=c(rho.min,rho.max),ylim=c(0,lambda.max))
    } else{
      contour(MASS::kde2d(theta.list[[i]][,1], theta.list[[i]][,4]),nlevels = 10,add=T,col=cols[i])
    }
  }
  points(truth[1],truth[4],col='red',pch=16)
  for(i in 1:n.theta){
    if(i==1){
      contour(MASS::kde2d(theta.list[[i]][,2], theta.list[[i]][,4]),nlevels = 10,col=cols[i],ylim=c(0,lambda.max),xlim=c(mu.min,mu.max))
    } else{
      contour(MASS::kde2d(theta.list[[i]][,2], theta.list[[i]][,4]),nlevels = 10,add=T,col=cols[i])
    }
  }
  points(truth[2],truth[4],col='red',pch=16)
  for(i in 1:n.theta){
    if(i==1){
      contour(MASS::kde2d(theta.list[[i]][,3]^2, theta.list[[i]][,4]),nlevels = 10,col=cols[i],ylim=c(0,lambda.max),xlim=c(var.min,var.max))
    } else{
      contour(MASS::kde2d(theta.list[[i]][,3]^2, theta.list[[i]][,4]),nlevels = 10,add=T,col=cols[i])
    }
  }
  points(truth[3],truth[4],col='red',pch=16)
  
  lambda.dens = list(n.theta)
  max = 0
  for(i in 1:n.theta){
    lambda.dens[[i]] = density(theta.list[[i]][,4],from=0,to=lambda.max)
    M = max(lambda.dens[[i]]$y)
    if(M>max)
      max = M
  }
  plot(lambda.x,lambda.prior,type='l',xlim=c(0,lambda.max),ylim=c(0,max),lty=2,xlab='',yaxt='n')
  for(i in 1:n.theta){
    lines(lambda.dens[[i]],col=cols[i])
  }
  abline(v=truth[4],col='red',lty=2)
  legend(inset=c(.05,.05),'topright',legend=expression(lambda),cex=1.5,bty='n')
}
