load(here::here("examples/publication/Figure_6_7_8/sfnf_4_data.RData"))
truth = c(NA,mean(sfnf_4$r),var(sfnf_4$r),nrow(sfnf_4)/30/30)
load(here::here("examples/publication/Figure_6_7_8/sfnf_4_posterior.RData"))

# Figure 6
png('examples/publication/Figure_6_7_8/sfnf_4.png',width=5,height=5,units = 'in',res = 1200)
fuel = list()
fuel$dat = list()
fuel$dat[[1]] = sfnf_4
fuel$reps = 1; fuel$dimX = 30; fuel$dimY = 30
plot_fuels(fuel)
dev.off()
# Figure 7
contour_pairs(list(theta),truth=truth,var.max=.5,labels = c('post'))

# Figure 8
par(usr=c(0,30,0,30))
set.seed(8)
n = 5
rho = runif(n,0,10)
mu = truncdist::rtrunc(n,'norm',a=0,b=3,mean=1.5,.5)
var = rgamma(n,1,scale=1e-3)
lambda = rgamma(n,1,1/truth[4])
prior_realization = vector(mode='list',length=n)
par(mfrow=c(3,5),mar=c(.5,.5,.5,.5))
for(i in 1:n){
  prior_realization[[i]] = gen_fuels(30,30,lambda[i],mu[i],sqrt(var[i]),NULL,NULL,rho[i],1,NULL,NULL,NULL,1,100,1)
  plot(x=c(),y=c(),xlim=c(0,30),ylim=c(0,30),xlab='X',ylab='Y',xaxt='n',yaxt='n',asp=1)
  for(j in 1:nrow(prior_realization[[i]]$dat[[1]])){
    plotrix::draw.circle(prior_realization[[i]]$dat[[1]]$X[j],prior_realization[[i]]$dat[[1]]$Y[j],prior_realization[[i]]$dat[[1]]$r[j],col='darkgreen')
  }
}
plot.new(); plot.new()
plot(x=c(),y=c(),xlim=c(0,30),ylim=c(0,30),xlab='X',ylab='Y',xaxt='n',yaxt='n',asp=1)
for(j in 1:nrow(fuel$dat[[1]])){
  plotrix::draw.circle(fuel$dat[[1]]$X[j],fuel$dat[[1]]$Y[j],fuel$dat[[1]]$r[j],col=adjustcolor('grey',alpha.f = .5))
}
plot.new(); plot.new()
rho = sample(theta[,1],n)
mu = sample(theta[,2],n)
sigma = sample(theta[,3],n)
lambda = sample(theta[,4],n)
post_realization = vector(mode='list',length=n)
for(i in 1:n){
  post_realization[[i]] = gen_fuels(30,30,lambda[i],mu[i],sigma[i],NULL,NULL,rho[i],1,NULL,NULL,NULL,1,100,i)
  plot(x=c(),y=c(),xlim=c(0,30),ylim=c(0,30),xlab='X',ylab='Y',xaxt='n',yaxt='n',asp=1)
  for(j in 1:nrow(post_realization[[i]]$dat[[1]])){
    plotrix::draw.circle(post_realization[[i]]$dat[[1]]$X[j],post_realization[[i]]$dat[[1]]$Y[j],post_realization[[i]]$dat[[1]]$r[j],col='darkgreen')
  }
}
#dev.off()