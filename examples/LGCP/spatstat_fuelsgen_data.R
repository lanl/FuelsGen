library(spatstat)
library(fuelsgen)

# fuelsgen::plot_fuels(data.cl.naip[[1]])
# reps = data.cl.naip[[1]]$reps
# n.obs = data.cl.naip[[1]]$nshrub
# ppp.list = vector(mode='list',length = reps)
# for(i in 1:reps)){
#   ppp.list[[i]] = ppp(x=data.cl.naip[[1]]$dat[[i]]$X,
#                       y=data.cl.naip[[1]]$dat[[i]]$Y,
#                       c(0,20),c(0,20))
# }
set.seed(2)
fuels = fuelsgen::gen_fuels(20,20,.2,1,.1,reps=10,heterogeneity = 5)
plot_fuels(fuels)
reps = fuels$reps
n.obs = numeric(reps)

ppp.list = vector(mode='list',length = reps)
for(i in 1:reps){
  ppp.list[[i]] = ppp(x=fuels$dat[[i]]$X,
                      y=fuels$dat[[i]]$Y,
                      c(0,20),c(0,20))
  n.obs[i] = nrow(fuels$dat[[i]])
}

G = hyperframe(shrubs = ppp.list)
# mppm function allows for fitting to multiple realizations but does not support log gaussian cox processes

# Fit the trend part of the model to your replicated point pattern data using mppm, for example m <- mppm(Y ~ X, data=G)
mod = mppm(shrubs ~ 1, G, interaction = Poisson())
summary(mod)
par(mfrow=c(3,4))
for(i in 1:12){
  plot(simulate(mod,1))
}
mppm.sims = simulate.mppm(mod,12)
plot(mppm.sims)


# Attempt to hack a LGCP model to use multiple realizations

# Extract the fitted intensities for each point pattern using predict.mppm
pred = predict(mod,type='lambda')

# For each point pattern, using the corresponding intensity obtained from the model, 
# compute the inhomogeneous K function using Kinhom (with argument ratio=TRUE)
K = lapply(1:reps, function(i) Kinhom(X=ppp.list[[i]],lambda=pred$cif[[i]],ratio=T))

# Combine the K functions using pool
Kpool = pool(K[[1]],K[[2]],K[[3]],K[[4]],K[[5]],K[[6]],K[[7]],K[[8]],K[[9]],K[[10]])

# Estimate the cluster parameters of the LGCP by applying lgcp.estK to the pooled K function.
cox.mod = lgcp.estK(Kpool, covmodel =  list('gauss'))
cox.mod$modelpar
mu_func = as.im(function(x,y){
  p = coef(mod)
  p[1]# + p[2]*x + p[3]*y + p[4]*x^2 + p[5]*x*y + p[6]*y^2
},W=owin(c(0,20),c(0,20)))

nsims = 1000
tmp = rLGCP(model='gauss',
            mu = coef(mod)[1], 
            var=cox.mod$par[1],scale=cox.mod$par[2],
      win = owin(c(0,20),c(0,20)),nsim = nsims)

n.sim = numeric(nsims)
for(i in 1:nsims){
  n.sim[i] = tmp[[i]]$n
}
hist(n.sim)
summary(n.sim)

plot(tmp)
# fuelsgen::plot_fuels(data.cl.naip[[1]])
fuelsgen::plot_fuels(fuels)
hist(n.sim,col=adjustcolor( "darkorange", alpha.f = 0.5),freq=F)
hist(n.obs,add=T,col=adjustcolor( "cornflowerblue", alpha.f = 0.5),freq=F)
