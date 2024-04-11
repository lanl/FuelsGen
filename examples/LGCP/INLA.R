m = 100
points = matrix(runif(m*2),m,2)
mesh = inla.mesh.create.helper( points=points, cutoff=0.05, offset=c(0.1,0.4), max.edge=c(0.05,0.5) )

reps = 10
fuels = fuelsgen::gen_fuels(20,20,.1,1,.1,heterogeneity = 5,reps=reps)
fuelsgen::plot_fuels(fuels)
points = matrix(nrow=0,ncol=3)
for(i in 1:reps){
  points = rbind(points,cbind(fuels$dat[[i]][1:2],i))
}
loc.bnd = matrix(c(0,0, 20,0, 20,20, 0,20), 4, 2, byrow=TRUE)
segm.bnd = inla.mesh.segment(loc.bnd)
mesh = inla.mesh.create.helper( points=points[,1:2], 
                                boundary=segm.bnd,
                                cutoff=1, offset=c(1,.5), max.edge=c(.1,1) )
par(mfrow=c(1,1))
plot(mesh)
points(points[,1],points[,2])

sigma0 = 1
range0 = 3
kappa0 = sqrt(8)/range0
tau0 = 1/(sqrt(4*pi)*kappa0*sigma0)
spde=inla.spde2.matern(mesh,B.tau=cbind(log(tau0),1,0), B.kappa=cbind(log(kappa0),0,1), theta.prior.prec=1)

# Q=inla.spde2.precision( spde,theta=c(0,0))
# x=as.vector(inla.qsample(n=1,Q))
# proj = inla.mesh.projector(mesh) 
# image(inla.mesh.project( proj,field=x))

# example with 2 replicates at same location
A=inla.spde.make.A( mesh, 
                    loc = as.matrix(points[,1:2]),            # 2d matrix of observed locations stack realizations
                    repl = points[,3])
                    #index=rep(1:m,times=2),# vector of row indices for each location (remove this for my data because locations are not shared)
                    #repl=rep(1:2,each=m) ) # vector indicating which replicate each observation location belongs to

# something is probably wrong with the effects
stk.sp = inla.stack(data = list(y = 1, e = 0),
                    effects=list( c(mesh.index,list(offset=1)), list(Intercept=1)),
                    A = list(A, 1), 
                    tag = 'sp')

Q=inla.spde.precision( spde,theta=c(0,0)) 
x=(inla.qsample(n=9,Q)) 
par(mfrow=c(3,3))
for(i in 1:9){
  image(inla.mesh.project( proj,field=x[,1]))
  points(fuels$dat[[i]][,1:2]/20)
}

covariate = rep(1,nrow(points))#rnorm(m*2) 
# y = 5 + covariate*2 + as.vector(A %*% x) + rnorm(m*2)*0.01

mesh.index=inla.spde.make.index( name="field", n.spde=mesh$n, n.repl=2)

st.est=inla.stack( data=list(y=y), A=list(A,1), 
                   effects=list( c(mesh.index,list(offset=1)), list(cov=covariate)), tag="est")

st.pred=inla.stack( data=list(y=NA), A=list(1), effects=list( c(mesh.index,list(offset=1))), tag="pred") 
# We can now join the estimation and prediction stack into a single stack,
stack = inla.stack(st.est,st.pred)

formula = y ~-1 + offset + cov + f(field, model=spde, replicate=field.repl) 
inla.result = inla(formula, data=inla.stack.data(stack), 
                   family="normal", control.predictor= list(A=inla.stack.A(stack), compute=TRUE))
result = inla.spde2.result( inla.result, "field", spde)
plot(result$"marginals.range.nominal"[[1]])

index=inla.stack.index( stack,"pred")$data 
image(inla.mesh.project(proj, inla.result$summary.linear.predictor$mean[ index[mesh.index$field.repl==1]]))
