library(spatstat.random)

library(spatstat) 
setwd("~/fuelsgen/examples/eucalyptus")
load("Eucalyptus sparsifolia Atlas 2012.RData") #Contains X and Y 
load("Quad1000.RData") #Contains quad 
ux = sort(unique(quad$X)) 
uy = sort(unique(quad$Y)) 
nx = length(ux) 
ny = length(uy) 
col.ref = match(quad$X, ux) 
row.ref = match(quad$Y, uy) 
all.vec = rep(NA, max(row.ref)*max(col.ref)) 
vec.ref = (col.ref- 1)*max(row.ref) + row.ref 
all.vec[vec.ref] = 1 
Sydney.mask = matrix(all.vec, max(row.ref), max(col.ref), dimnames = list(uy, ux)) 
Sydney.win = as.owin(im(Sydney.mask, xcol = ux, yrow = uy))

ppp.dat = ppp(X, Y, window = Sydney.win, check = FALSE) 
quads = ppp(quad$X, quad$Y, window = Sydney.win) 
Q = quadscheme(data = ppp.dat, dummy = quads, method = "grid", ntile = c(nx, ny), npix = c(nx, ny))

X.des = cbind(poly(quad$FC, quad$MNT, quad$MXT, quad$Rain, degree = 2, raw = TRUE), 
              poly(sqrt(quad$D.Main), sqrt(quad$D.Urb), degree = 2, raw = TRUE), quad$soil) 
int.list = list() 
for (i in 1:dim(X.des)[2]){ 
  all.vec = rep(NA, max(row.ref)*max(col.ref)) 
  vec.ref = (col.ref- 1)*max(row.ref) + row.ref 
  all.vec[vec.ref] = X.des[,i] 
  int.list[[i]] = im(matrix(all.vec, max(row.ref), max(col.ref), dimnames = list(uy, ux)), xcol = ux, yrow = uy) 
} 
names(int.list) = paste("V", 1:dim(X.des)[2], sep = "") 
pred.list = int.list 
set.0 = 15:19 #Variables to set to 0 
for (v in set.0){ 
  pred.list[[v]]$v = 0*pred.list[[v]]$v 
}

int.form = as.formula(paste("~", paste(names(int.list), collapse = "+"))) 
ft.int = ppm(Q, trend = as.formula(int.form), covariates = int.list)

tmp = simulate.ppm(ft.int,nsim = 10)
plot(tmp)
