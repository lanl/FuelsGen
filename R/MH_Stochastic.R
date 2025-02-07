Metro_Hastings_Stochastic = function (li_func, pars, prop_sigma = NULL, par_names = NULL, 
                                      iterations = 50000, burn_in = 1000, adapt_par = c(100, 20, 0.5, 0.75), quiet = FALSE, ...) 
{
  if (!is.finite(li_func(pars, ...))) 
    stop("Seed parameter values <pars> are not in the defined parameter space.  Try new starting values for <pars>.")
  if (is.null(par_names)) 
    par_names <- letters[1:length(pars)]
  if (!is.null(dim(prop_sigma))) {
    if ((dim(prop_sigma)[1] != length(pars) || dim(prop_sigma)[2] != 
         length(pars)) && !is.null(prop_sigma)) 
      stop("prop_sigma not of dimension length(pars) x length(pars)")
  }
  if (is.null(prop_sigma)) {
    if (length(pars) != 1) {
      fit <- optim(pars, li_func, control = list(fnscale = -1), 
                   hessian = TRUE, ...)
      fisher_info <- solve(-fit$hessian)
      prop_sigma <- sqrt(diag(fisher_info))
      prop_sigma <- diag(prop_sigma)
    }
    else {
      prop_sigma <- 1 + pars/2
    }
  }
  prop_sigma <- makePositiveDefinite(prop_sigma)
  mu <- pars
  pi_X <- li_func(pars, ...)
  k_X <- pars
  trace <- array(dim = c(iterations, length(pars)))
  deviance <- array(dim = iterations)
  announce <- floor(seq(iterations/100, iterations, length.out = 100))
  for (i in 1:iterations) {
    # cat(i,' ')
    k_Y <- mvrnorm(1, mu = k_X, Sigma = prop_sigma)
    # cat('prop:',round(k_Y,2),' ')
    pi_Y <- li_func(k_Y, ...)
    # cat('prop llh:',round(pi_Y,2),' ')
    a_X_Y = (pi_Y) - (pi_X)
    if (is.nan(a_X_Y)) 
      a_X_Y <- -Inf
    if (log(runif(1, 0, 1)) <= a_X_Y) {
      k_X = k_Y
      # cat('cur:',round(k_X,2),' ')
      pi_X = pi_Y
    } else{
      # if not accept recalc current state due to stochastic llh function
      # cat('cur:',round(k_X,2),' ')
      pi_X = li_func(k_X, ...)
    }
    # cat('cur llh:',round(pi_X,2),'\n ')
    trace[i, ] <- k_X
    deviance[i] <- (-2 * pi_X)
    if (i > adapt_par[1] && i%%adapt_par[2] == 0 && i < (adapt_par[4] * 
                                                         iterations)) {
      len <- floor(i * adapt_par[3]):i
      x <- trace[len, ]
      N <- length(len)
      p_sigma <- (N - 1) * var(x)/N
      p_sigma <- makePositiveDefinite(p_sigma)
      if (!(0 %in% p_sigma)) 
        prop_sigma <- p_sigma
    }
    if (!quiet && i %in% announce) 
      print(paste("updating: ", i/iterations * 100, "%", 
                  sep = ""))
    # if(i %% 10 == 0){save(trace,file='~/fuelsgen/examples/LGCP/tmp.RData')}
  }
  trace <- trace[burn_in:iterations, ]
  DIC <- NULL
  D_bar <- mean(deviance[burn_in:iterations])
  if (length(pars) > 1) {
    theta_bar <- sapply(1:length(pars), function(x) {
      mean(trace[, x])
    })
  }
  else theta_bar <- mean(trace)
  D_hat <- li_func(theta_bar, ...)
  pD <- D_bar - D_hat
  DIC <- D_hat + 2 * pD
  if (length(pars) > 1) 
    accept_rate <- length(unique(trace[, 1]))/(iterations - 
                                                 burn_in)
  else accept_rate <- length(unique(trace))/(iterations - burn_in)
  val <- list(trace = trace, prop_sigma = prop_sigma, par_names = par_names, 
              DIC = DIC, acceptance_rate = accept_rate)
  class(val) <- "MHposterior"
  return(val)
}