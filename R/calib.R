llh = function(y, theta, dimX, dimY, ldetS, Sinv, mets_info, gen_reps = 1, avg = 'mets', gen_parallel = F, mets_parallel = F)
{
  sim_fuels = gen_fuels(dimX, dimY,
                 theta[4],
                 theta[2], theta[3],
                 height = NULL, sd_height = 0,
                 theta[1], heterogeneity.scale = 1,
                 sp.cov.locs = NULL, sp.cov.vals = NULL, sp.cov.scale = NULL,
                 reps=gen_reps, GP.init.size = 100, seed=NULL, logis.scale=.217622, verbose=F, parallel = gen_parallel)

  metrics = get_mets(sim_fuels, mets_info, mets_parallel)
  
  if(avg == 'mets'){
    # average over the metrics and compute the likelihood once
    metrics = colMeans(metrics)
    llh = 0
    for(i in 1:nrow(y)){
      # loop over replicate observations
      llh = llh - 0.5 * (ldetS + ((y[i,,drop=F] - metrics) %*% Sinv %*% t(y[i,,drop=F] - metrics)))
    }
  } else if(avg == 'llh'){
    llh = numeric(gen_reps)
    for(j in 1:length(llh)){
      # loop over generated metrics
      llh[j] = 0
      for(i in 1:nrow(y)){
        # loop over replicate observations
        llh[j] = llh[j] - 0.5 * (ldetS + ((y[i,,drop=F] - metrics[j,]) %*% Sinv %*% t(y[i,,drop=F] - metrics[j,])))
      }
    }
  }
  
  return(mean(llh))
}

lprior = function(theta,params){
  if(is.null(params$rho_shape)){
    # flat prior but limited domain
    rho_prior = dunif(theta[1],0,params$rho_max)
  } else{
    # gamma prior with q95 set to be approximately rho_max
    rho_prior = dgamma(theta[1],params$rho_shape,params$rho_rate,log=T)
  }
  prior = rho_prior + 
    truncdist::dtrunc(theta[2],'norm',a=0,b=3,mean=params$mu_mean,sd=params$mu_sd,log=T) + 
    dgamma(1/theta[3]^2,params$prec_shape,params$prec_rate,log=T) +
    dgamma(theta[4],params$lambda_shape,params$lambda_rate,log=T)
  return(prior)
}

matsplitter=function(M, r, c)
{
  rg = (row(M)-1)%/%r+1
  cg = (col(M)-1)%/%c+1
  rci = (rg-1)*max(cg) + cg
  N = prod(dim(M))/r/c
  cv = unlist(lapply(1:N, function(x) M[rci==x]))
  dim(cv)=c(r,c,N)
  cv
}
moran_geary = function(afs, rq="rook", moran=T, geary=T)
{
  nm = length(afs)
  # previously
  # dd = dim(afs)[1]
  # ids = expand.grid(1:dd, 1:dd)
  # i think this was an error which caused matrix size bugs for rectanglular domains
  dd = dim(afs)
  ids = expand.grid(1:dd[1], 1:dd[2])
  dm = sqrt(plgp::distance(ids))

  if (rq == "rook"){
    wmat = ifelse(dm > 1, 0, 1)
  } else if(rq == "queen"){
    wmat = ifelse(dm > sqrt(2), 0, 1)
  }

  ret = c()
  # Moran's I
  if(moran){
    xminusxbar = c(afs - mean(afs))
    xijm = tidyr::expand_grid(xminusxbar, xminusxbar)
    xtx = matrix((xijm[,1] * xijm[,2])[,1], ncol=nm)
    mi = ( nm * sum(wmat * xtx) ) /
      ( sum(xminusxbar^2) * sum(wmat) )
    ret = c(ret,mi)
  }
  # Geary's C
  if(geary){
    af_vec = c(afs)
    xijg = tidyr::expand_grid(af_vec, af_vec)
    ximinusxj2mat = matrix(((xijg[,1] - xijg[,2])[,1])^2, ncol = nm)
    gc = ((nm - 1) * sum(wmat * ximinusxj2mat)) /
      (2 * sum(wmat) * sum((af_vec - mean(afs))^2))
    ret = c(ret,gc)
  }

  return(ret)
}

# function to remove all shrubs that are completely contained within another shrub
# not removing these rows makes logic for computing segments and breaks significantly more difficult
# this function should only be called from the transect function
remove_rows = function(pp_j,dim)
{
  remove = numeric(nrow(pp_j))
  for(i in 1:nrow(pp_j)){
    for(j in 1:nrow(pp_j)){
      if(i==j)
        next
      if(dim == 'X'){
        if(pp_j$Ypr[i] <= pp_j$Ypr[j] & pp_j$Ymr[i] >= pp_j$Ymr[j]){
          remove[i] = 1
          # we know we are removing row i so break out of j and go to next i
          break
        }
      } else{
        if(pp_j$Xpr[i] <= pp_j$Xpr[j] & pp_j$Xmr[i] >= pp_j$Xmr[j]){
          remove[i] = 1
          # we know we are removing row i so break out of j and go to next i
          break
        }
      }
    }
  }
  return(pp_j[!remove,])
}
transect = function(t_locs,dims,dat,dimX,dimY)
{
  n_transects = length(t_locs)
  segment_lengths = vector(mode='list',n_transects)
  break_lengths = vector(mode='list',n_transects)

  for(j in 1:n_transects){
    if(dims[j]=='X'){
      x_dist_from_transect = abs(dat$X - t_locs[j])
      x_within_r_from_transect = x_dist_from_transect <= dat$r

      # calculate segment lengths
      pp_j = dat[x_within_r_from_transect,]
      if(nrow(pp_j)>0){
        # pmax(0,pp_j$Y - pp_j$r) takes care of the case where the left boundary is < 0
        # pmin(pp_j$Y + pp_j$r,dimY) takes care of the case where the right boundary is > dimY
        pp_j$Ymr = pmax(0,pp_j$Y - pp_j$r)
        pp_j$Ypr = pmin(pp_j$Y + pp_j$r,dimY)
        # sort shrubs in order of left boundary
        pp_j = pp_j[order(pp_j$Ymr, pp_j$Ypr),]
        # remove shrubs that are completely contained within previous shrub
        pp_j = remove_rows(pp_j,dims[j])
        segment_lengths[[j]] = pp_j$Ypr[1] - pp_j$Ymr[1]
        if(pp_j$Ymr[1]>0){
          break_lengths[[j]] = c(pp_j$Ymr[1])
        } else{
          break_lengths[[j]] = numeric()
        }
        if(nrow(pp_j)>1){
          k = 1
          for(i in 2:nrow(pp_j)){
            # this is an oversimplification that does not account for segment length of the line bisecting the circle
            if(pp_j$Ypr[i-1] == dimY){
              # end if we reached the boundary
              break
            }
            # check if circle is overlapping previous
            if( pp_j$Ymr[i] <= pp_j$Ypr[i-1]){
              # update segment lengths to reflect total length of these connected circles
              segment_lengths[[j]][k] = segment_lengths[[j]][k] + (pp_j$Ypr[i] - pp_j$Ypr[i-1])
            } else if(pp_j$Ymr[i] > pp_j$Ypr[i-1]){
              segment_lengths[[j]] = c(segment_lengths[[j]], pp_j$Ypr[i] - pp_j$Ymr[i])
              break_lengths[[j]] = c(break_lengths[[j]],pp_j$Ymr[i] - pp_j$Ypr[i-1])
              k = k + 1
            }
          }
        }
        d_to_bound = dimY - pp_j$Ypr[nrow(pp_j)]
        if(d_to_bound>0){
          break_lengths[[j]] = c(break_lengths[[j]],d_to_bound)
        }
      } else{
        # no shrubs touch transect a.k.a. one break of length dimY and no segments
        segment_lengths[[j]] = 0
        break_lengths[[j]] = dimY
      }

    } else if(dims[j]=='Y'){
      y_dist_from_transect = abs(dat$Y - t_locs[j])
      y_within_r_from_transect = y_dist_from_transect <= dat$r

      # calculate segment lengths
      pp_j = dat[y_within_r_from_transect,]
      if(nrow(pp_j)>0){
        pp_j$Xmr = pmax(0,pp_j$X - pp_j$r)
        pp_j$Xpr = pmin(pp_j$X + pp_j$r,dimX)
        pp_j = pp_j[order(pp_j$Xmr),]
        # remove shrubs that are completely contained within previous shrub
        pp_j = remove_rows(pp_j,dims[j])
        segment_lengths[[j]] = pp_j$Xpr[1] - pp_j$Xmr[1]
        if(pp_j$Xmr[1]>0){
          break_lengths[[j]] = c(pp_j$Xmr[1])
        } else{
          break_lengths[[j]] = numeric()
        }
        if(nrow(pp_j)>1){
          k = 1
          for(i in 2:nrow(pp_j)){
            # this is an oversimplification that does not account for segment length of the line bisecting the circle
            if(pp_j$Xpr[i-1] == dimX){
              break
            }
            if(pp_j$Xpr[i] <= pp_j$Xpr[i-1]){
              # move to next as this shrub is completely contained in the previous shrub
              next
            }
            # check if cirle is overlapping previous
            if( pp_j$Xmr[i] <= pp_j$Xpr[i-1] ){
              # update segment lengths to reflect total length of these connected circles
              segment_lengths[[j]][k] = segment_lengths[[j]][k] + (pp_j$Xpr[i] - pp_j$Xpr[i-1])
            } else{
              segment_lengths[[j]] = c(segment_lengths[[j]],pp_j$Xpr[i] - pp_j$Xmr[i])
              break_lengths[[j]] = c(break_lengths[[j]],pp_j$Xmr[i] - pp_j$Xpr[i-1])
              k = k + 1
            }
          }
        }
        d_to_bound = dimX - pp_j$Xpr[nrow(pp_j)]
        if(d_to_bound>0){
          break_lengths[[j]] = c(break_lengths[[j]],d_to_bound)
        }
      } else{
        segment_lengths[[j]] = 0
        break_lengths[[j]] = dimX
      }
    } else{
      stop("dim must be one of 'X' or 'Y'")
    }
  }
  for(i in 1:n_transects){
    seg_sum = ifelse(length(segment_lengths[[i]])>0, sum(segment_lengths[[i]]), 0)
    break_sum = ifelse(length(break_lengths[[i]])>0, sum(break_lengths[[i]]), 0)
    max = ifelse(dims[i]=='X',dimY,dimX)
    if(abs(seg_sum + break_sum - max)>1e-4){
      #cat('Location: ',t_locs[i],'\nDim: ',dims[i],'\nSegments: ',round(segment_lengths[[i]],2),'\nBreaks: ',round(break_lengths[[i]],2),'\nSum: ',seg_sum+break_sum,'\n')
      warning('Error in transect ',i,': sum of segments and breaks != ',max,'\n')
    }
  }
  return(list(segments = segment_lengths,
              breaks = break_lengths))
}

get_mets = function(fuels, info, parallel = F)
{
  if(fuels$reps>1 & parallel){
    metrics = mets_parallel(fuels, info)
  } else{
    if(fuels$reps>1){
      metrics = c()
      for(i in 1:fuels$reps){
        metrics = rbind(metrics,mets(fuels$dat[[i]], info, fuels$dimX, fuels$dimY))
      }
    } else{
      metrics = mets(fuels$dat, info, fuels$dimX, fuels$dimY)
    }
  }
  return(metrics)
}

mets_parallel = function(fuels,info)
{
  #cl = parallel::makeCluster(min(fuels$reps,parallel::detectCores()))
  #doParallel::registerDoParallel(cl)
  metrics = foreach::foreach(i=1:fuels$reps,.combine = 'rbind') %dopar% mets(fuels$dat[[i]],info,fuels$dimX,fuels$dimY)
  #parallel::stopCluster(cl)
  return(metrics)
}

mets = function(dat, info, dimX, dimY)
{
  metrics = c()
  n = nrow(dat)
  if (n == 0){
    return(rep(0,11))
  }

  #--- NCC ---#
  if(info$ncc){
    ptm = proc.time()[3]
    r_mat =  matrix(rep(dat$r, n), ncol=n, byrow=T)
    r_pairs = r_mat + t(r_mat)
    cent_dist = sqrt(plgp::distance(dat[, c("X", "Y")]))
    incidence = matrix(as.numeric(cent_dist <= r_pairs), nrow=n, ncol=n, byrow = T)

    ig = igraph::graph_from_adjacency_matrix(incidence * (1 / cent_dist),
                                             weighted = T, diag = F, mode="undirected")
    #ds = igraph::distances(ig)
    #ds[is.infinite(ds)] = 0
    #max_d = max(ds)

    ncc.val = igraph::components(ig)$no
    ncc.time = proc.time()[3] - ptm
  }
  #--- N holes ---#
  if(info$nholes){
    ptm = proc.time()[3]
    if(!info$ncc){
      # these are normally computed in ncc so if we don't do ncc we must compute them here
      r_mat =  matrix(rep(dat$r, n), ncol=n, byrow=T)
      r_pairs = r_mat + t(r_mat)
      cent_dist = sqrt(plgp::distance(dat[, c("X", "Y")]))
    }
    # why a threshold of 5?
    homo = TDAstats::calculate_homology(cent_dist / r_pairs,
                                        format = "distmat", threshold = 5, return_df = T)
    #homo = TDApplied::PyH(cent_dist / r_pairs,
    #                      distance_mat = TRUE, thresh = 5, ripser=ripser)
    nholes.val = nrow(dplyr::filter(homo, dimension==1 & birth < 1 & death > 1))
    nholes.time = proc.time()[3] - ptm
  }
  #--- Gridded area ---#
  if(info$grid.area){
    ptm = proc.time()[3]
    dx = 1
    dy = 1
    stepsize = 0.25
    xstart = 0

    area_fractions = matrix(ncol=dimX / dx, nrow=dimY / dy)

    for (i in 1:(dimX / dx)) {
      ystart = 0
      for (j in 1:(dimY / dx)) {
        my_grid = pracma::meshgrid(seq(xstart, xstart + dx - stepsize, by=stepsize),
                                   seq(ystart, ystart + dy - stepsize, by=stepsize))
        mc_d = sqrt(plgp::distance(cbind(c(my_grid$X), c(my_grid$Y)), dat[,c("X", "Y")]))

        area_fractions[j,i] = mean(apply(mc_d < dat$r, 1, any))
        ystart = ystart + dy
      }
      xstart = xstart + dx
    }
    af.time = proc.time()[3] - ptm
    total_area = mean(area_fractions) * dimX * dimY
    empty_cells = sum(area_fractions == 0)
    full_cells = sum(area_fractions == 1)
    af_var = var(c(area_fractions))
  }

  #--- Moran I ---#
  if(info$moran | info$geary){
    ptm = proc.time()[3]
    af0 = ifelse(area_fractions == 0, 0, 1)
    if(info$rook){
      r1 = moran_geary(area_fractions, "rook", info$moran, info$geary)
      r01 = moran_geary(af0, "rook", info$moran, info$geary)
    }
    if(info$queen){
      q1 = moran_geary(area_fractions, "queen", info$moran, info$geary)
      q01 = moran_geary(af0, "queen", info$moran, info$geary)
    }

    af_half = matrix(apply(matsplitter(area_fractions, 2, 2), 3, mean),
                     ncol = dimX / (dx * 2), byrow = T)
    af_half0 = ifelse(af_half == 0, 0, 1)

    if(info$rook){
      r2 = moran_geary(af_half, "rook", info$moran, info$geary)
      r02 = moran_geary(af_half0, "rook", info$moran, info$geary)
    }
    if(info$queen){
      q2 = moran_geary(af_half, "queen", info$moran, info$geary)
      q02 = moran_geary(af_half0, "queen", info$moran, info$geary)
    }

    # what does this mean? These are 2-vectors so x[is.nan(x)] could be length 1 or length 2
    # If it's length 1, we are setting it to zero always and if its length 2 we set to c(0,1)
    # so this is saying in words "if 1 of them is NaN, make that one zero, and if both are Nan, set the moran to 0 and the geary to 1"
    # Need to ask Braden since it's not intuative. Try with commented out to see if NaN's are a problem,.
    #
    # this  might be saying set moran NaN's to 0 and Geary NaN's to 1
    if(info$rook){
      r01[is.nan(r01)] = c(0, 1)
      r02[is.nan(r02)] = c(0, 1)
    }
    if(info$queen){
      q1[is.nan(q1)] = 1
      q2[is.nan(q2)] = 1
      q01[is.nan(q01)] = 1
      q02[is.nan(q02)] = 1
    }
    moran.time = proc.time()[3] - ptm
  }
  #--- Perimeter ---#
  if(info$perim){
    ptm = proc.time()[3]
    r_sum = sum(dat$r)
    store = c(0, 0)

    for (i in 1:n){
      n_perm = max(c(1, round((dat$r[i] / r_sum) * 5000)))
      rads = seq(0, 2 * pi, length.out=(n_perm + 1))[1:n_perm]
      xy = rep(c(dat[i,1], dat[i,2]), each=n_perm)
      this = cbind(dat$r[i] * cos(rads) + dat[i,1], dat$r[i] * sin(rads) + dat[i,2])
      dis = sqrt(plgp::distance(this, dat[,c("X", "Y")]))
      perim_pts = sum(abs(matrixStats::rowMins(dis) - dat$r[i]) < 0.000001)

      store[1] = store[1] + n_perm
      store[2] = store[2] + perim_pts
    }

    pct_perm = store[2] / store[1]
    perimeter = pct_perm * sum(2 * dat$r * pi)
    perim.time = proc.time()[3] - ptm
  }
  #--- Transects ---#
  if(info$transect){
    ptm = proc.time()[3]
    if(is.null(info$dim.trans)){
      if(info$n.trans==1){
        info$dim.trans = 'X'
      } else if(info$n.trans>1){
        if(info$n.trans%%2==0){
          info$dim.trans = c(rep('X',as.integer(info$n.trans/2)),rep('Y',as.integer(info$n.trans/2)))
        } else{
          info$dim.trans = c(rep('X',as.integer(info$n.trans/2)+1),rep('Y',as.integer(info$n.trans/2)))
        }
      } else{
        stop('n.trans<1')
      }
    }
    if(is.null(info$loc.trans)){
      # evenly space the transects in each dimension
      if(info$n.trans==1){
        info$loc.trans = dimX/2
      } else{
        if(info$n.trans%%2==0){
          n.x = info$n.trans/2; n.y = info$n.trans/2
        } else{
          n.x = as.integer(info$n.trans/2) + 1; n.y = as.integer(info$n.trans/2)
        }
        x_start = dimX/(info$n.trans/2+1); x_end = dimX - x_start
        y_start = dimY/(info$n.trans/2+1); y_end = dimX - x_start
        info$loc.trans = c(seq(x_start,x_end,length.out=n.x),seq(y_start,y_end,length.out=n.y))
      }
    }
    t = transect(info$loc.trans,info$dim.trans,dat,dimX,dimY)
    segments = unlist(t$segments)
    breaks = unlist(t$breaks)
    transect.time = proc.time()[3] - ptm
  }

  if(info$ncc){
    metrics = c(metrics,ncc.val)
  }
  if(info$nholes){
    metrics = c(metrics,nholes.val)
  }
  if(info$grid.area){
    metrics = c(metrics, total_area, full_cells, empty_cells, af_var * 100)
  }
  if(info$moran | info$geary){
    if(info$rook){
      metrics = c(metrics,100*c(r1,r2,r01,r02))
    }
    if(info$queen){
      metrics = c(metrics,100*c(q1,q2,q01,q02))
    }
  }
  if(info$perim){
    metrics = c(metrics,perimeter)
  }
  if(info$transect){
    nseg = length(segments)
    nbreak = length(breaks)
    metrics = c(metrics,nseg,ifelse(nseg>0,mean(segments),0),ifelse(nseg>1,var(segments),0),
                     nbreak,ifelse(nbreak>0,mean(breaks),0),ifelse(nbreak>1,var(breaks),0))
  }
  # include emperical estimates of parameters as metrics
  metrics = c(metrics,2*n/dimX/dimY,mean(dat$r),sd(dat$r))
  # returns = list()
  # returns$mets = c()
  # returns$time = c()
  # returns$names = c()
  # if(info$ncc){
  #   returns$time = c(returns$time,ncc.time)
  #   returns$mets = c(returns$mets,ncc.val)
  #   returns$names = c(returns$names,'ncc')
  # }
  # if(info$nholes){
  #   returns$time = c(returns$time,nholes.time)
  #   returns$mets = c(returns$mets,nholes.val)
  #   returns$names = c(returns$names,'nholes')
  # }
  # if(info$grid.area){
  #   returns$time = c(returns$time,af.time)
  #   returns$mets = c(returns$mets, total_area, full_cells, empty_cells, af_var * 100)
  #   returns$names = c(returns$names,'total_area','full_cells','empty_cells','af_var')
  # }
  # if(info$moran | info$geary){
  #   returns$time = c(returns$time,moran.time)
  #   if(info$rook){
  #     returns$mets = c(returns$mets,100*c(r1,r2,r01,r02))
  #     returns$names = c(returns$names,'r1','r2','r01','r02')
  #   }
  #   if(info$queen){
  #     returns$mets = c(returns$mets,100*c(q1,q2,q01,q02))
  #     returns$names = c(returns$names,'q1','q2','q01','q02')
  #   }
  # }
  # if(info$perim){
  #   returns$time = c(returns$time,perim.time)
  #   returns$mets = c(returns$mets,perimeter)
  #   returns$names = c(returns$names,'perim')
  # }
  # if(info$transect){
  #   returns$time = c(returns$time,transect.time)
  #   nseg = length(segments)
  #   nbreak = length(breaks)
  #   returns$mets = c(returns$mets,nseg,ifelse(nseg>0,mean(segments),0),ifelse(nseg>1,var(segments),0),
  #                    nbreak,ifelse(nbreak>0,mean(breaks),0),ifelse(nbreak>1,var(breaks),0))
  #   returns$names = c(returns$names,'transects')
  # }
  # # include emperical estimates of parameters as metrics
  # returns$mets = c(returns$mets,2*n/dimX/dimY,mean(dat$r),sd(dat$r))
  # returns$names = c(returns$names,'lambda_est','mu_est','sigma_est')
  return(metrics)
}

est_cov = function(data.type = 'perfect',
                   data = NULL,
                   rho = NULL, mu = 1, sigma = .1, lambda = .5, dimX = 10, dimY = 10,
                   prior.samples = 100, reps = 100,
                   ncc = TRUE,
                   nholes = TRUE,
                   grid.area = TRUE,
                   moran = FALSE, geary = TRUE,
                   rook = FALSE, queen = TRUE,
                   perim = FALSE,
                   transect = TRUE, n.trans = 4, dim.trans = NULL, loc.trans = NULL)
{
  cl = parallel::makeCluster(min(reps,parallel::detectCores()))
  doParallel::registerDoParallel(cl)

  # Save information for calculating metrics
  mets_info = list()
  mets_info$ncc = ncc
  mets_info$nholes = nholes
  mets_info$grid.area = grid.area
  mets_info$moran = moran; mets_info$geary = geary
  mets_info$rook = rook; mets_info$queen = queen
  mets_info$perim = perim
  mets_info$transect = transect
  mets_info$n.trans = n.trans
  mets_info$dim.trans = dim.trans
  mets_info$loc.trans = loc.trans

  if(data.type == 'perfect'){
    if(is.null(rho)){
      stop('Enter rho to generate perfect data.')
    }
    theta = c(rho,mu,sigma,lambda)
    truth_realizations = foreach::foreach(i=1:reps) %dopar%
                gen_fuels(dimX, dimY,
                          theta[4],
                          theta[2], theta[3],
                          height = NULL, sd_height = 0,
                          theta[1], heterogeneity.scale = 1,
                          sp.cov.locs = NULL, sp.cov.vals = NULL, sp.cov.scale = NULL,
                          reps=1, GP.init.size = 100, seed=NULL, logis.scale=.217622, verbose=F)


    y_obs = colMeans(foreach::foreach(i=1:reps,.combine='rbind') %dopar%
                       mets(truth_realizations[[i]],mets_info)$mets)

    # estimate parameters and sample from informative priors for covariance
    mu_est = numeric(reps)
    sigma_est = numeric(reps)
    lambda_est = numeric(reps)
    for(i in 1:reps){
      mu_est[i] = mean(truth_realizations[[i]]$dat$r)
      sigma_est[i] = sd(truth_realizations[[i]]$dat$r)
      lambda_est[i] = 2*nrow(truth_realizations[[i]]$dat)/dimX/dimY
    }

    mu_est = mean(mu_est)
    sigma_est = mean(sigma_est)
    lambda_est = mean(lambda_est)

  } else if(data.type=='sfnf'){
    sfnf = list()
    sfnf$dat = read.csv(paste0(data.dir,'data.csv'),header=F)
    names(sfnf$dat) = c('X','Y','sd')
    sfnf$dat$r = 2*sfnf$dat$sd; sfnf$dat$sd = NULL
    sfnf$dat$X = sfnf$dat$X + 15
    sfnf$dat$Y = sfnf$dat$Y + 15
    ntree = nrow(sfnf$dat)
    sfnf$dimX = 30
    sfnf$dimY = 30
    y_obs = mets(sfnf,mets_info)$mets

    mu_est = mean(sfnf$dat$r)
    sigma_est = sd(sfnf$dat$r)
    lambda_est = 2*ntree/sfnf$dimX/sfnf$dimY

  } else if(data.type == 'single_sim'){
    dimX = data$dimX
    dimY = data$dimY
    ntree = nrow(data$dat)
    y_obs = mets(data,mets_info)$mets
    mu_est = mean(data$dat$r)
    sigma_est = sd(data$dat$r)
    lambda_est = ntree/dimX/dimY
    if(!is.null(rho)){
      rho_est = rho
    }
  } else{
    stop('supported data types are \'perfect\', \'sfnf\', and \'single_sim\'.')
  }

  mu_prior_samp = truncdist::rtrunc(prior.samples,'norm',a=0,mean=mu_est,sd=.1*mu_est)
  sigma_prior_samp = truncdist::rtrunc(prior.samples,'norm',a=0,mean=sigma_est,sd=.1*sigma_est)
  lambda_prior_samp = truncdist::rtrunc(prior.samples,'norm',a=0,mean=lambda_est,sd=.1*lambda_est)
  
  rho_max = max(dimX,dimY)/3
  if(is.null(rho)){
    # we don't have good information about rho, but we know that with the given kernel
    # point 3*rho apart are effectively uncorrelated, so we cannot infer larger rho's anyway
    rho_prior_samp = runif(prior.samples,0,rho_max)
    rho_est = mean(rho_prior_samp)
  } else{
    # given strong prior information about rho
    rho_prior_samp = truncdist::rtrunc(prior.samples,'norm',a=0,mean=rho_est,sd=.1*rho_est)
  }


  if(data.type=='sfnf'){
    # gamma parameters s.t. q95 ~~ rho_max
    #rho_shape = 1
    #rho_rate = .3
    #rho_prior = "gamma"
    rho_prior = "uniform"
    rho_shape = NULL
    rho_rate = NULL
  } else{
    rho_prior = "uniform"
    rho_shape = NULL
    rho_rate = NULL
  }
  # set bounds for each parameters

  # the lower bound on lambda is very important. Certain metrics calculations will fail if not enough trees are sampled.
  # The error we get with small lambda is "Point cloud must have at least 2 points and at least 2 dimensions." I need to find out
  # which metric calculation is causing this.
  # We choose to assume that the true lambda such that on average produce half the number of trees as observed are produced
  # i.e. lambda >= lambda_est/2. And similarly we assume that lambda <= lambda_est*2
  prior.lb = c(0,       0,   0, lambda_est/2)
  prior.ub = c(rho_max, 3, Inf, lambda_est*2)

  metrics = matrix(nrow=(prior.samples*reps),ncol=length(y_obs))
  for(i in 0:(prior.samples-1)){
    if(i>0 & i%%10==0){cat(i,' ')}
    j1 = i*reps + 1
    j2 = j1 + reps - 1
    theta = c(rho_prior_samp[i+1],mu_prior_samp[i+1],sigma_prior_samp[i+1],lambda_prior_samp[i+1])
    pp = foreach::foreach(k=1:reps) %dopar% gen_fuels(dimX, dimY,
                                                      theta[4],
                                                      theta[2], theta[3],
                                                      height = NULL, sd_height = 0,
                                                      theta[1], heterogeneity.scale = 1,
                                                      sp.cov.locs = NULL, sp.cov.vals = NULL, sp.cov.scale = NULL,
                                                      reps=1, GP.init.size = 100, seed=NULL, logis.scale=.217622, verbose=F)
    metrics[j1:j2,] = foreach::foreach(k=1:reps,.combine='rbind') %dopar%
      mets(pp[[k]],mets_info)$mets
  }

  Sigma = cov(metrics)
  Sinv = solve(Sigma + 1e-8*diag(dim(Sigma)[1]))
  ldetS = determinant(Sigma)$modulus

  # save everything needed to calculate priors
  prior_params = list()
  prior_params$lb = prior.lb; prior_params$ub = prior.ub
  prior_params$rho_prior = rho_prior
  prior_params$rho_max = rho_max
  prior_params$rho_shape = rho_shape
  prior_params$rho_rate = rho_rate
  prior_params$mu_prior = "trunc normal"
  prior_params$mu_mean = 1.5
  prior_params$mu_sd = .5
  prior_params$prec_prior = 'gamma'
  prior_params$prec_shape = 1
  prior_params$prec_rate = 1e-3
  prior_params$lambda_prior = 'gamma'
  prior_params$lambda_shape = 1
  prior_params$lambda_rate = 1/lambda_est
  
  parallel::stopCluster(cl)
  
  return(list('y_obs'=y_obs,
              'prior_params'=prior_params,
              'mets_info'=mets_info,
              'rho_est'=rho_est,
              'mu_est'=mu_est,
              'sigma_est'=sigma_est,
              'metrics'=metrics,
              'Sigma'=Sigma,
              'Sinv'=Sinv,
              'ldetS'=ldetS,
              'dimX'=dimX,
              'dimY'=dimY))
}

get_prior_info = function(fuel,metrics){
  if(fuel$reps==1){
    fuel$dat = list(fuel$dat)
  }
  for(i in 1:fuel$reps){
    lambda_est = numeric(fuel$reps)
    mu_est = numeric(fuel$reps)
    sigma_est = numeric(fuel$reps)
    for(i in 1:fuel$reps){
      lambda_est[i] = nrow(fuel$dat[[i]])/fuel$dimX/fuel$dimY
      mu_est[i] = mean(fuel$dat[[i]]$r)
      sigma_est[i] = sd(fuel$dat[[i]]$r)
    }
  }
  lambda_est = mean(lambda_est)
  mu_est = mean(mu_est)
  sigma_est = mean(sigma_est)
  rho_max = max(fuel$dimX,fuel$dimY)/3
  rho_est = rho_max/2
  
  # gamma parameters s.t. q95 ~~ rho_max
  #rho_prior = "gamma"
  rho_prior = "uniform"
  rho_shape = NULL # 1
  rho_rate = NULL # .3

  # set bounds for each parameters
  
  # the lower bound on lamda is very important. Certain metrics calculations will fail if not enough trees are sampled.
  # The error we get with small lamda is "Point cloud must have at least 2 points and at least 2 dimensions." I need to find out
  # which metric calculation is causing this.
  # We choose to assume that the true lambda such that on average produce half the number of trees as observed are produced
  # i.e. lambda >= lambda_est/2. And similarly we assume that lambda <= lambda_est*2
  prior.lb = c(0,       0,   0, lambda_est/2)
  prior.ub = c(rho_max, 3, Inf, lambda_est*2)
  
  Sigma = cov(metrics)
  Sinv = solve(Sigma + 1e-8*diag(dim(Sigma)[1]))
  ldetS = determinant(Sigma)$modulus
  
  # save everything needed to calculate priors
  prior_params = list()
  prior_params$lb = prior.lb; prior_params$ub = prior.ub 
  prior_params$rho_prior = rho_prior
  prior_params$rho_max = rho_max
  prior_params$rho_shape = rho_shape
  prior_params$rho_rate = rho_rate
  prior_params$mu_prior = "trunc normal"
  prior_params$mu_mean = 1.5
  prior_params$mu_sd = .5
  prior_params$prec_prior = 'gamma'
  prior_params$prec_shape = 1
  prior_params$prec_rate = 1e-3
  prior_params$lambda_prior = 'gamma'
  prior_params$lambda_shape = 1
  prior_params$lambda_rate = 1/lambda_est
  
  return(list('prior_params'=prior_params,
              'theta_est'=c(rho_est,mu_est,sigma_est,lambda_est),
              'Sinv'=Sinv,
              'ldetS'=ldetS))
}

get_mets_info = function(ncc=T,nholes=T,grid.area=T,moran=F,geary=T,rook=F,queen=T,perim=F,transect=F,n.trans=0)
{
  return(list('ncc'=ncc,
              'nholes'=nholes,
              'grid.area'=grid.area,
              'moran'=moran,
              'geary'=geary,
              'rook'=rook,
              'queen'=queen,
              'perim'=perim,
              'transect'=transect,
              'n.trans'=n.trans))
}

#' @title MCMC with joint proposal
#'
#' @description Addaptive MCMC as defined in Haario et al. 2001 - "An adaptive Metropolis algorithm"
#' @param data.dir directory with prior information file
#' @param n.samples total number of mcmc samples
#' @param n.burn number of samples to burn and use for addapting proposal covariance
#' @param update.every how often to update proposal covariance during burn in phase
#' @param load.theta load and add to previous mcmc samples
#' @param verbose frequency with which status updates are printed
#' @details Returns predictions at X.pred.orig
#' @export
#' @examples
#' # See examples folder for R markdown notebooks.
#'
mcmc = function(y_obs,fuel,prior,data.dir,n.samples=10000,n.burn=1000,
                gen_reps = 1, avg = 'llh', gen_parallel = F, mets_parallel = T,
                update.every=1,load.theta=0,make.cluster=T,parallel=T,verbose=0){
  cores=min(100,parallel::detectCores())
  if(make.cluster){
    cl = parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  }
  dimX = fuel$dimX
  dimY = fuel$dimY
  
  #prior.file = paste0(data.dir,'prior.RData')
  save.file = paste0(data.dir, 'mcmc.RData')
  #cat('prior file:',prior.file,'\n')
  #load(prior.file)
  cat('save file:',save.file,'\n')
  
  cat('n.samples:',n.samples,'\n')
  cat('n.burn:',n.burn,'\n')
  cat('load.theta:',load.theta,'\n')
  
  p.t = 4
  if(load.theta==0){
    theta = matrix(0,nrow=n.samples, ncol=p.t)
    accept = numeric(n.samples)
    theta.curr = prior$theta_est
    ll = llh(y_obs, theta.curr, dimX, dimY, prior$ldetS, prior$Sinv, mets_info, gen_reps, avg, gen_parallel, mets_parallel) + lprior(theta.curr,prior$prior_params)
    start = 1
  } else{
    cat('loading previous samples.\n')
    load(save.file)
    start = nrow(theta) + 1
    accept = c(accept,numeric(n.samples))
    theta.curr = theta[nrow(theta),]
    theta = rbind(theta,matrix(0,nrow=n.samples, ncol=p.t))
    ll = llh(y_obs, theta.curr, dimX, dimY, prior$ldetS, prior$Sinv, mets_info, gen_reps, avg, gen_parallel, mets_parallel) + lprior(theta.curr,prior$prior_params)
    theta = rbind(theta,matrix(0,nrow=n.samples, ncol=p.t))
  }
  
  prop.sd = .1*theta.curr # initial proposal sd is 10% of the data value
  if(!load.theta){
    prop.cov = diag(prop.sd^2)
  } else{
    prop.cov = cov(theta[1:(start-1),])
  }
  s.d = 2.4^2/p.t
  eps = sqrt(.Machine$double.eps)
  
  cat('parameter bounds:\n')
  cat('rho:',prior$prior_params$lb[1],prior$prior_params$ub[1],'\n')
  cat('mu:',prior$prior_params$lb[2],prior$prior_params$ub[2],'\n')
  cat('sigma:',prior$prior_params$lb[3],prior$prior_params$ub[3],'\n')
  cat('lambda:',prior$prior_params$lb[4],prior$prior_params$ub[4],'\n')
  cat('MCMC Start #--',format(Sys.time(), "%a %b %d %X"),'--#\n')
  
  for (i in start:(n.samples+start-1)){
    if(verbose>0){
      if(i %% verbose == 0){
        cat('MCMC iteration',i,'#--',format(Sys.time(), "%a %b %d %X"),'--#\n')
      }
    }
    prop = mvtnorm::rmvnorm(1,theta.curr,prop.cov)
    if(any(prop<prior$prior_params$lb) | any(prop>prior$prior_params$ub)){
      ll_prop = -Inf
    } else{
      ll_prop = llh(y_obs, prop, dimX, dimY, prior$ldetS, prior$Sinv, mets_info, gen_reps, avg, gen_parallel, mets_parallel) + lprior(prop,prior$prior_params)
    }
    ll = llh(y_obs, theta.curr, dimX, dimY, prior$ldetS, prior$Sinv, mets_info, gen_reps, avg, gen_parallel, mets_parallel) + lprior(theta.curr,prior$prior_params)
    
    # might help with mixing to recompute ll here since it's stochastic
    # This was necessary for my FlaGP stochastic likelihood mcmc
    # however, this might not be necessary because of the averaging
    if (runif(1) < exp(ll_prop - ll)){
      accept[i] = 1
      theta.curr = prop
      theta[i,] = prop
      ll = ll_prop
    } else {
      accept[i] = 0
      theta[i,] = theta.curr # reject -> reset theta
    }
    
    # for the first 100 iterations be greedy, updating the covriance matrix every time there is an accepted sample
    # and use only accepted samples for the estimate. We specify sum(accept>2) as the estimate of the covariance matrix
    # needs a few samples to work
    if(update.every>0 & i<100 & accept[i] & sum(accept>2)){
      # greedy in the beginning, update any time there is an accepted sample
      prop.cov = cov(theta[accept,]) + diag(eps,p.t)
    }
    if(update.every>0 & i<= n.burn & i %% update.every == 0){
      # if previous samples were passed to the function, use them for estimation
      prop.cov = cov(theta[1:i,]) + diag(eps,p.t)
      #C.store[as.integer(i/update.every),,] = C
    }
    
    if(i%%1000==0){
      theta.save = theta[1:i,]
      save(prop.cov,theta.save,file=save.file)
    }
  }
  
  cat('saving mcmc results to',save.file,'\n')
  old.samples = 1:(start-1)
  new.samples = (start+n.burn):(n.samples+start-1)
  theta = theta[c(old.samples,new.samples),]
  accept = accept[c(old.samples,new.samples)]
  save(prop.cov,theta,accept,file=save.file)
  if(make.cluster)
    parallel::stopCluster(cl)
  cat('done.')
  return(list('theta'=theta,
              'prop.cov'=prop.cov,
              'accept'=accept))
}

