library(patchwork)
library(ggplot2)
library(dplyr)

diagnose_margs_gpd <- function(margs, dfin, var_str="hs", exc_str='exc'){
  #' @export

  n.bstrp <- length(margs)
  dfin <- dfin
  gpd_sims <- array(1, c(n.bstrp, dim(dfin)[1]))*NA
  gpd_sims_exc <- array(1, c(n.bstrp, dim(dfin)[1]))*NA
  for (b in 1:n.bstrp){
    print(c("bootstrap nr:", b))
    m_gpd <- margs[[b]][['gpd']][[var_str]]
    m_ald <- margs[[b]][['thr']][[var_str]]
    # predict threshold
    thr <- predict(m_ald, newdata = dfin, type = "response")$location
    # predict exceedance
    gpd_param_sims <- predict(m_gpd, newdata = dfin, type= "response")
    scales <- gpd_param_sims$scale
    shapes <- gpd_param_sims$shape
    gpd_sims[b,] <- revd(length(scales), scale = scales, shape = shapes,
                         threshold = thr, type="GP")
    gpd_sims_exc[b,] <- revd(length(scales), scale = scales, shape = shapes,
                             threshold = 0, type="GP")
  }

  qs <- c(seq(.1,.99,.05),.995, .999)
  qs_gpd <- apply(gpd_sims,1,quantile,qs)
  median_gpd <- apply(qs_gpd,1,quantile,.5)
  mean_gpd <- apply(qs_gpd,1,mean)
  q99_gpd <- apply(qs_gpd,1,quantile,.99)
  q01_gpd <- apply(qs_gpd,1,quantile,.01)
  qs_obs <- quantile(dfin[[var_str]], qs)

  qs_gpd_exc <- apply(gpd_sims_exc,1,quantile,qs)
  median_gpd_exc <- apply(qs_gpd_exc,1,quantile,.5)
  mean_gpd_exc <- apply(qs_gpd_exc,1,mean)
  q99_gpd_exc <- apply(qs_gpd_exc,1,quantile,.99)
  q01_gpd_exc <- apply(qs_gpd_exc,1,quantile,.01)
  qs_obs_exc <- quantile(dfin[[exc_str]], qs)

  par(mfrow = c(1, 2))

  xlim <- c(0,max((q99_gpd) + 0.1*max(q99_gpd)))
  ylim <- xlim
  plot(qs_obs, median_gpd, xlim=xlim, ylim=ylim, xlab = '', ylab = '', xaxt='n', yaxt='n', main = "")
  par(new=TRUE)
  plot(qs_obs, q99_gpd, xlim=xlim, ylim=ylim, col='red', type='l',lwd = .5, xlab = '', ylab = '', xaxt='n', yaxt='n', main = "")
  par(new=TRUE)
  plot(qs_obs, q01_gpd, xlim=xlim, ylim=ylim, col='red', type='l',lwd = .5, xlab = 'observation quantiles', ylab = 'model quantiles', main='QQ-plot')
  # Adding a diagonal line
  abline(a = 0, b = 1, col = "gray", lwd = 2)

  xlim <- c(0,(max(q99_gpd_exc) + 0.1*max(q99_gpd_exc)))
  ylim <- xlim
  plot(qs_obs_exc, median_gpd_exc, xlim=xlim, ylim=ylim, xlab = '', ylab = '', xaxt='n', yaxt='n', main = "")
  par(new=TRUE)
  plot(qs_obs_exc, q99_gpd_exc, xlim=xlim, ylim=ylim, col='red', type='l',lwd = .5, xlab = '', ylab = '', xaxt='n', yaxt='n', main = "")
  par(new=TRUE)
  plot(qs_obs_exc, q01_gpd_exc, xlim=xlim, ylim=ylim, col='red', type='l',lwd = .5, xlab = 'observation quantiles', ylab = 'model quantiles', main='QQ-plot exceedances')
  # Adding a diagonal line
  abline(a = 0, b = 1, col = "gray", lwd = 2)
}

diagnose_margs_occ_rejection <- function(
                               marg, dfin, nbins, covarlst,
                               nr_of_years,
                               covar_mins,
                               covar_maxes,
                               RP=NULL, cslims=NULL,
                               condition=NULL,
                               grid_interval=NULL,
                               plot_ls=TRUE)
  {
  #' @export

  if (is.null(RP)){
    RP <- nr_of_years
  }
  nr_of_events <- dim(dfin)[1]
  dfobs <- dfin[,covarlst]
  dfcounts <- produce_storm_occurrences_rejection(nr_of_events, RP,
                                        nr_of_years,
                                        marg,
                                        covarlst,
                                        covar_mins, covar_maxes,
                                        dfin = dfin,
                                        condition = condition,
                                        grid_interval = grid_interval)

  # 2D plots
  if (plot_ls == TRUE){
    p1 <- ggplot(dfobs, aes(dfobs[[covarlst[1]]], dfobs[[covarlst[2]]])) +
    geom_bin2d(bins = nbins) +  # bins controls the number of bins
    scale_fill_gradient(low = "white", high = "blue", limits = cslims) +
    labs(title = "2D Histogram (data)",
         x = covarlst[[1]], y = covarlst[[2]])

    p2 <- ggplot(dfcounts, aes(dfcounts[[covarlst[1]]], dfcounts[[covarlst[2]]])) +
    geom_bin2d(bins = nbins) +  # bins controls the number of bins
    scale_fill_gradient(low = "white", high = "blue", limits = cslims) +
    labs(title = "2D Histogram (predicted)",
         x = covarlst[[1]], y = covarlst[[2]])
    p1 + p2
    return(list(p1 + p2, dim(dfcounts)[1]))
  } else {return(list(dim(dfcounts)[1]))}
}

diagnose_margs_occ_rejection_bstrp_1D <- function(margs, varstr, nbins, covarlst,
                                                  nr_of_years,
                                                  covar_mins,
                                                  covar_maxes,
                                                  RP=NULL, cslims=NULL,
                                                  grid_interval=NULL){
  #' export
  #'

  dimlst_rej <- NULL
  nbstrp <- length(margs)
  for (b in 1:nbstrp){
    print(b)
    tmpres <- diagnose_margs_occ_rejection(margs[[b]]$occ[[varstr]],
                                           margs[[b]]$gpd[[varstr]]$data,
                                           nbins,covarlst,nr_of_years,
                                           covar_mins = covar_mins,
                                           covar_maxes = covar_maxes,
                                           cslims = cslims, RP = RP,
                                           grid_interval = grid_interval)
    dimlst_rej[[b]] <- tmpres[[2]]
  }

  occ_from_data <- NULL
  for(b in 1:nbstrp){occ_from_data[[b]] <- dim(margs[[b]]$gpd[[varstr]]$data)[1]}

  par(mfrow = c(1, 1))
  ulim_y <- max(density(unlist(occ_from_data))$y)
  plot(density(unlist(dimlst_rej)),
       main = "", xlab = "", ylab = "", xaxt='n', yaxt='n',
       lty=1, col="brown", ylim=c(0,ulim_y))
  par(new=TRUE)
  plot(density(unlist(occ_from_data)),
       main = "", xlab = "", ylab = "",
       col="orange", lty=1, ylim=c(0,ulim_y))
  abline(v = mean(unlist(occ_from_data)), col = "orange", lwd = 1, lty = 2)
  abline(v = mean(unlist(dimlst_rej)), col = "brown", lwd = 1, lty = 1)
  legend("topright",
         legend = c("occurrences simulated","occurrences in data"),
         col = c("brown", "orange"),
         lty = c(1,3), # line type
         lwd = 1,        # line width
         bty = "n")      # box type ("n" = no box)
}

diagnose_margs_occ_rejection_bstrp_2D <- function(margs, varstr, nbins, covarlst,
                                                  nr_of_years,
                                                  covar_mins,
                                                  covar_maxes,
                                                  RP=NULL, cslims=NULL,
                                                  condition=NULL,
                                                  grid_interval=NULL,
                                                  plot_ls=TRUE){
  #' @export

  if (is.null(RP)){
    RP <- nr_of_years
  }
  dfcountslst <- NULL
  dfobslst <- NULL
  nbstrp <- length(margs)
  for (b in 1:nbstrp){
    dfin <- margs[[b]]$gpd[[varstr]]$data
    nr_of_events <- dim(dfin)[1]
    dfobslst[[b]] <- dfin[,covarlst]
    dfcounts_tmp <- produce_storm_occurrences_rejection(nr_of_events, RP,
                                                    nr_of_years,
                                                    margs[[b]]$occ[[varstr]],
                                                    covarlst,
                                                    covar_mins, covar_maxes,
                                                    dfin = dfin,
                                                    condition = condition,
                                                    grid_interval = grid_interval)
    dfcountslst[[b]] <- dfcounts_tmp
  }
  dfobs <- bind_rows(dfobslst)
  dfcounts <- bind_rows(dfcountslst)

  if (plot_ls == TRUE){
    p1 <- ggplot(dfobs, aes(dfobs[[covarlst[1]]], dfobs[[covarlst[2]]])) +
      geom_bin2d(bins = nbins) +  # bins controls the number of bins
      scale_fill_gradient(low = "white", high = "blue", limits = cslims) +
      labs(title = "2D Histogram (data)",
           x = covarlst[[1]], y = covarlst[[2]])

    p2 <- ggplot(dfcounts, aes(dfcounts[[covarlst[1]]], dfcounts[[covarlst[2]]])) +
      geom_bin2d(bins = nbins) +  # bins controls the number of bins
      scale_fill_gradient(low = "white", high = "blue", limits = cslims) +
      labs(title = "2D Histogram (predicted)",
           x = covarlst[[1]], y = covarlst[[2]])
    p1 + p2
    return(list(p1 + p2, dim(dfcounts)[1]))
  } else {return(list(dim(dfcounts)[1]))}
}

diagnose_margs_occ <- function(marg, dfin, nbins, covarlst,
                               nr_of_years,
                               covar_mins,
                               covar_maxes,
                               RP=NULL, cslims=NULL,
                               condition=NULL){
  #' @export

  if (is.null(RP)){
    RP <- nr_of_years
  }
  nr_of_events <- dim(dfin)[1]
  dfobs <- dfin[,covarlst]
  dfcounts <- produce_storm_occurrences(nr_of_events, RP,
                                        nr_of_years,
                                        marg,
                                        covarlst,
                                        covar_mins, covar_maxes,
                                        condition = condition)

  p1 <- ggplot(dfobs, aes(dfobs[[covarlst[1]]], dfobs[[covarlst[2]]])) +
    geom_bin2d(bins = nbins) +  # bins controls the number of bins
    scale_fill_gradient(low = "white", high = "blue", limits = cslims) +
    labs(title = "2D Histogram (data)",
         x = covarlst[[1]], y = covarlst[[2]])

  p2 <- ggplot(dfcounts, aes(dfcounts[[covarlst[1]]], dfcounts[[covarlst[2]]])) +
    geom_bin2d(bins = nbins) +  # bins controls the number of bins
    scale_fill_gradient(low = "white", high = "blue", limits = cslims) +
    labs(title = "2D Histogram (predicted)",
         x = covarlst[[1]], y = covarlst[[2]])

  p1 + p2

  return(list(p1 + p2, dim(dfcounts)[1]))
}

diagnose_maxds <- function(models_maxds_lst, X_var_str, Y_var_str, ulim=1.3, llim=-1.3){
  #' @export

  maxd_params_bstrp <- unfold_maxd_params_bstrp_v2(models_maxds_lst, X_var_str, Y_var_str)

  par(mfrow = c(2, 2))
  plot(maxd_thr, maxd_params_bstrp$alpha_cntr, type = "l",ylim=c(llim,ulim), main = "", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
  par(new=TRUE)
  plot(maxd_thr, maxd_params_bstrp$alpha_llim, type = "l",ylim=c(llim,ulim), col="red",  main = "", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
  par(new=TRUE)
  plot(maxd_thr, maxd_params_bstrp$alpha_ulim, type = "l",ylim=c(llim,ulim), col="red",  main = "", xlab = "maxd threshold quantile", ylab = "alpha")
  abline(h = 0, col = "gray", lwd = .8, lty = 1)
  abline(h = 1, col = "gray", lwd = .8, lty = 1)
  abline(h = -1, col = "gray", lwd = .8, lty = 1)
  abline(h = .5, col = "gray", lwd = .8, lty = 3)
  abline(h = -.5, col = "gray", lwd = .8, lty = 3)

  plot(maxd_thr, maxd_params_bstrp$beta_cntr, type = "l",ylim=c(llim,ulim), main = "", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
  par(new=TRUE)
  plot(maxd_thr, maxd_params_bstrp$beta_llim, type = "l",ylim=c(llim,ulim), col="red", main = "", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
  par(new=TRUE)
  plot(maxd_thr, maxd_params_bstrp$beta_ulim, type = "l",ylim=c(llim,ulim), col="red", main = "", xlab = "maxd threshold quantile", ylab = "beta")
  abline(h = 0, col = "gray", lwd = .8, lty = 1)
  abline(h = 1, col = "gray", lwd = .8, lty = 1)
  abline(h = -1, col = "gray", lwd = .8, lty = 1)
  abline(h = .5, col = "gray", lwd = .8, lty = 3)
  abline(h = -.5, col = "gray", lwd = .8, lty = 3)

  plot(maxd_thr, maxd_params_bstrp$mu_cntr, type = "l",ylim=c(llim,ulim), main = "", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
  par(new=TRUE)
  plot(maxd_thr, maxd_params_bstrp$mu_llim, type = "l",ylim=c(llim,ulim), col="red", main = "", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
  par(new=TRUE)
  plot(maxd_thr, maxd_params_bstrp$mu_ulim, type = "l",ylim=c(llim,ulim), col="red", main = "", xlab = "maxd threshold quantile", ylab = "mu")
  abline(h = 0, col = "gray", lwd = .8, lty = 1)
  abline(h = 1, col = "gray", lwd = .8, lty = 1)
  abline(h = -1, col = "gray", lwd = .8, lty = 1)
  abline(h = .5, col = "gray", lwd = .8, lty = 3)
  abline(h = -.5, col = "gray", lwd = .8, lty = 3)

  plot(maxd_thr, maxd_params_bstrp$sigma_cntr, type = "l",ylim=c(llim,ulim), main = "", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
  par(new=TRUE)
  plot(maxd_thr, maxd_params_bstrp$sigma_llim, type = "l",ylim=c(llim,ulim), col="red", main = "", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
  par(new=TRUE)
  plot(maxd_thr, maxd_params_bstrp$sigma_ulim, type = "l",ylim=c(llim,ulim), col="red", main = "", xlab = "maxd threshold quantile", ylab = "sigma")
  abline(h = 0, col = "gray", lwd = .8, lty = 1)
  abline(h = 1, col = "gray", lwd = .8, lty = 1)
  abline(h = -1, col = "gray", lwd = .8, lty = 1)
  abline(h = .5, col = "gray", lwd = .8, lty = 3)
  abline(h = -.5, col = "gray", lwd = .8, lty = 3)
}

diagnose_maxds_fitted <- function(maxd, X_str, Y_str, nsim){
  #' @export

  # plot Laplace margins
  X_all <- maxd[[X_str]][[Y_str]]$X_all
  X_fit <- maxd[[X_str]][[Y_str]]$X_fit
  Y_all <- maxd[[X_str]][[Y_str]]$Y_all
  Y_fit <- maxd[[X_str]][[Y_str]]$Y_fit
  xlim <- c(min(c(X_all, Y_all)),
            max(c(X_all, Y_all)))
  ylim <- xlim
  plot(X_all, Y_all,
       pch=20, col="black", xlim=xlim, ylim=ylim, cex=.5,
       main = "", xlab = "X on Laplace", ylab = "Y on Laplace")
  par(new=TRUE)
  plot(X_fit, Y_fit,
       pch=20, col="orange", xlim=xlim, ylim=ylim, cex=.5,
       main = "", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
  abline(h = 0, col = "gray", lwd = .5, lty = 1)
  abline(v = 0, col = "gray", lwd = .5, lty = 1)
  abline(v = quantile(X_all, maxd[[X_str]][[Y_str]]$thr),
         col = "orange", lwd = .5, lty = 1)

  # plot simulations from HT2004
  sims <- predict_maxd(maxd, X_str, Y_str, nsim)
  lp_X <- sims$Xsim
  lp_Y <- sims$Ysim
  xlim <- c(min(c(X_all,unlist(sims))),max(c(X_all,unlist(sims))))
  ylim <- xlim
  plot(X_all, Y_all,
       pch=20, col="black", xlim=xlim, ylim=ylim, cex=.5,
       main = "", xlab = "X on Laplace", ylab = "Y on Laplace")
  par(new=TRUE)
  plot(lp_X, lp_Y,
       pch=20, col="orange", xlim=xlim, ylim=ylim, cex=.5,
       main = "", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
  abline(h = 0, col = "gray", lwd = .5, lty = 1)
  abline(v = 0, col = "gray", lwd = .5, lty = 1)
  abline(v = quantile(X_all, maxd[[X_str]][[Y_str]]$thr),
         col = "orange", lwd = .5, lty = 1)

  # plot division when sim X or sim Y is largest (diagonal split)
  plot(X_all, Y_all,
       pch=20, col="black", xlim=xlim, ylim=ylim, cex=.5,
       main = "", xlab = "X on Laplace", ylab = "Y on Laplace")
  for (i in 1:length(lp_X)){
    par(new=TRUE)
    if (lp_X[i] > lp_Y[i]){
      plot(lp_X[i], lp_Y[i],
           pch=20, col="red", xlim=xlim, ylim=ylim, cex=.5,
           main = "", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
    } else {
      plot(lp_X[i], lp_Y[i],
           pch=20, col="cyan", xlim=xlim, ylim=ylim, cex=.5,
           main = "", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
    }
  }
  abline(h = 0, col = "gray", lwd = .5, lty = 1)
  abline(v = 0, col = "gray", lwd = .5, lty = 1)
  abline(v = quantile(X_all, maxd[[X_str]][[Y_str]]$thr),
         col = "orange", lwd = .5, lty = 1)
  abline(coef = c(0,1), col="orange")

  # plot division when sim X or sim Y is largest (diagonal split)
  # and simulated values from respective model
  sims <- predict_maxd(maxd, X_str, Y_str, nsim)
  lp_X <- sims$Xsim
  lp_Y <- sims$Ysim

  plot(X_all, Y_all,
       pch=20, col="black", xlim=xlim, ylim=ylim, cex=.5,
       main = "", xlab = "X on Laplace", ylab = "Y on Laplace")
  for (i in 1:length(lp_X)){
    par(new=TRUE)
    if (lp_X[i] > lp_Y[i]){
      plot(lp_X[i], lp_Y[i],
           pch=20, col="red", xlim=xlim, ylim=ylim, cex=.5,
           main = "", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
    }
  }
  sims <- predict_maxd(maxd, Y_str, X_str, nsim)
  lp_X <- sims$Xsim
  lp_Y <- sims$Ysim
  for (i in 1:length(lp_X)){
    par(new=TRUE)
    if (lp_X[i] > lp_Y[i]){
      plot(lp_Y[i], lp_X[i],
           pch=20, col="cyan", xlim=xlim, ylim=ylim, cex=.5,
           main = "", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
    }
  }

  abline(h = 0, col = "gray", lwd = .5, lty = 1)
  abline(v = 0, col = "gray", lwd = .5, lty = 1)

  abline(v = quantile(X_all, maxd[[X_str]][[Y_str]]$thr),
         col = "orange", lwd = .5, lty = 1)
  abline(h = quantile(X_all, maxd[[Y_str]][[X_str]]$thr),
         col = "orange", lwd = .5, lty = 1)
  abline(coef = c(0,1), col="orange")
}

diagnose_margs_preds_density <- function(margs_preds, varstr,
                                         n.bstr=NULL, bw=NULL,
                                         xlim=NULL, ylim=NULL,
                                         qlow=.05, qhigh=.95, qmed=.5,
                                         xlab=NULL){
  #' @export

  # if is.null(n.bstr) plot all combined
  if (is.null(n.bstr)){
    n.bstr <- seq(length(margs_preds))
  }

  maxvals <- NULL
  for (b in n.bstr){
    tmp <- margs_preds[[b]][[varstr]]$maxval
    maxvals[[b]] <- tmp[is.finite(tmp)]
  }

  if (is.null(xlim)){
    xlim <- c(min(unlist(maxvals)), max(unlist(maxvals)))
  }

  y_max <- array(dim=length(maxvals))
  for (b in 1:length(y_max)){
    y_max[b] <- max(density(maxvals[[b]],bw=bw,n=512)$y)
  }

  if (is.null(ylim)){
    ylim <- c(0,max(y_max))
  }

  if (is.null(bw)){
    bw <- .01*abs(xlim[1]-xlim[2])
  }

  # plot combined density with individual densities
  for (b in n.bstr){
    plot(density(maxvals[[b]], bw=bw, n=512,
                 from=xlim[1], to=xlim[2]),
         col="grey", lwd=.5, xlim=xlim, ylim=ylim,
         xlab = '', ylab = '', xaxt='n', yaxt='n', main = "")
    par(new=TRUE)
  }
  if (is.null(xlab)){
    xlab <- varstr
  }
  plot(density(unlist(maxvals), bw=bw, n=512,
               from=xlim[1], to=xlim[2]),
       xlim=xlim, ylim=ylim, main="", xlab=xlab)

  # plot density with uncertainty quantiles
  densities <- NULL
  for (b in n.bstr){
    densities[['d']][[b]] <- density(maxvals[[b]], bw=bw, n=512,
                                     from=xlim[1], to=xlim[2])
    densities[['x']][[b]] <- densities[['d']][[b]]$x
    densities[['y']][[b]] <- densities[['d']][[b]]$y
  }
  df_dens <- data.frame(densities$y)

  qlow_ts <- apply(df_dens,1,quantile,qlow)
  qhigh_ts <- apply(df_dens,1,quantile,qhigh)
  qmed_ts <- apply(df_dens,1,quantile,qmed)

  plot(densities$x[[1]], qlow_ts,col='blue', xlim=xlim, ylim=ylim, main="",
       lty=3, lwd=.5, type='l', xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
  par(new=TRUE)
  plot(densities$x[[1]], qhigh_ts,col='blue',xlim=xlim, ylim=ylim, main="",
       lty=3, lwd=.5, type='l', xlab = "", ylab = "", xaxt = 'n', yaxt = 'n',)
  par(new=TRUE)
  plot(densities$x[[1]], qmed_ts,xlim=xlim, ylim=ylim, main="",
       xlab=xlab, ylab="Density", lty=3, lwd=.5, type='l', col="blue")
  par(new=TRUE)
  plot(density(unlist(maxvals), bw=bw, n=512,
               from=xlim[1], to=xlim[2]),
       xlim=xlim, ylim=ylim, main="", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')

  print("summary maxes:")
  print(summary(unlist(maxvals)))
  return(list('maxvals'=maxvals, 'densities'=densities))
}

display_joint_densities <- function(preds_maxds,
                                    preds_maxds_LP,
                                    preds_margs_lst,
                                    X_str, Y_str,
                                    Xmax=NULL,Ymax=NULL){
  # use function to retrieve valids
  valids <- retrieve_valid_HT_samples(preds_maxds,
                                      preds_maxds_LP,
                                      preds_margs_lst,
                                      X_str, c(Y_str))

  library(hexbin)
  library(RColorBrewer)
  library(lattice)
  rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  r <- rf(32)
  # Create hexbin object and plot
  dfin <- data.frame(unlist(valids$X_valid[[Y_str]]), unlist(valids$Y_valid[[Y_str]]))
  colnames(dfin) <- c(X_str, Y_str)
  h <- hexbin(dfin)
  counts <- h@count
  formula <- as.formula(paste(Y_str, "~", X_str))
  # Define a custom panel function to add a vertical line
  custom_panel <- function(x, y, Xmax,Ymax, ...) {
    panel.hexbinplot(x, y, ...)
    panel.abline(v = Xmax, col = "gray", lwd = 1, lty = 1)
    panel.abline(h = Ymax, col = "gray", lwd = 1, lty = 1)
  }
  hexbinplot(formula, data=dfin,
             xbins=50, colramp=rf,
             mincnt=(min(counts)+1),
             maxcnt = max(counts),
             trans=log, inv=exp,
             panel = custom_panel,
             Xmax = Xmax, Ymax=Ymax)
             #ylim=c(6,12), xlim=c(5,10))
  # ylim=c(0,2.5), xlim=c(10,24)
  #if (is.numeric(Xmax)){
    #panel.abline(v = Xmax, col = "gray", lwd = .5, lty = 1)
  #}
}

plot_cvres <- function(cvres, limits=NULL, show_errors=NULL){
  #' Plots the results from the cross-validation procedure
  #' @export
  #'
  # print values of minimum
  print(c(cvres$thr[cvres$cost$costfct1==min(cvres$cost$costfct1)],
          cvres$thr[cvres$cost$costfct2==min(cvres$cost$costfct2)],
          cvres$thr[cvres$cost$costfct3==min(cvres$cost$costfct3)]))

  if (is.null(limits)){
    # plot mean cost-function
    par(mfrow = c(1,3))
    plot(cvres$thr, cvres$cost$costfct1)
    abline(v=cvres$thr[cvres$cost$costfct1==min(cvres$cost$costfct1)])
    plot(cvres$thr, cvres$cost$costfct2)
    abline(v=cvres$thr[cvres$cost$costfct2==min(cvres$cost$costfct2)])
    plot(cvres$thr, cvres$cost$costfct3)
    abline(v=cvres$thr[cvres$cost$costfct3==min(cvres$cost$costfct3)])

  } else {
    # Add the shaded region
    # Define the polygon coordinates
    x1 = limits[1]
    x2 = limits[2]
    polygon_x <- c(x1, x2, x2, x1)
    polygon_y <- c(-1000, -1000, 1000, 1000) # Extend to the plot limits vertically
    # plot mean cost-function
    par(mfrow = c(1,3))
    plot(cvres$thr, cvres$cost$costfct1)
    polygon(polygon_x, polygon_y, col = rgb(0.6, 0.6, 0.6, 0.3), border = NA)
    abline(v=cvres$thr[cvres$cost$costfct1==min(cvres$cost$costfct1)])
    plot(cvres$thr, cvres$cost$costfct2)
    polygon(polygon_x, polygon_y, col = rgb(0.6, 0.6, 0.6, 0.3), border = NA)
    abline(v=cvres$thr[cvres$cost$costfct2==min(cvres$cost$costfct2)])
    plot(cvres$thr, cvres$cost$costfct3)
    polygon(polygon_x, polygon_y, col = rgb(0.6, 0.6, 0.6, 0.3), border = NA)
    abline(v=cvres$thr[cvres$cost$costfct3==min(cvres$cost$costfct3)])
  }

  if (!is.null(show_errors)){
    # plot errors
    par(mfrow = c(2,2))
    boxplot(t(cvres$errors$bias),names=cvres$thr, main="BIAS")
    abline(h=0)
    boxplot(t(cvres$errors$mae),names=cvres$thr, main="MAE")
    abline(h=0)
    boxplot(t(cvres$errors$rse),names=cvres$thr, main="RSE")
    abline(h=0)
    boxplot(t(cvres$errors$drse),names=cvres$thr, main="dRSE")
    abline(h=0)
  }
}

visualize_storm_picking <- function(dfin, dfinall=NULL, xstr = 'dt', ystr='hs', sidx=NULL, eidx=NULL){
  #' function to visualize the outcome of the storm picking procedure
  #'
  #' @export
  #'

  if (is.null(sidx)){sidx <- 1}
  if (is.null(eidx)){eidx <- length(dfin[[ystr]])}

  if (!is.null(dfinall)){
    plot(dfinall[[xstr]][sidx:(eidx*3)], dfinall[[ystr]][sidx:(eidx*3)],
         pch=20, cex=.5,
         col = adjustcolor("gray", alpha.f = 0.5),
         xlim=c(dfin$storms[[xstr]][sidx],dfin$storms[[xstr]][eidx]), ylim=c(0,15),
         xlab = '', ylab = '', xaxt='n', yaxt='n', main = "")  # e.g. all hs from data
    par(new=TRUE)
  }
  plot(dfin$storms[[xstr]][sidx:(eidx*2)], dfin$storms$thr[sidx:(eidx*2)],
       pch=20, cex=.5,
       col = adjustcolor("red", alpha.f = 0.5),
       xlim=c(dfin$storms[[xstr]][sidx],dfin$storms[[xstr]][eidx]), ylim=c(0,15),
       xlab = '', ylab = '', xaxt='n', yaxt='n', main = "")  # non-stationary local threshold
  par(new=TRUE)
  plot(dfin$storms[[xstr]][sidx:(eidx*2)], dfin$storms$exc[sidx:(eidx*2)]+df_pick$storms$thr[sidx:(eidx*2)],
       pch=20, cex=.5,
       col = adjustcolor("black", alpha.f = 0.5),
       xlim=c(dfin$storms[[xstr]][sidx],dfin$storms[[xstr]][eidx]), ylim=c(0,15),
       xlab = '', ylab = '', xaxt='n', yaxt='n', main = "")  # exceedences
  par(new=TRUE)
  plot(dfin$pots[[xstr]][sidx:eidx], dfin$pots[[ystr]][sidx:eidx],
       pch=1, cex=1.,
       col = adjustcolor("orange", alpha.f = 0.5),
       xlim=c(dfin$storms[[xstr]][sidx],dfin$storms[[xstr]][eidx]), ylim=c(0,15),
       xlab = '', ylab = 'Hs [m]', main = "")  # peaks
}
