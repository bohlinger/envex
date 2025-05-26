library(patchwork)
library(ggplot2)

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

diagnose_margs_occ <- function(marg, dfin, nbins, covarlst, nr_of_years, cslims=NULL){
  #' @export

  nr_of_events <- dim(dfin)[1]
  dfobs <- dfin[,covarlst]
  dfcounts <- produce_storm_occurrences(nr_of_events,
                                        nr_of_years, nr_of_years,
                                        marg)

  print(c('max predicted counts:', max(dfcounts$counts)))

  if (is.null(cslims)){
    cslims <- c(0,max(dfcounts$counts))
  }

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
                                         qlow=.05, qhigh=.95, qmed=.5){
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
    y_max[b] <- max(density(maxvals[[b]])$y)
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
  plot(density(unlist(maxvals), bw=bw, n=512,
               from=xlim[1], to=xlim[2]),
       xlim=xlim, ylim=ylim, main="", xlab=varstr)

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

  plot(densities$x[[1]], qlow_ts,col='grey', xlim=xlim, ylim=ylim, main="",
       lty=1, lwd=.5, type='l', xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
  par(new=TRUE)
  plot(densities$x[[1]], qhigh_ts,col='grey',xlim=xlim, ylim=ylim, main="",
       lty=1, lwd=.5, type='l', xlab = "", ylab = "", xaxt = 'n', yaxt = 'n',)
  par(new=TRUE)
  plot(densities$x[[1]], qmed_ts,xlim=xlim, ylim=ylim, main="",
       xlab=varstr, ylab="Density", lty=1, lwd=.5, type='l')

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
  # ylim=c(0,2.5), xlim=c(10,24)
  #if (is.numeric(Xmax)){
    #panel.abline(v = Xmax, col = "gray", lwd = .5, lty = 1)
  #}
}
