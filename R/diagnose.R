library(patchwork)
library(ggplot2)
library(dplyr)

diagnose_margs_gpd <- function(margs, dfin, var_str = "hs", exc_str = "exc") {
  #' @export

  n.bstrp <- length(margs)
  dfin <- dfin
  gpd_sims <- array(1, c(n.bstrp, dim(dfin)[1])) * NA
  gpd_sims_exc <- array(1, c(n.bstrp, dim(dfin)[1])) * NA
  for (b in 1:n.bstrp) {
    print(c("bootstrap nr:", b))
    m_gpd <- margs[[b]][["gpd"]][[var_str]]
    m_ald <- margs[[b]][["thr"]][[var_str]]
    # predict threshold
    thr <- predict(m_ald, newdata = dfin, type = "response")$location
    # predict exceedance
    gpd_param_sims <- predict(m_gpd, newdata = dfin, type = "response")
    scales <- gpd_param_sims$scale
    shapes <- gpd_param_sims$shape
    gpd_sims_exc[b, ] <- revd(length(scales), scale = scales, shape = shapes,
                             threshold = 0, type = "GP")
    gpd_sims[b, ] <- gpd_sims_exc[b,] + thr
    # HERE
    #gpd_sims[b, ] <- revd(length(scales), scale = scales, shape = shapes,
    #                      threshold = thr, type = "GP")
    #gpd_sims_exc[b, ] <- revd(length(scales), scale = scales, shape = shapes,
    #                         threshold = 0, type = "GP")
  }

  qs <- c(seq(.1, .99, .05), .995, .999)
  qs_gpd <- apply(gpd_sims, 1, quantile, qs)
  median_gpd <- apply(qs_gpd, 1, quantile, .5)
  q99_gpd <- apply(qs_gpd, 1, quantile, .99)
  q01_gpd <- apply(qs_gpd, 1, quantile, .01)
  qs_obs <- quantile(dfin[[var_str]], qs)

  qs_gpd_exc <- apply(gpd_sims_exc, 1, quantile, qs)
  median_gpd_exc <- apply(qs_gpd_exc, 1, quantile, .5)
  q99_gpd_exc <- apply(qs_gpd_exc, 1, quantile, .99)
  q01_gpd_exc <- apply(qs_gpd_exc, 1, quantile, .01)
  qs_obs_exc <- quantile(dfin[[exc_str]], qs)

  par(mfrow = c(1, 2))

  xlim <- c(0, max((q99_gpd) + 0.1 * max(q99_gpd)))
  ylim <- xlim
  plot(qs_obs, median_gpd, xlim = xlim, ylim = ylim, xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = "")
  par(new = TRUE)
  plot(qs_obs, q99_gpd, xlim = xlim, ylim = ylim, col = "red", type = "l", lwd = .5, xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = "")
  par(new = TRUE)
  plot(qs_obs, q01_gpd, xlim = xlim, ylim = ylim, col = "red", type = "l", lwd = .5, xlab = "observation quantiles", ylab = "model quantiles", main = "QQ-plot")
  # Adding a diagonal line
  abline(a = 0, b = 1, col = "gray", lwd = 2)

  xlim <- c(0, (max(q99_gpd_exc) + 0.1 * max(q99_gpd_exc)))
  ylim <- xlim
  plot(qs_obs_exc, median_gpd_exc, xlim = xlim, ylim = ylim, xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = "")
  par(new = TRUE)
  plot(qs_obs_exc, q99_gpd_exc, xlim = xlim, ylim = ylim, col = "red", type = "l", lwd = .5, xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = "")
  par(new = TRUE)
  plot(qs_obs_exc, q01_gpd_exc, xlim = xlim, ylim = ylim, col = "red", type = "l", lwd = .5, xlab = "observation quantiles", ylab = "model quantiles", main = "QQ-plot exceedances")
  # Adding a diagonal line
  abline(a = 0, b = 1, col = "gray", lwd = 2)
}

diagnose_margs_occ_rejection_noplot <- function(
    marg, dfin, nbins, covarlst,
    nr_of_years,
    RP = NULL, cslims = NULL,
    condition = NULL,
    grid_interval = NULL,
    plot_ls = TRUE) {
  #' @export

  if (is.null(RP)) {
    RP <- nr_of_years
  }
  nr_of_events <- dim(dfin)[1]
  dfcounts <- produce_storm_occurrences_rejection(nr_of_events, RP,
                                                  nr_of_years,
                                                  marg,
                                                  covarlst,
                                                  dfin = dfin,
                                                  condition = condition,
                                                  grid_interval = grid_interval)
  return (dim(dfcounts)[1])
}

diagnose_margs_occ_rejection <- function(
                               marg, dfin, nbins, covarlst,
                               nr_of_years,
                               RP = NULL, cslims = NULL,
                               condition = NULL,
                               grid_interval = NULL,
                               plot_ls = TRUE)
  {
  #' @export

  if (is.null(RP)) {
    RP <- nr_of_years
  }
  nr_of_events <- dim(dfin)[1]
  dfobs <- dfin[, covarlst]
  dfcounts <- produce_storm_occurrences_rejection(nr_of_events, RP,
                                        nr_of_years,
                                        marg,
                                        covarlst,
                                        dfin = dfin,
                                        condition = condition,
                                        grid_interval = grid_interval)
  countnr <- dim(dfcounts)[1]
  # 2D plots
  if (plot_ls == TRUE) {
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
  } else {
    return(countnr)
    }
}

diagnose_margs_occ_rejection_bstrp_1D <- function(margs, varstr, nbins, covarlst,
                                                  nr_of_years,
                                                  RP = NULL, cslims = NULL,
                                                  grid_interval = NULL){
  #' export
  #'

  dimlst_rej <- NULL
  nbstrp <- length(margs)
  for (b in 1:nbstrp) {
    print(b)
    tmpres <- diagnose_margs_occ_rejection_noplot(margs[[b]]$occ[[varstr]],
                                           margs[[b]]$gpd[[varstr]]$data,
                                           nbins, covarlst, nr_of_years,
                                           cslims = cslims, RP = RP,
                                           grid_interval = grid_interval)
    dimlst_rej[[b]] <- tmpres
  }

  occ_from_data <- NULL
  for(b in 1:nbstrp){occ_from_data[[b]] <- dim(margs[[b]]$gpd[[varstr]]$data)[1]}

  par(mfrow = c(1, 1))
  ulim_y <- max(density(unlist(occ_from_data))$y)
  plot(density(unlist(dimlst_rej)),
       main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n",
       lty = 1, col = "brown", ylim = c(0, ulim_y))
  par(new = TRUE)
  plot(density(unlist(occ_from_data)),
       main = "", xlab = "", ylab = "",
       col="orange", lty = 1, ylim = c(0, ulim_y))
  abline(v = mean(unlist(occ_from_data)), col = "orange", lwd = 1, lty = 2)
  abline(v = mean(unlist(dimlst_rej)), col = "brown", lwd = 1, lty = 1)
  legend("topright",
         legend = c("occurrences simulated", "occurrences in data"),
         col = c("brown", "orange"),
         lty = c(1, 3), # line type
         lwd = 1,        # line width
         bty = "n")      # box type ("n" = no box)
}

diagnose_margs_occ_rejection_bstrp_2D <- function(margs, varstr, nbins, covarlst,
                                                  nr_of_years,
                                                  RP = NULL, cslims = NULL,
                                                  condition = NULL,
                                                  grid_interval = NULL,
                                                  plot_ls = TRUE) {
  #' @export

  if (is.null(RP)) {
    RP <- nr_of_years
  }
  dfcountslst <- NULL
  dfobslst <- NULL
  nbstrp <- length(margs)
  for (b in 1:nbstrp) {
    dfin <- margs[[b]]$gpd[[varstr]]$data
    nr_of_events <- dim(dfin)[1]
    dfobslst[[b]] <- dfin[, covarlst]
    dfcounts_tmp <- produce_storm_occurrences_rejection(nr_of_events, RP,
                                                    nr_of_years,
                                                    margs[[b]]$occ[[varstr]],
                                                    covarlst,
                                                    dfin = dfin,
                                                    condition = condition,
                                                    grid_interval = grid_interval)
    dfcountslst[[b]] <- dfcounts_tmp
  }
  dfobs <- bind_rows(dfobslst)
  dfcounts <- bind_rows(dfcountslst)

  if (plot_ls == TRUE) {
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
  } else {
    return(list(dim(dfcounts)[1]))
    }
}

diagnose_margs_occ_pois <- function(mpois, xstr='doy', ystr='Pdir', zstr='counts',
                                    colorbar_intervals=10){
  occ_df <- as.data.frame(mpois$data)

  # Compute fitted and observed values
  fitted_vals   <- exp(unlist(fitted(mpois)))
  observed_vals <- mpois$data[[zstr]]

  # Global range across both plots
  global_min <- min(c(fitted_vals, observed_vals))
  global_max <- max(c(fitted_vals, observed_vals))
  zlim       <- c(global_min, global_max)
  legend_labs <- round(seq(global_min, global_max, length.out = colorbar_intervals))

  par(mfrow = c(1, 2))

  # Simulated/fitted plot
  occ_df[[zstr]] <- fitted_vals
  image(
    x    = sort(unique(occ_df[[xstr]])),
    y    = sort(unique(occ_df[[ystr]])),
    z    = matrix(as.numeric(occ_df[[zstr]]), nrow = length(unique(occ_df[[xstr]]))),
    col  = heat.colors(colorbar_intervals),
    zlim = zlim,
    xlab = xstr,
    ylab = ystr,
    main = "Simulated 2D Grid Counts"
  )
  legend("topright",
         legend = legend_labs,
         fill   = heat.colors(colorbar_intervals),
         title  = "Counts",
         cex = .6)

  # Observed plot
  occ_df[[zstr]] <- observed_vals
  image(
    x    = sort(unique(occ_df[[xstr]])),
    y    = sort(unique(occ_df[[ystr]])),
    z    = matrix(as.numeric(occ_df[[zstr]]), nrow = length(unique(occ_df[[xstr]]))),
    col  = heat.colors(colorbar_intervals),
    zlim = zlim,
    xlab = xstr,
    ylab = ystr,
    main = "Observed 2D Grid Counts"
  )
  legend("topright",
         legend = legend_labs,
         fill   = heat.colors(colorbar_intervals),
         title  = "Counts",
         cex = .6)

  # Sanity check
  print('< Sanity check >')
  print('----------------')
  print(c('Nr of observed counts:', sum(observed_vals)))
  print(c('Nr of simulated counts:', sum(fitted_vals)))
}

#' diagnose_maxds <- function(models_maxds_lst, maxd_thr, X_var_str, Y_var_str, ylim=NULL) {
#'   #' @export
#'
#'   maxd_params_bstrp <- unfold_maxd_params_bstrp_v2(models_maxds_lst, X_var_str, Y_var_str)
#'
#'   # Define the parameters, labels, and data columns
#'   parameters <- c("alpha", "beta", "mu", "sigma")
#'   ylabels <- c("alpha", "beta", "mu", "sigma")
#'
#'   # Set up a 2x2 plotting layout
#'   par(mfrow = c(2, 2))
#'
#'   # Loop through each parameter and generate the plots
#'   for (i in seq_along(parameters)) {
#'     param <- parameters[i]
#'     ylabel <- ylabels[i]
#'
#'     # Extract the corresponding columns for center, lower, and upper limits
#'     p_cntr <- maxd_params_bstrp[[paste0(param, "_cntr")]]
#'     p_llim <- maxd_params_bstrp[[paste0(param, "_llim")]]
#'     p_ulim <- maxd_params_bstrp[[paste0(param, "_ulim")]]
#'
#'     # Remove NA values (if any)
#'     valid_indices <- complete.cases(p_cntr, p_llim, p_ulim)
#'     p_cntr <- p_cntr[valid_indices]
#'     p_llim <- p_llim[valid_indices]
#'     p_ulim <- p_ulim[valid_indices]
#'     maxd_thr_valid <- maxd_thr[valid_indices]
#'
#'     # Set y-axis range
#'     if (is.null(ylim)) {
#'       ylim <- c(-1.5, 1.5)
#'     }
#'
#'
#'     # Plot the center line
#'     plot(maxd_thr_valid, p_cntr, type = "l", ylim = ylim, main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
#'     par(new = TRUE)
#'
#'     # Plot the lower limit
#'     plot(maxd_thr_valid, p_llim, type = "l", ylim = ylim, col = "red", main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
#'     par(new = TRUE)
#'
#'     # Plot the upper limit
#'     plot(maxd_thr_valid, p_ulim, type = "l", ylim = ylim, col = "red", main = "", xlab = "maxd threshold quantile", ylab = ylabel)
#'
#'     # Add horizontal reference lines
#'     abline(h = 0, col = "gray", lwd = .8, lty = 1)
#'     abline(h = 1, col = "gray", lwd = .8, lty = 1)
#'     abline(h = -1, col = "gray", lwd = .8, lty = 1)
#'     abline(h = .5, col = "gray", lwd = .8, lty = 3)
#'     abline(h = -.5, col = "gray", lwd = .8, lty = 3)
#'   }
#' }

diagnose_maxds <- function(models_maxds_lst, maxd_thr, X_var_str, Y_var_str, ylim = NULL) {
  #' @export

  maxd_params_bstrp <- unfold_maxd_params_bstrp_v2(
    models_maxds_lst,
    X_var_str,
    Y_var_str
  )

  # Define the parameters, labels, and data columns
  parameters <- c("alpha", "beta", "mu", "sigma")
  ylabels <- c("alpha", "beta", "mu", "sigma")

  # Default y-axis limits
  if (is.null(ylim)) {
    ylim <- c(-1.5, 1.5)
  }

  # Set up a 2x2 plotting layout
  par(mfrow = c(2, 2))

  # Loop through each parameter and generate the plots
  for (i in seq_along(parameters)) {

    param <- parameters[i]
    ylabel <- ylabels[i]

    # Extract the corresponding columns
    p_cntr <- maxd_params_bstrp[[paste0(param, "_cntr")]]
    p_llim <- maxd_params_bstrp[[paste0(param, "_llim")]]
    p_ulim <- maxd_params_bstrp[[paste0(param, "_ulim")]]

    # Remove NA values
    valid_indices <- complete.cases(p_cntr, p_llim, p_ulim)
    p_cntr <- p_cntr[valid_indices]
    p_llim <- p_llim[valid_indices]
    p_ulim <- p_ulim[valid_indices]
    maxd_thr_valid <- maxd_thr[valid_indices]

    # Use supplied limits, but force alpha and sigma to start at 0
    ylim_use <- ylim
    if (param %in% c("alpha", "sigma")) {
      ylim_use[1] <- 0
    }

    # Plot center estimate
    plot(
      maxd_thr_valid, p_cntr,
      type = "l",
      ylim = ylim_use,
      main = "",
      xlab = "",
      ylab = "",
      xaxt = "n",
      yaxt = "n"
    )
    par(new = TRUE)

    # Plot lower confidence limit
    plot(
      maxd_thr_valid, p_llim,
      type = "l",
      ylim = ylim_use,
      col = "red",
      main = "",
      xlab = "",
      ylab = "",
      xaxt = "n",
      yaxt = "n"
    )
    par(new = TRUE)

    # Plot upper confidence limit
    plot(
      maxd_thr_valid, p_ulim,
      type = "l",
      ylim = ylim_use,
      col = "red",
      main = "",
      xlab = "maxd threshold quantile",
      ylab = ylabel
    )

    # Add horizontal reference lines
    abline(h = 0,   col = "gray", lwd = 0.8, lty = 1)
    abline(h = 1,   col = "gray", lwd = 0.8, lty = 1)
    abline(h = -1,  col = "gray", lwd = 0.8, lty = 1)
    abline(h = 0.5, col = "gray", lwd = 0.8, lty = 3)
    abline(h = -0.5, col = "gray", lwd = 0.8, lty = 3)
  }
}

diagnose_maxds_fitted <- function(maxd, X_str, Y_str, nsim) {
  #' @export

  par(mfrow = c(2, 2))
  # plot Laplace margins
  X_all <- maxd[[X_str]][[Y_str]]$X_all
  X_fit <- maxd[[X_str]][[Y_str]]$X_fit
  Y_all <- maxd[[X_str]][[Y_str]]$Y_all
  Y_fit <- maxd[[X_str]][[Y_str]]$Y_fit
  xlim <- c(min(c(X_all, Y_all)),
            max(c(X_all, Y_all)))
  ylim <- xlim
  plot(X_all, Y_all,
       pch = 20, col = "black", xlim = xlim, ylim = ylim, cex = .5,
       main = "", xlab = "X on Laplace", ylab = "Y on Laplace")
  par(new = TRUE)
  plot(X_fit, Y_fit,
       pch = 20, col = "orange", xlim = xlim, ylim = ylim, cex = .5,
       main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  abline(h = 0, col = "gray", lwd = .5, lty = 1)
  abline(v = 0, col = "gray", lwd = .5, lty = 1)
  abline(v = quantile(X_all, maxd[[X_str]][[Y_str]]$thr),
         col = "orange", lwd = .5, lty = 1)

  # plot simulations from HT2004
  sims <- predict_maxd(maxd, X_str, Y_str, nsim)
  lp_X <- sims$Xsim
  lp_Y <- sims$Ysim
  xlim <- c(min(c(X_all, unlist(sims))), max(c(X_all, unlist(sims))))
  ylim <- xlim
  plot(X_all, Y_all,
       pch = 20, col = "black", xlim = xlim, ylim = ylim, cex = .5,
       main = "", xlab = "X on Laplace", ylab = "Y on Laplace")
  par(new = TRUE)
  plot(lp_X, lp_Y,
       pch = 20, col = "orange", xlim = xlim, ylim = ylim, cex = .5,
       main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  abline(h = 0, col = "gray", lwd = .5, lty = 1)
  abline(v = 0, col = "gray", lwd = .5, lty = 1)
  abline(v = quantile(X_all, maxd[[X_str]][[Y_str]]$thr),
         col = "orange", lwd = .5, lty = 1)

  # plot division when sim X or sim Y is largest (diagonal split)
  plot(X_all, Y_all,
       pch = 20, col = "black", xlim = xlim, ylim = ylim, cex = .5,
       main = "", xlab = "X on Laplace", ylab = "Y on Laplace")
  for (i in 1:length(lp_X)) {
    par(new = TRUE)
    if (lp_X[i] > lp_Y[i]) {
      plot(lp_X[i], lp_Y[i],
           pch = 20, col = "red", xlim = xlim, ylim = ylim, cex = .5,
           main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    } else {
      plot(lp_X[i], lp_Y[i],
           pch = 20, col = "cyan", xlim = xlim, ylim = ylim, cex = .5,
           main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    }
  }
  abline(h = 0, col = "gray", lwd = .5, lty = 1)
  abline(v = 0, col = "gray", lwd = .5, lty = 1)
  abline(v = quantile(X_all, maxd[[X_str]][[Y_str]]$thr),
         col = "orange", lwd = .5, lty = 1)
  abline(coef = c(0, 1), col = "orange")

  # plot division when sim X or sim Y is largest (diagonal split)
  # and simulated values from respective model
  sims <- predict_maxd(maxd, X_str, Y_str, nsim)
  lp_X <- sims$Xsim
  lp_Y <- sims$Ysim

  plot(X_all, Y_all,
       pch = 20, col = "black", xlim = xlim, ylim = ylim, cex = .5,
       main = "", xlab = "X on Laplace", ylab = "Y on Laplace")
  for (i in 1:length(lp_X)) {
    par(new=TRUE)
    if (lp_X[i] > lp_Y[i]) {
      plot(lp_X[i], lp_Y[i],
           pch = 20, col = "red", xlim = xlim, ylim = ylim, cex = .5,
           main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    }
  }
  sims <- predict_maxd(maxd, Y_str, X_str, nsim)
  lp_X <- sims$Xsim
  lp_Y <- sims$Ysim
  for (i in 1:length(lp_X)) {
    par(new = TRUE)
    if (lp_X[i] > lp_Y[i]) {
      plot(lp_Y[i], lp_X[i],
           pch = 20, col = "cyan", xlim = xlim, ylim = ylim, cex = .5,
           main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    }
  }

  abline(h = 0, col = "gray", lwd = .5, lty = 1)
  abline(v = 0, col = "gray", lwd = .5, lty = 1)

  abline(v = quantile(X_all, maxd[[X_str]][[Y_str]]$thr),
         col = "orange", lwd = .5, lty = 1)
  abline(h = quantile(X_all, maxd[[Y_str]][[X_str]]$thr),
         col = "orange", lwd = .5, lty = 1)
  abline(coef = c(0, 1), col = "orange")
}

diagnose_margs_preds_density <- function(margs_preds, varstr,
                                         n.bstr = NULL, bw = NULL,
                                         xlim = NULL, ylim = NULL,
                                         qlow = .05, qhigh = .95, qmed = .5,
                                         xlab = NULL) {
  #' @export
  #'
  if (is.null(n.bstr)) {
    n.bstr <- seq(length(margs_preds))
  }

  maxvals <- NULL
  for (b in n.bstr) {
    tmp <- margs_preds[[b]][[varstr]]$maxval
    maxvals[[b]] <- tmp[is.finite(tmp)]
  }

  if (is.null(xlim)) {
    xlim <- c(min(unlist(maxvals), na.rm = TRUE),
              max(unlist(maxvals), na.rm = TRUE))
  }

  y_max <- array(dim = length(maxvals))
  for (b in 1:length(y_max)) {
    if (!is.na(maxvals[[b]][1])){
      y_max[b] <- max(density(maxvals[[b]], bw = bw, n = 512)$y)
    }
  }

  if (is.null(ylim)) {
    ylim <- c(0, max(y_max, na.rm=TRUE))
  }

  if (is.null(bw)) {
    bw <- .01 * abs(xlim[1] - xlim[2])
  }

  # plot combined density with individual densities
  for (b in n.bstr) {
    if (!is.na(maxvals[[b]][1])){
      plot(density(maxvals[[b]], bw = bw, n = 512,
                   from = xlim[1], to = xlim[2]),
           col = "grey", lwd = .5, xlim = xlim, ylim = ylim,
           xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = "")
      par(new = TRUE)
    }
  }
  if (is.null(xlab)) {
    xlab <- varstr
  }

  plot(density(na.omit(unlist(maxvals)), bw = bw, n = 512,
               from = xlim[1], to = xlim[2]),
       xlim = xlim, ylim = ylim, main = "", xlab = xlab)

  # plot density with uncertainty quantiles
  densities <- NULL
  for (b in n.bstr){
    if (!is.na(maxvals[[b]][1])){
      densities[["d"]][[b]] <- density(maxvals[[b]], bw = bw, n = 512,
                                       from = xlim[1], to = xlim[2])
      densities[["x"]][[b]] <- densities[["d"]][[b]]$x
      densities[["y"]][[b]] <- densities[["d"]][[b]]$y
    }
  }

  # filter out empty entries
  filtered_y <- Filter(function(x) !is.null(x) && length(x) > 0, densities$y)
  df_dens <- data.frame(filtered_y)

  qlow_ts <- apply(df_dens, 1, quantile, qlow)
  qhigh_ts <- apply(df_dens, 1, quantile, qhigh)
  qmed_ts <- apply(df_dens, 1, quantile, qmed)

  plot(densities$x[[1]], qlow_ts, col = "blue", xlim = xlim, ylim = ylim, main = "",
       lty = 3, lwd = .5, type = "l", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  par(new = TRUE)
  plot(densities$x[[1]], qhigh_ts, col = "blue", xlim = xlim, ylim = ylim, main = "",
       lty = 3, lwd = .5, type = "l", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  par(new = TRUE)
  plot(densities$x[[1]], qmed_ts, xlim = xlim, ylim = ylim, main = "",
       xlab = xlab, ylab = "Density", lty = 3, lwd = .5, type = "l", col = "blue")
  par(new = TRUE)
  plot(density(na.omit(unlist(maxvals)), bw = bw, n = 512,
               from = xlim[1], to = xlim[2]),
       xlim = xlim, ylim = ylim, main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n")

  print("summary maxes:")
  print(summary(na.omit(unlist(maxvals))))
  return(list("maxvals" = maxvals, "densities" = densities))
}

display_joint_densities <- function(preds_maxds,
                                    preds_maxds_LP,
                                    preds_margs_lst,
                                    X_str, Y_str,
                                    Xmax = NULL, Ymax = NULL) {
  # use function to retrieve valids
  valids <- retrieve_valid_HT_samples(preds_maxds,
                                      preds_maxds_LP,
                                      preds_margs_lst,
                                      X_str, c(Y_str))

  library(hexbin)
  library(RColorBrewer)
  library(lattice)
  rf <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  # Create hexbin object and plot
  dfin <- data.frame(unlist(valids$X_valid[[Y_str]]), unlist(valids$Y_valid[[Y_str]]))
  colnames(dfin) <- c(X_str, Y_str)
  h <- hexbin(dfin)
  counts <- h@count
  formula <- as.formula(paste(Y_str, "~", X_str))
  # Define a custom panel function to add a vertical line
  custom_panel <- function(x, y, Xmax, Ymax, ...) {
    panel.hexbinplot(x, y, ...)
    panel.abline(v = Xmax, col = "gray", lwd = 1, lty = 1)
    panel.abline(h = Ymax, col = "gray", lwd = 1, lty = 1)
  }
  hexbinplot(formula, data = dfin,
             xbins = 50, colramp = rf,
             mincnt = (min(counts) + 1),
             maxcnt = max(counts),
             trans = log, inv = exp,
             panel = custom_panel,
             Xmax = Xmax, Ymax = Ymax)
}

plot_cvres <- function(cvres, limits = NULL, show_errors = NULL) {
  #' Plots the results from the cross-validation procedure
  #' @export
  #'
  # print values of minimum
  print(c(cvres$thr[cvres$cost$costfct1 == min(cvres$cost$costfct1)],
          #cvres$thr[cvres$cost$costfct2 == min(cvres$cost$costfct2)],
          cvres$thr[cvres$cost$costfct3 == min(cvres$cost$costfct3)]))

  if (is.null(limits)) {
    # plot mean cost-function
    par(mfrow = c(1, 2))
    plot(cvres$thr, cvres$cost$costfct1, xlab = "Threshold", main = "MSE based", ylab = "Costfunction 1")
    abline(v = cvres$thr[cvres$cost$costfct1 == min(cvres$cost$costfct1)])
    splc1 <- smooth.spline(cvres$thr, cvres$cost$costfct1, df=5)
    lines(splc1, col = "red", lwd = 2)
    #plot(cvres$thr, cvres$cost$costfct2)
    #abline(v = cvres$thr[cvres$cost$costfct2 == min(cvres$cost$costfct2)])
    plot(cvres$thr, cvres$cost$costfct3, xlab = "Threshold", main = "MAE based", ylab = "Costfunction 2")
    abline(v = cvres$thr[cvres$cost$costfct3 == min(cvres$cost$costfct3)])
    splc3 <- smooth.spline(cvres$thr, cvres$cost$costfct3, df=5)
    lines(splc3, col = "red", lwd = 2)

  } else {
    # Add the shaded region
    # Define the polygon coordinates
    x1 = limits[1]
    x2 = limits[2]
    polygon_x <- c(x1, x2, x2, x1)
    polygon_y <- c(-1000, -1000, 1000, 1000) # Extend to the plot limits vertically
    # plot mean cost-function
    par(mfrow = c(1, 2))
    plot(cvres$thr, cvres$cost$costfct1, xlab = "Threshold", main = "MSE based", ylab = "Costfunction 1")
    polygon(polygon_x, polygon_y, col = rgb(0.6, 0.6, 0.6, 0.3), border = NA)
    abline(v = cvres$thr[cvres$cost$costfct1 == -min(cvres$cost$costfct1)])
    splc1 <- smooth.spline(cvres$thr, cvres$cost$costfct1, df=5)
    lines(splc1, col = "red", lwd = 2)
    #plot(cvres$thr, cvres$cost$costfct2)
    #polygon(polygon_x, polygon_y, col = rgb(0.6, 0.6, 0.6, 0.3), border = NA)
    #abline(v = cvres$thr[cvres$cost$costfct2 == min(cvres$cost$costfct2)])
    plot(cvres$thr, cvres$cost$costfct3, xlab = "Threshold", main = "MAE based", ylab = "Costfunction 2")
    polygon(polygon_x, polygon_y, col = rgb(0.6, 0.6, 0.6, 0.3), border = NA)
    abline(v = cvres$thr[cvres$cost$costfct3 == min(cvres$cost$costfct3)])
    splc3 <- smooth.spline(cvres$thr, cvres$cost$costfct3, df=5)
    lines(splc3, col = "red", lwd = 2)
  }

  if (!is.null(show_errors)) {
    # plot errors
    par(mfrow = c(2, 2))
    boxplot(t(cvres$errors$bias), names = cvres$thr, main = "BIAS")
    abline(h = 0)
    boxplot(t(cvres$errors$mae), names = cvres$thr, main = "MAE")
    abline(h = 0)
    boxplot(t(cvres$errors$rse), names = cvres$thr, main = "RSE")
    abline(h = 0)
    boxplot(t(cvres$errors$drse), names = cvres$thr, main = "dRSE")
    abline(h = 0)
  }
}

visualize_storm_picking <- function(dfin, dfinall = NULL, xstr = "dt",
                                    ystr = "hs", sidx = NULL, eidx = NULL,
                                    ylim = c(0, 15), xlab = "", ylab = "Hs [m]") {
  #' function to visualize the outcome of the storm picking procedure
  #'
  #' @export
  #'

  if (is.null(sidx)) {
    sidx <- 1
    }
  if (is.null(eidx)) {
    eidx <- length(dfin[[ystr]])
    }

  if (!is.null(dfinall)) {
    plot(dfinall[[xstr]][sidx:(eidx * 3)], dfinall[[ystr]][sidx:(eidx * 3)],
         pch = 20, cex = .5,
         col = adjustcolor("gray", alpha.f = 0.5),
         xlim = c(dfin$storms[[xstr]][sidx], dfin$storms[[xstr]][eidx]), ylim = ylim,
         xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = "")  # e.g. all hs from data
    par(new = TRUE)
  }
  plot(dfin$storms[[xstr]][sidx:(eidx * 2)], dfin$storms$thr[sidx:(eidx * 2)],
       pch = 20, cex = .5,
       col = adjustcolor("red", alpha.f = 0.5),
       xlim = c(dfin$storms[[xstr]][sidx], dfin$storms[[xstr]][eidx]), ylim = ylim,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = "")  # non-stationary local threshold
  par(new = TRUE)
  plot(dfin$storms[[xstr]][sidx:(eidx * 2)], dfin$storms$exc[sidx:(eidx * 2)] + dfin$storms$thr[sidx:(eidx * 2)],
       pch = 20, cex = .5,
       col = adjustcolor("black", alpha.f = 0.5),
       xlim = c(dfin$storms[[xstr]][sidx], dfin$storms[[xstr]][eidx]), ylim = ylim,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = "")  # exceedences
  par(new = TRUE)
  plot(dfin$pots[[xstr]][sidx:eidx], dfin$pots[[ystr]][sidx:eidx],
       pch = 1, cex = 1.,
       col = adjustcolor("orange", alpha.f = 0.6),
       xlim = c(dfin$storms[[xstr]][sidx], dfin$storms[[xstr]][eidx]), ylim = ylim,
       xlab = xlab, ylab = ylab, main = "")  # peaks
  legend("topright",
         legend = c("All data", "Threshold", "Exceedances", "Peaks"),
         col = c("gray", "red", "black", "orange"),
         pch = c(20, 20, 20, 1),
         pt.cex = c(0.5, 0.5, 0.5, 1),
         bty = "n")
}

vis_sim_storms <- function(res, xlim=c(0,15), ylim=c(0,25), storm_idx=1){

  df_storm_hist <- res$hist_storms[(res$hist_storms$pseudo_storm_idx==storm_idx), ]
  df_storm_sim <- res$sim_storms[(res$sim_storms$pseudo_storm_idx==storm_idx), ]

  par(mfrow = c(2, 2))
  # joint Hs/Tm02
  plot(res$hist_storms$tm2, res$hist_storms$hs, pch=20, xlim=xlim, ylim=ylim, cex=.3, xlab='', ylab='', main='', xaxt = 'n', yaxt = 'n', col='gray')
  par(new=TRUE)
  plot(df_storm_hist$tm2, df_storm_hist$hs, pch=20, xlim=xlim, ylim=ylim, type='o', cex=1, lty=1, lwd=1, xlab='Tm02 [s]', ylab='Hs [s]', main='', col='blue')
  par(new=TRUE)
  plot(df_storm_sim$tm2, df_storm_sim$hs, pch=20, xlim=xlim, ylim=ylim, type='o', cex=1, lty=1, lwd=1, xlab='Tm02 [s]', ylab='Hs [s]', main='', col='red')
  legend("topleft", legend = c("All hist", "Traj Hist", "Traj Sim"),
         col = c("grey", "blue", "red"), pch = c(20, 20, 20), lty = c(NA, 1, 1),
         bg = "white", bty = "n")
  mtext("a)", side = 3, line = -1.5, adj = 0.95, cex = 1)

  # Ts Hs
  plot(df_storm_hist$dt, df_storm_hist$hs, xlab='', ylab='Hs [m]', ylim=ylim, main='', type='o', lty=1, pch=20, col="blue")
  par(new=TRUE)
  plot(df_storm_hist$dt, df_storm_sim$hs, xlab='', ylab='Hs [m]', ylim=ylim, main='', type='o', lty=1, pch=20, col="red")
  mtext("b)", side = 3, line = -1.5, adj = 0.95, cex = 1)

  # Polar plot historic
  max_extent <- max(ylim)
  magnitudes <- df_storm_hist$hs  # Magnitudes remain unchanged
  # Rotate angles to align with the oceanographic convention (add 90°) and transofrm to radians
  angles_rotated <- (90 - df_storm_hist$Pdir) %% 360 * (pi/180)
  # Convert polar coordinates to Cartesian coordinates
  x <- magnitudes * cos(angles_rotated)
  y <- magnitudes * sin(angles_rotated)
  # Define the plotting limits
  lims <- c(-max_extent, max_extent)
  # Plot a reference circle
  plot(0, 0, xlim = lims, ylim = lims, type = "n", asp = 1, xlab = "Hs [X component]", ylab = "Hs [Y component]", main = "")
  symbols(0, 0, circles = .1, add = TRUE, inches = FALSE, lty = 2)  # Add a reference circle
  # Add radial gridlines
  rad_labels <- seq(0, 2 * pi, length.out = 13)  # 12 evenly spaced angles (30° apart)
  for (angle in rad_labels) {
    angle_rotated <- angle + pi / 2  # Rotate the gridlines
    segments(0, 0, cos(angle_rotated) * max_extent, sin(angle_rotated) * max_extent,
             lty = 3, col = "gray")
  }
  # Add concentric circles
  for (r in seq(2, max_extent, by = 2)) {
    symbols(0, 0, circles = r, add = TRUE, inches = FALSE, lty = 1, col = "gray", lwd=.5)
  }
  # Add the line plot
  lines(x, y, col = "blue", lty = 1)
  # Highlight the first and last points
  points(x[1], y[1], pch = 19, col = "blue", cex = 1)  # First point
  points(x[length(x)], y[length(y)], pch = 2, col = "blue", cex = 1)  # Last point
  # Add degree labels in oceanographic convention
  label_radius <- max_extent  # Place labels slightly outside the largest magnitude
  deg_labels <- c(0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330)  # Oceanographic convention
  rad_labels <- (90 - deg_labels) %% 360 * (pi/180)
  # Place degree labels
  text(label_radius * cos(rad_labels), label_radius * sin(rad_labels),
       labels = deg_labels, col = "red", cex = 0.8)

  # Polar plot simulated
  magnitudes <- df_storm_sim$hs  # Magnitudes remain unchanged
  # Rotate angles to align with the oceanographic convention (add 90°) and transofrm to radians
  angles_rotated <- (90 - df_storm_sim$Pdir) %% 360 * (pi/180)
  # Convert polar coordinates to Cartesian coordinates
  x <- magnitudes * cos(angles_rotated)
  y <- magnitudes * sin(angles_rotated)

  # Add the line plot
  lines(x, y, col = "red", lty = 1)
  # Highlight the first and last points
  points(x[1], y[1], pch = 19, col = "red", cex = 1)  # First point
  points(x[length(x)], y[length(y)], pch = 2, col = "red", cex = 1)  # Last point

  # Add a legend
  legend("topleft", legend = c("Hist", "Sim", "Start", "End"),
         col = c("blue", "red", "grey", "grey"), pch = c(NA, NA, 19, 2), lty = c(1, 1, NA, NA),
         bg = "white", bty = "n")
  mtext("c)", side = 3, line = -1.5, adj = 0.95, cex = 1)

  # Ts Tm2
  plot(df_storm_hist$dt, df_storm_hist$tm2, xlab='', ylab='Tm02 [s]', ylim=xlim, main='', type='o', lty=1, pch=20, col="blue")
  par(new=TRUE)
  plot(df_storm_hist$dt, df_storm_sim$tm2, xlab='', ylab='Tm02 [s]', ylim=xlim, main='', type='o', lty=1, pch=20, col="red")
  mtext("d)", side = 3, line = -1.5, adj = 0.95, cex = 1)
}

show_sim_pop <- function(res, xlim=c(0,15), ylim=c(0,25), xstr="tm2", ystr="hs", steepness=FALSE){
  # show population of all simulated RP storms
  par(mfrow = c(1, 1))

  # joint xstr/ystr (e.g. tm2/hs)
  for (i in 1:length(unique(res$sim_storms$pseudo_storm_idx))){
    indiv_storm <- subset(res$sim_storms, pseudo_storm_idx==i)
    plot(indiv_storm[[xstr]], indiv_storm[[ystr]], pch=20, xlim=xlim, ylim=ylim, type='o', cex=.3, lty=1, lwd=.2, xlab='', ylab='', main='', xaxt = 'n', yaxt = 'n',  col = adjustcolor("red", alpha = 0.1))
    par(new=TRUE)
  }
  plot(res$hist_storms[[xstr]], res$hist_storms[[ystr]], pch=20, xlim=xlim, ylim=ylim, cex=.3, xlab='Tm02 [s]', ylab='Hs [m]', main='', col = adjustcolor("black", alpha = 0.4) )

  for (i in 1:length(unique(res$sim_storms$pseudo_storm_idx))){
    indiv_storm <- subset(res$sim_storms, pseudo_storm_idx == i)
    par(new=TRUE)
    indiv_peak <- subset(indiv_storm, hs == max(indiv_storm[[ystr]]))
    #plot(indiv_peak[[xstr]], indiv_peak[[ystr]], pch=1, xlim=xlim, ylim=ylim, cex=.5, xlab='', ylab='', main='', xaxt = 'n', yaxt = 'n',  col = adjustcolor("orange", alpha.f = .5))
    plot(indiv_peak[[xstr]], indiv_peak[[ystr]], pch=21, xlim=xlim, ylim=ylim, cex=.5,
         xlab='', ylab='', main='', xaxt = 'n', yaxt = 'n',
         col = "grey30", bg = adjustcolor("darkorange", alpha.f = .8), lwd = 0.3)
  }
  for (i in 1:length(unique(res$hist_storms$pseudo_storm_idx))){
    indiv_storm <- subset(res$hist_storms, pseudo_storm_idx == i)
    par(new=TRUE)
    indiv_peak <- subset(indiv_storm, hs == max(indiv_storm[[ystr]]))
    plot(indiv_peak[[xstr]], indiv_peak[[ystr]], pch=21, xlim=xlim, ylim=ylim, cex=.5,
         xlab='', ylab='', main='', xaxt = 'n', yaxt = 'n',
         col = "grey30", bg = adjustcolor("lightblue", alpha.f = .8), lwd = 0.3)
  }

  if (steepness == TRUE){
    # Define the valid domains for each steepness limit function
    Tz1 <- seq(0.1, 6, 0.1)  # Domain for the first function
    HS1 <- 1/10 * 9.81 * Tz1^2 / (2 * pi)  # Values for the first function

    Tz2 <- seq(12.1, 20, 0.1)  # Domain for the second function
    HS2 <- 1/15 * 9.81 * Tz2^2 / (2 * pi)  # Values for the second function

    # Define the range for interpolation between the two domains
    Tz_interp <- seq(6.1, 12, 0.1)  # Intermediate range for interpolation

    # Interpolate linearly between the two functions
    HS_interp <- approx(
      x = c(Tz1[length(Tz1)], Tz2[1]),   # Boundary points for interpolation
      y = c(HS1[length(HS1)], HS2[1]),  # Corresponding values at the boundaries
      xout = Tz_interp                  # Points to interpolate
    )$y

    # Combine all the values
    Tz_combined <- c(Tz1, Tz_interp, Tz2)  # Combined domain
    HS_combined <- c(HS1, HS_interp, HS2) # Combined values

    lines(Tz_combined, HS_combined, lty=2)
  }
  legend(
    "topleft",
    legend = c("all historic storms","sample historic storm peaks","sample simulated storms","sample simulated peaks"),
    col = c("black", "grey30", "red", "grey30"),
    pt.bg = c(NA, adjustcolor("lightblue", alpha.f = .8), NA, adjustcolor("darkorange", alpha.f = .8)),
    pch = c(20, 21, 20, 21),
    lty = c(NA, NA, 1, NA),
    lwd = c(NA, NA, 1, NA),
    pt.cex = c(1, 1, 1, 1),
    bty = "n"
  )
}

plot_cvres_distr <- function(cvres, ylim=NULL) {
  #' Plots the results from the cross-validation procedure
  #' @export
  #'
  par(mfrow = c(1, 3))
  boxplot(cvres$all_score, names = cvres$thr, main = "all score", ylim=ylim)
  boxplot(cvres$rt_score, names = cvres$thr, main = "rt score", ylim=ylim)
  boxplot(cvres$crps_score, names = cvres$thr, main = "crps score", ylim=ylim)
}

plot_preds_2d_gg <- function(df_preds, xvar, yvar,
                             df_obs = NULL,
                             bins = 100,
                             density_levels = c(0.01, 0.05, 0.1),
                             xlim = NULL, ylim = NULL) {
  #' 2d visualisation of df_preds point cloud
  #' @export

  # Set axis limits from data if not provided
  if (is.null(xlim)) xlim <- range(df_preds[[xvar]])
  if (is.null(ylim)) ylim <- range(df_preds[[yvar]])

  # Expand limits to include observations if provided
  if (!is.null(df_obs)) {
    xlim <- range(c(xlim, df_obs[[xvar]]))
    ylim <- range(c(ylim, df_obs[[yvar]]))
  }

  # Compute 2D kernel density manually with MASS::kde2d
  dens <- MASS::kde2d(df_preds[[xvar]], df_preds[[yvar]], n = 200,
                      lims = c(xlim, ylim))

  # Convert to long dataframe for ggplot
  df_dens <- expand.grid(x = dens$x, y = dens$y)
  df_dens$z <- as.vector(dens$z)
  names(df_dens)[1:2] <- c(xvar, yvar)

  p <- ggplot() +
    geom_bin2d(data = df_preds,
               aes(x = .data[[xvar]], y = .data[[yvar]]),
               bins = bins) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    geom_contour(data = df_dens,
                 aes(x = .data[[xvar]], y = .data[[yvar]], z = z),
                 color = "red",
                 linewidth = 0.3,
                 breaks = density_levels) +
    geom_text_contour(data = df_dens,
                      aes(x = .data[[xvar]], y = .data[[yvar]], z = z),
                      color = "red",
                      size = 3,
                      breaks = density_levels,
                      skip = 0,
                      stroke = 0.2) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(x = xvar, y = yvar, fill = "count") +
    theme_bw()

  if (!is.null(df_obs)) {
    p <- p + geom_point(data = df_obs,
                        aes(x = .data[[xvar]], y = .data[[yvar]]),
                        color = "black", size = 0.3, alpha = 0.6)
  }

  p
}

# Update plot_preds_2d_gg to accept fixed color limits and hide colorbar
plot_preds_2d_gg_shared <- function(df_preds, xvar, yvar,
                                    df_obs = NULL,
                                    bins = 100,
                                    density_levels = c(0.01, 0.05, 0.1),
                                    xlim = NULL, ylim = NULL,
                                    clim = NULL,
                                    show_legend = TRUE) {
  #' @export

  if (is.null(xlim)) xlim <- range(df_preds[[xvar]])
  if (is.null(ylim)) ylim <- range(df_preds[[yvar]])

  if (!is.null(df_obs)) {
    xlim <- range(c(xlim, df_obs[[xvar]]))
    ylim <- range(c(ylim, df_obs[[yvar]]))
  }

  dens <- MASS::kde2d(df_preds[[xvar]], df_preds[[yvar]], n = 200,
                      lims = c(xlim, ylim))
  df_dens <- expand.grid(x = dens$x, y = dens$y)
  df_dens$z <- as.vector(dens$z)
  names(df_dens)[1:2] <- c(xvar, yvar)

  p <- ggplot() +
    geom_bin2d(data = df_preds,
               aes(x = .data[[xvar]], y = .data[[yvar]]),
               bins = bins) +
    scale_fill_gradient(low = "white", high = "steelblue",
                        limits = clim,
                        guide = if (show_legend) "colorbar" else "none") +
    geom_contour(data = df_dens,
                 aes(x = .data[[xvar]], y = .data[[yvar]], z = z),
                 color = "red", linewidth = 0.3,
                 breaks = density_levels) +
    geom_text_contour(data = df_dens,
                      aes(x = .data[[xvar]], y = .data[[yvar]], z = z),
                      color = "red", size = 3,
                      breaks = density_levels,
                      skip = 0, stroke = 0.2) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(x = xvar, y = yvar, fill = "count") +
    theme_bw()

  if (!is.null(df_obs)) {
    p <- p + geom_point(data = df_obs,
                        aes(x = .data[[xvar]], y = .data[[yvar]]),
                        color = "black", size = 0.3, alpha = 0.6)
  }

  p
}

prepare_shared_plot_params <- function(df_list, xvar, yvar, df_obs = NULL, bins = 100) {
  #' @export
  #'
  # Compute global axis limits across all dataframes
  xlim <- range(c(unlist(lapply(df_list, function(df) df[[xvar]])),
                  if (!is.null(df_obs)) df_obs[[xvar]]))
  ylim <- range(c(unlist(lapply(df_list, function(df) df[[yvar]])),
                  if (!is.null(df_obs)) df_obs[[yvar]]))

  # Compute global bin breaks
  x_breaks <- seq(xlim[1], xlim[2], length.out = bins + 1)
  y_breaks <- seq(ylim[1], ylim[2], length.out = bins + 1)

  # Compute global max count for shared color scale
  max_count <- max(sapply(df_list, function(df) {
    x_idx <- findInterval(df[[xvar]], x_breaks)
    y_idx <- findInterval(df[[yvar]], y_breaks)
    max(table(paste(x_idx, y_idx)))
  }))

  list(xlim  = xlim,
       ylim  = ylim,
       clim  = c(0, max_count),
       bins  = bins)
}

plot_ly_3d <- function(df_preds, xvar, yvar, zvar, df_obs = NULL,
                       xlim = NULL, ylim = NULL, zlim = NULL) {
  #' 3d visualisation of df_preds point cloud
  #' @export
  p <- plot_ly()
  p <- add_trace(p,
                 data = df_preds,
                 x = df_preds[[xvar]], y = df_preds[[yvar]], z = df_preds[[zvar]],
                 type = "scatter3d", mode = "markers",
                 name = "predictions",
                 marker = list(size = 2, color = "steelblue", opacity = 0.4)
  )
  if (!is.null(df_obs)) {
    p <- add_trace(p,
                   data = df_obs,
                   x = df_obs[[xvar]], y = df_obs[[yvar]], z = df_obs[[zvar]],
                   type = "scatter3d", mode = "markers",
                   name = "observed",
                   marker = list(size = 3, color = "black", opacity = 0.8)
    )
  }
  p <- layout(p,
              scene = list(
                xaxis = list(title = xvar, range = xlim),
                yaxis = list(title = yvar, range = ylim),
                zaxis = list(title = zvar, range = zlim)
              )
  )
  p
}
