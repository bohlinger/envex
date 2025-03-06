# Heffernan and Tawn conditional extremes model from 2004 as used in e.g.
# Ewans and Jonathan (2013) and Jonathan, Ewans, Randell (2013).

library(extRemes)
# library(extraDistr)

make_marginal_models_manual <- function(mexdata, mqu=.9, time_units=NULL,
                                        number_of_years=NULL){
  #' @export
  margs <- NULL
  if (is.null(time_units)){
    time_units <- sprintf("%.4f/year", length(mexdata$hs)/number_of_years)
  } else {time_units <- time_units}
  print("time_units are:")
  print(time_units)
  if (length(mqu)==length(colnames(mexdata))){
    print("use all mqus in order")
    print("not yet implemented")
  } else {
    print("only use 1mqu (i.e. mqu[1])")
    for (v in colnames(mexdata)){
      margs[[v]] <-fevd(mexdata[[v]], data = mexdata,
                        type = "GP",
                        threshold = quantile(mexdata[[v]], mqu[1]),
                        time.units = time_units,
                        use.phi = TRUE)
      margs[[v]][["mqu"]] <- mqu[1]
    }
  }
  return(margs)
}

make_marginal_models <- function(mexdata, mqu=.9,
                                 time_units=NULL,
                                 number_of_years=NULL){
  #' @description
  #' Function to create marginal models of considered joint variables.
  #' Input should be i.i.d. peaks.
  #' @param mexdata data of joint variables
  #' @param mqu threshold for marginal models as quantile. Could be an array.
  #'
  #' @return margs object of marginal models
  #'
  #' @examples
  #'
  #'
  #' @export

  margs <- NULL
  margs <- make_marginal_models_manual(mexdata=mexdata, mqu=mqu,
                                       number_of_years=number_of_years,
                                       time_units=time_units)
  return(margs)
}

transform_to_laplace <- function(dfin, margs, names_in=NULL,
                                 maxd=NULL, which.covariate=NULL){
  #' @export

  lp_margins <- NULL
  if (is.null(names_in)){
    names_in <- names(margs)
  }
  for (n in names_in){
    natural_scale <- dfin[[n]]
    Fn_ecdf <- ecdf(natural_scale)
    probs_ecdf <- Fn_ecdf(natural_scale[natural_scale<=margs[[n]]$threshold])
    probs_gpd <- pevd(natural_scale[natural_scale>margs[[n]]$threshold],
                      scale = exp(margs[[n]]$results$par[1]),
                      shape = margs[[n]]$results$par[2],
                      threshold = margs[[n]]$threshold,
                      type = c("GP"),
                      lower.tail = TRUE)
    # combine probs
    probs <- c(probs_ecdf,probs_gpd)

    # transform to Laplace
    if (is.null(maxd)){
      lp_margins[[n]] <- qlaplace(probs)
    } else{
      if (n==which.covariate){
        lp_margins[[n]] <- qlaplace(probs)
      } else{
        lp_margins[[n]] <- qlaplace(probs,
                                    maxd$dependence$coefficients['m'],
                                    maxd$dependence$coefficients['s'])
      }
    }
  }
  return(lp_margins)
}

transform_to_laplace_ndim_preds <- function(dfin, maxd=NULL,
                                            which.covariate=NULL){
  #' @export

  threshold <- maxd[[which.covariate]][[1]]$margs[[which.covariate]]$threshold
  natural_scale <- dfin[which.covariate][[1]]
  tmpdata <- natural_scale[natural_scale>threshold]
  probs_gpd <- pevd(tmpdata,
                    scale = exp(margs[[which.covariate]]$results$par[1]),
                    shape = margs[[which.covariate]]$results$par[2],
                    threshold = margs[[which.covariate]]$threshold,
                    type = c("GP"),
                    lower.tail = TRUE)
  # combine probs
  probs <- c(probs_gpd)

  # transform to Laplace
  lp_margins <- NULL
  lp_margins[[which.covariate]] <- qlaplace(probs)

  return(lp_margins)
}


## custom function for random sampling of element
runif_func <-function(n, min=1, max=10) sample(min:max, n, replace=T)

fit_dependence_model_ndim_manual <- function(margs, which.target,
                                             which.covariate, dqu,
                                             trace=0){
  #' @export
  #'

  # if which.target is NULL then fit HT with all variables despite covariate
  if (is.null(which.target)){
    varnames <- names(margs)
  } else {varnames = which.target}

  # subset dataset according to all thresholds

  # Loop through the conditions
  df_subset <- margs[[1]]$cov.data
  for (i in 1:length(varnames)) {
    # Subset the dataframe based on the current condition
    df_subset <- df_subset[df_subset[[varnames[i]]] > margs[[varnames[i]]]$threshold,]
  }
  #df_subset <- margs[[1]]$cov.data
  df_natural_scale <- df_subset

  maxd <- NULL
  for (t in varnames){
    if (t != which.covariate){
      # subset datasets to threshold
      X_natural_scale <- df_natural_scale[[which.covariate]]
      Y_natural_scale <- df_natural_scale[[t]]

      plot(X_natural_scale, Y_natural_scale)

      # transform to Laplace scale, use gpd above threshold, use ecdf below
      lp_margins <- transform_to_laplace(df_natural_scale, margs)
      X_laplace_scale <- lp_margins[[which.covariate]]
      Y_laplace_scale <- lp_margins[[t]]

      plot(X_laplace_scale, Y_laplace_scale)

      # get dependency threshold
      if (is.null(dqu)){dependency_threshold<-.5} else {dependency_threshold<-dqu}
      print("dependency_threshold:")
      print(quantile(X_laplace_scale, dependency_threshold))

      # get all data above threshold and fit model
      X_fit_natural <- df_natural_scale[[which.covariate]][X_laplace_scale>quantile(X_laplace_scale, dependency_threshold)]
      Y_fit_natural <- df_natural_scale[[t]][X_laplace_scale>quantile(X_laplace_scale, dependency_threshold)]
      X_fit <- X_laplace_scale[X_laplace_scale>quantile(X_laplace_scale, dependency_threshold)]
      Y_fit <- Y_laplace_scale[X_laplace_scale>quantile(X_laplace_scale, dependency_threshold)]

      plot(X_fit, Y_fit)

      # Nelder Mead does return with initial values for parameters
      o <- optim(c(a=1, b=1, m=1, s=1), x=X_fit, y=Y_fit, HT2004_mse, control = list(trace=trace))
      print("optimized parameters for HT2004 are:")
      print(o$par)

      # return model
      maxd[[which.covariate]][[t]][["dependence"]][["coefficients"]] <- o$par
      maxd[[which.covariate]][[t]][["dependence"]][["data"]][["X"]] <- X_fit
      maxd[[which.covariate]][[t]][["dependence"]][["data"]][[which.covariate]] <- X_fit_natural
      maxd[[which.covariate]][[t]][["dependence"]][["data"]][["Y"]] <- Y_fit
      maxd[[which.covariate]][[t]][["dependence"]][["data"]][[t]] <- Y_fit_natural
      maxd[[which.covariate]][[t]][['margs']] <- margs
    }
  }
  return(maxd)
}

compute_Z_from_residuals <- function(maxd, which.covariate){
  #' @export

  # create Z matrix dataframe with n x (p-1)
  models <- names(maxd[[which.covariate]])
  nr_models <- length(models)
  X <- maxd[[which.covariate]][[1]]$dependence$data$X
  nr_i <- length(X)

  # take X from an arbitrary object as it needs to be the same for all
  Z_M <- array(0, c(nr_i, nr_models))*NA
  for (m in 1:nr_models){
    if (models[m] != 'Z'){
      a <- maxd[[which.covariate]][[models[m]]]$dependence$coefficients[['a']]
      b <- maxd[[which.covariate]][[models[m]]]$dependence$coefficients[['b']]
      Y <- maxd[[which.covariate]][[models[m]]]$dependence$data$Y

      Z_M[,m] = (Y-a*X)/(X^b)
    }
  }
  # name columns of Z matrix
  Zdf <- data.frame(Z_M)
  colnames(Zdf) <- c(names(maxd[[which.covariate]]))
  return(Zdf)
}

make_dependence_model_ndim_manual <- function(margs, which.target=NULL,
                                              which.covariate="hs", dqu=NULL,
                                              trace=0){
  #' @export

  # dependence models
  maxd <- fit_dependence_model_ndim_manual(margs, which.target,
                                           which.covariate, dqu, trace=trace)
  # compute Z
  Zdf <- compute_Z_from_residuals(maxd, which.covariate)
  maxd[[which.covariate]][['Z']] <- Zdf
  return(maxd)
}

make_dependence_model_manual <- function(margs, which.target="tp",
                                         which.covariate="hs", dqu=NULL){
  #' @export

  maxd <- NULL
  # subset datasets to threshold
  threshold_covariate <- margs[[which.covariate]]$threshold
  threshold_target <- margs[[which.target]]$threshold
  df_natural_scale <- margs[[which.covariate]]$cov.data
  df_natural_scale_sub <- subset(df_natural_scale, df_natural_scale[[which.covariate]]>threshold_covariate)
  df_natural_scale_sub <- subset(df_natural_scale_sub, df_natural_scale_sub[[which.target]]>threshold_target)

  X_narural_scale <- df_natural_scale_sub[[which.covariate]]
  Y_narural_scale <- df_natural_scale_sub[[which.target]]

  # transform to Laplace scale, use gpd above threshold, use ecdf below
  lp_margins <- transform_to_laplace(df_natural_scale_sub, margs)
  #lp_margins <- transform_to_laplace(df_natural_scale)
  X_laplace_scale <- lp_margins[[which.covariate]]
  Y_laplace_scale <- lp_margins[[which.target]]

  # get dependency threshold
  if (is.null(dqu)){dependency_threshold<-.5} else {dependency_threshold<-dqu}
  print("dependency_threshold:")
  print(quantile(X_laplace_scale, dependency_threshold))

  # get all data above threshold and fit model
  X_fit_natural <- df_natural_scale_sub[[which.covariate]][X_laplace_scale>quantile(X_laplace_scale, dependency_threshold)]
  Y_fit_natural <- df_natural_scale_sub[[which.target]][X_laplace_scale>quantile(X_laplace_scale, dependency_threshold)]
  X_fit <- X_laplace_scale[X_laplace_scale>quantile(X_laplace_scale, dependency_threshold)]
  Y_fit <- Y_laplace_scale[X_laplace_scale>quantile(X_laplace_scale, dependency_threshold)]

  # Nelder Mead does return with initial values for parameters
  o <- optim(c(a=1, b=1, m=1, s=1), x=X_fit, y=Y_fit, HT2004_mse, control = list(trace=2))
  print("optimized parameters for HT2004 are:")
  print(o$par)

  # compute residuals Z
  Z <- (Y_fit-o$par[1]*X_fit)/X_fit**o$par[2]

  # return model
  maxd[["dependence"]][["coefficients"]] <- o$par
  maxd[["dependence"]][["Z"]] <- Z
  maxd[["dependence"]][["data"]][["X"]] <- X_fit
  maxd[["dependence"]][["data"]][[which.covariate]] <- X_fit_natural
  maxd[["dependence"]][["data"]][["Y"]] <- Y_fit
  maxd[["dependence"]][["data"]][[which.target]] <- Y_fit_natural
  maxd[['margs']] <- margs
  return(maxd)
}

HT2004_mse <- function(params, x, y) {
  alpha <- params[1]
  beta <- params[2]
  mu <- params[3]
  sigma <- params[4]

  if (sigma <= 0) return (1e8)
  s <- x^beta
  m <- alpha * x + mu * s
  s <- sigma * s
  res <- .5 * ((y - m) / s)^2
  res <- res + log(s)
  res <- sum(res)

  return(res)
}

make_dependence_model <- function(margs, which="hs", dqu=.9){
  #' @export
  maxd <- NULL
  maxd <- make_dependence_model_manual(margs, which=which, dqu=dqu)
  return(maxd)
}

simulate_from_HT2004_ndim_manual <- function(maxd,
                                             data_in=NULL,
                                             which.covariate=NULL){
  #' @export

  X <- transform_to_laplace_ndim_preds(data_in, maxd, which.covariate=which.covariate)
  Xtmp <- X[[which.covariate]]
  Xtmp <- Xtmp[is.finite(Xtmp)]
  X[[which.covariate]] <- Xtmp

  sims <- NULL
  sims[["data"]][["simulated"]][['X']] <- X[[which.covariate]]

  varnames <- names(maxd[[which.covariate]])
  for (n in varnames){
    if (n!='Z'){
      # sim Y on laplace scale
      Y_sim <- maxd[[which.covariate]][[n]]$dependence$coefficients[["a"]] * X[[which.covariate]] +
        X[[which.covariate]]^(maxd[[which.covariate]][[n]]$dependence$coefficients[["b"]]) *
        maxd[[which.covariate]][['Z']][[n]][runif_func(
          length(X[[which.covariate]]), min=1,
          max=length(maxd[[which.covariate]][['Z']][[n]]))]
      # return simulated samples
      print(length(Y_sim))
      print(length(X[[which.covariate]]))
      print(length(maxd[[which.covariate]][['Z']][[n]]))
      sims[["data"]][["simulated"]][['Y']][[n]] <- Y_sim
    }
  }
  # retain only finite values in X and Y's
  # tmp_df for subsetting
  #tmp_df <- sims$data$simulated$Y
  #tmp_df[['X']] <- sims$data$simulated$X
  #tmp_df <- data.frame(tmp_df)
  #for (n in varnames){
  #  tmp_df <- subset(tmp_df, is.finite(n))
  #}
  #sims$data$simulated$Y <- tmp_df[['Y']]
  #sims$data$simulated$X <- tmp_df[['X']]
  return(sims)
}

simulate_from_HT2004_manual <- function(maxd,
                                        data_in=NULL, data_in_lp=NULL,
                                        names_in=NULL){
  #' @export

  if (!is.null(data_in)){
    X <- transform_to_laplace(storm_peaks, maxd$margs,
                              names_in=names_in, maxd=maxd)
    X <- X[[names_in]]
  }
  else if (!is.null(data_in_lp)){X <- data_in_lp}
  else {X <- maxd$dependence$data$X}

  sims <- NULL

  # sim Y on laplace scale
  Y_sim <- maxd$dependence$coefficients[["a"]] * X +
    X^(maxd$dependence$coefficients[["b"]]) *
    maxd$dependence$Z[runif_func(length(X), min=1, max=length(maxd$dependence$Z))]

  # return simulated samples
  sims[["data"]][["simulated"]][["Y"]] <- Y_sim
  sims[["data"]][["simulated"]][["X"]] <- X
  return(sims)
}

transform_to_natural <- function(sims, maxd, which.target){
  #' @export

  p_lp <- plaplace(sims$data$simulated$Y,
                   maxd$dependence$coefficients["m"],
                   maxd$dependence$coefficients["s"])

  p_lp <- p_lp[p_lp<1]
  Y_natural <- qevd(p_lp,
                    scale = exp(maxd$margs[[which.target]]$results$par[1]),
                    shape = maxd$margs[[which.target]]$results$par[2],
                    threshold = maxd$margs[[which.target]]$threshold,
                    type = c("GP"))
  return(Y_natural)
}

transform_to_natural_ndim <- function(sims, maxd, which.target, which.covariate){
  #' @export

  if (which.target == which.covariate){
    p_lp <- plaplace(sims$data$simulated$X)
    p_lp <- p_lp[p_lp<1]
    Y_natural <- qevd(p_lp,
                      scale = exp(maxd[[which.covariate]][[1]]$margs[[which.covariate]]$results$par[1]),
                      shape = maxd[[which.covariate]][[1]]$margs[[which.covariate]]$results$par[2],
                      threshold = maxd[[which.covariate]][[1]]$margs[[which.covariate]]$threshold,
                      type = c("GP"))
  } else {
    p_lp <- plaplace(sims$data$simulated$Y[[which.target]])
    p_lp <- p_lp[p_lp<1]
    Y_natural <- qevd(p_lp,
                      scale = exp(maxd[[which.covariate]][[which.target]]$margs[[which.target]]$results$par[1]),
                      shape = maxd[[which.covariate]][[which.target]]$margs[[which.target]]$results$par[2],
                      threshold = maxd[[which.covariate]][[which.target]]$margs[[which.target]]$threshold,
                      type = c("GP"))
  }
  return(Y_natural)
}


simulate_from_HT2004 <- function(maxd, pqu, nsim){
  #' @export
  #'
  preds <- predict(maxd, which='hs', pqu = pqu, nsim=nsim, trace=10)
  return(preds)
}

compute_probs <- function(margs_thr, margs_gpd, dfin = NULL, list_var = NULL,
                          thr_str = 'thr', thr_ecdf_gpd_transistion_margin = 0.01){
  #' @export
  #'
  probs_ecdf <- NULL
  probs_gpd <- NULL
  probs <- NULL
  if (is.null(list_var)){
    list_var <- names(margs_gpd)
  }
  for (n in list_var){
    if (is.null(dfin)){
      dfin <- margs_gpd[[n]]$data[[n]]
    }

    # ecdf
    Fn_ecdf <- ecdf(dfin[[n]])
    probs_ecdf[[n]] <- Fn_ecdf(dfin[[n]])

    # threshold
    threshold <- predict(margs_thr[[n]], newdata = dfin, type= "response")$location

    # gpd
    gpd_params <- predict(margs_gpd[[n]], newdata = dfin, type= "response")
    scales <- as.vector(gpd_params$scale)
    shapes <- as.vector(gpd_params$shape)

    tmp <- array(0, c(length(shapes)))*NA
    for (i in 1:length(tmp)) {
      tmp[i] <- margs_thr[[n]]$tau + (1-margs_thr[[n]]$tau) *
        pgpd(dfin[[n]][i], mu=threshold[i], sigma=scales[i],
             xi=shapes[i], lower.tail = TRUE)
    }
    probs_gpd[[n]] <- tmp
    rm(tmp)

    # combine probs
    tmp <- probs_gpd[[n]]
    tmp[tmp<(margs_thr[[n]]$tau+thr_ecdf_gpd_transistion_margin)] <- probs_ecdf[[n]][tmp<(margs_thr[[n]]$tau+thr_ecdf_gpd_transistion_margin)]
    probs[[n]] <- unlist(tmp)
    rm(tmp)
  }
  return (list('probs'=probs, 'probs_ecdf'=probs_ecdf, 'probs_gpd'=probs_gpd))
}

compute_probs_single <- function(marg_thr, marg_gpd, dfin = NULL, varstr = NULL,
                                 thr_str = 'thr', thr_ecdf_gpd_transistion_margin = 0.01){
  #' @export
  #'
  probs_ecdf <- NULL
  probs_gpd <- NULL
  probs <- NULL

  if (is.null(dfin)){
    dfin <- marg_gpd$data[[varstr]]
  }

  # ecdf
  Fn_ecdf <- ecdf(dfin[[varstr]])
  probs_ecdf <- Fn_ecdf(dfin[[varstr]])

  # threshold
  threshold <- predict(marg_thr, newdata = dfin, type= "response")$location

  # gpd
  gpd_params <- predict(marg_gpd, newdata = dfin, type= "response")
  scales <- as.vector(gpd_params$scale)
  shapes <- as.vector(gpd_params$shape)

  tmp <- array(0, c(length(shapes)))*NA
  for (i in 1:length(tmp)) {
    tmp[i] <- marg_thr$tau + (1-marg_thr$tau) *
      pgpd(dfin[[varstr]][i], mu=threshold[i], sigma=scales[i],
           xi=shapes[i], lower.tail = TRUE)
  }
  probs_gpd <- tmp
  rm(tmp)

  # combine probs
  tmp <- probs_gpd
  tmp[tmp<(marg_thr$tau+thr_ecdf_gpd_transistion_margin)] <- probs_ecdf[tmp<(marg_thr$tau+thr_ecdf_gpd_transistion_margin)]
  probs <- unlist(tmp)
  rm(tmp)

  return (list('probs'=probs, 'probs_ecdf'=probs_ecdf, 'probs_gpd'=probs_gpd))
}

compute_lp_margs <- function(probs){
  #' @export
  #'
  lp_margs <- NULL
  for (n in names(probs)){
    lp_margs[[n]] <- qlaplace(probs[[n]])
  }
  return(lp_margs)
}

fit_HT2004 <- function(lp_margs, thr){
  #' @export
  #'

  maxd <- NULL
  for (i in 1:length(names(lp_margs))){
    list_of_vars <- names(lp_margs)[-i]
    Zs <- NULL

    X_fit_str <- names(lp_margs)[i]
    X_fit_tmp <- lp_margs[[X_fit_str]]
    for (n in list_of_vars){
      Y_fit_tmp <- lp_margs[[n]]
      Y_fit <- Y_fit_tmp[X_fit_tmp>quantile(X_fit_tmp,thr)]
      X_fit <- X_fit_tmp[X_fit_tmp>quantile(X_fit_tmp,thr)]
      # upper limit on a and b above 1 to allow for a chance to hit 1 assuming noisy data
      o <- optim(c(a=.5, b=.5, m=1, s=1), x=X_fit, y=Y_fit, HT2004_mse, method='L-BFGS-B', lower = c(0,-100,-100,0.001), upper = c(1.2,1.2,100,100))
      maxd[[names(lp_margs)[i]]][[n]][['params']] <- o
      maxd[[names(lp_margs)[i]]][[n]][['Z']] <- (Y_fit-o$par[1]*X_fit)/X_fit**o$par[2]
      maxd[[names(lp_margs)[i]]][[n]][['X_fit']] <- X_fit
      maxd[[names(lp_margs)[i]]][[n]][['Y_fit']] <- Y_fit
      maxd[[names(lp_margs)[i]]][[n]][['X_all']] <- X_fit_tmp
      maxd[[names(lp_margs)[i]]][[n]][['Y_all']] <- Y_fit_tmp
      maxd[[names(lp_margs)[i]]][[n]][['thr']] <- thr
    }
  }

  return (maxd)
}


compute_HT2004_parameter_bstrp <- function(dfin,
                                           model_fml_thr, model_fml_occ, model_fml_gpd,
                                           maxd_thr_lst, extr_thr,
                                           thr_str='thr', list_var=NULL,
                                           margs_thr_orig=NULL,
                                           margs_gpd_orig=NULL){
  #' @export
  #'
  if (is.null(list_var)){
    list_var <- names(model_fml_thr)
  }
  print(c('Considered variables:',list_var))
  maxds_vars <- NULL
  for (n in list_var){
    maxds_thr <- NULL
    for (t in 1:length(maxd_thr_lst)){
      maxds_bstr <- NULL
      for (i in 1:length(dfin)){
        ### Fitting marginals ###
        print(c("variable:",n,"threshold nr:",t,"bootstrap nr:",i))
        print("fit threshold model")
        margs_thr <- fit_marginal_models_thr(dfin = dfin[[i]],
                                             list_var = list_var,
                                             thr = extr_thr,
                                             model_fml = model_fml_thr,
                                             m_params = margs_thr_orig)

        print("subset to pots above threshold")
        data_sub <- subset_df(margs_thr, thr_str='thr', exc_str='exc')

        print("fit GPD model")
        margs_gpd <- fit_marginal_models_gpd(dfin = data_sub,
                                             list_var = list_var,
                                             model_fml = model_fml_gpd,
                                             m_params = margs_gpd_orig)

        print("fit occ model")
        #margs_occ <- fit_marginal_models_occ(dfin = data_sub,
        #                                     list_var = list_var,
        #                                     model_fml = model_fml_occ,
        #                                     nr_of_years = nr_of_years)
        margs_occ <- fit_marginal_models_occ_v2(dfin = data_sub,
                                                list_var = list_var,
                                                model_fml = model_fml_occ,
                                                nr_of_years = nr_of_years)

        probs <- compute_probs(margs_thr, margs_gpd, dfin = dfin[[i]],
                               list_var = list_var, thr_str = 'thr')[['probs']]

        lp_margs <- compute_lp_margs(probs)

        maxd <- fit_HT2004(lp_margs, thr = maxd_thr_lst[t])

        maxds_bstr[[i]] <- maxd
      }
      maxds_thr[[t]] <- maxds_bstr
    }
    maxds_vars[[n]] <- maxds_thr
  }
  return(maxds_vars)
}

unfold_maxd_params_bstrp <- function(maxds_vars, X_var_str, Y_var_str,
                                     ulim=.01,llim=.99,cntr=.5){
  #' @export
  #'
  alphas_cntr <- array(0, c(length(maxds_vars[[X_var_str]])))*NA
  alphas_llim <- array(0, c(length(maxds_vars[[X_var_str]])))*NA
  alphas_ulim <- array(0, c(length(maxds_vars[[X_var_str]])))*NA

  betas_cntr <- array(0, c(length(maxds_vars[[X_var_str]])))*NA
  betas_llim <- array(0, c(length(maxds_vars[[X_var_str]])))*NA
  betas_ulim <- array(0, c(length(maxds_vars[[X_var_str]])))*NA

  for (t in 1:length(maxds_vars[[X_var_str]])){
    alphas <- array(0, c(length(maxds_vars[[X_var_str]][[1]])))*NA
    betas <- array(0, c(length(maxds_vars[[X_var_str]][[1]])))*NA
    for (i in 1:length(alphas)){
      alphas[i] <- as.numeric(maxds_vars[[X_var_str]][[t]][[i]][[X_var_str]][[Y_var_str]]$params$par[1])
    }
    for (i in 1:length(betas)){
      betas[i] <- as.numeric(maxds_vars[[X_var_str]][[t]][[i]][[X_var_str]][[Y_var_str]]$params$par[2])
    }

    alphas_cntr[t] <- quantile(alphas, cntr)
    alphas_llim[t] <- quantile(alphas, llim)
    alphas_ulim[t] <- quantile(alphas, ulim)

    betas_cntr[t] <- quantile(betas, cntr)
    betas_llim[t] <- quantile(betas, llim)
    betas_ulim[t] <- quantile(betas, ulim)
  }
  return(list('alphas_llim'=alphas_llim, "alphas_cntr"=alphas_cntr, "alphas_ulim"=alphas_ulim,
              'betas_llim'= betas_llim, "betas_cntr"=betas_cntr, "betas_ulim"=betas_ulim))
}

unfold_maxd_params_bstrp_v2 <- function(maxds_vars, X_var_str, Y_var_str,
                                        ulim=.01,llim=.99,cntr=.5){
  #' @export
  #' models_maxds_thr[[t_maxd]][[b]]$hs$U10
  #'
  alphas_cntr <- array(0, c(length(maxds_vars)))*NA
  alphas_llim <- array(0, c(length(maxds_vars)))*NA
  alphas_ulim <- array(0, c(length(maxds_vars)))*NA

  betas_cntr <- array(0, c(length(maxds_vars)))*NA
  betas_llim <- array(0, c(length(maxds_vars)))*NA
  betas_ulim <- array(0, c(length(maxds_vars)))*NA

  mus_cntr <- array(0, c(length(maxds_vars)))*NA
  mus_llim <- array(0, c(length(maxds_vars)))*NA
  mus_ulim <- array(0, c(length(maxds_vars)))*NA

  sigmas_cntr <- array(0, c(length(maxds_vars)))*NA
  sigmas_llim <- array(0, c(length(maxds_vars)))*NA
  sigmas_ulim <- array(0, c(length(maxds_vars)))*NA

  for (t in 1:length(maxds_vars)){
    alphas <- array(0, c(length(maxds_vars[[t]])))*NA
    betas <- array(0, c(length(maxds_vars[[t]])))*NA
    mus <- array(0, c(length(maxds_vars[[t]])))*NA
    sigmas <- array(0, c(length(maxds_vars[[t]])))*NA
    for (b in 1:length(alphas)){
      alphas[b] <- as.numeric(maxds_vars[[t]][[b]][[X_var_str]][[Y_var_str]]$params$par[1])
      betas[b] <- as.numeric(maxds_vars[[t]][[b]][[X_var_str]][[Y_var_str]]$params$par[2])
      mus[b] <- as.numeric(maxds_vars[[t]][[b]][[X_var_str]][[Y_var_str]]$params$par[3])
      sigmas[b] <- as.numeric(maxds_vars[[t]][[b]][[X_var_str]][[Y_var_str]]$params$par[4])
    }

    alphas_cntr[t] <- quantile(alphas, cntr)
    alphas_llim[t] <- quantile(alphas, llim)
    alphas_ulim[t] <- quantile(alphas, ulim)

    betas_cntr[t] <- quantile(betas, cntr)
    betas_llim[t] <- quantile(betas, llim)
    betas_ulim[t] <- quantile(betas, ulim)

    mus_cntr[t] <- quantile(mus, cntr)
    mus_llim[t] <- quantile(mus, llim)
    mus_ulim[t] <- quantile(mus, ulim)

    sigmas_cntr[t] <- quantile(sigmas, cntr)
    sigmas_llim[t] <- quantile(sigmas, llim)
    sigmas_ulim[t] <- quantile(sigmas, ulim)
  }

  return(list('alpha_llim'=alphas_llim, "alpha_cntr"=alphas_cntr, "alpha_ulim"=alphas_ulim,
              'beta_llim'= betas_llim, "beta_cntr"=betas_cntr, "beta_ulim"=betas_ulim,
              'mu_llim'= mus_llim, "mu_cntr"=mus_cntr, "mu_ulim"=mus_ulim,
              'sigma_llim'= sigmas_llim, "sigma_cntr"=sigmas_cntr, "sigma_ulim"=sigmas_ulim))
}

fit_margs_bstrp <- function(dfin,
                            model_fml_thr, model_fml_occ, model_fml_gpd,
                            extr_thr, nr_of_years, thr_str='thr', list_var=NULL,
                            margs_thr_orig=NULL,
                            margs_gpd_orig=NULL,
                            nbstrp=NULL,
                            nquad=40){
  #' @export
  #'
  margs <- NULL

  if (is.null(list_var)){
    list_var <- names(model_fml_thr)
  }
  print(c('Considered variables:', list_var))

  if (is.null(nbstrp)){
    nbstrp <- length(dfin)
  }
  print(c('number of bootstraps is:', nbstrp))

  margs <- NULL
  for (i in 1:nbstrp){

    ### Fitting marginals ###

    # bootstrap extremal threshold within bounds
    if (length(extr_thr)>1){
      t <- runif(1, min = extr_thr[1], max = extr_thr[2])
    } else {t <- extr_thr}

    print(c("bootstrap nr:",i))
    print("fit threshold model")
    margs_thr <- fit_marginal_models_thr(dfin = dfin[[i]],
                                         list_var = list_var,
                                         thr = t,
                                         model_fml = model_fml_thr,
                                         m_params = margs_thr_orig)

    print("subset to pots above threshold")
    data_sub <- subset_df(margs_thr, thr_str='thr', exc_str='exc')

    print("fit GPD model")
    margs_gpd <- fit_marginal_models_gpd(dfin = data_sub,
                                         list_var = list_var,
                                         model_fml = model_fml_gpd,
                                         m_params = margs_gpd_orig)

    print("fit occ model")
    #margs_occ <- fit_marginal_models_occ(dfin = data_sub,
    #                                     list_var = list_var,
    #                                     model_fml = model_fml_occ,
    #                                     nr_of_years = nr_of_years)
    margs_occ <- fit_marginal_models_occ_v2(dfin = data_sub,
                                            list_var = list_var,
                                            model_fml = model_fml_occ,
                                            nr_of_years = nr_of_years,
                                            nquad = nquad)

    probs <- compute_probs(margs_thr, margs_gpd, dfin = dfin[[i]],
                           list_var = list_var, thr_str = 'thr')[['probs']]

    margs_tmp <- NULL
    margs_tmp[['thr']] <- margs_thr
    margs_tmp[['gpd']] <- margs_gpd
    margs_tmp[['occ']] <- margs_occ
    margs_tmp[['probs']] <- probs

    margs[[i]] <- margs_tmp
  }
  return(margs)
}

fit_models_bstrp <- function(dfin,
                             model_fml_thr, model_fml_occ, model_fml_gpd,
                             maxd_thr_lst, extr_thr,
                             thr_str='thr', list_var=NULL,
                             margs_thr_orig=NULL,
                             margs_gpd_orig=NULL){
  #' @export
  #'
  if (is.null(list_var)){
    list_var <- names(model_fml_thr)
  }
  print(c('Considered variables:',list_var))
  maxds_vars <- NULL
  margs_vars <- NULL
  for (n in list_var){
    maxds_thr <- NULL
    margs_thr <- NULL
    for (t in 1:length(maxd_thr_lst)){
      maxds_bstr <- NULL
      margs_bstr <- NULL
      for (i in 1:length(dfin)){
        ### Fitting marginals ###
        print(c("variable:",n,"threshold nr:",t,"bootstrap nr:",i))
        print("fit threshold model")
        margs_thr <- fit_marginal_models_thr(dfin = dfin[[i]],
                                             list_var = list_var,
                                             thr = extr_thr,
                                             model_fml = model_fml_thr,
                                             m_params = margs_thr_orig)

        print("subset to pots above threshold")
        data_sub <- subset_df(margs_thr, thr_str='thr', exc_str='exc')

        print("fit GPD model")
        margs_gpd <- fit_marginal_models_gpd(dfin = data_sub,
                                             list_var = list_var,
                                             model_fml = model_fml_gpd,
                                             m_params = margs_gpd_orig)

        print("fit occ model")
        #margs_occ <- fit_marginal_models_occ(dfin = data_sub,
        #                                     list_var = list_var,
        #                                     model_fml = model_fml_occ,
        #                                     nr_of_years = nr_of_years)
        margs_occ <- fit_marginal_models_occ_v2(dfin = data_sub,
                                                list_var = list_var,
                                                model_fml = model_fml_occ,
                                                nr_of_years = nr_of_years)

        margs <- NULL
        margs[['thr']] <- margs_thr
        margs[['gpd']] <- margs_gpd
        margs[['occ']] <- margs_occ

        probs <- compute_probs(margs_thr, margs_gpd, dfin = dfin[[i]],
                               list_var = list_var, thr_str = 'thr')[['probs']]

        lp_margs <- compute_lp_margs(probs)

        maxd <- fit_HT2004(lp_margs, thr = maxd_thr_lst[t])

        maxds_bstr[[i]] <- maxd
        margs_bstr[[i]] <- margs
      }
      maxds_thr[[t]] <- maxds_bstr
      margs_thr[[t]] <- margs_bstr
    }
    maxds_vars[[n]] <- maxds_thr
    margs_vars[[n]] <- margs_thr
  }
  return(list("maxds"=maxds_vars,"margs"=margs_vars))
}

fit_maxds_bstrp <- function(dfin, margs,
                            maxd_thr, nbstrp=NULL){
  #' @export
  #'

  if (is.null(nbstrp)){
    nbstrp <- length(dfin)
  }

  if (length(maxd_thr)>1){
    t <- runif(1, min = maxd_thr[1], max = maxd_thr[2])
  } else {t <- maxd_thr}

  maxds <- NULL
  for (b in 1:nbstrp){
    ### Fitting marginals ###
    print(c("bootstrap nr:",b))
    lp_margs <- compute_lp_margs(models_margs[[b]]$probs)
    maxd <- fit_HT2004(lp_margs, thr = t)
    maxds[[b]] <- maxd
  }

  return(maxds)
}

predict_maxd <- function(maxd, varstr_X, varstr_Y, nsim){
  tmp_alpha <- maxd[[varstr_X]][[varstr_Y]]$params$par[1]
  tmp_beta <- maxd[[varstr_X]][[varstr_Y]]$params$par[2]
  tmp_Z <- maxd[[varstr_X]][[varstr_Y]][['Z']]
  lp_X <- maxd[[varstr_X]][[varstr_Y]]$X_fit

  Xsim <- array(0, nsim)*NA
  Ysim <- array(0, nsim)*NA
  for (i in 1:nsim){
    Xsim[i] <- lp_X[runif_func(1, min=1, max=length(lp_X))]
    Ysim[i] <- tmp_alpha * Xsim[i] + Xsim[i]**(tmp_beta)* tmp_Z[runif_func(1, min=1, max=length(tmp_Z))]
  }

  return(list("Xsim"=Xsim,"Ysim"=Ysim))
}

predict_from_HT2004_models <- function(margs, maxds, nr_of_years, RP,
                                       var_lst, nbstrp, nmc=1,
                                       preds=NULL, condition=NULL){
  #' @export

  if (is.null(var_lst)){
    var_lst <- names(maxds[[1]])
  }

  lp_Y_bstrp <- NULL
  lp_X_bstrp <- NULL
  X_bstrp <- NULL
  if (is.null(preds)){
    # predict values from non-stationary marginal model if not supplied
    preds <- predict_margs(margs, nr_of_years, RP,
                           nmc = nmc,
                           var_lst = var_lst,
                           nbstrp = nbstrp,
                           condition = condition)
  }
  for (b in 1:nbstrp){
    lp_Y_lst <- NULL
    lp_X_lst <- NULL
    X_lst <- NULL
    for (n in 1:length(var_lst)){
      lp_Y_lst_tmp <- NULL
      tmplst <- var_lst[-n]
      lp_X <- qlaplace(preds[[b]][[var_lst[n]]]$prob)
      X_sim <- preds[[b]][[var_lst[n]]]$maxval

      print("predict from HT2004")
      for (m in tmplst){
        tmp_alpha <- maxds[[b]][[var_lst[n]]][[m]]$params$par[1]
        tmp_beta <- maxds[[b]][[var_lst[n]]][[m]]$params$par[2]
        tmp_Z <- maxds[[b]][[var_lst[n]]][[m]][['Z']]
        lp_Y <- tmp_alpha * lp_X + lp_X**(tmp_beta) * tmp_Z[runif_func(length(lp_X), min=1, max=length(tmp_Z))]
        lp_Y_lst_tmp[[m]] <- as.numeric(lp_Y)
      }
      lp_Y_lst[[var_lst[n]]] <- lp_Y_lst_tmp
      lp_X_lst[[var_lst[n]]] <- lp_X
      X_lst[[var_lst[n]]] <- X_sim
    }
    lp_Y_bstrp[[b]] <- lp_Y_lst
    lp_X_bstrp[[b]] <- lp_X_lst
    X_bstrp[[b]] <- X_lst
  }
  return(list("lp_Y" = lp_Y_bstrp, "lp_X" = lp_X_bstrp, "X" = X_bstrp))
}

convert_HT2004_preds_to_original_space <- function(margs, maxds,
                                                   maxd_preds, marg_preds,
                                                   var_lst = NULL,
                                                   nbstrp = NULL){
  #' @export

  if (is.null(var_lst)){
    var_lst <- names(maxds[[1]])
  }

  if (is.null(nbstrp)){
    nbstrp <- length(marg_preds)
  }

  Y_bstrp <- NULL
  for (b in 1:nbstrp){
    Y_lst <- NULL
    X_lst <- NULL
    for (n in 1:length(var_lst)){
      Y_lst_tmp <- NULL
      tmplst <- var_lst[-n]
      for (m in tmplst){
        lp_prob <- plaplace(maxd_preds$lp_Y[[b]][[var_lst[n]]][[m]])
        thr_sim <- marg_preds[[1]][[m]]$thr
        scale_sim <- marg_preds[[1]][[m]]$scale
        shape_sim <- marg_preds[[1]][[m]]$shape
        Y_lst_tmp[[m]] <- qgpd(lp_prob, mu = thr_sim, sigma = scale_sim,
                               xi = shape_sim, lower.tail = TRUE, log.p = FALSE)
      }
      Y_lst[[var_lst[n]]] <- Y_lst_tmp
    }
    Y_bstrp[[b]] <- Y_lst
  }
  return(Y_bstrp)
}

predict_models_bstrp <- function(margs, maxds, nr_of_years, RP, var_lst=NULL, nbstrp=NULL){
  #' @export
  #'

  if (is.null(var_lst)){
    var_lst <- names(models$maxds)
  }

  if (is.null(nbstrp)){
    nbstrp <- length(margs)
  }

  max_vals <- NULL
  max_lp_X <- NULL
  max_lp_Y <- NULL
  accepted <- NULL
  rejected <- NULL
  for (b in 1:nbstrp){ # number of bootstraps
    for (n in 1:length(var_lst)){

      # predict from PP: produce storms occurrences at correct rate
      print("produce storms")

      tmp <- margs[[b]]$thr[[n]]$data[[var_lst[n]]] - fitted(margs[[b]]$thr[[var_lst[n]]])$location
      nr_of_events <- length(tmp[tmp>0])
      rm(tmp)

      df_storm_cov <- produce_storm_occurrences(nr_of_events, RP, nr_of_years, margs[[b]]$occ[[var_lst[n]]])

      # predict from ALD
      print("predict threshold")
      thr <- predict(margs[[b]]$thr[[var_lst[n]]], newdata = df_storm_cov, type = "response")$location

      # predict from GPD
      print("predict exceedences")
      gpd_param_sims <- predict(margs[[b]]$gpd[[var_lst[n]]], newdata = df_storm_cov, type= "response")
      scales <- gpd_param_sims$scale
      shapes <- gpd_param_sims$shape
      gpd_sims <- revd(length(scales), scale = scales, shape = shapes, threshold = thr, type="GP")

      max_idx <- which(gpd_sims==max(gpd_sims))

      # Save maximum
      print("save max and max_idx")
      max_val <- gpd_sims[max_idx]
      max_scale <- scales[max_idx]
      max_shape <- shapes[max_idx]
      thr_max <- thr[max_idx]
      df_pred <- df_storm_cov[max_idx,]

      # reduce tmplst
      tmplst <- var_lst[-n]
      print(tmplst)

      # predict from HT2004 for remaining response variables
      print("predict from HT2004")
      prob_X <- margs[[b]]$thr[[var_lst[n]]]$tau + (1-margs[[b]]$thr[[var_lst[n]]]$tau) *
        pgpd(max_val, mu=thr_max, sigma=max_scale,
             xi=max_shape, lower.tail = TRUE)
      lp_X <- qlaplace(prob_X)

      lp_Y <- NULL
      Y_sim <- NULL
      for (m in tmplst){
        tmp_alpha <- maxds[[b]][[var_lst[n]]][[m]]$params$par[1]
        tmp_beta <- maxds[[b]][[var_lst[n]]][[m]]$params$par[2]
        tmp_Z <- maxds[[b]][[var_lst[n]]][[m]][['Z']]
        lp_Y[[m]] <- tmp_alpha * lp_X
        + lp_X**(tmp_beta) * tmp_Z[runif_func(1, min=1, max=length(tmp_Z))]
        lp_prob <- plaplace(lp_Y[[m]])
        # predict threshold for max
        thr_sim <- predict(margs[[b]]$thr[[m]], newdata = df_pred, type= "response")$location
        # predict gpd params
        gpd_param_sims <- predict(margs[[b]]$gpd[[m]], newdata = df_pred, type= "response")
        scale_sim <- gpd_param_sims$scale
        shape_sim <- gpd_param_sims$shape
        Y_sim[[m]] <- qgpd(lp_prob, mu = thr_sim, sigma = scale_sim, xi = shape_sim, lower.tail = TRUE, log.p = FALSE)

        max_lp_Y[[var_lst[n]]][[m]][[b]] <- as.numeric(lp_Y[[m]])
        max_lp_X[[var_lst[n]]][[m]][[b]] <- as.numeric(lp_X)

        max_vals[[var_lst[n]]][[m]][[b]] <- as.numeric(Y_sim[[m]])
      }

      max_vals[[var_lst[n]]][[var_lst[n]]][[b]] <- max_val
      #if (lp_X > max(unlist(lp_Y))){
      #  for (m in names(Y_sim)){
      #    max_vals[[var_lst[n]]][[m]][[b]] <- as.numeric(Y_sim[[m]])
      #  }
      #}
    }
  }

  return(list('max_vals'=max_vals, 'max_lp_Y'=max_lp_Y, 'max_lp_X'=max_lp_X))
}

predict_marg <- function(margs, nr_of_years, RP, varstr = NULL, nmc = 1,
                         condition = NULL){
  #' @export
  #'

  df_max <- data.frame(matrix(ncol = 5, nrow = nmc))
  colnames(df_max) <- c("maxval", "scale", "shape", "thr", "prob")

  predictor_names <- margs$thr[[varstr]]$predictor.names
  df_covs <- data.frame(matrix(ncol = length(predictor_names),
                               nrow = nmc))
  colnames(df_covs) <- predictor_names

  # predict from PP: produce storms occurrences at correct rate
  print("produce storms")

  tmp <- margs$thr[[varstr]]$data[[varstr]] - fitted(margs$thr[[varstr]])$location
  nr_of_events <- length(tmp[tmp>0])
  rm(tmp)

  df_storm_cov <- produce_storm_occurrences(nr_of_events, RP, nr_of_years,
                                            margs$occ[[varstr]],
                                            condition = condition)

  # predict from ALD
  print("predict threshold")
  thr <- predict(margs$thr[[varstr]], newdata = df_storm_cov,
                 type = "response")$location

  # predict from GPD
  print("predict exceedences")
  gpd_param_sims <- predict(margs$gpd[[varstr]],
                            newdata = df_storm_cov, type= "response")
  scales <- gpd_param_sims$scale
  shapes <- gpd_param_sims$shape

  for (i in 1:nmc){
    gpd_sims <- revd(length(scales), scale = scales, shape = shapes,
                     threshold = thr, type="GP")
    max_idx <- which(gpd_sims==max(gpd_sims))

    # Save maximum
    max_val <- gpd_sims[max_idx]
    max_scale <- scales[max_idx]
    max_shape <- shapes[max_idx]
    thr_max <- thr[max_idx]
    df_pred_cov <- df_storm_cov[max_idx,]

    gpd_prob <- pevd(max_val, scale = max_scale, shape = max_shape,
                     threshold = thr_max, type = "GP", lower.tail = TRUE)

    # store in output field
    df_max$maxval[i] <- max_val
    df_max$scale[i] <- max_scale
    df_max$shape[i] <- max_shape
    df_max$thr[i] <- thr_max
    df_max$prob[i] <- gpd_prob

    # store covariates
    for (n in predictor_names){
      df_covs[[n]][[i]] <- df_pred_cov[[n]]
    }
  }

  # combine dfs
  df_out <- cbind(df_max, df_covs)

  return(df_out)
}

predict_margs <- function(margs, nr_of_years, RP, nmc = 1,
                          var_lst=NULL, nbstrp=NULL, condition=NULL){
  #' @export
  #'

  if (is.null(var_lst)){
    var_lst <- names(margs[[1]]$thr)
  }

  if (is.null(nbstrp)){
    nbstrp <- length(margs)
  }

  preds_lst <- NULL
  for (b in 1:nbstrp){ # number of bootstraps
    print(c("number of boostraps (predict_margs):", b))
    preds_tmp_lst <- NULL
    for (n in 1:length(var_lst)){
      preds <- predict_marg(margs[[b]], nr_of_years, RP,
                            varstr = var_lst[n], nmc = nmc,
                            condition = condition)
      preds_tmp_lst[[var_lst[n]]] <- preds
    }
    preds_lst[[b]] <- preds_tmp_lst
  }
  return(preds_lst)
}

predict_margs_bstrp <- function(margs, nr_of_years,
                                RP, var_lst=NULL, nbstrp=NULL,
                                condition=NULL){
  #' @export
  #'

  if (is.null(var_lst)){
    var_lst <- names(models$maxds)
  }

  if (is.null(nbstrp)){
    nbstrp <- length(margs)
  }

  df_max <- data.frame(matrix(ncol = 5, nrow = nbstrp))
  colnames(df_max) <- c("maxval", "scale", "shape", "thr", "prob")
  df_max_vals <- NULL
  for (n in var_lst){
    df_max_vals[[n]] <- df_max
  }

  df_covs <- data.frame(matrix(ncol = length(margs[[1]]$thr[[1]]$predictor.names),
                               nrow = nbstrp))
  colnames(df_covs) <- margs[[1]]$thr[[1]]$predictor.names
  df_max_covs <- NULL
  for (n in var_lst){
    df_max_covs[[n]] <- df_covs
  }

  max_vals <- NULL
  for (b in 1:nbstrp){ # number of bootstraps
    for (n in 1:length(var_lst)){

      # predict from PP: produce storms occurrences at correct rate
      print("produce storms")

      tmp <- margs[[b]]$thr[[n]]$data[[var_lst[n]]] - fitted(margs[[b]]$thr[[var_lst[n]]])$location
      nr_of_events <- length(tmp[tmp>0])
      rm(tmp)

      df_storm_cov <- produce_storm_occurrences(nr_of_events, RP, nr_of_years,
                                                margs[[b]]$occ[[var_lst[n]]],
                                                condition = condition)

      # predict from ALD
      print("predict threshold")
      thr <- predict(margs[[b]]$thr[[var_lst[n]]], newdata = df_storm_cov,
                     type = "response")$location

      # predict from GPD
      print("predict exceedences")
      gpd_param_sims <- predict(margs[[b]]$gpd[[var_lst[n]]],
                                newdata = df_storm_cov, type= "response")
      scales <- gpd_param_sims$scale
      shapes <- gpd_param_sims$shape
      gpd_sims <- revd(length(scales), scale = scales, shape = shapes,
                       threshold = thr, type="GP")
      max_idx <- which(gpd_sims==max(gpd_sims))

      # Save maximum
      print("save max and max_idx")
      max_val <- gpd_sims[max_idx]
      max_scale <- scales[max_idx]
      max_shape <- shapes[max_idx]
      thr_max <- thr[max_idx]
      df_pred <- df_storm_cov[max_idx,]

      gpd_prob <- pevd(max_val, scale = max_scale, shape = max_shape, threshold = thr_max,
                       type = "GP", lower.tail = TRUE)

      # store in output field
      df_max_vals[[var_lst[n]]][b,1] <- max_val
      df_max_vals[[var_lst[n]]][b,2] <- max_scale
      df_max_vals[[var_lst[n]]][b,3] <- max_shape
      df_max_vals[[var_lst[n]]][b,4] <- thr_max
      df_max_vals[[var_lst[n]]][b,5] <- gpd_prob

      # store covariates
      df_max_covs[[var_lst[n]]][b,] <- df_pred[1,]
    }
  }
  return(list("max_vals"=df_max_vals, "max_covs"=df_max_covs))
}

retrieve_valid_HT_samples <- function(preds_HT, preds_HT_lp, preds_margs,
                                      varstr_X, varstr_Y, nbstrp=NULL){
  #' @export

  # color the valid ones and add
  # preds_HT <- HT_sims
  # preds_HT_lp <- preds_HT
  # preds_margs <- preds_margs_lst

  if (is.null(nbstrp)){
    nbstrp <- length(preds_margs)
  }

  X_lp_valid <- NULL
  Y_lp_valid <- NULL
  X_lp_invalid <- NULL
  Y_lp_invalid <- NULL

  Y_valid <- NULL
  Y_invalid <- NULL
  X_valid <- NULL
  X_invalid <- NULL

  for (n in varstr_Y){
    X_lp_valid_b <- NULL
    Y_lp_valid_b <- NULL
    X_lp_invalid_b <- NULL
    Y_lp_invalid_b <- NULL

    Y_valid_b <- NULL
    Y_invalid_b <- NULL
    X_valid_b <- NULL
    X_invalid_b <- NULL
    for (b in 1:nbstrp){
      X_lp_valid_tmp <- NULL
      Y_lp_valid_tmp <- NULL
      X_lp_invalid_tmp <- NULL
      Y_lp_invalid_tmp <- NULL

      Y_valid_tmp <- NULL
      X_valid_tmp <- NULL
      Y_invalid_tmp <- NULL
      X_invalid_tmp <- NULL
      for (i in 1:length(preds_HT_lp$lp_X[[b]][[varstr_X]])){
        Y_maxes <- NULL
        for (m in length(varstr_Y)){
          Y_maxes[[m]] <- max(preds_HT_lp$lp_Y[[b]][[varstr_X]][[m]][i])
        }
        if (preds_HT_lp$lp_X[[b]][[varstr_X]][i]>max(unlist(Y_maxes))){
          X_lp_valid_tmp[[i]] <- preds_HT_lp$lp_X[[b]][[varstr_X]][i]
          Y_lp_valid_tmp[[i]] <- preds_HT_lp$lp_Y[[b]][[varstr_X]][[n]][i]
          Y_valid_tmp[[i]] <- preds_HT[[b]][[varstr_X]][[n]][i]
          X_valid_tmp[[i]] <- preds_margs[[b]][[varstr_X]]$maxval[i]
        } else{
          X_lp_invalid_tmp[[i]] <- preds_HT_lp$lp_X[[b]][[varstr_X]][i]
          Y_lp_invalid_tmp[[i]] <- preds_HT_lp$lp_Y[[b]][[varstr_X]][[n]][i]
          Y_invalid_tmp[[i]] <- HT_sims[[b]][[varstr_X]][[n]][i]
          X_invalid_tmp[[i]] <- preds_margs[[b]][[varstr_X]]$maxval[i]
        }
      }
      X_lp_valid_b[[b]] <- X_lp_valid_tmp
      Y_lp_valid_b[[b]] <- Y_lp_valid_tmp
      X_lp_invalid_b[[b]] <- X_lp_invalid_tmp
      Y_lp_invalid_b[[b]] <- Y_lp_invalid_tmp

      Y_valid_b[[b]] <- Y_valid_tmp
      Y_invalid_b[[b]] <- Y_invalid_tmp
      X_valid_b[[b]] <- X_valid_tmp
      X_invalid_b[[b]] <- X_invalid_tmp
    }
    X_lp_valid[[n]] <- X_lp_valid_b
    Y_lp_valid[[n]] <- Y_lp_valid_b
    X_lp_invalid[[n]] <- X_lp_invalid_b
    Y_lp_invalid[[n]] <- Y_lp_invalid_b

    Y_valid[[n]] <- Y_valid_b
    Y_invalid[[n]] <- Y_invalid_b
    X_valid[[n]] <- X_valid_b
    X_invalid[[n]] <- X_invalid_b
  }
  return(list('X_valid'=X_valid, 'Y_valid'=Y_valid, 'X_lp_valid'=X_lp_valid, 'Y_lp_valid'=Y_lp_valid))
}
