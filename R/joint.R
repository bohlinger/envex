# Heffernan and Tawn conditional extremes model from 2004 as used in e.g.
# Ewans and Jonathan (2013) and Jonathan, Ewans, Randell (2013).

library(extRemes)
library(extraDistr)

transform_to_laplace <- function(dfin, margs, names_in = NULL,
                                 maxd = NULL, covarstr = NULL) {
  #' @export

  lp_margins <- NULL
  if (is.null(names_in)) {
    names_in <- names(margs)
  }
  for (n in names_in){
    natural_scale <- dfin[[n]]
    Fn_ecdf <- ecdf(natural_scale)
    probs_ecdf <- Fn_ecdf(natural_scale[natural_scale <= margs[[n]]$threshold])
    probs_gpd <- pevd(natural_scale[natural_scale > margs[[n]]$threshold],
                      scale = exp(margs[[n]]$results$par[1]),
                      shape = margs[[n]]$results$par[2],
                      threshold = margs[[n]]$threshold,
                      type = c("GP"),
                      lower.tail = TRUE)
    # combine probs
    probs <- c(probs_ecdf, probs_gpd)

    # transform to Laplace
    if (is.null(maxd)) {
      lp_margins[[n]] <- qlaplace(probs)
    } else {
      if (n == covarstr) {
        lp_margins[[n]] <- qlaplace(probs)
      } else {
        lp_margins[[n]] <- qlaplace(probs,
                                    maxd$dependence$coefficients["m"],
                                    maxd$dependence$coefficients["s"])
      }
    }
  }
  return(lp_margins)
}


## custom function for random sampling of element
runif_func <- function(n, min = 1, max = 10) sample(min:max, n, replace = TRUE)

compute_Z_from_residuals <- function(maxd, covarstr) {
  #' @export

  # create Z matrix dataframe with n x (p-1)
  models <- names(maxd[[covarstr]])
  nr_models <- length(models)
  X <- maxd[[covarstr]][[1]]$dependence$data$X
  nr_i <- length(X)

  # take X from an arbitrary object as it needs to be the same for all
  Z_M <- array(0, c(nr_i, nr_models)) * NA
  for (m in 1:nr_models) {
    if (models[m] != "Z") {
      a <- maxd[[covarstr]][[models[m]]]$dependence$coefficients[["a"]]
      b <- maxd[[covarstr]][[models[m]]]$dependence$coefficients[["b"]]
      Y <- maxd[[covarstr]][[models[m]]]$dependence$data$Y

      Z_M[, m] <- (Y - a * X) / (X^b)
    }
  }
  # name columns of Z matrix
  Zdf <- data.frame(Z_M)
  colnames(Zdf) <- c(names(maxd[[covarstr]]))
  return(Zdf)
}

HT2004_mse <- function(params, x, y) {
  alpha <- params[1]
  beta <- params[2]
  mu <- params[3]
  sigma <- params[4]

  if (sigma <= 0) return(1e8)
  s <- x^beta
  m <- alpha * x + mu * s
  s <- sigma * s
  res <- .5 * ((y - m) / s)^2
  res <- res + log(s)
  res <- sum(res)

  return(res)
}

transform_to_natural <- function(sims, maxd, targetstr) {
  #' @export

  p_lp <- plaplace(sims$data$simulated$Y,
                   maxd$dependence$coefficients["m"],
                   maxd$dependence$coefficients["s"])

  p_lp <- p_lp[p_lp < 1]
  Y_natural <- qevd(p_lp,
                    scale = exp(maxd$margs[[targetstr]]$results$par[1]),
                    shape = maxd$margs[[targetstr]]$results$par[2],
                    threshold = maxd$margs[[targetstr]]$threshold,
                    type = c("GP"))
  return(Y_natural)
}

simulate_from_HT2004 <- function(maxd, pqu, nsim) {
  #' @export
  #'
  preds <- predict(maxd, which = "hs", pqu = pqu, nsim = nsim, trace = 10)
  return(preds)
}

compute_lp_margs <- function(probs) {
  #' @export
  #'
  lp_margs <- NULL
  for (n in names(probs)) {
    lp_margs[[n]] <- qlaplace(probs[[n]])
  }
  return(lp_margs)
}

compute_probs <- function(margs_thr, margs_gpd, dfin = NULL, list_var = NULL,
                          thr_str = "thr", thr_ecdf_gpd_transition_margin = 0.01) {
  #' @export
  #'
  probs_ecdf <- NULL
  probs_gpd <- NULL
  probs <- NULL
  if (is.null(list_var)) {
    list_var <- names(margs_gpd)
  }
  for (n in list_var) {
    if (is.null(dfin)) {
      dfin <- margs_gpd[[n]]$data[[n]]
    }

    # ecdf
    Fn_ecdf <- ecdf(dfin[[n]])
    probs_ecdf[[n]] <- Fn_ecdf(dfin[[n]])

    # threshold
    threshold <- predict(margs_thr[[n]], newdata = dfin, type = "response")$location

    # gpd
    gpd_params <- predict(margs_gpd[[n]], newdata = dfin, type = "response")
    scales <- as.vector(gpd_params$scale)
    shapes <- as.vector(gpd_params$shape)

    tmp <- array(0, c(length(shapes))) * NA
    for (i in seq_along(tmp)) {
      #tmp[i] <- margs_thr[[n]]$tau + (1 - margs_thr[[n]]$tau) *
      #  pgpd(dfin[[n]][i], mu = threshold[i], sigma = scales[i],
      #       xi = shapes[i], lower.tail = TRUE)
      tmp[i] <- margs_thr[[n]]$likdata$args$tau + (1 - margs_thr[[n]]$likdata$args$tau) *
        pgpd(dfin[[n]][i], mu = threshold[i], sigma = scales[i],
             xi = shapes[i], lower.tail = TRUE)
    }
    probs_gpd[[n]] <- tmp
    rm(tmp)

    # combine probs
    tmp <- probs_gpd[[n]]
    #tmp[tmp < (margs_thr[[n]]$tau + thr_ecdf_gpd_transition_margin)] <- probs_ecdf[[n]][tmp < (margs_thr[[n]]$tau + thr_ecdf_gpd_transition_margin)]
    tmp[tmp < (margs_thr[[n]]$likdata$args$tau + thr_ecdf_gpd_transition_margin)] <- probs_ecdf[[n]][tmp < (margs_thr[[n]]$likdata$args$tau + thr_ecdf_gpd_transition_margin)]
    probs[[n]] <- unlist(tmp)
    rm(tmp)
  }
  return(list("probs" = probs, "probs_ecdf" = probs_ecdf, "probs_gpd" = probs_gpd))
}

fit_HT2004 <- function(lp_margs, thr, var_lst = NULL) {
  #' @export
  #'

  if (is.null(var_lst)){
    var_lst <- names(lp_margs)
  }

  maxd <- NULL
  for (i in seq_along(var_lst)) {
    list_of_vars <- var_lst[-i]

    X_fit_str <- var_lst[i]
    X_fit_tmp <- lp_margs[[X_fit_str]]
    for (n in list_of_vars){
      Y_fit_tmp <- lp_margs[[n]]
      Y_fit <- Y_fit_tmp[X_fit_tmp > quantile(X_fit_tmp, thr)]
      X_fit <- X_fit_tmp[X_fit_tmp > quantile(X_fit_tmp, thr)]
      # upper limit on a and b above 1 to allow for a chance to hit 1 assuming noisy data
      o <- optim(c(a = .5, b = .5, m = 1, s = 1), x = X_fit, y = Y_fit,
                 HT2004_mse, method = "L-BFGS-B",
                 lower = c(0, -100, -100, 0.001),
                 upper = c(1.2, 1.2, 100, 100))
      maxd[[X_fit_str]][[n]][["params"]]  <- o
      maxd[[X_fit_str]][[n]][["Z"]]       <- (Y_fit - o$par[1] * X_fit) / X_fit**o$par[2]
      maxd[[X_fit_str]][[n]][["X_fit"]]   <- X_fit
      maxd[[X_fit_str]][[n]][["Y_fit"]]   <- Y_fit
      maxd[[X_fit_str]][[n]][["X_all"]]   <- X_fit_tmp
      maxd[[X_fit_str]][[n]][["Y_all"]]   <- Y_fit_tmp
      maxd[[X_fit_str]][[n]][["thr"]]     <- thr
    }
  }

  return(maxd)
}

unfold_maxd_params_bstrp <- function(maxds_vars, X_var_str, Y_var_str,
                                     ulim = .01, llim = .99, cntr = .5) {
  #' @export
  #'
  alphas_cntr <- array(0, c(length(maxds_vars[[X_var_str]]))) * NA
  alphas_llim <- array(0, c(length(maxds_vars[[X_var_str]]))) * NA
  alphas_ulim <- array(0, c(length(maxds_vars[[X_var_str]]))) * NA

  betas_cntr <- array(0, c(length(maxds_vars[[X_var_str]]))) * NA
  betas_llim <- array(0, c(length(maxds_vars[[X_var_str]]))) * NA
  betas_ulim <- array(0, c(length(maxds_vars[[X_var_str]]))) * NA

  for (t in seq_along(maxds_vars[[X_var_str]])) {
    alphas <- array(0, c(length(maxds_vars[[X_var_str]][[1]]))) * NA
    betas <- array(0, c(length(maxds_vars[[X_var_str]][[1]]))) * NA
    for (i in seq_along(alphas)) {
      alphas[i] <- as.numeric(maxds_vars[[X_var_str]][[t]][[i]][[X_var_str]][[Y_var_str]]$params$par[1])
    }
    for (i in seq_along(betas)) {
      betas[i] <- as.numeric(maxds_vars[[X_var_str]][[t]][[i]][[X_var_str]][[Y_var_str]]$params$par[2])
    }

    alphas_cntr[t] <- quantile(alphas, cntr)
    alphas_llim[t] <- quantile(alphas, llim)
    alphas_ulim[t] <- quantile(alphas, ulim)

    betas_cntr[t] <- quantile(betas, cntr)
    betas_llim[t] <- quantile(betas, llim)
    betas_ulim[t] <- quantile(betas, ulim)
  }
  return(list("alphas_llim" = alphas_llim, "alphas_cntr" = alphas_cntr,
              "alphas_ulim" = alphas_ulim, "betas_llim" = betas_llim,
              "betas_cntr" = betas_cntr, "betas_ulim" = betas_ulim))
}

unfold_maxd_params_bstrp_v2 <- function(maxds, X_var_str, Y_var_str,
                                        llim = .01, ulim = .99, cntr = .5) {
  #' @export
  #' models_maxds_thr[[t_maxd]][[b]]$hs$U10
  #'
  alphas_cntr <- array(0, c(length(maxds))) * NA
  alphas_llim <- array(0, c(length(maxds))) * NA
  alphas_ulim <- array(0, c(length(maxds))) * NA

  betas_cntr <- array(0, c(length(maxds))) * NA
  betas_llim <- array(0, c(length(maxds))) * NA
  betas_ulim <- array(0, c(length(maxds))) * NA

  mus_cntr <- array(0, c(length(maxds))) * NA
  mus_llim <- array(0, c(length(maxds))) * NA
  mus_ulim <- array(0, c(length(maxds))) * NA

  sigmas_cntr <- array(0, c(length(maxds))) * NA
  sigmas_llim <- array(0, c(length(maxds))) * NA
  sigmas_ulim <- array(0, c(length(maxds))) * NA

  for (t in seq_along(maxds)) {
    alphas <- array(0, c(length(maxds[[t]]))) * NA
    betas <- array(0, c(length(maxds[[t]]))) * NA
    mus <- array(0, c(length(maxds[[t]]))) * NA
    sigmas <- array(0, c(length(maxds[[t]]))) * NA
    for (b in seq_along(alphas)) {
      alphas[b] <- as.numeric(maxds[[t]][[b]][[X_var_str]][[Y_var_str]]$params$par[1])
      betas[b] <- as.numeric(maxds[[t]][[b]][[X_var_str]][[Y_var_str]]$params$par[2])
      mus[b] <- as.numeric(maxds[[t]][[b]][[X_var_str]][[Y_var_str]]$params$par[3])
      sigmas[b] <- as.numeric(maxds[[t]][[b]][[X_var_str]][[Y_var_str]]$params$par[4])
    }

    alphas_cntr[t] <- quantile(alphas, cntr, na.rm = TRUE)
    alphas_llim[t] <- quantile(alphas, llim, na.rm = TRUE)
    alphas_ulim[t] <- quantile(alphas, ulim, na.rm = TRUE)

    betas_cntr[t] <- quantile(betas, cntr, na.rm = TRUE)
    betas_llim[t] <- quantile(betas, llim, na.rm = TRUE)
    betas_ulim[t] <- quantile(betas, ulim, na.rm = TRUE)

    mus_cntr[t] <- quantile(mus, cntr, na.rm = TRUE)
    mus_llim[t] <- quantile(mus, llim, na.rm = TRUE)
    mus_ulim[t] <- quantile(mus, ulim, na.rm = TRUE)

    sigmas_cntr[t] <- quantile(sigmas, cntr, na.rm = TRUE)
    sigmas_llim[t] <- quantile(sigmas, llim, na.rm = TRUE)
    sigmas_ulim[t] <- quantile(sigmas, ulim, na.rm = TRUE)
  }

  return(list("alpha_llim" = alphas_llim, "alpha_cntr" = alphas_cntr, "alpha_ulim" = alphas_ulim,
              "beta_llim" = betas_llim, "beta_cntr" = betas_cntr, "beta_ulim" = betas_ulim,
              "mu_llim" = mus_llim, "mu_cntr" = mus_cntr, "mu_ulim" = mus_ulim,
              "sigma_llim" = sigmas_llim, "sigma_cntr" = sigmas_cntr, "sigma_ulim" = sigmas_ulim))
}

fit_margs_bstrp <- function(dfin,
                            model_fml_thr, model_fml_occ, model_fml_gpd,
                            extr_thr, nr_of_years,
                            thr_str = "thr",
                            list_var = NULL,
                            list_covar = NULL,
                            margs_thr_orig = NULL,
                            margs_gpd_orig = NULL,
                            margs_occ_orig = NULL,
                            nbstrp = NULL,
                            nquad = 40,
                            knots = NULL,
                            mids = NULL,
                            breaks = NULL,
                            interval = NULL,
                            nodes = NULL,
                            node_str_lst = NULL,
                            family_gpd = 'gpd2',
                            family_occ = 'pois',
                            grid_in = NULL) {
  #' @export
  #'

  if (is.null(list_var)) {
    list_var <- names(model_fml_thr)
  }
  print(c("Considered variables:", list_var))

  if (is.null(nbstrp)) {
    nbstrp <- length(dfin)
  }
  print(c("number of bootstraps is:", nbstrp))

  margs <- NULL
  for (i in 1:nbstrp) {

    ### Fitting marginals ###

    print(c("bootstrap nr:", i))
    print("fit threshold model")

    margs_thr <- fit_marginal_models_thr(dfin = dfin[[i]],
                                         list_var = list_var,
                                         thr = extr_thr,
                                         model_fml = model_fml_thr,
                                         m_params = margs_thr_orig,
                                         knots = knots)

    print("subset to pots above threshold")
    data_sub <- subset_df(margs_thr, thr_str = "thr", exc_str = "exc")

    print("fit GPD model")
    margs_gpd <- fit_marginal_models_gpd(dfin = data_sub,
                                         list_var = list_var,
                                         model_fml = model_fml_gpd,
                                         m_params = margs_gpd_orig,
                                         knots = knots,
                                         family = family_gpd)

    print("fit occ model")
    if (family_occ == 'pois'){
      margs_occ <- fit_marginal_models_pois(data_sub, model_fml_occ,
                                            list_var = list_var,
                                            list_covar = list_covar,
                                            knots = knots,
                                            breaks = breaks,
                                            m_params = margs_occ_orig)
    } else if(family_occ == 'pp'){
      margs_occ <- fit_marginal_models_occ(dfin = data_sub,
                                           list_var = list_var,
                                           model_fml = model_fml_occ,
                                           nr_of_years = nr_of_years,
                                           nquad = nquad,
                                           knots = knots,
                                           node_str_lst = node_str_lst,
                                           nodes = nodes,
                                           interval = interval)
    }

    probs <- compute_probs(margs_thr, margs_gpd, dfin = dfin[[i]],
                           list_var = list_var, thr_str = "thr")[["probs"]]

    margs_tmp <- NULL
    margs_tmp[["thr"]] <- margs_thr
    margs_tmp[["gpd"]] <- margs_gpd
    margs_tmp[["occ"]] <- margs_occ
    margs_tmp[["probs"]] <- probs

    margs[[i]] <- margs_tmp
  }
  return(margs)
}

fit_maxds_bstrp <- function(dfin, margs, maxd_thr, nbstrp = NULL, var_lst = NULL) {
  #' @export
  #'

  if (is.null(nbstrp)) {
    nbstrp <- length(dfin)
  }

  if (length(maxd_thr) > 1) {
    t <- runif(1, min = maxd_thr[1], max = maxd_thr[2])
  } else {
    t <- maxd_thr
  }

  maxds <- NULL
  for (b in 1:nbstrp) {
    # print(c("bootstrap nr:",b))
    lp_margs <- compute_lp_margs(margs[[b]]$probs)
    maxd <- fit_HT2004(lp_margs, thr = t, var_lst = var_lst)
    maxds[[b]] <- maxd
  }

  return(maxds)
}

predict_maxd <- function(maxd, varstr_X, varstr_Y, nsim) {
  tmp_alpha <- maxd[[varstr_X]][[varstr_Y]]$params$par[1]
  tmp_beta <- maxd[[varstr_X]][[varstr_Y]]$params$par[2]
  tmp_Z <- maxd[[varstr_X]][[varstr_Y]][["Z"]]
  lp_X <- maxd[[varstr_X]][[varstr_Y]]$X_fit

  Xsim <- array(0, nsim) * NA
  Ysim <- array(0, nsim) * NA
  for (i in 1:nsim) {
    Xsim[i] <- lp_X[runif_func(1, min = 1, max = length(lp_X))]
    Ysim[i] <- tmp_alpha * Xsim[i] + Xsim[i]**(tmp_beta) * tmp_Z[runif_func(1, min = 1, max = length(tmp_Z))]
  }

  return(list("Xsim" = Xsim, "Ysim" = Ysim))
}

predict_from_HT2004_models <- function(margs, maxds,
                                       var_lst, nbstrp,
                                       preds = NULL,
                                       covar_lst = NULL){
  #' @export

  if (is.null(var_lst)) {
    var_lst <- names(maxds[[1]])
  }

  if (is.null(covar_lst)) {
    covar_lst <- margs[[1]][[1]][[1]]$predictor.names
  }

  lp_Y_bstrp <- NULL
  lp_X_bstrp <- NULL
  X_bstrp <- NULL
  covar_bstrp <- NULL

  for (b in 1:nbstrp){
    lp_Y_lst <- NULL
    lp_X_lst <- NULL
    X_lst <- NULL
    covar <- NULL

    for (n in seq_along(var_lst)) {
      lp_Y_lst_tmp <- NULL
      tmplst <- var_lst[-n]
      lp_X <- qlaplace(preds[[b]][[var_lst[n]]]$prob)
      X_sim <- preds[[b]][[var_lst[n]]]$maxval

      print("predict from HT2004")
      for (m in tmplst) {
        tmp_alpha <- maxds[[b]][[var_lst[n]]][[m]]$params$par[1]
        tmp_beta <- maxds[[b]][[var_lst[n]]][[m]]$params$par[2]
        tmp_mu <- maxds[[b]][[var_lst[n]]][[m]]$params$par[3]
        tmp_sig <- maxds[[b]][[var_lst[n]]][[m]]$params$par[4]

        tmp_Z <- maxds[[b]][[var_lst[n]]][[m]][["Z"]]
        lp_Y <- tmp_alpha * lp_X + lp_X**(tmp_beta) *
                tmp_Z[runif_func(length(lp_X), min = 1, max = length(tmp_Z))]
        lp_Y_lst_tmp[[m]] <- as.numeric(lp_Y)
      }

      lp_Y_lst[[var_lst[n]]] <- lp_Y_lst_tmp
      lp_X_lst[[var_lst[n]]] <- lp_X
      X_lst[[var_lst[n]]] <- X_sim
      for (nc in seq_along(covar_lst)){
        covar[[var_lst[n]]][[covar_lst[nc]]] <- preds[[b]][[var_lst[n]]][[covar_lst[nc]]]
      }
    }

    lp_Y_bstrp[[b]] <- lp_Y_lst
    lp_X_bstrp[[b]] <- lp_X_lst
    X_bstrp[[b]] <- X_lst
    covar_bstrp[[b]] <- covar
  }

  return(list("lp_Y" = lp_Y_bstrp, "lp_X" = lp_X_bstrp, "X" = X_bstrp, "covar" = covar_bstrp))
}

convert_HT2004_preds_to_original_space <- function(maxds,
                                                   maxd_preds, marg_preds,
                                                   var_lst = NULL,
                                                   nbstrp = NULL) {
  #' @export

  if (is.null(var_lst)) {
    var_lst <- names(maxds[[1]])
  }

  if (is.null(nbstrp)) {
    nbstrp <- length(marg_preds)
  }

  Y_bstrp <- NULL
  for (b in 1:nbstrp) {
    Y_lst <- NULL
    covar <- NULL
    for (n in 1:length(var_lst)) {
      Y_lst_tmp <- NULL
      tmplst <- var_lst[-n]
      for (m in tmplst){
        lp_prob <- plaplace(maxd_preds$lp_Y[[b]][[var_lst[n]]][[m]])
        thr_sim <- marg_preds[[b]][[m]]$thr
        scale_sim <- marg_preds[[b]][[m]]$scale
        shape_sim <- marg_preds[[b]][[m]]$shape
        Y_lst_tmp[[m]] <- qgpd(lp_prob, mu = thr_sim, sigma = scale_sim,
                               xi = shape_sim, lower.tail = TRUE, log.p = FALSE)
      }
      Y_lst_tmp[[var_lst[n]]] <- maxd_preds$X[[b]][[var_lst[n]]]
      Y_lst[[var_lst[n]]] <- Y_lst_tmp
      covar[[var_lst[n]]] <- maxd_preds$covar[[b]][[var_lst[n]]]
    }
    Y_bstrp[[b]] <- Y_lst
    Y_bstrp[[b]][["covar"]] <- covar
  }
  return(Y_bstrp)
}

reorganize_HT2004_preds <- function(maxds,
                                    maxd_preds, marg_preds,
                                    var_lst = NULL,
                                    nbstrp = NULL) {
  #' @export

  if (is.null(var_lst)) {
    var_lst <- names(maxds[[1]])
  }

  if (is.null(nbstrp)) {
    nbstrp <- length(marg_preds)
  }

  Y_bstrp <- NULL
  for (b in 1:nbstrp) {
    Y_lst <- NULL
    covar <- NULL
    for (n in 1:length(var_lst)) {
      Y_lst_tmp <- NULL
      tmplst <- var_lst[-n]
      for (m in tmplst){
        Y_lst_tmp[[m]] <- maxd_preds$lp_Y[[b]][[var_lst[n]]][[m]]
      }
      Y_lst_tmp[[var_lst[n]]] <- maxd_preds$lp_X[[b]][[var_lst[n]]]
      Y_lst[[var_lst[n]]] <- Y_lst_tmp
      covar[[var_lst[n]]] <- maxd_preds$covar[[b]][[var_lst[n]]]
    }
    Y_bstrp[[b]] <- Y_lst
    Y_bstrp[[b]][["covar"]] <- covar
  }
  return(Y_bstrp)
}

filter_preds_maxds <- function(preds, var_lst, preds2 = NULL) {

  result1 <- vector("list", length(preds))
  result2 <- if (!is.null(preds2)) vector("list", length(preds)) else NULL

  for (b in seq_along(preds)) {
    result1[[b]] <- setNames(vector("list", length(var_lst)), var_lst)
    if (!is.null(preds2)) result2[[b]] <- setNames(vector("list", length(var_lst)), var_lst)

    for (xi in var_lst) {

      Xi     <- preds[[b]][[xi]]
      Xi_ref <- Xi[[xi]]                        # diagonal reference

      # find indices where ALL Xj <= Xi_ref
      keep_idx <- Reduce("&", lapply(Xi, function(Xj) Xj <= Xi_ref))

      # filter covariates once, reuse for both outputs
      covar_filtered <- lapply(preds[[b]]$covar[[xi]], function(cov) cov[keep_idx])

      # 1) filter preds and attach covariates
      result1[[b]][[xi]] <- lapply(Xi, function(Xj) Xj[keep_idx])
      result1[[b]][[xi]]$covar <- covar_filtered

      # 2) filter preds2 and attach same covariates
      if (!is.null(preds2)) {
        result2[[b]][[xi]] <- lapply(preds2[[b]][[xi]], function(Xj) Xj[keep_idx])
        result2[[b]][[xi]]$covar <- covar_filtered
      }
    }
  }

  # Build output
  if (!is.null(preds2)) {
    list(preds_filtered  = result1,
         preds2_filtered = result2)
  } else {
    list(preds_filtered = result1)
  }
}

#filter_preds_maxds <- function(preds, var_lst, preds2 = NULL) {
#  #' @export
#  result1 <- vector("list", length(preds))
#  result2 <- if (!is.null(preds2)) vector("list", length(preds)) else NULL
#  for (b in seq_along(preds)) {
#    result1[[b]] <- setNames(lapply(var_lst, function(xi) {
#      Xi <- preds[[b]][[xi]]
#      Xi_ref <- Xi[[xi]]
#      # find indices where ALL Xj <= Xi_ref
#      keep_idx <- Reduce("&", lapply(Xi, function(Xj) Xj <= Xi_ref))
#      # apply filter to preds
#      filtered1 <- lapply(Xi, function(Xj) Xj[keep_idx])
#      # apply same filter to preds2 if provided
#      if (!is.null(preds2)) {
#        result2[[b]][[xi]] <<- lapply(preds2[[b]][[xi]], function(Xj) Xj[keep_idx])
#      }
#      filtered1
#    }), var_lst)
#  }
#  if (!is.null(preds2)) {
#    result2 <- lapply(result2, function(b) setNames(b, var_lst))
#    list(preds_filtered = result1, preds2_filtered = result2)
#  } else {
#    result1
#  }
#}

# Helper: extract one dataframe from a given bootstrap and xi block
extract_df <- function(preds2_filtered, b, xi, var_lst) {
  xi_data <- preds2_filtered[[b]][[xi]]
  as.data.frame(setNames(
    lapply(var_lst, function(xj) xi_data[[xj]]),
    var_lst
  ))
}

#prepare_preds_df <- function(preds2_filtered, var_lst, xi_select = NULL) {
#
#  # Default to all xi if not specified
#  if (is.null(xi_select)) xi_select <- var_lst
#
#  # Combine all bootstraps and selected xi into one dataframe
#  do.call(rbind, lapply(xi_select, function(xi)
#    do.call(rbind, lapply(seq_along(preds2_filtered), function(b)
#      extract_df(preds2_filtered, b, xi, var_lst)
#    ))
#  ))
#}

prepare_preds_df <- function(preds_filtered, var_lst, xi_select = NULL) {

  # Default to all xi if not specified
  if (is.null(xi_select)) xi_select <- var_lst

  # Combine all bootstraps and selected xi into one dataframe
  do.call(rbind, lapply(xi_select, function(xi)
    do.call(rbind, lapply(seq_along(preds_filtered), function(b) {

      xi_data <- preds_filtered[[b]][[xi]]

      # extract the main variables (excluding covar)
      df_main <- as.data.frame(setNames(
        lapply(var_lst, function(xj) xi_data[[xj]]),
        var_lst
      ))

      # extract covariates if present and cbind to main df
      if (!is.null(xi_data$covar)) {
        df_covar <- as.data.frame(xi_data$covar)
        df_main  <- cbind(df_main, df_covar)
      }

      df_main
    }))
  ))
}

predict_marg <- function(margs, nr_of_years, RP, varstr = NULL, nmc = 1,
                         covarstr_lst,
                         condition = NULL, grid_interval = NULL,
                         dfin = NULL,
                         family_occ = 'pois',
                         breaks = NULL) {
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
  nr_of_events <- length(tmp[tmp > 0])
  rm(tmp)

  if (family_occ == 'pois'){
    df_storm_cov <- produce_storm_occurrences_pois(margs$occ[[varstr]],
                                                   breaks, covarstr_lst,
                                                   RP, nr_of_years)
  } else if (family_occ == 'pp'){
    df_storm_cov <- produce_storm_occurrences_rejection(nr_of_events, RP, nr_of_years,
                                                        margs$occ[[varstr]],
                                                        covarstr_lst,
                                                        dfin = dfin,
                                                        condition = condition,
                                                        grid_interval = grid_interval)
  }

  if (dim(df_storm_cov)[1] > 0) {
    # predict from ALD
    print("predict threshold")
    thr <- predict(margs$thr[[varstr]], newdata = df_storm_cov,
                   type = "response")$location

    # predict from GPD
    print("predict exceedences")
    gpd_param_sims <- predict(margs$gpd[[varstr]],
                              newdata = df_storm_cov, type = "response")
    scales <- gpd_param_sims$scale
    shapes <- gpd_param_sims$shape

    for (i in 1:nmc){
      gpd_sims <- rgpd(length(scales), mu = 0, sigma = scales, xi = shapes)
      res_sims <- rgpd(length(scales), mu = thr, sigma = scales, xi = shapes)

      max_idx <- which(res_sims == max(res_sims))

      # Save maximum
      max_val <- res_sims[max_idx]
      max_scale <- scales[max_idx]
      max_shape <- shapes[max_idx]
      thr_max <- thr[max_idx]

      if (length(predictor_names) == 1) {
        df_pred_cov <- NULL
        df_pred_cov[[predictor_names[1]]] <- df_storm_cov[max_idx, ]
        df_pred_cov <- as.data.frame(df_pred_cov)
      } else {
        df_pred_cov <- df_storm_cov[max_idx, ]
      }

      gpd_prob <- pevd(res_sims[max_idx], scale = max_scale, shape = max_shape,
                       threshold = thr_max, type = "GP", lower.tail = TRUE)

      # store in output field
      df_max$maxval[i] <- max_val
      df_max$scale[i] <- max_scale
      df_max$shape[i] <- max_shape
      df_max$thr[i] <- thr_max
      df_max$prob[i] <- gpd_prob

      # store covariates
      for (n in predictor_names) {
        df_covs[[n]][[i]] <- df_pred_cov[[n]]
      }
    }
  }

  # combine dfs
  df_out <- cbind(df_max, df_covs)

  return(df_out)
}

predict_margs <- function(margs, nr_of_years, RP, nmc = 1,
                          var_lst = NULL, nbstrp = NULL,
                          covarstr_lst,
                          condition = NULL, grid_interval = NULL,
                          family_occ = 'pois',
                          breaks = NULL) {
  #' @export
  #'

  if (is.null(var_lst)) {
    var_lst <- names(margs[[1]]$thr)
  }

  if (is.null(nbstrp)) {
    nbstrp <- length(margs)
  }

  preds_lst <- NULL
  for (b in 1:nbstrp){ # number of bootstraps
    print(c("number of boostraps (predict_margs):", b))
    preds_tmp_lst <- NULL
    for (n in seq_along(var_lst)){
      dfin <- margs[[b]]$gpd[[var_lst[n]]]$data
      preds <- predict_marg(margs[[b]], nr_of_years, RP,
                            varstr = var_lst[n], nmc = nmc,
                            covarstr_lst,
                            condition = condition,
                            grid_interval = grid_interval,
                            dfin = dfin,
                            family_occ = family_occ,
                            breaks = breaks)
      preds_tmp_lst[[var_lst[n]]] <- preds
    }
    preds_lst[[b]] <- preds_tmp_lst
  }
  return(preds_lst)
}

retrieve_valid_HT_samples <- function(preds_HT, preds_HT_lp, preds_margs,
                                      varstr_X, varstr_Y, nbstrp = NULL,
                                      region = "default", models_maxds = NULL) {
  #' @export

  # region: relates to region in LP space where different HT2004 models are
  #         valid/appropriate. X is always assumed to be the conditioning
  #         variable and Y the conditioned variable. By default only the
  #         valid/appropriate region is output, this can be extended to
  #         including the entire rectangle where X greater than a threshold.
  #         This is the case when not filtering for valids at all, keyword "all".
  #         Additionally this can be reduced to only the triangle where one
  #         model is above threshold for both X and Y. To cover this entire
  #         square where X and Y both are extreme this script needs to be run
  #         twice swapping X and Y and consolidating the results, the keyword
  #         for this is "both".

  if (is.null(nbstrp)) {
    nbstrp <- length(preds_margs)
  }

  Y_lp_valid <- NULL
  Y_lp_invalid <- NULL
  X_lp_valid <- NULL
  X_lp_invalid <- NULL

  Y_valid <- NULL
  Y_invalid <- NULL
  X_valid <- NULL
  X_invalid <- NULL

  covar_valid <- NULL
  covar_invalid <- NULL

  if (region == "default") {
    print("Greater than threshold on X, only valids")
    for (n in varstr_Y){
      X_lp_valid_b <- NULL
      Y_lp_valid_b <- NULL
      X_lp_invalid_b <- NULL
      Y_lp_invalid_b <- NULL

      Y_valid_b <- NULL
      Y_invalid_b <- NULL
      X_valid_b <- NULL
      X_invalid_b <- NULL

      covar_valid_b <- NULL
      covar_invalid_b <- NULL

      for (b in 1:nbstrp){
        lp_Y_maxes <- NULL
        for (m in length(varstr_Y)){
          lp_Y_maxes[[m]] <- unlist(preds_HT_lp$lp_Y[[b]][[varstr_X]][[m]])
        }
        lp_Y_maxes <- unlist(lp_Y_maxes)

        covar_names <- names(preds_HT[[b]]$covar[[n]])

        X_lp_valid_tmp <- preds_HT_lp$lp_X[[b]][[varstr_X]]
        Y_lp_valid_tmp <- preds_HT_lp$lp_Y[[b]][[varstr_X]][[n]]

        accept_idx <- which(X_lp_valid_tmp > Y_lp_valid_tmp)
        reject_idx <- which(X_lp_valid_tmp < Y_lp_valid_tmp)

        X_lp_valid_tmp <- X_lp_valid_tmp[accept_idx]
        Y_lp_valid_tmp <- Y_lp_valid_tmp[accept_idx]
        X_lp_invalid_tmp <- X_lp_valid_tmp[reject_idx]
        Y_lp_invalid_tmp <- Y_lp_valid_tmp[reject_idx]

        Y_valid_tmp <- preds_HT[[b]][[varstr_X]][[n]][accept_idx]
        X_valid_tmp <- preds_margs[[b]][[varstr_X]]$maxval[accept_idx]
        Y_invalid_tmp <- preds_HT[[b]][[varstr_X]][[n]][reject_idx]
        X_invalid_tmp <- preds_margs[[b]][[varstr_X]]$maxval[reject_idx]

        covar_valid_tmp <- NULL
        covar_invalid_tmp <- NULL
        for (cn in covar_names){
          covar_valid_tmp[[cn]] <- preds_HT[[b]][["covar"]][[n]][[cn]][accept_idx]
          covar_invalid_tmp[[cn]] <- preds_HT[[b]][["covar"]][[n]][[cn]][reject_idx]
        }

        X_lp_valid_b[[b]] <- X_lp_valid_tmp
        Y_lp_valid_b[[b]] <- Y_lp_valid_tmp
        X_lp_invalid_b[[b]] <- X_lp_invalid_tmp
        Y_lp_invalid_b[[b]] <- Y_lp_invalid_tmp

        Y_valid_b[[b]] <- Y_valid_tmp
        Y_invalid_b[[b]] <- Y_invalid_tmp
        X_valid_b[[b]] <- X_valid_tmp
        X_invalid_b[[b]] <- X_invalid_tmp

        covar_valid_b[[b]] <- covar_valid_tmp
        covar_invalid_b[[b]] <- covar_invalid_tmp
      }
      X_lp_valid[[n]] <- X_lp_valid_b
      Y_lp_valid[[n]] <- Y_lp_valid_b
      X_lp_invalid[[n]] <- X_lp_invalid_b
      Y_lp_invalid[[n]] <- Y_lp_invalid_b

      Y_valid[[n]] <- Y_valid_b
      Y_invalid[[n]] <- Y_invalid_b
      X_valid[[n]] <- X_valid_b
      X_invalid[[n]] <- X_invalid_b

      covar_valid[[n]] <- covar_valid_b
      covar_invalid[[n]] <- covar_invalid_b
    }
  } else if (region == "all") {
    X_lp_valid_b <- NULL
    Y_lp_valid_b <- NULL
    X_lp_invalid_b <- NULL
    Y_lp_invalid_b <- NULL

    Y_valid_b <- NULL
    Y_invalid_b <- NULL
    X_valid_b <- NULL
    X_invalid_b <- NULL

    covar_valid_b <- NULL
    covar_invalid_b <- NULL

    for (b in 1:nbstrp) {
      lp_Y_maxes <- NULL
      for (m in length(varstr_Y)) {
        lp_Y_maxes[[m]] <- unlist(preds_HT_lp$lp_Y[[b]][[varstr_X]][[m]])
      }
      lp_Y_maxes <- unlist(lp_Y_maxes)

      covar_names <- names(preds_HT[[b]]$covar[[n]])

      X_lp_valid_tmp <- preds_HT_lp$lp_X[[b]][[varstr_X]]
      Y_lp_valid_tmp <- preds_HT_lp$lp_Y[[b]][[varstr_X]][[n]]
      X_lp_invalid_tmp <- NULL
      Y_lp_invalid_tmp <- NULL

      Y_valid_tmp <- preds_HT[[b]][[varstr_X]][[n]]
      X_valid_tmp <- preds_margs[[b]][[varstr_X]]$maxval
      Y_invalid_tmp <- NULL
      X_invalid_tmp <- NULL

      covar_valid_tmp <- NULL
      covar_invalid_tmp <- NULL
      for (cn in covar_names) {
        covar_valid_tmp[[cn]] <- preds_HT[[b]][["covar"]][[n]][[cn]]
      }

      X_lp_valid_b[[b]] <- X_lp_valid_tmp
      Y_lp_valid_b[[b]] <- Y_lp_valid_tmp
      X_lp_invalid_b[[b]] <- X_lp_invalid_tmp
      Y_lp_invalid_b[[b]] <- Y_lp_invalid_tmp

      Y_valid_b[[b]] <- Y_valid_tmp
      Y_invalid_b[[b]] <- Y_invalid_tmp
      X_valid_b[[b]] <- X_valid_tmp
      X_invalid_b[[b]] <- X_invalid_tmp

      covar_valid_b[[b]] <- covar_valid_tmp
      covar_invalid_b[[b]] <- covar_invalid_tmp
    }
    X_lp_valid[[n]] <- X_lp_valid_b
    Y_lp_valid[[n]] <- Y_lp_valid_b
    X_lp_invalid[[n]] <- X_lp_invalid_b
    Y_lp_invalid[[n]] <- Y_lp_invalid_b

    Y_valid[[n]] <- Y_valid_b
    Y_invalid[[n]] <- Y_invalid_b
    X_valid[[n]] <- X_valid_b
    X_invalid[[n]] <- X_invalid_b

    covar_valid[[n]] <- covar_valid_b
    covar_invalid[[n]] <- covar_invalid_b
  }
  return(list("X_valid" = X_valid, "Y_valid" = Y_valid,
              "X_lp_valid" = X_lp_valid, "Y_lp_valid" = Y_lp_valid,
              "covar_valid" = covar_valid))
}

retrieve_valid_HT_samples_v2 <- function(preds_HT, preds_HT_lp, preds_margs,
                                         varstr_X, varstr_Y, nbstrp = NULL,
                                         region = "default", models_maxds = NULL) {
  #' @export

  if (is.null(nbstrp)) nbstrp <- length(preds_margs)

  # --- helper: process one bootstrap replicate for one Y variable ----------
  .process_one_b <- function(b, n) {
    # seq_along, so every element of varstr_Y is used
    lp_Y_maxes <- unlist(lapply(seq_along(varstr_Y), function(m)
      unlist(preds_HT_lp$lp_Y[[b]][[varstr_X]][[m]])))

    lp_X   <- preds_HT_lp$lp_X[[b]][[varstr_X]]
    lp_Y   <- preds_HT_lp$lp_Y[[b]][[varstr_X]][[n]]

    # derive both index sets from the *original* lp_X before subsetting
    accept_idx <- which(lp_X >  lp_Y)
    reject_idx <- which(lp_X <= lp_Y)   # use <= for a clean partition

    covar_names <- names(preds_HT[[b]]$covar[[n]])

    .subset_covar <- function(idx)
      lapply(setNames(covar_names, covar_names),
             function(cn) preds_HT[[b]][["covar"]][[n]][[cn]][idx])

    list(
      X_lp_valid   = lp_X[accept_idx],
      Y_lp_valid   = lp_Y[accept_idx],
      X_lp_invalid = lp_X[reject_idx],   # index into originals
      Y_lp_invalid = lp_Y[reject_idx],
      Y_valid       = preds_HT[[b]][[varstr_X]][[n]][accept_idx],
      X_valid       = preds_margs[[b]][[varstr_X]]$maxval[accept_idx],
      Y_invalid     = preds_HT[[b]][[varstr_X]][[n]][reject_idx],
      X_invalid     = preds_margs[[b]][[varstr_X]]$maxval[reject_idx],
      covar_valid   = .subset_covar(accept_idx),
      covar_invalid = .subset_covar(reject_idx)
    )
  }

  # --- helper: process all bootstrap replicates for one Y variable ----------
  .process_one_n <- function(n) {
    bs <- lapply(seq_len(nbstrp), .process_one_b, n = n)
    # Transpose: list-of-replicates → named list-of-fields
    lapply(setNames(nm = names(bs[[1]])), function(field)
      lapply(bs, `[[`, field))
  }

  # --- main: iterate over Y variables ---------------------------------------
  message("Greater than threshold on X, only valids")
  results <- lapply(setNames(nm = varstr_Y), .process_one_n)

  # --- reshape output to match original return structure --------------------
  fields_to_return <- c("X_valid", "Y_valid", "X_lp_valid", "Y_lp_valid", "covar_valid")
  out <- lapply(setNames(nm = fields_to_return), function(field)
    lapply(results, `[[`, field))

  out
}

consolidate_HT2004_preds <- function(preds_maxds){
  var_lst <- names(preds_maxds[[1]])
  preds_cons <- NULL
  for (b in seq_len(preds_maxds)){
    for (n in var_lst){
      tmp_lst <- var_lst[-n]
      for (m in tmp_lst){
       preds_cons[[m]][[n]][[b]] <- preds_maxds[[b]][[n]][[m]]
      }
    }
  tmp_lst <- var_lst[-n]}
}
