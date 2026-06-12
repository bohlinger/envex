# imports
library(ppgam)
library(evgam)
library(extRemes)
library(dplyr)
library(rlang)
library(MASS)

# adopted Ben's version to data_sub_orig[[varstr]]
# rounding function
round_to_nearest <- function(x, interval, offset) {
  rounded <- (round(x / interval) * interval) + offset
  wrapped <- rounded %% 360
  return(wrapped)
}

generate_prediction_grid <- function(var_breaks, var_lst, nr_sim){
  # Validate inputs
  stopifnot(all(var_lst %in% names(var_breaks)))

  # For each variable, sample uniformly between the innermost valid breaks
  pred_list <- lapply(var_lst, function(var) {
    breaks <- var_breaks[[var]]
    runif(nr_sim, min = breaks[2], max = breaks[length(breaks) - 1])
  })
  names(pred_list) <- var_lst

  return(as.data.frame(pred_list))
}

produce_storm_occurrences_pois <- function(model, breaks, covar_lst, RP, nr_of_years) {

  nr_per_year <- nrow(model$data) / nr_of_years
  nr_sim      <- rpois(1, nr_per_year * RP)

  newdfin <- generate_prediction_grid(breaks, covar_lst, nr_sim)

  preds_lambda <- exp(predict(model, newdata = newdfin)$location)
  newdfin$counts <- rpois(n = length(preds_lambda), lambda = preds_lambda)

  # Unfold: repeat each row according to its count
  #df_unfolded <- newdfin[rep(seq_len(nrow(newdfin)), times = newdfin$counts), covar_lst]
  df_unfolded <- newdfin[rep(seq_len(nrow(newdfin)), times = newdfin$counts), covar_lst, drop = FALSE]
  rownames(df_unfolded) <- NULL

  return(df_unfolded)
}

produce_storm_occurrences_rejection <- function(nr_of_events, RP,
                                      nr_of_years,
                                      model_nr_of_events,
                                      covarstr_lst,
                                      dfin = dfin,
                                      condition = NULL,
                                      grid_interval = NULL) {
  # simulate storm peaks from poisson
  print(c("nr_of_events:", nr_of_events))

  # probing covariate space on defined grid to cover everything efficiently
  # and hence retrieve the best lambda_max
  covar_grid_sizes <- NULL
  covar_grid <- NULL

  for (n in covarstr_lst) {
    if (is.null(grid_interval[[n]])) {
      grid_interval <- 1
    }
    vargrid <- seq(min(model_nr_of_events$knots[[n]]),
                   max(model_nr_of_events$knots[[n]]),
                   grid_interval[[n]])
    covar_grid[[n]] <- vargrid
    covar_grid_sizes[[n]] <- length(vargrid)
  }
  nr_of_grid_cells <- prod(unlist(covar_grid_sizes))
  mv_grid <- expand.grid(covar_grid)

  newdf <- as.data.frame(array(NA, dim = c((nr_of_grid_cells), (length(covarstr_lst) + 1))))
  colnames(newdf) <- c("counts", covarstr_lst)

  for (n in covarstr_lst) {
    newdf[[n]] <- mv_grid[[n]]
  }
  pred_lambdas <- predict(model_nr_of_events, newdata = newdf, type = "response")

  # get max of predicted lambdas for rejection probability
  lambda_max <- max(pred_lambdas)
  print(c("sum(predcounts)", sum(pred_lambdas)))
  print(c("lambda_max", lambda_max))

  # establish proposal_N, number of random samples from the covariate space
  proposal_N <- rpois(1, nr_of_events * lambda_max * RP)
  print(c("proposal_N:", proposal_N))

  if (proposal_N == 0){
    print("proposal_N is 0, please consider a lower extreme threshold")
  }
  # draw random values for each covariate according to their EDF
  # later replaced by rejection sampling to be more general
  dfobs <- dfin[, covarstr_lst]

  samples <- NULL
  if (length(covarstr_lst) == 1) {
    samples[[covarstr_lst[1]]] = sample(dfobs, proposal_N, replace = TRUE)
  } else if (length(covarstr_lst) > 1) {
    for (c in covarstr_lst){
      samples[[c]] = sample(dfobs[[c]], proposal_N, replace = TRUE)
    }
  }

  pred_data <- as.data.frame(samples)
  new_pred_lambda <- predict(model_nr_of_events, newdata = pred_data, type = "response")

  counts_lst <- NULL
  keep_lst <- NULL
  if (proposal_N > 0) {
    for (i in 1:proposal_N) {
      keep <- new_pred_lambda[i] / lambda_max > runif(1)
      if (keep == TRUE) {
        counts_lst[[i]] <- 1
        keep_lst[[i]] <- TRUE
      } else {
        keep_lst[[i]] <- FALSE
      }
    }

    counts_lst <- unlist(counts_lst)
    keep_unlist <- unlist(keep_lst)

    if (length(covarstr_lst) == 1) {
      dfcounts <- NULL
      dfcounts[[covarstr_lst[1]]] <- pred_data[keep_unlist, ]
      dfcounts <- as.data.frame(dfcounts)
    } else {
      dfcounts <- pred_data[keep_unlist, ]
      dfcounts$counts <- counts_lst
    }

    # filter according to covariate constraints
    if (!is.null(condition)) {
      print("filtering for condition:")
      print(c("  ", condition))
      condition_expr <- parse_expr(condition)
      dfcounts <- dfcounts %>% filter(!!condition_expr)
    }

    print(c("nr_of_sims:", dim(dfcounts)[1]))

  } else {
    dfcounts <- numeric()
    dim(dfcounts) <- c(0, 0)
  }
  return(dfcounts)
}

get_ann_max <- function(df_all, list_of_years, var_str, year_str) {
  #' @param df_all dataframe to use with all info included
  #' @param list_of_years vector of all years like seq(1980,2021)
  #' @param var_str string of variable of interest
  #' @param year_str string of variable for years
  #'
  #' @return ann_max vector of maxima per year
  #'
  #' @export
  ann_max <- array(0, length(seq_along(list_of_years))) * NA
  for (y in list_of_years) {
    df_sub <- subset(df_all, df_all[[year_str]] == y)[[var_str]]
    ann_max[y - (list_of_years[1] - 1)] <- max(df_sub)
  }
  return(ann_max)
}

calc_Hs_RV <- function(df_all, marginal_Hs_model, RP, Hs_str, year_str) {
  #' @param df_all dataframe to use with all info included
  #' @param marginal_Hs_model string of method e.g. "block"
  #' @param RP Return period e.g. 100
  #' @param Hs_str string of Hs
  #' @param year_str string of variable for years
  #'
  #' @return Hs_RV Return value of Hs for given RP
  #'
  #' @export
  if (marginal_Hs_model == "block") {
    # get annual max
    list_of_years <- seq(min(df_all[[year_str]]), max(df_all[[year_str]]))
    Hs_block_max <- get_ann_max(df_all, list_of_years, Hs_str, year_str)
    # fit GEV
    model_stat <- fevd(Hs_block_max)  # from extRemes package
    Hs_RV <- return.level(model_stat, return.period = c(RP),
                          do.ci = FALSE, alpha = 0.02, make.plot = FALSE)
  }
  return(Hs_RV)
}

get_Hs_bins_and_fit_ln <- function(df_all, Hs_str, T_str, sbinval = NULL,
                                   ebinval = NULL, binsize = 1) {
  #' @param df_all dataframe to use with all info included
  #' @param Hs_str string of Hs variable
  #' @param Tp_str string of T variable
  #'
  #' @return vectors of ln_means and ln_stds of T for Hs
  #'
  #' @export
  #'
  # bin hs in bins of 1m and perform a log-norm fit to all periods (e.g. Tp, Tm01, Tm02, ...) in one bin

  if (is.null(sbinval)) {
    sbinval <- 1
  }
  if (is.null(ebinval)) {
    ebinval <- as.integer(max(df_all[[Hs_str]]))
  }
  seqtmp <- seq(sbinval, ebinval, binsize)
  ln_means <- array(length(seqtmp)) * NA
  ln_stds <- array(length(seqtmp)) * NA
  for (i in seq_along(seqtmp)) {
    print(seqtmp[i])
    df_sub <- subset(df_all, df_all[[Hs_str]] >= (seqtmp[i] - binsize) & df_all[[Hs_str]] < (seqtmp[i]))
    fit_wb <- fitdistr(df_sub[[T_str]], "lognormal")
    ln_means[i] <- fit_wb$estimate[1]
    ln_stds[i] <- fit_wb$estimate[2]
  }
  return(list(ln_means, ln_stds))
}

fitMe <- function(params, X, Y) {
  a0 <- params[1]
  a1 <- params[2]
  a2 <- params[3]

  error <- sum((Y - (a0 + a1 * X**a2))^2)
  return(error)
}

calc_Tp_RV <- function(ln_means, Hs_RV,
                       sbinval = NULL, ebinval = NULL, binsize = 1) {
  #' @param ln_means log-normal fits of T variable for Hs bins
  #' @param Hs_RV Return value of Hs
  #'
  #' @return T_joint_RV Return value of T variable
  #'
  #' @export

  seqtmp <- seq(sbinval, ebinval, binsize)
  X <- seqtmp - (binsize / 2)
  Y <- ln_means

  o <- optim(c(a0 = 1, a1 = 1, a2 = 1), X = X, Y = Y, fitMe)
  # o$par are the optimal parameters

  T_joint_RV <- exp(o$par[1] + o$par[2] * Hs_RV**o$par[3])
  return(T_joint_RV)
}

fit_RV_from_joint_distr_ref_DNVRP <- function(df_all, marginal_Hs_model, RP,
                                              Hs_str, T_str, year_str) {
  #' @param df_all dataframe to use with all info included
  #' @param marginal_Hs_model string of method e.g. "block"
  #' @param RP Return period e.g. 100
  #' @param Hs_str string of Hs variable
  #' @param T_str string of T variable
  #' @param year_str string of variable for years
  #'
  #' @return list list of RP return value of Hs and joint T variable
  #'
  #' @export

  # 1) Fit marginal Hs and retrieve RL
  Hs_RV <- calc_Hs_RV(df_all, marginal_Hs_model, RP, Hs_str, year_str)
  # 2) Bin hs in bins of e.g. 1m and perform a log-norm fit to all periods (e.g. Tp, Tm02, ...) in one bin
  #    Caution: last bin may be underpopulated and should be treated with caution (possibly ignored) for a subsequent fit
  ln_means <- get_Hs_bins_and_fit_ln(df_all, Hs_str, T_str)[[1]]
  # 3) Optimize function for mu from DNV-BP and retrieve T_joint_RV
  T_joint_RV <- calc_Tp_RV(ln_means, Hs_RV)
  #
  return(list(Hs_RV, T_joint_RV))
}

fit_marginal_models_thr <- function(dfin, thr, model_fml,
                                    list_var = NULL, m_params = NULL,
                                    knots = NULL) {
  #' @export
  #'
  if (is.null(list_var)) {
    list_var <- names(model_fml)
  }
  margs <- NULL
  for (n in list_var){
    print(c("fit threshold model for", n))
    # bootstrap extremal threshold within bounds
    if (length(thr[[n]]) > 1) {
      t <- round(runif(1, min = thr[[n]][1], max = thr[[n]][2]), 2)
    } else {
      t <- thr[[n]]
      }

    if (is.null(m_params)) {
      margs[[n]] <- evgam(model_fml[[n]], dfin, family = "ald",
                          ald.args = list(tau = t), knots = knots)
    } else {
      margs[[n]] <- evgam(model_fml[[n]], dfin, family = "ald",
                          ald.args = list(tau = t), sp = m_params[[n]]$sp,
                          knots = knots)
    }
  }
  return (margs)
}


fit_marginal_models_gpd <- function(dfin, model_fml, list_var = NULL,
                                    family = "gpd2", m_params = NULL,
                                    knots = NULL, trace = 0) {
  #' @export
  #'

  margs <- NULL
  if (is.null(list_var)) {
    list_var <- names(model_fml)
  }
  for (n in list_var) {
    print(c("fit gpd model for", n))
    if (is.null(m_params)) {
      margs[[n]] <- evgam(model_fml[[n]], dfin[[n]], family = family,
                          trace = trace, knots = knots)
    } else {
      margs[[n]] <- evgam(model_fml[[n]], dfin[[n]], family = family,
                          sp = m_params[[n]]$sp, trace = trace,
                          knots = knots)
    }
  }
  return (margs)
}

fit_marginal_models_occ <- function(dfin, model_fml, nr_of_years,
                                    list_var = NULL, nquad = 48, m_params = NULL,
                                    knots = NULL, node_str_lst = NULL,
                                    nodes = NULL, interval = NULL) {
  #' @export
  #'
  margs <- NULL
  if (is.null(list_var)) {
    list_var <- names(model_fml)
  }

  # define weights (wts) and nodes given knots (mids and breaks)
  # if wts are not given equal weighting is assumed
  if (is.null(nodes)) {
    nodes <- NULL
    for (n in names(knots)){
      if (is.null(interval[[n]])) {
        interval[[n]] <- 1
      }
      if (n %in% node_str_lst) {
        print("computing weights and creating nodes")
        mids <- seq(min(knots[[n]]), max(knots[[n]]), interval[[n]])
        llims <- mids - interval[[n]] / 2
        breaks <- c(llims, llims[length(llims)] + interval[[n]])
        tmphist <- hist(dfin[[list_var[1]]][[n]], breaks = breaks, plot = FALSE)
        wts <- tmphist$counts / sum(tmphist$counts)
        nodes[[n]] = cbind(mids, wts)
        # occurrences only depend on prime-variable and same for all variables
        # which is why I just use [[list_var[1]]
      }
    }
  }

  for (n in list_var) {
    print(c("fit occurrence model for", n))

    if (is.null(m_params)) {
      margs[[n]] <- ppgam(model_fml[[n]],
                          data = dfin[[n]],
                          weights = nr_of_years * dim(dfin[[1]])[1],
                          nquad = nquad,
                          knots = knots[[n]],
                          nodes = nodes[[n]])
    } else {
      margs[[n]] <- ppgam(model_fml[[n]],
                          data = dfin[[n]],
                          weights = nr_of_years * dim(dfin[[1]])[1],
                          sp = m_params[[n]]$sp,
                          nquad = nquad,
                          knots = knots[[n]],
                          nodes = nodes[[n]])
    }
    margs[[n]][["nodes"]] <- nodes
    margs[[n]][["knots"]] <- knots
  }
  return (margs)
}

prepare_evgam_pois_grid <- function(dfin, var_breaks, variables) {

    # Validate inputs
    stopifnot(all(variables %in% names(dfin)))
    stopifnot(all(variables %in% names(var_breaks)))

    # Bin each variable
    for (var in variables) {
      breaks <- var_breaks[[var]]
      bin_col <- paste0(var, "_bin")
      dfin[[bin_col]] <- cut(dfin[[var]], breaks = breaks, include.lowest = TRUE)
    }

    # Build table of counts over all bin columns
    bin_cols <- paste0(variables, "_bin")
    grid_counts <- table(dfin[, bin_cols, drop = FALSE])

    # Convert to dataframe
    grid_df <- as.data.frame(grid_counts)
    names(grid_df) <- c(bin_cols, "counts")

    # Compute midpoints for each variable
    mids_list <- lapply(variables, function(var) {
      breaks <- var_breaks[[var]]
      (head(breaks, -1) + tail(breaks, -1)) / 2
    })
    names(mids_list) <- variables

    # Replace bin factor columns with numeric midpoints
    for (var in variables) {
      bin_col <- paste0(var, "_bin")
      grid_df[[var]] <- mids_list[[var]][as.integer(grid_df[[bin_col]])]
      grid_df[[bin_col]] <- NULL
    }

    # Reorder: variables first, then counts
    grid_df <- grid_df[, c(variables, "counts")]

    return(grid_df)
  }

fit_marginal_models_pois <- function(data_sub, model_fml,
                                     list_var = NULL, m_params = NULL,
                                     knots = NULL, list_covar = NULL,
                                     breaks = NULL) {
  #' @export
  #'
  margs <- NULL
  if (is.null(list_var)) {
    list_var <- names(model_fml)
  }

  for (n in list_var) {
    print(c("fit occurrence model for", n))

    # prepare gridded dataset
    grid_in <- prepare_evgam_pois_grid(
      dfin      = data_sub[[n]],
      var_breaks = breaks,
      variables  = list_covar
    )

    if (is.null(m_params)) {
      margs[[n]] <- evgam(model_fml[[n]], grid_in, family = 'poisson',
                          knots = knots)
    } else {
      margs[[n]] <- evgam(model_fml[[n]], grid_in, family = 'poisson',
                          knots = knots, sp = m_params[[n]]$sp)
    }
    margs[[n]][["knots"]] <- knots
  }
  return (margs)
}

subset_df <- function(margs, thr_str = "thr", exc_str = "exc"){
  #' @export
  #'
  list_var <- names(margs)
  data_sub <- NULL
  for (n in list_var){
    print(c("subset dataset for", n))
    dfin <- margs[[n]]$data
    dfin[[thr_str]] <- fitted(margs[[n]])$location
    dfin[[exc_str]] <- dfin[[n]] - dfin[[thr_str]]
    data_sub[[n]] <- subset(dfin, dfin[[exc_str]] > 0)
  }
  return (data_sub)
}


peak_picking <- function(dfin, thr_model_fml,
                         thr = .5, var_str = "hs", thr_str = "thr",
                         exc_str = "exc", time_str = "time",
                         decorrelation_time_scale = 2, idx = FALSE,
                         knots = NULL){
  #' @export


  if (is.null(thr_model_fml)) {
    print("apply constant threshold to data")
    dfin[[thr_str]] <- array(1, length(dfin[[var_str]])) * thr
    dfin[[exc_str]] <- dfin[[var_str]] - dfin[[thr_str]]
  } else {
    print("apply threshold model to data")
    m_ald <- evgam(thr_model_fml, dfin, family = "ald",
                   ald.args = list(tau = thr), knots = knots)

    dfin[[thr_str]] <- fitted(m_ald)$location
    dfin[[exc_str]] <- dfin[[var_str]] - dfin[[thr_str]]
  }

  print("label storms")
  dfin_labeled <- label_storms_variable_thr(dfin, exc_str = exc_str)

  print("combine and relabel storm peaks that are too close")
  # 1. round of peak finding
  df_pots_orig <- find_storm_peaks(dfin_labeled, var_str = var_str,
                                   time_str = time_str, exc_var = exc_str)

  # 2. check if storm peaks too close, if true combine
  if (isTRUE(idx)) {
    storm_idx_list <- get_list_of_close_storms_idx(df_pots_orig, var_str, time_str,
                                                   decorrelation_time_scale)
  } else{
    storm_idx_list <- get_list_of_close_storms(df_pots_orig, var_str, time_str,
                                               decorrelation_time_scale)
  }

  dfin_labeled_storms <- subset(dfin_labeled, dfin_labeled$storm_idx > 0)
  dfin_labeled_combined <- combine_storms(dfin_labeled_storms, storm_idx_list)

  # 3. find peaks again for new combined storms
  print("find peaks for combined storms")
  df_pots_orig_combined <- find_storm_peaks(dfin_labeled_combined,
                                            var_str = var_str,
                                            time_str = time_str)

  df_pots_orig_combined[[thr_str]] <- NULL
  df_pots_orig_combined[[exc_str]] <- NULL

  return (list("pots" = df_pots_orig_combined,
               "storms" = dfin_labeled_combined,
               "marg" = m_ald))
}
