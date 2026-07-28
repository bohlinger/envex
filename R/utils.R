create_settings_object <- function(settings) {
  #' Initialize a settings object to pass to simulate_storm_peaks function
  #'
  #' @param settings list of kwargs
  #' @return settings a settings object
  #' @examples
  #' fmla_ald <- hs ~ te(doy, Pdir, bs=c("cc","cc"), k=c(6,6))
  #' fmla_gpd_pots <- list(exc ~ te(doy, Pdir, k=c(6,6), bs=c("cc","cc")),
  #'                       ~ te(doy, Pdir, k=c(6,6), bs=c("cc","cc")))
  #' settings_list <- list(N_bstrp = 2, var_str = "hs", time_str = "time",
  #'                      fmla_gpd, fmla_ald,
  #'                      q_thr = .9, storm_thr = 2, nr_of_years = 42,
  #'                      RP = 100, min_nr_days = 2)
  #' settings <- create_settings_object(settings=settings_list)
  #' @export

  class(settings) <- "Settings"

  return(settings)
}

runif_func <- function(n, min = 0, max = 1) sample(min:max, n, replace = TRUE)

adjust_dirs_to_start_dir <- function(d, start_dir) {
  #' @export
  #'
  d_new <- (d - start_dir) %% 360
  return(d_new)
}

filter_dir_sector <- function(df, dir_str, start_dir, sector_width) {
  #'
  #' @param df A dataframe.
  #' @param dir_str string of direction variable
  #' @param start_dir direction to start sector from (clockwise)
  #' @param sector_width width in degrees

  #' @return df A filtered dataframe.
  #'
  #' @examples
  #' df_sector <- check_dir_sector(ekofisk, "Pdir", 270, 180)
  #'
  #' @export

  df$tmp_var <- adjust_dirs_to_start_dir(df[[dir_str]], start_dir)
  df_dir_sector <- subset(df, df$tmp_var < sector_width)
  df_dir_sector$tmp_var <- NULL

  return(df_dir_sector)
}

filter_sims_for_dir_sector <- function(dirs, start_dir, sector_width) {
  #'
  #' filters directions to fiven sector
  #'
  #' @param dirs An array of directions
  #' @param start_dir direction to start sector from (clockwise)
  #' @param sector_width width in degrees

  #' @return idxs indices of matching directions for sector
  #'
  #' @examples
  #' idxs <- check_dir_sector(dirs, 270, 180)
  #'
  #' @export

  tmp_var <- adjust_dirs_to_start_dir(dirs, start_dir)
  idxs <- which(tmp_var < sector_width)

  return(idxs)
}

var_to_zscore <- function(df, var_str) {
  #' transforms variable to zscore
  #'
  #' @param df A dataframe.
  #' @return df A dataframe.
  #'
  #' @examples
  #' df <- var_to_zscore(ekofisk, "Pdir")
  #'
  #' @export

  vals <- df[[var_str]]
  tmp_mean <- mean(vals)
  tmp_std <- sd(vals)
  vals_zscore <- (vals - tmp_mean) / tmp_std
  var_zscore_str <- paste(var_str, "_", "zscore", sep = "")
  df[[var_zscore_str]] <- vals_zscore

  return(df)
}

transform_var_to_rank_df <- function(df, var_str, max_val) {
  #' transforms variable to ranks
  #'
  #' @param df A dataframe.
  #' @return df A dataframe.
  #'
  #' @examples
  #' df <- transform_var_to_rank(ekofisk, "Pdir", 360)
  #'
  #' @export

  vals <- df[[var_str]]
  vals_trans <- transform_var_to_rank_array(vals, max_val)
  var_trans_str <- paste(var_str, "_", "trans", sep = "")
  df[[var_trans_str]] <- vals_trans
  return(df)
}

transform_var_to_rank_array <- function(vals, max_val) {
  #' transforms variable to ranks
  #'
  #' @param vals input values
  #' @param max_val maximum value fro transformation
  #'
  #' @return vals_trans transformed output values
  #'
  #' @examples
  #' vals_trans <- transform_var_to_rank_array(vals, 360)
  #'
  #' @export

  vals_ranked <- rank(vals)
  vals_trans <- array(0, c(length(vals))) * NA
  for (i in seq_along(vals)){
    vals_trans[i] <- (max_val / length(vals)) * (vals_ranked[i] - 1)
  }
  return(vals_trans)
}

transform_var_to_rank <- function(vals, max_val) {
  #' transforms variable to ranks
  #'
  #' @param vals input values
  #' @param max_val maximum value fro transformation
  #'
  #' @return vals_trans transformed output values
  #'
  #' @examples
  #' vals_trans <- transform_var_to_rank_array(vals, 360)
  #'
  #' @export

  vals_ranked <- rank(vals)
  vals_trans <- array(0, c(length(vals))) * NA
  for (i in seq_along(vals)){
    vals_trans[i] <- (max_val / length(vals)) * (vals_ranked[i] - 1)
  }
  return(vals_trans)
}

unfold_counts <- function(dfin, covarstr_lst) {
  #' Unfolds dataframe with covariates for binned count data
  #'
  #' @param df A dataframe.
  #' @return df A dataframe.
  #'
  #' @examples
  #' df <- data.frame(counts=seq(10), cov1=seq(10), cov2=seq(10))
  #' unfolded <- unfold_counts(df)
  #'
  #' @export

  covars <- list()
  # possibly outsource for loop in apply function and import here
  for (s in covarstr_lst) {
    covar <- list()
    for (i in seq_len(dim(dfin)[1])){
      covar[[i]] <- rep(dfin[i, 2], dfin[i, 1])
    }
    covars[[s]] <- unlist(covar)
  }
  covars <- as.data.frame(covars)
  colnames(covars) <- covarstr_lst
  covars[["counts"]] <- array(1, dim(covars)[1])
  return(covars)
}

make_data_frame <- function(varstr_lst, var_lst) {
  # make dataframe
  df <- data.frame(var_lst[[1]])
  colnames(df) <- c(varstr_lst[[1]])
  for (i in 2:length(varstr_lst)) {
    # Add a new column to the dataframe
    df[[varstr_lst[[i]]]] <- var_lst[[i]]
  }
  return(df)
}

closest_even <- function(x) {

  if (x %% 2 == 0) {
    # If x is already even, it's the closest even number
    return(x)
  } else {
    # If x is odd, find the closest even number
    return(round(x / 2) * 2)
  }
}

setup_pp_wts <- function(dfin, varstr_lst, bounds = NULL, res = NULL,
                         spltype = NULL, nquad = 256) {
  #' setting up weights for point process and ppgam
  #'
  #' @param dfin input dataframe
  #' @param varstr_lst list of variables to consider
  #' @param bounds boundaries for each variable
  #' @param res resolution for discretizing
  #' @param spltype for example 'cc', 'cs', ...
  #' @nquad number of quadrature points
  #'
  #' @return list of nodes (contain mids and wts), knots, and nquad
  #'
  #' @examples
  #' setup_pp_res <- setup_pp_wts(data_sub_orig$y, varstr_lst=c('doy','Pdir'),
  #'                              bounds=list('doy'=c(0,366), 'Pdir'=c(0,360)),
  #'                              res=c(36,36), spltype=c('cc','cc'), nquad=256)
  #'
  #' @export
  #'

  intervals <- list()
  llims <- list()
  mids <- list()
  breaks <- list()
  histres <- list()
  nodes <- list()
  wts <- list()
  knots <- list()

  for (i in seq_along(varstr_lst)){
    intervals[[varstr_lst[i]]] <- seq(bounds[[varstr_lst[i]]][1], bounds[[varstr_lst[i]]][2], length.out = res[i])
    width <- intervals[[varstr_lst[i]]][2] - intervals[[varstr_lst[i]]][1]
    mids[[varstr_lst[i]]] <- intervals[[varstr_lst[i]]][1:(length(intervals[[varstr_lst[i]]]) - 1)] + width / 2
    llims[[varstr_lst[i]]] <- mids[[varstr_lst[i]]] - width / 2
    breaks[[varstr_lst[i]]] <- c(llims[[varstr_lst[i]]], llims[[varstr_lst[i]]][length(llims[[varstr_lst[i]]])] + width)

    if (spltype[i] == "cc") {
      # if cyclic
      knots[[varstr_lst[i]]] <- c(breaks[[varstr_lst[i]]][1], breaks[[varstr_lst[i]]][length(breaks[[varstr_lst[i]]])])
    } else {
      # place knots at breaks or evenly or according to customized locations
      print("Not yet available")
    }

    # compute histograms for wts
    histres[[varstr_lst[i]]] <- hist(dfin[[varstr_lst[i]]], breaks = breaks[[varstr_lst[i]]], plot = FALSE)
    # compute wts
    wts[[varstr_lst[i]]] <- histres[[varstr_lst[i]]]$counts / sum(histres[[varstr_lst[i]]]$counts)

    # define nodes
    nodes[[varstr_lst[i]]] <- cbind(mids[[varstr_lst[i]]], wts[[varstr_lst[i]]])
  }
  return(list(nodes = nodes, knots = knots, nquad = nquad))
}

calc_q <- function(maxvallst, p) {
  #' Calculates the return value estimates with a given probability
  #' based on exp(-1) or 1-1/T such as q2 (uses A max), q2N (uses AN max),
  #' q3 (uses A max), q4 (uses AN max).
  #' - only q2 and q2N are based on averages of realizations of distributions
  #' of uncertain parameters.
  #' - q3 and q4 are applied to the combined uncertain distribution.
  #' For more info see Jonathan et al. (2021).
  #'
  #' @export
  #'
  for (i in seq_along(maxvallst)) {
    rvs[[i]] <- quantile(unlist(maxvallst[[i]]), p)
  }
  return(rvs)
}

divide_data_into_k <- function(dfin, k) {
  #' divides a dataframe into k datasets
  #' @param dfin input dataframe
  #' @param k integer nbumber of k
  #'
  #' @return dfout_lst returns a list of dataframes, each with a subset
  #'
  #' @export
  idx <- seq_len(dim(dfin)[1])
  chunks <- split(idx, cut(seq_along(idx), k, labels = FALSE))
  dfout_lst <- NULL
  for (i in 1:k) {
    dfout_lst[[i]] <- dfin[chunks[[k]], ]
  }
  return(dfout_lst)
}


cross_validation <- function(dfin, nr_cv, extr_thr_lst, varstr,
                             model_fml_thr, model_fml_gpd, knots,
                             tailfrac = .1) {
  #' perform a k-fold cross-validation scheme
  #'
  #' @export

  mean_sample_size_lst <- NULL

  chunked_df <- divide_data_into_k(dfin, nr_cv)

  cv_bias_thr_lst <- NULL
  cv_mae_thr_lst <- NULL
  cv_serror_thr_lst <- NULL
  cv_dserror_thr_lst <- NULL
  for (j in seq_along(extr_thr_lst)) {
    crossres <- try({
      extr_thr <- list()
      extr_thr[[varstr]] <- extr_thr_lst[j]
      cv_bias_lst <- NULL
      cv_mae_lst <- NULL
      cv_serror_lst <- NULL
      cv_dserror_lst <- NULL

      mean_sample_size_tmp <- NULL

      for (i in 1:nr_cv){
        print("###")
        print(c("extr_thr:", extr_thr))
        print(c("cross-validation step:", i))
        print("###")

        chunk_train_lst <- NULL
        for (c in 1:nr_cv){
          if (c != i) {
            chunk_train_lst[[c]] <- chunked_df[[c]]
          }
        }
        chunk_train_df <- do.call(rbind, chunk_train_lst)
        chunk_test_df <- chunked_df[[i]]

        margs_thr_orig <- NULL
        if (i == 1) {
          print("fit threshold model")
          margs_thr_orig[[varstr]] <- evgam(model_fml_thr[[varstr]], chunk_train_df,
                                            family = "ald",
                                            ald.args = list(tau = extr_thr[[varstr]]),
                                            knots = knots)
        } else {
          print("fit threshold model")
          margs_thr_orig[[varstr]] <- evgam(model_fml_thr[[varstr]], chunk_train_df,
                                            family = "ald",
                                            ald.args = list(tau = extr_thr[[varstr]]),
                                            knots = knots,
                                            sp = margs_thr_orig[[varstr]]$sp)
        }

        data_sub_orig <- subset_df(margs_thr_orig, thr_str = "thr", exc_str = "exc")

        margs_gpd_orig <- fit_marginal_models_gpd(dfin = data_sub_orig,
                                                  model_fml = model_fml_gpd,
                                                  list_var = varstr,
                                                  knots = knots)

        df_sorted <- chunk_test_df[order(chunk_test_df[[varstr]], decreasing = TRUE), ][1:as.integer(tailfrac * dim(chunk_test_df)[1]), ]

        mean_sample_size_tmp[[i]] <- dim(df_sorted)[1]
        ### predict explicitly for the left out row ###
        # predict threshold
        thr_pred <- predict(margs_thr_orig[[varstr]], newdata = df_sorted, type = "response")$location
        # predict exceedance
        gpd_param_preds <- predict(margs_gpd_orig[[varstr]], newdata = df_sorted, type = "response")
        scales <- gpd_param_preds$scale
        shapes <- gpd_param_preds$shape
        gpd_pred <- revd(length(scales), scale = scales, shape = shapes, threshold = thr_pred, type = "GP")

        # compute error
        cv_bias_tmp <- (gpd_pred) - df_sorted[[varstr]]
        cv_serror_tmp <- ((gpd_pred) - df_sorted[[varstr]])**2
        cv_bias_lst[[i]] <- mean(unlist(cv_bias_tmp))
        cv_mae_lst[[i]] <- mean(unlist(abs(cv_bias_tmp)))
        cv_serror_lst[[i]] <- sqrt(mean(unlist(cv_serror_tmp)))
        cv_dserror_lst[[i]] <- sqrt((mean(unlist(cv_serror_tmp))) - cv_bias_lst[[i]]**2)
      }
      cv_bias_thr_lst[[j]] <- unlist(cv_bias_lst)
      cv_mae_thr_lst[[j]] <- unlist(cv_mae_lst)
      cv_serror_thr_lst[[j]] <- unlist(cv_serror_lst)
      cv_dserror_thr_lst[[j]] <- unlist(cv_dserror_lst)
      mean_sample_size_lst[[j]] <- mean_sample_size_tmp
    }, silent = TRUE) # Suppress error messages
    if (inherits(crossres, "try-error")) {
      cv_bias_thr_lst[[j]] <- NA
      cv_mae_thr_lst[[j]] <- NA
      cv_serror_thr_lst[[j]] <- NA
      cv_dserror_thr_lst[[j]] <- NA
    }
  }
  cv_bias_df <- as.data.frame(do.call(rbind, cv_bias_thr_lst))
  cv_mae_df <- as.data.frame(do.call(rbind, cv_mae_thr_lst))
  cv_serror_df <- as.data.frame(do.call(rbind, cv_serror_thr_lst))
  cv_dserror_df <- as.data.frame(do.call(rbind, cv_dserror_thr_lst))

  cv_mean_bias <- apply(cv_bias_df, 1, mean, na.rm = TRUE)
  cv_mean_mae <- apply(cv_mae_df, 1, mean, na.rm = TRUE)
  cv_mean_serror <- apply(cv_serror_df, 1, mean, na.rm = TRUE)
  cv_mean_dserror <- apply(cv_dserror_df, 1, mean, na.rm = TRUE)

  costfunc1 <- ((cv_mean_bias)^2 + cv_mean_serror)
  costfunc2 <- ((cv_mean_bias)^2 + cv_mean_dserror)
  costfunc3 <- (abs(cv_mean_bias) + cv_mean_mae)

  mask <- !is.na(costfunc1)
  extr_thr_lst <- extr_thr_lst[!is.na(costfunc1)]

  costfunc1 <- costfunc1[!is.na(costfunc1)]
  costfunc2 <- costfunc2[!is.na(costfunc2)]
  costfunc3 <- costfunc3[!is.na(costfunc3)]

  cost <- NULL
  cost[["costfct1"]] <- costfunc1
  cost[["costfct2"]] <- costfunc2
  cost[["costfct3"]] <- costfunc3

  errors <- NULL
  errors[["bias"]] <- cv_bias_df[mask, ]
  errors[["mae"]] <- cv_mae_df[mask, ]
  errors[["rse"]] <- cv_serror_df[mask, ]
  errors[["drse"]] <- cv_dserror_df[mask, ]

  mean_sample_size <- mean(unlist(mean_sample_size_lst))
  return(list("cost" = cost, "errors" = errors, "thr" = extr_thr_lst,
              "mean_sample_size" = mean_sample_size, "tailfrac" = tailfrac))
}

simple_prediction_thr <- function(dfin, extr_thr, varstr,
                                  model_fml_thr, knots) {
  #' @export
  print("fit threshold model")
  margs_thr_orig <- NULL
  margs_thr_orig[[varstr]] <- evgam(model_fml_thr[[varstr]], dfin,
                          family = "ald",
                          ald.args = list(tau = extr_thr),
                          knots = knots)

  data_sub_orig <- subset_df(margs_thr_orig, thr_str = "thr", exc_str = "exc")

  # predict threshold
  thr_pred <- predict(margs_thr_orig[[varstr]], newdata = data_sub_orig[[varstr]],
                      type = "response")$location

  return(list('thr_pred'=thr_pred, 'data_sub'=data_sub_orig))
}

simple_prediction_gpd <- function(dfin, varstr,
                                  model_fml_gpd, knots) {
  #' @export
  # fit gpd
  margs_gpd_orig <- fit_marginal_models_gpd(dfin = dfin,
                                            model_fml = model_fml_gpd,
                                            list_var = varstr,
                                            knots = knots)
  # predict exceedance)
  gpd_param_preds <- predict(margs_gpd_orig[[varstr]], newdata = dfin[[varstr]], type = "response")
  scales <- gpd_param_preds$scale
  shapes <- gpd_param_preds$shape
  gpd_pred <- revd(as.numeric(length(scales)), scale = scales, shape = shapes, threshold = 0, type = "GP")

  return(list('gpd_pred'=gpd_pred))
}

cross_validation_distr_simple <- function(dfin, extr_thr, varstr,
                                          model_fml_thr, model_fml_gpd,
                                          knots) {

  res_thr <- simple_prediction_thr(dfin, extr_thr, varstr,
                                   model_fml_thr, knots)

  res_gpd <- simple_prediction_gpd(res_thr$data_sub, varstr,
                                   model_fml_gpd, knots)

  var_pred <- res_thr$thr_pred + res_gpd$gpd_pred

  cv_score_lst <- crps_samples(var_pred, res_thr$data_sub[[varstr]][[varstr]])

  return(list('score'=cv_score_lst, 'gpd_pred'=var_pred, 'data'=res_thr$data_sub[[varstr]][[varstr]]))
}

cross_validation_distr <- function(dfin, nr_cv, extr_thr_lst, varstr,
                                   model_fml_thr, model_fml_gpd,
                                   knots, tail_cut=.1) {

  chunked_df <- divide_data_into_k(dfin, nr_cv)

  all_score_lst <- NULL
  rt_score_lst <- NULL
  crps_score_lst <- NULL

  for (i in 1:length(extr_thr_lst)){
    extr_thr <- extr_thr_lst[i]

    all_score <- NULL
    rt_score <- NULL
    crps_score <- NULL

    for (j in 1:nr_cv){
      print("###")
      print(c("extr_thr:", extr_thr))
      print(c("cross-validation step:", j))
      print("###")

      chunk_train_lst <- NULL
      res_extr_lst_tmp <- NULL
      for (c in 1:nr_cv){
        if (c != j) {
          chunk_train_lst[[c]] <- chunked_df[[c]]
        }
      }
      chunk_train_df <- do.call(rbind, chunk_train_lst)
      chunk_test_df <- chunked_df[[j]]

      cvres_distr <- cross_validation_distr_simple(dfin = dfin,
                                                   extr_thr = extr_thr_lst[i],
                                                   varstr = varstr,
                                                   model_fml_thr = model_fml_thr,
                                                   model_fml_gpd = model_fml_gpd,
                                                   knots = knots)

      cvres_ts <- tail_score(cvres_distr$gpd_pred, cvres_distr$data, tail_cut = tail_cut)

      all_score[[j]] <- cvres_ts$score
      rt_score[[j]] <- cvres_ts$right_score
      crps_score[[j]] <- cvres_distr$score
    }
    all_score_lst[[i]] <- unlist(all_score)
    rt_score_lst[[i]] <- unlist(rt_score)
    crps_score_lst[[i]] <- unlist(crps_score)
  }
  return(list('all_score'=all_score_lst, 'rt_score'=rt_score_lst, 'crps_score'=crps_score_lst, 'thr'=extr_thr_lst))
}

# ============================================================
# Distribution Comparison Scores
# - crps_samples():     CRPS between two empirical distributions
# - tail_score():       Tail similarity score (Anderson-Darling-based)
# ============================================================


# ------------------------------------------------------------
# CRPS for two empirical distributions (sample vs sample)
#
# Uses the energy-score identity:
#   CRPS(F, G) = E|X - Y| - 0.5 * E|X - X'| - 0.5 * E|Y - Y'|
#
# Args:
#   x        : numeric vector — "forecast" / reference samples
#   y        : numeric vector — "observation" / target samples
#   method   : "pwm"  → O(n log n), recommended (default)
#              "exact"→ O(n²),  brute-force pairwise; use for small n
#
# Returns: scalar CRPS value (lower = more similar; 0 = identical)
# ------------------------------------------------------------
crps_samples <- function(x, y, method = "pwm") {
  #' @export

  stopifnot(is.numeric(x), is.numeric(y),
            as.numeric(length(x)) >= 2, as.numeric(length(y)) >= 2)

  method <- match.arg(method, c("pwm", "exact"))

  # Internal helper: E|X - X'| via sorted PWM (O(n log n))
  # Identity: E|X-X'| = 2 * sum_i x_(i) * (2i - n - 1) / (n*(n-1))
  mean_abs_diff_pwm <- function(v) {
    n <- as.numeric(length(v))
    v <- sort(v)
    i <- seq_len(n)
    2 * sum(v * (2 * i - n - 1)) / (n * (n - 1))
  }

  # Internal helper: E|X - X'| brute force (O(n²))
  mean_abs_diff_exact <- function(v) {
    mean(abs(outer(v, v, "-")))
  }

  if (method == "pwm") {
    # E|X - Y|: combine and use sorted trick across joint samples
    # (exact cross-term via sorting both arrays)
    cross_term <- mean_abs_cross(x, y)
    self_x     <- mean_abs_diff_pwm(x)
    self_y     <- mean_abs_diff_pwm(y)
  } else {
    cross_term <- mean(abs(outer(x, y, "-")))
    self_x     <- mean_abs_diff_exact(x)
    self_y     <- mean_abs_diff_exact(y)
  }

  cross_term - 0.5 * self_x - 0.5 * self_y
}

# O(n log n) cross-term E|X - Y| using sorted merge
mean_abs_cross <- function(x, y) {
  #' @export
  nx <- as.numeric(length(x)); ny <- as.numeric(length(y))
  x  <- sort(x);   y  <- sort(y)

  # Use the identity:
  #   sum_{i,j} |x_i - y_j| = sum_i x_i*(2*rank_in_merged(x_i) - nx - 1)
  #                          + sum_j y_j*(2*rank_in_merged(y_j) - ny - 1)
  # where rank_in_merged counts position among all (nx+ny) values.
  combined <- c(x, y)
  labels   <- c(rep(1L, nx), rep(2L, ny))
  ord      <- order(combined)
  combined <- combined[ord]
  labels   <- labels[ord]

  # Running counts of how many x and y values have been seen so far
  cx <- cumsum(labels == 1L)   # # of x values <= combined[k]
  cy <- cumsum(labels == 2L)   # # of y values <= combined[k]

  total <- 0
  for (k in seq_along(combined)) {
    if (labels[k] == 1L) {
      # x value: contributes x_k * (cy[k] below) - x_k * (ny - cy[k] above)
      # = x_k * (2*cy[k] - ny)  ... but cy[k] here counts y values < or = x_k
      total <- total + combined[k] * (2 * cy[k] - ny)
    } else {
      total <- total + combined[k] * (2 * cx[k] - nx)
    }
  }
  total / (nx * ny)
}


# ------------------------------------------------------------
# Tail Similarity Score
#
# Quantifies how similar the *tails* of two distributions are.
# Strategy: weighted KS-type statistic that up-weights tail regions.
#
# Two components are combined:
#   1. Left tail score  — focuses on quantiles [0, tail_cut]
#   2. Right tail score — focuses on quantiles [1-tail_cut, 1]
#
# Each component is the weighted L2 distance between the two ECDFs,
# with weight w(u) = 1 / (u * (1-u))  (Anderson-Darling weight),
# which diverges at 0 and 1, giving heavy emphasis to extremes.
#
# The final score is normalised to [0, 1]:  0 = identical tails, 1 = max diff.
#
# Args:
#   x         : numeric vector (distribution 1)
#   y         : numeric vector (distribution 2)
#   tail_cut  : quantile threshold defining "tail" region (default 0.10)
#   n_grid    : number of quantile grid points (default 500)
#   normalize : if TRUE (default), returns score in [0,1] via
#               comparison to a reference worst-case
#
# Returns: list with
#   $score       overall tail similarity (0 = identical, 1 = maximally different)
#   $left_score  left-tail component
#   $right_score right-tail component
#   $detail      data.frame of quantile grid and per-point contributions
# ------------------------------------------------------------
tail_score <- function(x, y,
                       tail_cut  = 0.10,
                       n_grid    = 500,
                       normalize = TRUE) {
  #' @export

  stopifnot(is.numeric(x), is.numeric(y),
            as.numeric(length(x)) >= 10, as.numeric(length(y)) >= 10,
            tail_cut > 0, tail_cut < 0.5)

  # Empirical quantile functions
  qx <- function(p) quantile(x, probs = p, type = 8)
  qy <- function(p) quantile(y, probs = p, type = 8)

  # Anderson-Darling weight (clamped to avoid Inf at boundaries)
  ad_weight <- function(u, eps = 1e-4) {
    u <- pmax(eps, pmin(1 - eps, u))
    1 / (u * (1 - u))
  }

  # Weighted L2 distance on a quantile grid restricted to [lo, hi]
  tail_distance <- function(lo, hi) {
    probs <- seq(lo, hi, length.out = n_grid)
    w     <- ad_weight(probs)
    w     <- w / sum(w)                     # normalise weights
    dx    <- qx(probs) - qy(probs)
    sum(w * dx^2)
  }

  left_raw  <- tail_distance(0.001,         tail_cut)
  right_raw <- tail_distance(1 - tail_cut,  0.999)
  combined  <- 0.5 * (left_raw + right_raw)

  if (normalize) {
    # Reference: compare x to its own mirror image shifted by 2*sd(x)
    # This gives a stable, data-adaptive worst-case baseline
    ref   <- c(x, 2 * mean(x) - x + 2 * sd(x))   # reflected + shifted copy
    ref_q <- function(p) quantile(ref, probs = p, type = 8)

    ref_dist <- function(lo, hi) {
      probs <- seq(lo, hi, length.out = n_grid)
      w     <- ad_weight(probs)
      w     <- w / sum(w)
      dx    <- qx(probs) - ref_q(probs)
      sum(w * dx^2)
    }
    ref_left  <- ref_dist(0.001,        tail_cut)
    ref_right <- ref_dist(1 - tail_cut, 0.999)
    ref_comb  <- 0.5 * (ref_left + ref_right)

    left_score  <- min(1, left_raw  / ref_left)
    right_score <- min(1, right_raw / ref_right)
    score       <- min(1, combined  / ref_comb)
  } else {
    left_score  <- left_raw
    right_score <- right_raw
    score       <- combined
  }

  # Detail data.frame for plotting / inspection
  probs  <- seq(0.001, 0.999, length.out = n_grid)
  is_tail <- (probs <= tail_cut) | (probs >= 1 - tail_cut)
  detail <- data.frame(
    quantile     = probs,
    q_x          = qx(probs),
    q_y          = qy(probs),
    diff         = qx(probs) - qy(probs),
    ad_weight    = ad_weight(probs),
    in_tail      = is_tail
  )

  list(
    score       = score,
    left_score  = left_score,
    right_score = right_score,
    detail      = detail
  )
}


# ------------------------------------------------------------
# Convenience: compare a list of distributions pairwise
#
# Args:
#   dists  : named list of numeric vectors
#   metric : "crps", "tail", or "both"
#   ...    : extra args forwarded to crps_samples() or tail_score()
#
# Returns: matrix (or list of matrices) of pairwise scores
# ------------------------------------------------------------
pairwise_scores <- function(dists, metric = "both", ...) {
  #' @export

  stopifnot(is.list(dists), as.numeric(length(dists)) >= 2)
  nms <- names(dists)
  if (is.null(nms)) nms <- paste0("D", seq_along(dists))
  n <- as.numeric(length(dists))

  make_mat <- function() {
    m <- matrix(0, n, n, dimnames = list(nms, nms))
    m
  }

  compute <- function(fn) {
    m <- make_mat()
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i != j) m[i, j] <- fn(dists[[i]], dists[[j]])
      }
    }
    m
  }

  if (metric == "crps") {
    return(compute(function(a, b) crps_samples(a, b, ...)))
  }
  if (metric == "tail") {
    return(compute(function(a, b) tail_score(a, b, ...)$score))
  }
  # both
  list(
    crps = compute(function(a, b) crps_samples(a, b, ...)),
    tail = compute(function(a, b) tail_score(a, b, ...)$score)
  )
}

#strip_margs <- function(margs_in) {
#  return margs_out
#}
