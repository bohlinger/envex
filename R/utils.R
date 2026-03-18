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

  costfunc1 <- (abs(cv_mean_bias) + cv_mean_serror)
  costfunc2 <- (abs(cv_mean_bias) + cv_mean_dserror)
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
