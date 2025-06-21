create_settings_object <- function(settings){
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

runif_func <-function(n, min=0, max=1) sample(min:max, n, replace=T)

adjust_dirs_to_start_dir <- function(d, start_dir){
  #' @export
  #'
  d_new <- (d-start_dir)%%360
  return(d_new)
}

filter_dir_sector <- function(df, dir_str, start_dir, sector_width){
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

filter_sims_for_dir_sector <- function(dirs, start_dir, sector_width){
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

var_to_zscore <- function(df, var_str){
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
  vals_zscore <- (vals - tmp_mean)/tmp_std
  var_zscore_str <- paste(var_str, "_", "trans", sep = "")
  df[[var_zscore_str]] <- vals_zscore
  return(df)
}

transform_var_to_rank_df <- function(df, var_str, max_val){
  #' transforms variable to ranks
  #'
  #' @param df A dataframe.
  #' @return df A dataframe.
  #'
  #' @examples
  #' df <- transform_var_to_rank(ekofisk, "Pdir", 360)
  #'
  #' @export

  #vals <- sort(df[[var_str]])
  #vals_ranked <- rank(vals)
  vals <- df[[var_str]]
  vals_trans <- transform_var_to_rank_array(vals, max_val)
  var_trans_str <- paste(var_str, "_", "trans", sep = "")
  df[[var_trans_str]] <- vals_trans
  return(df)
}

transform_var_to_rank_array <- function(vals, max_val){
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
  for (i in 1:length(vals)){
    vals_trans[i] <- (max_val/length(vals)) * (vals_ranked[i] -1)
  }
  return(vals_trans)
}

transform_var_to_rank <- function(vals, max_val){
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
  for (i in 1:length(vals)){
    vals_trans[i] <- (max_val/length(vals)) * (vals_ranked[i] -1)
  }
  return(vals_trans)
}

unfold_counts <- function(dfin, covarstr_lst){
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
  for (s in covarstr_lst){
    covar <- list()
    for (i in 1:dim(dfin)[1]){
      covar[[i]] <- rep(dfin[i,2], dfin[i,1])
    }
    covars[[s]] <- unlist(covar)
  }
  covars <- as.data.frame(covars)
  colnames(covars) <- covarstr_lst
  covars[['counts']] <- array(1,dim(covars)[1])
  return(covars)
}

make_data_frame <- function(varstr_lst, var_lst){
  # make dataframe
  df <- data.frame(var_lst[[1]])
  colnames(df) <- c(varstr_lst[[1]])
  for (i in 2:length(varstr_lst)){
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
    return(round(x/2) * 2)
  }
}

#simulate.ppgam <- function(object, nsim = 1e3, seed = NULL, newdata,
#                           type = "link", ...) {
#  #' Add simulate function to ppgam
#  #'
#  #'@param object ppgam object
#  #'
#  #'@return X draws from posterior for distribution parameters
#  #'
#  #'@examples
#  #'sim_counts <- simulate.ppgam(m1, nsim=10, newdata = newd, type = 'response')
#  #'
#  #'@export

#  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
#    runif(1) # initialize the RNG if necessary
#  if(is.null(seed)) {
#    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
#  } else {
#    R.seed <- get(".Random.seed", envir = .GlobalEnv)
#    set.seed(seed)
#    RNGstate <- structure(seed, kind = as.list(RNGkind()))
#    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
#  }
#  family <- object$family
#  V.type <- "Vp"
#  B <- evgam:::.pivchol_rmvn(nsim, object$coefficients, object[[V.type]])
#  X <- mgcv:::predict.gam(object, newdata, type = "lpmatrix")
#  X <- X %*% B
#  if (type == "response")
#    X <- exp(X)
#  return(X)
#}

filter_output <- function(df){
  # filter results dataframe
  return(0)
}
