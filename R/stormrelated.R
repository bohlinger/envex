define_storms <- function(df, var_str, storm_thr, exc_var='excess'){
  #' index all storms, no storms has index 0
  #'
  #' @param df A dataframe.
  #' @param var_str variable str for target variable.
  #' @param storm_thr threshold for storm definition.
  #' @return df A dataframe.
  #'
  #' @examples
  #' df_storms <- define_storms(ekofisk[1:2000,], "hs", 2)
  #'
  #' @export

  df$storm_thr <- array(storm_thr, c(length(df[[var_str]])))
  df$storm_ts <- df[[var_str]] - df$storm_thr

  df$storm_idx <- array(0, c(length(df[[var_str]])))

  # detect blocks and index the storms and index all storms between the 0
  idx <- 1
  for (i in 1:length(df$storm_idx)){
    if (df[[exc_var]][i]>0){
      if (i>1){
        if (df[[exc_var]][(i-1)]<0){
          idx <- idx +1
        }
      }
      df$storm_idx[i]<-idx
    }
  }
  return(df)
}

#label_storms_variable_thr <- function(df, exc_str){
#  #' index all storms, no storms has index 0
#  #'
#  #' @param df A dataframe.
#  #' @param exc_str variable str for exceedance
#  #' @return df A dataframe.
#  #'
#  #' @export
#
#  df$storm_idx <- array(0, c(length(df[[exc_str]])))
#
#  # detect blocks and index the storms and index all storms between the 0
#  idx <- 1
#  for (i in 1:length(df$storm_idx)){
#    if (df[[exc_str]][i]>0){
#      if (i>1){
#        if (df[[exc_str]][(i-1)]<0){
#          idx <- idx +1
#        }
#      }
#      df$storm_idx[i]<-idx
#    }
#  }
#  return(df)
#}

label_storms_variable_thr <- function(df, exc_str) {
  #' Index all storms, no storms have index 0
  #'
  #' @param df A dataframe.
  #' @param exc_str Variable string for exceedance.
  #' @return df A dataframe.
  #'
  #' @export

  # Create a logical vector indicating exceedances
  exceedance <- df[[exc_str]] > 0

  # Use rle to identify runs of consecutive TRUE values
  rle_obj <- rle(exceedance)

  # Create a vector of storm indices
  storm_idx <- rep(0, length(exceedance)) # Initialize with zeros

  # Assign storm indices where exceedance is TRUE
  storm_ids <- cumsum(rle_obj$values) * rle_obj$values
  storm_idx <- inverse.rle(list(lengths = rle_obj$lengths, values = storm_ids))

  # Add the storm indices to the dataframe
  df$storm_idx <- storm_idx

  return(df)
}

find_storm_peaks_and_decluster <- function(df, var_str, time_str,
                                           min_nr_days=1, exc_var='excess',
                                           idx_str='storm_idx'){
  #' Decluster storm events by ensuring a minimum time distance to pass between
  #' storm peaks
  #'
  #' @param df A dataframe.
  #' @param var_str variable str for target variable.
  #' @param time_str variable str for time variable.
  #' @param min_nr_days decorrelation time scale in days.
  #' @return pots A dataframe of declustered POTs
  #'
  #' @examples
  #' df_storms <- define_storms(ekofisk[1:2000,], "hs", 2)
  #' df_pots <- find_storm_peaks(df_storms,
  #'                             var_str = "hs", time_str = "time",
  #'                             min_nr_days = 10)
  #'@export

  # rm all negative storm_ts ("non-storms")
  df_non_neg <- subset(df,df[[exc_var]]>0)

  # split according to storm_idx
  newdf_pos_storms <- split(df_non_neg, df_non_neg[[idx_str]])
  df_pots_clustered <- do.call(rbind, lapply(newdf_pos_storms,
                                             function(x) x[which.max(x[[exc_var]]),]))

  ## sort pots to later check time difference between pots
  df_pots_sorted <- df_pots_clustered[order(df_pots_clustered[[idx_str]]),]

  # decluster pots based on time difference of peaks given min_nr_days
  dt_ex <- as.POSIXct(df_pots_sorted[[time_str]], format="%Y-%m-%d %H:%M:%S", tz="UTC")
  date_diffs <- difftime(df_pots_sorted[[time_str]][2:length(dt_ex)],
                         df_pots_sorted[[time_str]][1:length(dt_ex)-1], units="days")
  df_date_diffs <- data.frame(date_diffs)
  df_date_diffs <- cbind(df_date_diffs, seq(1,length(date_diffs)))
  colnames(df_date_diffs) <- c('diffs','idx')

  ## determine idx with a shorter time difference than min_nr_days and disregard them
  ## difference needs to be larger than 0 to not exclude double counts due to bootstrapping
  idx_tdiff <- which(df_date_diffs$diffs<min_nr_days & df_date_diffs$diffs > 0)
  tmp_lst <- list()
  for (idx in idx_tdiff){
    tmpl <- c(idx,(idx + 1))
    # mark the minim to neglect later, only keep maximum of neighboring pots
    tmpidx <- tmpl[which(df_pots_sorted[[var_str]][tmpl] == min(df_pots_sorted[[var_str]][tmpl]))]
    tmp_lst <- append(tmp_lst, tmpidx)}
  df_pots_declustered <- df_pots_sorted

  ## neglect minima if necessary
  if (length(tmp_lst)>0){
    print("no minimum detected")
    df_pots_declustered[unlist(tmp_lst),][[var_str]] <- NA
    df_pots_declustered <- na.omit(df_pots_declustered)
  }

  return(df_pots_declustered)
}

find_storm_peaks <- function(dfin, var_str, time_str,
                             exc_var='exc', idx_str='storm_idx'){
  #' Decluster storm events by ensuring a minimum time distance to pass between
  #' storm peaks
  #'
  #' @param dfin A dataframe.
  #' @param var_str variable str for target variable.
  #' @param time_str variable str for time variable.
  #' @param min_nr_days decorrelation time scale in days.
  #' @return pots A dataframe of declustered POTs
  #'
  #' @examples
  #' df_storms <- define_storms(ekofisk[1:2000,], "hs", 2)
  #' df_pots <- find_storm_peaks(df_storms, var_str = "hs", time_str = "time")
  #'@export

  # rm all negative storm_ts ("non-storms")
  df_non_neg <- subset(dfin, dfin[[exc_var]] > 0)

  # split according to storm_idx
  newdf_pos_storms <- split(df_non_neg, df_non_neg[[idx_str]])
  df_pots_clustered <- do.call(rbind, lapply(newdf_pos_storms,
                                             function(x) x[which.max(x[[exc_var]]), ]))

  ## sort pots to later check time difference between pots
  df_pots_sorted <- df_pots_clustered[order(df_pots_clustered[[idx_str]]), ]
  return(df_pots_sorted)
}

estimate_nr_of_peaks <- function(df, m_ald, settings){
  #' Estimate nr of peaks in df given min number of days
  # '
  #' @param df A dataframe.
  #' @param m_ald evgam ald threshold model object
  #' @param settings an object containing:
  #' @param var_str variable str for target variable.
  #' @param time_str variable str for time variable.
  #' @param min_nr_days decorrelation time scale in days.
  #'
  #' @return pots A dataframe of declustered POTs
  #'
  #' @examples
  #' df_storms <- define_storms(ekofisk[1:2000,], "hs", 2)
  #' nr_of_peaks <- estimate_nr_of_peaks(df_storms, m_ald, settings)
  #'
  #' @export


  df$thr <- predict(m_ald, nsim = 1,  newdata = df, type='response')$location
  df_tmp <- df
  df_tmp$exc <- df_tmp[[settings$var_str]] - df_tmp$thr
  df_tmp <- subset(df_tmp, exc > 0)
  df_pots <- find_storm_peaks(df_tmp,
                              settings$var_str, settings$time_str,
                              settings$min_nr_days)
  return(dim(df_pots)[1])
}

estimate_nr_of_peaks_sector <- function(df, m_ald, settings){
  #' Estimate nr of peaks in df given min number of days
  # '
  #' @param df A dataframe.
  #' @param m_ald evgam ald threshold model object
  #' @param settings an object containing:
  #' @param var_str variable str for target variable.
  #' @param time_str variable str for time variable.
  #' @param min_nr_days decorrelation time scale in days.
  #'
  #' @return pots A dataframe of declustered POTs
  #'
  #' @examples
  #' df_storms <- define_storms(ekofisk[1:2000,], "hs", 2)
  #' nr_of_peaks_sector <- estimate_nr_of_peaks_sector(df_storms, m_ald, settings)
  #'
  #' @export


  df$thr <- predict(m_ald, nsim = 1,  newdata = df, type='response')$location
  df_tmp <- df
  df_tmp$exc <- df_tmp[[settings$var_str]] - df_tmp$thr
  df_tmp <- subset(df_tmp, exc > 0)
  df_pots <- find_storm_peaks(df_tmp,
                              settings$var_str, settings$time_str,
                              settings$min_nr_days)
  # filter for sector if given
  df_pots_sector <- filter_dir_sector(df_pots, "Pdir", settings$start_dir, settings$sector_width)
  return(dim(df_pots_sector)[1])
}

decluster_exceedances <- function(dfin, var_str, time_str, min_nr_days){
  #'
  #' Decluster exceedences to only retain declustered peaks
  #' where storms are separated by min_nr_days
  #'
  #' @param df A dataframe.
  #' @param var_str string of variable of interest
  #' @param time_str string of time variable
  #' @param min_nr_of_days minimum number of days for keeping storms apart

  #' @return df_pots a data frame of declustered POTs
  #'
  #' @examples
  #' df_pots <- decluster_exceedances(df_data, var_str = "hs",
  #'                                  time_str = "time", min_nr_days = 1)
  #'
  #' @export
  #'

  # create blocks from which only the respective maximum is retained
  dt_ex <- as.POSIXct(dfin[[time_str]], format="%Y-%m-%d %H:%M:%S", tz="UTC")
  date_diffs <- difftime(dfin[[time_str]][2:length(dt_ex)],dfin[[time_str]][1:length(dt_ex)-1], units="days")
  df_date_diffs <- data.frame(date_diffs)

  df_date_diffs <- cbind(df_date_diffs, seq(1,length(date_diffs)))
  colnames(df_date_diffs) <- c('diffs','idx')

  # define blocks by difference of minimal number of days
  mask <- df_date_diffs$diffs >= min_nr_days
  blck_idx <- df_date_diffs$idx[mask]
  end_idx <- c(blck_idx,length(dfin[[time_str]]))
  start_idx <- c(1,end_idx[1:length(end_idx)-1]+1)

  column_names <- names(dfin)
  newdf <- as.data.frame(matrix(ncol = length(column_names), nrow = length(blck_idx)))
  names(newdf) <- column_names

  # create df of POTs + covariates and get POTs (peaks from blocks)

  for (x in 1:length(blck_idx)) {
    tmpval <- max(dfin[[var_str]][start_idx[x]: end_idx[x]])
    tmpmask <- dfin[[var_str]]==tmpval
    tmpidx <- which(tmpmask[start_idx[x]:end_idx[x]]==TRUE)
    if (length(tmpidx)>1){
      print(start_idx[x]); print(end_idx[x]);
      #print("Storm is considered multiple times in a row, choosing the first one.");
      #tmpidx <- tmpidx[1]}
      print("Storm is considered multiple times in a row, choosing a random one.");
      tmpidx <- tmpidx[runif_func(1, min = 1, max = length(tmpidx))]}
    for (name in column_names){
      newdf[[name]][x] <- dfin[[name]][start_idx[x]:end_idx[x]][tmpidx]
    }
  }
  return(newdf)
}

get_list_of_close_storms <- function(dfin, var_str, time_str, min_nr_days){
  #'
  #' Decluster exceedences to only retain declustered peaks
  #' where storms are separated by min_nr_days
  #'
  #' @param df A dataframe.
  #' @param var_str string of variable of interest
  #' @param time_str string of time variable
  #' @param min_nr_of_days minimum number of days for keeping storms apart

  #' @return df_pots a data frame of declustered POTs
  #'
  #' @examples
  #'
  #' @export
  #'

  # create blocks from which only the respective maximum is retained
  dt_ex <- as.POSIXct(dfin[[time_str]], format="%Y-%m-%d %H:%M:%S", tz="UTC")
  date_diffs <- difftime(dfin[[time_str]][2:length(dt_ex)],dfin[[time_str]][1:length(dt_ex)-1], units="days")
  df_date_diffs <- data.frame(date_diffs)

  df_date_diffs <- cbind(df_date_diffs, seq(1,length(date_diffs)))
  colnames(df_date_diffs) <- c('diffs','idx')

  # define blocks by difference of minimal number of days
  mask <- df_date_diffs$diffs >= min_nr_days
  blck_idx <- df_date_diffs$idx[mask]
  end_idx <- c(blck_idx,length(dfin[[time_str]]))
  start_idx <- c(1,end_idx[1:length(end_idx)-1]+1)

  storm_idx_list = NULL
  for (x in 1:length(blck_idx)) {
    # create new df of storm peaks with new indices
    storm_idx_list[[x]] <- dfin[['storm_idx']][start_idx[x]: end_idx[x]]
  }
  return(storm_idx_list)
}

#get_list_of_close_storms_idx <- function(dfin, var_str, time_str, stepsize){
#  #'
#  #' Decluster exceedences to only retain declustered peaks
#  #' where storms are separated by min_nr_days
#  #'
#  #' @param df A dataframe.
#  #' @param var_str string of variable of interest
#  #' @param time_str string of time variable
#  #' @param stepsize minimum number of steps for keeping storms apart
#  #'
#  #' @export
#  #'
#
#  # create blocks from which only the respective maximum is retained
#  date_diffs <- dfin[[time_str]][2:dim(dfin)[1]]-dfin[[time_str]][1:(dim(dfin)[1]-1)]
#  df_date_diffs <- data.frame(date_diffs)
#
#  df_date_diffs <- cbind(df_date_diffs, seq(1,length(date_diffs)))
#  colnames(df_date_diffs) <- c('diffs','idx')
#
#  # define blocks by difference of minimal number of days
#  mask <- df_date_diffs$diffs >= stepsize
#  blck_idx <- df_date_diffs$idx[mask]
#  end_idx <- c(blck_idx,length(dfin[[time_str]]))
#  start_idx <- c(1,end_idx[1:length(end_idx)-1]+1)
#
#  storm_idx_list = NULL
#  #newdf <- dfin
#  for (x in 1:length(blck_idx)) {
#    # create new df of storm peaks with new indices
#    storm_idx_list[[x]] <- dfin[['storm_idx']][start_idx[x]: end_idx[x]]
#  }
#  return(storm_idx_list)
#}

get_list_of_close_storms_idx <- function(dfin, var_str, time_str, stepsize) {
  #'
  #' Decluster exceedances to only retain declustered peaks
  #' where storms are separated by a minimum number of steps
  #'
  #' @param dfin A dataframe.
  #' @param var_str String of variable of interest (not used in the function but kept for compatibility).
  #' @param time_str String of time variable.
  #' @param stepsize Minimum number of steps for keeping storms apart.
  #' @return A list of storm indices grouped into blocks.
  #'
  #' @export

  # Calculate the differences in time between consecutive rows
  date_diffs <- diff(dfin[[time_str]])

  # Identify the indices where the time difference is greater than or equal to stepsize
  block_boundaries <- which(date_diffs >= stepsize)

  # Define the start and end indices for each block
  start_idx <- c(1, block_boundaries + 1)
  end_idx <- c(block_boundaries, nrow(dfin))

  # Use mapply to extract the storm indices for each block
  storm_idx_list <- mapply(function(start, end) {
    dfin[['storm_idx']][start:end]
  }, start_idx, end_idx, SIMPLIFY = FALSE)

  return(storm_idx_list)
}

#combine_storms <- function(dfin, storm_idx_list, idx_str="storm_idx"){
#  dfnew <- dfin[0,]
#  colnames(dfnew) <- names(dfin)
#  for (i in 1:length(storm_idx_list)){
#    subdf <- subset(dfin,
#                    (dfin[[idx_str]]>=min(storm_idx_list[[i]]) & dfin[[idx_str]]<=max(storm_idx_list[[i]])))
#    subdf[[idx_str]] <- i
#    dfnew <- rbind(dfnew, subdf)
#  }
#  return(dfnew)
#}

combine_storms <- function(dfin, storm_idx_list, idx_str = "storm_idx") {
  #' Combine storms into a single dataframe with updated storm indices
  #'
  #' @param dfin A dataframe.
  #' @param storm_idx_list A list of storm indices to combine.
  #' @param idx_str Column name for storm indices.
  #' @return A dataframe with combined storms and updated indices.
  #'
  #' @export

  # Create a vector to store the new storm indices
  new_storm_idx <- numeric(nrow(dfin))

  # Assign a new storm index for each block in storm_idx_list
  for (i in seq_along(storm_idx_list)) {
    storm_indices <- unlist(storm_idx_list[[i]])
    new_storm_idx[dfin[[idx_str]] %in% storm_indices] <- i
  }

  # Add the new storm index to the dataframe
  dfin[[idx_str]] <- new_storm_idx

  # Filter rows that are part of the new storms
  dfin <- dfin[new_storm_idx > 0, ]

  return(dfin)
}

declust <- function(x, u, r = 1) {
  #' @export
  idu <- x > u
  exc <- x[idu]
  n_exc <- sum(idu)
  s <- which(idu)
  clst <- rep(1, n_exc)
  tu <- diff(s)
  if (n_exc > 1)
    clst[2:n_exc] <- 1 + cumsum(tu > r)
  out <- rep(NA, length(x))
  out[idu] <- clst
  out
}

bootstrap_all <- function(df){
  #' Bootstrap storm time series without breaking storms.
  #' Uses storm_idx variable in df for identifying storms and re-shuffling.
  #'
  #' @param df A data frame.
  #' @return df A bootstrapped un-ordered data frame.
  #'
  #' @examples
  #' bstrp_df <- bootstrap_storms(ekofisk)
  #'
  #' @export

  print("resample time series: storms and non-storms")

  d1 <- declust(df$storm_idx, .5)
  d2 <- declust(-df$storm_idx, -.5)
  d1[is.na(d1)] <- na.omit(d2) + .5

  df$storm_idx2 <- d1
  df_spl <- split(df, df$storm_idx2)
  n_boot <- length(df_spl)

  newdf_spl <- df_spl[sample(n_boot, n_boot, replace = TRUE)]
  for (i in 1:n_boot) newdf_spl[[i]]$storm_idx2 <- i
  newdf <- do.call(rbind, newdf_spl)

  return(newdf)
}

bootstrap_storms <- function(dfin, exc_var='exc'){
  #' Bootstrap storm time series without breaking storms.
  #' Uses storm_idx variable in df for identifying storms and re-shuffling.
  #'
  #' @param dfin A data frame.
  #' @return newdf A bootstrapped un-ordered data frame.
  #'
  #' @examples
  #' bstrp_df <- bootstrap_storms(ekofisk, exc_var='exc')
  #'
  #' @export

  print("resample storms")
  # remove valleys := zeros
  df_storms <- subset(dfin, dfin[[exc_var]]>0)
  df_spl <- split(df_storms, df_storms$storm_idx)
  n_boot <- length(df_spl)

  newdf_spl <- df_spl[sample(n_boot, n_boot, replace = TRUE)]
  for (i in 1:n_boot) newdf_spl[[i]]$storm_idx <- i
  newdf <- do.call(rbind, newdf_spl)

  return(newdf)
}

run_bootstrap_storms <- function(df, group_col, n_boot) {

  # Pre-compute once outside the loop
  unique_groups  <- unique(df[[group_col]])
  group_row_idx  <- split(seq_len(nrow(df)), df[[group_col]])

  # Pre-allocate results list
  boot_samples <- vector("list", n_boot)

  for (i in seq_len(n_boot)) {
    sampled_groups <- sample(unique_groups, size = length(unique_groups), replace = TRUE)
    row_idx        <- unlist(group_row_idx[as.character(sampled_groups)], use.names = FALSE)
    boot_samples[[i]] <- df[row_idx, ]
  }

  return(boot_samples)
}

run_bootstrap_storms_pots <- function(df, group_col, n_boot, max_var) {

  # Pre-compute once outside the loop
  unique_groups <- unique(df[[group_col]])
  group_row_idx <- split(seq_len(nrow(df)), df[[group_col]])

  # Pre-allocate results lists
  boot_samples  <- vector("list", n_boot)
  boot_max      <- vector("list", n_boot)

  for (i in seq_len(n_boot)) {
    sampled_groups <- sample(unique_groups, size = length(unique_groups), replace = TRUE)
    row_idx        <- unlist(group_row_idx[as.character(sampled_groups)], use.names = FALSE)
    boot_df        <- df[row_idx, ]

    # Extract rows corresponding to the max of max_var within each group
    max_vals <- tapply(boot_df[[max_var]], boot_df[[group_col]], max)
    max_idx  <- which(boot_df[[max_var]] == max_vals[as.character(boot_df[[group_col]])])
    boot_max_df <- boot_df[max_idx, ]

    # Identify and remove groups with missing values in max_var
    groups_with_na <- unique(boot_max_df[[group_col]][is.na(boot_max_df[[max_var]])])

    if (length(groups_with_na) > 0) {
      print("NA removed for bootstrap nr:")
      print(i)
      boot_max_df <- boot_max_df[!boot_max_df[[group_col]] %in% groups_with_na, ]
      boot_df     <- boot_df[!boot_df[[group_col]] %in% groups_with_na, ]
    }

    boot_samples[[i]] <- boot_df
    boot_max[[i]]     <- boot_max_df
  }

  return(list(
    boot_samples = boot_samples,
    boot_max     = boot_max
  ))
}

# total euclidian distance
total_eucdist_fct <- function(vardists_df){
  sum_df <- apply(vardists_df,1,sum)
  sqrt_df <- sqrt(sum_df)
  return(sqrt_df)
}

# circular distance
circ_diff <- function(a, b, L = 360) {
  ((a - b + L/2) %% L) - L/2
}

find_idx_closest_storm_peak <- function(sim, hists, idxs, varlst=NULL,
                                        dist_fct_lst=NULL, dist_circ_L_lst=NULL,
                                        sidx=1, eidx=10){
  #' @param simPeakMV multivariate simulated peak of Hs
  #' @param histPeaksMV multivariate historic peaks of Hs to match
  #' @param varlst list of matching variables, e.g. hs, s_tm1
  #' @param sidx start idx
  #' @param eidx end idx
  #'
  #' @return Array of closest storm peak idx
  #'
  #' @export

  hists_scaled <- scale(hists)
  scaled_center <- as.data.frame(t(attr(hists_scaled, "scaled:center")))
  scaled_scale <- as.data.frame(t(attr(hists_scaled, "scaled:scale")))
  hists_scaled_df <- as.data.frame(hists_scaled)

  if (is.null(varlst)){
    varlst <- names(sim)
  }

  vardists <- NULL
  for (n in varlst) {
    sim_scaled <- scale(
      sim[[n]],
      center = scaled_center[[n]],
      scale = scaled_scale[[n]]
    )
    if (dist_fct_lst[[n]]=='circ'){
      L = dist_circ_L_lst[[n]]
      vardists[[n]] <- circ_diff(hists_scaled_df[[n]], sim_scaled, L = L)**2
    } else {
      vardists[[n]] <- (hists_scaled_df[[n]] - sim_scaled)**2
    }
  }

  # make dataframe from list
  vardists_df <- as.data.frame(vardists)

  # Compute Euclidean distances
  distances <- apply(vardists_df, 1, function(row) {
    sqrt(sum(row))
  })

  # create indices
  idx_SP <- order(distances)[runif_func(1, min=sidx, max=eidx)] # random pick

  return(idx_SP)
}

find_idx_closest_storm_peaks <- function(simPeaksMV, histPeaksMV, varlst=NULL,
                                         dist_fct_lst=NULL, dist_circ_L_lst=NULL,
                                         sidx=1, eidx=10){
  #' @param simPeakMV multivariate simulated peak of Hs
  #' @param histPeaksMV multivariate historic peaks of Hs to match
  #' @param varlst list of matching variables, e.g. hs, s_tm1
  #' @param sidx start idx
  #' @param eidx end idx
  #'
  #' @return Array of closest storm peak idx
  #'
  #' @export

  idx_SP_vec <- NULL
  for (i in 1:length(simPeaksMV[[varlst[1]]])){
    simPeakMV <- simPeaksMV[i, ]
    idx_SP_vec[[i]] <- find_idx_closest_storm_peak(simPeakMV, histPeaksMV,
                                       idxs = histPeaksMV[['storm_idx']],
                                       varlst = varlst,
                                       dist_fct_lst = dist_fct_lst,
                                       dist_circ_L_lst = dist_circ_L_lst,
                                       sidx = sidx, eidx = eidx)
  }
  return(unlist(idx_SP_vec))
}

add_angle <- function(a, b, L=360) {
  (a + b) %% L
}

anchor_sim_storms <- function(df_joint_pots, df_hist_pots, df_hist_storms, matched_storm_idx){
  # scale and rotate storms

  var_lst <- c("hs", "tm2", "Pdir", "doy", "storm_idx")

  sim_storm_df <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(sim_storm_df) <- var_lst

  hist_storm_df <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(hist_storm_df) <- var_lst

  for (i in 1:length(matched_storm_idx)){
    tmplst <- NULL

    df_pots_filtered <- df_pick$pots[df_pick$pots$storm_idx %in% matched_storm_idx[i], ]
    df_storms_filtered <- df_pick$storms[df_pick$storms$storm_idx %in% matched_storm_idx[i], ]

    scaling <- get_constant_scaling(df_joint_pots$hs[i], df_pots_filtered$hs)
    dirdiff <- circ_diff(df_joint_pots$Pdir[i], df_pots_filtered$Pdir)

    rotated_angle <- add_angle(df_storms_filtered$Pdir, dirdiff, L=360)
    scaled_storm <- scale_storm(df_storms_filtered$hs, df_storms_filtered$tm2, scaling)

    doydiff <- circ_diff(df_joint_pots$doy[i], df_pots_filtered$doy, L=365)
    adjusted_doy <- add_angle(df_storms_filtered$doy, doydiff, L=365)

    tmplst[['hs']] <- scaled_storm[[1]]
    tmplst[['tm2']] <- scaled_storm[[2]]
    tmplst[['Pdir']] <- rotated_angle
    tmplst[['doy']] <- adjusted_doy

    storm_idx_tmp_array <- array(data = 1, dim = c(length(df_storms_filtered$hs)))
    tmplst[['storm_idx']] <- storm_idx_tmp_array * matched_storm_idx[i]
    tmplst[['pseudo_storm_idx']] <- storm_idx_tmp_array * i

    df_tmp <- as.data.frame(tmplst)
    sim_storm_df <- rbind(sim_storm_df, df_tmp)

    hist_storm <- df_storms_filtered[var_lst]
    hist_storm["pseudo_storm_idx"] <- storm_idx_tmp_array * i
    hist_storm["dt"] <- df_storms_filtered$dt

    hist_storm_df <- rbind(hist_storm_df, hist_storm)
  }
  return (list(sim_storms = sim_storm_df, hist_storms = hist_storm_df))
}
