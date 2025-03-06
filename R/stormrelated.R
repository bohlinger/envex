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

label_storms_variable_thr <- function(df, exc_str){
  #' index all storms, no storms has index 0
  #'
  #' @param df A dataframe.
  #' @param exc_str variable str for exceedance
  #' @return df A dataframe.
  #'
  #' @examples
  #' df_storms <- define_storms(ekofisk[1:2000,], "exc")
  #'
  #' @export

  df$storm_idx <- array(0, c(length(df[[exc_str]])))

  # detect blocks and index the storms and index all storms between the 0
  idx <- 1
  for (i in 1:length(df$storm_idx)){
    if (df[[exc_str]][i]>0){
      if (i>1){
        if (df[[exc_str]][(i-1)]<0){
          idx <- idx +1
        }
      }
      df$storm_idx[i]<-idx
    }
  }
  return(df)
}

find_storm_peaks_and_decluster <- function(df, var_str, time_str,
                                           min_nr_days=2, exc_var='excess',
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
  idx_tdiff <- which(df_date_diffs$diffs<min_nr_days & df_date_diffs$diffs>0)
  tmp_lst <- list()
  for (idx in idx_tdiff){
    tmpl <- c(idx,(idx+1))
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
  df_non_neg <- subset(dfin,dfin[[exc_var]]>0)

  # split according to storm_idx
  newdf_pos_storms <- split(df_non_neg, df_non_neg[[idx_str]])
  df_pots_clustered <- do.call(rbind, lapply(newdf_pos_storms,
                                             function(x) x[which.max(x[[exc_var]]),]))

  ## sort pots to later check time difference between pots
  df_pots_sorted <- df_pots_clustered[order(df_pots_clustered[[idx_str]]),]
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

get_list_of_close_storms_idx <- function(dfin, var_str, time_str, stepsize){
  #'
  #' Decluster exceedences to only retain declustered peaks
  #' where storms are separated by min_nr_days
  #'
  #' @param df A dataframe.
  #' @param var_str string of variable of interest
  #' @param time_str string of time variable
  #' @param stepsize minimum number of steps for keeping storms apart

  #' @return df_pots a data frame of declustered POTs
  #'
  #' @examples
  #' df_pots <- decluster_exceedances(df_data, var_str = "hs",
  #'                                  time_str = "time", stepsize = 1)
  #'
  #' @export
  #'

  # create blocks from which only the respective maximum is retained
  date_diffs <- dfin[[time_str]][2:dim(dfin)[1]]-dfin[[time_str]][1:(dim(dfin)[1]-1)]
  df_date_diffs <- data.frame(date_diffs)

  df_date_diffs <- cbind(df_date_diffs, seq(1,length(date_diffs)))
  colnames(df_date_diffs) <- c('diffs','idx')

  # define blocks by difference of minimal number of days
  mask <- df_date_diffs$diffs >= stepsize
  blck_idx <- df_date_diffs$idx[mask]
  end_idx <- c(blck_idx,length(dfin[[time_str]]))
  start_idx <- c(1,end_idx[1:length(end_idx)-1]+1)

  storm_idx_list = NULL
  #newdf <- dfin
  for (x in 1:length(blck_idx)) {
    # create new df of storm peaks with new indices
    storm_idx_list[[x]] <- dfin[['storm_idx']][start_idx[x]: end_idx[x]]
  }
  return(storm_idx_list)
}

combine_storms <- function(dfin, storm_idx_list, idx_str="storm_idx"){
  dfnew <- dfin[0,]
  colnames(dfnew) <- names(dfin)
  for (i in 1:length(storm_idx_list)){
    subdf <- subset(dfin,
                    (dfin[[idx_str]]>=min(storm_idx_list[[i]]) & dfin[[idx_str]]<=max(storm_idx_list[[i]])))
    subdf[[idx_str]] <- i
    dfnew <- rbind(dfnew, subdf)
  }
  return(dfnew)
}

combine_storms_including_valleys <- function(dfin, storm_idx_list, idx_str="storm_idx"){
  dfnew <- dfin[0,]
  colnames(dfnew) <- names(dfin)
  for (i in 1:length(storm_idx_list)){
    subdf <- subset(dfin,
                    (dfin[[idx_str]]>=min(storm_idx_list[[i]]) & dfin[[idx_str]]<=max(storm_idx_list[[i]])))
    subdf[[idx_str]] <- i
    dfnew <- rbind(dfnew, subdf)
  }
  return(dfnew)
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
