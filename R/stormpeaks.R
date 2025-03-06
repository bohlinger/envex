# imports
library(ppgam)
library(evgam)
library(extRemes)
library(lubridate)
library(dplyr)
library(rlang)
library(MASS)

produce_storm_occurrences <- function(nr_of_events, RP,
                                      nr_of_years,
                                      model_nr_of_events,
                                      condition=NULL){
  #' Produce storm occurrences by producing the appropriate number
  #' of coinciding covariates, i.e. direction and day of year.
  #'
  #' @param nr_of_events nr of events in dataset (integer)
  #' @param RP return period (RP) (integer)
  #' @param nr_of_years nr of years the dataset covers (integer)
  #' @param model_nr_of_events ppgam model object for occurrences given the covariates "doy" and "Pdir"
  #' @return df dataset of unfolded set of covariates ready to be used for GPD or other model
  #'
  #' @examples
  #' storms <- produce_storm_occurrences(nr_of_events, RP, model_nr_of_events)
  #'
  #' @export

  # simulate storm peaks from poisson
  print(c('nr_of_events:', nr_of_events))
  lambda <- nr_of_events*RP
  nr_of_events_sim <- rpois(1, lambda=lambda)

  rdirs <- runif(nr_of_events_sim, min = 1, max = 360)
  rdays <- runif(nr_of_events_sim, min = 1, max = 365.25)

  newdf <- data.frame(cbind(rdays, rdirs))
  colnames(newdf) <- c('doy', 'Pdir')
  predcounts <- predict(model_nr_of_events, newdata = newdf, type = 'response')

  simcounts <- rpois(length(predcounts), lambda=predcounts)
  dftmp <- data.frame(counts=simcounts, cov1=rdays, cov2=rdirs)

  print(c("max predicted counts",max(predcounts)))
  print(c("max simulated counts",max(simcounts)))

  dfcounts_tmp <- subset(dftmp, counts>0)
  dfcounts_tmp_1 <- subset(dfcounts_tmp, counts==1)
  dfcounts_tmp_l1 <- subset(dfcounts_tmp, counts>1)
  if (dim(dfcounts_tmp_l1)[1]>0){
    dfcounts_1 <- unfold_counts(dfcounts_tmp_l1)
    dfcounts <- rbind(dfcounts_1,dfcounts_tmp_1[,c('cov1','cov2')])
  } else {
    dfcounts <- dfcounts_tmp_1[,c('cov1','cov2')]
  }

  colnames(dfcounts) <- c('doy', 'Pdir')
  rm(dftmp)

  # filter according to covariate constraints
  if (!is.null(condition)){
    condition_expr <- parse_expr(condition)
    dfcounts <- dfcounts %>% filter(!!condition_expr)
    nr_of_events_sim <- dim(dfcounts)[1]
  }

  print(c('nr_of_events_sim:', dim(dfcounts)[1]))
  return(dfcounts)
}

#' produce_storm_occurrences <- function(nr_of_events, RP,
#'                                       nr_of_years, model_nr_of_events,
#'                                       condition=NULL){
#'   #' Produce storm occurrences by producing the appropriate number
#'   #' of coinciding covariates, i.e. direction and day of year.
#'   #'
#'   #' @param nr_of_events nr of events in dataset (integer)
#'   #' @param RP return period (RP) (integer)
#'   #' @param nr_of_years nr of years the dataset covers (integer)
#'   #' @param model_nr_of_events ppgam model object for occurrences given the covariates "doy" and "Pdir"
#'   #' @return df dataset of unfolded set of covariates ready to be used for GPD or other model
#'   #'
#'   #' @examples
#'   #' storms <- produce_storm_occurrences(nr_of_events, RP, nr_of_years, model_nr_of_events)
#'   #'
#'   #' @export
#'
#'   # simulate nr of events from poisson
#'   nr_of_events_sim <- rpois(1, lambda=nr_of_events/nr_of_years*RP)
#'   print('nr_of_events_sim:')
#'   print(nr_of_events_sim)
#'
#'   # sample the covariate space
#'   rdirs <- runif(nr_of_events_sim, min = 1, max = 360)
#'   rdays <- runif(nr_of_events_sim, min = 1, max = 365.25)
#'
#'   newdf <- data.frame(cbind(rdays, rdirs))
#'   colnames(newdf) <- c('doy', 'Pdir')
#'   pred_count <- predict(model_nr_of_events, newdata = newdf, type = 'response')
#'
#'   # convert mean counts
#'   simcounts <- rpois(length(pred_count), lambda=pred_count/nr_of_years*RP)
#'   dftmp <- data.frame(counts=simcounts, cov1=rdays, cov2=rdirs)
#'
#'   # unfold the counts together with covariates
#'   # NOTE: the number of the unfolded matrix is being shorted according
#'   # to the mean occurrence to not produce too many values
#'   dfcounts <- unfold_counts(dftmp[1:(floor(nr_of_events_sim/max(1,mean(dftmp$counts)))),])
#'   colnames(dfcounts) <- c('doy', 'Pdir')
#'   rm(dftmp)
#'
#'   # filter according to covariate constraints
#'   if (!is.null(condition)){
#'     condition_expr <- parse_expr(condition)
#'     dfcounts <- dfcounts %>% filter(!!condition_expr)
#'     nr_of_events_sim <- dim(dfcounts)[1]
#'   }
#'
#'   idxs <- runif_func(nr_of_events_sim, min = 1, max = dim(dfcounts)[1])
#'   return(dfcounts[idxs,])
#' }

simulate_storm_peaks_boostrap <- function(df_storms, m_gpd, m_ald){
  #'
  #' @param df_storms A data frame unfolded according to count including with covariates.
  #' @param m_gpd evgam model object for GPD distribution.
  #' @param m_ald evgam model object for ALD distribution.
  #'
  #' @return list of simulated peaks, thresholds, and exceedances.
  #'
  #' @examples
  #' res <- estimate_storm_peaks(df_sims_nr_years, m_gpd, m_ald)
  #'
  #' @export

  # init output fields
  ald_sim <- array(0, c(dim(df_storms)[1]))*NA
  gpd_sim <- array(0, c(dim(df_storms)[1]))*NA
  peak_sim <- array(0, c(dim(df_storms)[1]))*NA
  counter <- 1

  gpd_param_preds <- predict(m_gpd, newdata = df_storms, type= "response")
  scale <- gpd_param_preds$scale
  shape <- gpd_param_preds$shape

  # threshold
  tmp_ald_sim <- predict(m_ald, newdata = df_storms, type='response')
  ald_sim <- tmp_ald_sim$location
  thr <- tmp_ald_sim$location

  # simulate exceedance given predicted parameters given covariates
  #gpd_sim <- revd(length(scale), scale = scale, shape = shape, threshold = 0, type="GP")
  exc <- revd(length(scale), scale = scale, shape = shape, threshold = 0, type="GP")

  # simulate peak directly
  peaks <- revd(length(scale), scale = scale, shape = shape, threshold = thr, type="GP")

  # combine threshold on excess
  #peaks_sim <- ald_sim + gpd_sim
  comb <- thr + exc

  # probabilities
  probs <- array(0, c(length(scale)))*NA
  for (i in 1:length(scale)){
    probs[i] <- m_ald$tau+(1-m_ald$tau)*pgpd(peaks[i], mu=thr[i],
                                             sigma=scale[i], xi=shape[i],
                                             lower.tail = TRUE)
  }

  #return(list(peaks_sim, ald_sim, gpd_sim))
  return(list("peaks"=peaks, "exc"=exc, "comb"=comb, "probs"=probs, "thr"=thr, "scale"=scale, "shape"=shape))
}

simulate_storm_peaks_bayes <- function(df, m_gpd, m_ald, nsims){
  #'
  #' @param df A data frame unfolded according to count including with covariates.
  #' @param m_gpd evgam model object for GPD distribution.
  #' @param m_ald evgam model object for ALD distribution.
  #' @param nsims number of posterior simulations.
  #'
  #' @return list of simulated peaks, thresholds, and exceedances.
  #'
  #' @examples
  #' res <- simulate_storm_peaks_bayes(df, m_gpd, m_ald, nsims)
  #'
  #' @export

  # print('simulate from gpd')
  gpd_param_sims <- simulate(m_gpd, nsim = nsims, newdata = df, type= "response")
  scales <- gpd_param_sims$scale
  shapes <- gpd_param_sims$shape

  # exceedance
  gpd_sims <- revd(dim(df)[1], scale = scales, shape = shapes, threshold = 0, type="GP")
  gpd_sims <- array(1, c(length(gpd_sims),1))*gpd_sims

  # threshold
  ald_sims_tmp <- simulate(m_ald, nsim = nsims, newdata = df, type='response')
  ald_sims <- ald_sims_tmp$location

  # combine threshold and excess
  peak_sims <- ald_sims + gpd_sims

  # store peak probabilities for later transition to laplace space
  peak_probs <- array(0, c(length(gpd_sims),1))*NA
  tau = m_ald$tau
  for (i in 1:length(gpd_sims)){
    peak_probs[i] <- (1-(tau*(1-tau)*pevd(gpd_sims[i], threshold = ald_sims[i],
                                          scale=scales[i], shape=shapes[i], type='GP',
                                          lower.tail = TRUE)))
  }

  return(list(peak_sims, ald_sims, gpd_sims, peak_probs))
}

# simulate_storm_peaks function
simulate_storm_peaks_bootstrap_wrapper <- function(settings,
                                                   df,
                                                   m_ald=NULL,
                                                   transform_covar=NULL,
                                                   trace=0){
  #' @param df data with covariates
  #' @param m_ald threshold model if existing for speed-up
  #' @param transform_covar transforms covariate into rank based variable for optimal spline knot distribution
  #' @param trace degree for debugging, 0 for none
  #'
  #' @param settings settings object with necessary params as follows:
  #' @param N_bstrp number of bootstrap samples.
  #' @param fmla_gpd formula for GPD model.
  #' @param fmla_ald formula for ALD model.
  #' @param var_str string of variable of interest.
  #' @param q_thr threshold for ALD threshold model as quantile.
  #' @param storm_thr threshold for storm definition in [m].
  #' @param RP return period in years.
  #' @param min_nr_days decorrelation time in days for declusering.
  #' @param nr_of_years number of years represented by dataset.
  #'
  #' @return list of RP return levels and ALD threshold model.
  #'
  #' @examples
  #' fmla_ald <- hs ~ te(doy, Pdir, bs=c("cc","cc"), k=c(6,6))
  #' fmla_gpd_pots <- list(exc ~ te(doy, Pdir, k=c(6,6), bs=c("cc","cc")),
  #'                       ~ te(doy, Pdir, k=c(6,6), bs=c("cc","cc")))
  #' settings_list <- list(N_bstrp = 2, var_str = "hs", time_str = "time",
  #'                       fmla_gpd = fmla_gpd, fmla_ald = fmla_ald, fmla_pp = fmla_pp,
  #'                       q_thr = .9, storm_thr = 2, nr_of_years = 42,
  #'                       RP = 100, min_nr_days = 2)
  #' settings <- create_settings_object(settings=settings_list)
  #' res <- simulate_storm_peaks_bootstrap_wrapper(settings, df = ekofisk)
  #'
  #' @export

  # initialize array for maxima
  RP_max_peaks <- array(0, c(settings$N_bstrp))*NA

  # define storms for bootstrapping
  df_storms <- define_storms(df, settings$var_str, settings$storm_thr)

  for (b in 1:settings$N_bstrp){
    print("boostrap nr:")
    print(b)

    # boostrap original data by alternating storm and no-storm
    bstrp_df <- bootstrap_storms(df_storms)

    # add time as numeric
    bstrp_df$dt <- as.numeric(as.POSIXct(bstrp_df$time,
                                         format="%Y-%m-%d %H:%M:%S",
                                         tz="UTC"))

    # ------------------------------------------------------------------------ #
    # Fit threshold model
    m_ald <- fit_threshold_model(settings=settings,
                                 df=bstrp_df,
                                 transform_covar = transform_covar,
                                 trace=trace,
                                 m_ald=m_ald)
    # ------------------------------------------------------------------------ #

    # exclude all storm peaks that are below gam ALD threshold
    bstrp_df$thr <- fitted(m_ald)$location
    bstrp_df$exc <- bstrp_df[[settings$var_str]] - bstrp_df$thr
    df_tmp <- subset(bstrp_df, exc > 0)

    # find peaks and decluster
    print("find peaks and decluster")
    df_pots <- find_storm_peaks(df_tmp,
                                settings$var_str, settings$time_str,
                                settings$min_nr_days)
    rm(df_tmp)

    # fit occurrence model
    m_occ <- ppgam(settings$fmla_pp,
                   data = df_pots,
                   weights = settings$nr_of_years*settings$nr_of_years,
                   trace = trace)
    #weights = settings$nr_of_years*settings$nr_of_years/settings$RP,

    # produce storms
    df_sims_nr_years <- produce_storm_occurrences(dim(df_pots)[1],
                                                  settings$RP,
                                                  settings$nr_of_years,
                                                  m_occ)

    # fit GPD model to exceedances
    if (!is.null(m_gpd)){
      print("initialize gpd model fitting with previous parameters for speed-up")
      m_gpd <- try(evgam(settings$fmla_gpd,
                         data=df_pots,
                         family="gpd",
                         trace = trace,
                         sp=m_gpd$sp))
    } else {
      m_gpd <- try(evgam(settings$fmla_gpd,
                         data=df_pots,
                         family="gpd",
                         trace = trace))
    }

    if (inherits(m_gpd, 'try-error')) {
      print("##### try-error #####")
      next}

    # simulate peaks
    print("Simulate peaks")

    res <- simulate_storm_peaks_boostrap(df_sims_nr_years, m_gpd, m_ald)

    # store maxima
    A_res <- res[[1]]
    RP_max_peaks[b] <- max(A_res)
  }
  return(list(RP_max_peaks, m_ald, m_gpd, m_occ))
}

simulate_storm_peaks_bayes_wrapper <- function(settings,
                                               df,
                                               m_ald=NULL,
                                               transform_covar=NULL,
                                               trace=0){
  #' @param df data with covariates
  #' @param m_ald threshold model if existing for speed-up
  #' @param transform_covar transforms covariate into rank based variable for optimal spline knot distribution
  #' @param trace degree for debugging, 0 for none
  #'
  #' @param settings settings object with necessary params as follows:
  #' @param fmla_gpd formula for GPD model.
  #' @param fmla_ald formula for ALD model.
  #' @param fmla_pp formula for PP model.
  #' @param var_str string of variable of interest.
  #' @param q_thr threshold for ALD threshold model as quantile.
  #' @param min_nr_days decorrelation time in days for declusering.
  #' @param nr_of_years number of years represented by dataset.
  #' @param nsims_thr nr of threshold simulations
  #' @param nsims_years nr of years to simulate
  #' @param nsims_pois nr of occurrence simulations
  #' @param nsims_gpd nr of GPD simulations
  #'
  #' @return list of RP return levels and ALD threshold model.
  #'
  #' @examples
  #' fmla_ald <- hs ~ te(doy, Pdir, bs=c("cc","cc"), k=c(6,6))
  #' fmla_gpd_pots <- list(exc ~ te(doy, Pdir, k=c(6,6), bs=c("cc","cc")),
  #'                       ~ te(doy, Pdir, k=c(6,6), bs=c("cc","cc")))
  #' settings_list <- list(nsims_thr = 2,
  #'                       nsims_pois = 10,
  #'                       nsims_gpd = 1,
  #'                       nsims_years = 100,
  #'                       var_str = "hs", time_str = "time",
  #'                       fmla_gpd = fmla_gpd,
  #'                       fmla_ald = fmla_ald,
  #'                       fmla_pp, = fmla_pp,
  #'                       q_thr = .9,
  #'                       nr_of_years = 42, RP = 100,
  #'                       min_nr_days = 2)
  #' settings <- create_settings_object(settings=settings_list)
  #' res <- simulate_storm_peaks_bayes_wrapper(settings, df = ekofisk)
  #' @export

  # -------------------------------------------------------------------------- #
  # Fit threshold model
  m_ald <- fit_threshold_model(settings=settings,
                               df=df,
                               transform_covar = transform_covar,
                               trace=trace,
                               m_ald=m_ald)
  # -------------------------------------------------------------------------- #

  ## decluster exceedances and obtain POTs with location mean prediction

  # temporary data frame
  df_tmp <- df

  df_tmp[['threshold']] <- fitted(m_ald)$location

  df_tmp[['excess']] <- df_tmp[[settings$var_str]] - df_tmp$threshold

  # subset dataset for use
  df_tmp <- subset(df_tmp, excess > 0)

  # decluster and retrieve POTs
  df_pots_mean <- decluster_exceedances(df_tmp,
                                        var_str = settings$var_str,
                                        time_str = settings$time_str,
                                        min_nr_days = settings$min_nr_days)
  # delete temporary data frame
  rm(df_tmp)

  ## determine number of storms in Dataset
  nsims_storms_real <- dim(df_pots_mean)[1]
  print("nsims_storms_real")
  print(nsims_storms_real)

  ## Determine number of storms to be simulated representative of nsims_years
  nsims_storms <- nsims_storms_real/settings$nr_of_years*settings$nsims_years
  print("nsims_storms")
  print(nsims_storms)

  ## initialize final array of return levels

  # Initialize an empty list of arrays covering all accompanying variables
  array_list <- list()

  # Number of arrays to create
  array_names <- colnames(df_pots_mean)
  num_arrays <- length(array_names)

  # Create a loop
  for (i in 1:num_arrays) {
    # Create an array
    array <- array(0, c(settings$nsims_thr, settings$nsims_pois, nsims_storms))*NA

    # Add the array to the list
    array_list[[array_names[i]]] <- array
  }

  # add peaks, threshold, exceedences

  array_list[["peaks"]] <- array
  array_list[["thresholds"]] <- array
  array_list[["exceedences"]] <- array
  array_list[["probabilities"]] <- array

  #A_peaks <- array(0, c(settings$nsims_thr, settings$nsims_pois, nsims_storms))*NA
  #A_thr <- array(0, c(settings$nsims_thr, settings$nsims_pois, nsims_storms))*NA
  #A_exc <- array(0, c(settings$nsims_thr, settings$nsims_pois, nsims_storms))*NA

  A_doy <- array(0, c(settings$nsims_thr, settings$nsims_pois, nsims_storms))*NA
  A_Pdir <- array(0, c(settings$nsims_thr, settings$nsims_pois, nsims_storms))*NA

  #A_real_peaks <- array(0, c(settings$nsims_thr, nsims_storms_real))*NA
  array_list[["real_peaks"]] <- array(0, c(settings$nsims_thr, nsims_storms_real))*NA
  #A_real_exc <- array(0, c(settings$nsims_thr, nsims_storms_real))*NA

  for (n_thr in 1:settings$nsims_thr) {
    print('Dataset from threshold simulation Nr.')
    print(n_thr)

    # find peaks and decluster
    print("simulate threshold")

    # temporary data frame
    df_tmp <- df
    df_tmp$thr <- simulate(m_ald, nsim = 1,
                           newdata = df_tmp,
                           type='response')$location
    df_tmp$exc <- df_tmp[[settings$var_str]] - df_tmp$thr
    # subset dataset for use
    df_tmp <- subset(df_tmp, exc > 0)

    # decluster and retrieve POTs
    print("find peaks, decluster and retrieve POTs")
    df_pots <- decluster_exceedances(df_tmp,
                                     var_str = settings$var_str,
                                     time_str = settings$time_str,
                                     min_nr_days = settings$min_nr_days)
    # delete temporary data frame
    rm(df_tmp)

    # fit occurrence model
    print("fit occurrence model")
    m_occ <- ppgam(settings$fmla_pp,
                   data = df_pots,
                   weights = settings$nr_of_years*settings$nr_of_years)
    #weights = settings$nr_of_years*settings$nr_of_years/settings$nsims_years)

    ## produce storms
    print("produce storms")
    nr_of_events_sim <- rpois(1, lambda=nsims_storms)

    # draw dirs and day of year for subsequent simulations
    rdirs <- runif(nr_of_events_sim, min = 1, max = 360)
    rdays <- runif(nr_of_events_sim, min = 1, max = 365.25)

    # simulate mean number of occurrences
    df_tmp <- data.frame(cbind(rdays,rdirs))
    colnames(df_tmp) <- c('doy','Pdir')
    sim_counts_ppgam <- simulate.ppgam(m_occ, nsim=settings$nsims_pois,
                                       newdata = df_tmp, type = 'response')
    rm(df_tmp)

    print("Enter Pois-loop")
    for (pc in 1:settings$nsims_pois){
      # print("simulate covariate data set")
      # simulate integer number of occurrences to be used for creation of covariate data set
      sim_counts <- sim_counts_ppgam[,pc]
      sim_counts_pois <- rpois(length(sim_counts),
                               lambda=sim_counts * settings$nsims_years)
      #lambda=sim_counts/settings$nr_of_years*settings$nsims_years)

      df_tmp <- data.frame(counts=sim_counts_pois, cov1=rdays, cov2=rdirs)
      #
      mean_multiplicator <- mean(df_tmp$counts)
      df_tmp_red <- df_tmp[1:(round(length(df_tmp$counts)/floor(mean_multiplicator))+1),]
      df_sim_counts_pois <- unfold_counts(df_tmp_red)
      #
      # create covariate data set
      #df_sim_counts_pois <- unfold_counts(df_tmp)
      colnames(df_sim_counts_pois) <- c('doy', 'Pdir')

      # delete temporary data frame
      rm(df_tmp)
      rm(df_tmp_red)

      # draw at random events matching the number of storms peaks given years to be simulated
      idxs <- runif_func(nr_of_events_sim, min = 1, max = dim(df_sim_counts_pois)[1])
      df_sims_nr_years <- df_sim_counts_pois[idxs,]

      # fit exceedance model
      # print("fit exceedance model")
      m_gpd <- evgam(settings$fmla_gpd, data=df_pots, family="gpd", trace=trace)

      df_pred <- df_sims_nr_years
      # simulate exceedances, thresholds, and compute peaks for df_sims_nr_years
      # !! Caution, if covariates where transformed this needs to be forwarded into here !!
      res <- simulate_storm_peaks_bayes(df_pred, m_gpd, m_ald,
                                        nsims = settings$nsims_gpd)

      A_peaks_res <- res[[1]]  # get peaks
      A_thr_res <- res[[2]]  # get thresholds
      A_exc_res <- res[[3]]  # get exceedances
      A_peaks_probs <- res[[4]]  # get probabilities

      # draw at random nsims_storms nr of simulated values
      idx1 <- runif_func(nsims_storms, min=1, max=dim(A_peaks_res)[1])
      array_list[["peaks"]][n_thr,pc,] <- t(A_peaks_res[idx1,])
      array_list[["thresholds"]][n_thr,pc,] <- t(A_thr_res[idx1,])
      array_list[["exceedences"]][n_thr,pc,] <- t(A_exc_res[idx1,])
      array_list[["probabilities"]][n_thr,pc,] <- t(A_peaks_probs[idx1,])
      A_doy[n_thr,pc,] <- df_sims_nr_years$doy[idx1]
      A_Pdir[n_thr,pc,] <- df_sims_nr_years$Pdir[idx1]
    }

    idx2 <- runif_func(nsims_storms_real, min=1, max=length(df_pots[[settings$var_str]]))
    array_list[["real_peaks"]][n_thr,] <- df_pots[[settings$var_str]][idx2]
  }

  nr_peaks <- nsims_storms_real
  average_annual_nr_of_storms <- nr_peaks/settings$nr_of_years
  length_all <- settings$nsims_thr*settings$nsims_pois*settings$nsims_gpd

  return(list("peaks"=array_list[["peaks"]],
              "thresholds"=array_list[["thresholds"]],
              "exceedences"=array_list[["exceedences"]],
              "probabilities"=array_list[["probabilities"]],
              "df_pots"=df_pots,
              "nsims_storms"=nsims_storms, "nsims_storms_real"=nsims_storms_real,
              "real_peaks"=array_list[["real_peaks"]],
              "cov_dir"=A_Pdir, "cov_doy"=A_doy,
              "average_annual_nr_of_storms"=average_annual_nr_of_storms,
              "m_ald"=m_ald, "m_gpd"=m_gpd, "m_occ"=m_occ))
}

simulate_storm_peaks <- function(settings,
                                 df,
                                 m_ald=NULL,
                                 strategy="Bayes",  # Bootstrap
                                 transform_covar=NULL,
                                 trace=0){
  #' @param df data with covariates
  #' @param m_ald threshold model if existing for speed-up
  #' @param strategy Bayes or Bootstrap
  #' @param transform_covar transforms covariate into rank based variable for optimal spline knot distribution
  #' @param trace degree for debugging, 0 for none
  #'
  #' @param settings settings object with necessary params as follows:
  #' @param N_bstrp number of bootstrap samples.
  #' @param fmla_gpd formula for GPD model.
  #' @param fmla_ald formula for ALD model.
  #' @param var_str string of variable of interest.
  #' @param q_thr threshold for ALD threshold model as quantile.
  #' @param storm_thr threshold for storm definition in [m].
  #' @param RP return period in years.
  #' @param min_nr_days decorrelation time in days for declusering.
  #' @param nr_of_years number of years represented by dataset.
  #' @param nsims_thr nr of threshold simulations
  #' @param nsims_years nr of years to simulate
  #' @param nsims_pois nr of occurrence simulations
  #' @param nsims_gpd nr of GPD simulations
  #'
  #' @return list output from subfcts
  #'
  #' @examples
  #' fmla_ald <- hs ~ te(doy, Pdir, bs=c("cc","cc"), k=c(6,6))
  #' fmla_gpd_pots <- list(exc ~ te(doy, Pdir, k=c(6,6), bs=c("cc","cc")),
  #'                       ~ te(doy, Pdir, k=c(6,6), bs=c("cc","cc")))
  #' settings_list <- list(N_bstrp = 2, var_str = "hs", time_str = "time",
  #'                       fmla_gpd = fmla_gpd, fmla_ald = fmla_ald,
  #'                       q_thr = .9, storm_thr = 2, nr_of_years = 42,
  #'                       RP = 100, min_nr_days = 2)
  #' settings <- create_settings_object(settings=settings_list)
  #' RVs <- simulate_storm_peaks(settings, df = ekofisk, strategy="Bootstrap")
  #'
  #' @export

  strategy_names <- c('Bootstrap', 'Bayes')
  strategy_list <- list(simulate_storm_peaks_bootstrap_wrapper,
                        simulate_storm_peaks_bayes_wrapper)

  strategy_idx <- which(strategy_names == strategy)

  res <- strategy_list[[strategy_idx]](settings, df = df,
                                       m_ald = m_ald,
                                       transform_covar = transform_covar,
                                       trace = trace)

  return(res)
}

get_ann_max <- function(df_all, list_of_years, var_str, year_str){
  #' @param df_all dataframe to use with all info included
  #' @param list_of_years vector of all years like seq(1980,2021)
  #' @param var_str string of variable of interest
  #' @param year_str string of variable for years
  #'
  #' @return ann_max vector of maxima per year
  #'
  #' @export
  ann_max <- array(0, length(seq(length(list_of_years))))*NA
  for (y in list_of_years) {
    df_sub <- subset(df_all, df_all[[year_str]] == y)[[var_str]]
    ann_max[y-(list_of_years[1]-1)] <- max(df_sub)
  }
  return(ann_max)
}

calc_Hs_RV <- function(df_all, marginal_Hs_model, RP, Hs_str, year_str){
  #' @param df_all dataframe to use with all info included
  #' @param marginal_Hs_model string of method e.g. "block"
  #' @param RP Return period e.g. 100
  #' @param Hs_str string of Hs
  #' @param year_str string of variable for years
  #'
  #' @return Hs_RV Return value of Hs for given RP
  #'
  #' @export
  if (marginal_Hs_model=="block"){
    # get annual max
    list_of_years <- seq(min(df_all[[year_str]]), max(df_all[[year_str]]))
    Hs_block_max <- get_ann_max(df_all, list_of_years, Hs_str, year_str)
    # fit GEV
    model_stat <- fevd(Hs_block_max)  # from extRemes package
    Hs_RV <- return.level(model_stat, return.period = c(RP),
                          do.ci=FALSE, alpha=0.02, make.plot=FALSE)
  }
  return(Hs_RV)
}

get_Hs_bins_and_fit_ln <- function(df_all, Hs_str, T_str){
  #' @param df_all dataframe to use with all info included
  #' @param Hs_str string of Hs variable
  #' @param Tp_str string of T variable
  #'
  #' @return vectors of ln_means and ln_stds of T for Hs
  #'
  #' @export
  # bin hs in bins of 1m and perform a log-norm fit to all periods (e.g. Tp, Tm01, Tm02, ...) in one bin
  array_size <- as.integer(max(df_all[[Hs_str]]))
  ln_means <- array(0, array_size)*NA
  ln_stds <- array(0, array_size)*NA
  for (i in 1:array_size){
    df_sub <- subset(df_all, df_all[[Hs_str]]>=(i-1) & df_all[[Hs_str]]<(i))
    fit_wb <- fitdistr(df_sub[[T_str]], 'lognormal')
    ln_means[i] <- fit_wb$estimate[1]
    ln_stds[i] <- fit_wb$estimate[2]
  }
  return(list(ln_means, ln_stds))
}

fitMe <- function(params, X, Y) {
  a0 <- params[1]
  a1 <- params[2]
  a2 <- params[3]

  #points(X, a0+a1*X**a2, type="l")
  #Sys.sleep(0.1) # slow it down to watch progress

  error <- sum((Y - (a0+a1*X**a2))^2)
  #print(error)
  return(error)
}

calc_Tp_RV <- function(ln_means, Hs_RV){
  #' @param ln_means log-normal fits of T variable for Hs bins
  #' @param Hs_RV Return value of Hs
  #'
  #' @return T_joint_RV Return value of T variable
  #'
  #' @export

  array_size <- length(ln_means)
  X <- seq(array_size)-.5
  Y <- ln_means[1:array_size]

  o <- optim(c(a0=1, a1=1, a2=1), X=X, Y=Y, fitMe)
  # o$par are the optimal parameters

  T_joint_RV <- exp(o$par[1]+o$par[2]*Hs_RV**o$par[3])
  return(T_joint_RV)
}

fit_RV_from_joint_distr_ref_DNVRP <- function(df_all, marginal_Hs_model, RP,
                                              Hs_str, T_str, year_str){
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


fit_threshold_model <- function(settings, type='ald', df=df,
                                m_ald=NULL, trace=0,
                                transform_covar = NULL){
  print("Fit threshold model:")

  df_ald <- df

  if (!is.null(m_ald)){
    print("initialize threshold model fitting with previous parameters for speed-up")
    m_ald <- evgam(settings$fmla_ald, df_ald, family="ald",
                   ald.args=list(tau=settings$q_thr),
                   sp = m_ald$sp, trace=trace)
  } else {
    m_ald <- evgam(settings$fmla_ald, df_ald, family="ald",
                   ald.args=list(tau=settings$q_thr),
                   trace=trace)
  }

  return(m_ald)
}

get_probability <- function(dfin, m_ald, m_gpd, var_str = 'hs', nsim=1){
  #' @export

  #ald_thr <- simulate(m_ald, nsim = nsim, newdata = df_pots, type='response')$location
  ald_thr <- predict(m_ald, newdata = dfin, type='response')$location

  dfin[['threshold']] <- ald_thr
  dfin[['exc']] <- dfin[[var_str]]-ald_thr

  dfgpd <- subset(dfin,exc>0)
  dfecdf <- subset(dfin,exc<=0)

  #gpd_params <- simulate(m_gpd, nsim = nsim, newdata = df_pots, type= "response")
  gpd_params <- predict(m_gpd, newdata = dfgpd, type= "response")
  scales <- as.vector(gpd_params$scale)
  shapes <- as.vector(gpd_params$shape)
  thr <- dfgpd[['threshold']]

  # convert to probabilities using gpd
  probs_gpd <- array(0, c(length(scales)))*NA
  tau = q_thr
  for (i in 1:length(scales)){
    #probs_gpd[i] <- tau+(1-tau)*pgpd(dfgpd[[var_str]][i], mu=thr[i] ,sigma=scales[i], xi=shapes[i],lower.tail = TRUE)
    #probs_gpd[i] <- tau+(1-tau)*pgpd(dfgpd[['exc']][i], mu=thr[i] ,sigma=scales[i], xi=shapes[i],lower.tail = TRUE)
    probs_gpd[i] <- tau+(1-tau)*pgpd(dfgpd[["exc"]][i], mu=0 ,sigma=scales[i], xi=shapes[i],lower.tail = TRUE)
  }

  # convert to probabilities using ecdf
  Fn_ecdf <- ecdf(dfin[[var_str]])
  probs_ecdf <- Fn_ecdf(dfecdf[[var_str]])

  # combine
  probs <- c(probs_ecdf, probs_gpd)
  #probs <- c(probs_gpd)
  return (probs)
}


fit_marginal_models_thr <- function(dfin, thr, model_fml, list_var=NULL, m_params=NULL){
  #' @export
  #'
  if (is.null(list_var)){
    list_var <- names(model_fml)
  }
  margs <- NULL
  for (n in list_var){
    print(c("fit threshold model for",n))
    if (is.null(m_params)){
      margs[[n]] <- evgam(model_fml[[n]], dfin, family="ald", ald.args=list(tau=thr))
    } else {
      margs[[n]] <- evgam(model_fml[[n]], dfin, family="ald", ald.args=list(tau=thr), sp=m_params[[n]]$sp)
    }
  }
  return (margs)
}


fit_marginal_models_gpd <- function(dfin, model_fml, list_var=NULL, m_params=NULL){
  #' @export
  #'
  margs <- NULL
  if (is.null(list_var)){
    list_var <- names(model_fml)
  }
  for (n in list_var){
    print(c("fit gpd model for",n))
    if (is.null(m_params)){
      margs[[n]] <- evgam(model_fml[[n]], dfin[[n]], family="gpd")
    } else {
      margs[[n]] <- evgam(model_fml[[n]], dfin[[n]], family="gpd", sp=m_params[[n]]$sp)
    }
  }
  return (margs)
}

fit_marginal_models_occ <- function(dfin, model_fml, nr_of_years, list_var=NULL){
  #' @export
  #'
  margs <- NULL
  if (is.null(list_var)){
    list_var <- names(model_fml)
  }
  for (n in list_var){
    print(c("fit occurrence model for",n))
    margs[[n]] <- ppgam(model_fml[[n]], data = dfin[[n]], weights = nr_of_years*nr_of_years)
  }
  return (margs)
}

fit_marginal_models_occ_v2_ben <- function(dfin, model_fml, nr_of_years, list_var=NULL, sp0 = -1, sp = sp){
  #' @export
  #'
  margs <- NULL
  if (is.null(list_var)){
    list_var <- names(model_fml)
  }
  for (n in list_var){
    print(c("fit occurrence model for",n))
    margs[[n]] <- ppgam(model_fml[[n]], data = dfin[[n]], weights = nr_of_years*dim(dfin[[1]])[1], sp0 = sp0, sp = sp)
  }
  return (margs)
}

fit_marginal_models_occ_v2 <- function(dfin, model_fml, nr_of_years, list_var=NULL, nquad=20){
  #' @export
  #'
  margs <- NULL
  if (is.null(list_var)){
    list_var <- names(model_fml)
  }
  for (n in list_var){
    print(c("fit occurrence model for",n))
    margs[[n]] <- ppgam(model_fml[[n]], data = dfin[[n]], weights = nr_of_years*dim(dfin[[1]])[1], nquad=nquad)
  }
  return (margs)
}

subset_df <- function(margs, thr_str='thr', exc_str='exc'){
  #' @export
  #'
  list_var <- names(margs)
  data_sub <- NULL
  for (n in list_var){
    print(c("subset dataset for",n))
    dfin <- margs[[n]]$data
    dfin[[thr_str]] <- fitted(margs[[n]])$location
    dfin[[exc_str]] <- dfin[[n]] - dfin[[thr_str]]
    data_sub[[n]] <- subset(dfin, dfin[[exc_str]] > 0)
  }
  return (data_sub)
}


peak_picking <- function(dfin, lst_vars, thr_model_fml,
                         thr=.5, var_str='hs', thr_str='thr',
                         exc_str='exc', time_str='time',
                         decorrelation_time_scale=2){
  #' @export

  print('apply threshold model over all data')
  m_ald <- evgam(thr_model_fml, dfin, family="ald",
                 ald.args=list(tau=thr))
  q_tmp <- fitted(m_ald)$location

  dfin[[thr_str]] <- fitted(m_ald)$location
  dfin[[exc_str]] <- dfin[[var_str]] - dfin[[thr_str]]

  print('label storms')
  dfin_labeled <- label_storms_variable_thr(dfin, exc_str = exc_str)

  print('reduce dataset to storms')
  lst_vars <- c(lst_vars, exc_str, thr_str)
  dfin_labeled_red <- dfin_labeled[lst_vars]
  colnames(dfin_labeled_red) <- lst_vars

  print('combine and relabel storm peaks that are too close')
  # 1. round of peak finding
  df_pots_orig <- find_storm_peaks(dfin_labeled_red, var_str = var_str,
                                   time_str = time_str, exc_var = exc_str)

  # 2. check if storm peaks too close, if true combine
  storm_idx_list <- get_list_of_close_storms(df_pots_orig, var_str, time_str,
                                             decorrelation_time_scale)

  dfin_labeled_storms <- subset(dfin_labeled_red, dfin_labeled_red$storm_idx>0)
  dfin_labeled_combined <- combine_storms(dfin_labeled_storms, storm_idx_list)

  # 3. find peaks again for new combined storms
  print('find peaks for combined storms')
  df_pots_orig_combined <- find_storm_peaks(dfin_labeled_combined,
                                            var_str = var_str,
                                            time_str = time_str)

  df_pots_orig_combined[[thr_str]] <- NULL
  df_pots_orig_combined[[exc_str]] <- NULL

  return (list("pots"=df_pots_orig_combined, "storms"=dfin_labeled_combined))
}

peak_picking_idx <- function(dfin, lst_vars, thr_model_fml,
                             thr=.5, var_str='hs', thr_str='thr',
                             exc_str='exc', time_str='time',
                             decorrelation_time_scale=(2*24)){
  #' @export

  print('apply threshold model over all data')
  m_ald <- evgam(thr_model_fml, dfin, family="ald",
                 ald.args=list(tau=thr))
  q_tmp <- fitted(m_ald)$location

  dfin[[thr_str]] <- fitted(m_ald)$location
  dfin[[exc_str]] <- dfin[[var_str]] - dfin[[thr_str]]

  print('label storms')
  dfin_labeled <- label_storms_variable_thr(dfin, exc_str = exc_str)

  print('reduce dataset to storms')
  lst_vars <- c(lst_vars, exc_str, thr_str)
  dfin_labeled_red <- dfin_labeled[lst_vars]
  colnames(dfin_labeled_red) <- lst_vars

  print('combine and relabel storm peaks that are too close')
  # 1. round of peak finding
  df_pots_orig <- find_storm_peaks(dfin_labeled_red, var_str = var_str,
                                   time_str = time_str, exc_var = exc_str)

  # 2. check if storm peaks too close, if true combine
  storm_idx_list <- get_list_of_close_storms_idx(df_pots_orig, var_str, time_str,
                                                 decorrelation_time_scale)

  dfin_labeled_storms <- subset(dfin_labeled_red, dfin_labeled_red$storm_idx>0)
  dfin_labeled_combined <- combine_storms(dfin_labeled_storms, storm_idx_list)

  # 3. find peaks again for new combined storms
  print('find peaks for combined storms')
  df_pots_orig_combined <- find_storm_peaks(dfin_labeled_combined,
                                            var_str = var_str,
                                            time_str = time_str)

  df_pots_orig_combined[[thr_str]] <- NULL
  df_pots_orig_combined[[exc_str]] <- NULL

  return (list("pots"=df_pots_orig_combined, "storms"=dfin_labeled_combined))
}
