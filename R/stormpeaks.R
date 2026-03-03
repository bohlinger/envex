# imports
library(ppgam)
library(evgam)
library(extRemes)
#library(lubridate)
library(dplyr)
library(rlang)
library(MASS)

# applying Ben's version to data_sub_orig$y
# rounding function
round_to_nearest <- function(x,interval,offset) {
  rounded <- (round(x / interval) * interval) + offset
  wrapped <- rounded %% 360
  return(wrapped)
}

produce_storm_occurrences_rejection <- function(nr_of_events, RP,
                                      nr_of_years,
                                      model_nr_of_events,
                                      covarstr_lst,
                                      covar_mins, covar_maxes,
                                      dfin = dfin,
                                      condition = NULL,
                                      grid_interval = NULL)
{
  # simulate storm peaks from poisson
  print(c('nr_of_events:', nr_of_events))

  # probing covariate space on defined grid to cover everything efficiently
  # and hence retrieve the best lambda_max
  covar_grid_sizes <- NULL
  covar_grid <- NULL
  for (i in 1:length(covarstr_lst)){
    vargrid <- seq(covar_mins[i], covar_maxes[i], grid_interval[i])
    covar_grid[[covarstr_lst[i]]] <- vargrid
    covar_grid_sizes[[covarstr_lst[i]]] <- length(vargrid)
  }
  nr_of_grid_cells <- prod(unlist(covar_grid_sizes))
  mv_grid <- expand.grid(covar_grid)

  newdf <- as.data.frame(array(NA,dim=c((nr_of_grid_cells),(length(covarstr_lst)+1))))
  colnames(newdf) <- c('counts', covarstr_lst)

  for (i in 1:length(covarstr_lst)){
    newdf[[covarstr_lst[i]]] <-mv_grid[[covarstr_lst[i]]]
  }
  pred_lambdas <- predict(model_nr_of_events, newdata = newdf, type = 'response')

  # get max of predicted lambdas for rejection probability
  lambda_max <- max(pred_lambdas)
  print(c('sum(predcounts)',sum(pred_lambdas)))
  print(c('lambda_max',lambda_max))

  # establish proposal_N, number of random samples from the covariate space
  proposal_N <- rpois(1, nr_of_events * lambda_max * RP)
  print(c('proposal_N:',proposal_N))

  # draw random values for each covariate according to their EDF
  # later replaced by rejection sampling to be more general
  dfobs <- dfin[,covarstr_lst]

  samples = NULL
  if (length(covarstr_lst)==1){
    samples[[covarstr_lst[1]]] = sample(dfobs, proposal_N, replace = TRUE)
  }
  else if (length(covarstr_lst)>1){
    for (c in covarstr_lst){
      samples[[c]] = sample(dfobs[[c]], proposal_N, replace = TRUE)
    }
  }

  pred_data <- as.data.frame(samples)
  new_pred_lambda <- predict(model_nr_of_events, newdata = pred_data, type = 'response')

  counts_lst <- NULL
  keep_lst <- NULL
  for (i in 1:proposal_N){
    keep <- new_pred_lambda[i] / lambda_max > runif(1)
    if (keep==TRUE){
      counts_lst[[i]] <- 1
      keep_lst[[i]] <- TRUE
    } else{
      keep_lst[[i]] <- FALSE
    }
  }

  counts_lst <- unlist(counts_lst)
  keep_unlist <- unlist(keep_lst)

  if (length(covarstr_lst)==1){
    dfcounts <- NULL
    dfcounts[[covarstr_lst[1]]] <- pred_data[keep_unlist,]
    dfcounts <- as.data.frame(dfcounts)
  } else{
    dfcounts <- pred_data[keep_unlist,]
    dfcounts$counts <- counts_lst
  }

  # filter according to covariate constraints
  if (!is.null(condition)){
    print("filtering for condition:")
    print(c("  ", condition))
    condition_expr <- parse_expr(condition)
    dfcounts <- dfcounts %>% filter(!!condition_expr)
    nr_of_events_sim <- dim(dfcounts)[1]
  }

  print(c('nr_of_sims:',dim(dfcounts)[1]))

  return(dfcounts)
}

produce_storm_occurrences <- function(nr_of_events, RP,
                                      nr_of_years,
                                      model_nr_of_events,
                                      covarstr_lst,
                                      covar_mins, covar_maxes,
                                      condition=NULL){
  #' Produce storm occurrences by producing the appropriate number
  #' of coinciding covariates, i.e. direction and day of year.
  #'
  #' @param nr_of_events nr of events in dataset (integer)
  #' @param RP return period (RP) (integer)
  #' @param nr_of_years nr of years the dataset covers (integer)
  #' @param model_nr_of_events ppgam model object for occurrences given the covariates
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

  newdf <- as.data.frame(array(NA,dim=c(nr_of_events_sim,(length(covarstr_lst)+1))))
  colnames(newdf) <- c('counts', covarstr_lst)

  for (i in 1:length(covarstr_lst)){
    covar_samples <- runif(nr_of_events_sim, min = covar_mins[i], max = covar_maxes[i])
    newdf[[covarstr_lst[i]]] <- covar_samples
  }

  predcounts <- predict(model_nr_of_events, newdata = newdf, type = 'response')

  simcounts <- rpois(length(predcounts), lambda=predcounts)
  dftmp <- newdf
  dftmp[['counts']] <- simcounts

  print(c("simulate total counts of", sum(simcounts)))

  print(c("max predicted counts per covar value",max(predcounts)))
  print(c("max simulated counts per covar value",max(simcounts)))

  dfcounts_tmp <- subset(dftmp, counts>0)
  dfcounts_tmp_1 <- subset(dfcounts_tmp, counts==1)
  dfcounts_tmp_l1 <- subset(dfcounts_tmp, counts>1)

  if (dim(dfcounts_tmp_l1)[1]>0){
    dfcounts_1 <- unfold_counts(dfcounts_tmp_l1, covarstr_lst)
    dfcounts <- rbind(dfcounts_1,dfcounts_tmp_1)
  } else {
    dfcounts <- dfcounts_tmp_1[,covarstr_lst]
  }
  rm(dftmp)

  # filter according to covariate constraints
  if (!is.null(condition)){
    print("filtering for condition:")
    print(c("  ", condition))
    condition_expr <- parse_expr(condition)
    dfcounts <- dfcounts %>% filter(!!condition_expr)
    nr_of_events_sim <- dim(dfcounts)[1]
  }

  print(c('nr_of_events_sim:', dim(dfcounts)[1]))
  return(dfcounts)
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

get_Hs_bins_and_fit_ln <- function(df_all, Hs_str, T_str, sbinval=NULL, ebinval=NULL, binsize=1){
  #' @param df_all dataframe to use with all info included
  #' @param Hs_str string of Hs variable
  #' @param Tp_str string of T variable
  #'
  #' @return vectors of ln_means and ln_stds of T for Hs
  #'
  #' @export
  # bin hs in bins of 1m and perform a log-norm fit to all periods (e.g. Tp, Tm01, Tm02, ...) in one bin
  if (is.null(sbinval)){
    sbinval=1
  }
  if (is.null(ebinval)){
    ebinval=as.integer(max(df_all[[Hs_str]]))
  }
  seqtmp <- seq(sbinval,ebinval,binsize)
  ln_means <- array(length(seqtmp))*NA
  ln_stds <- array(length(seqtmp))*NA
  for (i in 1:length(seqtmp)){
    print(seqtmp[i])
    df_sub <- subset(df_all, df_all[[Hs_str]]>=(seqtmp[i]-binsize) & df_all[[Hs_str]]<(seqtmp[i]))
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

  error <- sum((Y - (a0+a1*X**a2))^2)
  return(error)
}

calc_Tp_RV <- function(ln_means, Hs_RV, sbinval=NULL, ebinval=NULL, binsize=1){
  #' @param ln_means log-normal fits of T variable for Hs bins
  #' @param Hs_RV Return value of Hs
  #'
  #' @return T_joint_RV Return value of T variable
  #'
  #' @export

  if (is.null(sbinval)){
    sbinval=1
  }
  if (is.null(ebinval)){
    ebinval=as.integer(max(df_all[[Hs_str]]))
  }
  seqtmp <- seq(sbinval,ebinval,binsize)
  X <- seqtmp-(binsize/2)
  Y <- ln_means

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

get_probability <- function(dfin, m_ald, m_gpd, var_str = 'hs', nsim=1){
  #' @export

  ald_thr <- predict(m_ald, newdata = dfin, type='response')$location

  dfin[['threshold']] <- ald_thr
  dfin[['exc']] <- dfin[[var_str]]-ald_thr

  dfgpd <- subset(dfin,exc>0)
  dfecdf <- subset(dfin,exc<=0)

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
  return (probs)
}


fit_marginal_models_thr <- function(dfin, thr, model_fml, list_var=NULL, m_params=NULL, knots=NULL){
  #' @export
  #'
  if (is.null(list_var)){
    list_var <- names(model_fml)
  }
  margs <- NULL
  for (n in list_var){
    print(c("fit threshold model for",n))
    # bootstrap extremal threshold within bounds
    if (length(thr[[n]])>1){
      t <- round(runif(1, min = thr[[n]][1], max = thr[[n]][2]),2)
    } else {t <- thr[[n]]}

    if (is.null(m_params)){
      margs[[n]] <- evgam(model_fml[[n]], dfin, family="ald", ald.args=list(tau=t), knots=knots)
    } else {
      margs[[n]] <- evgam(model_fml[[n]], dfin, family="ald",
                          ald.args=list(tau=t), sp=m_params[[n]]$sp,
                          knots=knots)
    }
  }
  return (margs)
}


fit_marginal_models_gpd <- function(dfin, model_fml, list_var=NULL, family='gpd2',
                                    m_params=NULL, knots=NULL, trace=0){
  #' @export
  #'
  margs <- NULL
  if (is.null(list_var)){
    list_var <- names(model_fml)
  }
  for (n in list_var){
    print(c("fit gpd model for",n))
    if (is.null(m_params)){
      margs[[n]] <- evgam(model_fml[[n]], dfin[[n]], family=family,
                          trace=trace, knots=knots)
    } else {
      margs[[n]] <- evgam(model_fml[[n]], dfin[[n]], family=family,
                          sp=m_params[[n]]$sp, trace=trace,
                          knots=knots)
    }
  }
  return (margs)
}

fit_marginal_models_occ <- function(dfin, model_fml, nr_of_years,
                                    list_var=NULL, nquad=20, m_params=NULL,
                                    knots=NULL, nodes=NULL){
  #' @export
  #'
  margs <- NULL
  if (is.null(list_var)){
    list_var <- names(model_fml)
  }

  for (n in list_var){
    print(c("fit occurrence model for",n))

    if (is.null(m_params)){
      margs[[n]] <- ppgam(model_fml[[n]],
                          data = dfin[[n]],
                          weights = nr_of_years*dim(dfin[[1]])[1],
                          nquad = nquad,
                          knots = knots,
                          nodes = nodes)
    } else {
      margs[[n]] <- ppgam(model_fml[[n]],
                          data = dfin[[n]],
                          weights = nr_of_years*dim(dfin[[1]])[1],
                          sp = m_params[[n]]$sp,
                          nquad = nquad,
                          knots = knots,
                          nodes = nodes)
    }
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
                         decorrelation_time_scale=2, idx=FALSE){
  #' @export

  if (is.null(thr_model_fml)){
    print('apply constant threshold to all data')
    dfin[[thr_str]] <- array(1, length(dfin[[var_str]]))*thr
    dfin[[exc_str]] <- dfin[[var_str]] - dfin[[thr_str]]
  } else {
    print('apply threshold model to all data')
    m_ald <- evgam(thr_model_fml, dfin, family="ald",
                   ald.args=list(tau=thr))
    q_tmp <- fitted(m_ald)$location

    dfin[[thr_str]] <- fitted(m_ald)$location
    dfin[[exc_str]] <- dfin[[var_str]] - dfin[[thr_str]]
  }

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
  if (isTRUE(idx)){
    storm_idx_list <- get_list_of_close_storms_idx(df_pots_orig, var_str, time_str,
                                                   decorrelation_time_scale)
  } else{
    storm_idx_list <- get_list_of_close_storms(df_pots_orig, var_str, time_str,
                                               decorrelation_time_scale)
  }

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
