Hmax_distr_from_traj_forristall_2008 <- function(A_hs, tdelta=3600, llim=0,
                                                 ulim=50, integr_step=.1,
                                                 tm2_const=10, multiplicator=8,
                                                 storm_base_pct=0,
                                                 verbose=FALSE){
  #' @param A_hs Array of hs
  #' @param tdelta Time-step in seconds
  #' @param llim Low limit of integration
  #' @param ulim Upper limit of integration
  #' @param inter_step Integration step
  #' @param tm2_const Constant tm2 value
  #' @param multiplicator Parameter for storm duration
  #' @param storm_base_pct Storm base percentage of HsPeak e.g. 0.7
  #' @param verbose If wanted set TRUE
  #'
  #' @return A_Hmax array of max individual waves during storm
  #'
  #' @examples
  #' A_Hmax <- Hmax_distr_from_traj_forristall_2008(A_rl)
  #'
  #' @export

  A_Hmax <- array(1, length(A_hs))*NA
  for (i in 1:length(A_hs)){
    trajs <- triangular_storm_traj(A_hs[i], tm2_const, multiplicator, tdelta, storm_base_pct)
    hs_storm <- trajs[[1]]
    tm2_storm <- trajs[[2]]
    A_Hmax[i] <- storm_trajectory_Hmax(hs_storm, tm2_storm,tdelta)[1] # 1=mode, 2=mean
  }
  return(A_Hmax)
}

storm_trajectory_Hmax <- function(hs_traj, tm02_traj, tdelta=3600, llim=0,
                                  ulim=50, integr_step=.1, dist="forristall",
                                  verbose=FALSE){

  #' @param hs_traj Array of hs resembling the storm trajectory in hs-space
  #' @param tm02_traj Array tm02 hs resembling the storm trajectory in tm02-space
  #' @param tdelta Time-step in seconds
  #' @param llim Low limit of integration
  #' @param ulim Upper limit of integration
  #' @param inter_step Integration step
  #' @param verbose If wanted set TRUE
  #'
  #' @return list Returns expected and most likely Hmax
  #'
  #' @examples
  #' A_Hmax <- storm_trajectory_Hmax(hs_storm, tm2_storm, tdelta)[1]
  #'
  #' @export

  x_lst <- seq(0,ulim,integr_step)
  P_i_lst <- array(0, length(hs_traj))*NA
  P_s_lst <- array(0, length(x_lst))*NA

  for (j in 1:length(x_lst)) {
    for (i in 1:length(hs_traj)) {
      # Compute P based on Hmax: x
      Hs <- hs_traj[i]

      # P(H>=x) one draw
      if (dist == "rayleigh"){
        P_H <- rayleigh_Hs(Hs, x_lst[j])
      } else if (dist == "forristall") {
        P_H <- forristall_Hs(Hs, x_lst[j])
      } else if (dist == "prevosto") {
        P_H <- prevosto_Hs_2nd(Hs, x_lst[j])
      } else {
        print("dist is not valid!")
      }
      # Compute N based on Tm02
      N_i <- 1/tm02_traj[i]*tdelta

      # get probability for all draws (N) from P_H
      P_i <- 1-(1-P_H)**N_i
      P_i_lst[i] <- 1-P_i
    }
    P_s <- 1-prod(P_i_lst)
    P_s_lst[j] <- P_s
  }

  deriv <- diff(1-P_s_lst)
  x_lst_diff <- x_lst[2:length(x_lst)]

  # find max of deriv
  mode_of_deriv <- x_lst_diff[deriv==max(deriv)]
  E_of_deriv <- sum(deriv*x_lst_diff)
  if (verbose==TRUE){
    print(c('Most likely Hmax:', mode_of_deriv))
    print(c('Expected Hmax:', E_of_deriv))
  }
  return(list("Hmax_mode"=mode_of_deriv,"Hmax_E"=E_of_deriv,"deriv"=deriv,"x_lst_diff"=x_lst_diff))
}

forristall_Hs <- function(Hs,H,a=0.681,b=2.126){

  #' @param Hs Significant wave height
  #' @param H Individual wave height
  #'
  #' @return P Probability of H/Hs
  #'
  #' @export

  P <- exp(-(H/(a*Hs))**b)
  return(P)
}

prevosto_Hs_2nd <- function(Hs,H,a=2.13,b=8.42){
  # see e.g. Feld et al 2015 with
  # distribution parameter values from Prevosto et al. 2000

  #' @param Hs Significant wave height
  #' @param H Individual wave height
  #'
  #' @return P Probability of H/Hs
  #'
  #' @export

  P <- exp(-(1/b)*(H/(Hs/4))**a)
  return(P)
}

rayleigh_Hs <- function(Hs,H){

  #' @param Hs Significant wave height
  #' @param H Individual wave height
  #'
  #' @return P Probability of H/Hs
  #'
  #' @export

  P <- exp(-(2*H**2/Hs**2))
  return(P)
}

Hmax_rayleigh_Hs <- function(P,Hs) {

  #' @param Hs Significant wave height
  #' @param P Probability of H/Hs
  #'
  #' @return Hmax
  #'
  #' @export

  # e.g. Holthuijsen p.69
  # m0 <- Hs_to_m0(Hs)
  # Hmax <- sqrt(-8*m0*log(1-P))
  Hmax <- sqrt(log(1/P)/2)*Hs
  return(Hmax)
}

#sample_example <- function(n, Hs, k){
#  return(sqrt(log(1/runif(n))/2)*Hs)
#}

Hmax_forristall_Hs <- function(P,Hs,a=0.681,b=2.126) {

  #' @param Hs Significant wave height
  #' @param P Probability of H/Hs
  #'
  #' @return Hmax
  #'
  #' @export

  Hmax <- ((log(1/P))**b)*a
  return(Hmax)
}

triangular_storm_traj <- function(HsPeak, tm2_const=10, multiplicator=8,
                                  tdelta=3600, storm_base_pct=.7) {
  #' @param HsPeak Peak Hs
  #' @param tm2_const Constant mean period
  #' @param tdelta Time-step in seconds
  #' @param multiplicator Parameter for storm duration
  #' @param storm_base_pct Storm base percentage of HsPeak e.g. 0.7
  #'
  #' @return list Relevant Hs and Tm
  #'
  #' @export

  # triangular storm model
  # model conditional on HsPeak with duration D = multiplicator*HsPeak with unit hours
  # multiplicator is 8 according literature
  # Tm2 is constant with 10s

  D_storm <- closest_even(HsPeak*multiplicator)
  sinalpha <- HsPeak/(D_storm/2)
  Hs_storm_tmp <- array(1, D_storm/2)*NA
  tm2_storm <- array(1, (D_storm-1))*10
  for (j in seq(D_storm/2)){
    Hs_storm_tmp[j] <- sinalpha*j # in hours
  }
  Hs_storm <- c(Hs_storm_tmp, rev(Hs_storm_tmp)[2:(length(Hs_storm_tmp))])
  relevant_Hs <- Hs_storm[Hs_storm>(storm_base_pct*HsPeak)]
  relevant_tm2 <- tm2_storm[Hs_storm>(storm_base_pct*HsPeak)]
  return(list(relevant_Hs, relevant_tm2))
}

power_storm_traj <- function(Peak, B=NULL, lam=NULL, xPeak=NULL) {
  #'
  #' @return list Relevant Hs
  #'
  #' @export

  # lam = 1 <- triangular
  # lam > 1 round
  # lam < 1 approaching delta fct

  A <- Peak

  if (is.null(B)){
    B <- 150*A**(-.25)
    #B <- 50*A**(-.25)
  }

  if (is.null(lam)){
    lam <- 1
  }

  Bhalf <- as.integer(B/2)
  t <- seq(-Bhalf,Bhalf)

  # Hs <- A*(1-(2*abs(t-t0)/B)**lam)
  y <- A*(1-(2*abs(t)/B)**lam)
  t <- t + Bhalf
  if (is.null(xPeak)){
    x <- t
  } else {
    x <- t - Bhalf + xPeak*24
    x2 <- t
  }
  y <- y[x>0]
  t <- t[x>0]
  x <- x[x>0]

  # returns hourly ts of one storm
  return(list('y'=y,'x'=x, 't'=t))
}

power_storm_trajs <- function(Peaks, B=NULL, lam=1, xPeaks=NULL, slen=NULL){
  # establish final matrix containing all single storm ts
  dimensions <- c(length(Peaks), slen)
  dummy_data <- runif(prod(dimensions))*NA
  ss_array <- array(dummy_data, dim = dimensions)
  for (i in 1:length(Peaks)) {
    ss <- power_storm_traj(Peaks[i], B=B, lam=lam, xPeak=xPeaks[i])
    ss_array[i,ss$x[1]:ss$x[length(ss$x)]] <- ss$y
  }
  return(ss_array)
}

power_storm_traj_composite_max <- function(Peaks, B=NULL, lam=1, xPeaks=NULL, slen=NULL){
  #' @export
  #'
  dimensions_tmp <- c(length(Peaks), slen)
  dummy_data_tmp <- runif(prod(dimensions_tmp))*NA
  ss_array_tmp <- array(dummy_data_tmp, dim = dimensions_tmp)

  print(dim(ss_array_tmp))

  ss1 <- power_storm_traj(Peaks[1], B=B, lam=lam, xPeak=xPeaks[1])

  for (i in 2:length(Peaks)) {
    #print(c('nr peaks:',i))
    ss <- power_storm_traj(Peaks[i], B=B, lam=lam, xPeak=xPeaks[i])
    ss_array_tmp[i,(ss$x[1]:(ss$x[length(ss$x)]))] <- ss$y

    #ss_array_tmp[1,] <- apply(ss_array_tmp, 2, max, na.rm=TRUE)
  }
  return(ss_array_tmp)
}

compute_steepness <- function(Hm0, Tm01){
  #' @param Hm0 Significant wave height
  #' @param Tm01 Mean period
  #'
  #' @return s Steepness
  #'
  #' @export

  # computes wave steepness from Tm02/02/p
  # e.g. Holthuijsen p.58
  # instead of Tm01 it could be Tm02

  g <- 9.81
  s <- Hm0 / (g*Tm01^{2}/(2*pi))
  return(s)
}

compute_T_from_steepness <- function(hs,s,g=9.81){
  #' @param hs Significant wave height
  #' @param s Steepness
  #'
  #' @return T Period used for steepness
  #'
  #' @export

  # Computes mean T from Hs and wave steepness
  # Could be Tm01/Tm02 depending on the definition used for steepness

  T <- sqrt((hs/g*2*pi/s))
  return(T)
}

find_idxs_closest_storm_peaks_hs_tm <- function(simPeakHs, simPeakTm,
                                          histPeaksHs, histPeaksTm,
                                          sidx=1, eidx=10){
  #' @param simPeakHs simulated peak of Hs
  #' @param simPeakTm simulated peak of Tm
  #' @param histPeaksHs historic peaks of Hs to match
  #' @param histPeaksTm historic peaks of Tm to match
  #' @param sidx start idx
  #' @param eidx end idx
  #'
  #' @return Array of closest storm peak idx
  #'
  #' @export

  distsHs <- abs(histPeaksHs-simPeakHs)
  distsTm <- abs(histPeaksTm-simPeakTm)
  dists <- sqrt(distsHs**2+distsTm**2)
  return(order(dists)[sidx:eidx])
}

find_idxs_closest_storm_peaks_hs <- function(simPeakHs, histPeaksHs,
                                             sidx=1, eidx=10){
  #' @param simPeakHs simulated peak of Hs
  #' @param histPeaksHs historic peaks of Hs to match
  #' @param sidx start idx
  #' @param eidx end idx
  #'
  #' @return Array of closest storm peak idx
  #'
  #' @export

  dists <- abs(histPeaksHs-simPeakHs)
  return(order(dists)[sidx:eidx])
}

find_idx_closest_storm_peaks <- function(simPeaksMV, histPeaksMV, varlst=NULL,
                                         sidx=1, eidx=10){
  #' @param simPeaksMV multivariate simulated peak of Hs
  #' @param histPeaksMV multivariate historic peaks of Hs to match
  #' @param varlst list of matching variables
  #' @param sidx start idx
  #' @param eidx end idx
  #'
  #' @return Array of closest storm peak idx
  #'
  #' @export

  dists <- NULL
  if (is.null(varlst)){
    varlst <- names(simPeaksMV)
  }

  for (n in 1:length(varlst)) {
    dists[[varlst[n]]] <- abs(histPeaksMV[[varlst[n]]]-simPeakHs[[varlst[n]]])
  }

  # multidim Pythagoras
  sumdist <- sqrt(unlist(dists))/length(unlist(dists))

  distsHs <- abs(histPeaksHs-simPeakHs)
  distsTm <- abs(histPeaksTm-simPeakTm)
  dists <- sqrt(distsHs**2+distsTm**2)
  return(order(dists)[sidx:eidx])
}

get_constant_scaling <- function(simPeak, histPeaks){
  #' @param simPeak simulated peak
  #' @param histPeaks historic peaks to match
  #'
  #' @return scaling scaling factor
  #'
  #' @export

  scaling <- simPeak/histPeaks
  return(scaling)
}

scale_storm <- function(stormTrajHs, stormTrajTm, scaling){
  #' @param simPeak simulated peak
  #' @param histPeaks historic peaks to match
  #'
  #' @return scaling scaling factor
  #'
  #' @export

  # this function scales the simulated storm according to the closest observed one

  stormTrajHs_scaled <- stormTrajHs * scaling
  stormTrajHsSteep <- compute_steepness(stormTrajHs, stormTrajTm)
  stormTrajTm_scaled <- compute_tm1_or_tm2_from_steepness(stormTrajHs_scaled,
                                                          stormTrajHsSteep,
                                                          g=9.81)
  return(list(stormTrajHs_scaled, stormTrajTm_scaled))
}

find_storm_match <- function(simPeakHs, simPeakTm=0,
                             sidx=1, eidx=10,
                             histPeaks, var_str,
                             storm_df){
  #' @param simPeak simulated peak of Hs
  #' @param simPeak simulated peak of Tm
  #' @param sidx start idx
  #' @param eidx end idx
  #' @param hist_df historic dataframe
  #' @param storm_thr threshold for storm definition for integration
  #' @param var_str string of variable of interest
  #'
  #' @return stormTraj_scaled scaled storm trajectory
  #'
  #' @export

  # This function scales the simulated storm according to the closest observed one
  # get closest peak indices
  if (simPeakTm==0){
    storm_peak_idxs <- find_idxs_closest_storm_peaks_hs(simPeakHs, histPeaks$hs,
                                                        sidx=sidx, eidx=eidx)
  } else {
    storm_peak_idxs <- find_idxs_closest_storm_peaks(simPeakHs, simPeakTm,
                                                     histPeaks$hs, histPeaks$tm2,
                                                     sidx=sidx, eidx=eidx)
  }
  # draw storm at random
  random_storm_idx <- runif_func(1, min=1, max=length(storm_peak_idxs))
  # get idx for matching random peak
  idx_peak <- storm_peak_idxs[random_storm_idx]
  # get idx for matching storm
  idx_storm <- histPeaks$storm_idx[idx_peak]
  # isolate storm trajectory
  stormTraj <- subset(storm_df, storm_df$storm_idx==idx_storm)
  # scale storm
  scaling <-  get_constant_scaling(simPeakHs, histPeaks$hs[idx_peak])
  stormTraj_scaled <- scale_storm(stormTraj$hs, stormTraj$tm2, scaling)
  return(stormTraj_scaled)
}

Hmax_distr_from_hist_matching_v0 <- function(simPeaksHs, simPeakTm=0,
                                          sidx=1, eidx=10,
                                          hist_df, storm_thr, var_str,
                                          tdelta=3600){
  #' @param simPeaksHs Array of simulated peaks of Hs
  #' @param simPeakTm Simulated peak of Tm
  #' @param sidx Start idx
  #' @param eidx End idx
  #' @param hist_df Historic dataframe
  #' @param storm_thr Threshold for storm definition for integration
  #' @param var_str String of variable of interest
  #' @param tdelta Time-step in seconds
  #'
  #' @return A_Hmax Array of Hmax'es
  #'
  #' @export

  # define storm peaks
  storm_df <- define_storms(hist_df, var_str, storm_thr)
  histPeaks <- find_storm_peaks(storm_df, "hs", "time", min_nr_days = 1)

  A_Hmax <- array(1, length(simPeaksHs))*NA
  for (p in 1:length(simPeaksHs)) {
    trajs <- find_storm_match(simPeaksHs[p], simPeakTm = simPeakTm,
                              histPeaks = histPeaks, var_str = var_str,
                              storm_df = storm_df)
    hs_storm <- trajs[[1]]
    tm2_storm <- trajs[[2]]
    A_Hmax[p] <- storm_trajectory_Hmax(hs_storm, tm2_storm, tdelta)[1] # 1=mode, 2=mean
  }
  return(A_Hmax)
}

Hmax_distr_from_hist_matching  <- function(simPeaksMV, histPeaksMV, histStormsMV,
                                           idxmatches, tdelta=3600, dist="forristall"){
  #' @param simPeaksMV Array of simulated peaks of Hs
  #' @param histPeaksMV Array of simulated peaks of Hs
  #' @param histStormsMV Array of simulated peaks of Hs
  #' @param idxmatches
  #' @param tdelta Time-step in seconds
  #'
  #' @return A_Hmax Array of Hmax'es
  #'
  #' @export

  # define storm peaks

  A_Hmax <- array(1, length(idxmatches))*NA
  for (i in 1:length(idxmatches)) {
    print(i)
    # isolate storm trajectory
    print("get storm_idx")
    storm_idx <- histPeaksMV$storm_idx[idxmatches[i]]
    print(storm_idx)
    print("get storm traj")
    stormTraj <- subset(histStormsMV, histStormsMV$storm_idx==storm_idx)
    # scale storm
    print("scale storm")
    scaling <-  get_constant_scaling(simPeaksMV$hs[i], histPeaksMV$hs[idxmatches[i]])
    stormTraj_scaled <- stormTraj$hs * scaling
    hs_storm <- stormTraj_scaled
    print("compute T")
    tm_storm <- compute_T_from_steepness(stormTraj_scaled, stormTraj$s_tm1)
    print(hs_storm)
    print(tm_storm)
    print("compute Hmax for storm")
    A_Hmax[i] <- storm_trajectory_Hmax(hs_storm, tm_storm, tdelta, dist=dist)[1] # 1=mode, 2=mean
  }
  return(A_Hmax)
}
