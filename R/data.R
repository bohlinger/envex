#' Data to showcase modelling of sea state extremes depending smoothly on covariates
#'
#' Contains various bulk sea state variables plus time from the Ekofisk oil and gas field
#'
#' @format A data frame with 22 variables and 385704 rows:
#'  \describe{
#'    \item{time}{"character"}
#'    \item{Pdir}{peak direction}
#'    \item{fpI}{peak frequency (interpolated)}
#'    \item{hs}{significant wave height combined (m0 based)}
#'    \item{hs_sea}{significant wave height (wind)}
#'    \item{hs_swell}{significant wave height (swell)}
#'    \item{thq}{mean wave direction combined}
#'    \item{thq_sea}{mean wave direction (wind)}
#'    \item{thq_swell}{mean wave direction (swell)}
#'    \item{tm1}{mean period (m1)}
#'    \item{tm2}{mean period (m2)}
#'    \item{tmp}{mean period (m-1)}
#'    \item{tp}{peak period combined interpolated}
#'    \item{tp_sea}{peak period (wind)}
#'    \item{tp_swell}{peak period (swell)}
#'    \item{yy}{year}
#'    \item{mm}{month}
#'    \item{doy}{day of year}
#'    \item{dt}{dt <- as.POSIXct(ekofisk$time, format="%Y-%m-%d %H:%M:%S", tz="UTC")}
#'    \item{nt}{as.numeric(dt)}
#'    \item{U10}{wind speed at 10m height}
#'    \item{zeta}{storm surge, tides removed}
#'    }
#' @source {NORA3 wave hindcast data available under:
#'          https://thredds.met.no/thredds/catalog/windsurfer/mywavewam3km_files/catalog.html.
#'          https://thredds.met.no/thredds/projects/stormrisk.html
#'          Time series nearest- neighbour interpolated to location of Ekofisk oil and gas field.}
#' @examples
#' data(ekofisk_wave_surge)
"ekofisk_wave_surge"
