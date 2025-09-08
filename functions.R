lc.uvot <- function(file.path, flux.units){
  #' @title Reads in uvotsource output files and makes a light curve data frame
  #' @description This function will read in all of the output files produced by
  #' the uvotsource HEASoft tool that are stored in the given file path.
  #' @param file.path System file path to the directory that contains the output
  #' files from uvotsource.
  #' @param flux.units Light curve flux units. Options are: rate (counts/second)
  #' , flux (ergs/s/cm^2/A), and mJy (milliJansky).
  #' @returns A data frame with columns for observation ID, UVOT telescope
  #' filter, time (MET), time error (MET), flux, and flux error.
  #' @export

  working.path <- getwd()

  ### Determine what flux and error columns to grab based on user input
  if (flux.units == "rate"){
    val <- "CORR_RATE"
    err <- "CORR_RATE_ERR"
  } else if (flux.units == "flux"){
    val <- "FLUX_AA"
    err <- "FLUX_AA_ERR"
  } else if (flux.units == "mJy"){
    val <- "FLUX_HZ"
    err <- "FLUX_HZ_ERR"
  } else {
    cat("***flux units not recognized!\nplease enter 'rate', 'flux', or 'Jy'\n")
    break
  }

  ### Read in all of the file names, splitting in OBSID and filter information
  setwd(file.path)
  all.files <- list.files(pattern = ".fits")
  df.files <- data.frame(FILENAME = rep(NA, length.out = length(all.files)))
  df.files$OBSID <- NA
  df.files$FILTER <- NA
  for (i in 1:length(all.files)){
    df.files$FILENAME[i] <- all.files[i]
    df.files$OBSID[i] <- strsplit(all.files[i], "_")[[1]][1]
    df.files$FILTER[i] <- strsplit(strsplit(all.files[i], "_")[[1]][2], ".fits")[[1]][1]
  }
  rm(all.files)

  ### Create the full UVOT light curve
  obs.ids <- unique(df.files$OBSID)
  lc <- data.frame(OBSID = NA,
                   FILTER = NA,
                   TIME = NA_real_, TIME.ERR = NA_real_,
                   RATE = NA_real_, RATE.ERR = NA_real_)

  # for (k in 1:length(obs.ids)){
  for (k in 1:10){
    cat(paste("Working on observation:\t", k, " / ", length(obs.ids)), "\r")
    tmp.files <- df.files[which(df.files$OBSID == obs.ids[k]),]

    if (length(tmp.files$OBSID) == 0){
      next
    }

    for (j in 1:length(tmp.files$FILENAME)){
      file.tmp <- FITSio::readFITS(tmp.files$FILENAME[j])
      idx.val <- which(file.tmp$colNames == val)
      idx.err <- which(file.tmp$colNames == err)
      lc <- rbind(lc, data.frame(OBSID = tmp.files$OBSID[j],
                                 FILTER = tmp.files$FILTER[j],
                                 TIME = file.tmp$col[[1]], TIME.ERR = file.tmp$col[[4]]-file.tmp$col[[3]],
                                 RATE = file.tmp$col[[idx.val]], RATE.ERR = file.tmp$col[[idx.err]]))
    }
    rm(list = ls(pattern = "tmp"))
  }
  lc <- lc[-1,]
  lc <- lc[order(lc$TIME),]

  setwd(working.path)
  return(lc)
}


remove.dropouts <- function(light.curve, minimum.separation = NULL, tolerance = 0){
  #' @title Remove drop-out points from UVOT light curves
  #' @description Takes an input UVOT light curve and removes the so-called
  #' "drop-out" points according to the method of Edelson et al. (2015).
  #' @param light.curve Light curve data frame for 1 of the UVOT filters (should
  #' have column names: TIME, TIME.ERR, RATE, RATE.ERR).
  #' @param minimum.separation Minimum time separation required between two
  #' consecutive light curve points.
  #' @param tolerance Threshold level (in terms of sigma) to compare with
  #' maximum positive deviation.
  #' @return Data frame of the light curve without drop-outs
  #' @export

  if (is.null(minimum.separation) == T){
    cat("***please supply a minimum separation between data points!\n")
    break
  }

  light.curve$FLAG <- NA
  light.curve$DF <- 0

  for (i in 2:(length(light.curve$TIME)-1)){
    t.sep <- light.curve$TIME[i+1]-light.curve$TIME[i-1]
    if (t.sep <= minimum.separation){
      light.curve$FLAG[i] <- "G"
      light.curve$DF[i] <- (light.curve$RATE[i] - 0.5*(light.curve$RATE[i+1] + light.curve$RATE[i-1]))/light.curve$RATE.ERR[i]
    } else if (t.sep > minimum.separation){
      light.curve$FLAG[i] <- "B"
    }
  }
  light.curve$FLAG[which(is.na(light.curve$FLAG))] <- "G"
  light.curve$DF[which(light.curve$FLAG == "B")] <- NA
  cut.val <- max(sd(na.omit(light.curve$DF[which(light.curve$FLAG == "G")]))*tolerance, max(na.omit(light.curve$DF)))
  bad.idx <- which(light.curve$DF < -cut.val)

  if (length(bad.idx) > 0){
    light.curve <- light.curve[-bad.idx,]
    while(length(bad.idx) > 0){
      for (i in 2:(length(light.curve$TIME)-1)){
        t.sep <- light.curve$TIME[i+1]-light.curve$TIME[i-1]
        if (t.sep <= minimum.separation){
          light.curve$FLAG[i] <- "G"
          light.curve$DF[i] <- (light.curve$RATE[i] - 0.5*(light.curve$RATE[i+1] + light.curve$RATE[i-1]))/light.curve$RATE.ERR[i]
        } else if (t.sep > minimum.separation){
          light.curve$FLAG[i] <- "B"
        }
      }
      cut.val <- max(sd(na.omit(light.curve$DF[which(light.curve$FLAG == "G")]))*tolerance, max(na.omit(light.curve$DF)))
      bad.idx <- which(light.curve$DF < -cut.val)
      if (length(bad.idx) > 0){
        light.curve <- light.curve[-bad.idx,]
      } else if (length(bad.idx) == 0){
        break
      }
    }
  }

  return(light.curve[,c("TIME", "TIME.ERR", "RATE", "RATE.ERR")])
}


detrend.SavitzkyGolay <- function(light.curve, filter.width = NULL, polynomial.order = 1, output = NULL){
  #' @title De-trend a light curve using a Savitzky-Golay filter
  #' @description Using a Savitzky-Golay filter from the 'signal' package with
  #' the user-defined filter width and polynomial order, de-trend an input light
  #' curve.
  #' @param light.curve A light curve data frame (cols: TIME, TIME.ERR, RATE,
  #' RATE.ERR).
  #' @param filter.width The number of consecutive points over which to smooth.
  #' @param polynomial.order The order of polynomial to fit over the points.
  #' @param output Either "SG" to output the trend line or "detrended" to output
  #' the original RATE - the trend line rate.
  #' @export

  SG.rate <- signal::sgolayfilt(light.curve$RATE, p = polynomial.order, n = filter.width)

  if (is.null(output) == TRUE){
    cat("***you need to specify an output format as either 'SG' or 'detrended'!", "\n")
    stopifnot(TRUE)
  } else if (output == "SG"){
    df <- data.frame(TIME = light.curve$TIME, TIME.ERR = light.curve$TIME.ERR,
                     RATE = SG.rate)
  } else if (output == "detrended"){
    detrended.rate <- light.curve$RATE - SG.rate
    df <- data.frame(TIME = light.curve$TIME, TIME.ERR = light.curve$TIME.ERR,
                     RATE = detrended.rate, RATE.ERR = light.curve$RATE.ERR)
  }

  return(df)
}


calc.sf <- function(light.curve, min.counts = 10, var.signal = 1, var.noise = 0){
  #' @title Compute the structure function of a light curve
  #' @description Compute the structure function of a light curve using the
  #' method of Collier and Peterson (2001) and Gallo et al. (2018).
  #' @param light.curve light curve data with columns of TIME, RATE, and RATE.ERR (time must be zeroed at first index!)
  #' @return Data frame of the structure function
  #' @export

  dtau <- round(median(diff(light.curve$TIME)), digits = 2)
  tau.bins <- seq(from = 0, to = max(light.curve$TIME), by = dtau)
  sq.diff <- matrix(NA_real_, nrow = length(light.curve$TIME), ncol = length(light.curve$TIME))
  delta.t <- matrix(NA_real_, nrow = length(light.curve$TIME), ncol = length(light.curve$TIME))
  for (i in 1:length(light.curve$TIME)){
    for (j in 1:length(light.curve$TIME)){
      sq.diff[i,j] <- (light.curve$RATE[j] - light.curve$RATE[i])^2
      delta.t[i,j] <- light.curve$TIME[j] - light.curve$TIME[i]
    }
  }

  sf <- data.frame(tau = rep(NA_real_, length.out = length(tau.bins)-1),
                   min = rep(NA_real_, length.out = length(tau.bins)-1),
                   max = rep(NA_real_, length.out = length(tau.bins)-1),
                   val = rep(NA_real_, length.out = length(tau.bins)-1),
                   err = rep(NA_real_, length.out = length(tau.bins)-1),
                   num = rep(NA_real_, length.out = length(tau.bins)-1))

  for (i in 1:(length(tau.bins)-1)){
    bin.idx <- which(delta.t > tau.bins[i] & delta.t <= tau.bins[i+1])
    sf$num[i] <- length(bin.idx)
    sf$min[i] <- tau.bins[i]
    sf$max[i] <- tau.bins[i+1]
    sf$tau[i] <- 10^((log10(tau.bins[i])+log10(tau.bins[i+1]))/2)
    sf$val[i] <- (sum(sq.diff[bin.idx])/sf$num[i] - var.noise) / var.signal
    sf$err[i] <- sqrt(sum((sq.diff[bin.idx]-sf$val[i])^2)/sf$num[i]) / sqrt(sf$num[i]/2) / var.signal
  }

  sf <- sf[which(sf$num >= min.counts),]
  return(sf)
}


calc.iccf <- function(lightcurve.1, lightcurve.2, delta.tau, max.lag = NA){
  #' @title Compute the interpolation cross-correlation function
  #' @description Takes the input data frames of two light curves and computes the interpolation cross-correlation function (ICCF) as in Gaskell & Peterson (1987).
  #' @param lightcurve.1 data frame of first light curve
  #' @param lightcurve.2 data frame of second light curve
  #' @param delta.tau time step to sample lags over
  #' @param max.lag maximum lag to calculate
  #' @return Data frame of the ICCF
  #' @export

  ### Extract the light curve time and count rates
  x.1 <- lightcurve.1$TIME - lightcurve.1$TIME[1]
  y.1 <- lightcurve.1$RATE
  x.2 <- lightcurve.2$TIME - lightcurve.1$TIME[1]
  y.2 <- lightcurve.2$RATE

  ### Set up the lag values based on input delta.tau and max.lag
  if (is.na(max.lag) == F){
    if(max.lag%%delta.tau > 1){
      max.lag <- round(max.lag/delta.tau)*delta.tau
    }
  } else if (is.na(max.lag) == T){
    max.lag <- round(max(x.1)/delta.tau)*delta.tau
  }
  tau.vals <- seq(from = -max.lag, to = max.lag, by = delta.tau)

  ### Set up the dual ICCF arrays and compute
  iccf.vals.12 <- rep(0, length.out = length(tau.vals))
  iccf.vals.21 <- rep(0, length.out = length(tau.vals))
  iccf.vals <- rep(0, length.out = length(tau.vals))

  for (i in 1:length(tau.vals)){
    interp.y.2 <- approx(x.2, y.2, x.1+tau.vals[i])$y
    bad.idx.12 <- which(is.na(interp.y.2))
    if (length(bad.idx.12) > 0){
      iccf.vals.12[i] <- sum((y.1[-bad.idx.12]-mean(y.1[-bad.idx.12]))*(interp.y.2[-bad.idx.12]-mean(interp.y.2[-bad.idx.12])))/(sd(y.1[-bad.idx.12]-mean(y.1[-bad.idx.12]))*sd(interp.y.2[-bad.idx.12]-mean(interp.y.2[-bad.idx.12])))/(length(y.1[-bad.idx.12])-1)
    } else if (length(bad.idx.12) == 0){
      iccf.vals.12[i] <- sum((y.1-mean(y.1))*(interp.y.2-mean(interp.y.2)))/(sd(y.1-mean(y.1))*sd(interp.y.2-mean(interp.y.2)))/(length(y.1)-1)
    }

    interp.y.1 <- approx(x.1, y.1, x.2-tau.vals[i])$y
    bad.idx.21 <- which(is.na(interp.y.1))
    if (length(bad.idx.21) > 0){
      iccf.vals.21[i] <- sum((y.2[-bad.idx.21]-mean(y.2[-bad.idx.21]))*(interp.y.1[-bad.idx.21]-mean(interp.y.1[-bad.idx.21])))/(sd(y.2[-bad.idx.21]-mean(y.2[-bad.idx.21]))*sd(interp.y.1[-bad.idx.21]-mean(interp.y.1[-bad.idx.21])))/(length(y.2[-bad.idx.21])-1)
    } else if (length(bad.idx.21) == 0){
      iccf.vals.21[i] <- sum((y.2-mean(y.2))*(interp.y.1-mean(interp.y.1)))/(sd(y.2-mean(y.2))*sd(interp.y.1-mean(interp.y.1)))/(length(y.2)-1)
    }

    iccf.vals[i] <- mean(c(iccf.vals.12[i], iccf.vals.21[i]))
  }

  return(data.frame(tau = tau.vals, iccf = iccf.vals))
}
