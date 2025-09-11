library(foreach)
library(parallel)
library(doParallel)
library(optimx)


#' lc.uvot <- function(file.path, flux.units){
#'   #' @title Reads in uvotsource output files and makes a light curve data frame
#'   #' @description This function will read in all of the output files produced by
#'   #' the uvotsource HEASoft tool that are stored in the given file path.
#'   #' @param file.path System file path to the directory that contains the output
#'   #' files from uvotsource.
#'   #' @param flux.units Light curve flux units. Options are: rate (counts/second)
#'   #' , flux (ergs/s/cm^2/A), and mJy (milliJansky).
#'   #' @returns A data frame with columns for observation ID, UVOT telescope
#'   #' filter, time (MET), time error (MET), flux, and flux error.
#'   #' @export
#'
#'   working.path <- getwd()
#'
#'   ### Determine what flux and error columns to grab based on user input
#'   if (flux.units == "rate"){
#'     val <- "CORR_RATE"
#'     err <- "CORR_RATE_ERR"
#'   } else if (flux.units == "flux"){
#'     val <- "FLUX_AA"
#'     err <- "FLUX_AA_ERR"
#'   } else if (flux.units == "mJy"){
#'     val <- "FLUX_HZ"
#'     err <- "FLUX_HZ_ERR"
#'   } else {
#'     cat("***flux units not recognized!\nplease enter 'rate', 'flux', or 'Jy'\n")
#'     break
#'   }
#'
#'   ### Read in all of the file names, splitting in OBSID and filter information
#'   setwd(file.path)
#'   all.files <- list.files(pattern = ".fits")
#'   df.files <- data.frame(FILENAME = rep(NA, length.out = length(all.files)))
#'   df.files$OBSID <- NA
#'   df.files$FILTER <- NA
#'   for (i in 1:length(all.files)){
#'     df.files$FILENAME[i] <- all.files[i]
#'     df.files$OBSID[i] <- strsplit(all.files[i], "_")[[1]][1]
#'     df.files$FILTER[i] <- strsplit(strsplit(all.files[i], "_")[[1]][2], ".fits")[[1]][1]
#'   }
#'   rm(all.files)
#'
#'   ### Create the full UVOT light curve
#'   obs.ids <- unique(df.files$OBSID)
#'   lc <- data.frame(OBSID = NA,
#'                    FILTER = NA,
#'                    TIME = NA_real_, TIME.ERR = NA_real_,
#'                    RATE = NA_real_, RATE.ERR = NA_real_)
#'
#'   for (k in 1:length(obs.ids)){
#'   # for (k in 1:50){
#'     cat(paste("Working on observation:\t", k, " / ", length(obs.ids)), "\r")
#'     tmp.files <- df.files[which(df.files$OBSID == obs.ids[k]),]
#'
#'     if (length(tmp.files$OBSID) == 0){
#'       next
#'     }
#'
#'     for (j in 1:length(tmp.files$FILENAME)){
#'       file.tmp <- FITSio::readFITS(tmp.files$FILENAME[j])
#'       idx.val <- which(file.tmp$colNames == val)
#'       idx.err <- which(file.tmp$colNames == err)
#'       lc <- rbind(lc, data.frame(OBSID = tmp.files$OBSID[j],
#'                                  FILTER = tmp.files$FILTER[j],
#'                                  TIME = file.tmp$col[[1]], TIME.ERR = file.tmp$col[[4]]-file.tmp$col[[3]],
#'                                  RATE = file.tmp$col[[idx.val]], RATE.ERR = file.tmp$col[[idx.err]]))
#'     }
#'     rm(list = ls(pattern = "tmp"))
#'   }
#'   lc <- lc[-1,]
#'   lc <- lc[order(lc$TIME),]
#'
#'   setwd(working.path)
#'   return(lc)
#' }


lc.uvot <- function(file.path){
  #' @title Reads in uvotsource output files and makes a light curve data frame
  #' @description This function will read in all of the output files produced by
  #' the uvotsource HEASoft tool that are stored in the given file path.
  #' @param file.path System file path to the directory that contains the output
  #' files from uvotsource.
  #' @returns A data frame with columns for observation ID, UVOT telescope
  #' filter, time (MET), time error (MET), flux, and flux error.
  #' @export

  # Read in all of the file names, splitting in OBSID and filter information
  working.path <- getwd()
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

  # Create the full UVOT light curve
  obs.ids <- unique(df.files$OBSID)
  lc <- data.frame(OBSID = NA,
                   FILTER = NA,
                   TIME = NA_real_, TIME.ERR = NA_real_,
                   RATE = NA_real_, RATE.ERR = NA_real_,
                   FLUX.AA = NA_real_, FLUX.AA.ERR = NA_real_,
                   FLUX.HZ = NA_real_, FLUX.HZ.ERR = NA_real_
                   )

  for (k in 1:length(obs.ids)){
    cat(paste("Working on observation:\t", k, " / ", length(obs.ids)), "\r")
    tmp.files <- df.files[which(df.files$OBSID == obs.ids[k]),]

    if (length(tmp.files$OBSID) == 0){
      next
    }

    for (j in 1:length(tmp.files$FILENAME)){
      file.tmp <- FITSio::readFITS(tmp.files$FILENAME[j])
      lc <- rbind(lc, data.frame(OBSID = tmp.files$OBSID[j],
                                 FILTER = tmp.files$FILTER[j],
                                 TIME = file.tmp$col[[1]], TIME.ERR = file.tmp$col[[4]]-file.tmp$col[[3]],
                                 RATE = file.tmp$col[[which(file.tmp$colNames == "CORR_RATE")]], RATE.ERR = file.tmp$col[[which(file.tmp$colNames == "CORR_RATE_ERR")]],
                                 FLUX.AA = file.tmp$col[[which(file.tmp$colNames == "FLUX_AA")]], FLUX.AA.ERR = file.tmp$col[[which(file.tmp$colNames == "FLUX_AA_ERR")]],
                                 FLUX.HZ = file.tmp$col[[which(file.tmp$colNames == "FLUX_HZ")]], FLUX.HZ.ERR = file.tmp$col[[which(file.tmp$colNames == "FLUX_HZ_ERR")]]
                                 )
                  )
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


calc.fvar <- function(light.curve, method = "Vaughan"){
  #' @title Compute fractional variability of light curve
  #' @description Takes the input data frame light curve and compute the fractional variability as in Vaughan et al. (2003), MNRAS, 345, 1271-1284.
  #' @param light.curve data frame of light curve
  #' @param method method used to compute fractional variability
  #' @return Data frame of fractional variability and error.
  #' @export

  y.vals <- light.curve$RATE
  y.errs <- light.curve$RATE.ERR
  s.sq <- sum((y.vals - mean(y.vals))^2) / (length(y.vals)-1)

  if (method == "Edelson"){
    ### Edelson+2002
    f.var <- sqrt(s.sq-mean(y.errs^2))/mean(y.vals)
    f.var.err <- (1/f.var) * sqrt(1/(2*length(y.vals))) * (s.sq/(mean(y.vals)^2))
  } else if (method == "Vaughan"){
    ### Vaughan+2003
    sig.sq <- sum(y.errs^2) / length(y.errs)
    f.var <- sqrt((s.sq - sig.sq) / mean(y.vals)^2)
    f.var.err <- sqrt((sqrt(1 / (2*length(y.vals))) * sig.sq / (mean(y.vals)^2 * f.var))^2 + (sqrt(sig.sq/length(y.vals)) * 1/mean(y.vals))^2)
  }

  return(data.frame(FVAR=f.var, FVAR.ERR=f.var.err))
}


detrend.SavitzkyGolay <- function(light.curve, filter.width = NULL, polynomial.order = 1, output = 'detrended'){
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

  if (length(which(sf$val[-1] < 0)) > 0 | length(which(sf$val[-1] < sf$err[-1])) > 0){
    cat("***be careful with this SF, it has bad values!", "\n")
  }

  return(sf)
}


broken.power.law <- function(params, x){
  #' @title Broken power law model
  #' @description A broken power law model for fitting structure functions
  #' @param params list of parameters for the broken power law model (1: scale,
  #' 2: slope, 3: break time)
  #' @param x time values
  #' @return model y (i.e structure function) values
  #' @export

  x.prebreak <- x[x <= params[3]]
  x.postbreak <- x[x > params[3]]
  y.prebreak <- (10^params[1])*(x.prebreak^params[2])
  y.postbreak <- rep((10^params[1])*(x.prebreak[length(x.prebreak)]^params[2]), length.out = length(x.postbreak))
  return(c(y.prebreak, y.postbreak))
}


calc.chisq <- function(model.values, y, e){
  return(sum( (y - model.values)^2 / e^2 ))
}


chisq.brknpwrlw <- function(par, data){
  #' @title Chi-squared calculator
  #' @description Computes the chi-squared value of a model evaluated against some data
  #' @param par list of parameters for the model; note that this MUST be the first argument!
  #' @param data the data to be used with columns for x, y, and e (error)
  #' @return chi-squared value
  #' @export

  return(calc.chisq(broken.power.law(par, data$x), data$y, data$e))
}


fit.sf.brknpwrlw <- function(sf, start.params, tchar.vals){
  min.chisq <- Inf
  for (tchar in tchar.vals){
    fit.results <- optim(par = start.params,
                         fn = chisq.brknpwrlw,
                         data = data.frame(x = sf$tau, y = sf$val, e = sf$err)
                         )
    if (fit.results$value < min.chisq){
      min.chisq <- fit.results$value
      min.params <- fit.results$par
    }
  }
  return(min.params)
}


calc.iccf <- function(light.curve.1, light.curve.2, delta.tau, max.lag = NA){
  #' @title Compute the interpolation cross-correlation function
  #' @description Takes the input data frames of two light curves and computes the interpolation cross-correlation function (ICCF) as in Gaskell & Peterson (1987).
  #' @param light.curve.1 data frame of first light curve
  #' @param light.curve.2 data frame of second light curve
  #' @param delta.tau time step to sample lags over
  #' @param max.lag maximum lag to calculate
  #' @return Data frame of the ICCF
  #' @export

  ### Extract the light curve time and count rates
  x.1 <- light.curve.1$TIME - light.curve.1$TIME[1]
  y.1 <- light.curve.1$RATE
  x.2 <- light.curve.2$TIME - light.curve.1$TIME[1]
  y.2 <- light.curve.2$RATE

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


TK95.simulate.lightcurve <- function(beta, bins = 1024, length = 100000, scale.factor = 1, shift.factor = 0, error.scale = 0.1) {
  bins <- 2*bins
  time <- seq(1,2*length, length.out = bins)
  fourier.frequencies <- seq(1,bins)/length
  step.one <- rnorm(bins)
  step.two <- (1/fourier.frequencies)^(beta/2.0)
  step.three <- step.one*step.two
  step.four <- fft(step.three, inverse = TRUE)
  step.five <- Re(step.four)
  simulated.lc <- data.frame(TIME = time, TIME.ERR = time, RATE = step.five, RATE.ERR = step.five)
  simulated.lc$TIME <- simulated.lc$TIME - simulated.lc$TIME[1]
  simulated.lc$TIME.ERR <- (simulated.lc$TIME[2]-simulated.lc$TIME[1])/2
  simulated.lc <- subset(simulated.lc, simulated.lc$TIME < length)
  simulated.lc$RATE <- simulated.lc$RATE - mean(simulated.lc$RATE)
  simulated.lc$RATE <- simulated.lc$RATE / max(abs(simulated.lc$RATE))
  simulated.lc$RATE <- simulated.lc$RATE * scale.factor
  simulated.lc$RATE <- simulated.lc$RATE + shift.factor
  simulated.lc$RATE.ERR <- mean(simulated.lc$RATE) * error.scale
  return(simulated.lc)
}


impose.gappy.cadence <- function(parent, child){
  idx <- NA
  for (i in 1:length(parent$TIME)){
    idx <- c(idx, which.min(abs(child$TIME - parent$TIME[i])))
  }
  idx <- idx[-1]
  gappy.child <- child[idx,]
  gappy.child <- gappy.child[which(duplicated(gappy.child) == F),]
  return(gappy.child)
}


iccf.sims.inparallel <- function(iccf, light.curve.1, light.curve.2, delta.tau, max.lag, beta, length, bins, n.sims, n.cores,
                                 detrend = FALSE, filter.width = NULL, polynomial.order = 1){
  cl <- makeCluster(n.cores[1])
  registerDoParallel(cl)
  all.sim.iccfs.parallel <- foreach(i=1:n.sims,
                                    .combine=cbind,
                                    .export = c("calc.iccf",
                                                "detrend.SavitzkyGolay",
                                                "impose.gappy.cadence",
                                                "TK95.simulate.lightcurve"
                                                )
                                    ) %dopar% {
    slc <- TK95.simulate.lightcurve(beta, length = length, bins = bins, scale.factor = var(light.curve.1$RATE), shift.factor = mean(light.curve.1$RATE), error.scale = mean(light.curve.1$RATE.ERR/light.curve.1$RATE))
    slc <- impose.gappy.cadence(light.curve.1, slc)
    if (detrend == TRUE & is.null(filter.width) == FALSE){
      slc <- detrend.SavitzkyGolay(slc, filter.width = filter.width, polynomial.order = polynomial.order, output = 'detrended')
    }
    sim.iccf <- calc.iccf(slc, light.curve.2, delta.tau = delta.tau, max.lag = max.lag)
    all.sim.iccfs <- sim.iccf$iccf

    all.sim.iccfs
  }
  stopCluster(cl)
  all.sim.iccfs <- all.sim.iccfs.parallel
  rm(all.sim.iccfs.parallel)

  ### Compile results into iccf data frame
  length(which(is.nan(all.sim.iccfs) == T))
  iccf$p99 <- NA_real_
  for (i in 1:length(iccf$tau)){
    iccf$p99[i] <- quantile(all.sim.iccfs[i,], 0.99, na.rm = T)
  }

  return(iccf)
}


iccf.centroid.inparallel <- function(iccf, light.curve.1, light.curve.2, delta.tau, max.lag, peak.width = NULL, n.sims, n.cores){
  if (is.null(peak.width) == TRUE){
    peak.width <- max.lag
  }

  cl <- makeCluster(n.cores[1])
  registerDoParallel(cl)

  centroids.parallel <- foreach(i=1:n.sims,
                                .combine=rbind,
                                .export = c("calc.iccf")
                                ) %dopar% {
    base.tmp <- light.curve.1
    compare.tmp <- light.curve.2

    ### perform flux randomisation using error bars as st.dev. with mean = 0
    base.tmp$RATE <- base.tmp$RATE + rnorm(n = length(base.tmp$RATE.ERR), mean = 0, sd = base.tmp$RATE.ERR)
    compare.tmp$RATE <- compare.tmp$RATE + rnorm(n = length(compare.tmp$RATE.ERR), mean = 0, sd = compare.tmp$RATE.ERR)

    ### resample the light curves keeping only unique times
    base.rand.idx <- floor(runif(length(base.tmp$TIME), min=1, max=length(base.tmp$TIME)+1))
    compare.rand.idx <- floor(runif(length(compare.tmp$TIME), min=1, max=length(compare.tmp$TIME)+1))
    base.tmp <- base.tmp[unique(sort(base.rand.idx)),]
    compare.tmp <- compare.tmp[unique(sort(compare.rand.idx)),]

    ### compute the new iccf
    tmp.iccf <- calc.iccf(base.tmp, compare.tmp, delta.tau = delta.tau, max.lag = max.lag)
    tmp.iccf <- tmp.iccf[which(is.nan(iccf$iccf) == F),]

    ### extract tmp.iccf time bins where it is > 99% contour of REAL iccf, then fit
    tmp.peak.idxs <- which(tmp.iccf$iccf > iccf$p99 & abs(tmp.iccf$tau) <= peak.width)
    tmp.iccf.peak <- tmp.iccf[tmp.peak.idxs,]
    tmp.iccf.peak <- tmp.iccf.peak[which(is.na(tmp.iccf.peak$iccf) == F),]
    centroids <- weighted.mean(tmp.iccf.peak$tau, w = tmp.iccf.peak$iccf)

    centroids
  }
  stopCluster(cl)
  centroids <- unname(centroids.parallel[,1])
  rm(centroids.parallel)
  centroids <- centroids[which(is.nan(centroids) == F)]

  return(centroids)
}


iccf.pipeline <- function(light.curve.1, light.curve.2, delta.tau, max.lag = NA, peak.width = NULL,
                          beta, length, bins,
                          n.sims, n.cores,
                          detrend = FALSE, filter.width = NULL, polynomial.order = 1){

  iccf <- calc.iccf(light.curve.1, light.curve.2, delta.tau = delta.tau, max.lag = max.lag)
  iccf <- iccf.sims.inparallel(iccf = iccf, light.curve.1 = light.curve.1, light.curve.2 = light.curve.2,
                               delta.tau = delta.tau, max.lag = max.lag,
                               beta = beta, length = length, bins = bins,
                               n.sims = n.sims, n.cores = n.cores,
                               detrend = detrend, filter.width = filter.width, polynomial.order = polynomial.order)
  centroid <- iccf.centroid.inparallel(iccf = iccf, light.curve.1 = light.curve.1, light.curve.2 = light.curve.2,
                                       delta.tau = delta.tau, max.lag = max.lag, peak.width = peak.width,
                                       n.sims = n.sims, n.cores = n.cores)
  return(list("iccf" = iccf, "centroid" = centroid))
}


deredden <- function(light.curve, wavelength, EBmV, RV = 3.1, update = F){
  #' @title De-redden optical/UV data
  #' @description This function reads in a data-frame of wavelengths (in micrometers) and optical/UV magnitudes to de-redden them according to Cardelli+1989
  #' @param df data frame with columns 'wavelength' (in micrometers) and 'flux' (counts per second)
  #' @param EBmV foreground Galactic extinction in V-band (i.e. from NASA/IPAC Extragalactic Database, NED)
  #' @param RV extinction law (default = 3.1 for Milky Way extinction)
  #' @param update boolean (TRUE or FALSE) whether or not to include the O'Donnell (1994) correction (default = F)
  #' @return Returns the input data frame with an additional 'dered' column of de-reddened values
  #' @export

  AV <- EBmV*RV
  x <- 1/(wavelength/1e4)
  if (x >= 1.1 & x < 3.3){
    y <- x-1.82
    if (update == T){
      acoeff <- 1 + 0.104*y - 0.609*y^2 + 0.701*y^3 + 1.137*y^4 - 1.718*y^5 - 0.827*y^6 + 1.647*y^5 - 0.505*y^8
      bcoeff <- 1.952*y + 2.908*y^2 - 3.989*y^3 - 7.985*y^4 + 11.102*y^5 + 5.941*y^6 - 10.805*y^7 + 3.347*y^8
    } else if (update == F){
      acoeff <- 1 + 0.17699*y - 0.50447*y^2 - 0.02427*y^3 + 0.72085*y^4 + 0.01979*y^5 - 0.77530*y^6 + 0.32999*y^7
      bcoeff <- 1.41338*y + 2.28305*y^2 + 1.07233*y^3 - 5.38434*y^4 - 0.62251*y^5 + 5.30260*y^6 - 2.09002*y^7
    }
  } else if (x >= 3.3 & x < 5.9){
    acoeff <- 1.752 - 0.316*x - 0.104/((x-4.67)^2 + 0.341)
    bcoeff <- -3.090 + 1.825*x + 1.206/((x-4.62)^2 + 0.263)
  } else if (x >= 5.9 & x < 8){
    acoeff <- 1.752 - 0.316*x - 0.104/((x-4.67)^2 + 0.341) - 0.04473*(x-5.9)^2 - 0.009779*(x-5.9)^3
    bcoeff <- -3.090 + 1.825*x + 1.206/((x-4.62)^2 + 0.263) + 0.2130*(x-5.9)^2 + 0.1207*(x-5.9)^3
  } else if (x >= 8 & x < 10){
    acoeff <- -1.073 - 0.628*(x-8) + 0.137*(x-8)^2 - 0.070*(x-8)^3
    bcoeff <- 13.670 + 4.257*(x-8) - 0.420*(x-8)^2 + 0.374*(x-8)^3
  } else if (x >= 10){
    cat("...wavelength not covered by curve (too short)!\n")
    break
  } else if (x < 1.1){
    cat("...wavelength not covered by curve (too long)!\n")
    break
  } else {
    cat("...check that the wavelength units are micrometers!\n")
    break
  }
  Alambda <- AV*(acoeff + bcoeff/RV)
  light.curve$RATE <- light.curve$RATE*10^(0.4*Alambda)
  light.curve$RATE.ERR <- light.curve$RATE.ERR*10^(0.4*Alambda)

  return(light.curve)
}


normalize.lightcurve <- function(light.curve){
  #' @title Normalize an input light curve
  #' @description Subtract the RATE mean and divide by RATE standard deviation
  #' @param light.curve data frame of the light curve with at least columns of RATE and RATE.ERR
  #' @return Data frame of the normlized light curve
  #' @export

  avg <- mean(light.curve$RATE)
  std <- sd(light.curve$RATE)
  light.curve$RATE <- (light.curve$RATE-avg)/std
  light.curve$RATE.ERR <- light.curve$RATE.ERR/std
  return( light.curve )
}


closest.idx <- function(light.curve, approx.times){
  idx <- NA_real_
  for (i in 1:length(light.curve$TIME)){
    idx <- c(idx, which.min(abs(approx.times-light.curve$TIME[i])))
  }
  idx <- unique(idx[-1])
  return( idx )
}


filter.times <- function(light.curve, approx.times){
  idx <- NA_real_
  for (i in 1:length(approx.times)){
    idx <- c(idx, which.min(abs(light.curve$TIME-approx.times[i])))
  }
  idx <- unique(idx[-1])
  light.curve <- light.curve[idx,]
  return( light.curve )
}


chisq.xt <- function(par, data, xt){
  w2.model <- par[1]+par[2]*xt
  chi.w2 <- calc.chisq(w2.model, data$w2, data$w2.err)

  m2.model <- par[3]+par[4]*xt
  chi.m2 <- calc.chisq(m2.model, data$m2, data$m2.err)

  w1.model <- par[5]+par[6]*xt
  chi.w1 <- calc.chisq(w1.model, data$w1, data$w1.err)

  u.model <- par[7]+par[8]*xt
  chi.u <- calc.chisq(u.model, data$u, data$u.err)

  b.model <- par[9]+par[10]*xt
  chi.b <- calc.chisq(b.model, data$b, data$b.err)

  v.model <- par[11]+par[12]*xt
  chi.v <- calc.chisq(v.model, data$v, data$v.err)

  return( sum(c(chi.w2, chi.m2, chi.w1, chi.u, chi.b, chi.v)) )
}


fit.xt <- function(data, start.params, xt){
  fit.results <- optimx(par = start.params,
                        fn = chisq.xt,
                        data = data,
                        xt = xt,
                        method = c("BFGS")
                        )

  pars <- list("w2.avg" = fit.results$p1, "w2.rms" = fit.results$p2,
               "m2.avg" = fit.results$p3, "m2.rms" = fit.results$p4,
               "w1.avg" = fit.results$p5, "w1.rms" = fit.results$p6,
               "u.avg" = fit.results$p7, "u.rms" = fit.results$p8,
               "b.avg" = fit.results$p9, "b.rms" = fit.results$p10,
               "v.avg" = fit.results$p11, "v.rms" = fit.results$p12,
               "chisq" = fit.results$value, "red.chisq" = fit.results$value/(length(xt)*6-2*6)
               )

  return(pars)
}


generate.variability.spectrum <- function(fit.par, xt){
  w2.gal <- -fit.par$w2.avg/fit.par$w2.rms
  w2.min <- fit.par$w2.avg + fit.par$w2.rms*min(xt)
  w2.max <- fit.par$w2.avg + fit.par$w2.rms*max(xt)
  w2.dif <- w2.max - w2.min

  m2.gal <- fit.par$m2.avg + fit.par$m2.rms*w2.gal
  m2.min <- fit.par$m2.avg + fit.par$m2.rms*min(xt)
  m2.max <- fit.par$m2.avg + fit.par$m2.rms*max(xt)
  m2.dif <- m2.max - m2.min

  w1.gal <- fit.par$w1.avg + fit.par$w1.rms*w2.gal
  w1.min <- fit.par$w1.avg + fit.par$w1.rms*min(xt)
  w1.max <- fit.par$w1.avg + fit.par$w1.rms*max(xt)
  w1.dif <- w1.max - w1.min

  u.gal <- fit.par$u.avg + fit.par$u.rms*w2.gal
  u.min <- fit.par$u.avg + fit.par$u.rms*min(xt)
  u.max <- fit.par$u.avg + fit.par$u.rms*max(xt)
  u.dif <- u.max - u.min

  b.gal <- fit.par$b.avg + fit.par$b.rms*w2.gal
  b.min <- fit.par$b.avg + fit.par$b.rms*min(xt)
  b.max <- fit.par$b.avg + fit.par$b.rms*max(xt)
  b.dif <- b.max - b.min

  v.gal <- fit.par$v.avg + fit.par$v.rms*w2.gal
  v.min <- fit.par$v.avg + fit.par$v.rms*min(xt)
  v.max <- fit.par$v.avg + fit.par$v.rms*max(xt)
  v.dif <- v.max - v.min

  return(data.frame(WAVELENGTH = c(1928, 2246, 2600, 3465, 4392, 5468),
                    AVG = c(fit.par$w2.avg, fit.par$m2.avg, fit.par$w1.avg, fit.par$u.avg, fit.par$b.avg, fit.par$v.avg),
                    RMS = c(fit.par$w2.rms, fit.par$m2.rms, fit.par$w1.rms, fit.par$u.rms, fit.par$b.rms, fit.par$v.rms),
                    AGN = c(w2.dif, m2.dif, w1.dif, u.dif, b.dif, v.dif),
                    GAL = c(NA_real_, m2.gal, w1.gal, u.gal, b.gal, v.gal)
                    )
         )
}

