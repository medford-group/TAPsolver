#' Wavelet Smoothing
#'
#' Smooth the \code{pulseVector} by an integer size of the filter length.
#'
#' Savitzky-Golay performs smoothing by using a least squares fit on small subsections of the data to a polynomial curve.
#' As such, the algorithm requires a choice of window sizes of the subsections of the data.  As the window size increases, the pulse becomes more smooth.
#' This algorithm should be used minimally due to the decrease of the magnitude of the pulse.
#'
#' This is a wrapper function on the savgol function from the package pracma.
#' Options on the savgol function allow for the filter order to be managed e.g. 2 = quadratic filter and 4 = quartic (currently using 2).
#' The derivative order may be changed as well e.g. 0 = smoothing and 1 = first derivative and so on.
#' The function is currently setup so that it does quadratic smoothing.
#'
#' @param TAPexperiment A TAP object of pulse values.
#' @param gasName Value
#' @return A vector of smoothed pulse values.
#' @examples
#'
#' TAPexperiment = readTAP("data/pumpProbe.tdms")
#' TAPexperiment$AMU40$options$baselineStart = 3.9
#' TAPexperiment$AMU40$options$baselineStart = 4
#' TAPexperiment$AMU40$options$baselineType = "mean"
#'
#' TAPexperiment = baselineCorrect(TAPexperiment, "AMU40")
#'
#' TAPexperiment$AMU40$options$waveSmoothing = 8
#' TAPexperiment = waveTAP(TAPexperiment, "AMU40")
#' plot(TAPexperiment$AMU40$matrix[[5]])
#'
#'
#' @export waveTAP

waveTAP = function(TAPexperiment, gasName){
  for(i in 1:length(gasName)){
    TAPobj = TAPexperiment[[gasName[i]]]

    if(TAPobj$options$waveSmoothing == 0)
      next

    helperOptions = list()
    helperOptions$waveSmoothing = TAPobj$options$waveSmoothing

    waveHelper = function(x,helperOptions){
      originalLength = length(x)
      if(is.integer(log2(length(x))) == F){
        x = c(x, rep(x[length(x)], 2^(ceiling(log2(length(x)))) - length(x) ))
      }
      waveletDecomp = wavethresh::wd(family = "DaubExPhase",data = x, filter.number = helperOptions$waveSmoothing)
      waveletThresholdedTune = wavethresh::threshold(waveletDecomp, policy="universal", dev = wavethresh::madmad, type = "soft",
                                                     levels = wavethresh::nlevelsWT(waveletDecomp) - 1, return.threshold = TRUE, verbose = F)
      waveletThreshold = wavethresh::threshold(waveletDecomp, policy = "manual", value= waveletThresholdedTune, type = "soft",
                                               levels = 0:(wavethresh::nlevelsWT(waveletDecomp) - 1))
      return(wavethresh::wr(waveletThreshold)[1:originalLength])
    }


    TAPobj$matrix = lapply(TAPobj$matrix, FUN = waveHelper, helperOptions = helperOptions)


    TAPexperiment[[gasName[i]]]$matrix = TAPobj$matrix
    tempMat = fluxToDataFrame(TAPexperiment, gasName[i])

    for(j in 1:TAPexperiment$parameters$experiment$numTime){
      tempPulse = unlist(tempMat[j, ],F,F)
      originalLength = length(tempPulse)
      if(is.integer(log2(length(tempPulse))) == F){
        tempPulse = c(tempPulse, rep(tempPulse[length(tempPulse)], 2^(ceiling(log2(length(tempPulse)))) - length(tempPulse) ))
      }
      waveletDecomp = wavethresh::wd(family = "DaubExPhase",data = tempPulse, filter.number = TAPobj$options$waveSmoothing)
      waveletThresholdedTune = wavethresh::threshold(waveletDecomp, policy="universal", dev = wavethresh::madmad, type = "soft",
                                                     levels = wavethresh::nlevelsWT(waveletDecomp) - 1, return.threshold = TRUE, verbose = F)
      waveletThreshold = wavethresh::threshold(waveletDecomp, policy = "manual", value= waveletThresholdedTune, type = "soft",
                                               levels = 0:(wavethresh::nlevelsWT(waveletDecomp) - 1))
      tempMat[j, ] = wavethresh::wr(waveletThreshold)[1:originalLength]
    }

    TAPexperiment[[gasName[i]]]$matrix = as.list(tempMat)

    #TAPexperiment[[TAPobj$options$name]] = TAPobj
  }
  return(TAPexperiment)
}

