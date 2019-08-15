#' Reaction Rate
#'
#' Blah
#'
#' blah
#'
#' @param TAPexperiment A set of pre-processing parameters based on the Y-Procedure
#' @param gasName Blah
#' @param method Blah
#' @return Blah
#' @examples
#' TAPexperiment = readTAP("data/pumpProbe.tdms")
#' TAPexperiment$AMU40$options$baselineStart = 3.9
#' TAPexperiment$AMU40$options$baselineStart = 4
#' TAPexperiment$AMU40$options$baselineType = "mean"
#'
#' TAPexperiment = baselineCorrect(TAPexperiment, "AMU40")
#'
#' TAPexperiment$AMU40$options$isProduct = T
#'
#' TAPexperiment = reactionRate(TAPexperiment, "AMU40")
#' plot(TAPexperiment$AMU40rate$matrix[[5]])
#'
#' TAPexperiment = reactionRate(TAPexperiment, "AMU40", method = "g")
#' plot(TAPexperiment$AMU40rate$matrix[[5]])
#'@export reactionRate
#'
reactionRate = function(TAPexperiment, gasName, method = "y"){

  timeMinPos = which.min(abs(TAPexperiment$parameters$experiment$time - TAPexperiment$parameters$experiment$timeStart))
  timeMaxPos = which.min(abs(TAPexperiment$parameters$experiment$time - TAPexperiment$parameters$experiment$timeEnd))
  timeIndex = timeMinPos:timeMaxPos
  currentTime = TAPexperiment$parameters$experiment$time[timeIndex]

  for(i in 1:length(gasName)){
    if(is.null(TAPexperiment[[gasName[i]]]$rateScalar)){
      TAPexperiment = initYproc(TAPexperiment, gasName[i], method)
    }
    TAPobj = TAPexperiment[[gasName[i]]]
    TAPobj$options$name = paste0(TAPobj$options$name, "rate")

    if(is.na(TAPobj$options$isProduct))
      TAPobj$options$isProduct = F

    if(!(TAPobj$options$isProduct)){
      for(j in 1:TAPexperiment$parameters$experiment$numPulses){
        TAPobj$matrix[[j]] = TAPexperiment[[TAPobj$options$inert]]$matrix[[j]] - TAPobj$matrix[[j]]
      }

      TAPexperiment[[TAPobj$options$name]] = TAPobj
      TAPexperiment = moments(TAPexperiment, TAPobj$options$name)
      TAPobj = TAPexperiment[[TAPobj$options$name]]
    }

    if(method == "y"){
      for(j in 1:TAPexperiment$parameters$experiment$numPulses){
        tempPulse = TAPobj$matrix[[j]][timeIndex]

        result = Re(fft(fft(tempPulse) * TAPobj$rateScalar, inverse = T)) * TAPobj$reactionRateScalar[j]
        if(length(tempPulse) != TAPexperiment$parameters$experiment$numTime){
          trailingMedian = rep(median(result[round(.9*length(result)):length(result)]), TAPexperiment$parameters$experiment$numTime - length(result))
          result = c(result, trailingMedian)
        }
        TAPobj$matrix[[j]] = result
      }
    }else{
      rateScalar = 1 - 1.5 * (TAPexperiment$parameters$reactor$inertZone2Length / TAPexperiment$parameters$reactor$reactorLength)^2
      basicShape =  3 / 2

      for(j in 1:TAPexperiment$parameters$experiment$numPulses){
        tempPulse = TAPobj$matrix[[j]][timeIndex] / TAPobj$moments[[j]]$M0

        tempMean = mean(tempPulse * currentTime * TAPexperiment$parameters$experiment$timeEnd) * pi / 4
        if((mean(tempPulse * currentTime * TAPexperiment$parameters$experiment$timeEnd) < 0)){
          # use median if mean is below zero due to noise
          tempMean = currentTime[TAPcodeV2::robustMode(tempPulse)] * 3 * pi / 4
        }
        pulseGamma = pgamma(currentTime, shape = basicShape, scale = tempMean * 2/3)
        gasGamma = pgamma(currentTime, shape = rateScalar * basicShape, scale = tempMean * 2/3)
        gasRatio = gasGamma / pulseGamma

        result = c(0, diff(cumsum(tempPulse) * gasRatio)) * TAPobj$moments[[j]]$M0 * TAPobj$reactionRateScalar[j]
        result[is.na(result)] = 0

        if(length(tempPulse) != TAPexperiment$parameters$experiment$numTime){
          trailingMedian = rep(median(result[round(.9*length(result)):length(result)]), TAPexperiment$parameters$experiment$numTime - length(result))
          result = c(result, trailingMedian)
        }
        TAPobj$matrix[[j]] = result
      }
    }

    TAPexperiment[[TAPobj$options$name]] = TAPobj
  }
  return(TAPexperiment)
}
