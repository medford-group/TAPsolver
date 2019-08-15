#' Gas Concentration
#'
#' Optimization of the Y-Procedure Smoothing value
#'
#' blah
#'
#' @param TAPobj A set of pre-processing parameters based on the Y-Procedure
#' @param method Blah
#' @return Blah
#' @examples
#' TAPexperiment = readTAP("data/pumpProbe.tdms")
#' TAPexperiment$AMU40$options$baselineStart = 3.9
#' TAPexperiment$AMU40$options$baselineStart = 4
#' TAPexperiment$AMU40$options$baselineType = "mean"
#'
#' TAPexperiment = baselineCorrect(TAPexperiment, "AMU40")
#' TAPexperiment = gasCon(TAPexperiment, "AMU40")
#' plot(TAPexperiment$AMU40gasCon$matrix[[5]])
#'
#' TAPexperiment = gasCon(TAPexperiment, "AMU40", method = "g")
#' plot(TAPexperiment$AMU40gasCon$matrix[[5]])
#'
#'@export gasCon

gasCon = function(TAPexperiment, gasName, method = "y"){
  timeMinPos = which.min(abs(TAPexperiment$parameters$experiment$time - TAPexperiment$parameters$experiment$timeStart))
  timeMaxPos = which.min(abs(TAPexperiment$parameters$experiment$time - TAPexperiment$parameters$experiment$timeEnd))
  timeIndex = timeMinPos:timeMaxPos
  currentTime = TAPexperiment$parameters$experiment$time[timeIndex]


  for(i in 1:length(gasName)){
    if(is.null(TAPexperiment[[gasName[i]]]$gasScalar)){
      TAPexperiment = initYproc(TAPexperiment, gasName[i], method)
    }
    TAPobj = TAPexperiment[[gasName[i]]]

    if(method == "y"){

      for(j in 1:TAPexperiment$parameters$experiment$numPulses){
        result = Re(fft(fft(TAPobj$matrix[[j]][timeIndex]) * TAPobj$gasScalar, inverse = T)) * TAPobj$gasConScalar[j]

        if(length(result) != TAPexperiment$parameters$experiment$numTime){
          trailingMedian = rep(median(result[round(.9*length(result)):length(result)]), TAPexperiment$parameters$experiment$numTime - length(result))
          result = c(result, trailingMedian)
        }
        TAPobj$matrix[[j]] = result
      }

    }else{
      gasScalar = 1 - (TAPexperiment$parameters$reactor$inertZone2Length / TAPexperiment$parameters$reactor$reactorLength)^2 / 6
      basicShape =  3 / 2

      for(j in 1:TAPexperiment$parameters$experiment$numPulses){
        tempPulse = TAPobj$matrix[[j]][timeIndex] / TAPobj$moments[[j]]$M0

        tempMean = mean(tempPulse * currentTime * TAPexperiment$parameters$experiment$timeEnd) * pi / 4
        if(mean(tempPulse * currentTime * TAPexperiment$parameters$experiment$timeEnd) < 0){
          # use median if mean is below zero due to noise
          tempMean = currentTime[which.max(tempPulse)] * 3 * pi / 4
        }
        pulseGamma = pgamma(currentTime, shape = basicShape, scale = tempMean * 2/3)
        gasGamma = pgamma(currentTime, shape = gasScalar * basicShape, scale = tempMean * 2/3)
        gasRatio = gasGamma / pulseGamma

        result = c(0, diff(cumsum(tempPulse) * gasRatio)) * TAPobj$moments[[j]]$M0 * TAPobj$gasConScalar[j]
        result[is.na(result)] = 0

        if(length(result) != TAPexperiment$parameters$experiment$numTime){
          trailingMedian = rep(median(result[round(.9*length(result)):length(result)]), TAPexperiment$parameters$experiment$numTime - length(result))
          result = c(result, trailingMedian)
        }
        TAPobj$matrix[[j]] = result
      }
    }
    TAPobj$options$name = paste0(gasName[i], "gasCon")
    TAPexperiment[[TAPobj$options$name]] = TAPobj
  }

  return(TAPexperiment)
}








