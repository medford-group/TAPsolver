#' Moments of an TAP Pulse
#'
#' Function for determining the moments of a pulse.
#'
#'The \eqn{n^{th}} moment is defined by:
#'\deqn{ \int t^n pulse(t) dt}
#'
#'
#'
#' @param TAPobj A TAP pulse or series of pulses.
#' @return A vector of moments based on the pulse.
#' @examples
#'
#' TAPexperiment = readTAP("data/pumpProbe.tdms")
#' TAPexperiment$AMU40$options$baselineStart = 3.9
#' TAPexperiment$AMU40$options$baselineStart = 4
#' TAPexperiment$AMU40$options$baselineType = "mean"
#'
#' TAPexperiment = baselineCorrect(TAPexperiment, "AMU40")
#' TAPexperiment = moments(TAPexperiment, "AMU40")
#' df = momentsToDataFrame(TAPexperiment, "AMU40")
#' plot(df$M1)
#'
#'@export moments
moments = function(TAPexperiment, gasName){

  momentsHelper = function(x, helperOptions){
    result = list("max" = 0, "baselineShift" = 0, "M0" = 0,
                  "M1" = 0, "M2" = 0, "M0norm" = 0, "M1norm" = 0,
                  "M2norm" = 0, "M1M0" = 0, "M2M0" = 0, "mode" = 0,
                  "knudsenRatio" = 0,
                  "diffusion" = 0)
    result$max = max(x)
    result$baselineShift = mean(x[baseIndex])
    result$M0 = pracma::trapz(x[timeIndex]) * helperOptions$maxTime / helperOptions$lenTime
    result$M1 = pracma::trapz((x * helperOptions$time)[timeIndex]) * helperOptions$maxTime / helperOptions$lenTime
    result$M2 = pracma::trapz((x * helperOptions$time^2)[timeIndex]) * helperOptions$maxTime / helperOptions$lenTime
    result$M0norm = result$M0 / result$M0
    result$M1norm = result$M1 / result$M0
    result$M2norm = result$M2 / result$M0
    result$M1M0 = result$M1 / result$M0
    result$M2M0 = result$M2 / result$M0
    result$mode = tempTime[which.max(x[helperOptions$timeIndex])]
    result$knudsenRatio = result$M1norm / result$mode

    if(helperOptions$diffusionCoef != 0){
      result$diffusion = helperOptions$diffusionCoef
    }else{
      result$diffusion = helperOptions$diffusionScalar * result$M0 / (2 * result$M1)
    }

    return(result)
  }


  tempParams = TAPexperiment$reactor$reactorParams
  timeMinPos = which.min(abs(TAPexperiment$parameters$experiment$time - TAPexperiment$parameters$experiment$timeStart))
  timeMaxPos = which.min(abs(TAPexperiment$parameters$experiment$time - TAPexperiment$parameters$experiment$timeEnd))

  timeIndex = timeMinPos:timeMaxPos
  tempTime = TAPexperiment$parameters$experiment$time[timeMinPos:timeMaxPos]

  helperOptions = list()
  helperOptions$time = tempTime
  helperOptions$timeIndex = timeIndex

  helperOptions$maxTime = max(tempTime)
  helperOptions$lenTime = length(tempTime)
  helperOptions$diffusionCoef = TAPexperiment$parameters$experiment$diffusionCoef
  helperOptions$diffusionScalar = TAPexperiment$parameters$reactor$bedPorosity * TAPexperiment$parameters$reactor$reactorLength^2

  for(i in 1:length(gasName)){
    TAPobj = TAPexperiment[[gasName[i]]]
    baseMinPos = which.min(abs(TAPexperiment$parameters$experiment$time - TAPobj$options$baselineStart))
    baseMaxPos = which.min(abs(TAPexperiment$parameters$experiment$time - TAPobj$options$baselineEnd))
    baseIndex = baseMinPos:baseMaxPos
    helperOptions$baseIndex = baseIndex

    TAPobj$moments = lapply(TAPobj$matrix, FUN = momentsHelper, helperOptions = helperOptions)
    for( j in 1:TAPexperiment$parameters$experiment$numPulses){
      TAPobj$moments[[j]]$temperature = TAPobj$temperature[j] #TAPexperiment$parameters$experiment$temperature[[j]]
    }

    TAPexperiment[[TAPobj$options$name]] = TAPobj

  }


  return(TAPexperiment)
}
