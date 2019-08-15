#' Standard Diffusion Curve
#'
#' Create a standard diffusion curve based on the inert gas to determine if the procedure is under the Knudsen diffusion regime.
#'
#'
#' \deqn{ \frac{1}{2} \pi M0/M1 \sum_{j=0}^{2000} (-1)^j (2j + 1) exp(-0.5 \pi^2 (j + 0.5)^2 M0/M1 )}
#'
#' @param diffusion Estimated diffusion coefficient.
#' @param tempTime A vector of time values.
#' @param L The length of the reactor.
#' @param bedPorosity Bed Porosity of the reactor.
#' @return A data frame of the pulse and the standard diffusion values based on time.
#' @examples
#' TAPexperiment = readTAP("data/pumpProbe.tdms")
#' TAPexperiment$AMU40$options$baselineStart = 3.9
#' TAPexperiment$AMU40$options$baselineStart = 4
#' TAPexperiment$AMU40$options$baselineType = "mean"
#'
#' TAPexperiment = baselineCorrect(TAPexperiment, "AMU40")
#' TAPexperiment = SDC(TAPexperiment, "AMU40")
#' plot(TAPexperiment$AMU40SDC$matrix[[5]])
#' @export SDC
SDC = function(TAPexperiment, gasName){

  timeMinPos = which.min(abs(TAPexperiment$parameters$experiment$time - TAPexperiment$parameters$experiment$timeStart))
  timeMaxPos = which.min(abs(TAPexperiment$parameters$experiment$time - TAPexperiment$parameters$experiment$timeEnd))

  timeIndex = timeMinPos:timeMaxPos
  tempTime = TAPexperiment$parameters$experiment$time[timeMinPos:timeMaxPos]

  if(timeMinPos != 1){
    timePre = rep(0,timeMinPos)[-1]
  }
  if(timeMaxPos != TAPexperiment$parameters$experiment$numTime){
    timePost = rep(0, (TAPexperiment$parameters$experiment$numTime - timeMaxPos))
  }

  tempTime = tempTime - min(tempTime)
  tempTime = tempTime[-1]
  iVals = 1:(length(tempTime))
  tempTime = matrix(tempTime, 1, length(tempTime))

  reactorLength = TAPexperiment$parameters$reactor$reactorLength
  nTerms = 2000
  jTerms = 0: nTerms
  negativeIteration = rep(1,(nTerms + 1))
  negativeIteration[jTerms %% 2 == 0] = -1
  negativeIteration = negativeIteration * -1
  negWithJ = matrix(negativeIteration * (2 * jTerms + 1), 1, (nTerms + 1))


  for(i in 1:length(gasName)){
    TAPobj = TAPexperiment[[gasName[i]]]

    if(is.null(TAPobj$moments)){
      TAPexperiment = TAPcodeV2::moments(TAPexperiment, TAPobj$options$name)
      TAPobj = TAPexperiment[[TAPobj$options$name]]
    }

    for(j in 1:TAPexperiment$parameters$experiment$numPulses){
      mRatio = TAPexperiment$parameters$reactor$bedPorosity * reactorLength^2 / (2 * TAPobj$moments[[j]]$diffusion)

      exponentTerms = as.matrix(-.5 * pi^2 * (jTerms + .5)^2 / mRatio)
      exponentWithTime = negWithJ %*% exp(exponentTerms %*% tempTime)
      sdc = pi * as.numeric(exponentWithTime) / (mRatio * 2)
      sdc = c(sdc,0)

      if(timeMinPos != 1){
        sdc = c(timePre, sdc)
      }
      if(timeMaxPos != TAPexperiment$parameters$experiment$numTime){
        sdc = c(sdc, timePost)
      }
      TAPobj$matrix[[j]] = sdc
    }

    TAPobj$options$name = paste0(TAPobj$options$name, "SDC")
    TAPexperiment[[TAPobj$options$name]] = TAPobj
    TAPexperiment = TAPcodeV2::moments(TAPexperiment, TAPobj$options$name)
  }

  return(TAPexperiment)
}
