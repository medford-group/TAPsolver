#' Pulse Average
#'
#' Pulse Average
#'
#' blah
#'
#' @param TAPexperiment A set of pre-processing parameters based on the Y-Procedure
#' @param gasName Blah
#' @return Blah
#' @examples
#' TAPexperiment = readTAP("data/pumpProbe.tdms")
#'
#' TAPexperiment$AMU40$options$baselineStart = 3.9
#' TAPexperiment$AMU40$options$baselineStart = 4
#' TAPexperiment$AMU40$options$baselineType = "mean"
#'
#' TAPexperiment = baselineCorrect(TAPexperiment, "AMU40")
#' TAPexperiment = pulseAverage(TAPexperiment, "AMU40", 1:10)
#'
#' plot(TAPexperiment$AMU40avg$matrix[[1]])
#'
#'@export pulseAverage
pulseAverage = function(TAPexperiment, gasName, pulseRange = NULL){

  for(i in 1:length(gasName)){
    TAPobj = TAPexperiment[[gasName[i]]]

    if(is.null(TAPobj$moments)){
      TAPexperiment = moments(TAPexperiment, gasName[i])
      TAPobj = TAPexperiment[[gasName[i]]]
    }


    if(is.null(pulseRange))
      pulseRange = 1:TAPexperiment$parameters$experiment$numPulses

    tempMean = fluxToDataFrame(TAPexperiment, gasName[i])[,pulseRange]
    TAPobj$matrix = list("avg" = apply(tempMean, 1, mean))

    tempMoments = momentsToDataFrame(TAPexperiment, gasName[i])[pulseRange,]
    TAPobj$moments = as.list(apply(tempMoments, 2, mean))

    TAPobj$options$name = paste0(gasName[i], "avg")
    TAPexperiment[[TAPobj$options$name]] = TAPobj
  }

  return(TAPexperiment)
}








