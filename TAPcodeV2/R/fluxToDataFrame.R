#' Flux to Data frame
#'
#' Helper function for subsetting data
#'
#'
#' @param Array An array.
#' @param index An index of values.
#' @return The TAPexperiment with applied function.
#' @examples
#' TAPexperiment = readTAP("data/pumpProbe.tdms")
#'
#' TAPexperiment$AMU40$options$baselineStart = 3.9
#' TAPexperiment$AMU40$options$baselineStart = 4
#' TAPexperiment$AMU40$options$baselineType = "mean"
#'
#' TAPexperiment = baselineCorrect(TAPexperiment, "AMU40")
#'
#' # plot a baseline corrected pulse
#' plot(TAPexperiment$AMU40$matrix[, 5])
#' @export fluxToDataFrame
#'

fluxToDataFrame = function(TAPexperiment, gasName){
  result = as.data.frame((matrix(unlist(unlist(TAPexperiment[[gasName]]$matrix,F,F)), ncol = TAPexperiment$parameters$experiment$numPulses, byrow = F )))
  names(result) = names(TAPexperiment[[gasName]]$matrix)
  return(result)
}
