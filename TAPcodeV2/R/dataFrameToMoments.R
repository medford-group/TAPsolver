#' Data Frame to Moments
#'
#' Helper function for subsetting data
#'
#'
#' @param df An array.
#' @param momentsObj An index of values.
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
#' @export dataFrameToMoments
#'

dataFrameToMoments = function(df, momentsObj){
  for(i in 1:length(momentsObj)){
    momentsObj[[i]] = as.list(df[i,])
  }
  return(momentsObj)
}



