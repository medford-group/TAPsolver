#' Robust Mode Calculation
#'
#' Roughly based on "On a fast, robust estimator of the mode: Comparisons to other robust estimators with applications"
#'
#' blah
#'
#' @param tempPulse a temporary Pulse
#' @return Blah
#' @examples
#' data("pulseData")
#'
#'
#'@export robustMode

robustMode = function(tempPulse){
  tempIndex = 1:length(tempPulse)
  maxVals = rep(0, 8)
  for(j in 1:length(maxVals)){
    modNum = j
    subIndex = c(1, which(tempIndex%%modNum == 0), length(tempPulse))
    tempMedian = rep(0, (length(subIndex) - 1))
    for(i in 2:length(subIndex)){
      tempMedian[i-1] = median(tempPulse[subIndex[i - 1] : subIndex[i]])
    }
    maxVals[j] = which.max(tempMedian) * j
  }
  result = floor(median(maxVals))
  if(result > length(tempPulse))
    result = round(length(tempPulse)/2)
  return(result)
}
