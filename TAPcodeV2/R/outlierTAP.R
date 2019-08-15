#' Outlier Removal
#'
#' Partition the pulse with an exponential decay from the mode.  Any points above or below 2 standard deviations are set to the median of the partition.
#'
#' @param TAPexperiment A TAP object of pulse values.
#' @param gasName Value
#' @return A vector of smoothed pulse values.
#' @examples
#' TAPexperiment = readTAP("data/pumpProbe.tdms")
#'
#' TAPexperiment$AMU40$options$baselineStart = 3.9
#' TAPexperiment$AMU40$options$baselineStart = 4
#' TAPexperiment$AMU40$options$baselineType = "mean"
#'
#' TAPexperiment = baselineCorrect(TAPexperiment, "AMU40")
#' TAPexperiment = outlierTAP(TAPexperiment, "AMU40")
#'
#' plot(TAPexperiment$AMU40$matrix[[5]])
#' @export outlierTAP

outlierTAP = function(TAPexperiment, gasName){

  helperOptions = list()
  helperOptions$time = TAPexperiment$parameters$experiment$time
  helperOptions$maxTime = max(helperOptions$time)

  outlierHelper = function(x, helperOptions){
    tempMode = helperOptions$time[robustMode(x)]

    testSpacing = c(exp(seq(0,log(max(helperOptions$time) + 1),by = .05)) - 1 + tempMode, helperOptions$maxTime)
    negativeTestSpacing = c(0,sort(abs(exp(seq(0,log(tempMode + 1),by = .05)) - 1 - tempMode)))
    testSpacing = c(negativeTestSpacing[-length(negativeTestSpacing)], testSpacing)

    testIndex = rep(0,length(testSpacing))
    for(k in 1:length(testIndex)){
      testIndex[k] = which.min(abs(helperOptions$time - testSpacing[k]))
    }

    withoutOutliers = x
    for(k in 1:(length(testIndex)-1)){
      subPulse = x[testIndex[k]:testIndex[k+1]]
      subMedian = median(subPulse)
      subSD = sd(subPulse)
      subPulse[abs(subPulse) > (subMedian + 1.5 * subSD)] = subMedian
      withoutOutliers[testIndex[k]:testIndex[k+1]] = subPulse
    }
    return(withoutOutliers)
  }

  for(i in 1:length(gasName)){
    TAPobj = TAPexperiment[[gasName[i]]]
    TAPobj$matrix = lapply(TAPobj$matrix, outlierHelper, helperOptions = helperOptions)

    TAPexperiment[[TAPobj$options$name]] = TAPobj
  }
  return(TAPexperiment)
}
