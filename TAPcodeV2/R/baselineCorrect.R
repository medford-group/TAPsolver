#' Baseline correction on a TAP object
#'
#' Correct the baseline of a pulse towards zero
#'
#' The TAP pulse mass spec results may sometimes not reach zero appropriately when no product is produced or when all the reactant is consumed.
#' In this case, it is nessecary to perform a baseline correction of the respective pulse.  This is slightly different from typical baseline correction because it is only necessary
#' to shift the pulse vertically such that the baseline portion is at zero.  This method requires a user to know at what time segment the baseline is not changing, i.e., the slope of the pulse
#' remains constant or flat.  Once baseline has been determined, this function will shift the whole pulse vertically by subtracting a choice of either the min or mean of the baseline time segment.
#'
#'
#' @param TAPexperiment A collection of TAP objects.
#' @param gasName The gas in which to apply the function on.
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
#' plot(TAPexperiment$AMU40$matrix[[5]])
#' @export baselineCorrect
#'

baselineCorrect = function(TAPexperiment, gasName){

  baselineHelper = function(x, helperOptions){
    result = switch(helperOptions$Type,
                    "mean" = x - mean(x[helperOptions$index]),
                    "min" = x - min(x[helperOptions$index]),
                    "none" = x
                    )
    return(result)
  }
  for(i in 1:length(gasName)){
    if(TAPexperiment[[gasName[i]]]$options$baselineType == "none")
      next
    baselineMinPos = which.min(abs(TAPexperiment$parameters$experiment$time - TAPexperiment[[gasName[i]]]$options$baselineStart))
    baselineMaxPos = which.min(abs(TAPexperiment$parameters$experiment$time - TAPexperiment[[gasName[i]]]$options$baselineEnd))
    baselineIndex = baselineMinPos : baselineMaxPos
    baselineOptions = list("Type" = TAPexperiment[[gasName[i]]]$options$baselineType, "index" = baselineIndex)

    #cl = snow::makeCluster(TAPexperiment$parameters$experiment$numCores)
    #TAPexperiment[[gasName[i]]]$matrix = snow::parLapply(cl, TAPexperiment[[gasName[i]]]$matrix, baselineHelper, baselineOptions = baselineOptions)
    TAPexperiment[[gasName[i]]]$matrix = lapply(TAPexperiment[[gasName[i]]]$matrix, baselineHelper, helperOptions = baselineOptions)
    #snow::stopCluster(cl)
    #TAPexperiment[[gasName[i]]] = TAPexperiment[[gasName[i]]]
  }
  return(TAPexperiment)
}

