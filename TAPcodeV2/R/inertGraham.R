#' Inert Graham Calculation
#'
#' Scale inert flux by Graham's Law
#'
#' This funciton is for manipulation of the inert flux to appear as if it has a different mass.  This allows the inert gas to behave as the reactant gas
#' without any interaction on the catalyst.  This function scales the magnitude and time by the square root of the ratio of masses.  For example, if the
#' inert gas was Argon (AMU of 40) and the reactant was Oxygen (AMU 32), the inert will be scaled by \eqn{\sqrt{40/32}}.
#'
#'
#' @param TAPobj The inert gas flux
#' @param newMass Transformed mass
#' @return A matrix or vector of transformed inert flux
#'
#' @examples
#' TAPexperiment = readTAP("data/pumpProbe.tdms")
#' TAPexperiment$AMU40$options$baselineStart = 3.9
#' TAPexperiment$AMU40$options$baselineStart = 4
#' TAPexperiment$AMU40$options$baselineType = "mean"
#'
#' TAPexperiment = baselineCorrect(TAPexperiment, "AMU40")
#' TAPexperiment = inertGraham(TAPexperiment, "AMU40", 20)
#' plot(TAPexperiment[["AMU40"]]$matrix[[5]])
#' lines(TAPexperiment[["AMU40to20"]]$matrix[[5]])
#'
#' @export inertGraham
#'
inertGraham = function(TAPexperiment, gasName, newMass){

  inertHelper = function(x, helperOptions){
    tempFun = approxfun(helperOptions$time, x)

    if(helperOptions$grahamsConstant >= 1){
      x = helperOptions$grahamsConstant * c(tempFun( seq(helperOptions$minTime, helperOptions$maxTime,length.out = helperOptions$roundPt) ), rep(tempFun(helperOptions$maxTime), length(helperOptions$time) - helperOptions$roundPt))
    }else{
      x = tempFun( seq(helperOptions$minTime, helperOptions$maxTime,length.out = length(helperOptions$time)))  *  helperOptions$grahamsConstant
    }
    return(x)
  }


  for(i in 1:length(gasName)){
    TAPobj = TAPexperiment[[gasName[i]]]
    minTime = min(TAPexperiment$parameters$experiment$time)
    maxTime = max(TAPexperiment$parameters$experiment$time)
    grahamsConstant = sqrt(TAPexperiment[[gasName[i]]]$options$amu / newMass)
    roundPt = round(length(TAPexperiment$parameters$experiment$time) / grahamsConstant)
    helperOptions = list("minTime"=minTime, "maxTime" = maxTime,
                         "roundPt" = roundPt, "grahamsConstant" = grahamsConstant,
                         "time" = TAPexperiment$parameters$experiment$time)

    TAPobj$matrix = lapply(TAPobj$matrix, FUN = inertHelper, helperOptions = helperOptions)

    TAPobj$options$name = paste0(gasName, "to", newMass)
    TAPexperiment[[TAPobj$options$name]] = TAPobj
  }
  return(TAPexperiment)
}
