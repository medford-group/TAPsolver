#' Set Diffusion Coefficient based on moments
#'
#' Blah
#'
#' blah
#'
#' @param TAPexperiment A set of pre-processing parameters based on the Y-Procedure
#' @param gasName Blah
#' @return Blah
#' @examples
#' data("pulseData")
#'
#'
#'@export setDiffusion

setDiffusion = function(TAPexperiment, gasName){
  TAPobj = TAPexperiment[[gasName]]
  inertObj = TAPexperiment[[TAPobj$options$inert]]
  reactorParams = TAPexperiment$parameters$reactor

  reactorLength = reactorParams$reactorLength

  if(is.null(inertObj$moments)){
    TAPexperiment[[TAPobj$options$inert]] = inertObj
    TAPexperiment = TAPcodeV2::moments(TAPexperiment, TAPobj$options$inert)
    inertObj = TAPexperiment[[TAPobj$options$inert]]
  }

  inertMoments = momentsToDataFrame(inertObj$moments)
  inertObj$options$diffusionCoef = mean(reactorParams$bedPorosity * reactorLength^2 * inertMoments$M0 / (2 * inertMoments$M1))
  TAPobj$options$diffusionCoef = inertObj$options$diffusionCoef * sqrt(inertObj$options$amu / TAPobj$options$amu)

  TAPexperiment[[gasName]] = TAPobj
  TAPexperiment[[TAPobj$options$inert]]
  return(TAPexperiment)
}
