#' TAP Flux Constructor
#'
#' Generates prefilled Flux parameters for options.  Does not include prefilled flux
#' @param none No value necessary
#' @return Prefilled Flux parameters
#' @examples
#' TAPstruct = constructFlux()
#' names(TAPstruct)
#' @export constructFlux

constructFlux = function(none = NULL){
  TAPobj = list(
    "matrix" = NULL, # A list of flux values per the number of pulses
    "options" = NULL, # A list of optional parameters for
    "moments" = NULL # A list of moments per the number of pulses
  )

  TAPobj$options = list(
    "name" = "AMU40", # the name associated with the gas species
    "amu" = 40, # the amu
    "gain" = 9, # gain setting in measurement
    "baselineStart" = 0, # start of the baseline
    "baselineEnd" = 2, # end of the baseline
    "baselineType" = "mean", # baseline type
    "inert" = "AMU40", # associated inert
    "isProduct" = TRUE, # if the gas is a product or inert. Used in calculating rates
    "calibrationCoef" = 1, # the mass spec calibration coefficient
    "savGolSmoothing" = 0, # the amount of smoothing to be used in SavitzkyGolay
    "yProcSmoothing" = 0, # Amount of smoothing for Y-Procedure.  If zero, smoothing parameter will be fit
    "waveSmoothing" = 0 # The amount of smoothing using wavelet thresholding
  )
  return(TAPobj)
}
