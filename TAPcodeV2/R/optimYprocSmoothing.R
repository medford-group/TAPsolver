#' Optimization of the Y-Procedure Smoothing value
#'
#' Optimization of the Y-Procedure Smoothing value
#'
#' blah
#' TAPgasConVars = list("tempPulse" = tempPulse, "omega" = omega, "timeStep" = timeStep, "gasScalar" = gasScalar, "gamma3" = gamma3, "gasConScalar" = gasConScalar)
#'
#'
#' @param TAPgasConVars A set of pre-processing parameters based on the Y-Procedure
#' @return A vector of moments based on the pulse.
#' @examples
#' data("pulseData")
#'
#'
#'@export optimYprocSmoothing

optimYprocSmoothing = function(TAPgasConVars){

  # Needed parameters resulting from the Y-Procedure pre-processing
  # TAPgasConVars = list("tempPulse" = tempPulse, "omega" = omega, "timeStep" = timeStep, "gasScalar" = gasScalar, "gamma3" = gamma3, "gasConScalar" = gasConScalar)
  optimizeSmoothing = function(smoothing, TAPgasConVars){
    smoothScalar = exp(-TAPgasConVars$omega^2 * TAPgasConVars$timeStep^2 * smoothing^2 / 2) #TAPobj$options$yProcSmoothing
    gasScalar = TAPgasConVars$gasScalar * smoothScalar / TAPgasConVars$gamma3
    gasResult = Re(fft(fft(TAPgasConVars$tempPulse) * gasScalar, inverse = T))
    gasResult = gasResult / pracma::trapz(gasResult)
    gasFlux = TAPgasConVars$tempPulse / pracma::trapz(TAPgasConVars$tempPulse)

    varResult = var(gasResult[round(.95 * length(TAPgasConVars$tempPulse)):length(TAPgasConVars$tempPulse)])
    varFlux = var(gasFlux[round(.95 * length(TAPgasConVars$tempPulse)):length(TAPgasConVars$tempPulse)])
    return(abs(1 - varResult / varFlux))
  }

  result = optim(1, fn = optimizeSmoothing,lower = 0, upper = 6, TAPgasConVars = TAPgasConVars, method = "L-BFGS-B")
  return(result$par)
}
