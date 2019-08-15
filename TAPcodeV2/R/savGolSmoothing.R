#' Savitzky-Golay Smoothing
#'
#' Smooth the \code{pulseVector} by an integer size of the filter length.
#'
#' Savitzky-Golay performs smoothing by using a least squares fit on small subsections of the data to a polynomial curve.
#' As such, the algorithm requires a choice of window sizes of the subsections of the data.  As the window size increases, the pulse becomes more smooth.
#' This algorithm should be used minimally due to the decrease of the magnitude of the pulse.
#'
#' This is a wrapper function on the savgol function from the package pracma.
#' Options on the savgol function allow for the filter order to be managed e.g. 2 = quadratic filter and 4 = quartic (currently using 2).
#' The derivative order may be changed as well e.g. 0 = smoothing and 1 = first derivative and so on.
#' The function is currently setup so that it does quadratic smoothing.
#'
#' @param TAPexperiment A TAP object of pulse values.
#' @param gasName Value
#' @return A vector of smoothed pulse values.
#' @examples
#' TAPexperiment = readTAP("data/pumpProbe.tdms")
#' TAPexperiment$AMU40$options$baselineStart = 3.9
#' TAPexperiment$AMU40$options$baselineStart = 4
#' TAPexperiment$AMU40$options$baselineType = "mean"
#'
#' TAPexperiment = baselineCorrect(TAPexperiment, "AMU40")
#'
#' TAPexperiment$AMU40$options$savGolSmoothing = 15
#' TAPexperiment = savGol(TAPexperiment, "AMU40")
#' plot(TAPexperiment$AMU40$matrix[[5]])
#' @export savGol

#### code For Savitzky-Golay Smoothing
savGol = function(TAPexperiment, gasName){

  for(i in 1:length(gasName)){
    TAPobj = TAPexperiment[[gasName[i]]]

    if(TAPobj$options$savGolSmoothing == 0)
      next

    # have to use odd integer
    TAPobj$options$savGolSmoothing = abs(round(TAPobj$options$savGolSmoothing))

    if((TAPobj$options$savGolSmoothing %% 2) == 0){
      TAPobj$options$savGolSmoothing = TAPobj$options$savGolSmoothing + 1
    }
    TAPobj$matrix = lapply(TAPobj$matrix, FUN = pracma::savgol, fl = TAPobj$options$savGolSmoothing, forder = 2, dorder =0)
    TAPexperiment[[TAPobj$options$name]] = TAPobj
  }

  return(TAPexperiment)
}
