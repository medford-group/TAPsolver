#' Pulse Normalization
#'
#' This function normalizes the plot based on \code{normalization}.  Typically, this function is used to compare the inert gas with the Standard Diffusion Curve.
#' Options for normalization include "none", "height", and "area".
#'
#' For height normalization:
#' \deqn{\frac{pulse(t)}{max(pulse(t))}}
#'
#' For area normalization:
#' \deqn{\frac{pulse(t)}{\int pulse(t)dt}}
#'
#'
#'
#' @param pulseVector A vector of pulse values.
#' @param timeVector A vector of time values.
#' @param normalization A string that indicates the normalization type.  Options: "none", "height", "area".
#' @return A normalized pulse vector.
#' @examples
#' data("pulseData")
#'
#' pulse = pulseData$`32AMU`$pulses[,33]
#' timeVector = pulseData$`32AMU`$time
#' pulseNorm(pulse, timeVector, "area")
#'
#'
#'
#' @export pulseNorm
pulseNorm = function(pulseVector, timeVector, normalization = "none"){

  if(is.null(dim(pulseVector))){
    result = switch(normalization,
                    "none" = pulseVector,
                    "height" = pulseVector/max(pulseVector, na.rm = T),
                    "area" = pulseVector/TAPcode::moments(pulseVector, timeVector,0))
  }else{
    if(normalization == "none"){
      divisionNum = 1
    }else if(normalization == "height"){
      divisionNum = rep(1, dim(pulseVector)[2])
      for(i in 1:dim(pulseVector)[2]){
        divisionNum[i] = max(pulseVector[,i], na.rm = T)
      }
    }else if(normalization == "area"){
      divisionNum = TAPcode::moments(pulseVector, timeVector, 0)
    }
    result = sweep(pulseVector, 2, divisionNum, FUN = "/")
  }

  return(result)
}

