#' Init Y-Procedure
#'
#' Optimization of the Y-Procedure Smoothing value
#'
#' blah
#'
#' @param TAPexperiment A set of pre-processing parameters based on the Y-Procedure
#' @param gasName Blah
#' @param method Blah
#' @return Blah
#' @examples
#' TAPexperiment = readTAP("data/pumpProbe.tdms")
#' TAPexperiment$AMU40$options$baselineStart = 3.9
#' TAPexperiment$AMU40$options$baselineStart = 4
#' TAPexperiment$AMU40$options$baselineType = "mean"
#'
#' TAPexperiment = baselineCorrect(TAPexperiment, "AMU40")
#' TAPexperiment$AMU40$options$inert = "AMU40"
#' TAPexperiment = initYproc(TAPexperiment, "AMU40")
#' plot(TAPexperiment$AMU40$gasScalar)
#'
#'@export initYproc

initYproc = function(TAPexperiment, gasName, method = "y"){

  for(i in 1:length(gasName)){
    reactorParams = TAPexperiment$parameters$reactor
    TAPobj = TAPexperiment[[gasName[i]]]
    inertObj = TAPexperiment[[TAPobj$options$inert]]

    timeMinPos = which.min(abs(TAPexperiment$parameters$experiment$time - TAPexperiment$parameters$experiment$timeStart))
    timeMaxPos = which.min(abs(TAPexperiment$parameters$experiment$time - TAPexperiment$parameters$experiment$timeEnd))
    timeIndex = timeMinPos:timeMaxPos

    currentTime = TAPexperiment$parameters$experiment$time[timeIndex]
    # Ensure that time starts at 0
    if(currentTime[1] == 0){
      timeStep = currentTime[2]
    }else{
      currentTime = currentTime - min(currentTime)
      timeStep = currentTime[2]
    }
    timeLen = length(currentTime)


    #TAPobj$reactorLength = reactorLength

    #TAPexperiment[[TAPobj$options$inert]] = inertObj
    #TAPexperiment[[TAPobj$options$Name]] = TAPobj

    if(is.null(TAPobj$moments)){
      TAPexperiment = moments(TAPexperiment, TAPobj$options$name)
      TAPobj = TAPexperiment[[TAPobj$options$name]]
    }

    if(is.null(inertObj$moments)){
      TAPexperiment = moments(TAPexperiment, TAPobj$options$inert)
      inertObj = TAPexperiment[[TAPobj$options$inert]]
    }
    inertMoments = momentsToDataFrame(TAPexperiment, TAPobj$options$inert)

    inertDiffusion = inertMoments$diffusion * sqrt(inertObj$options$amu / TAPobj$options$amu)

    TAPobj$gasConScalar = reactorParams$catalystBedLength / (inertDiffusion * reactorParams$crossSectionalArea * reactorParams$reactorLength) * inertMoments$M0 * reactorParams$molPerM0Inert / timeLen
    TAPobj$reactionRateScalar = inertMoments$M0 * reactorParams$molPerM0Inert / timeLen / reactorParams$catalystWeight

    if(method == "y"){
      k1 = seq(0, (timeLen / 2))
      k2 = sort(seq(-timeLen / 2 + 1)[-c(1, 2)])
      k = c(k1,k2)
      omega = 2 * pi * k / (timeLen * timeStep)
      omega[1] = 1.0e-10
      TAPobj$omega = omega

      gamma1 = mean(inertDiffusion) / reactorParams$inertZone1Length
      gamma3 = mean(inertDiffusion) / reactorParams$inertZone2Length
      tau1 = reactorParams$bedPorosity * reactorParams$inertZone1Length^2 / mean(inertDiffusion)
      tau3 = reactorParams$bedPorosity * reactorParams$inertZone2Length^2 / mean(inertDiffusion)
      iwt1 = sqrt(1i * omega * tau1)
      iwt3 = sqrt(1i * omega * tau3)
      gasScalar = sinh(iwt3) / iwt3
      gasScalar[1] = 1

      if(!is.numeric(TAPobj$options$yProcSmoothing)){
        TAPgasConVars = list("tempPulse" = TAPobj$matrix[[5]], "omega" = omega, "timeStep" = timeStep, "gasScalar" = gasScalar, "gamma3" = gamma3, "gasConScalar" = TAPobj$gasConScalar[5])
        TAPobj$options$yProcSmoothing = TAPcodeV2::optimYprocSmoothing(TAPgasConVars)
      }
      smoothScalar = exp(-TAPobj$omega^2 * timeStep^2 * TAPobj$options$yProcSmoothing^2 / 2)

      TAPobj$gasScalar = gasScalar * smoothScalar / gamma3
      TAPobj$rateScalar = (cosh(iwt3) + sqrt(tau1 * gamma1^2 / (tau3 * gamma3^2)) * sinh(iwt1) * sinh(iwt3) / cosh(iwt1)) * smoothScalar
    }

    TAPexperiment[[TAPobj$options$inert]] = inertObj
    TAPexperiment[[TAPobj$options$name]] = TAPobj
  }


  return(TAPexperiment)
}
