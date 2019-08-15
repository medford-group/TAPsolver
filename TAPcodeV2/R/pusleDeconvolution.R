#' Pulse Deconvolution
#'
#' Pulse Deconvolution
#'
#' blah
#'
#' @param TAPexperiment A set of collected gas species.
#' @param gasName Name of the Gas to be processed.
#' @return Blah
#' @examples
#' TAPexperiment = readTAP("data/pumpProbe.tdms")
#'
#' TAPexperiment$AMU28$options$baselineStart = 3.9
#' TAPexperiment$AMU28$options$baselineStart = 4
#' TAPexperiment$AMU28$options$baselineType = "mean"
#'
#' TAPexperiment = baselineCorrect(TAPexperiment, "AMU28")
#'
#' TAPexperiment$parameters$experiment$pumpProbeSpacing = 0.1
#'
#' TAPexperiment = pulseDeconvolution(TAPexperiment, "AMU28")
#'
#' # plot a baseline corrected pulse
#' plot(TAPexperiment$AMU28$matrix[[5]])
#' lines(TAPexperiment$AMU28pre$matrix[[5]], col = "green")
#' lines(TAPexperiment$AMU28post$matrix[[5]], col = "red")
#'
#'@export pulseDeconvolution
pulseDeconvolution = function(TAPexperiment, gasName){


  for(i in 1:length(gasName)){
    TAPobj = TAPexperiment[[gasName[i]]]

    if(is.null(TAPexperiment$parameters$experiment$pumpProbeSpacing))
      stop("Must include a time for the start of the second pulse.")

    pulseDelayPosition = 1:(which.min(abs(TAPexperiment$parameters$experiment$time - TAPexperiment$parameters$experiment$pumpProbeSpacing)))

    preTAPobj = TAPobj
    preTAPobj$options$timeEnd = TAPexperiment$parameters$experiment$pumpProbeSpacing
    preTAPobj$options$Name = paste0(TAPobj$options$name, "pre")
    TAPexperiment[[paste0(TAPobj$options$name, "pre")]] = preTAPobj
    TAPexperiment = TAPcodeV2::moments(TAPexperiment, preTAPobj$options$name)
    preTAPobj = TAPexperiment[[preTAPobj$options$name]]

    preTime = TAPexperiment$parameters$experiment$time[pulseDelayPosition]

    postPulses = TAPobj$matrix
    prePulses = TAPobj$matrix


    for(j in 1:length(preTAPobj$matrix)){
      normPreFlux = preTAPobj$matrix[[j]][pulseDelayPosition] / preTAPobj$moments[[j]]$M0

      # tempMean = mean(normPreFlux * preTAPobj$time * max(preTAPobj$time)) * pi / 4
      # if( mean(normPreFlux * preTAPobj$time * max(preTAPobj$time)) < 0){
      #   tempMean = preTAPobj$time[which.max(normPreFlux)] * 3 * pi /4
      # }

      tempMean = preTime[which.max(normPreFlux)] * 3 * pi /4
      basicShape = 3/2

      preFluxTail = density(rgamma(1000000, scale =  tempMean * 2/3, shape = basicShape), n = (TAPexperiment$parameters$experiment$numTime - length(pulseDelayPosition)), from = max(preTime), to = TAPexperiment$parameters$experiment$timeEnd)$y
      prePulses[[j]] = c(normPreFlux, preFluxTail* normPreFlux[length(normPreFlux)] / (preFluxTail[1] + 1e-10)) * preTAPobj$moments[[j]]$M0
      postPulses[[j]] = TAPobj$matrix[[j]] - prePulses[[j]]
    }

    preTAPobj$matrix = prePulses
    preTAPobj$options$name = paste0(TAPobj$options$name, "pre")
    #preTAPobj$time = TAPobj$time
    postTAPobj = TAPobj
    postTAPobj$matrix = postPulses
    postTAPobj$options$name = paste0(TAPobj$options$name, "post")

    TAPexperiment[[preTAPobj$options$name]] = preTAPobj
    TAPexperiment[[postTAPobj$options$name]] = postTAPobj

    # if(!is.null(TAPobj$moments)){
    #   TAPexperiment = TAPcodeV2::moments(TAPexperiment, preTAPobj$options$name)
    #   TAPexperiment = TAPcodeV2::moments(TAPexperiment, postTAPobj$options$name)
    # }
  }
  return(TAPexperiment)
}



