#' Execute Series of functions
#'
#' @param TAPexperiment Data to be applied
#' @return A collection of TAP objects that describe the experiment performed in the reactor.
#' @examples
#'
#'
#' @export buildTAPScript


buildTAPScript = function(TAPexperiment){

  scriptList = TAPexperiment$reactor$scriptList
  scriptList = scriptList[-1] # may have to delete
  for(i in 1:length(scriptList)){
    firstSplit = strsplit(scriptList[i], "[*]")[[1]]
    functionPerformed = firstSplit[1]

    TAPexperiment = switch(functionPerformed,
                           "Changed Inert and Product of " = {
                             AMUs = gsub("'" , "", strsplit(firstSplit[2], "[,]")[[1]])
                             for(i in 1:length(AMUs)){
                               TAPexperiment[[AMUs[i]]]$options$inert = firstSplit[3]
                               TAPexperiment[[AMUs[i]]]$options$isProduct = firstSplit[4]
                               TAPexperiment[[AMUs[i]]]$options$calibrationCoef = firstSplit[5]
                               TAPexperiment[[AMUs[i]]]$pulses = TAPexperiment[[AMUs[i]]]$pulses * firstSplit[5]
                               TAPexperiment[[AMUs[i]]]$options$diffusionCoef = firstSplit[6]
                             }
                           },
                           "Changed Baseline of " = {
                             AMUs = gsub("'" , "", strsplit(firstSplit[2], "[,]")[[1]])
                             for(i in 1:length(AMUs)){
                               TAPexperiment[[AMUs[i]]]$options$baselineStart = firstSplit[3]
                               TAPexperiment[[AMUs[i]]]$options$baselineEnd = firstSplit[4]
                               TAPexperiment[[AMUs[i]]]$options$baselineType = firstSplit[5]
                               TAPexperiment[[AMUs[i]]]$options$timeStart = firstSplit[6]
                               TAPexperiment[[AMUs[i]]]$options$timeEnd = firstSplit[7]
                             }
                             TAPexperiment = TAPcodeV2::baselineCorrect(TAPexperiment, AMUs)
                           },
                           "Applied SavitskyGolay Smoothing to " = {
                             AMUs = gsub("'" , "", strsplit(firstSplit[2], "[,]")[[1]])
                             for(i in 1:length(AMUs)){
                               TAPexperiment[[AMUs]]$options$savGolSmoothing = as.numeric(firstSplit[3])
                             }
                             TAPexperiment = TAPcodeV2::savGol(TAPexperiment, AMUs)
                           },
                           "Updated Reactor Settings of Inert Zone 1 length to " = {
                             tempParams = TAPexperiment$reactor$reactorParams
                             tempParams$inertZone1Length = firstSplit[2]
                             tempParams$inertZone2Length = firstSplit[3]
                             tempParams$catalystBedLength = firstSplit[4]
                             tempParams$crossSectionalArea = firstSplit[5]
                             tempParams$catalystWeight = firstSplit[6]
                             tempParams$bedPorosity = firstSplit[7]
                             tempParams$molPerM0Inert = firstSplit[8]
                             tempParams$pumpProbeSpacing = firstSplit[9]

                             TAPexperiment$reactor$reactorParams = tempParams
                           },
                           "Calculated Gas Concentration of " = {
                             AMUs = gsub("'" , "", strsplit(firstSplit[2], "[,]")[[1]])
                             for(i in 1:length(AMUs)){
                               if(is.numeric(firstSplit[3])){
                                 TAPexperiment[[AMUs[i]]]$options$yProcSmoothing = firstSplit[3]
                               }else{
                                 TAPexperiment[[AMUs[i]]]$options$yProcSmoothing = NA
                               }
                             }
                             TAPexperiment = TAPcodeV2::gasCon(TAPexperiment, AMUs, method = firstSplit[4])
                           },
                           "Calculated Reaction Rate of " = {
                             AMUs = gsub("'" , "", strsplit(firstSplit[2], "[,]")[[1]])
                             for(i in 1:length(AMUs)){
                               if(is.numeric(firstSplit[3])){
                                 TAPexperiment[[AMUs[i]]]$options$yProcSmoothing = firstSplit[3]
                               }else{
                                 TAPexperiment[[AMUs[i]]]$options$yProcSmoothing = NA
                               }
                             }
                             TAPexperiment = TAPcodeV2::reactionRate(TAPexperiment, AMUs, method = firstSplit[4])
                           },
                           "Applied Pulse Deconvolution to " = {
                             AMUs = gsub("'" , "", strsplit(firstSplit[2], "[,]")[[1]])
                             TAPexperiment = TAPcodeV2::pulseDeconvolution(TAPexperiment, AMUs)
                           },
                           "Applied GrahamInert transform from " = {
                             AMUs = gsub("'" , "", strsplit(firstSplit[2], "[,]")[[1]])
                             TAPexperiment = TAPcodeV2::inertGraham(TAPexperiment, AMUs, firstSplit[3])
                           },
                           "Calculated moments of " = {
                             AMUs = gsub("'" , "", strsplit(firstSplit[2], "[,]")[[1]])
                             TAPexperiment = TAPcodeV2::moments(TAPexperiment, AMUs)
                           },
                           "Calculated reactivities of " = {
                             AMUs = gsub("'" , "", strsplit(firstSplit[2], "[,]")[[1]])
                             productAMUs = gsub("'" , "", strsplit(firstSplit[3], "[,]")[[1]])
                             TAPexperiment = TAPcodeV2::reactivities(TAPexperiment, AMUs, productAMUs)
                           },
                           "Calculated average of " = {
                             AMUs = gsub("'" , "", strsplit(firstSplit[2], "[,]")[[1]])
                             pulseRange = as.numeric(gsub("'" , "", strsplit(firstSplit[3], "[,]")[[1]]))
                             TAPexperiment = pulseAverage(TAPexperiment, AMUs, pulseRange = pulseRange)
                           },
                           "Scaled " = {
                             TAPobj1 = TAPexperiment[[firstSplit[2]]]
                             TAPobj2 = TAPexperiment[[firstSplit[5]]]
                             TAPobjNew = TAPobj1
                             TAPobjNew$pulses = switch(firstSplit[4],
                                                       "-" = (TAPobj1$pulses * firstSplit[3]) - (TAPobj2$pulses * firstSplit[6]),
                                                       "+" = (TAPobj1$pulses * firstSplit[3]) + (TAPobj2$pulses * firstSplit[6]),
                                                       "*" = (TAPobj1$pulses * firstSplit[3]) * (TAPobj2$pulses * firstSplit[6]),
                                                       "/" = (TAPobj1$pulses * firstSplit[3]) / (TAPobj2$pulses * firstSplit[6])
                             )

                             TAPobjNew$temperature = (TAPobj1$temperature + TAPobj2$temperature) / 2
                             TAPobjNew$options$Name = paste0(firstSplit[2], firstSplit[4], firstSplit[5])
                             TAPexperiment[[TAPobjNew$options$Name]] = TAPobjNew
                           },
                           "Calculated SDC of " = {
                             AMUs = gsub("'" , "", strsplit(firstSplit[2], "[,]")[[1]])
                             TAPexperiment = TAPcodeV2::SDC(TAPexperiment, gasName = AMUs)
                           },
                           "Transformed " = {
                             TAPobj = TAPexperiment[[firstSplit[2]]]
                             momentInfo = TAPexperiment[[firstSplit[5]]]$moments[[firstSplit[4]]]
                             TAPobj$pulses = switch(firstSplit[3],
                                                    "-" = sweep(TAPobj$pulses, 2, STATS = momentInfo, FUN = "-"),
                                                    "+" = sweep(TAPobj$pulses, 2, STATS = momentInfo, FUN = "+"),
                                                    "*" = sweep(TAPobj$pulses, 2, STATS = momentInfo, FUN = "*"),
                                                    "/" = sweep(TAPobj$pulses, 2, STATS = momentInfo, FUN = "/")
                             )
                             TAPobj$options$Name = paste0(firstSplit[2], firstSplit[3], firstSplit[5], firstSplit[4])
                             TAPexperiment[[TAPobj$options$Name]] = TAPobj
                           },
                           "Applied Wavelet Smoothing to " = {
                             AMUs = gsub("'" , "", strsplit(firstSplit[2], "[,]")[[1]])
                             for(i in 1:length(AMUs)){
                               TAPexperiment[[AMUs]]$options$waveSmoothing = as.numeric(firstSplit[3])
                             }
                             TAPexperiment = TAPcodeV2::waveTAP(TAPexperiment, AMUs)
                           },
                           "Applied outlier removal to " = {
                             AMUs = gsub("'" , "", strsplit(firstSplit[2], "[,]")[[1]])
                             TAPexperiment = TAPcodeV2::outlierTAP(TAPexperiment, AMUs)
                           }
    )
  }
  return(TAPexperiment)
}





