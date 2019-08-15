#' Write Excel TAP file
#'
#' Write an Excel .xlsx of the original TDMS file.
#'
#' This function is based on R package rio.
#'
#' @param dataPath This is a string path to an .xlsx file.
#' @return A collection of TAP objects that describe the experiment performed in the reactor.
#' @examples
#' TAPexperiment = readTAP("data/pumpProbe.tdms")
#' TAPexperiment$AMU40$options$baselineStart = 3.9
#' TAPexperiment$AMU40$options$baselineStart = 4
#' TAPexperiment$AMU40$options$baselineType = "mean"
#'
#' TAPexperiment = baselineCorrect(TAPexperiment, "AMU40")
#' TAPexperiment = moments(TAPexperiment, "AMU40")
#' writeTAP(TAPexperiment, "data/pumpProbe.xlsx")
#' @export writeTAP

writeTAP = function(TAPexperiment, dataPath){

  sheetList = list()
  TAPnames = names(TAPexperiment)
  # reactorParams is manipulated and needs to be adjusted in the DF

  reactorDF = data.frame("Key" = names(TAPexperiment$parameters$reactor),
                         "Value" = unlist(TAPexperiment$parameters$reactor))
  sheetList[["reactor"]] = reactorDF

  experimentDF = data.frame("Key" = names(TAPexperiment$parameters$experiment),
                            "Value" = names(TAPexperiment$parameters$experiment))
  experimentDF$Value = as.character(experimentDF$Value)
  for(i in 1:dim(experimentDF)[1]){
    tempStr = TAPexperiment$parameters$experiment[[i]]
    if(length(tempStr) > 1){
      tempStr = paste(tempStr, collapse = ",")
    }
    experimentDF$Value[i] = tempStr
  }
  sheetList[["experiment"]] = experimentDF

  simulationDF = data.frame("Key" = names(TAPexperiment$parameters$simulation),
                            "Value" = names(TAPexperiment$parameters$simulation))
  simulationDF$Value = as.character(simulationDF$Value)
  for(i in 1:dim(simulationDF)[1]){
    tempStr = as.character(TAPexperiment$parameters$simulation[[i]])
    if(length(tempStr) > 1){
      tempStr = paste(tempStr, collapse = ",")
    }
    simulationDF$Value[i] = tempStr
  }
  sheetList[["simulation"]] = simulationDF

  sheetList[["analysis"]] = TAPexperiment$parameters$analysis


  for(i in 1:length(TAPexperiment)){
    if(TAPnames[i] == "parameters")
      next
    TAPobj = TAPexperiment[[i]]
    pulseNames = names(TAPobj$matrix)
    tempParamNames = names(TAPobj$options)
    tempParamValue = as.character(unlist(TAPobj$options,F,F))
    if(length(tempParamNames) < TAPexperiment$parameters$experiment$numTime){
      tempParamNames = c(tempParamNames, rep("", TAPexperiment$parameters$experiment$numTime - length(tempParamNames)))
      tempParamValue = c(tempParamValue, rep("", TAPexperiment$parameters$experiment$numTime - length(tempParamValue)))
    }
    tempDF = cbind(tempParamNames, tempParamValue, TAPexperiment$parameters$experiment$time, fluxToDataFrame(TAPexperiment, TAPnames[i]))
    names(tempDF) = c("Key", "Value", "Time", pulseNames)
    sheetList[[TAPnames[i]]] = tempDF
  }
  rio::export(sheetList, dataPath)
}

