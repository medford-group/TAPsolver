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
#' writeMoments(TAPexperiment, "data/pumpProbeMoments.xlsx")
#' @export writeMoments

writeMoments = function(TAPexperiment, dataPath){
  sheetList = list()
  TAPnames = names(TAPexperiment)

  # reactorParams is manipulated and needs to be adjusted in the DF
  #TAPexperiment$parameters$reactor = unlist(TAPexperiment$reactor$reactorParams,F,F)
  sheetList[["parameters"]] = as.data.frame(TAPexperiment$parameters$reactor)
  sheetList[["secondaryInfo"]] = ""
  sheetList[["metaData"]] = ""


  for(i in 1:length(TAPexperiment)){
    if(TAPnames[i] == "parameters")
      next
    TAPobj = TAPexperiment[[i]]
    if(is.null(TAPobj$moments))
      next
    sheetList[[TAPnames[i]]] = momentsToDataFrame(TAPexperiment, names(TAPexperiment)[i])
  }
  rio::export(sheetList, dataPath)

}

