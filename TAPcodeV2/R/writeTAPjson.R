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
#' writeTAPjson(TAPexperiment, "data/pumpProbe.json")
#' @export writeTAPjson

writeTAPjson = function(TAPexperiment, dataPath){
  tempJson = RJSONIO::toJSON(TAPexperiment)
  write(tempJson, dataPath)
  #jsonlite::write_json(TAPexperiment, dataPath)
}
