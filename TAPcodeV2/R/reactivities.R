#' Reactivities for TAP
#'
#' Reactivities for TAP
#'
#' blah
#'
#' @param TAPexperiment A set of pre-processing parameters based on the Y-Procedure
#' @param reactantGasName Blah
#' @param productGasName Blah
#' @return Blah
#' @examples
#' TAPexperiment = readTAP("data/pumpProbe.tdms")
#' TAPexperiment$AMU40$options$baselineStart = 3.9
#' TAPexperiment$AMU40$options$baselineStart = 4
#' TAPexperiment$AMU40$options$baselineType = "mean"
#'
#' TAPexperiment = baselineCorrect(TAPexperiment, "AMU40")
#' TAPexperiment = reactivities(TAPexperiment, "AMU40", "AMU40")
#' df = momentsToDataFrame(TAPexperiment, "AMU40")
#' plot(df$R0)
#'
#'@export reactivities
reactivities = function(TAPexperiment, reactantGasName, productGasName){

  # Reactor parameters
  tempParams = TAPexperiment$parameters$reactor
  lengthOfInertZone = tempParams$inertZone1Length + tempParams$inertZone2Length

  reactantObj = TAPexperiment[[reactantGasName]]
  if(is.null(reactantObj$moments)){
    TAPexperiment = TAPcodeV2::moments(TAPexperiment, reactantGasName)
    reactantObj = TAPexperiment[[reactantGasName]]
  }

  inertObj = TAPexperiment[[reactantObj$options$inert]]
  if(is.null(inertObj$moments)){
    TAPexperiment = TAPcodeV2::moments(TAPexperiment, reactantObj$options$inert)
    inertObj = TAPexperiment[[reactantObj$options$inert]]
  }

  inertMoments = momentsToDataFrame(TAPexperiment, reactantObj$options$inert)
  reactMoments = momentsToDataFrame(TAPexperiment, reactantObj$options$name)
  # Calculating the residence time for the inert and reactant
  inertMoments$resTime = tempParams$bedPorosity * lengthOfInertZone^2 / (2 * inertMoments$diffusion)
  reactMoments$diffusion = inertMoments$diffusion * sqrt(inertObj$options$amu / reactantObj$options$amu)
  reactMoments$resTime = tempParams$catalystBedLength * lengthOfInertZone / ( 2 * reactMoments$diffusion)

  # Calculating the reactivities for the reactant
  reactMoments$R0 = -1 / reactMoments$resTime + 1 / (reactMoments$resTime * reactMoments$M0norm)
  reactMoments$R1 = -2 * inertMoments$resTime / (3 * reactMoments$resTime) - inertMoments$resTime /
    (3 * reactMoments$resTime * reactMoments$M0norm) + reactMoments$M1norm / (reactMoments$resTime * reactMoments$M0norm^2)
  reactMoments$R2 = 4 * inertMoments$resTime^2/(45 * reactMoments$resTime) + 7 * inertMoments$resTime^2 / (90 * reactMoments$resTime * reactMoments$M0norm) -
    inertMoments$resTime * reactMoments$M1norm / (3 * reactMoments$resTime * reactMoments$M0norm^2) - reactMoments$M2norm /
    (2 * reactMoments$resTime * reactMoments$M0norm^2) + reactMoments$M1norm^2 / (reactMoments$resTime * reactMoments$M0norm^3)

  reactantObj$moments = dataFrameToMoments(reactMoments, reactantObj$moments)
  TAPexperiment[[reactantGasName]] = reactantObj
  TAPexperiment[[reactantObj$options$inert]] = inertObj

  if(length(productGasName) == 0)
    return(TAPexperiment)

  for(i in 1:length(productGasName)){
    productObj = TAPexperiment[[productGasName[i]]]
    if(is.null(productObj$moments)){
      TAPexperiment = TAPcodeV2::moments(TAPexperiment, productGasName[i])
      productObj = TAPexperiment[[productGasName[i]]]
    }

    productMoments = momentsToDataFrame(TAPexperiment, productObj$options$name)

    inertMoments$diffusion * sqrt(inertObj$options$amu / productObj$options$amu)
    productMoments$resTime = tempParams$catalystBedLength * lengthOfInertZone / ( 2 * productMoments$diffusion)
    productMoments$R0 = productMoments$M0norm / (productMoments$resTime * reactMoments$M0norm)
    productMoments$R1 = productMoments$R0 * (inertMoments$resTime / 12 * (8 * reactMoments$M0norm + 3 + 9 * reactMoments$diffusion / productMoments$diffusion) +
                                     productMoments$resTime * reactMoments$M0norm * reactMoments$R1 - productMoments$M1M0)
    productMoments$R2 = productMoments$R0 / 2 * (productMoments$M2M0 - 19 * reactMoments$diffusion^2 * inertMoments$resTime / (16 * productMoments$diffusion) *
                                         (productMoments$R0 * (inertMoments$resTime * (3 + 8 * reactMoments$M0norm) + 12 * reactMoments$M0norm *
                                                           reactMoments$R1 * productMoments$resTime) - 12 * productMoments$R1) +
                                         productMoments$R1 / (6 * productMoments$R0) * ((3 + 8 * reactMoments$M0norm) * inertMoments$resTime + 12 * reactMoments$M0norm *
                                                                                reactMoments$R1 * productMoments$resTime) - inertMoments$resTime^2 *
                                         (5 / 48 + reactMoments$M0norm / 45 * (23 + 40 * reactMoments$M0norm)) -
                                         reactMoments$M0norm / 6 * (3 + 16 * reactMoments$M0norm) * reactMoments$R1 * productMoments$resTime * inertMoments$resTime -
                                         2 * reactMoments$M0norm * productMoments$resTime * (reactMoments$M0norm * reactMoments$R1^2 * productMoments$resTime - reactMoments$R2))
    productObj$moments = dataFrameToMoments(productMoments, productObj$moments)
    TAPexperiment[[productGasName[i]]] = productObj
  }

  return(TAPexperiment)
}
