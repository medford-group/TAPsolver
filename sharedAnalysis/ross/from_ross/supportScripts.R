library(TAPcodeV2)
library(hierNet)
library(glmnet)
library(corrr)

# build TAP object

solverToObj = function(solverDF){
  TAPexperiment = list()
  
  
  reactorDF = data.frame("reactorParameters" = c("inertZone1Length", "inertZone2Length",
                                                 "catalystBedLength", "crossSectionalArea",
                                                 "catalystWeight", "bedPorosity",
                                                 "molPerM0Inert", "pumpProbeSpacing"),
                         "value" = c(0.014, 0.014, 1e-3, 1.8e-05, 1, 0.53, 1.63E-09, 0),
                         "unit" = c("m", "m", "m", "m^2", "g", "none", "mol", "s"))
  reactorParams = as.list(reactorDF$value)
  names(reactorParams) = reactorDF$reactorParameters
  TAPexperiment[["reactor"]] = list("reactorParams" = reactorParams,"reactorDF" = reactorDF)
  
  # CO2
  currentList = list()
  # Assign AMU
  currentList$AMU = 44
  currentList$time = solverDF$timeVector
  currentList$temperature = 30
  currentList$pulses = data.frame("1" = solverDF$co2, "2" = solverDF$co2)
  
  tempOptions = list()
  tempOptions[["Name"]] = paste0("AMU", currentList$AMU)
  tempOptions[["AMU"]] = currentList$AMU
  tempOptions[["baselineStart"]] = 0
  tempOptions[["baselineEnd"]] = .5
  tempOptions[["baselineType"]] = "mean"
  tempOptions[["timeStart"]] = 0
  tempOptions[["timeEnd"]] = 6
  tempOptions[["inert"]] = "AMU40"
  tempOptions[["isProduct"]] = T
  tempOptions[["calibrationCoef"]] = 1
  tempOptions[["diffusionCoef"]] = 0
  tempOptions[["savGolSmoothing"]] = 0
  tempOptions[["yProcSmoothing"]] = 0
  tempOptions[["waveSmoothing"]] = 0
  currentList$options = tempOptions
  TAPexperiment[[currentList$options$Name]] = currentList
  
  # CO
  currentList = list()
  # Assign AMU
  currentList$AMU = 28
  currentList$time = solverDF$timeVector
  currentList$temperature = 30
  currentList$pulses = data.frame("1" = solverDF$co, "2" = solverDF$co)
  
  tempOptions = list()
  tempOptions[["Name"]] = paste0("AMU", currentList$AMU)
  tempOptions[["AMU"]] = currentList$AMU
  tempOptions[["baselineStart"]] = 0
  tempOptions[["baselineEnd"]] = .5
  tempOptions[["baselineType"]] = "mean"
  tempOptions[["timeStart"]] = 0
  tempOptions[["timeEnd"]] = 6
  tempOptions[["inert"]] = "AMU40"
  tempOptions[["isProduct"]] = F
  tempOptions[["calibrationCoef"]] = 1
  tempOptions[["diffusionCoef"]] = 0
  tempOptions[["savGolSmoothing"]] = 0
  tempOptions[["yProcSmoothing"]] = 0
  tempOptions[["waveSmoothing"]] = 0
  currentList$options = tempOptions
  TAPexperiment[[currentList$options$Name]] = currentList
  
  # Ar
  currentList = list()
  # Assign AMU
  currentList$AMU = 40
  currentList$time = solverDF$timeVector
  currentList$temperature = 30
  currentList$pulses = data.frame("1" = solverDF$inert, "2" = solverDF$inert)
  
  tempOptions = list()
  tempOptions[["Name"]] = paste0("AMU", currentList$AMU)
  tempOptions[["AMU"]] = currentList$AMU
  tempOptions[["baselineStart"]] = 0
  tempOptions[["baselineEnd"]] = .5
  tempOptions[["baselineType"]] = "mean"
  tempOptions[["timeStart"]] = 0
  tempOptions[["timeEnd"]] = 6
  tempOptions[["inert"]] = "AMU40"
  tempOptions[["isProduct"]] = T
  tempOptions[["calibrationCoef"]] = 1
  tempOptions[["diffusionCoef"]] = 0
  tempOptions[["savGolSmoothing"]] = 0
  tempOptions[["yProcSmoothing"]] = 0
  tempOptions[["waveSmoothing"]] = 0
  currentList$options = tempOptions
  TAPexperiment[[currentList$options$Name]] = currentList
  
  TAPexperiment = TAPcodeV2::inertGraham(TAPexperiment, "AMU40", 28)
  TAPexperiment = TAPcodeV2::inertGraham(TAPexperiment, "AMU44", 28)
  TAPexperiment$AMU28$options$inert = "AMU40to28"
  TAPexperiment = gasCon(TAPexperiment, c("AMU28", "AMU44"), method = "g" )
  TAPexperiment = TAPcodeV2::reactionRate(TAPexperiment, c("AMU28", "AMU44","AMU44to28"), method = "g" )
  
  TAPexperiment = TAPcodeV2::moments(TAPexperiment, c("AMU28gasCon", "AMU44gasCon", "AMU28rate", "AMU44rate", "AMU44to28rate"))
  return(TAPexperiment)
}



solverToObj2 = function(dataPath){

  co2 = read.csv(paste0(dataPath,"CO2.csv"), header = F)
  co = read.csv(paste0(dataPath,"CO.csv"), header = F)
  inert = read.csv(paste0(dataPath,"Inert.csv"), header = F)
  timeVector = inert[,1]
  timeRange = 1:which.min(abs(timeVector - 3))
  
  co2 = co2[timeRange, ]
  co = co[timeRange, ]
  inert = inert[timeRange, ]
  timeVector = timeVector[timeRange]
  
  
  TAPexperiment = list()
  
  reactorDF = data.frame("reactorParameters" = c("inertZone1Length", "inertZone2Length",
                                                 "catalystBedLength", "crossSectionalArea",
                                                 "catalystWeight", "bedPorosity",
                                                 "molPerM0Inert", "pumpProbeSpacing"),
                         "value" = c(0.02894, 0.02894, 1.6e-3, 2e-03, 1, 0.4, 1.63E-09, 0),
                         "unit" = c("m", "m", "m", "m^2", "g", "none", "mol", "s"))
  reactorParams = as.list(reactorDF$value)
  names(reactorParams) = reactorDF$reactorParameters
  TAPexperiment[["reactor"]] = list("reactorParams" = reactorParams,"reactorDF" = reactorDF)
  
  # CO2
  currentList = list()
  # Assign AMU
  currentList$AMU = 44
  currentList$time = timeVector
  currentList$temperature = 30
  currentList$pulses = co2[,-1]
  
  tempOptions = list()
  tempOptions[["Name"]] = paste0("AMU", currentList$AMU)
  tempOptions[["AMU"]] = currentList$AMU
  tempOptions[["baselineStart"]] = 0
  tempOptions[["baselineEnd"]] = .5
  tempOptions[["baselineType"]] = "mean"
  tempOptions[["timeStart"]] = 0
  tempOptions[["timeEnd"]] = 6
  tempOptions[["inert"]] = "AMU40"
  tempOptions[["isProduct"]] = T
  tempOptions[["calibrationCoef"]] = 1
  tempOptions[["diffusionCoef"]] = 0
  tempOptions[["savGolSmoothing"]] = 0
  tempOptions[["yProcSmoothing"]] = 0
  tempOptions[["waveSmoothing"]] = 0
  currentList$options = tempOptions
  TAPexperiment[[currentList$options$Name]] = currentList
  
  # CO
  currentList = list()
  # Assign AMU
  currentList$AMU = 28
  currentList$time = timeVector
  currentList$temperature = 30
  currentList$pulses = co[, -1]
  
  tempOptions = list()
  tempOptions[["Name"]] = paste0("AMU", currentList$AMU)
  tempOptions[["AMU"]] = currentList$AMU
  tempOptions[["baselineStart"]] = 0
  tempOptions[["baselineEnd"]] = .5
  tempOptions[["baselineType"]] = "mean"
  tempOptions[["timeStart"]] = 0
  tempOptions[["timeEnd"]] = 6
  tempOptions[["inert"]] = "AMU40"
  tempOptions[["isProduct"]] = F
  tempOptions[["calibrationCoef"]] = 1
  tempOptions[["diffusionCoef"]] = 0
  tempOptions[["savGolSmoothing"]] = 0
  tempOptions[["yProcSmoothing"]] = 0
  tempOptions[["waveSmoothing"]] = 0
  currentList$options = tempOptions
  TAPexperiment[[currentList$options$Name]] = currentList
  
  # Ar
  currentList = list()
  # Assign AMU
  currentList$AMU = 40
  currentList$time = timeVector
  currentList$temperature = 30
  currentList$pulses = inert[, -1]
  
  tempOptions = list()
  tempOptions[["Name"]] = paste0("AMU", currentList$AMU)
  tempOptions[["AMU"]] = currentList$AMU
  tempOptions[["baselineStart"]] = 0
  tempOptions[["baselineEnd"]] = .5
  tempOptions[["baselineType"]] = "mean"
  tempOptions[["timeStart"]] = 0
  tempOptions[["timeEnd"]] = 6
  tempOptions[["inert"]] = "AMU40"
  tempOptions[["isProduct"]] = T
  tempOptions[["calibrationCoef"]] = 1
  tempOptions[["diffusionCoef"]] = 0
  tempOptions[["savGolSmoothing"]] = 0
  tempOptions[["yProcSmoothing"]] = 0
  tempOptions[["waveSmoothing"]] = 0
  currentList$options = tempOptions
  TAPexperiment[[currentList$options$Name]] = currentList
  
  TAPexperiment = TAPcodeV2::inertGraham(TAPexperiment, "AMU40", 28)
  TAPexperiment = TAPcodeV2::inertGraham(TAPexperiment, "AMU44", 28)
  TAPexperiment$AMU28$options$inert = "AMU40to28"
  TAPexperiment = gasCon(TAPexperiment, c("AMU28", "AMU44"), method = "g" )
  TAPexperiment = TAPcodeV2::reactionRate(TAPexperiment, c("AMU28", "AMU44","AMU44to28"), method = "g" )
  
  TAPexperiment = TAPcodeV2::moments(TAPexperiment, c("AMU28gasCon", "AMU44gasCon", "AMU28rate", "AMU44rate", "AMU44to28rate"))
  return(TAPexperiment)
}





readElucFolder = function(dataPath){
  co2 = read.csv(paste0(dataPath,"CO2.csv"), header = F)
  co = read.csv(paste0(dataPath,"CO.csv"), header = F)
  inert = read.csv(paste0(dataPath,"Inert.csv"), header = F)
  timeVector = inert$V1
  co2 = co2$V2
  co = co$V2
  inert = inert$V2
  
  # cutting to .5 seconds to make it more clear
  timeRange = 1:which.min(abs(timeVector - .3))
  co2 = co2[timeRange]
  co = co[timeRange]
  inert = inert[timeRange]
  timeVector = timeVector[timeRange]
  
  # mixup with the labeled gas species
  tempCO = co
  tempCO2 = co2
  tempInert = inert
  
  inert = tempCO2
  co2 = tempCO
  co = tempInert
  
  solverDF = data.frame("co" = co, "co2" = co2, "inert" = inert, "timeVector" = timeVector)
  return(solverDF)
}

plotGas = function(TAPexperiment, pulseNum = 1, Type = ""){
  # if(length(pulseNum) > 1){
  # 
  # }else{
  # 
  # }
  
  plotDF = data.frame("flux" = c(TAPexperiment$AMU44$pulses[,pulseNum], TAPexperiment$AMU28$pulses[,pulseNum], TAPexperiment$AMU40$pulses[,pulseNum]),
                      "Time" = rep(TAPexperiment$AMU40$time, 3),
                      "Type" = rep(c("CO2", "CO", "Inert"), each = length(TAPexperiment$AMU40$time)))
  m = ggplot(plotDF, aes(x = Time, y = flux, color = Type)) + geom_line(size = 1) +
    ggplot2::theme(text = ggplot2::element_text(size = 24),
                   legend.position = c(.8,.8),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   panel.grid.minor = element_line(colour = "black", linetype = "dashed", size = .05),
                   panel.grid.major = element_line(colour = "black", linetype = "dashed", size = .1),
                   panel.border = element_rect(fill = "NA"),
                   legend.key.width = unit(1.5,"cm"),
                   legend.key = element_rect(fill = "NA")
    ) +
    ggplot2::ggtitle(paste0(Type, "Outlet")) +
    ggplot2::xlab("Time (s)") +
    ggplot2::ylab("Flux") +
    labs(color='Type')
  m
  
}

plotRates = function(TAPexperiment,pulseNum, Type = ""){
  plotDF = data.frame("flux" = c(TAPexperiment$AMU28rate$pulses[,pulseNum], 
                                 TAPexperiment$AMU44rate$pulses[,pulseNum]),
                      "Time" = rep(TAPexperiment$AMU44rate$time, 2),
                      "Type" = rep(c("CO", "CO2"), each = length(TAPexperiment$AMU44rate$time)))
  m = ggplot(plotDF, aes(x = Time, y = flux, color = Type)) + geom_line(size = 1) +
    ggplot2::theme(text = ggplot2::element_text(size = 24),
                   legend.position = c(.8,.8),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   panel.grid.minor = element_line(colour = "black", linetype = "dashed", size = .05),
                   panel.grid.major = element_line(colour = "black", linetype = "dashed", size = .1),
                   panel.border = element_rect(fill = "NA"),
                   legend.key.width = unit(1.5,"cm"),
                   legend.key = element_rect(fill = "NA")
    ) +
    ggplot2::ggtitle(paste0(Type, " Rates")) +
    ggplot2::xlab("Time (s)") +
    ggplot2::ylab("Rate") +
    labs(color='Type')
  m
}


createPredMat = function(TAPexperiment, normalization = T, logTransform = F){
  if(normalization){
    df = data.frame("r28" = TAPexperiment$AMU28rate$pulses[,1] / TAPexperiment$AMU28rate$moments$M0[1],
                    "g28" = TAPexperiment$AMU28gasCon$pulses[,1] / TAPexperiment$AMU28gasCon$moments$M0[1],
                    "r44" = TAPexperiment$AMU44rate$pulses[,1] / TAPexperiment$AMU44rate$moments$M0[1])
    df$s28 = cumsum(df$r28)
    df$timeVector = TAPexperiment$AMU44gasCon$time
    df$g28r28 = df$g28 * df$r28
    df$g28s28 = df$g28 * df$s28
    df$r28s28 = df$r28 * df$s28
    df$g28t = df$g28 * df$timeVector
    df$r28t = df$r28 * df$timeVector
    df$s28t = df$s28 * df$timeVector
    df$g28r28s28 = df$g28 * df$s28 * df$r28
    df$g28r28t = df$g28 * df$r28 * df$timeVector
    df$g28s28t = df$g28 * df$s28 * df$timeVector
    df$r28s28t = df$r28 * df$s28 * df$timeVector
    df$g28r28s28t = df$g28 * df$r28 * df$s28 * df$timeVector
  }else{
    df = data.frame("r28" = TAPexperiment$AMU28rate$pulses[,1],
                    "g28" = TAPexperiment$AMU28gasCon$pulses[,1],
                    "r44" = TAPexperiment$AMU44rate$pulses[,1])
    df$s28 = cumsum(df$r28)
    df$timeVector = TAPexperiment$AMU44gasCon$time
    df$g28r28 = df$g28 * df$r28
    df$g28s28 = df$g28 * df$s28
    df$r28s28 = df$r28 * df$s28
    df$g28t = df$g28 * df$timeVector
    df$r28t = df$r28 * df$timeVector
    df$s28t = df$s28 * df$timeVector
    df$g28r28s28 = df$g28 * df$s28 * df$r28
    df$g28r28t = df$g28 * df$r28 * df$timeVector
    df$g28s28t = df$g28 * df$s28 * df$timeVector
    df$r28s28t = df$r28 * df$s28 * df$timeVector
    df$g28r28s28t = df$g28 * df$r28 * df$s28 * df$timeVector
  }
  
  if(logTransform){
    df = data.frame("r28" = log(TAPexperiment$AMU28rate$pulses[,1] - min(TAPexperiment$AMU28rate$pulses[,1]) + 1e-10),
                    "g28" = log(TAPexperiment$AMU28gasCon$pulses[,1] - min(TAPexperiment$AMU28gasCon$pulses[,1]) + 1e-10),
                    "r44" = log(TAPexperiment$AMU44rate$pulses[,1] - min(TAPexperiment$AMU44rate$pulses[,1]) + 1e-10))
    df$s28 = cumsum(df$r28)
    df$timeVector = TAPexperiment$AMU44gasCon$time
    df$g28r28 = df$g28 * df$r28
    df$g28s28 = df$g28 * df$s28
    df$r28s28 = df$r28 * df$s28
    df$g28t = df$g28 * df$timeVector
    df$r28t = df$r28 * df$timeVector
    df$s28t = df$s28 * df$timeVector
    df$g28r28s28 = df$g28 * df$s28 * df$r28
    df$g28r28t = df$g28 * df$r28 * df$timeVector
    df$g28s28t = df$g28 * df$s28 * df$timeVector
    df$r28s28t = df$r28 * df$s28 * df$timeVector
    df$g28r28s28t = df$g28 * df$r28 * df$s28 * df$timeVector
    
  }
  return(df)
}


glmnetTAP = function(TAPexperiment, alpha = .5, logTransform = F, normalization = T){
  
  predMat = createPredMat(TAPexperiment, logTransform = logTransform, normalization = normalization)
  response = predMat$r44
  timeVector = predMat$timeVector
  
  predMat = scale(predMat[,-which(names(predMat) %in% c("r44") )])
  #predMat = (cbind(predMat$rate28, predMat$gasS28, predMat$gasS28 * predMat$rate28))
  
  flag = T
  iteration = 0
  fit = cv.glmnet(predMat, response, alpha = alpha, intercept = F)
  while(flag){
    if(fit$lambda.min != min(fit$lambda)){
      break
    }
    newLamb = fit$lambda / 2
    fit = cv.glmnet(predMat, response, alpha = alpha, lambda = newLamb, intercept = F)
    iteration = iteration + 1
    #print(iteration)
  }
  beta = fit$glmnet.fit$beta[,which(fit$lambda == fit$lambda.min)]
  yHat = predict(fit, newx = predMat, s = "lambda.min")
  result = list("beta" = beta, "yHat" = yHat, "response" = response, "timeVector"= timeVector,
                "residuals" = response - yHat)
}


plotFit = function(fitResults, Type = ""){
  plotDF = data.frame("flux" = c(fitResults$response, 
                                 fitResults$yHat),
                      "Time" = rep(fitResults$timeVector, 2),
                      "Type" = rep(c("response", "predicted"), each = length(fitResults$timeVector)))
  m = ggplot(plotDF, aes(x = Time, y = flux, color = Type)) + geom_line(size = 1, alpha = .4) +
    ggplot2::theme(text = ggplot2::element_text(size = 24),
                   legend.position = c(.8,.8),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   panel.grid.minor = element_line(colour = "black", linetype = "dashed", size = .05),
                   panel.grid.major = element_line(colour = "black", linetype = "dashed", size = .1),
                   panel.border = element_rect(fill = "NA"),
                   legend.key.width = unit(1.5,"cm"),
                   legend.key = element_rect(fill = "NA")
    ) +
    ggplot2::ggtitle(paste0(Type, " Predict CO2")) +
    ggplot2::xlab("Time (s)") +
    ggplot2::ylab("Rate") +
    labs(color='Type')
  m
}



hierNetTAP = function(TAPexperiment, diagVal = T, strong = T, logTransform = F, normalization = F){
  predMat = createPredMat(TAPexperiment, logTransform = logTransform, normalization = normalization)
  response = predMat$r44
  timeVector = predMat$timeVector
  predMat = scale(predMat[,which(names(predMat) %in% c("g28", "r28", "s28", "timeVector") ) ])
  
  fitPath = hierNet.path(predMat,response,strong = strong, diagonal = diagVal)
  fitcv = hierNet.cv(fitPath,predMat,response)
  fit2=hierNet(predMat,response,lam=fitcv$lamhat, strong = strong, diagonal = diagVal)
  
  yHat=predict(fit2,predMat)
  intEffects = fit2$th
  colnames(intEffects) = colnames(predMat)
  rownames(intEffects) = colnames(predMat)
  result = list("intEffects" = intEffects, "posCoefs" = fit2$bp, "negCoefs" = fit2$bn, "yHat" = yHat,
                "response" = response, "timeVector" = timeVector, "residuals" = response - yHat)
  return(result)
}


rmse = function(fitResults){
  sqrt(mean((fitResults$residuals)^2))
}



covTAP = function(TAPexperiment, corVal = F, logTransform = F, normalization = T){
  predMat = createPredMat(TAPexperiment, logTransform = logTransform, normalization = normalization)
  #response = predMat$rate28
  timeVector = predMat$timeVector
  #predMat = predMat[, which(names(predMat) %in% c("rate28", "rate44", "gas28", "surface44", "surface28") )]
  predMat = predMat[, which(names(predMat) %in% c("r28", "g28", "s28", "r44", "timeVector") )]
  
  if(corVal){
    tempCor = cor(predMat)  
  }else{
    tempCor = cov(predMat)
  }
  
  m = corrr::network_plot(cor(predMat), min_cor = 0)
  result = list("matrix" = tempCor, "plot" = m)
  return(result)
}



