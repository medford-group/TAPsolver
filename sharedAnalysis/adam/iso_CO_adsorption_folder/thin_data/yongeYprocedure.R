setwd("C:/Users/ayong/Documents/tap_solver/TAPsolver/sharedAnalysis/adam/iso_CO_adsorption_folder/flux_data/")
knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)
library(TAPcodeV2)
source("./supportScripts.R")

dataPath = "C:/Users/ayong/Documents/tap_solver/TAPsolver/sharedAnalysis/adam/iso_CO_adsorption_folder/flux_data/"
solverDFE = readElucFolder(dataPath)
TAPexperimentE = solverToObj(solverDFE)
#print(TAPexperimentE$AMU28rate$time)
#predMatE = createPredMat(TAPexperimentE)
#print(TAPexperimentE$AMU44rate$pulses$X2)
#print(TAPexperimentE$AMU44gasCon$pulses$X2)
print(TAPexperimentE$AMU28rate$pulses$X2)
print(TAPexperimentE$AMU28gasCon$pulses$X2)
#print(TAPexperimentE$AMU44rate$reactionRateScalar)
print(TAPexperimentE$AMU28rate$reactionRateScalar)

#print(typeof(TAPexperimentE$AMU44gasCon$pulses$X2))
#print(TAPexperimentE$AMU44gasCon$pulses$X2)
print(TAPexperimentE$AMU28gasCon$gasConScalar)

#print(TAPexperimentE$AMU40to28$)
print(TAPexperimentE$AMU40to28$pulses$X1)
print(TAPexperimentE$AMU40to28$pulses$X2)
write.table(TAPexperimentE$AMU28rate$pulses$X2,row.names=FALSE,col.names=FALSE,'g_CO.csv')
write.table(TAPexperimentE$AMU28gasCon$pulses$X2,row.names=FALSE,col.names=FALSE,'g_con_CO.csv')

plotGas(TAPexperimentE,1, "ER ")
plotRates(TAPexperimentE,1 , "ER ")

covInfo = covTAP(TAPexperimentE, corVal = F, normalization = T)
covInfo$plot

fit = lm(r44 ~ g28 * r28 * s28 * timeVector - 1, data = predMatE)
fitResults = data.frame("yHat" = fit$fitted.values, "response" = predMatE$r44, "timeVector" = predMatE$timeVector)
plotFit(fitResults, "ER OLS")
summary(fit)
