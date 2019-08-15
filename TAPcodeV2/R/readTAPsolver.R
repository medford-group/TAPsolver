#' Read TAPsolver Output
#'
#' Reads the output of the TAPsolver and turns it into the required TAPstructure.
#'
#' This function is based on R package RJSONIO.
#'
#' @param dataPath This is a string path to the output folder of the TAPsolver Experiment
#' @return A collection of TAP objects that describe the experiment performed in the reactor.
#' @examples
#'
#' print("none")
#'
#' @export readTAPsolver

readTAPsolver = function(dataPath){
  inputFile = RJSONIO::readJSONStream(paste0(dataPath, "inputFileRoss.txt"))
  TAPsolverStruct = list("parameters" = NULL)
  TAPsolverStruct$parameters = list(
    "reactor" = NULL,
    "experiment" = NULL,
    "simulation" = NULL,
    "analysis" = NULL
  )

  tempCatalystLen = inputFile$reactor_kinetics_input$`Reactor Length` * inputFile$reactor_kinetics_input$`Catalyst Fraction`
  tempInertZone = (inputFile$reactor_kinetics_input$`Reactor Length` - tempCatalystLen) / 2

  TAPsolverStruct$parameters$reactor = list(
    "inertZone1Length" = tempInertZone,
    "inertZone2Length" = tempInertZone,
    "catalystBedLength" = tempCatalystLen,
    "reactorLength" = inputFile$reactor_kinetics_input$`Reactor Length`,
    "crossSectionalArea" = pi * inputFile$reactor_kinetics_input$`Reactor Radius`^2,
    "catalystWeight" = 1,
    "bedPorosity" = inputFile$reactor_kinetics_input$`Void Fraction Catalyst`,
    "molPerM0Inert" = inputFile$reactor_kinetics_input$`Reference Pulse Size`
  )

  pulseStartTimes = as.numeric(strsplit(inputFile$reactor_kinetics_input$`Pulse Time`,",")[[1]])
  TAPsolverStruct$parameters$experiment = list(
    "numTime" = inputFile$reactor_kinetics_input$`Time Steps`,
    "numPulses" = inputFile$reactor_kinetics_input$`Number of Pulses`,
    "time" = seq(as.numeric(strsplit(inputFile$reactor_kinetics_input$`Pulse Time`,",")[[1]][1]), inputFile$reactor_kinetics_input$`Pulse Duration`, inputFile$reactor_kinetics_input$`Time Steps`),
    "timeStart" = pulseStartTimes[1],
    "timeEnd" = inputFile$reactor_kinetics_input$`Pulse Duration`,
    "timeIndex" = c(1:inputFile$reactor_kinetics_input$`Time Steps`),
    "pumpProbeSpacing" = pulseStartTimes[1] - pulseStartTimes[2],
    "temperature" = rep(inputFile$reactor_kinetics_input$`Reactor Temperature`, inputFile$reactor_kinetics_input$`Number of Pulses`),
    "diffusionCoef" = inputFile$reactor_kinetics_input$`Reference Diffusion Inert`,
    "inert" = "Inert",
    "reactant" = NULL,
    "products" = NULL,
    "inert2" = NULL,
    "reactant2" = NULL,
    "products2" = NULL,
    "numCores" = 4
  )

  tempKnames = names(inputFile$kinetic_parameters)
  maxK = 1
  kBackward = rep(0,20)
  kForward = rep(0,20)
  for(i in 1:length(tempKnames)){
    lastK = as.numeric(substr(tempKnames[i],3,nchar(tempKnames))) + 1
    forwardOrBack = substr(tempKnames[i],2,2)
    if(forwardOrBack == "f"){
      kForward[lastK] = inputFile$kinetic_parameters[i]
    }else{
      kBackward[lastK] = inputFile$kinetic_parameters[i]
    }
    if(lastK > maxK){
      maxK = lastK
    }
  }
  kBackward = kBackward[1:maxK]
  kForward = kForward[1:maxK]

  TAPsolverStruct$parameters$simulation = list(
    "meshSize" = inputFile$reactor_kinetics_input$`Mesh Size`,
    "outputFolderName" = inputFile$reactor_kinetics_input$`Output Folder Name`,
    "experimentalFolderName" = inputFile$reactor_kinetics_input$`Experimental Data Folder`,
    "noise" = inputFile$reactor_kinetics_input$Noise,
    "reactorType" = inputFile$reactor_kinetics_input$`Reactor Type`,
    "theta" = inputFile$reactor_kinetics_input$Theta,
    "solverMethod" = inputFile$reactor_kinetics_input$`Solver Method`,
    "storeOutletFlux" = inputFile$reactor_kinetics_input$`Store Outlet Flux`,
    "storeGraph" = inputFile$reactor_kinetics_input$`Store Graph`,
    "displayExperimentalData" = inputFile$reactor_kinetics_input$`Display Experimental Data`,
    "displayGraph" = inputFile$reactor_kinetics_input$`Display Graph`,
    "sensitivityAnalysis" = inputFile$reactor_kinetics_input$`Sensitivity Analysis`,
    "fitParameters" = inputFile$reactor_kinetics_input$`Fit Parameters`,
    "optimizationMethod" = inputFile$reactor_kinetics_input$`Optimization Method`,
    "objectivePoints" = inputFile$reactor_kinetics_input$`Objective Points`,
    "rrmAnalysis" = inputFile$reactor_kinetics_input$`RRM Analysis`,
    "mkmAnalysis" = inputFile$reactor_kinetics_input$`MKM Analysis`,
    "petalPlots" = inputFile$reactor_kinetics_input$`Petal Plots`,
    "reactionsTest" = inputFile$reactor_kinetics_input$reactions_test, # this form is very explicit, see TAPsolver documentation
    # gas species requirements.  For example, gas A and an Inert value
    "massList" = inputFile$reactor_kinetics_input$`Mass List`, # first value is the mass of A, last value is the inert. Total length of the vector depends on the number of gas species.
    "pulseTime" = inputFile$reactor_kinetics_input$`Pulse Time`, # pulseTiming of each gas
    "pulseRatio" = inputFile$reactor_kinetics_input$`Pulse Ratio`, # balance of the gas injected into the reactor
    # surface Requirements
    "initialSurfaceComposition" = inputFile$reactor_kinetics_input$`Initial Surface Composition`, # depends on the number of sites. For example A* and *.  Similar to the inert example in the gas species, the open sites are always placed last (* has 1000 sites)
    # Elementary step requirements
    "kForward" = kForward, # length depends on the length of reactionsTest
    "kBackward" = kBackward # length depends on the length of reactionsTest
  )

  # construct flux data
  subDataPath = paste0(dataPath, "flux_data/")
  tempFileNames = list.files(subDataPath)
  for(i in 1:length(tempFileNames)){
    # currently there is no connection to the actual mass or initial surface constribution. May want to spend some time to build a class in python for TAPsolver.
    tempData = read.csv(paste0(subDataPath, tempFileNames[i]), header = F)
    TAPobj = TAPcodeV2::constructFlux()
    for(j in 2:dim(tempData)[2]){
      TAPobj$matrix[[j-1]] = tempData[,j]
    }
    tempName = substr(tempFileNames[i], 1, (nchar(tempFileNames[i]) - 4))
    TAPobj$options$name = tempName
    TAPsolverStruct[[TAPobj$options$name]] = TAPobj
  }
  # construct TAPobj for each different csv in folder
  subDataPath = paste0(dataPath, "thin_data/")
  tempFileNames = list.files(subDataPath)
  for(i in 1:length(tempFileNames)){
    # currently there is no connection to the actual mass or initial surface constribution. May want to spend some time to build a class in python for TAPsolver.
    tempData = read.csv(paste0(subDataPath, tempFileNames[i]), header = F)
    TAPobj = TAPcodeV2::constructFlux()
    for(j in 2:dim(tempData)[2]){
      TAPobj$matrix[[j-1]] = tempData[,j]
    }
    tempName = substr(tempFileNames[i], 1, (nchar(tempFileNames[i]) - 4))
    TAPobj$options$name = tempName
    TAPsolverStruct[[TAPobj$options$name]] = TAPobj
  }

  return(TAPsolverStruct)
}

