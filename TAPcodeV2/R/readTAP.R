#' Read Excel TDMS
#'
#' Reads in an Excel .xlsx of the original TDMS file.  This requires that data has been previously converted to .xlsx by the tdms addon in Excel.  This is the function that will be called when an user uploads data.
#'
#' This function is based on R package rio.
#'
#' @param dataPath This is a string path to an .tdms file.
#' @return A collection of TAP objects that describe the experiment performed in the reactor.
#' @examples
#' TAPexperiment = readTAP("data/pumpProbe.tdms")
#' names(TAPexperiment)
#' @export readTAP

readTAP = function(dataPath){
  # constructor for TAP parameters
  TAPstruct = list("parameters" = NULL)
  TAPstruct$parameters = list(
    "reactor" = NULL,
    "experiment" = NULL,
    "simulation" = NULL,
    "analysis" = ""
  )
  TAPstruct$parameters$reactor = list(
    "inertZone1Length" = 0.02,
    "inertZone2Length" = 0.02,
    "catalystBedLength" = 0.00075,
    "reactorLength" = .04075,
    "crossSectionalArea" = 1.14009E-05,
    "catalystWeight" = 1,
    "bedPorosity" = 0.4,
    "molPerM0Inert" = 1.63E-09
  )
  TAPstruct$parameters$experiment = list(
    "numTime" = 100,
    "numPulses" = 1,
    "time" = seq(0,2,length.out = 100),
    "timeStart" = 0,
    "timeEnd" = 2,
    "timeIndex" = (1:100),
    "pumpProbeSpacing" = 0,
    "temperature" = 23,
    "diffusionCoef" = 0,
    "inert" = "AMU40",
    "reactant" = "AMU40",
    "products" = "AMU40",
    "inert2" = "AMU40",
    "reactant2" = "AMU40",
    "products2" = "AMU40",
    "numCores" = 4
  )
  TAPstruct$parameters$simulation = list(
    "meshSize" = 360,
    "outputFolderName" = "results",
    "experimentalFolderName" = "none",
    "noise" = "FALSE",
    "reactorType" = "tap",
    "theta" = 1,
    "solverMethod" = "None",
    "storeOutletFlux" = "TRUE",
    "storeGraph" = "FALSE",
    "displayExperimentalData" = "FALSE",
    "displayGraph" = "FALSE",
    "sensitivityAnalysis" = "FALSE",
    "fitParameters" = "FALSE",
    "optimizationMethod" = "CG",
    "objectivePoints" = 1,
    "rrmAnalysis" = "FALSE",
    "mkmAnalysis" = "TRUE",
    "petalPlots" = "FALSE",
    "reactionsTest" = "A + * -> A*", # this form is very explicit, see TAPsolver documentation
    # gas species requirements.  For example, gas A and an Inert value
    "massList" = c(40,40), # first value is the mass of A, last value is the inert. Total length of the vector depends on the number of gas species.
    "pulseTime" = c(0,0), # pulseTiming of each gas
    "pulseRatio" = c(1,1), # balance of the gas injected into the reactor
    # surface Requirements
    "initialSurfaceComposition" = c(0,1000), # depends on the number of sites. For example A* and *.  Similar to the inert example in the gas species, the open sites are always placed last (* has 1000 sites)
    # Elementary step requirements
    "kForward" = c(200), # length depends on the length of reactionsTest
    "kBackward" = c(0) # length depends on the length of reactionsTest
  )
  if(is.null(dataPath)){
    return(TAPstruct)
  }

  if(substr(dataPath,(nchar(dataPath)-3), nchar(dataPath)) == "tdms"){
    sheetPos = 3
    currentName = NULL

    library(reticulate)
    #use_python("/opt/anaconda/bin/python3.7")
    #use_virtualenv("r-reticulate")
    nptdms = reticulate::import("nptdms")
    tempFile = nptdms$TdmsFile(dataPath)
    groupNames = nptdms$TdmsFile$groups(tempFile)

    temperatureForAll = list()
    tempTime = NULL
    tempNumPulses = list()
    nonIntegerAMUvalues = tempFile$object(groupNames[2])
    nonIntegerAMUvalues = nonIntegerAMUvalues$as_dataframe()

    # need information from flux
    tempIndex = 1
    for(i in sheetPos:length(groupNames)){
      attempt = tempFile$object(groupNames[i])
      attempt = attempt$as_dataframe()

      # Assign Time, only need to do this once
      if(i == sheetPos){
        tempTime = attempt$Time
        tempTime[length(tempTime)] = tempTime[2] + tempTime[length(tempTime) - 1]
        tempTime = tempTime[-1]
      }

      # Collect temperature information for Gas Species
      temperatureForAll[[tempIndex]] = as.numeric(unlist(attempt[1,],F,F)[-(1:3)])

      # Flux and options information
      TAPobj = TAPcodeV2::constructFlux()
      TAPobj$matrix = as.list(attempt[-1,-(1:3)])
      TAPobj$temperature = as.numeric(unlist(attempt[1,],F,F)[-(1:3)])
      #TAPobj$options$name = paste0("AMU",attempt$Value[1]) # Need to change this due to nonInteger amu values
      TAPobj$options$name = paste0("AMU",nonIntegerAMUvalues$Value[i - sheetPos + 1])
      TAPobj$options$amu = as.numeric(attempt$Value[1])
      TAPobj$options$gain = as.numeric(attempt$Value[2])
      TAPobj$options$baselineStart = min(tempTime)
      TAPobj$options$baselineEnd = max(tempTime)
      TAPobj$options$inert = paste0("AMU",attempt$Value[1])

      # store the total number of pulses in case the experiment was cut short of a complete set
      tempNumPulses[[i]] = length(TAPobj$matrix)

      TAPstruct[[TAPobj$options$name]] = TAPobj

      tempIndex = tempIndex + 1
    }
  }else if(substr(dataPath,(nchar(dataPath)-3), nchar(dataPath)) == "xlsx"){
    result = list()
    sheetPos = 5
    currentName = NULL

    # build reactor information
    reactorParams = rio::import(dataPath, which = 1)
    for(i in 1:length(reactorParams$Value)){
      TAPstruct$parameters$reactor[[i]] = reactorParams$Value[i]
    }
    # build experiment information
    experimentParams = rio::import(dataPath, which = 2)
    for(i in 1:length(experimentParams$Value)){
      if(grepl(",",experimentParams$Value[i])){
        # account for string vectors
        tempValue = unlist(strsplit(experimentParams$Value[i], split = ","))
        if(sum(is.na(as.numeric(tempValue))) != 0){
          TAPstruct$parameters$experiment[[i]] = tempValue
        }else{
          TAPstruct$parameters$experiment[[i]] = as.numeric(tempValue)
        }
      }else{
        TAPstruct$parameters$experiment[[i]] = experimentParams$Value[i]
      }
    }
    # build simulation information
    simulationParams = rio::import(dataPath, which = 3)
    for(i in 1:length(simulationParams$Value)){
      if(grepl(",",simulationParams$Value[i])){
        # account for string vectors
        tempValue = unlist(strsplit(simulationParams$Value[i], split = ","))
        if(sum(is.na(as.numeric(tempValue))) != 0){
          TAPstruct$parameters$experiment[[i]] = tempValue
        }else{
          TAPstruct$parameters$experiment[[i]] = as.numeric(tempValue)
        }
      }else{
        TAPstruct$parameters$experiment[[i]] = simulationParams$Value[i]
      }
    }
    # build analysis information
    TAPstruct$parameters$analysis = rio::import(dataPath, which = 4)

    temperatureForAll = list()
    tempTime = NULL
    tempNumPulses = list()

    while(sheetPos != 0){
      attempt = try(rio::import(dataPath, which = sheetPos), silent = T)
      if(is.data.frame(attempt)){
        TAPobj = TAPcodeV2::constructFlux()

        TAPobj$options$name = attempt$Value[1]
        TAPobj$options$amu = as.numeric(attempt$Value[2])
        TAPobj$options$gain = as.numeric(attempt$Value[3])
        TAPobj$options$baselineStart = as.numeric(attempt$Value[4])
        TAPobj$options$baselineEnd = as.numeric(attempt$Value[5])
        TAPobj$options$baselineType = attempt$Value[6]
        TAPobj$options$inert = attempt$Value[7]
        if(grepl("T", toupper(attempt$Value[8]))){
          TAPobj$options$isProduct = TRUE
        }else{
          TAPobj$options$isProduct = FALSE
        }
        TAPobj$options$calibrationCoef = as.numeric(attempt$Value[9])
        TAPobj$options$savGolSmoothing = as.numeric(attempt$Value[10])
        TAPobj$options$yProcSmoothing = as.numeric(attempt$Value[11])
        TAPobj$options$waveSmoothing = as.numeric(attempt$Value[12])


        # Assign Time values
        if(i == sheetPos){
          tempTime = attempt$Time
          timeLen = length(tempTime)
          tempTime[timeLen] = tempTime[2] + tempTime[timeLen - 1]
          tempTime = tempTime[-1]
        }
        temperatureForAll[[i]] = as.numeric(unlist(attempt[1,],F,F)[-(1:3)])
        TAPobj$matrix = as.list(attempt[-1,-(1:3)])
        tempNumPulses[[i]] = length(TAPobj$matrix)

        TAPstruct[[TAPobj$options$name]] = TAPobj
        sheetPos = sheetPos + 1
      }else{
        sheetPos = 0
      }
    }
  }

  for(i in 2:length(TAPstruct)){
    tempNumPulses = unlist(tempNumPulses)
    if(tempNumPulses[i-1] != min(tempNumPulses)){
      TAPstruct[[i]]$matrix = TAPstruct[[i]]$matrix[1:min(tempNumPulses)]
      temperatureForAll[[i-1]] = temperatureForAll[[i-1]][1:min(tempNumPulses)]
    }
  }
  temperatureMean = matrix(unlist(temperatureForAll), nrow = (length(TAPstruct) - 1), byrow = FALSE)

  # Update Experiment
  TAPstruct$parameters$experiment$numTime = length(TAPstruct[[2]]$matrix[[1]])
  TAPstruct$parameters$experiment$numPulses = length(TAPstruct[[2]]$matrix)
  TAPstruct$parameters$experiment$time = tempTime
  TAPstruct$parameters$experiment$timeStart = min(tempTime)
  TAPstruct$parameters$experiment$timeEnd = max(tempTime)
  TAPstruct$parameters$experiment$timeIndex = 1:length(tempTime)
  TAPstruct$parameters$experiment$temperature = apply(temperatureMean,2,mean)
  TAPstruct$parameters$experiment$inert = TAPstruct[[2]]$options$name
  TAPstruct$parameters$experiment$reactant = TAPstruct[[2]]$options$name
  TAPstruct$parameters$experiment$products = TAPstruct[[2]]$options$name
  TAPstruct$parameters$experiment$inert2 = TAPstruct[[2]]$options$name
  TAPstruct$parameters$experiment$reactant2 = TAPstruct[[2]]$options$name
  TAPstruct$parameters$experiment$products2 = TAPstruct[[2]]$options$name
  return(TAPstruct)
}
