library(methods)
library(glmnet)
library(TAPcode)

#filePath = "/home/adam/research_medford/python_code/tap_code/csv_input_sim/eley_eluc_folder/"#!!!

#Start of the simplest RRM script

#list of reactants generated from the tap simulator
#number of gasses and number of reactants fed

simpleRRM = function(filePath,reactants,feed){
  as.list(strsplit(reactants, ","))
  listOfFiles = list.files(filePath)#!!!
  # Process the inert (read and store the data)
  inert = read.csv(paste0(filePath, "inert.csv"), header = F)#!!!
  timeVector = inert[,1]
  inert = inert[,-1]
  if(is.null(dim(inert))){
    inert = as.data.frame(inert)
  }
  
  # Process the reactants (read and store the data)
  reactant = read.csv(paste0(filePath, "reactant.csv"), header = F)#!!!
  reactant = reactant[,-1]
  if(is.null(dim(reactant))){
    reactant = as.data.frame(reactant)
  }
  
  # Process the products (read and store the data)
  product = list()#will go through folder and look for files with Pro in them #!!!
  findProduct = grep("pro",substr(listOfFiles,1,3))
  for(i in 1:length(findProduct)){
    tempProd = read.csv( paste0(filePath, listOfFiles[findProduct[i]]), header = F)[,-1]
    if(is.null(dim(tempProd))){
      tempProd = as.data.frame(tempProd)
    }
    product[[i]] = tempProd
  }
  
  # Y-Procedure parameters (read in reactor / species parameters for the y-procedure & RRM)
  yParameters = read.csv(paste0(filePath, "parameters_opt.csv")) #!!!
  
  # Requirements
  inertZone1 = (yParameters$reactorLength - yParameters$catalystLength) / 2
  bedPorosity = .4
  lCat = yParameters$catalystLength # length of catalyst
  crossArea = .001140092 # cross sectional area
  molecules = 1.63e-09 * 6.022e23
  reactorLength = yParameters$reactorLength
  
  # Diffusion Calculations 
  inertResidenceTime = mean(inert[,1] * timeVector * max(timeVector))
  diffusion =  bedPorosity * reactorLength^2 / (inertResidenceTime)
  reactantDiffusion = diffusion * sqrt(yParameters$inertAMU / yParameters$reactantAMU)
  productDiffusion = list()
  productLocation = grep("pro", names(yParameters))
  for(i in 1:length(productLocation) ){
    productDiffusion[[i]] =  diffusion * sqrt(yParameters$inertAMU / yParameters[1,productLocation[i]]) 
  }
  
  resultMatrix = matrix( rep(0,length(timeVector) * (1 + length(product))),  length(timeVector), (1 + length(product)) )

  if(dim(inert)[2] == 100){# if the number of pulses is equal to 1?
    # gas concentration and interaction terms with product rate 
    resultMatrix = matrix( rep(0,length(timeVector) * (1 + length(product))),  length(timeVector), (1 + length(product)) )
    
    reactantResults = TAPcode::Yclust(inert = inert, reactant = reactant, timeVector = timeVector, 
                                      molecules = molecules, crossArea = crossArea, 
                                      catalystZone = lCat, bedPorosity = bedPorosity, 
                                      diffusion = reactantDiffusion, catalystWeight = 1, 
                                      inertZone1 = inertZone1, smoothing = 3, isProduct = F)
    
    reactantRate = reactantResults$reactionRate[,1]
    
    resultMatrix[,1] = reactantResults$gasCon[,1]
    
  
    for(i in 1:length(product)){
      productResults = TAPcode::Yclust(inert = inert, reactant = product[[i]], timeVector = timeVector, 
                                       molecules = molecules, crossArea = crossArea, 
                                       catalystZone = lCat, bedPorosity = bedPorosity, 
                                       diffusion = productDiffusion[[i]], catalystWeight = 1, 
                                       inertZone1 = inertZone1, smoothing = 3, isProduct = T)
      resultMatrix[,(i+1)] = productResults$reactionRate[,1] * resultMatrix[,1]
    }
    
    # No variable selection yet, just gas concentration of reactant and the interactions of the product rate and reactant gas concentrations
    fit = glmnet::glmnet(resultMatrix, reactantRate, alpha = 0, lambda = 1e-3)
    # lambda also may be done through cross validation ?cv.glmnet
    result = as.numeric(fit$beta) # gas concentration, product1, product2 ...
    
  }else{#If there is more than one pulse
    #rep() = total number of cells, timevector = nrow, (1 + length(product)) = ncol

    resultMatrix = matrix( rep(0,length(timeVector) * (1 + length(product))),  length(timeVector), (1 + length(product)) )

    reactantResults = TAPcode::Yclust(inert = inert[1], reactant = reactant, timeVector = timeVector, 
                                      molecules = molecules, crossArea = crossArea, 
                                      catalystZone = lCat, bedPorosity = bedPorosity, 
                                      diffusion = reactantDiffusion, catalystWeight = 1, 
                                      inertZone1 = inertZone1, smoothing = 3, isProduct = F)
        
    reactantRate = reactantResults$reactionRate
    reactantGas = reactantResults$gasCon
    reactantSurf = reactantResults$surfaceCon
   
    productInteractions = list()
    for(i in 1:length(product)){
      productResults = TAPcode::Yclust(inert = inert, reactant = product[[i]], timeVector = timeVector, 
                                       molecules = molecules, crossArea = crossArea, 
                                       catalystZone = lCat, bedPorosity = bedPorosity, 
                                       diffusion = productDiffusion[[i]], catalystWeight = 1, 
                                       inertZone1 = inertZone1, smoothing = 3, isProduct = T)
      productInteractions[[i]] = productResults$reactionRate * reactantRate
    }
    
    
    resultMatrix[,1] = reactantResults$gasCon[,1]
    for(i in 1:length(product)){
      productResults = TAPcode::Yclust(inert = inert, reactant = product[[i]], timeVector = timeVector, 
                                       molecules = molecules, crossArea = crossArea, 
                                       catalystZone = lCat, bedPorosity = bedPorosity, 
                                       diffusion = productDiffusion[[i]], catalystWeight = 1, 
                                       inertZone1 = inertZone1, smoothing = 3, isProduct = T)
      
    }
    
    resultMatrix = matrix(rep(0,dim(inert)[1] * (1 + length(product)) ), dim(inert)[1], (1 + length(product)))
    # gas concentration and interaction terms with product rate 
    tempMatrix = matrix( rep(0,dim(inert)[2] * (1 + length(product))),  dim(inert)[2], (1 + length(product)) ) #fill a new matrix with zeros
    # Number of pulses to processed -> print(dim(inert)[2])
    # Number of gas species -> print((1 + length(product)))
    
    #Construct new matrix here!
    for(i in 1:dim(inert)[1]){# for each time step, do process
      response = reactantRate[i,]
      #print(tempMatrix)
      tempMatrix[,1] <- t(reactantGas[i,])#tempMatrix[,2]#1234567123456712345671:1234567123456712345675
      #print(tempMatrix)
      for( j in 1:length(product)){
        tempMatrix[,(j+1)] = t(productInteractions[[j]][i,])#t(productInteractions[[j]][i,])
      }
      print(tempMatrix[,1])
      #print(response)
      
      fit = glmnet::glmnet(tempMatrix, t(response), alpha = 0, lambda = 1e-3)
      resultMatrix[i,] = as.numeric(fit$beta)
      print(fit$beta)
      print(coef(fit,s=0.9))
      #print(resultMatrix)
      #print(coef(fit,s=0.9))
      #print(fit)
      #print(tempMatrix)
      Sys.sleep(1)
     # quit()
    }
    
    # No variable selection yet, just gas concentration of reactant and the interactions of the product rate and reactant gas concentrations
    # lambda also may be done through cross validation ?cv.glmnet
    result = resultMatrix # gas concentration, product1, product2 ...
  }
  

  write.csv(result, paste0(filePath, "rrmResults.csv"))
  return(result)
}
args <- commandArgs(trailingOnly = TRUE)
print(args[1])
test = simpleRRM(args[1],args[2],args[3])#filePath



