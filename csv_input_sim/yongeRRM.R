library(methods)
library(glmnet)
library(TAPcode)

simpleRRM = function(filePath,reactants,feed){
  
  ### Read in outlet flux data based on a list of reactants and products specified and fed to the script
  listOfFiles = as.list(strsplit(reactants, ',')[[1]])
  inert = read.csv(paste0(filePath, "flux_data/Inert.csv"), header = F)#!!!
  
  ### Store time and inert data 
  timeVector = inert[,1]
  inert = inert[,-1]
  if(is.null(dim(inert))){
    inert = as.data.frame(inert)
  }
  
  ### Used later to setup the linear model correctly
  x_sum<-seq(0,length(listOfFiles)-1,1)
  result2 <- sum((length(listOfFiles)-x_sum))
  total <- result2+length(listOfFiles)*length(listOfFiles)+length(listOfFiles)*length(listOfFiles)
  param_matrix <- array(rep(0, length(listOfFiles)*(total+1)*dim(inert)[2]), dim=c(length(listOfFiles), 1+total, dim(inert)[2]))
  
  ### for each pulse found in the inert data file, run the process
  for (k_pulse in 1:dim(inert)[2]){
  
    ### Generate the vector/names and store the data for all gas species alongside each name
    all_vect = c()
    my_names = c()
    for(i in 1:(length(listOfFiles))){#strtoi(feed)
      my_names = c(my_names,listOfFiles[i])
      all_vect = c(all_vect,read.csv(paste0(filePath, paste0('flux_data/',paste(listOfFiles[i],'.csv',sep=''))), header = F)[k_pulse+1])
    }
    names(all_vect) = my_names
    
    ### Read in the parameters used for the simulator for the yprocedure
    yParameters = read.csv(paste0(filePath, "RRM_params.csv")) #!!!
    
    ### Convert the parameters for use in the yprocedure
    inertZone1 = (yParameters$length_reac - yParameters$length_reac * yParameters$factor) / 2
    bedPorosity = yParameters$void_inert  #!!!
    lCat = yParameters$length_reac * yParameters$factor # length of catalyst
    crossArea = 3.14159*(lCat)^2 #.001140092 # cross sectional area #!!!
    molecules = yParameters$Inert_pulse_size
    reactorLength = yParameters$length_reac
    mass_list_temp = toString(yParameters$mass_list)
    mass_list = as.list(strsplit(mass_list_temp, ',')[[1]])
    
    ### Calculate the diffusion coefficients 
    inertResidenceTime = mean(inert[,1] * timeVector * max(timeVector))
    diffusion =  bedPorosity * reactorLength^2 / (inertResidenceTime)

    allDiffusion = c()
    my_names = c()
    for(i in 1:(length(listOfFiles))){#strtoi(feed)
      my_names = c(my_names,listOfFiles[i])
      allDiffusion = c(allDiffusion, diffusion * sqrt(as.double(mass_list[[length(mass_list)]]) / as.double(mass_list[[i]]  )))
    }
    
    names(allDiffusion) = my_names
    
    ### generate 'dictionaries' to store the y-procedure data
    gas_vect = matrix(nrow=dim(inert[k_pulse])[1], ncol = length(listOfFiles))
    surf_vect = matrix(nrow=dim(inert[k_pulse])[1], ncol = length(listOfFiles))
    rate_vect = matrix(nrow=dim(inert[k_pulse])[1], ncol = length(listOfFiles))
    
    for(i in 1:(length(listOfFiles))){
      
      ### If it is a product or reactant
      if(i < strtoi(feed)){
        new_step = T
      }
      else{
        new_step = F
      }
      reactantResults = TAPcode::Yclust(inert = inert[k_pulse], reactant = matrix( all_vect[[listOfFiles[[i]]]] , length(all_vect[[listOfFiles[[i]]]]), 1 ), timeVector = timeVector, 
                                        molecules = molecules, crossArea = crossArea, 
                                        catalystZone = lCat, bedPorosity = bedPorosity, 
                                        diffusion = allDiffusion[[listOfFiles[[i]]]], catalystWeight = 1, 
                                        inertZone1 = inertZone1, smoothing = 3, isProduct = new_step)
      
      gas_vect[,i] = data.matrix(reactantResults$gasCon[,1])
      surf_vect[,i] = data.matrix(reactantResults$surfaceCon[,1])
      rate_vect[,i] = data.matrix(reactantResults$reactionRate[,1])

    }
    
    ### generate empty matricies with rate data and c1, c2, u1, u2, c1u1 ... data
    resultMatrix = matrix(rep(0,dim(inert[k_pulse])[1] * total ), dim(inert[k_pulse])[1], total)
    tempMatrix = matrix( rep(0,dim(inert[k_pulse])[2] * total),  dim(inert[k_pulse])[2], total )
    
    for(i in 1:length(listOfFiles)){# for each time step, do process
      resultMatrix[,i] = gas_vect[,i]
    }
    for(i in 1:length(listOfFiles)){# for each time step, do process
      resultMatrix[,i+length(listOfFiles)] = surf_vect[,i]
    }
    
    for(i in 1:(length(listOfFiles))){# for each time step, do process
      for(j in 1:length(listOfFiles)){
        resultMatrix[,j+(i-1)*length(listOfFiles)+2*length(listOfFiles)] = gas_vect[,i]*surf_vect[,j]
      }
    }
    
    for(i in 1:(length(listOfFiles))){
      for(j in 1:(length(listOfFiles)-(i-1))){
        resultMatrix[,j+(i-1)*length(listOfFiles)+2*length(listOfFiles)+length(listOfFiles)*length(listOfFiles)] = surf_vect[,i]*surf_vect[,j+(i-1)]
      }
    }
    
    ### for each gas (except the inert), use glmnet to fit the parameters  
    for(j_p in 1:(length(listOfFiles))){
      fit = glmnet::glmnet(resultMatrix, rate_vect[,i], alpha = 0, lambda = 1e-3)
      param_matrix[j_p,,k_pulse] = coef(fit,s=c(0.1))[,1]
    }
    write.csv(resultMatrix, paste0(paste0(filePath,'RRM_results/'), paste0(k_pulse,"_y_proc.csv")))
  }
  for(j_p in 1:(length(listOfFiles))){
    write.csv(t(param_matrix[j_p,,]), paste0(paste0(filePath,'RRM_results/'), paste0(listOfFiles[[j_p]],"_reactivities.csv")))
  }
  
  }
args <- commandArgs(trailingOnly = TRUE)
test = simpleRRM(args[1],args[2],args[3])
quit()
  
  
############################################################################################################################################
###################################     Extra code -- could be useful later, but can ignore for now      ###################################
############################################################################################################################################

#Generate the vector/names and store the data for each reactant
#react_vect = c()
#my_names = c()
#for(i in 1:(strtoi(feed))){#strtoi(feed)
#  my_names = c(my_names,listOfFiles[i])
#  react_vect = c(react_vect,read.csv(paste0(filePath, paste0('flux_data/',paste(listOfFiles[i],'.csv',sep=''))), header = F)[k_pulse+1])
#  #print(paste(listOfFiles[i],'.csv',sep=''))
#}
#names(react_vect) = my_names

#Generate the vector/names and store the data for each product
#prod_vect = c()
#my_names = c()
#for(i in (strtoi(feed)+1):(length(listOfFiles))){#strtoi(feed)
#  my_names = c(my_names,listOfFiles[i])
#  prod_vect = c(prod_vect,read.csv(paste0(filePath, paste0('flux_data/',paste(listOfFiles[i],'.csv',sep=''))), header = F)[k_pulse+1])
#  #print(paste(listOfFiles[i],'.csv',sep=''))
#}
#names(prod_vect) = my_names

#reactantDiffusion = c()
#my_names = c()
#for(i in 1:(strtoi(feed))){
#  my_names = c(my_names,listOfFiles[i])
#  reactantDiffusion = c(reactantDiffusion, diffusion * sqrt(as.double(mass_list[[length(mass_list)]]) / as.double(mass_list[[i]]  )))
#  #print(paste(listOfFiles[i],'.csv',sep=''))
#}
#names(reactantDiffusion) = my_names

#productDiffusion = c()
#my_names = c()
#for(i in (strtoi(feed)+1):(length(listOfFiles))){#strtoi(feed)
#  my_names = c(my_names,listOfFiles[i])
#  productDiffusion = c(productDiffusion, diffusion * sqrt(as.double(mass_list[[length(mass_list)]]) / as.double(mass_list[[i]]  )))
#  #print(paste(listOfFiles[i],'.csv',sep=''))
#}
#names(productDiffusion) = my_names
