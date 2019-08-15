#' Plot Pulses
#'
#' Plots pulses in 2D, 3D waterfall or sequence based on users choice.
#'
#' This function combines pulses in a data frame that can be used graphically.
#' Several inputs are based on the user, for example, normalization of the pulses, color ramping and axis labels.
#' Options on the type of plot include "2D", "3D", "Sequence" and "Heatmap"
#'
#'
#' @param xPulses A list of data frames containing the pulses to be plotted.
#' @param xTime A list of vectors of time values for xPulses.
#' @param yPulses An optional data frame for values to be plotted against xPulses.
#' @param yTime An optional time for yPulses.  This must be present if yPulses are used.
#' @param normalization A string that indicates the normalization type.  Options: "none", "height", "area".
#' @param newXlab The x axis label.
#' @param newYlab The y axis label.
#' @param colRamp A boolean value if the pulses are colored by a single color or matlab like gradient.
#' @param showLegend A boolean value to display the legend.
#' @param plotType A character describing the plot type.  Options include: "2D", "3D", "Sequence" and "Heatmap"
#' @return A plot of the pulses.
#' @examples
#' data("pulseData")
#'
#' xPulses = pulseData$`32AMU`$pulses
#' xTime = pulseData$`32AMU`$time
#' plotPulse(xPulses, xTime)
#'
#'
#'
#'
#'
#'
#' @seealso
#' \code{\link{pulseNorm}}
#' @export plotPulseTesting

plotPulseTesting = function(TAPexperiment, xNames, yNames, timeIndex, pulseIndex, yIndex = NULL, xNorm = "none",yNorm = "none", xTrans = "none", yTrans="none",  newXlab = "", newYlab = "", colRamp = F, showLegend = F, plotType = "2D"){
  lineWidth = 2
  cexSize = 2

  if(xNames == "time"){
    tempXaxis = matrix(rep(TAPexperiment[[2]]$time, dim(TAPexperiment[[2]]$pulses)[2]), length(TAPexperiment[[2]]$time), dim(TAPexperiment[[2]]$pulses)[2])
  }else{
    tempXaxis = TAPexperiment[[xNames]]$pulses
  }

  if(xTrans == "log"){
    tempXaxis[tempXaxis <= 0] = 1e-5
    tempXaxis = log(tempXaxis)
  }else if(xTrans == "log10"){
    tempXaxis[tempXaxis <= 0] = 1e-5
    tempXaxis = log10(tempXaxis)
  }else if(xTrans == "sqrt"){
    tempXaxis[tempXaxis <= 0] = 1e-5
    tempXaxis = sqrt(tempXaxis)
  }

  if(xNorm == "height"){
    tempValues = apply(tempXaxis,2,max)
    tempXaxis = sweep(tempXaxis, 2, tempValues, FUN = "/")
  }else if(xNorm == "area"){
    if(xNames != "time"){
      tempXaxis = apply(tempXaxis, 2, TAPexperiment[[xNames]]$moments$M0, FUN = "/")
    }
  }

  for(i in 1:length(yNames)){
    TAPobj = TAPexperiment[[yNames[i]]]
    tempData = TAPobj$pulses
    if(yTrans == "log"){
      tempData[tempData <= 0] = 1e-5
      tempData = log(tempData)
    }else if(yTrans == "log10"){
      tempData[tempData <= 0] = 1e-5
      tempData = log10(tempData)
    }else if(yTrans == "sqrt"){
      tempData[tempData <= 0] = 1e-5
      tempData = sqrt(tempData)
    }

    if(yNorm == "height"){
      tempValues = apply(tempData,2,max)
      tempData = sweep(tempData, 2, tempValues, FUN = "/")
    }else if(yNorm == "area"){
      tempData = apply(tempData, 2, TAPobj$moments$M0, FUN = "/")
    }

    TAPobj$pulses = tempData
    TAPexperiment[[TAPobj$options$Name]] = TAPobj
  }








  pulseNames = unlist(lapply(xPulses,colnames),F,F)


  pulseIndex = 1
  timeStartPos = which.min(abs((TAPobjects[[1]])$time - timeStart))
  timeEndPos = which.min(abs((TAPobjects[[1]])$time - timeEnd))
  timeRange = timeStartPos:timeEndPos
  minTimes = rep(0, length(TAPobjects))
  maxTimes = rep(0, length(TAPobjects))
  minPulses = rep(0, length(TAPobjects))
  maxPulses = rep(0, length(TAPobjects))

  for(i in 1:length(TAPobjects)){
    TAPobj =  TAPobjects[[i]]

    if(length(pulseIndex) == 1){
      tempPulses = TAPobj$pulses[timeRange, c(pulseIndex, 1)]
      tempMax = TAPobj$moments$max[c(pulseIndex, 1)]
      tempM0 = TAPobj$moments$M0[c(pulseIndex, 1)]
    }else{
      tempPulses = TAPobj$pulses[timeRange, pulseIndex]
      tempMax = TAPobj$moments$max[pulseIndex]
      tempM0 = TAPobj$moments$M0[pulseIndex]
    }
    if(transformation == "height"){
      tempPulses = sweep(tempPulses, 2, tempMax, FUN = "/")
    }else if(transformation == "area"){
      tempPulses = sweep(tempPulses, 2, tempM0, FUN = "/")
    }else if(transformation == "log10"){
      tempPulses[tempPulses <= 0 ] = 1e-10
      tempPulses =  log10(tempPulses)
    }else if(transformation == "log"){
      tempPulses[tempPulses <= 0 ] = 1e-10
      tempPulses =  log(tempPulses)
    }else if(transformation == "sqrt"){
      tempPulses[tempPulses < 0 ] = 1e-10
      tempPulses =  sqrt(tempPulses)
    }
    minTimes[i] = TAPobj$options$timeStart
    maxTimes[i] = TAPobj$options$timeEnd
    minPulses[i] = min(TAPobj$pulses)
    maxPulses[i] = max(TAPobj$pulses)
    (TAPobjects[[i]])$pulses = tempPulses
  }

  ### colors
  if(plotType != "Heatmap"){
    if(colRamp){
      pal = grDevices::rainbow( length(TAPobjects) * length(pulseIndex) + 1)
    }else{
      palfun1 = colorRampPalette(c("royalblue", "paleturquoise"))
      palfun2 = colorRampPalette(c("seagreen", "springgreen"))
      palfun3 = colorRampPalette(c("gold4","yellow"))
      palfun4 = colorRampPalette(c("darkred", "pink"))
      palTracker = 1
      pal = "0"
      for(i in 1:length(TAPobjects)){
        tempPalFun = switch(palTracker,
                            palfun1,
                            palfun2,
                            palfun3,
                            palfun4)
        pal = c(pal,tempPalFun(pulseIndex))

        if(palTracker == 4){
          palTracker = 1
        }else{
          palTracker = palTracker + 1
        }
      }
      pal = pal[-1]
    }
  }





  if(plotType == "Heatmap"){
    catMat = NULL
    for(i in 1:length(TAPobjects)){
      catMat = rbind(catMat, t(as.matrix((TAPobjects[[i]])$pulses)))
    }
    heatmap(catMat, Rowv = NA, Colv = NA, scale = "none", xlab = newXlab, ylab = newYlab, col = colorRamps::matlab.like(256))
  }else{
    #tempColor = pal
    if(is.null(yIndex)){
      minXaxis = min(minTimes)
      maxXaxis = max(maxTimes)
      minYaxis = min(minPulses)
      maxYaxis = max(maxPulses)
    }else{
      minXaxis = minPulses[yIndex]
      maxXaxis = maxPulses[yIndex]
      minYaxis = min(minPulses[-yIndex])
      maxYaxis = max(maxPulses[-yIndex])
    }

    if(is.null(yIndex)){
      numYpoints = 0
      for(i in 1:length(TAPobjects)){
        tempDim = dim((TAPobjects[[i]])$pulses)
        numYpoints = numYpoints + tempDim[1] * tempDim[2]
      }
    }else{
      numYpoints = sum(apply((TAPobjects[[yIndex]])$pulses,2,max)) * length(TAPobjects)
    }



    if(showLegend){
      legScale = .15
    }else{
      legScale = 0
    }

    if(plotType != "Sequence"){
      plot(c(minXaxis,maxXaxis + (maxXaxis - minXaxis)*legScale), c(minYaxis,maxYaxis),
           type = "n", xlab = newXlab, ylab = newYlab, cex.lab = cexSize, cex.axis = cexSize, cex.main = cexSize, cex.sub = cexSize)
    }else{
      plot(c(1,numYpoints +  numYpoints * legScale), c(minYaxis,maxYaxis),
           type = "n", xlab = newXlab, ylab = newYlab, cex.lab = cexSize, cex.axis = cexSize, cex.main = cexSize, cex.sub = cexSize)
    }


    if(is.null(yIndex)){
      palFlag = 1
      timeFlag = 1
      if(plotType == "3D"){
        for(i in 1:length(TAPobjects)){
          TAPobj = TAPobjects[[i]]
          for(j in 1:dim(TAPobj$pulses)[2]){
            subPal = colorRampPalette(c(pal[palFlag], pal[palFlag + 1]))
            lines((TAPobj$time[timeRange] / 2 + (palFlag)/length(pal) * maxXaxis / 2  ),
                   (TAPobj$pulses[,j] / 2  + (palFlag)/length(pal) * maxYaxis/2 ),
                   col =   subPal(length(timeRange)), lwd = lineWidth)
            palFlag = palFlag + 1
          }
        }
      }else if(plotType == "Sequence"){
        for(i in 1:length(TAPobjects)){
          TAPobj = TAPobjects[[i]]
          for(j in 1:dim(xPulses[[i]])[2]){
            lines( seq(timeFlag,(timeFlag + length(xTime[[i]])-1)), xPulses[[i]][,j], col = pal[palFlag], lwd = lineWidth)
            timeFlag = timeFlag + length(xTime[[i]])
            palFlag = palFlag + 1
          }
        }
      }else{
        for(i in 1:length(xPulses)){
          for(j in 1:dim(xPulses[[i]])[2]){
            lines(xTime[[i]], xPulses[[i]][,j], col = pal[palFlag], lwd = lineWidth)
            palFlag = palFlag + 1
          }
        }
      }

    }else{
      palFlag = 1
      timeFlag = 0
      if(plotType == "3D"){
        for(i in 1:length(xPulses)){
          for(j in 1:dim(xPulses[[i]])[2]){
            minLen = min(c(length(yPulses[,j]), length(xPulses[[i]][,j])))
            lines( (yPulses[1:minLen,j]/2 + (palFlag)/length(pal) * maxXaxis/2  ),
                   (xPulses[[i]][1:minLen,j]/2  + (palFlag)/length(pal) * maxYaxis/2 ), col = pal[palFlag], lwd = lineWidth)
            palFlag = palFlag + 1
          }
        }
      }else if(plotType == "Sequence"){
        for(i in 1:length(xPulses)){
          for(j in 1:dim(xPulses[[i]])[2]){
            minLen = min(c(length(yPulses[,j]), length(xPulses[[i]][,j])))
            lines(yPulses[1:minLen,j] + timeFlag  , xPulses[[i]][1:minLen,j], col = pal[palFlag], lwd = lineWidth)
            timeFlag = timeFlag + max(yPulses[,j])
            palFlag = palFlag + 1
          }
        }
      }else{
        for(i in 1:length(xPulses)){
          for(j in 1:dim(xPulses[[i]])[2]){
            minLen = min(c(length(yPulses[,j]), length(xPulses[[i]][,j])))
            lines(yPulses[1:minLen,j], xPulses[[i]][1:minLen,j], col = pal[palFlag], lwd = lineWidth)
            palFlag = palFlag + 1
          }
        }
      }
    }

    if(plotType == "seq" ){
      legMax = numYpoints
    }else{
      legMax = maxXaxis
    }

    if(showLegend){
      legend("topright",
             pulseNames,
             lty=rep(1,length(pulseNames)),
             lwd=rep(2.5,length(pulseNames)),
             col=pal)
    }
  }




}
